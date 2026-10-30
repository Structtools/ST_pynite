// A pass-through between Python and PolyFrame.Connector. It holds no
// engineering knowledge and must not acquire any: the model arrives as JSON
// that Python authored, gets deserialised into PolyFrameDataModel, calculated,
// and serialised back. Anything that decides what a frame *is* belongs in
// tools/polyframe/translate.py, where it can be tested without a licence.
//
//     dotnet run --project tools/polyframe/bridge -- <job.json> <results.json>
//
// The job file:
//
//     {
//       "installationDirectory": "C:\\Program Files\\PolyStruc\\PolyFrame",
//       "cases": [ { "caseId": "steel-ipe300-portal-fixed-udl-wind",
//                    "model": { ...PolyFrameDataModel... },
//                    "pdfPath": "optional\\report.pdf",
//                    "pdfOnly": false } ]
//     }
//
// `pdfOnly` skips Calculate and only writes the report: PolyFrame prints the
// buckling lengths in the PDF and nowhere in the API, so exploring them needs
// no forces.
//
// The result file mirrors it, one entry per case, each either a result or an
// error. A case that throws does not abort the batch — a licence hiccup on
// case 4 should not cost the three captures that already succeeded.
//
// IMPORTANT: PolyFrame is started once per OPERATION, never once per batch. A
// connection works for its first call only; see docs/POLYFRAME_DISCOVERY.md.

using System.Text.Json;
using System.Text.Json.Serialization;
using PolyFrame.Connector;
using PolyFrame.DataModels;
using PolyFrame.DataModels.Results;

namespace PolyFrameBridge;

internal sealed class JobCase
{
    public string CaseId { get; set; } = "";
    public JsonElement Model { get; set; }
    public string? PdfPath { get; set; }
    public bool PdfOnly { get; set; }
}

internal sealed class Job
{
    public string InstallationDirectory { get; set; } = "";
    public List<JobCase> Cases { get; set; } = new();
}

internal sealed class CaseOutcome
{
    public string CaseId { get; set; } = "";
    public bool Calculated { get; set; }
    public string? Error { get; set; }
    public PolyFrameDataModelResult? Result { get; set; }
    public string? PdfPath { get; set; }
}

internal sealed class BatchOutcome
{
    public string BridgeVersion { get; set; } = "";
    public string ConnectorVersion { get; set; } = "";
    public List<CaseOutcome> Cases { get; set; } = new();
}

internal static class Program
{
    // Python authors readable JSON: PascalCase property names and enums as
    // their C# names, which the string-enum converter reads back.
    private static readonly JsonSerializerOptions Options = new()
    {
        WriteIndented = true,
        PropertyNameCaseInsensitive = true,
        DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull,
        Converters = { new JsonStringEnumConverter() },
        // The other PolyStruc programs return Infinity for a utilisation whose
        // capacity is zero. Without this the whole batch fails to serialise and
        // every capture in it is lost. Writing the literal keeps the run intact;
        // tools/polyframe/extract.py then drops non-finite values rather than
        // pretending they are comparable.
        NumberHandling = JsonNumberHandling.AllowNamedFloatingPointLiterals,
    };

    private static int Main(string[] args)
    {
        if (args.Length != 2)
        {
            Console.Error.WriteLine("usage: PolyFrameBridge <job.json> <results.json>");
            return 2;
        }

        Job? job;
        try
        {
            job = JsonSerializer.Deserialize<Job>(File.ReadAllText(args[0]), Options);
        }
        catch (Exception exc)
        {
            Console.Error.WriteLine($"could not read job file '{args[0]}': {exc.Message}");
            return 2;
        }

        if (job is null || job.Cases.Count == 0)
        {
            Console.Error.WriteLine("job file contains no cases");
            return 2;
        }

        if (!Directory.Exists(job.InstallationDirectory))
        {
            Console.Error.WriteLine(
                $"PolyFrame installation directory not found: '{job.InstallationDirectory}'");
            return 2;
        }

        var batch = new BatchOutcome
        {
            BridgeVersion = "1",
            ConnectorVersion =
                typeof(PolyFrameConnection).Assembly.GetName().Version?.ToString() ?? "unknown",
        };

        // One connection per *operation*, not per batch. PolyFrame 1.1.11 breaks
        // after its first call on a connection: a second Calculate returns null
        // and a SavePDF after a Calculate throws "LoadCombination not found in
        // results" (docs/POLYFRAME_DISCOVERY.md). A fresh application per
        // operation costs 2-6 s each and is the only reliable pattern.
        foreach (var jobCase in job.Cases)
        {
            batch.Cases.Add(RunOne(job.InstallationDirectory, jobCase));
        }

        File.WriteAllText(args[1], JsonSerializer.Serialize(batch, Options));

        var failed = batch.Cases.Count(c => !c.Calculated);
        Console.WriteLine($"{batch.Cases.Count - failed}/{batch.Cases.Count} cases calculated");
        // A per-case failure is data, not a crash: the result file records it
        // and the Python side decides what to do. Exit 0 so a partial batch is
        // still readable.
        return 0;
    }

    private static CaseOutcome RunOne(string installationDirectory, JobCase jobCase)
    {
        var outcome = new CaseOutcome { CaseId = jobCase.CaseId };
        try
        {
            if (!jobCase.PdfOnly)
            {
                var result = WithConnection(installationDirectory, jobCase, (polyFrame, model) =>
                    polyFrame.Calculate(model));
                if (result is null)
                {
                    // PolyFrame returns null rather than throwing when it cannot
                    // calculate, so this is a normal outcome to record.
                    outcome.Error = "Calculate returned null";
                    return outcome;
                }

                outcome.Result = result;
                outcome.Calculated = true;
            }

            if (!string.IsNullOrWhiteSpace(jobCase.PdfPath))
            {
                // The PDF is the only place the buckling lengths and every
                // Eurocode intermediate exist, so it is worth a second
                // application start.
                Directory.CreateDirectory(Path.GetDirectoryName(jobCase.PdfPath)!);
                var saved = WithConnection(installationDirectory, jobCase, (polyFrame, model) =>
                    polyFrame.SavePDF(model, jobCase.PdfPath!));
                if (saved)
                {
                    outcome.PdfPath = jobCase.PdfPath;
                    outcome.Calculated |= jobCase.PdfOnly;
                }
            }
        }
        catch (Exception exc)
        {
            outcome.Error = $"{exc.GetType().Name}: {exc.Message}";
        }

        return outcome;
    }

    private static T WithConnection<T>(
        string installationDirectory,
        JobCase jobCase,
        Func<PolyFrameConnection, PolyFrameDataModel, T> operation)
    {
        var model = JsonSerializer.Deserialize<PolyFrameDataModel>(
            jobCase.Model.GetRawText(), Options)
            ?? throw new InvalidOperationException("model JSON deserialised to null");

        PolyFrameConnection? polyFrame = null;
        try
        {
            polyFrame = new PolyFrameConnection(installationDirectory);
            if (!polyFrame.IsConnected)
            {
                throw new InvalidOperationException("connected to PolyFrame but IsConnected is false");
            }

            return operation(polyFrame, model);
        }
        finally
        {
            // Leaving a WPF process behind would make the next run connect to a
            // stale instance. CloseApplication is what actually stops it;
            // Dispose alone does not.
            try { polyFrame?.CloseApplication(); } catch { /* already gone */ }
            polyFrame?.Dispose();
        }
    }
}
