"""The one place that knows PolyFrame is an external Windows program.

Writes a job file, runs the C# bridge in ``bridge/`` with ``dotnet run`` and reads
the result file back. One application start per operation (``Calculate``, then
``SavePDF``) because a PolyFrame connection works for its first call only; a case
costs about 15 s. Nothing in the test suite imports this module - see
``test_polyframe_buckling.py::test_tests_never_start_polyframe``.
"""

import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

DEFAULT_INSTALLATION_DIRECTORY = Path(
    os.environ.get("POLYFRAME_INSTALL_DIR", r"C:\Program Files\PolyStruc\PolyFrame")
)
BRIDGE_PROJECT = Path(__file__).resolve().parent / "bridge" / "PolyFrameBridge.csproj"


class PolyFrameUnavailableError(RuntimeError):
    """PolyFrame or the .NET SDK is not usable on this machine."""


@dataclass
class CaseJob:
    case_id: str
    model: dict[str, Any]
    pdf_path: Path | None = None


@dataclass
class CaseOutcome:
    case_id: str
    calculated: bool
    error: str | None
    result: dict[str, Any] | None
    pdf_path: str | None


@dataclass
class RunOutcome:
    connector_version: str
    cases: list[CaseOutcome]
    stdout: str = ""
    stderr: str = ""
    exit_code: int = 0

    def outcome_for(self, case_id: str) -> CaseOutcome:
        for case in self.cases:
            if case.case_id == case_id:
                return case
        raise KeyError(case_id)


def _dotnet() -> str:
    found = shutil.which("dotnet")
    if found is None:
        raise PolyFrameUnavailableError(
            "the .NET SDK is not on PATH; the bridge targets net9.0."
        )
    return found


def check_available(installation_directory: Path = DEFAULT_INSTALLATION_DIRECTORY) -> None:
    if sys.platform != "win32":
        raise PolyFrameUnavailableError(
            f"PolyFrame is a Windows desktop application; this is {sys.platform}."
        )
    _dotnet()
    if not installation_directory.is_dir():
        raise PolyFrameUnavailableError(
            f"PolyFrame is not installed at '{installation_directory}'. Install it "
            "from https://www.PolyStruc.com/Download or set POLYFRAME_INSTALL_DIR."
        )


def run_bridge(
    jobs: list[CaseJob],
    work_dir: Path,
    installation_directory: Path = DEFAULT_INSTALLATION_DIRECTORY,
    timeout_seconds: int = 3600,
) -> RunOutcome:
    check_available(installation_directory)
    work_dir.mkdir(parents=True, exist_ok=True)
    job_path = work_dir / "job.json"
    result_path = work_dir / "results.json"
    if result_path.exists():
        result_path.unlink()

    job_path.write_text(
        json.dumps(
            {
                "installationDirectory": str(installation_directory),
                "cases": [
                    {
                        "caseId": job.case_id,
                        "model": job.model,
                        **({"pdfPath": str(job.pdf_path)} if job.pdf_path else {}),
                    }
                    for job in jobs
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    completed = subprocess.run(  # nosec B603
        [
            _dotnet(), "run", "--project", str(BRIDGE_PROJECT), "-c", "Release", "--",
            str(job_path), str(result_path),
        ],
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
        check=False,
    )
    if not result_path.exists():
        raise PolyFrameUnavailableError(
            "the bridge wrote no result file.\n"
            f"exit code: {completed.returncode}\n"
            f"stdout: {completed.stdout.strip()}\nstderr: {completed.stderr.strip()}"
        )

    payload = json.loads(result_path.read_text(encoding="utf-8"))
    cases = [
        CaseOutcome(
            case_id=c["CaseId"],
            calculated=bool(c.get("Calculated")),
            error=c.get("Error"),
            result=c.get("Result"),
            pdf_path=c.get("PdfPath"),
        )
        for c in payload.get("Cases", [])
    ]
    return RunOutcome(
        connector_version=payload.get("ConnectorVersion", "unknown"),
        cases=cases,
        stdout=completed.stdout,
        stderr=completed.stderr,
        exit_code=completed.returncode,
    )
