# PolyFrame buckling comparison

Compares PyNite's linear buckling analysis (`FEModel3D.analyze_buckling`, this
branch) with PolyFrame 1.1.11's **general method** (EN 1993-1-1 6.3.4), which
reports for a planar steel frame

- the **in-plane** critical load factor `alpha_cr` (EN 1993-1-1 5.2.1), and
- the **out-of-plane** critical load factor `alpha_cr,op` for lateral and
  lateral-torsional buckling (6.3.4), with the governing element.

PyNite computes the first with `analyze_buckling(plane='XY')` and the second with
`analyze_buckling(plane='XY', out_of_plane=True)`. The margin is **±5 %** on the
factor. Results: [`docs/POLYFRAME_BUCKLING_COMPARISON.md`](../../docs/POLYFRAME_BUCKLING_COMPARISON.md).

**How to read the gate.** "N passed, 4 xfailed" is not "all within ±5 %". The four
xfails are cases where PyNite is *above* PolyFrame by more than 5 % (the unconservative
direction) and a `KnownDeviation` with a stated cause bounds how far: 10 % on three
fixed-base portals, and 100 % on the top-flange-loaded portal, which is a watch, not a
gate. Bounds are one-sided: PyNite falling below PolyFrame beyond 5 % fails. An
independent third reference (`ltb_reference.py`, below) is what adjudicates them.

**Licence.** The corpus was captured under a student licence; see `corpus/README.md`
for what that means before project use or publication.

The FEM Design references that the `*_comparison.py` scripts in `Testing/` carry
are gated the same way by `Testing/test_fem_design_buckling.py` (in-plane only;
FEM Design's out-of-plane factors were not captured).

## Layout

| File | Role |
|---|---|
| `cases.py` | the catalogue: one `FrameCase` per frame, rendered into both tools |
| `sections.py` | IPE200 / IPE300 / HE200A properties (SI) and S355 |
| `translate.py` | `FrameCase` -> PolyFrame data model (JSON for the bridge) |
| `pynite_side.py` | `FrameCase` -> `FEModel3D`, and the three buckling runs |
| `bridge/` | C# pass-through to `PolyFrame.Connector` (one application start per operation) |
| `runner.py` | runs the bridge with `dotnet run`; the only module that starts PolyFrame |
| `pdf_report.py` | reads `alpha_cr`, `alpha_cr,op`, the governing element and the utilisations off the PDF |
| `ltb_reference.py` | independent Rayleigh-Ritz lateral-torsional buckling solver (Vlasov theory, true warping energy); tested in `../test_ltb_reference.py` |
| `reference_checks.py` | the governing member of each open case in that solver, against both tools |
| `capture.py`, `corpus.py` | bridge outcome + PDF -> `corpus/<case>.json` |
| `summary.py` | writes the results document |
| `__main__.py` | the CLI |
| `corpus/` | the captured records; `corpus/reports/` the PDFs they were read from (gitignored) |
| `../test_polyframe_buckling.py` | the gate |

## Running

```bash
uv run --no-sync python -m pytest Testing/test_polyframe_buckling.py Testing/test_fem_design_buckling.py
```

pytest never starts PolyFrame (a test checks that). Capturing needs Windows, the
PolyFrame application at `C:\Program Files\PolyStruc\PolyFrame` (or
`POLYFRAME_INSTALL_DIR`), the .NET 9 SDK, a licence signed in once by hand, and
the `polyframe` dependency group (PyMuPDF, Pillow, RapidOCR):

```bash
uv sync --group polyframe
uv run --no-sync python -m Testing.polyframe list                 # catalogue and capture state
uv run --no-sync python -m Testing.polyframe refresh              # recapture everything (~1 min per case)
uv run --no-sync python -m Testing.polyframe refresh --case portal-ipe300-pinned-100-wind
uv run --no-sync python -m Testing.polyframe reparse              # re-read saved PDFs, no PolyFrame
uv run --no-sync python -m Testing.polyframe pynite               # PyNite's factors
uv run --no-sync python -m Testing.polyframe summary              # docs/POLYFRAME_BUCKLING_COMPARISON.md
```

Review the diff of `corpus/` after a `refresh`: every changed number is a changed
reference value. A record carries the SHA-256 of the PolyFrame model it was captured
from, and the test fails if the catalogue entry no longer hashes to it.

## Where PolyFrame's numbers come from

PolyFrame's API returns sectional forces, reactions and aggregate utilisations only.
Both critical load factors exist in the PDF report and nowhere else:

- `alpha_cr` is on the **text layer** when it is above 10 ("The in-plane critical
  load factor 17.2 is larger than 10"). Below 10 PolyFrame runs a second-order
  analysis with an imperfection and the factor is inside a rendered formula image.
- `alpha_cr,op` is always inside images: the formula `alpha_cr,op = 5.63` in the
  stability check and the captions of the two buckling-shape plots in the appendix
  ("Critical load factor: alpha_cr = 5.63").

`pdf_report.py` renders each image region at 300 and 500 dpi with a white margin
and reads it with RapidOCR, which reads the LaTeX-style glyphs reliably; the formula
and the two plot captions are cross-checked and a disagreement is written to the
record's `messages`. The record says which source each value came from (`text`,
`ocr`, or `hand` for a value typed into `corpus/hand_read.json` from the report).

## What is being compared, and the modelling assumptions

- **Material and sections.** S355, E = 210 GPa, G = E / 2.6. Sections by name in
  PolyFrame; tabulated catalogue values in PyNite (`sections.py`), including the
  warping constant and the EN 1993-1-1 6.2.6 shear areas. Both tools use
  Timoshenko members.
- **Supports.** In plane, PolyFrame's three DOFs. Out of plane, PolyFrame's general
  method restrains the lateral displacement at **every node** and nothing else
  (its report prints a "Stability Supports" table with FZ at every node and the
  MX / MY columns empty, fixed bases included). PyNite does the same: DZ restrained
  at every node, RX and RY free. A single straight member would be a twist
  mechanism under that rule, so such cases get fork supports (twist restrained at
  the supported nodes).
- **Warping.** PyNite's 6-DOF element has no warping DOF; `warping='equivalent_torsion'`
  uses `J + pi^2 E Iw / (G L_w^2)` with the member's **warping length** `L_w`
  (`FEModel3D.def_warping_length`), by default the member length, which assumes the
  twist is held at both member ends. PolyFrame holds no twist at a support, so a column
  on a base, or a span of a continuous beam at an intermediate support, is twist-free
  at one end. `pynite_side.warping_length_factor` declares `2 L` for such members (a
  member end counts as twist-restrained when a non-collinear member meets it there)
  and `math.inf` for a member held at neither end; the test gates on this variant
  (`VARIANT = "restraints"`), and the summary also reports the member-length and
  no-warping values. Without the rule, columns come out 13-24 % above PolyFrame and a
  two-span beam 39 % above; with it, within 2 %. PyNite now warns when a member on the
  default has an end where neither a support nor another member holds the twist.
  **Status of the `2 L` rule:** it was found by matching PolyFrame. The independent
  solver shows it is not textbook physics: for pure torsional buckling a member held in
  twist at one end with free warping has *no* warping contribution (the critical mode
  is a uniform twist rate), and `L_T = 2 L` belongs to warping fixed at the held end.
  Under bending the true contribution lies between the two, and `2 L` is a working
  approximation for it whose accuracy the reference table below measures per case.
- **Single members need a stub.** PolyFrame's general method refuses a structure whose
  elements are all collinear, so every column and beam case carries an unloaded 1 m
  stub at a right angle. The stub's tip is a node and therefore a lateral restraint,
  which makes it a fork for the member end it hangs from, in both tools.
- **Load height.** A vertical load's attack point `Top` / `Bottom` becomes
  `load_height = +h/2` / `-h/2` in PyNite (EN 1993-1-1 Annex F `z_g` sign).
- **Subdivision.** 8 elements per member.

## Findings (2026-10-02, PolyFrame 1.1.11)

The numbers are in `docs/POLYFRAME_BUCKLING_COMPARISON.md` (including the independent
reference table); deviations are `KnownDeviation` entries on the cases they concern,
so they gate.

1. **In plane, the two tools agree to 0.2 % on every case** (15 frames, factors from
   5.7 to 225), and PyNite agrees with FEM Design on 58 modes within 5 %. Both PolyStruc
   and PyNite include shear deformation.
2. **Out of plane, 18 of 22 cases agree within 5 %** once PyNite's warping
   approximation uses the twist-restraint spacing (above). Pure-compression columns
   agree to 0.1 %, simply supported beams under UDL, point load and top/bottom flange
   loads to 0-4 %, pinned-base portals, the gable and the two-bay frame to 2 %.
3. **The independent reference sides with both tools on pinned-base columns and with
   PyNite on the load height.** The Rayleigh-Ritz Vlasov solver (`ltb_reference.py`),
   fed each governing member's first-order forces, reproduces PolyFrame to 0.1-0.6 % and
   PyNite to 0.8-2.7 % on the four pinned-base columns and the gable column, and the
   two-span beam to 1.5 % (17.23 against 16.96 / 16.97; the member-length value 23.5 is
   wrong). PolyFrame therefore computes an exact member solution for a column with zero
   base moment, and PyNite's `2 L` warping length approximates it to about 2 %.
4. **Fixed-base portals: accepted as safe-side in both tools (decision 2026-10-02).**
   PolyFrame is 6.5-7.6 % below PyNite on the three fixed-base frames. The isolated
   member reference gives 6.46 on the fixed-100 column (PyNite 6.31, PolyFrame 5.88)
   and 4.52 / 6.87 on the rafters (PyNite 4.08 / 6.39, PolyFrame 3.79 / 6.00): both
   tools below it. The column carries a fixing moment at a node that holds no twist,
   a configuration that exists only because PolyFrame's model holds no twist at a
   fixed base; a real fixed base does, and PyNite with the twist held at the bases
   gives 15.1 / 4.50 / 6.92. So the twist-free numbers of both tools are conservative
   for a real fixed base, PolyFrame's more so. Bound 10 %, PyNite high, kept as a watch.
5. **PolyFrame's load-height term inside a frame is the outlier.** On a simply supported
   beam both tools and the reference move alpha_cr,op by -22 % for a top-flange load.
   On the portal rafter the reference gives 2.655, PyNite 2.64 (-0.5 %), PolyFrame 1.44
   (-46 %); a bottom-flange probe raised PolyFrame only 11 %. Bound 100 %, a watch on a
   PolyFrame finding.
6. **Mesh.** Refining from 8 to 32 elements per member moves the four open cases by
   -0.02 % to -1.9 % and the simply supported beam by -1.7 %; the gaps are not
   discretisation.
7. **PolyFrame's "out-of-plane buckling mode: E<n>" is not the critical member.** It
   names a column on the UDL portals while the factor it prints equals PyNite's beam
   mode. Not gated.
8. **The general method cannot model a laterally free end**: every node is a lateral
   restraint, so a cantilever column is pinned-pinned out of plane (its in-plane factor
   is the cantilever's, 53.4 against 53.42). Its 20 kNm fixing moment sits at a
   twist-free node, where both tools agree with each other (5.69 / 5.68) and the member
   reference does not (3.37); see 4.
9. **PyNite's equivalent-torsion warping length is now explicit.**
   `FEModel3D.def_warping_length` sets it per member, the analysis warns when a member
   on the default has a twist-free end, the docstrings no longer claim conservatism,
   and `Testing/test_buckling.py` holds the closed-form torsional case (twist held at
   one end, warping fixed: `L_T = 2 L`). The harness uses the library argument.
