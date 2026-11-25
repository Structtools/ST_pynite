# Board ruling: next steps for PyNite buckling verification

## Follow-up, same day: what was done, what remains

**Done (the board's Must list, plus the Should items that need no licence or external tool).**

| Item | Where | Result |
|---|---|---|
| Per-member warping length in the library, no false conservatism claim, warning on a twist-free end left at the default | `Pynite/FEModel3D.py` (`def_warping_length`, `_twist_held_at`), `Pynite/Member3D.py`, `docs/source/buckling.rst` | harness moved onto the library argument; all PolyFrame numbers unchanged |
| Closed-form torsional test with one twist-free end, and a test that the default fails there | `Testing/test_buckling.py` | 36 passed |
| Mesh rerun at 16 and 32 elements of the four open cases | `Testing/polyframe/README.md` finding 6 | -0.02 to -1.9 %, gaps are not discretisation |
| One-sided gates and the 100 % watch stated as such | `Testing/test_polyframe_buckling.py`, `cases.py`, README | a known deviation now only excuses PyNite being high |
| Licence note on the student-licence corpus | `Testing/polyframe/corpus/README.md` | recapture conditions stated |
| Independent third reference | `Testing/polyframe/ltb_reference.py` (Rayleigh-Ritz Vlasov solver), `test_ltb_reference.py` (7 closed forms), `reference_checks.py`, table in `docs/POLYFRAME_BUCKLING_COMPARISON.md` | see below |

**What the independent reference settled.**

- The `2 L` rule under bending: the Nay Sayer's decider, `beam-ipe300-two-span`, gives 17.23 in the reference against 16.96 (PolyFrame) and 16.97 (PyNite with the rule); the member-length value 23.5 is wrong. Pinned-base portal columns and the gable column: reference within 0.1-0.6 % of PolyFrame and 0.8-2.7 % of PyNite. The rule is a working approximation, not textbook physics: for pure torsional buckling a member held in twist at one end with free warping has *no* warping contribution, and `L_T = 2 L` belongs to warping fixed at the held end. Kill criterion not triggered.
- The top-flange portal: reference 2.655, PyNite 2.64 (-0.5 %), PolyFrame 1.44 (-46 %). PolyFrame is the outlier. Kill criterion ("hand M_cr sides with PolyFrame") not triggered.
- Rafters at the shear centre: the reference (rigid-fork member) sits 5-10 % above PyNite and 6-16 % above PolyFrame; both tools are conservative, PyNite less so, and the rigid-fork upper bound cannot split a 7 % difference between them.

**Not solved, and what has to happen before it can be.**

1. *Fixed-base portals, PolyFrame 6.5-7.6 % below PyNite.* **Decided (user, 2026-10-02): accepted, both tools are on the safe side.** The column carries a fixing moment at a node that holds no twist, a configuration that exists only because PolyFrame's general method holds no twist at a fixed base (stability supports are FZ only). A real fixed base holds the twist, and PyNite with the twist held at the bases gives 15.1 / 4.50 / 6.92 against the twist-free 6.31 / 4.08 / 6.39 of the three frames; both tools also sit below the rigid-fork member reference on the rafters (4.52 / 6.87). The twist-free model is therefore conservative for a real fixed base in both tools, PolyFrame's more so. The 10 % one-sided bound stays as a watch. Out-of-plane buckling can be represented on this basis.
2. *Commercial recapture of the 22 PolyFrame cases.* Needs the company's commercial PolyFrame sign-in on this machine; 25 minutes once signed in. Until then, out-of-plane PolyFrame numbers stay out of project documents.
3. *FEM Design out of plane.* Whether FEM Design's bar element carries a warping DOF decides whether it can act as a reference; a manual check by the engineer, then four captures if it can.
4. *Insurer.* A written answer from the PI insurer on in-house modified open-source FE tools before any verification note is cited to a third-party check. Parallel, not blocking.
5. *Project volume.* The CFO's threshold (5 to 15 frame-buckling checks a year) decides whether the narrow B slice (pinned bases, restraints and loads at the shear centre, rafter restraints at L/2 and L/4) is worth its half day. Only the office knows the number.
6. *Merge.* The branch is ready for the user to merge; nothing has been committed from this session.

---


Date 2026-10-02. Eight-seat board (blind round, one rebuttal, neutral director) on: what to do after the PolyFrame buckling harness, whether to add PolyBeam/PolyColumn testing, and how far the results justify real-world use. Evidence base: `docs/POLYFRAME_BUCKLING_COMPARISON.md`, `Testing/polyframe/README.md`, `Testing/test_fem_design_buckling.py`.

## Ruling (director, unedited)

**Ruling:** Merge only after C-lite (explicit per-member warping length, no silent default, corrected docstring) is in the library and the two closed-form checks of the 2L rule pass; then settle the four open out-of-plane cases with an independent reference before any new PolyFrame cases; do not build a PolyBeam/PolyColumn harness (option A).

**Confidence:** medium. In-plane facts are strong and two-sourced; the out-of-plane ordering rests on a rule found after seeing the reference, and on four unresolved cases where PyNite is the unconservative tool.

**Why:**
1. Nay Sayer, conceded by CFO, Risk and Customer: the rule "warping length = twist-restraint spacing" is applied by the harness (`pynite_side.warping_length_factor`), not by PyNite, while FEModel3D.py claims conservatism and hard-wires the member length. A verification that depends on the verifier's patch does not verify the library. Fix is 0.5 agent-day.
2. Eye Doctor, conceded by Yes Man and Nay Sayer: all four xfails have PyNite higher, and the symmetric ±5 % policy hides it. More PolyFrame captures cannot adjudicate a two-tool disagreement; a third reference can. Hand M_cr (SN003 C1, C2, z_g) is decisive on the top-flange portal (gap is a factor 1.8); a shell model settles the fixed-base 7 %.
3. Purist and Nay Sayer, confirmed in repo: `test_buckling.py` already covers M_cr against Timoshenko and Gere and the three-factor formula. PolyBeam's M_cr is that same formula; comparing to it measures the formula twice through a 775-line OCR reader. A adds no information.
4. Yes Man and Long-Term Thinker, uncontested: in-plane α_cr agrees with two independent commercial tools (PolyFrame 15/15 within 0.2 %, FEM Design 58/58 modes within 5 %). That is project grade for planar steel frames.
5. Eye Doctor: twist is interpolated linearly, so a localised double-curvature mode over 8 elements can sit several percent high, the direction and size of the fixed-base gap. A 16/32-element rerun of the four cases takes minutes.

**Moved in rebuttal:** every seat moved toward C-lite before merge (Yes Man, Long-Term, CFO, Risk from "merge now"). Yes Man dropped A and withdrew "explained and bounded". Nay Sayer withdrew "curve fitting" to "physically defensible, not independently shown". Risk downgraded the student licence from a merge blocker to a staged condition. Customer withdrew FEM Design as the out-of-plane referee.

**Dissent on record (Customer, Yes Man):** PolyBeam/PolyColumn produce the numbers a Danish checker already accepts, and the four open cases (fixed-base halls with cranes, top-flange rafter loads) are the weekly project, roughly 5 to 8 of 15 to 25 steel frames a year. Rejecting A removes the one procedurally unquestioned reference, and the path to usable out-of-plane results is 1.5 to 2 agent-days.

**Conditions and kill criteria:**
- Before merge: C-lite in library; harness moved onto the new argument; a test that fails on the member-length default; closed-form N_cr,T test (twist held at one end, L_T = 2L); README states three gates sit at 10 % one-sided and one case at 100 % is a watch, not a gate; corpus README states the licence the PDFs were captured under.
- Before out-of-plane results enter a project document: recapture the 22 cases under a commercial licence (25 min).
- Kill: an independent 7-DOF or shell check of `beam-ipe300-two-span` near 23.5 rather than 16.1 to 17.8 reopens all out-of-plane gates. Hand M_cr on the top-flange portal siding with PolyFrame stops out-of-plane use until PyNite's load-height term is fixed. A 32-element rerun moving fixed-base portals by more than 3 % means the default mesh is inadequate. No real project needing α_cr,op within two quarters: B waits.

**Open questions:** shell model of one twist-free column and the three fixed-base portals (engineer, about half a day); hand M_cr for the top-flange portal (engineer, 1 h); whether FEM Design's bar carries a warping DOF; licence type of the PolyBeam/PolyColumn installs; PI insurer's position on in-house modified FE tools.

**Confidence for real-world use:**
- In plane: high. Planar steel frames, doubly symmetric rolled sections, 2 to 5 members, 1.0 G + 1.0 W. Usable now for the α_cr ≥ 10 sway check (EN 1993-1-1 5.2.1) and N_cr into the engineer's own 6.3 checks. Cite the FEM Design leg until the PolyFrame corpus is recaptured.
- Out of plane: low to medium, envelope-bound. Usable for pinned-base frames, columns and simply supported beams with loads at the shear centre, with the twist-restraint spacing set explicitly. Fixed-base portals: divide α_cr,op by 1.10 as an interim margin or confirm by hand. Destabilising loads inside frames: not usable without an independent M_cr on file. Not verified: timber, haunches, tapers, intermediate purlin restraints, mono-symmetric sections, second-order forces, imperfections.

**Scope:**
- Must (now): C-lite; closed-form N_cr,T 2L test; mesh rerun at 16/32 elements of the four open cases; README scope and one-sided gate statements; licence note in corpus; merge.
- Should (this quarter): hand M_cr on the top-flange portal; shell model of one twist-free column and the fixed-base portals; commercial recapture; one B slice inside PolyFrame's valid domain (pinned bases, restraints and loads at shear centre, rafter restraints at L/2 and L/4); one-page dated verification note with commit hash, versions, mode table, envelope.
- Won't (now): PolyBeam/PolyColumn harness (A); warping DOF; utilisation or χ output; tapered, haunched or timber cases; upstream PR until rebased and C-lite is part of the buckling feature.
- Main design decision: explicit warping length (twist-restraint spacing) as a required per-member input to the equivalent-torsion approximation, enforced by test. Rejected: a 7-DOF warping element (5 to 10 agent-days, no customer waiting).
- Measure after: out-of-plane cases within 5 % one-sided per reference; open cases with PyNite on the unconservative side (target zero or explained); real projects per quarter where a PyNite α_cr,op replaced or pre-checked a PolyFrame run.

## Board minutes

| Seat | Round 1 | Round 2 | Core argument |
|---|---|---|---|
| Nay Sayer | Conditions: C first, oppose A | Conditions (same) | Harness-side rule is not library verification; all four misses are unconservative |
| Yes Man | Conditions: D, C, B, A | Conditions: C-lite, recapture, merge; A dropped | In-plane two-sourced and project grade; "explained" was too strong |
| Eye Doctor | Conditions: C-lite, narrow B | Conditions: C-lite, mesh rerun, third reference | Sign of the gaps; shell model is the decider; hand M_cr settles top-flange |
| Purist | Conditions: D plus C-lite, cut A | Conditions: same plus third reference before fixed-base use | Closed forms beat OCR of PolyBeam; add the 2L torsional test |
| Long-Term Thinker | Conditions: D checkpoint, then C, B | Conditions: C-lite before merge | A false conservatism claim must not reach main or an upstream PR |
| Customer | Conditions: C-lite, B, A later | Conditions (same); signs in-plane on FEM Design leg only | Fixed-base and top-flange cases are the weekly project |
| CFO | Conditions: merge, B; oppose A | Conditions: C-lite, third reference, then B | 0.5 + 0.5-1 + 0.5 agent-days to usable out-of-plane |
| Risk Officer | Conditions: C-lite, then merge; licence recapture | Conditions; licence staged, not blocking | Student licence; one-sided gate; insurer parallel |
