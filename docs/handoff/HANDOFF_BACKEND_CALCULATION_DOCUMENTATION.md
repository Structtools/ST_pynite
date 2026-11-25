# Handoff to the backend: how out-of-plane buckling is computed (for the calculation documentation)

Date 2026-10-02. Source: ST_pynite branch `out-of-plane-buckling-2D-frames`,
`Pynite/Analysis.py`, `Pynite/FEModel3D.py`, `Pynite/Member3D.py`, `docs/source/buckling.rst`,
verification in `docs/POLYFRAME_BUCKLING_COMPARISON.md` and `Testing/polyframe/README.md`.
This is the text the calculation documentation can be written from; the numbers are the
verified ones.

## 1. Method

Linear (eigenvalue) buckling of the whole frame: `[K - lambda Kg] phi = 0`, where `K` is the
elastic stiffness of the frame and `Kg` the geometric stiffness assembled from the first-order
member forces under the load combination. The smallest positive `lambda` is the critical load
factor `alpha_cr`: the loads of the combination multiplied by `alpha_cr` reach the elastic
critical load (EN 1993-1-1 5.2.1(3) for in plane; 6.3.4(3) for `alpha_cr,op`).

- **Elements.** 12-DOF Timoshenko beam elements (shear deformation included; EN 1993-1-1
  does not require it for steel, PolyFrame and FEM Design include it). Each member is
  subdivided into 8 elements for the eigenproblem (refining to 32 changes the verified
  factors by less than 2 %).
- **Geometric stiffness.** The axial-force terms (flexural and torsional buckling) and
  the bending-moment, shear and torque terms (lateral-torsional buckling and
  beam-column interaction), derived from the second-order work of the initial stresses
  for a doubly symmetric section with a linear moment variation over each element
  (`Derivations/Geometric Stiffness Matrix - Moments and Torque.py`). A transverse member
  load applied at height `z_g` from the shear centre adds `-|p| z_g` to the twist stiffness
  (EN 1993-1-1 Annex F / NCCI SN003 sign of `z_g`: positive towards the shear centre,
  destabilising).
- **Two analyses per frame.** In plane: the out-of-plane DOFs are held at every node and
  the eigenproblem gives `alpha_cr`. Out of plane: the frame's own out-of-plane supports
  are used for the static pre-solve and the in-plane DOFs are held in the eigenproblem,
  which gives `alpha_cr,op` for lateral and lateral-torsional buckling (6.3.4). For a planar
  frame loaded in its plane the two sets of modes are uncoupled and their union is the 3D
  result.
- **Warping.** The element has no warping degree of freedom. The warping stiffness is
  approximated by the equivalent torsion constant `J_eff = J + pi^2 E Iw / (G L_w^2)`, with
  `L_w` the member's **warping length**: the distance between the points where its twist is
  held. `L_w` is the member length between two fork supports or frame corners, `2 L` for a
  member whose twist is held at one end only (a column on a base that does not hold its
  twist; each span of a continuous beam at an intermediate support), and no warping
  stiffness for a member held at neither end. This is exact for a member in uniform bending
  between fork supports and an approximation otherwise; against an independent Vlasov
  solution it is within 2 % for columns with zero base moment and 3 to 6 % conservative for
  members in non-uniform bending.

## 2. Boundary conditions (the conventions the numbers are verified for)

- Every node of the frame holds the lateral (out-of-plane) displacement. This is
  PolyFrame's general-method convention (its report prints "Stability Supports: FZ" at
  every node) and the basis of the comparison.
- A support holds the lateral displacement and **not** the twist of the member, unless the
  user says the detail does. This is conservative: for the three verified fixed-base frames,
  holding the twist raised `alpha_cr,op` from 6.31 / 4.08 / 6.39 to 15.1 / 4.50 / 6.92.
- Frame corners hold the twist of both members meeting there (the other member's bending
  stiffness does it; stiffness ratios of 50 to 300 justify treating it as a fork).
- Warping is free everywhere.

## 3. Verification (PolyFrame 1.1.11 general method, FEM Design 2024, independent Vlasov solver)

| Quantity | Reference | Cases | Result |
|---|---|---|---|
| `alpha_cr` in plane | PolyFrame | 15 frames, factors 5.7 to 225 | all within 0.2 % |
| in-plane modes | FEM Design | 16 scenarios, 58 modes | all within 5 %, worst 4.5 % |
| `alpha_cr,op` | PolyFrame | 22 frames (columns, beams, portals, gable, two-bay) | 18 within 5 %; 3 fixed-base portals PolyFrame 6.5 to 7.6 % lower; 1 top-flange portal PolyFrame 46 % lower |
| `alpha_cr,op` | Rayleigh-Ritz Vlasov member solver (`Testing/polyframe/ltb_reference.py`) | governing member of 15 cases | pinned-base columns and gable 0.8 to 2.7 %; two-span beam 1.5 %; simply supported beams 0.5 to 1.2 %; top-flange portal rafter 0.5 % |

Readings of the two open clusters:

- **Fixed-base portals.** Both tools lie below the member reference and far below a base
  that holds the twist; both are on the safe side, PolyFrame more so. Accepted as such on
  2026-10-02.
- **Top-flange load on a portal rafter.** The independent reference gives 2.655, PyNite
  2.64, PolyFrame 1.44. PolyFrame's load-height response in frames is asymmetric (a
  bottom-flange load raised it only 11 %) and is the outlier.

## 4. Verified envelope (to be stated in the documentation)

Planar steel frames of doubly symmetric rolled I-sections (IPE200, IPE300, HE200A tested),
S355, 2 to 5 members, lateral restraint at every node, twist held at frame corners and not
at supports, loads at the shear centre or on a flange, one load combination of the form
1.0 G + 1.0 W, first-order forces. Not verified: intermediate restraints along a member,
haunches, tapered members, timber, mono-symmetric sections, second-order forces,
imperfections, SLS, fire.

## 5. Interpretation for the engineer of record

`alpha_cr` and `alpha_cr,op` are elastic critical factors. The design check is the
engineer's (EN 1993-1-1 6.3.4 general method, or 6.3.1 to 6.3.3 with `N_cr = alpha_cr N_Ed`
and `M_cr = alpha_cr,op M_Ed`), with imperfections and second-order effects where
`alpha_cr < 10` (5.2.1(3), 5.3.2). The engine does not apply partial factors and does not
produce a utilisation.

## 6. References

EN 1993-1-1:2005 5.2.1, 5.2.2, 5.3.2, 6.3.1 to 6.3.4; ENV 1993-1-1 Annex F and NCCI SN003
(three-factor `M_cr`, `C1`, `C2`, `z_g`); Timoshenko & Gere, *Theory of Elastic Stability*,
ch. 2 and 6; Trahair, *Flexural-Torsional Buckling of Structures*; PolyStruc PolyFrame 1.1.11
reports (`Testing/polyframe/corpus/reports/`); FEM Design 2024 critical parameters
(`Testing/*_comparison.py`).
