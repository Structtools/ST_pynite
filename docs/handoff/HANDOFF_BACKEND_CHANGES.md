# Handoff to the backend: changes required to carry out-of-plane buckling into the Eurocode checks

Date 2026-10-02. Companion of `HANDOFF_BACKEND_CALCULATION_DOCUMENTATION.md`. Written against
what the structural_backend harness recorded of the current backend (`docs/POLYFRAME_*.md`,
`tools/polyframe/translate.py` there): steel and timber member stability checks are skipped
for `entityName: "frame"`, buckling lengths come from `segments` as `beta * L`, the FE solve
is Bernoulli-Euler, supports are the 2D presets `fixed / pinned / roller-x / roller-z`.
Changes are required; they are listed in the order they unblock each other.

## A. Engine integration

1. **Run two eigenvalue analyses per ULS combination** on the frame model:
   `analyze_buckling(combo, plane='XY')` for `alpha_cr` and
   `analyze_buckling(combo, plane='XY', out_of_plane=True, warping='equivalent_torsion')`
   for `alpha_cr,op`. Keep `elements_per_member=8`. Store `alpha_cr`, `alpha_cr,op`, the
   governing member (largest out-of-plane displacement in mode 1) and the mode shapes.
2. **Sections** must carry `Iw` and the EN 1993-1-1 6.2.6 shear areas (`Avz` web, `Avy`
   flanges) into `add_section(..., Asy=Avy, Asz=Avz, Iw=Iw)`. The library has `Iw`; the
   shear areas are derived from `h, b, tw, tf, r` as in `Testing/polyframe/sections.py`.
3. **Member orientation.** Members are rotated 90 degrees so the Eurocode strong axis bends
   in the frame plane (the convention in `Testing/polyframe/pynite_side.py` and the PyNite
   frame tests).
4. **Out-of-plane supports.** Every frame node: hold `DZ`. Supports: hold `DZ`; hold `RX`
   and `RY` only when the user's support says it holds twist (new input, default off).
   This is the verified convention; see B.3.
5. **Warping length per member.** Call `def_warping_length(member, L_w)` with `L_w` from the
   twist restraints: member length when both ends are held (frame corner, fork, twist-holding
   support), `2 L` when one end is held, `math.inf` when none. The rule and a reference
   implementation are `warping_length_factor` in `Testing/polyframe/pynite_side.py`. PyNite
   warns when a member is left at the default with a twist-free end; treat that warning as a
   validation error of the input, not as noise.
6. **Load height.** Pass `load_height = +h/2` for a gravity load on the top flange and
   `-h/2` on the bottom flange to `add_member_dist_load` / `add_member_pt_load` (sign of
   Annex F `z_g`). Loads at the shear centre: 0. This needs the per-load position input (B.2).
7. **Shear deformation.** The backend FE is Bernoulli-Euler; PyNite's buckling analysis is
   Timoshenko (as PolyFrame and FEM Design are). Either accept the two solves differing by
   the shear effect (0.2 to 15 % on forces for stocky IPE300 frames, 76 % on a GL28h base
   moment, per the earlier harness) or move the static solve to PyNite as well. EN 1993-1-1
   does not require shear deformation for steel; EN 1995-1-1 2.2.3 arguably does for timber.

## B. Input model

1. **Lateral restraints along a member**: positions and flange (`top | bottom | both`,
   `holdsTwist`). These become extra nodes with `DZ` held (and `RX`/`RY` for a twist
   restraint) in the engine model, and they split the member's warping length. Replaces the
   out-of-plane role of `segments.beta` for `weak`, `top`, `bottom` once the eigenvalue path
   is used; keep `beta` for the in-plane strong-axis check only if 6.3.1 is still applied
   that way (see C.2).
2. **Load position per load** (`centroid | top | bottom`), replacing the single
   `globalElements.loadPosition: "centroid"`.
3. **Support twist flag** (`holdsTwist`, default false) on the support presets. A real
   fixed base holds twist; the default stays off because it is the conservative and
   verified convention.
4. **Laterally free nodes** (`laterallyFree: true`, default false) for cantilever tips and
   unbraced ridges; the default follows PolyFrame.

## C. Eurocode checks (EN 1993-1-1), steel

The frame member stability checks are currently skipped. With `alpha_cr` and
`alpha_cr,op` available, two routes are open; the first is the one PolyFrame's general
method uses and was verified against.

1. **General method, 6.3.4.** Per combination: `alpha_ult,k` = the minimum load amplifier
   for the in-plane cross-section resistance of the critical section (first-order forces,
   with second-order amplification when `alpha_cr < 10`, see C.3);
   `lambda_op = sqrt(alpha_ult,k / alpha_cr,op)`; `chi_op = min(chi, chi_LT)` from the
   buckling curves 6.3.1.2 and 6.3.2.2 at `lambda_op` (DK NA values of `alpha_LT` and
   `lambda_LT,0`); check `chi_op * alpha_ult,k / gamma_M1 >= 1`. Report the four
   intermediates as PolyFrame does. 6.3.4 is permitted for "structural components ... with
   mono-symmetric cross-sections"; for doubly symmetric sections it is the usual
   application. Its DK NA status must be confirmed by the engineer owning the checks.
2. **Member method, 6.3.1 to 6.3.3, with eigenvalue-derived lengths.** For each compression
   member `N_cr = alpha_cr N_Ed` gives `L_cr,y = pi sqrt(E I_y / N_cr)` for the in-plane
   check (this replaces `beta * L`, which the earlier harness showed is 2 to 4 times too
   short on sway frames when `beta = 1`). Out of plane, `M_cr = alpha_cr,op M_Ed` for the
   governing member gives `lambda_LT` directly (6.3.2.2) without the three-factor formula;
   `N_cr,z` and `N_cr,T` for 6.3.1 follow from the member's own restraint spacing. 6.3.3
   interaction with Annex A or B factors as now. Note the earlier harness's open findings on
   the existing 6.3.3 implementation (division by `gamma_M1`, `C_my` on sway columns) must be
   closed before this route is trusted.
3. **Second-order effects in plane, 5.2.1 and 5.2.2.** When `alpha_cr < 10` (3 for plastic
   analysis), amplify the sway effects by `1 / (1 - 1 / alpha_cr)` (5.2.2(5)B, valid for
   `alpha_cr >= 3`) or run a second-order analysis, and include the sway imperfection of
   5.3.2 (or the equivalent horizontal forces). PolyFrame uses the first critical mode as
   the imperfection shape (5.3.2(11)); the mode shape is available from the engine.
4. **Partial factors.** DK NA `gamma_M0 = 1.10`, `gamma_M1 = 1.20` (as read back from
   PolyFrame's DK reports), applied in the backend, never in the engine.

## D. Timber, EN 1995-1-1

Not verified; the engine runs the same eigenvalue analysis for a rectangular glulam
section if `Iw` is set (for a solid rectangle `Iw` is small and may be taken as zero, with
`def_warping_length(..., math.inf)` to say so). 6.3.3 lateral-torsional check with
`sigma_m,crit = M_cr / W_y` from `alpha_cr,op`, and `k_crit`; 6.3.2 with `lambda_rel`
from `N_cr`. Keep the timber frame stability check skipped until a timber comparison exists.

## E. Validation to add on the backend side

- Re-run the 22 catalogue frames through the backend path and assert the same `alpha_cr`
  and `alpha_cr,op` as `Testing/polyframe` (they are in `Testing/polyframe/corpus/` and the
  comparison document), so the integration cannot drift from the verified engine.
- A test that a frame with a twist-free end left at the default warping length is rejected
  with PyNite's warning.
- A test of the 6.3.4 intermediates against PolyFrame's printed `alpha_ult,k`, `lambda_op`,
  `chi_op` and utilisation for the pinned-base portals (in the corpus records and PDFs).

## F. What is not required

No PolyBeam or PolyColumn harness (the board ruled it adds nothing a closed form does not).
No warping DOF in the engine. No utilisation output from the engine. No change to the
in-plane statics beyond the shear-deformation decision in A.7.
