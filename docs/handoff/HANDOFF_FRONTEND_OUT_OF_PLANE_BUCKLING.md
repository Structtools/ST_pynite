# Handoff to the frontend: out-of-plane buckling of 2D frames

Date 2026-10-02. From the ST_pynite buckling work (branch `out-of-plane-buckling-2D-frames`).
Companion documents: `HANDOFF_BACKEND_CALCULATION_DOCUMENTATION.md` (how the numbers are
computed) and `HANDOFF_BACKEND_CHANGES.md` (what the backend must change). Verification:
`docs/POLYFRAME_BUCKLING_COMPARISON.md`.

## What the engine will deliver per frame and load combination

| Result | Meaning | Status |
|---|---|---|
| `alpha_cr` (in plane) | smallest load factor at which the frame buckles in its plane (EN 1993-1-1 5.2.1). Below 10 the frame is sway-sensitive and second-order effects apply. | verified: 0.2 % against PolyFrame on 15 frames, 5 % against FEM Design on 58 modes |
| `alpha_cr,op` (out of plane) | smallest load factor for lateral and lateral-torsional buckling of the frame as a whole (EN 1993-1-1 6.3.4) | verified: 18 of 22 frames within 5 % of PolyFrame, the rest on the safe side; see "Limits" |
| governing member | the member whose displacement dominates the first out-of-plane mode | PyNite's own; PolyFrame's label was found unreliable |
| mode shape | normalised displacements, for drawing | available for both planes |

Nothing else changes in the static results. The engine does not deliver utilisations or
buckling lengths for out-of-plane buckling; the backend derives what the Eurocode checks need
(see the backend handoffs).

## What the frontend must collect

Out-of-plane buckling is decided by restraints that a 2D drawing does not show. Each item
below changes `alpha_cr,op` directly; defaults are stated so the engine can run without
them, but the user must be able to see and change every one.

1. **Lateral restraint at nodes.** Default: every node of the frame holds the lateral
   (out-of-plane) displacement. That is PolyFrame's convention and the verified envelope.
   A node that is *not* laterally held (a free cantilever tip, an unbraced ridge) must be
   markable as such. Note PolyFrame cannot represent that; our engine can.
2. **Intermediate lateral restraints along a member** (purlins, side rails, sheeting). Input:
   positions along the member and which flange they hold (top, bottom, both). Today the
   payload carries `segments` with a `beta` per axis (`stiff`, `weak`, `top`, `bottom`);
   that stays for the member checks, but the engine needs the **positions**, not a beta.
   Suggested: a list of restraint positions per member with `flange: top | bottom | both`
   and `holdsTwist: bool`. Not yet in the verified envelope; see "Limits".
3. **Twist restraint at supports.** Default: a support holds the lateral displacement but
   **not** the twist of the member (PolyFrame's convention, conservative). A real fixed
   base, a stiff base plate, or a fork support holds the twist: offer a checkbox
   "holds twist" on supports (default off). Holding it raised `alpha_cr,op` by 10 to 140 %
   on the verified fixed-base frames, so it is worth asking.
4. **Load position on the section.** Per line load and point load: `centroid` (shear
   centre), `top flange`, `bottom flange`. Today `globalElements.loadPosition` carries
   `"centroid"`; it must become per load and offer the two flanges. A gravity load on the
   top flange lowered `alpha_cr,op` by 22 % on a simply supported beam and 34 % on a
   portal rafter; on the bottom flange it raised it 30 %.
5. **Section warping constant `Iw`** must be present for every steel section used
   (it is in the section library already). Timber sections: not in the verified envelope.

## What to show

- `alpha_cr` and `alpha_cr,op` per ULS combination, the governing combination, the
  governing member, and the two mode shapes (front view for in plane, top and side view
  for out of plane, as PolyFrame draws them).
- The sway flag: `alpha_cr < 10` means second-order effects must be included
  (EN 1993-1-1 5.2.1(3)); display it as a state, not only as a number.
- A **verified-envelope notice** on the out-of-plane result (text supplied by the backend):
  planar steel frames of rolled I-sections, lateral restraint at every node, loads at the
  stated position, 1.0 G + 1.0 W type combinations. Outside it the result is shown with a
  warning, not hidden.
- Where the engine warns that a member's twist is free at an end (the engine does this
  itself), surface the warning next to that member: it means the user should confirm the
  support's twist restraint.

## Limits the user must be told (as of 2026-10-02)

- Fixed bases are modelled as holding no twist unless the user says otherwise. That is
  on the safe side: both PolyFrame and our engine give 6 to 140 % lower factors than a
  base that holds twist. Showing the "holds twist" option is the way to recover that
  margin when the detail justifies it.
- Our engine and PolyFrame agree on top-flange loads for simply supported beams (both
  -22 %). On a portal rafter PolyFrame drops 64 % where our engine and an independent
  reference drop 34 %; a user comparing with PolyFrame will see that difference. Ours is
  the one the reference supports.
- Not yet verified: intermediate lateral restraints along a member, haunches, tapered
  members, timber, mono-symmetric sections, load combinations with partial factors other
  than 1.0, second-order forces. Results for those must carry the warning above.
- PolyFrame refuses straight-only structures and laterally free ends in its general
  method. Our engine handles both; do not expect PolyFrame to confirm them.
