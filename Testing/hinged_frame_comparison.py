r"""
hinged_frame_comparison.py
==========================
Single-bay portal frame with a hinge (charnier) at the top-left corner
(connection between left column and beam at node B).

Usage
-----
    uv run python Testing/hinged_frame_comparison.py

Structure topology
------------------
Single-bay portal frame in the XY-plane (Z=0).
Columns are vertical (along Y). Beam is horizontal (along X).
Column bases are fixed. Hinge at node B (beam i-end).
Frame is restrained from out-of-plane moment.

    B --o--------- C
    |               |
    | Col1          | Col2
    |               |
    A               D
  (base)          (base)

    o = hinge (moment release at beam i-end)

    H      = 3.0 m  (column height)
    B_SPAN = 4.0 m  (beam span)
    Section: IPE 300 (all members)
    Steel: S235

Scenarios
---------
    1  Vertical UDL on beam: 10 kN/m downward
    2  Horizontal point loads at B and C, same direction (-X): 10 kN each
    3  Distributed loads on columns pointing away from frame:
       Col1: -10 kN/m (leftward), Col2: +10 kN/m (rightward)
    4  Combined: vertical UDL 10 kN/m on beam + horizontal point loads
       at B (+10 kN, rightward/inward) and C (-10 kN, leftward/inward)
    5  Same-direction column loads: both -10 kN/m leftward (-X)
    6  Vertical point loads at B and C: -10 kN downward each
    7  Asymmetric sway: single -10 kN point load at C (leftward, -X)
"""

import sys
import os
import numpy as np

# Make sure the project root is on the path when running from the Testing/ dir
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from Pynite import FEModel3D

# ============================================================================
# SCENARIO SELECTOR
# ============================================================================
#
#  0  Manual         use the PARAMETERS block below as-is
#  1  Gravity UDL    10 kN/m downward on beam
#  2  Horizontal P   10 kN point loads at B & C in -X direction
#  3  Column wind    10 kN/m on columns pointing away from frame
#  4  Combined       10 kN/m UDL on beam + 10 kN inward point loads at B & C
#  5  Same-dir col  10 kN/m on both columns in -X direction
#  6  Vertical P    10 kN downward point loads at B & C
#  7  Asymmetric    single 10 kN point load at C in -X direction
#
SCENARIO = int(os.getenv("SCENARIO", "1"))

# ============================================================================
# PARAMETERS -- used when SCENARIO = 0, or as defaults overridden by presets
# ============================================================================

# --- Geometry (metres) ------------------------------------------------------
H = 3.0          # column height
B_SPAN = 4.0     # beam span (horizontal distance between columns)

# Number of FEM sub-elements per member (more = better accuracy, slower)
N_ELEM = 8

# --- Material: S235 steel (SI units: N, m, Pa) -----------------------------
E = 210e9        # Young's modulus
G = 80.769e9     # Shear modulus
NU = 0.3         # Poisson's ratio
RHO = 0.0        # density (set to 7850 for self-weight)

# --- Cross-section: IPE 300 ------------------------------------------------
SECTIONS = {
    "IPE300": dict(A=53.8e-4, Iy=8356e-8, Iz=604e-8, J=20.1e-8,
                   Asy=32.10e-4, Asz=25.67e-4),
}

COLUMN_SEC = "IPE300"
BEAM_SEC = "IPE300"

# --- Boundary conditions ----------------------------------------------------
FIXED_BASE = (True, True, True, True, True, True)

# --- Distributed load on beam (N/m, global FY) -----------------------------
BEAM_LOAD = -10e3   # negative = downward

# --- Point loads at column-top nodes (N, global FX) -------------------------
LOAD_B_FX = 0.0
LOAD_C_FX = 0.0

# --- Point loads at column-top nodes (N, global FY) -------------------------
LOAD_B_FY = 0.0
LOAD_C_FY = 0.0

# --- Distributed loads on columns (N/m, global FX) -------------------------
COL_LOAD_1_FX = 0.0   # Col1 (A->B)
COL_LOAD_2_FX = 0.0   # Col2 (D->C)

# --- Lateral bracing / 2D constraint ----------------------------------------
LATERAL_BRACE_NODES = True

# --- Beam element formulation -----------------------------------------------
USE_TIMOSHENKO = True

# --- Buckling analysis settings --------------------------------------------
NUM_MODES = 5
COMBO_NAME = "Combo 1"

# ============================================================================
# END OF PARAMETERS
# ============================================================================


# ---------------------------------------------------------------------------
# Scenario presets
# ---------------------------------------------------------------------------

_PRESETS = {
    1: dict(
        name="Vertical UDL on beam: 10 kN/m downward",
        BEAM_LOAD=-10e3,
        LOAD_B_FX=0.0,
        LOAD_C_FX=0.0,
        LOAD_B_FY=0.0,
        LOAD_C_FY=0.0,
        COL_LOAD_1_FX=0.0,
        COL_LOAD_2_FX=0.0,
    ),
    2: dict(
        name="Horizontal point loads at B & C, same direction (-X): 10 kN each",
        BEAM_LOAD=0.0,
        LOAD_B_FX=-10e3,
        LOAD_C_FX=-10e3,
        LOAD_B_FY=0.0,
        LOAD_C_FY=0.0,
        COL_LOAD_1_FX=0.0,
        COL_LOAD_2_FX=0.0,
    ),
    3: dict(
        name="Column loads pointing away from frame: Col1 -10 kN/m, Col2 +10 kN/m",
        BEAM_LOAD=0.0,
        LOAD_B_FX=0.0,
        LOAD_C_FX=0.0,
        LOAD_B_FY=0.0,
        LOAD_C_FY=0.0,
        COL_LOAD_1_FX=-10e3,
        COL_LOAD_2_FX=10e3,
    ),
    4: dict(
        name="Combined: 10 kN/m UDL on beam + 10 kN inward point loads at B & C",
        BEAM_LOAD=-10e3,
        LOAD_B_FX=10e3,
        LOAD_C_FX=-10e3,
        LOAD_B_FY=0.0,
        LOAD_C_FY=0.0,
        COL_LOAD_1_FX=0.0,
        COL_LOAD_2_FX=0.0,
    ),
    5: dict(
        name="Same-direction column loads: both +10 kN/m rightward (+X)",
        BEAM_LOAD=0.0,
        LOAD_B_FX=0.0,
        LOAD_C_FX=0.0,
        LOAD_B_FY=0.0,
        LOAD_C_FY=0.0,
        COL_LOAD_1_FX=10e3,
        COL_LOAD_2_FX=10e3,
    ),
    6: dict(
        name="Vertical point loads at B and C: -10 kN downward each",
        BEAM_LOAD=0.0,
        LOAD_B_FX=0.0,
        LOAD_C_FX=0.0,
        LOAD_B_FY=-10e3,
        LOAD_C_FY=-10e3,
        COL_LOAD_1_FX=0.0,
        COL_LOAD_2_FX=0.0,
    ),
    7: dict(
        name="Asymmetric sway: single -10 kN point load at C (leftward, -X)",
        BEAM_LOAD=0.0,
        LOAD_B_FX=0.0,
        LOAD_C_FX=-10e3,
        LOAD_B_FY=0.0,
        LOAD_C_FY=0.0,
        COL_LOAD_1_FX=0.0,
        COL_LOAD_2_FX=0.0,
    ),
}


def _apply_scenario():
    """Overwrite module globals with the selected scenario preset."""
    if SCENARIO == 0:
        return
    if SCENARIO not in _PRESETS:
        raise ValueError(
            f"SCENARIO={SCENARIO} is not defined.  "
            f"Valid choices: 0 (manual) or {sorted(_PRESETS)}"
        )
    g = globals()
    for key, val in _PRESETS[SCENARIO].items():
        if key == "name":
            continue
        g[key] = val


# ---------------------------------------------------------------------------
# Model builder
# ---------------------------------------------------------------------------


def _add_column(model, name, base_node, top_node, mat, sec, n_elem, col_height, base_x):
    """Add a vertical column with n_elem intermediate nodes."""
    h = col_height / n_elem
    for i in range(1, n_elem):
        nn = f"_{name}_int{i}"
        model.add_node(nn, base_x, i * h, 0.0)
    model.def_support(base_node, *FIXED_BASE)
    model.add_member(name, base_node, top_node, mat, sec, rotation=90)


def _add_beam(model, name, node_i, node_j, mat, sec, n_elem, xi, yi, xj, yj):
    """Add a beam with n_elem intermediate nodes."""
    for i in range(1, n_elem):
        t = i / n_elem
        nn = f"_{name}_int{i}"
        model.add_node(nn, xi + t * (xj - xi), yi + t * (yj - yi), 0.0)
    model.add_member(name, node_i, node_j, mat, sec, rotation=90)


def build_model():
    model = FEModel3D()

    model.add_material("Steel", E, G, NU, RHO)
    for sec_name, props in SECTIONS.items():
        if USE_TIMOSHENKO:
            model.add_section(
                sec_name, props["A"], props["Iy"], props["Iz"], props["J"],
                Asy=props.get("Asy"), Asz=props.get("Asz"),
            )
        else:
            model.add_section(
                sec_name, props["A"], props["Iy"], props["Iz"], props["J"],
            )

    # --- Nodes ---------------------------------------------------------------
    # Column bases
    model.add_node("A", 0.0, 0.0, 0.0)         # left column base
    model.add_node("D", B_SPAN, 0.0, 0.0)       # right column base

    # Column tops / beam ends
    model.add_node("B", 0.0, H, 0.0)            # top-left  (hinge here)
    model.add_node("C", B_SPAN, H, 0.0)          # top-right

    # --- Members -------------------------------------------------------------
    _add_column(model, "Col1", "A", "B", "Steel", COLUMN_SEC, N_ELEM, H, 0.0)
    _add_column(model, "Col2", "D", "C", "Steel", COLUMN_SEC, N_ELEM, H, B_SPAN)

    _add_beam(model, "Beam1", "B", "C", "Steel", BEAM_SEC, N_ELEM,
              0.0, H, B_SPAN, H)

    # --- Hinge at top-left corner (beam i-end at node B) ---------------------
    # Release in-plane moment (Ryi) and out-of-plane moment (Rzi) at i-end.
    model.def_releases("Beam1", Ryi=True, Rzi=True)

    # --- Distributed load on beam (global FY) --------------------------------
    if BEAM_LOAD != 0.0:
        model.add_member_dist_load("Beam1", "FY", BEAM_LOAD, BEAM_LOAD)

    # --- Point loads at column-top nodes (global FX) -------------------------
    if LOAD_B_FX != 0.0:
        model.add_node_load("B", "FX", LOAD_B_FX)
    if LOAD_C_FX != 0.0:
        model.add_node_load("C", "FX", LOAD_C_FX)

    # --- Point loads at column-top nodes (global FY) -------------------------
    if LOAD_B_FY != 0.0:
        model.add_node_load("B", "FY", LOAD_B_FY)
    if LOAD_C_FY != 0.0:
        model.add_node_load("C", "FY", LOAD_C_FY)

    # --- Distributed loads on columns (global FX) ----------------------------
    if COL_LOAD_1_FX != 0.0:
        model.add_member_dist_load("Col1", "FX", COL_LOAD_1_FX, COL_LOAD_1_FX)
    if COL_LOAD_2_FX != 0.0:
        model.add_member_dist_load("Col2", "FX", COL_LOAD_2_FX, COL_LOAD_2_FX)

    # --- Lateral bracing / 2D constraint --------------------------------------
    if LATERAL_BRACE_NODES:
        for node_name in model.nodes:
            if node_name not in ("A", "D"):  # bases already fully supported
                model.def_support(node_name, False, False, True, True, True, False)

    return model


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

SEP = "-" * 64


def _print_section_info():
    s = SECTIONS[COLUMN_SEC]
    print(f"  Section (all members): {COLUMN_SEC}")
    print(f"    A  = {s['A'] * 1e4:.2f} cm^2")
    print(f"    Iy = {s['Iy'] * 1e8:.0f} cm^4  (strong axis)")
    print(f"    Iz = {s['Iz'] * 1e8:.0f} cm^4  (weak axis)")
    print(f"    J  = {s['J'] * 1e8:.2f} cm^4")
    if USE_TIMOSHENKO and s.get('Asy') and s.get('Asz'):
        print(f"    Asy= {s['Asy'] * 1e4:.2f} cm^2  (shear area, weak-axis bending)")
        print(f"    Asz= {s['Asz'] * 1e4:.2f} cm^2  (shear area, strong-axis bending)")
    beam_type = "Timoshenko" if USE_TIMOSHENKO else "Euler-Bernoulli"
    print(f"  Beam formulation : {beam_type}")
    print()


def _print_loads():
    printed = False
    if BEAM_LOAD != 0.0:
        print(f"  Beam1 (B->C) : w = {BEAM_LOAD / 1e3:+.2f} kN/m  (distributed FY)")
        printed = True
    for node, F in [("B", LOAD_B_FX), ("C", LOAD_C_FX)]:
        if F != 0.0:
            print(f"  Node {node}       : FX = {F / 1e3:+.1f} kN  (point load)")
            printed = True
    for node, F in [("B", LOAD_B_FY), ("C", LOAD_C_FY)]:
        if F != 0.0:
            print(f"  Node {node}       : FY = {F / 1e3:+.1f} kN  (point load)")
            printed = True
    for col, w in [("Col1", COL_LOAD_1_FX), ("Col2", COL_LOAD_2_FX)]:
        if w != 0.0:
            print(f"  {col:14s}: w = {w / 1e3:+.2f} kN/m  (distributed FX)")
            printed = True
    if not printed:
        print("  (no loads applied)")
    brace = (
        "Yes -- 2D constraint (DZ, RX, RY restrained at all non-base nodes)"
        if LATERAL_BRACE_NODES
        else "No -- full 3D analysis"
    )
    print(f"  Lateral bracing : {brace}")
    print()


# ---------------------------------------------------------------------------
# FEM Design reference values (critical parameters)
# ---------------------------------------------------------------------------

FEM_DESIGN_REF = {
    1: [358.320, 1721.181, 1874.093, 3408.881],
    2: [7751.024, 18331.246, 28128.972],
    3: [59770.732, 140032.503, 212360.038],
    4: [337.393, 1158.908, 1859.044, 2074.372],
    5: [3728.163, 9953.064, 15721.756, 17457.187],
    6: [719.489, 3342.843, 3957.773, 7295.908],
    7: [4025.712, 11628.951, 16075.308, 21983.368],
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    _apply_scenario()

    print("=" * 64)
    print("  PyNite -- Hinged Portal Frame Buckling Analysis")
    print("=" * 64)
    print()
    if SCENARIO != 0:
        print(f"  Scenario {SCENARIO}: {_PRESETS[SCENARIO]['name']}")
        print()
    print(f"  Frame geometry : H = {H:.2f} m, B_SPAN = {B_SPAN:.2f} m")
    print(f"  Base condition : Fixed")
    print(f"  Hinge          : top-left corner (beam i-end at node B)")
    print(f"  Sub-elements   : {N_ELEM} per member")
    print()
    _print_section_info()
    print("  Applied loads")
    _print_loads()

    # --- Build and run -------------------------------------------------------
    model = build_model()

    print("  Running buckling analysis ...")
    results = model.analyze_buckling(combo_name=COMBO_NAME, num_modes=NUM_MODES)
    print()

    lams = results.load_multipliers
    n_modes = len(lams)

    # --- Load multipliers table ----------------------------------------------
    print(SEP)
    print(f"  {'Mode':>4}   {'Lambda_cr':>10}   Notes")
    print(SEP)
    for i, lam in enumerate(lams):
        note = "  <-- critical (lowest)" if i == 0 else ""
        print(f"  {i + 1:4d}   {lam:10.4f}{note}")
    print(SEP)
    print()

    # --- Effective lengths per member and mode --------------------------------
    all_members = ["Col1", "Col2", "Beam1"]
    ref_lengths = {"Col1": H, "Col2": H, "Beam1": B_SPAN}

    sec = SECTIONS[COLUMN_SEC]
    EIy = E * sec["Iy"]
    EIz = E * sec["Iz"]

    print(SEP)
    print(f"  Effective buckling lengths  (H = {H:.2f} m, B_SPAN = {B_SPAN:.2f} m)")
    print(SEP)
    print(
        f"  {'Member':8s}  {'Mode':>4}  {'Plane':>10}"
        f"  {'L_cr[m]':>8}  {'L_cr/L':>7}  {'N_cr[kN]':>9}"
    )
    print(SEP)

    for mem_name in all_members:
        member = model.members[mem_name]
        L_ref = ref_lengths[mem_name]

        try:
            N_Ed = member.axial(x=0.0, combo_name=COMBO_NAME)
        except Exception:
            N_Ed = float("nan")

        for mode_idx in range(n_modes):
            lam = lams[mode_idx]

            for plane_label, plane_char in (("y (strong)", "y"), ("z (weak)", "z")):
                try:
                    L_cr = results.effective_length(
                        mem_name, mode=mode_idx, plane=plane_char
                    )
                except Exception:
                    L_cr = float("nan")

                N_cr_kN = lam * abs(N_Ed) / 1e3

                print(
                    f"  {mem_name:8s}  {mode_idx + 1:4d}  {plane_label:>10}"
                    f"  {L_cr:8.3f}  {L_cr / L_ref:7.3f}  {N_cr_kN:9.1f}"
                )

        print()

    # --- Euler reference (pinned-pinned) -------------------------------------
    print(SEP)
    print("  Euler reference  (pinned-pinned, N_cr = pi^2*EI/L^2)")
    print(SEP)
    N_Ey_col = np.pi**2 * EIy / H**2 / 1e3
    N_Ez_col = np.pi**2 * EIz / H**2 / 1e3
    print(f"  Column  (L = H = {H:.2f} m):")
    print(f"    Strong axis (Iy):  N_cr = {N_Ey_col:.1f} kN")
    print(f"    Weak axis  (Iz):  N_cr = {N_Ez_col:.1f} kN")
    N_Ey_beam = np.pi**2 * EIy / B_SPAN**2 / 1e3
    N_Ez_beam = np.pi**2 * EIz / B_SPAN**2 / 1e3
    print(f"  Beam    (L = B_SPAN = {B_SPAN:.2f} m):")
    print(f"    Strong axis (Iy):  N_cr = {N_Ey_beam:.1f} kN")
    print(f"    Weak axis  (Iz):  N_cr = {N_Ez_beam:.1f} kN")
    print()

    # --- Member axial forces -------------------------------------------------
    print(SEP)
    print("  Member axial forces from static pre-solve")
    print(SEP)
    print(f"  {'Member':8s}  {'N_Ed [kN]':>12}  {'Sign'}")
    print(SEP)
    for mem_name in all_members:
        member = model.members[mem_name]
        try:
            N = member.axial(x=0.0, combo_name=COMBO_NAME)
            sign = "compression" if N > 0 else "tension"
            print(f"  {mem_name:8s}  {N / 1e3:12.2f}  {sign}")
        except Exception as exc:
            print(f"  {mem_name:8s}  {'N/A':>12}  ({exc})")
    print(SEP)
    print()

    # --- Beam & column internal forces ---------------------------------------
    print(SEP)
    print("  Member internal forces  (N_Ed at i-end, M_max along member)")
    print(SEP)
    print(f"  {'Member':8s}  {'N_Ed [kN]':>12}  {'M_max [kNm]':>12}  {'Sign_N'}")
    print(SEP)
    for mem_name in all_members:
        m = model.members[mem_name]
        try:
            N = m.axial(x=0.0, combo_name=COMBO_NAME)
            sign = "compr." if N > 0 else "tension"
            x_pts = [i * m.L() / 20 for i in range(21)]
            M_max = max(abs(m.moment("My", x, combo_name=COMBO_NAME)) for x in x_pts)
            print(f"  {mem_name:8s}  {N / 1e3:12.2f}  {M_max / 1e3:12.2f}  {sign}")
        except Exception as exc:
            print(f"  {mem_name:8s}  {'N/A':>12}  {'N/A':>12}  ({exc})")
    print(SEP)
    print()
    print("  Done.  Compare Lambda_cr and L_cr/L with reference output.")

    # --- FEM Design comparison -----------------------------------------------
    if SCENARIO in FEM_DESIGN_REF:
        ref = FEM_DESIGN_REF[SCENARIO]
        print()
        print(SEP)
        print("  FEM Design comparison")
        print(SEP)
        print(
            f"  {'Mode':>4}   {'PyNite':>10}   {'FEM Design':>10}"
            f"   {'Diff [%]':>8}   Status"
        )
        print(SEP)
        n_compare = min(len(lams), len(ref))
        all_pass = True
        for i in range(n_compare):
            diff_pct = (lams[i] - ref[i]) / ref[i] * 100
            status = "OK" if abs(diff_pct) < 5.0 else "CHECK"
            if abs(diff_pct) >= 5.0:
                all_pass = False
            print(
                f"  {i + 1:4d}   {lams[i]:10.3f}   {ref[i]:10.3f}"
                f"   {diff_pct:+8.2f}   {status}"
            )
        print(SEP)
        if all_pass:
            print("  RESULT: ALL MODES WITHIN 5% TOLERANCE")
        else:
            print("  RESULT: SOME MODES EXCEED 5% TOLERANCE")
        print()


if __name__ == "__main__":
    main()
