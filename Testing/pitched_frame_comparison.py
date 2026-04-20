"""
pitched_frame_comparison.py
===========================
Interactive boilerplate for comparing PyNite buckling results against
PolyFrame, FEM Design, or any other reference.

Set SCENARIO (see list below) and run — or set SCENARIO = 0 and edit the
PARAMETERS block manually.

Usage
-----
    uv run python Testing/pitched_frame_comparison.py

Structure topology
------------------
Single-bay pitched-roof (gable) portal frame in the global XY-plane (Z=0).
Columns are vertical (along Y).  Rafters are inclined from eave to ridge.
Column bases can be pinned or fixed.  Eave connections are rigid.

              R (ridge)
             / \\
            /   \\
     Raft1 /     \\ Raft2
          /       \\
         B ─ ─ ─ ─ C        (eave level, Y = H)
         |         |
   Col1  |         |  Col2
         |         |
         A         D
       (base)    (base)

    B_SPAN = 5.0 m  (horizontal distance A to D)
    H      = 3.0 m  (column height, A to B / D to C)
    RIDGE_H = 1.0 m (ridge rise above eave level)

    Ridge node R is at (B_SPAN/2, H + RIDGE_H, 0).
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
#  1  Gravity UDL    fixed bases, 10 kN/m UDL on rafters (global FY)
#  2  Midspan P      fixed bases, midspan point loads on rafters (global FY)
#  3  Eave FX loads  fixed bases, FX point loads at eave nodes B & C
#  4  Wind-like      fixed bases, column FX UDL + asymmetric rafter FY UDLs
#  5  Asym rafter    fixed bases, asymmetric rafter FY UDLs (no horiz. load)
#
SCENARIO = int(os.getenv("SCENARIO", "4"))

# ============================================================================
# PARAMETERS -- used when SCENARIO = 0, or as defaults overridden by presets
# ============================================================================

# --- Geometry (metres) ------------------------------------------------------
H = 3.0          # column height (eave level)
B_SPAN = 5.0     # horizontal span (distance between column bases)
RIDGE_H = 1.0    # ridge rise above eave level

# Number of FEM sub-elements per member (more = better accuracy, slower)
N_ELEM = 8

# --- Material: steel (SI units: N, m, Pa) -----------------------------------
E = 210e9    # Young's modulus
G = 80.769e9     # Shear modulus
NU = 0.3     # Poisson's ratio
RHO = 0.0    # density -- set to 7850 if you want self-weight included

# --- Cross-sections ---------------------------------------------------------
# Keys: (A [m^2], Iy [m^4], Iz [m^4], J [m^4], Asy [m^2], Asz [m^2])
# Asy = shear area for y-dir shear (bending about z / weak-axis) = flange area
# Asz = shear area for z-dir shear (bending about y / strong-axis) = web area
SECTIONS = {
    "IPE200": dict(A=28.48e-4, Iy=1943e-8, Iz=142e-8, J=7.0e-8, Asy=17.0e-4, Asz=14.02e-4),
    "IPE300": dict(A=53.8e-4, Iy=8356e-8, Iz=604e-8, J=20.1e-8, Asy=32.10e-4, Asz=25.67e-4),
    "HEA200": dict(A=53.8e-4, Iy=3692e-8, Iz=1336e-8, J=21.1e-8, Asy=40.0e-4, Asz=18.05e-4),
}

# Which section for each member type
COLUMN_SEC = "IPE200"
RAFTER_SEC = "IPE200"

# --- Boundary conditions ----------------------------------------------------
# True = restrained, False = free.  Order: (DX, DY, DZ, RX, RY, RZ)
FIXED_BASE = (True, True, True, True, True, True)    # fully fixed
PINNED_BASE = (True, True, True, True, True, False)   # pinned (RZ free in-plane)

BASE_SUPPORT = FIXED_BASE  # <-- switch to PINNED_BASE to try pinned bases

# --- Distributed loads on rafters (N/m) -------------------------------------
# Global FY direction.  Negative = downward.  Set to 0.0 to skip.
RAFTER_LOAD_1 = -10e3   # N/m on Raft1 (B->R, left rafter)
RAFTER_LOAD_2 = -10e3   # N/m on Raft2 (R->C, right rafter)

# --- Optional rafter midspan point loads (N) --------------------------------
# Global FY; negative = downward. Applied at x = L_rafter/2 from member i-end.
RAFTER_MID_LOAD_1 = 0.0  # Raft1 (B->R)
RAFTER_MID_LOAD_2 = 0.0  # Raft2 (R->C)

# --- Optional point loads at eave nodes (N, global FY) ----------------------
# Set to 0.0 to skip.
LOAD_B = 0.0   # left eave  (node B)
LOAD_C = 0.0   # right eave (node C)
LOAD_R = 0.0   # ridge      (node R)

# --- Optional point loads at nodes in global FX (N) -------------------------
LOAD_B_FX = 0.0   # node B
LOAD_C_FX = 0.0   # node C
LOAD_R_FX = 0.0   # node R

# --- Optional distributed loads on columns in global FX (N/m) --------------
COL_LOAD_1_FX = 0.0   # Col1 (A->B)
COL_LOAD_2_FX = 0.0   # Col2 (D->C)

# --- Lateral bracing / 2D constraint ----------------------------------------
# True  = restrain out-of-plane DOFs (DZ, RX, RY) at all nodes
#         forces pure in-plane buckling — use for braced frames
# False = full 3D analysis including lateral-torsional buckling
LATERAL_BRACE_NODES = True

# --- Beam element formulation -----------------------------------------------
# True  = Timoshenko beam (includes shear deformation via Asy/Asz)
# False = Euler-Bernoulli (classical, no shear deformation)
USE_TIMOSHENKO = True

# --- Buckling analysis settings --------------------------------------------
NUM_MODES = 5          # number of buckling modes to compute
COMBO_NAME = "Combo 1"

# ============================================================================
# END OF PARAMETERS
# ============================================================================

# ---------------------------------------------------------------------------
# FEM Design reference results (critical parameters / load multipliers)
# ---------------------------------------------------------------------------
# Each key is a scenario number; values are lists of lambda_cr per mode.
# Source: FEM Design 2024, IPE200 all members, E=210 GPa, G=80.769 GPa,
#         rigid line supports (out-of-plane restrained), member-length loads.

_FEM_DESIGN_REF = {
    1: [95.414, 236.049, 354.718, 465.246],
    2: [255.666, 589.676, 939.767, 1238.565],
    3: [2071.518, 2172.698],
    4: [980.953, 1172.103],
    5: [669.486, 769.573],
}


# ---------------------------------------------------------------------------
# Scenario presets
# ---------------------------------------------------------------------------

_PRESETS = {
    1: dict(
        name="Pitched frame, fixed bases, 10 kN/m UDL on rafters (gravity)",
        H=3.0,
        B_SPAN=5.0,
        RIDGE_H=1.0,
        COLUMN_SEC="IPE200",
        RAFTER_SEC="IPE200",
        _base="FIXED",
        RAFTER_LOAD_1=-10e3,
        RAFTER_LOAD_2=-10e3,
        RAFTER_MID_LOAD_1=0.0,
        RAFTER_MID_LOAD_2=0.0,
        LOAD_B=0.0,
        LOAD_C=0.0,
        LOAD_R=0.0,
        LOAD_B_FX=0.0,
        LOAD_C_FX=0.0,
        LOAD_R_FX=0.0,
        COL_LOAD_1_FX=0.0,
        COL_LOAD_2_FX=0.0,
        LATERAL_BRACE_NODES=True,
        N_ELEM=8,
        NUM_MODES=5,
    ),
    2: dict(
        name="Pitched frame, fixed bases, midspan point loads on rafters (global FY)",
        H=3.0,
        B_SPAN=5.0,
        RIDGE_H=1.0,
        COLUMN_SEC="IPE200",
        RAFTER_SEC="IPE200",
        _base="FIXED",
        RAFTER_LOAD_1=0.0,
        RAFTER_LOAD_2=0.0,
        RAFTER_MID_LOAD_1=-10e3,
        RAFTER_MID_LOAD_2=-10e3,
        LOAD_B=0.0,
        LOAD_C=0.0,
        LOAD_R=0.0,
        LOAD_B_FX=0.0,
        LOAD_C_FX=0.0,
        LOAD_R_FX=0.0,
        COL_LOAD_1_FX=0.0,
        COL_LOAD_2_FX=0.0,
        LATERAL_BRACE_NODES=True,
        N_ELEM=8,
        NUM_MODES=5,
    ),
    3: dict(
        name="Pitched frame, fixed bases, FX point loads at eave nodes",
        H=3.0,
        B_SPAN=5.0,
        RIDGE_H=1.0,
        COLUMN_SEC="IPE200",
        RAFTER_SEC="IPE200",
        _base="FIXED",
        RAFTER_LOAD_1=0.0,
        RAFTER_LOAD_2=0.0,
        RAFTER_MID_LOAD_1=0.0,
        RAFTER_MID_LOAD_2=0.0,
        LOAD_B=0.0,
        LOAD_C=0.0,
        LOAD_R=0.0,
        LOAD_B_FX=10e3,
        LOAD_C_FX=10e3,
        LOAD_R_FX=0.0,
        COL_LOAD_1_FX=0.0,
        COL_LOAD_2_FX=0.0,
        LATERAL_BRACE_NODES=True,
        N_ELEM=8,
        NUM_MODES=5,
    ),
    4: dict(
        name="Pitched frame, wind-like: column FX UDL + asymmetric rafter FY UDLs",
        H=3.0,
        B_SPAN=5.0,
        RIDGE_H=1.0,
        COLUMN_SEC="IPE200",
        RAFTER_SEC="IPE200",
        _base="FIXED",
        RAFTER_LOAD_1=-10e3,
        RAFTER_LOAD_2=+10e3,
        RAFTER_MID_LOAD_1=0.0,
        RAFTER_MID_LOAD_2=0.0,
        LOAD_B=0.0,
        LOAD_C=0.0,
        LOAD_R=0.0,
        LOAD_B_FX=0.0,
        LOAD_C_FX=0.0,
        LOAD_R_FX=0.0,
        COL_LOAD_1_FX=10e3,
        COL_LOAD_2_FX=10e3,
        LATERAL_BRACE_NODES=True,
        N_ELEM=8,
        NUM_MODES=5,
    ),
    5: dict(
        name="Pitched frame, asymmetric rafter FY UDLs (no horizontal load)",
        H=3.0,
        B_SPAN=5.0,
        RIDGE_H=1.0,
        COLUMN_SEC="IPE200",
        RAFTER_SEC="IPE200",
        _base="FIXED",
        RAFTER_LOAD_1=-10e3,
        RAFTER_LOAD_2=+10e3,
        RAFTER_MID_LOAD_1=0.0,
        RAFTER_MID_LOAD_2=0.0,
        LOAD_B=0.0,
        LOAD_C=0.0,
        LOAD_R=0.0,
        LOAD_B_FX=0.0,
        LOAD_C_FX=0.0,
        LOAD_R_FX=0.0,
        COL_LOAD_1_FX=0.0,
        COL_LOAD_2_FX=0.0,
        LATERAL_BRACE_NODES=True,
        N_ELEM=8,
        NUM_MODES=5,
    ),
    6: dict(
        name="Pitched frame, wind-like: asymmetric rafter FY UDLs",
        H=3.0,
        B_SPAN=5.0,
        RIDGE_H=1.0,
        COLUMN_SEC="IPE200",
        RAFTER_SEC="IPE200",
        _base="FIXED",
        RAFTER_LOAD_1=0.0,
        RAFTER_LOAD_2=0.0,
        RAFTER_MID_LOAD_1=0.0,
        RAFTER_MID_LOAD_2=0.0,
        LOAD_B=0.0,
        LOAD_C=0.0,
        LOAD_R=0.0,
        LOAD_B_FX=0.0,
        LOAD_C_FX=0.0,
        LOAD_R_FX=0.0,
        COL_LOAD_1_FX=10e3,
        COL_LOAD_2_FX=10e3,
        LATERAL_BRACE_NODES=True,
        N_ELEM=8,
        NUM_MODES=5,
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
        if key == "_base":
            g["BASE_SUPPORT"] = FIXED_BASE if val == "FIXED" else PINNED_BASE
        else:
            g[key] = val


# ---------------------------------------------------------------------------
# Model builder
# ---------------------------------------------------------------------------


def _add_column(
    model, name, base_node, top_node, mat, sec, n_elem, col_height, base_x, base_support
):
    """Add a vertical column with n_elem intermediate nodes."""
    h = col_height / n_elem
    for i in range(1, n_elem):
        nn = f"_{name}_int{i}"
        model.add_node(nn, base_x, i * h, 0.0)
    model.def_support(base_node, *base_support)
    model.add_member(name, base_node, top_node, mat, sec, rotation=90)


def _add_beam(
    model, name, node_i, node_j, mat, sec, n_elem, xi, yi, xj, yj
):
    """Add a beam with n_elem intermediate nodes.

    Intermediate nodes are placed by linear interpolation between
    (xi, yi, 0) and (xj, yj, 0).  Works for both horizontal and
    inclined members.
    """
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

    # Ridge coordinates
    rx = B_SPAN / 2.0
    ry = H + RIDGE_H

    # Column base nodes
    model.add_node("A", 0.0, 0.0, 0.0)
    model.add_node("D", B_SPAN, 0.0, 0.0)

    # Eave nodes (column tops)
    model.add_node("B", 0.0, H, 0.0)
    model.add_node("C", B_SPAN, H, 0.0)

    # Ridge node
    model.add_node("R", rx, ry, 0.0)

    # Columns
    _add_column(
        model, "Col1", "A", "B", "Steel", COLUMN_SEC, N_ELEM, H, 0.0, BASE_SUPPORT
    )
    _add_column(
        model, "Col2", "D", "C", "Steel", COLUMN_SEC, N_ELEM, H, B_SPAN, BASE_SUPPORT
    )

    # Rafters
    _add_beam(
        model, "Raft1", "B", "R", "Steel", RAFTER_SEC, N_ELEM,
        0.0, H, rx, ry,
    )
    _add_beam(
        model, "Raft2", "R", "C", "Steel", RAFTER_SEC, N_ELEM,
        rx, ry, B_SPAN, H,
    )

    # --- Distributed rafter loads (global FY) --------------------------------
    if RAFTER_LOAD_1 != 0.0:
        model.add_member_dist_load("Raft1", "FY", RAFTER_LOAD_1, RAFTER_LOAD_1)
    if RAFTER_LOAD_2 != 0.0:
        model.add_member_dist_load("Raft2", "FY", RAFTER_LOAD_2, RAFTER_LOAD_2)

    # --- Rafter midspan point loads (global FY) ------------------------------
    raft1_L = np.sqrt((rx - 0.0) ** 2 + (ry - H) ** 2)
    raft2_L = np.sqrt((B_SPAN - rx) ** 2 + (H - ry) ** 2)
    if RAFTER_MID_LOAD_1 != 0.0:
        model.add_member_pt_load("Raft1", "FY", RAFTER_MID_LOAD_1, raft1_L / 2)
    if RAFTER_MID_LOAD_2 != 0.0:
        model.add_member_pt_load("Raft2", "FY", RAFTER_MID_LOAD_2, raft2_L / 2)

    # --- Optional node point loads (global FY) -------------------------------
    if LOAD_B != 0.0:
        model.add_node_load("B", "FY", LOAD_B)
    if LOAD_C != 0.0:
        model.add_node_load("C", "FY", LOAD_C)
    if LOAD_R != 0.0:
        model.add_node_load("R", "FY", LOAD_R)

    # --- Optional node point loads (global FX) -------------------------------
    if LOAD_B_FX != 0.0:
        model.add_node_load("B", "FX", LOAD_B_FX)
    if LOAD_C_FX != 0.0:
        model.add_node_load("C", "FX", LOAD_C_FX)
    if LOAD_R_FX != 0.0:
        model.add_node_load("R", "FX", LOAD_R_FX)

    # --- Optional distributed column loads (global FX) -----------------------
    if COL_LOAD_1_FX != 0.0:
        model.add_member_dist_load("Col1", "FX", COL_LOAD_1_FX, COL_LOAD_1_FX)
    if COL_LOAD_2_FX != 0.0:
        model.add_member_dist_load("Col2", "FX", COL_LOAD_2_FX, COL_LOAD_2_FX)

    # --- Lateral bracing / 2D constraint --------------------------------------
    # When True, restrain out-of-plane DOFs (DZ, RX, RY) at ALL non-base
    # nodes — forces pure in-plane buckling (DX, DY, RZ active).
    # When False, full 3D analysis including lateral-torsional buckling.
    if LATERAL_BRACE_NODES:
        for node_name in model.nodes:
            if node_name not in ("A", "D"):  # bases already fully supported
                model.def_support(node_name, False, False, True, True, True, False)

    return model


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

SEP = "-" * 64


def _rafter_length():
    """Compute the geometric length of one rafter (symmetric frame)."""
    rx = B_SPAN / 2.0
    return np.sqrt(rx**2 + RIDGE_H**2)


def _print_section_info():
    for role, sec_name in (("Column", COLUMN_SEC), ("Rafter", RAFTER_SEC)):
        s = SECTIONS[sec_name]
        print(f"  {role} section : {sec_name}")
        print(f"    A  = {s['A'] * 1e4:.2f} cm^2")
        print(f"    Iy = {s['Iy'] * 1e8:.0f} cm^4  (strong axis)")
        print(f"    Iz = {s['Iz'] * 1e8:.0f} cm^4  (weak axis)")
        print(f"    J  = {s['J'] * 1e8:.2f} cm^4")
        if USE_TIMOSHENKO and s.get('Asy') and s.get('Asz'):
            print(f"    Asy= {s['Asy'] * 1e4:.2f} cm^2  (shear area, weak-axis bending)")
            print(f"    Asz= {s['Asz'] * 1e4:.2f} cm^2  (shear area, strong-axis bending)")
    beam_type = "Timoshenko (shear deformation included)" if USE_TIMOSHENKO else "Euler-Bernoulli (no shear deformation)"
    print(f"  Beam formulation : {beam_type}")
    print()


def _print_loads():
    printed = False
    if RAFTER_LOAD_1 != 0.0:
        print(
            f"  Raft1 (B->R) : w = {RAFTER_LOAD_1 / 1e3:+.2f} kN/m  (distributed FY)"
        )
        printed = True
    if RAFTER_LOAD_2 != 0.0:
        print(
            f"  Raft2 (R->C) : w = {RAFTER_LOAD_2 / 1e3:+.2f} kN/m  (distributed FY)"
        )
        printed = True
    if RAFTER_MID_LOAD_1 != 0.0:
        print(
            f"  Raft1 (B->R) : P = {RAFTER_MID_LOAD_1 / 1e3:+.1f} kN"
            f" @ midspan (global FY)"
        )
        printed = True
    if RAFTER_MID_LOAD_2 != 0.0:
        print(
            f"  Raft2 (R->C) : P = {RAFTER_MID_LOAD_2 / 1e3:+.1f} kN"
            f" @ midspan (global FY)"
        )
        printed = True
    for node, F in [("B", LOAD_B), ("C", LOAD_C), ("R", LOAD_R)]:
        if F != 0.0:
            print(f"  Node {node}       : FY = {F / 1e3:+.1f} kN  (point load)")
            printed = True
    for node, F in [("B", LOAD_B_FX), ("C", LOAD_C_FX), ("R", LOAD_R_FX)]:
        if F != 0.0:
            print(f"  Node {node}       : FX = {F / 1e3:+.1f} kN  (point load)")
            printed = True
    for col, w in [("Col1", COL_LOAD_1_FX), ("Col2", COL_LOAD_2_FX)]:
        if w != 0.0:
            print(f"  {col:14s}: w = {w / 1e3:+.2f} kN/m  (distributed FX)")
            printed = True
    if not printed:
        print("  (no loads applied)")
    brace = (
        "Yes — 2D constraint (DZ, RX, RY restrained at all nodes)"
        if LATERAL_BRACE_NODES
        else "No — full 3D analysis"
    )
    print(f"  Lateral bracing : {brace}")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    _apply_scenario()

    L_raft = _rafter_length()
    slope_deg = np.degrees(np.arctan2(RIDGE_H, B_SPAN / 2.0))

    print("=" * 64)
    print("  PyNite -- Pitched-Roof Portal Frame Buckling Analysis")
    print("=" * 64)
    print()
    if SCENARIO != 0:
        print(f"  Scenario {SCENARIO}: {_PRESETS[SCENARIO]['name']}")
        print()
    print(f"  Frame geometry : H={H:.2f} m, B_SPAN={B_SPAN:.2f} m, RIDGE_H={RIDGE_H:.2f} m")
    print(f"  Ridge at       : ({B_SPAN / 2:.2f}, {H + RIDGE_H:.2f}, 0)")
    print(f"  Rafter length  : {L_raft:.3f} m   slope = {slope_deg:.1f} deg")
    bc_label = "Fixed" if BASE_SUPPORT == FIXED_BASE else "Pinned"
    print(f"  Base condition : {bc_label}")
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
    all_members = ["Col1", "Col2", "Raft1", "Raft2"]
    ref_lengths = {
        "Col1": H,
        "Col2": H,
        "Raft1": L_raft,
        "Raft2": L_raft,
    }

    sec_col = SECTIONS[COLUMN_SEC]
    EIy_col = E * sec_col["Iy"]
    EIz_col = E * sec_col["Iz"]

    print(SEP)
    print(
        f"  Effective buckling lengths  (H = {H:.2f} m, L_raft = {L_raft:.3f} m)"
    )
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
    N_Ey = np.pi**2 * EIy_col / H**2 / 1e3
    N_Ez = np.pi**2 * EIz_col / H**2 / 1e3
    print(f"  Column  (L=H={H:.2f} m):")
    print(f"    Strong axis (Iy):  N_cr = {N_Ey:.1f} kN   L_eff/H = 1.00")
    print(f"    Weak axis  (Iz):  N_cr = {N_Ez:.1f} kN   L_eff/H = 1.00")
    sec_raft = SECTIONS[RAFTER_SEC]
    EIy_raft = E * sec_raft["Iy"]
    EIz_raft = E * sec_raft["Iz"]
    N_Ey_r = np.pi**2 * EIy_raft / L_raft**2 / 1e3
    N_Ez_r = np.pi**2 * EIz_raft / L_raft**2 / 1e3
    print(f"  Rafter  (L=L_raft={L_raft:.3f} m):")
    print(f"    Strong axis (Iy):  N_cr = {N_Ey_r:.1f} kN   L_eff/L = 1.00")
    print(f"    Weak axis  (Iz):  N_cr = {N_Ez_r:.1f} kN   L_eff/L = 1.00")
    print()

    # --- Column axial forces -------------------------------------------------
    col_names = ["Col1", "Col2"]
    print(SEP)
    print("  Column axial forces from static pre-solve")
    print(SEP)
    print(f"  {'Member':8s}  {'N_Ed [kN]':>12}  {'Sign'}")
    print(SEP)
    for col in col_names:
        member = model.members[col]
        try:
            N = member.axial(x=0.0, combo_name=COMBO_NAME)
            sign = "compression" if N > 0 else "tension"
            print(f"  {col:8s}  {N / 1e3:12.2f}  {sign}")
        except Exception as exc:
            print(f"  {col:8s}  {'N/A':>12}  ({exc})")
    print(SEP)
    print()

    # --- Rafter internal forces ----------------------------------------------
    rafter_names = ["Raft1", "Raft2"]
    print(SEP)
    print("  Rafter internal forces  (compare with PolyFrame / FEM Design)")
    print(SEP)
    print(f"  {'Member':8s}  {'N_Ed [kN]':>12}  {'M_max [kNm]':>12}  {'Sign_N'}")
    print(SEP)
    for rname in rafter_names:
        rm = model.members[rname]
        try:
            N = rm.axial(x=0.0, combo_name=COMBO_NAME)
            sign = "compr." if N > 0 else "tension"
            x_pts = [i * rm.L() / 20 for i in range(21)]
            M_max = max(abs(rm.moment("My", x, combo_name=COMBO_NAME)) for x in x_pts)
            print(f"  {rname:8s}  {N / 1e3:12.2f}  {M_max / 1e3:12.2f}  {sign}")
        except Exception as exc:
            print(f"  {rname:8s}  {'N/A':>12}  {'N/A':>12}  ({exc})")
    print(SEP)
    print()
    print("  Done.  Compare Lambda_cr and L_cr/L with PolyFrame / FEM Design output.")

    # --- FEM Design comparison ------------------------------------------------
    fem_ref = _FEM_DESIGN_REF.get(SCENARIO)
    if fem_ref:
        print()
        print(SEP)
        print("  FEM Design comparison")
        print(SEP)
        print(f"  {'Mode':>4}  {'PyNite':>10}  {'FEM Design':>10}  {'Diff':>8}")
        print(SEP)
        for i, lam in enumerate(lams):
            if i < len(fem_ref):
                ref = fem_ref[i]
                diff = (lam - ref) / ref * 100
                print(f"  {i + 1:4d}  {lam:10.3f}  {ref:10.3f}  {diff:+7.1f}%")
            else:
                print(f"  {i + 1:4d}  {lam:10.3f}  {'n/a':>10}")
        print(SEP)


if __name__ == "__main__":
    main()
