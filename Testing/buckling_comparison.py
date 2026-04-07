"""
buckling_comparison.py
======================
Interactive boilerplate for comparing PyNite buckling results against
PolyFrame, FEM Design, or any other reference.

Set SCENARIO (see list below) and run — or set SCENARIO = 0 and edit the
PARAMETERS block manually.

Usage
-----
    uv run python Testing/buckling_comparison.py

Structure topology
------------------
Two-bay / one-storey portal frame in the global XY-plane (Z=0).
All columns are vertical (along Y).  Beams are horizontal (along X).
Column bases can be pinned or fixed.  Beam-column connections are rigid.

    B ---------- C ---------- D
    |            |            |
    | Col1       | Col2       | Col3
    |            |            |
    A            E            F
   (base)      (base)       (base)

Loads: uniform distributed load (FY) on beams, with optional point loads
       at column-top nodes.

To model a single-bay frame set B2 = 0.0 (or pick a single-bay scenario).
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
#  1  PolyFrame ref  2-bay, fixed, 10 kN/m UDL  (default PolyFrame example)
#  2  PolyFrame+brc  same + DZ bracing at beam nodes -> in-plane alpha_cr
#  3  1-bay fixed    1-bay, fixed bases, 10 kN/m UDL
#  4  1-bay pinned   1-bay, pinned bases, 10 kN/m UDL
#  5  2-bay fixed P  2-bay, fixed bases, 100 kN column-top point loads
#  6  2-bay pinned P 2-bay, pinned bases, 100 kN column-top point loads
#  7  Taller fixed   2-bay, fixed, H=5 m, 10 kN/m UDL
#  8  IPE300 cols    2-bay, fixed, IPE300 columns / IPE200 beams, 20 kN/m UDL
#
#  --- Isolated column hand-check (compare FEM λ_cr with Euler formula) ---
#  9  Pin-pin col    pinned-pinned, K=1.0, L_cr = H
# 10  Fix-fix col    fixed-fixed no sway, K=0.5, L_cr = 0.5H
# 11  Cantilever     fixed-free, K=2.0, L_cr = 2H
# 12  Fix-pin col    fixed-pinned no sway, K=0.699, L_cr = 0.699H
#
#  --- Beam LTB hand-check (M_cr from EN 1993-1-1 Annex F formula) -------
# 13  SS beam M_cr   simply supported IPE200 beam, UDL, step-by-step M_cr
# 14  Midspan FY P   point loads at midspan of each beam segment (global FY)
# 15  Top FX loads   point loads at column-top nodes (global FX)
# 16  Wind-like      outer-column FX UDL + inverted beam FY UDLs
#
SCENARIO = int(os.getenv("SCENARIO", "2"))

# ============================================================================
# PARAMETERS -- used when SCENARIO = 0, or as defaults overridden by presets
# ============================================================================

# --- Geometry (metres) ------------------------------------------------------
H = 3.0  # column height
B1 = 5.0  # span of left bay  (between Col1 and Col2)
B2 = 5.0  # span of right bay (between Col2 and Col3)
# set B2 = 0.0 for a single-bay frame

# Number of FEM sub-elements per member (more = better accuracy, slower)
N_ELEM = 8

# --- Material: steel (SI units: N, m, Pa) -----------------------------------
E = 200e9  # Young's modulus
G = 80e9  # Shear modulus
NU = 0.3  # Poisson's ratio
RHO = 0.0  # density -- set to 7850 if you want self-weight included

# --- Cross-sections ---------------------------------------------------------
# Add more entries here and assign them per member below.
# Keys: (A [m^2], Iy [m^4], Iz [m^4], J [m^4])
SECTIONS = {
    "IPE200": dict(A=28.5e-4, Iy=1943e-8, Iz=142e-8, J=7.0e-8),
    "IPE300": dict(A=53.8e-4, Iy=8356e-8, Iz=604e-8, J=20.1e-8),
    "HEA200": dict(A=53.8e-4, Iy=3692e-8, Iz=1336e-8, J=21.1e-8),
}

# Which section for each member type
COLUMN_SEC = "IPE200"
BEAM_SEC = "IPE200"

# --- Boundary conditions ----------------------------------------------------
# True = restrained, False = free.  Order: (DX, DY, DZ, RX, RY, RZ)
FIXED_BASE = (True, True, True, True, True, True)  # fully fixed
PINNED_BASE = (True, True, True, True, True, False)  # pinned (RZ free in-plane)

BASE_SUPPORT = FIXED_BASE  # <-- switch to PINNED_BASE to try pinned bases

# --- Distributed loads on beams (N/m) ---------------------------------------
# Negative = downward.  Set to 0.0 to skip.
BEAM_LOAD_1 = -10e3  # N/m on Beam1 (left bay,  B->C)
BEAM_LOAD_2 = -10e3  # N/m on Beam2 (right bay, C->D) -- ignored if B2 = 0

# --- Optional point loads at column-top nodes (N) ---------------------------
# Set to 0.0 to skip.  These are IN ADDITION to the beam distributed loads.
LOAD_B = 0.0  # left column top   (node B)
LOAD_C = 0.0  # middle column top (node C)
LOAD_D = 0.0  # right column top  (node D) -- ignored if B2 = 0

# --- Optional beam point loads at midspan (N) -------------------------------
# Global FY; negative = downward. Applied at x = L/2 from member i-end.
BEAM_MID_LOAD_1 = 0.0  # Beam1 (B->C)
BEAM_MID_LOAD_2 = 0.0  # Beam2 (C->D) -- ignored if B2 = 0

# --- Optional point loads at column-top nodes in global FX (N) -------------
LOAD_B_FX = 0.0  # node B
LOAD_C_FX = 0.0  # node C
LOAD_D_FX = 0.0  # node D -- ignored if B2 = 0

# --- Optional distributed loads on columns in global FX (N/m) --------------
# Useful for wind-style loading. Positive is +X direction.
COL_LOAD_1_FX = 0.0  # Col1 (A->B)
COL_LOAD_2_FX = 0.0  # Col2 (E->C)
COL_LOAD_3_FX = 0.0  # Col3 (F->D) -- ignored if B2 = 0

# --- Lateral bracing at beam-level nodes ------------------------------------
# True  = restrain DZ (global Z) at nodes B, C, D
#         Matches PolyFrame's "Stability Supports: FZ" at beam level.
#         Prevents lateral sway of the whole frame; isolates LTB of beams.
# False = no lateral restraint -- full 3D eigenvalue of the whole frame.
#         This is the more complete structural analysis (recommended).
LATERAL_BRACE_BEAM_NODES = False

# --- Buckling analysis settings --------------------------------------------
NUM_MODES = 5  # number of buckling modes to compute
COMBO_NAME = "Combo 1"

# ============================================================================
# END OF PARAMETERS
# ============================================================================


# ---------------------------------------------------------------------------
# Scenario presets
# ---------------------------------------------------------------------------

_PRESETS = {
    1: dict(
        name="2-bay, fixed bases, 10 kN/m UDL on beams  [PolyFrame reference]",
        H=3.0,
        B1=5.0,
        B2=5.0,
        COLUMN_SEC="IPE200",
        BEAM_SEC="IPE200",
        _base="FIXED",
        BEAM_LOAD_1=-10e3,
        BEAM_LOAD_2=-10e3,
        LOAD_B=0.0,
        LOAD_C=0.0,
        LOAD_D=0.0,
        LATERAL_BRACE_BEAM_NODES=False,
        N_ELEM=8,
        NUM_MODES=5,
    ),
    2: dict(
        name="2-bay, fixed bases, midspan point loads on beams (global FY)",
        H=3.0,
        B1=5.0,
        B2=5.0,
        COLUMN_SEC="IPE200",
        BEAM_SEC="IPE200",
        _base="FIXED",
        BEAM_LOAD_1=0.0,
        BEAM_LOAD_2=0.0,
        BEAM_MID_LOAD_1=-10e3,
        BEAM_MID_LOAD_2=-10e3,
        LOAD_B=0.0,
        LOAD_C=0.0,
        LOAD_D=0.0,
        LOAD_B_FX=0.0,
        LOAD_C_FX=0.0,
        LOAD_D_FX=0.0,
        COL_LOAD_1_FX=0.0,
        COL_LOAD_2_FX=0.0,
        COL_LOAD_3_FX=0.0,
        LATERAL_BRACE_BEAM_NODES=False,
        N_ELEM=8,
        NUM_MODES=5,
    ),
    3: dict(
        name="2-bay, fixed bases, top-of-column point loads in global FX",
        H=3.0,
        B1=5.0,
        B2=5.0,
        COLUMN_SEC="IPE200",
        BEAM_SEC="IPE200",
        _base="FIXED",
        BEAM_LOAD_1=0.0,
        BEAM_LOAD_2=0.0,
        BEAM_MID_LOAD_1=0.0,
        BEAM_MID_LOAD_2=0.0,
        LOAD_B=0.0,
        LOAD_C=0.0,
        LOAD_D=0.0,
        LOAD_B_FX=10e3,
        LOAD_C_FX=10e3,
        LOAD_D_FX=10e3,
        COL_LOAD_1_FX=0.0,
        COL_LOAD_2_FX=0.0,
        COL_LOAD_3_FX=0.0,
        LATERAL_BRACE_BEAM_NODES=False,
        N_ELEM=8,
        NUM_MODES=5,
    ),
    4: dict(
        name="2-bay wind-like: outer-column FX UDL + inverted beam FY UDLs",
        H=3.0,
        B1=5.0,
        B2=5.0,
        COLUMN_SEC="IPE200",
        BEAM_SEC="IPE200",
        _base="FIXED",
        BEAM_LOAD_1=-10e3,
        BEAM_LOAD_2=+10e3,
        BEAM_MID_LOAD_1=0.0,
        BEAM_MID_LOAD_2=0.0,
        LOAD_B=0.0,
        LOAD_C=0.0,
        LOAD_D=0.0,
        LOAD_B_FX=0.0,
        LOAD_C_FX=0.0,
        LOAD_D_FX=0.0,
        COL_LOAD_1_FX=10e3,
        COL_LOAD_2_FX=0.0,
        COL_LOAD_3_FX=10e3,
        LATERAL_BRACE_BEAM_NODES=False,
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
    """Add a column with n_elem intermediate nodes."""
    h = col_height / n_elem
    for i in range(1, n_elem):
        nn = f"_{name}_int{i}"
        model.add_node(nn, base_x, i * h, 0.0)
    model.def_support(base_node, *base_support)
    model.add_member(name, base_node, top_node, mat, sec)


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
    model.add_member(name, node_i, node_j, mat, sec)


def build_model():
    model = FEModel3D()
    two_bay = B2 > 0.0

    model.add_material("Steel", E, G, NU, RHO)
    for sec_name, props in SECTIONS.items():
        model.add_section(sec_name, props["A"], props["Iy"], props["Iz"], props["J"])

    # Column base nodes
    model.add_node("A", 0.0, 0.0, 0.0)
    model.add_node("E", B1, 0.0, 0.0)
    if two_bay:
        model.add_node("F", B1 + B2, 0.0, 0.0)

    # Column top / beam nodes
    model.add_node("B", 0.0, H, 0.0)
    model.add_node("C", B1, H, 0.0)
    if two_bay:
        model.add_node("D", B1 + B2, H, 0.0)

    # Columns
    _add_column(
        model, "Col1", "A", "B", "Steel", COLUMN_SEC, N_ELEM, H, 0.0, BASE_SUPPORT
    )
    _add_column(
        model, "Col2", "E", "C", "Steel", COLUMN_SEC, N_ELEM, H, B1, BASE_SUPPORT
    )
    if two_bay:
        _add_column(
            model,
            "Col3",
            "F",
            "D",
            "Steel",
            COLUMN_SEC,
            N_ELEM,
            H,
            B1 + B2,
            BASE_SUPPORT,
        )

    # Beams
    _add_beam(model, "Beam1", "B", "C", "Steel", BEAM_SEC, N_ELEM, 0.0, H, B1, H)
    if two_bay:
        _add_beam(model, "Beam2", "C", "D", "Steel", BEAM_SEC, N_ELEM, B1, H, B1 + B2, H)

    # --- Distributed beam loads ----------------------------------------------
    if BEAM_LOAD_1 != 0.0:
        model.add_member_dist_load("Beam1", "FY", BEAM_LOAD_1, BEAM_LOAD_1)
    if two_bay and BEAM_LOAD_2 != 0.0:
        model.add_member_dist_load("Beam2", "FY", BEAM_LOAD_2, BEAM_LOAD_2)

    # --- Beam midspan point loads (global FY) -------------------------------
    if BEAM_MID_LOAD_1 != 0.0:
        model.add_member_pt_load("Beam1", "FY", BEAM_MID_LOAD_1, B1 / 2)
    if two_bay and BEAM_MID_LOAD_2 != 0.0:
        model.add_member_pt_load("Beam2", "FY", BEAM_MID_LOAD_2, B2 / 2)

    # --- Optional column-top point loads -------------------------------------
    if LOAD_B != 0.0:
        model.add_node_load("B", "FY", LOAD_B)
    if LOAD_C != 0.0:
        model.add_node_load("C", "FY", LOAD_C)
    if two_bay and LOAD_D != 0.0:
        model.add_node_load("D", "FY", LOAD_D)

    # --- Optional column-top point loads in global FX -----------------------
    if LOAD_B_FX != 0.0:
        model.add_node_load("B", "FX", LOAD_B_FX)
    if LOAD_C_FX != 0.0:
        model.add_node_load("C", "FX", LOAD_C_FX)
    if two_bay and LOAD_D_FX != 0.0:
        model.add_node_load("D", "FX", LOAD_D_FX)

    # --- Optional distributed column loads in global FX ---------------------
    if COL_LOAD_1_FX != 0.0:
        model.add_member_dist_load("Col1", "FX", COL_LOAD_1_FX, COL_LOAD_1_FX)
    if COL_LOAD_2_FX != 0.0:
        model.add_member_dist_load("Col2", "FX", COL_LOAD_2_FX, COL_LOAD_2_FX)
    if two_bay and COL_LOAD_3_FX != 0.0:
        model.add_member_dist_load("Col3", "FX", COL_LOAD_3_FX, COL_LOAD_3_FX)

    # --- Lateral bracing (PolyFrame-style FZ stability supports) -------------
    if LATERAL_BRACE_BEAM_NODES:
        brace_nodes = ["B", "C"] + (["D"] if two_bay else [])
        for node in brace_nodes:
            model.def_support(node, False, False, True, False, False, False)

    return model, two_bay


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

SEP = "-" * 64


def _print_section_info():
    for role, sec_name in (("Column", COLUMN_SEC), ("Beam", BEAM_SEC)):
        s = SECTIONS[sec_name]
        print(f"  {role} section : {sec_name}")
        print(f"    A  = {s['A'] * 1e4:.2f} cm^2")
        print(f"    Iy = {s['Iy'] * 1e8:.0f} cm^4  (strong axis)")
        print(f"    Iz = {s['Iz'] * 1e8:.0f} cm^4  (weak axis)")
        print(f"    J  = {s['J'] * 1e8:.2f} cm^4")
    print()


def _print_loads(two_bay):
    printed = False
    if BEAM_LOAD_1 != 0.0:
        print(f"  Beam1 (B->C) : w = {BEAM_LOAD_1 / 1e3:+.2f} kN/m  (distributed FY)")
        printed = True
    if two_bay and BEAM_LOAD_2 != 0.0:
        print(f"  Beam2 (C->D) : w = {BEAM_LOAD_2 / 1e3:+.2f} kN/m  (distributed FY)")
        printed = True
    if BEAM_MID_LOAD_1 != 0.0:
        print(
            f"  Beam1 (B->C) : P = {BEAM_MID_LOAD_1 / 1e3:+.1f} kN"
            f" @ midspan (global FY)"
        )
        printed = True
    if two_bay and BEAM_MID_LOAD_2 != 0.0:
        print(
            f"  Beam2 (C->D) : P = {BEAM_MID_LOAD_2 / 1e3:+.1f} kN"
            f" @ midspan (global FY)"
        )
        printed = True
    pt_loads = [("B", LOAD_B), ("C", LOAD_C)] + ([("D", LOAD_D)] if two_bay else [])
    for node, F in pt_loads:
        if F != 0.0:
            print(f"  Node {node}       : FY = {F / 1e3:+.1f} kN  (point load)")
            printed = True
    pt_loads_fx = [("B", LOAD_B_FX), ("C", LOAD_C_FX)] + (
        [("D", LOAD_D_FX)] if two_bay else []
    )
    for node, F in pt_loads_fx:
        if F != 0.0:
            print(f"  Node {node}       : FX = {F / 1e3:+.1f} kN  (point load)")
            printed = True

    col_fx = [("Col1", COL_LOAD_1_FX), ("Col2", COL_LOAD_2_FX)] + (
        [("Col3", COL_LOAD_3_FX)] if two_bay else []
    )
    for col, w in col_fx:
        if w != 0.0:
            print(f"  {col:14s}: w = {w / 1e3:+.2f} kN/m  (distributed FX)")
            printed = True
    if not printed:
        print("  (no loads applied)")
    brace = (
        "Yes (DZ restrained at B, C" + (", D" if two_bay else "") + ")"
        if LATERAL_BRACE_BEAM_NODES
        else "No"
    )
    print(f"  Lateral brace at beam nodes : {brace}")
    print()


# ---------------------------------------------------------------------------
# Isolated-column runner  (scenarios 9-12)
# ---------------------------------------------------------------------------


def _build_isolated_column(p):
    model = FEModel3D()
    sec = SECTIONS[p["section"]]
    model.add_material("Steel", E, G, NU, 0.0)
    model.add_section(p["section"], sec["A"], sec["Iy"], sec["Iz"], sec["J"])
    n, H = p["N_ELEM"], p["H"]
    dh = H / n
    for i in range(n + 1):
        model.add_node(f"N{i}", 0.0, i * dh, 0.0)
    model.def_support("N0", *p["bc_base"])
    if any(p["bc_top"]):
        model.def_support(f"N{n}", *p["bc_top"])
    model.add_member("Col", "N0", f"N{n}", "Steel", p["section"])
    model.add_node_load(f"N{n}", "FY", p["P"])
    return model


def main_column(p):
    sec = SECTIONS[p["section"]]
    H, P, K = p["H"], p["P"], p["K_theory"]
    EIy = E * sec["Iy"]
    EIz = E * sec["Iz"]
    L_cr_th = K * H
    N_cr_z = np.pi**2 * EIz / L_cr_th**2
    N_cr_y = np.pi**2 * EIy / L_cr_th**2

    print("=" * 64)
    print("  PyNite -- Isolated Column Buckling Verification")
    print("=" * 64)
    print()
    print(f"  Scenario {SCENARIO}: {p['name']}")
    print()
    print(f"  Column height : H  = {H:.2f} m")
    print(f"  Section       : {p['section']}")
    print(f"  Load          : P  = {abs(P) / 1e3:.0f} kN  (compression)")
    print(
        f"  E*Iy = {EIy / 1e3:.0f} kNm^2  (strong)    "
        f"E*Iz = {EIz / 1e3:.0f} kNm^2  (weak)"
    )
    print()

    print(SEP)
    print(f"  Theoretical reference   K = {K:.3f}  ->  L_cr = K*H = {L_cr_th:.3f} m")
    print(SEP)
    print(
        f"  N_cr (weak,  Iz) = pi^2*{EIz / 1e3:.0f}/{L_cr_th:.3f}^2"
        f" = {N_cr_z / 1e3:.2f} kN   lam_cr = {N_cr_z / abs(P):.4f}"
    )
    print(
        f"  N_cr (strong, Iy) = pi^2*{EIy / 1e3:.0f}/{L_cr_th:.3f}^2"
        f" = {N_cr_y / 1e3:.2f} kN   lam_cr = {N_cr_y / abs(P):.4f}"
    )
    print(SEP)
    print()

    model = _build_isolated_column(p)
    print("  Running buckling analysis (num_modes=4) ...")
    results = model.analyze_buckling(num_modes=4)
    lams = results.load_multipliers
    print()

    print(SEP)
    print(
        f"  {'Mode':>4}   {'lam_cr':>8}   {'N_cr [kN]':>10}  "
        f" {'L_cr,z [m]':>10}   {'K_z':>6}   "
        f"{'L_cr,y [m]':>10}   {'K_y':>6}"
    )
    print(SEP)
    for i, lam in enumerate(lams):
        N_cr = lam * abs(P)
        try:
            Lz = results.effective_length("Col", mode=i, plane="z")
        except Exception:
            Lz = float("nan")
        try:
            Ly = results.effective_length("Col", mode=i, plane="y")
        except Exception:
            Ly = float("nan")
        tag = "  <-- mode 1" if i == 0 else ""
        print(
            f"  {i + 1:4d}   {lam:8.4f}   {N_cr / 1e3:10.2f}"
            f"   {Lz:10.3f}   {Lz / H:6.3f}"
            f"   {Ly:10.3f}   {Ly / H:6.3f}{tag}"
        )
    print(SEP)
    print()

    if len(lams) > 0:
        Lz0 = results.effective_length("Col", mode=0, plane="z")
        err = abs(Lz0 - L_cr_th) / L_cr_th * 100
        print(
            f"  FEM L_cr,z (mode 1) = {Lz0:.4f} m   theory = {L_cr_th:.4f} m"
            f"   error = {err:.2f}%"
        )
    print()
    print("  Done.")


# ---------------------------------------------------------------------------
# Beam M_cr runner  (scenario 13)
# ---------------------------------------------------------------------------


def _mcr_annex_f(E, Iz, G, J, Iw, L, C1, k, kw):
    """EN 1993-1-1 Annex F three-factor M_cr formula."""
    kL = k * L
    pi2 = np.pi**2
    return (
        C1
        * (pi2 * E * Iz / kL**2)
        * np.sqrt((kw / k) ** 2 * Iw / Iz + kL**2 * G * J / (pi2 * E * Iz))
    )


def main_beam_mcr(p):
    sec = SECTIONS[p["section"]]
    L, w, C1, k, kw, Iw = p["L"], p["w"], p["C1"], p["k"], p["kw"], p["Iw"]
    Iz, J = sec["Iz"], sec["J"]
    EIz, GJ = E * Iz, G * J

    kL = k * L
    pi2 = np.pi**2

    M_max = w * L**2 / 8  # simply-supported UDL max moment
    M_cr_noIw = _mcr_annex_f(E, Iz, G, J, 0.0, L, C1, k, kw)
    M_cr_Iw = _mcr_annex_f(E, Iz, G, J, Iw, L, C1, k, kw)

    term_warp = (kw / k) ** 2 * Iw / Iz  # m²
    term_sv = kL**2 * GJ / (pi2 * EIz)  # m²
    coeff = pi2 * EIz / kL**2  # N  (the leading coefficient)

    print("=" * 64)
    print("  PyNite -- Beam M_cr  (EN 1993-1-1 Annex F formula)")
    print("=" * 64)
    print()
    print(f"  Scenario {SCENARIO}: {p['name']}")
    print()
    print(f"  Span    : L  = {L:.1f} m         Section : {p['section']}")
    print(f"  UDL     : w  = {w / 1e3:.1f} kN/m")
    print(f"  C1      : {C1:.3f}  (1.132 = UDL on simply supported beam)")
    print(f"  k / kw  : {k:.1f} / {kw:.1f}   (both ends free)")
    print()

    print(SEP)
    print("  Section properties (from model SECTIONS dict)")
    print(SEP)
    print(f"  Iz  = {Iz * 1e8:.0f} cm^4     E*Iz = {EIz / 1e3:.0f} kNm^2")
    print(f"  J   = {J * 1e8:.2f} cm^4    G*J  = {GJ / 1e3:.3f} kNm^2")
    print(f"  Iw  = {Iw * 1e12:.0f} cm^6   (warping constant, set in preset)")
    print(SEP)
    print()

    print(SEP)
    print("  Static moment  (hand formula: M = |w|*L^2/8)")
    print(SEP)
    print(f"  M_max = {abs(w) / 1e3:.1f}*{L:.1f}^2/8 = {abs(M_max) / 1e3:.4f} kNm")
    print(SEP)
    print()

    print(SEP)
    print("  Annex F step-by-step  (verify each line by hand)")
    print(SEP)
    print(f"  Leading coeff : C1*pi^2*E*Iz/(kL)^2")
    print(f"                = {C1:.3f}*pi^2*{EIz / 1e3:.0f}/{kL:.1f}^2")
    print(f"                = {C1 * coeff / 1e3:.4f} kN")
    print()
    print(f"  Terms under sqrt  [units = m^2]:")
    print(f"    Warping   (kw/k)^2*Iw/Iz    = {term_warp:.6e} m^2")
    print(f"    St.Venant (kL)^2*GJ/(pi^2*EIz) = {term_sv:.6e} m^2")
    print(
        f"    Sum (no Iw)                  = {term_sv:.6e} m^2   sqrt = {np.sqrt(term_sv):.5f} m"
    )
    print(
        f"    Sum (with Iw)                = {term_warp + term_sv:.6e} m^2   sqrt = {np.sqrt(term_warp + term_sv):.5f} m"
    )
    print()
    print(
        f"  M_cr (Iw = 0) = {C1 * coeff / 1e3:.4f} x {np.sqrt(term_sv):.5f} m = {abs(M_cr_noIw) / 1e3:.4f} kNm"
    )
    print(
        f"  M_cr (with Iw)= {C1 * coeff / 1e3:.4f} x {np.sqrt(term_warp + term_sv):.5f} m = {abs(M_cr_Iw) / 1e3:.4f} kNm"
    )
    print(SEP)
    print()

    eta_noIw = abs(M_max) / abs(M_cr_noIw)
    eta_Iw = abs(M_max) / abs(M_cr_Iw)
    pass_noIw = "OK    (<1.0)" if eta_noIw <= 1.0 else "FAILS (>1.0)"
    pass_Iw = "OK    (<1.0)" if eta_Iw <= 1.0 else "FAILS (>1.0)"

    print(SEP)
    print("  Utilization  M_max / M_cr  (no partial factors)")
    print(SEP)
    print(
        f"  Without Iw:  {abs(M_max) / 1e3:.4f} / {abs(M_cr_noIw) / 1e3:.4f}"
        f" = {eta_noIw:.3f}  {pass_noIw}"
    )
    print(
        f"  With Iw:     {abs(M_max) / 1e3:.4f} / {abs(M_cr_Iw) / 1e3:.4f}"
        f" = {eta_Iw:.3f}  {pass_Iw}"
    )
    print(SEP)
    print()
    print("  Note: PyNite uses 6-DOF elements (no warping DOF).")
    print("  LTB cannot be extracted from the eigenvalue; use the formula above.")
    print()
    print("  Done.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    _apply_scenario()
    two_bay = B2 > 0.0

    print("=" * 64)
    print("  PyNite -- Linear Buckling Analysis")
    print("=" * 64)
    print()
    if SCENARIO != 0:
        print(f"  Scenario {SCENARIO}: {_PRESETS[SCENARIO]['name']}")
        print()
    bay_str = f"  B1={B1:.2f} m, B2={B2:.2f} m" if two_bay else f"  B={B1:.2f} m"
    print(f"  Frame geometry : H={H:.2f} m,{bay_str}")
    bc_label = "Fixed" if BASE_SUPPORT == FIXED_BASE else "Pinned"
    print(f"  Base condition : {bc_label}")
    print(f"  Sub-elements   : {N_ELEM} per member")
    print()
    _print_section_info()
    print("  Applied loads")
    _print_loads(two_bay)

    # --- Build and run -------------------------------------------------------
    model, two_bay = build_model()

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
    col_names = ["Col1", "Col2"] + (["Col3"] if two_bay else [])
    sec_col = SECTIONS[COLUMN_SEC]
    EIy_col = E * sec_col["Iy"]
    EIz_col = E * sec_col["Iz"]

    print(SEP)
    print(f"  Effective buckling lengths  (column height H = {H:.2f} m)")
    print(SEP)
    print(
        f"  {'Member':8s}  {'Mode':>4}  {'Plane':>10}"
        f"  {'L_cr[m]':>8}  {'L_cr/H':>7}  {'N_cr[kN]':>9}"
    )
    print(SEP)

    for col in col_names:
        member = model.members[col]

        try:
            N_Ed = member.axial(x=0.0, combo_name=COMBO_NAME)
        except Exception:
            N_Ed = float("nan")

        for mode_idx in range(n_modes):
            lam = lams[mode_idx]

            for plane_label, plane_char in (("y (strong)", "y"), ("z (weak)", "z")):
                try:
                    L_cr = results.effective_length(
                        col, mode=mode_idx, plane=plane_char
                    )
                except Exception:
                    L_cr = float("nan")

                N_cr_kN = lam * abs(N_Ed) / 1e3

                print(
                    f"  {col:8s}  {mode_idx + 1:4d}  {plane_label:>10}"
                    f"  {L_cr:8.3f}  {L_cr / H:7.3f}  {N_cr_kN:9.1f}"
                )

        print()

    # --- Single-column Euler reference ---------------------------------------
    print(SEP)
    print("  Single-column Euler reference  (pinned-pinned, N_cr = pi^2*EI/H^2)")
    print(SEP)
    N_Ey = np.pi**2 * EIy_col / H**2 / 1e3
    N_Ez = np.pi**2 * EIz_col / H**2 / 1e3
    print(f"  Strong axis (Iy):  N_cr = {N_Ey:.1f} kN   L_eff/H = 1.00")
    print(f"  Weak axis  (Iz):  N_cr = {N_Ez:.1f} kN   L_eff/H = 1.00")
    print()

    # --- Column axial forces -------------------------------------------------
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

    # --- Beam internal forces ------------------------------------------------
    beam_names = ["Beam1"] + (["Beam2"] if two_bay else [])
    print(SEP)
    print("  Beam internal forces  (compare with PolyFrame N_Ed, M_Ed)")
    print(SEP)
    print(f"  {'Member':8s}  {'N_Ed [kN]':>12}  {'M_max [kNm]':>12}  {'Sign_N'}")
    print(SEP)
    for bname in beam_names:
        bm = model.members[bname]
        try:
            N = bm.axial(x=0.0, combo_name=COMBO_NAME)
            sign = "compr." if N > 0 else "tension"
            x_pts = [i * bm.L() / 20 for i in range(21)]
            M_max = max(abs(bm.moment("Mz", x, combo_name=COMBO_NAME)) for x in x_pts)
            print(f"  {bname:8s}  {N / 1e3:12.2f}  {M_max / 1e3:12.2f}  {sign}")
        except Exception as exc:
            print(f"  {bname:8s}  {'N/A':>12}  {'N/A':>12}  ({exc})")
    print(SEP)
    print()
    print("  Done.  Compare Lambda_cr and L_cr/H with PolyFrame / FEM Design output.")


if __name__ == "__main__":
    main()
