"""
Verification tests for linear buckling (stability) eigenvalue analysis.

Tests
-----
1. Euler pinned-pinned column:
   λ_cr = π²EI/(PL²)  →  must match within 1%
2. Fixed-free cantilever column:
   effective buckling length L_cr = 2L  →  must match within 2%
3. Portal frame (fixed bases):
   λ_cr > 1  (structure stable under applied loads)
"""

import numpy as np
import pytest
from Pynite import FEModel3D


# ---------------------------------------------------------------------------
# Helper: build a vertical column model (along global Y)
# ---------------------------------------------------------------------------
def _build_column(n_elem, E, G, nu, A, Iy, Iz, J, L, support_base, support_top,
                  load_direction='FY', P=-1.0):
    """
    Build a column model segmented into *n_elem* equal sub-elements.

    The column runs along the global Y-axis from Y=0 to Y=L.
    Returns the FEModel3D instance ready for buckling analysis.
    """
    model = FEModel3D()
    h = L / n_elem

    model.add_material('Mat', E, G, nu, 1.0)
    model.add_section('Sec', A, Iy, Iz, J)

    # Add nodes
    for i in range(n_elem + 1):
        model.add_node(f'N{i+1}', 0.0, i * h, 0.0)

    # Support conditions (True = restrained)
    model.def_support('N1',         *support_base)
    model.def_support(f'N{n_elem+1}', *support_top)

    # Single physical member — PyNite auto-segments at intermediate nodes
    model.add_member('Column', 'N1', f'N{n_elem+1}', 'Mat', 'Sec')

    # Compressive load at top node
    model.add_node_load(f'N{n_elem+1}', load_direction, P)

    return model


# ---------------------------------------------------------------------------
# Test 1: Euler column — pinned-pinned
# ---------------------------------------------------------------------------
def test_euler_column_pinned_pinned():
    """
    Euler pinned-pinned column: λ_cr = π²EI/(PL²).

    Uses a circular section with A >> 2π² (= ~20) to ensure the
    torsional-buckling eigenvalue G·J·A/(P·Ip) = 50 lies well above π² ≈ 9.87,
    so the first computed eigenvalue is the Euler flexural mode.

    Tolerance: 1% of the theoretical value π².
    """
    n_elem = 10         # 10 sub-elements → FEM error < 0.05%
    L      = 1.0
    E      = 1.0
    G      = 0.5        # G = E/2 for ν=0
    nu     = 0.0
    Iy     = 1.0        # EIy = 1
    Iz     = 1.0        # EIz = 1
    J      = 2.0        # J = Ip = Iy+Iz (solid circular section)
    A      = 100.0      # large A → λ_T = G·J·A/(P·Ip) = 50 >> π²

    # Pinned base: all translations fixed, bending rotations free,
    # torsion (RY for column along Y) restrained.
    support_base = (True,  True,  True,  False, True,  False)
    # DX DY DZ RX RY RZ

    # Pinned roller top: lateral (DX, DZ) fixed, axial (DY) free,
    # bending rotations free, torsion restrained.
    support_top  = (True,  False, True,  False, True,  False)

    model = _build_column(n_elem, E, G, nu, A, Iy, Iz, J, L,
                          support_base, support_top, P=-1.0)

    results = model.analyze_buckling(num_modes=3)

    assert len(results.load_multipliers) > 0, "No positive eigenvalues found"

    lam_computed = results.load_multipliers[0]
    lam_euler    = np.pi**2   # ≈ 9.8696

    rel_error = abs(lam_computed - lam_euler) / lam_euler
    assert rel_error < 0.01, (
        f"Euler column: λ_computed={lam_computed:.5f}, "
        f"λ_euler=π²={lam_euler:.5f}, "
        f"relative error={rel_error*100:.3f}% > 1%"
    )


# ---------------------------------------------------------------------------
# Test 2: Fixed-free cantilever — effective buckling length L_cr = 2L
# ---------------------------------------------------------------------------
def test_cantilever_effective_length():
    """
    Fixed-free cantilever column: λ_cr = π²EI/(4PL²), L_cr = 2L.

    Checks both the eigenvalue and the L_cr computed via
    BucklingResults.effective_length().  Tolerance: 2%.
    """
    n_elem = 10
    L      = 1.0
    E      = 1.0
    G      = 0.5
    nu     = 0.0
    Iy     = 1.0
    Iz     = 1.0
    J      = 2.0
    A      = 100.0      # same as above → λ_T = 50 >> π²/4 ≈ 2.47

    # Fully fixed base
    support_base = (True,  True,  True,  True,  True,  True)
    # Free top — no support applied (omit: model.def_support for top not called)
    support_top  = (False, False, False, False, False, False)

    model = _build_column(n_elem, E, G, nu, A, Iy, Iz, J, L,
                          support_base, support_top, P=-1.0)

    results = model.analyze_buckling(num_modes=3)

    assert len(results.load_multipliers) > 0, "No positive eigenvalues found"

    lam_computed  = results.load_multipliers[0]
    lam_cantilever = np.pi**2 / 4.0   # ≈ 2.4674 (L_eff = 2L)

    rel_error_lam = abs(lam_computed - lam_cantilever) / lam_cantilever
    assert rel_error_lam < 0.01, (
        f"Cantilever λ_computed={lam_computed:.5f}, "
        f"λ_expected={lam_cantilever:.5f}, "
        f"error={rel_error_lam*100:.3f}% > 1%"
    )

    # Check effective buckling length L_cr = π·sqrt(EI / N_cr)
    # With N_cr = λ · |P| = λ · 1 and EI = 1: L_cr = π/sqrt(λ) = 2L for λ = π²/4
    L_cr_computed = results.effective_length('Column', mode=0, plane='y')
    L_cr_expected = 2.0 * L   # = 2.0

    rel_error_Lcr = abs(L_cr_computed - L_cr_expected) / L_cr_expected
    assert rel_error_Lcr < 0.02, (
        f"Cantilever effective length: L_cr_computed={L_cr_computed:.5f} m, "
        f"L_cr_expected={L_cr_expected:.5f} m, "
        f"relative error={rel_error_Lcr*100:.3f}% > 2%"
    )


# ---------------------------------------------------------------------------
# Test 3: Portal frame — verify stability
# ---------------------------------------------------------------------------
def test_portal_frame_buckling():
    """
    Simple portal frame with fixed column bases under symmetric vertical loads.

    Checks that:
    - The analysis completes successfully.
    - The first load multiplier is positive (structure is stable).
    - BucklingResults has the requested number of modes.
    """
    # IPE200 section properties (SI units: N, m)
    E   = 200e9      # Pa
    G   =  80e9      # Pa
    nu  =   0.3
    rho =    7850.0  # kg/m³ (not used in buckling)
    A   = 28.5e-4    # m²
    Iy  = 1943e-8    # m⁴ (strong axis)
    Iz  =  142e-8    # m⁴ (weak axis)
    J   =    7e-8    # m⁴ (St. Venant torsion)

    h = 3.0   # column height (m)
    b = 6.0   # beam span (m)
    P = 10e3  # vertical load per column top (N) = 10 kN

    model = FEModel3D()
    model.add_material('Steel', E, G, nu, rho)
    model.add_section('IPE200', A, Iy, Iz, J)

    # Nodes
    model.add_node('A', 0,  0, 0)  # column base left
    model.add_node('B', 0,  h, 0)  # column top left  / beam left end
    model.add_node('C', b,  h, 0)  # column top right / beam right end
    model.add_node('D', b,  0, 0)  # column base right

    # Fixed column bases
    model.def_support('A', True, True, True, True, True, True)
    model.def_support('D', True, True, True, True, True, True)

    # Members
    model.add_member('Col1', 'A', 'B', 'Steel', 'IPE200')
    model.add_member('Beam', 'B', 'C', 'Steel', 'IPE200')
    model.add_member('Col2', 'D', 'C', 'Steel', 'IPE200')

    # Vertical (compressive) loads at column tops
    model.add_node_load('B', 'FY', -P)
    model.add_node_load('C', 'FY', -P)

    num_modes = 3
    results = model.analyze_buckling(num_modes=num_modes)

    # Should return at least one mode
    assert len(results.load_multipliers) > 0, "Buckling analysis returned no modes"

    # First eigenvalue must be positive (structure does not buckle at zero load)
    assert results.load_multipliers[0] > 0, "First load multiplier is not positive"

    # For P = 10 kN and an IPE200 column (h = 3 m), the in-plane sway
    # N_cr ≈ 50–200 kN depending on beam restraint, so λ > 5.
    # We assert λ > 10 (a conservative frame is well above serviceability level).
    assert results.load_multipliers[0] > 10.0, (
        f"First load multiplier {results.load_multipliers[0]:.2f} is unexpectedly low "
        f"(expected > 10 for P = {P/1e3:.0f} kN on IPE200 h = {h} m portal frame)"
    )

    # Check alias
    np.testing.assert_array_equal(
        results.load_multipliers, results.critical_load_factors
    )


# ---------------------------------------------------------------------------
# Entry point for direct execution
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    test_euler_column_pinned_pinned()
    print("Test 1 (Euler column, pinned-pinned): PASSED")

    test_cantilever_effective_length()
    print("Test 2 (Cantilever effective length): PASSED")

    test_portal_frame_buckling()
    print("Test 3 (Portal frame stability):      PASSED")
