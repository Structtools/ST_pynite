"""
Verification tests for linear buckling (stability) eigenvalue analysis.

In-plane (flexural) buckling
----------------------------
1. Euler pinned-pinned column: λ_cr = π²EI/(PL²), within 1%.
2. Fixed-free cantilever column: effective length L_cr = 2L, within 2%.
3. Portal frame (fixed bases): λ_cr > 1 and the analysis runs.
4. Portal frame with a rigid beam: column effective length factors K = 1.0 (fixed bases) and
   K = 2.0 (pinned bases), within 1%.

Out-of-plane buckling
---------------------
5. Column with unequal second moments of area: the 3D analysis finds the weak-axis Euler load,
   labels the mode out-of-plane, and the in-plane / out-of-plane analyses each return their own
   Euler load.
6. Lateral-torsional buckling without warping, against the closed-form solutions of Timoshenko &
   Gere, "Theory of Elastic Stability", Chapter 6 (all with Iw = 0):
   - simply supported beam, uniform moment:      M_cr = (π/L)·sqrt(E·Iy·G·J)
   - simply supported beam, central point load:  P_cr = 16.94·sqrt(E·Iy·G·J)/L²
   - simply supported beam, uniform load:        q_cr = 28.3·sqrt(E·Iy·G·J)/L³
   - cantilever, tip load:                        P_cr = 4.013·sqrt(E·Iy·G·J)/L²
   - laterally and torsionally clamped ends, uniform moment: M_cr = (2π/L)·sqrt(E·Iy·G·J)
7. Equivalent torsion constant approximation of warping: exact for the fork-supported beam in
   uniform bending, M_cr = (π/L)·sqrt(E·Iy·(G·J + π²·E·Iw/L²)).
8. Consistency: a beam whose second half is rotated 90° about its axis (with Iy and Iz swapped)
   gives the same critical load as the unrotated beam. This exercises the My terms and the
   local-to-global transformation of the moment terms.
9. Planar frame: the 3D modes are the union of the in-plane and out-of-plane modes, and for an
   unbraced IPE200 portal frame lateral-torsional buckling governs by a wide margin. The
   axial-force-only matrix misses it.

Units: the beam tests are non-dimensional (E = G = 1, L = 1). The frame tests use N and m.
"""

import numpy as np
import pytest
from Pynite import FEModel3D

# Closed-form lateral-torsional buckling coefficients for Iw = 0 (Timoshenko & Gere, Chapter 6)
GAMMA_CENTRAL_POINT_LOAD = 16.94   # P_cr L² / sqrt(EI GJ), simply supported, load at the shear centre
GAMMA_UNIFORM_LOAD = 28.3          # q_cr L³ / sqrt(EI GJ), simply supported, load at the shear centre
GAMMA_CANTILEVER_TIP_LOAD = 4.013  # P_cr L² / sqrt(EI GJ), cantilever, load at the shear centre

# IPE200 (SI units: N, m). Iy is the strong axis in Eurocode notation; members are rotated by 90°
# in the frame models so that the strong axis bends in the XY plane.
E_STEEL = 210e9
G_STEEL = 80.77e9
IPE200 = dict(A=28.48e-4, Iy=1943e-8, Iz=142e-8, J=7.0e-8, Iw=12990e-12)

# Fork support DOFs for a beam along X: lateral displacement and twist restrained, bending
# rotations free. The i-end also carries the axial restraint.
FORK_I = (True, True, True, True, False, False)
FORK_J = (False, True, True, True, False, False)


# ---------------------------------------------------------------------------
# Helpers
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


def _beam(L=1.0, Iy=1.0, Iz=10.0, J=1.0, Iw=None, A=1e3):
    """A single-member beam along X with E = G = 1, lateral (weak) axis Iy and strong axis Iz."""
    model = FEModel3D()
    model.add_material('Mat', 1.0, 1.0, 0.0, 0.0)
    model.add_section('Sec', A, Iy, Iz, J, Iw=Iw)
    model.add_node('N0', 0.0, 0.0, 0.0)
    model.add_node('N1', L, 0.0, 0.0)
    model.add_member('Beam', 'N0', 'N1', 'Mat', 'Sec')
    return model


def _uniform_moment_beam(support_i=FORK_I, support_j=FORK_J, **kwargs):
    """A beam under equal and opposite end moments of 1 about Z (uniform in-plane moment)."""
    model = _beam(**kwargs)
    model.def_support('N0', *support_i)
    model.def_support('N1', *support_j)
    model.add_node_load('N0', 'MZ', 1.0)
    model.add_node_load('N1', 'MZ', -1.0)
    return model


def _portal_frame(brace_eaves=False, w=-10e3):
    """A 5 m x 3 m IPE200 portal frame in the XY plane with fixed bases, strong axis in-plane."""
    model = FEModel3D()
    model.add_material('Steel', E_STEEL, G_STEEL, 0.3, 0.0)
    model.add_section('IPE200', **IPE200)

    model.add_node('A', 0.0, 0.0, 0.0)
    model.add_node('B', 0.0, 3.0, 0.0)
    model.add_node('C', 5.0, 3.0, 0.0)
    model.add_node('D', 5.0, 0.0, 0.0)

    # Rotating the members by 90° puts the Eurocode strong axis (Iy) in the plane of the frame
    model.add_member('Col1', 'A', 'B', 'Steel', 'IPE200', rotation=90)
    model.add_member('Beam', 'B', 'C', 'Steel', 'IPE200', rotation=90)
    model.add_member('Col2', 'D', 'C', 'Steel', 'IPE200', rotation=90)

    model.def_support('A', True, True, True, True, True, True)
    model.def_support('D', True, True, True, True, True, True)

    if brace_eaves:
        # Lateral and torsional restraint at the eaves
        model.def_support('B', False, False, True, True, False, False)
        model.def_support('C', False, False, True, True, False, False)

    model.add_member_dist_load('Beam', 'FY', w, w)

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
# Test 4: In-plane sway buckling of a portal frame with a rigid beam
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('fixed_base, K', [(True, 1.0), (False, 2.0)])
def test_rigid_beam_portal_effective_length_factor(fixed_base, K):
    """
    A sway portal frame whose beam is rigid compared to its columns has the textbook column
    effective length factors K = 1.0 (fixed bases) and K = 2.0 (pinned bases).

    Checks the critical load factor and the derived buckling length, within 1%.
    """
    E, I, h = 1.0, 1.0, 1.0

    model = FEModel3D()
    model.add_material('Mat', E, E/2.6, 0.3, 0.0)
    model.add_section('Col', 1e3, I, I, 2*I)
    model.add_section('Rigid', 1e3, 1e5*I, 1e5*I, 2e5*I)

    model.add_node('A', 0.0, 0.0, 0.0)
    model.add_node('B', 0.0, h, 0.0)
    model.add_node('C', 2.0, h, 0.0)
    model.add_node('D', 2.0, 0.0, 0.0)

    model.add_member('Col1', 'A', 'B', 'Mat', 'Col')
    model.add_member('Beam', 'B', 'C', 'Mat', 'Rigid')
    model.add_member('Col2', 'D', 'C', 'Mat', 'Col')

    model.def_support('A', True, True, True, True, True, fixed_base)
    model.def_support('D', True, True, True, True, True, fixed_base)

    model.add_node_load('B', 'FY', -1.0)
    model.add_node_load('C', 'FY', -1.0)

    results = model.analyze_buckling(num_modes=2, plane='XY')

    N_cr = np.pi**2*E*I/(K*h)**2
    assert results.load_multipliers[0] == pytest.approx(N_cr, rel=0.01)

    # In-plane bending of a vertical member is about its local z-axis
    assert results.effective_length('Col1', mode=0, plane='z') == pytest.approx(K*h, rel=0.01)

    # The first-order member forces of the analysed combination are available on the model
    assert abs(model.members['Col1'].axial(0.0, 'Combo 1')) == pytest.approx(1.0, rel=1e-6)


# ---------------------------------------------------------------------------
# Test 5: Weak-axis (out-of-plane) flexural buckling of a column
# ---------------------------------------------------------------------------
def test_weak_axis_column_buckles_out_of_plane():
    """
    A pin-ended column in the XY plane with Iy = 1 (bending out of plane, in Z) and Iz = 4
    (bending in plane). The 3D analysis must find the weak-axis Euler load π² first and label it
    out-of-plane; the in-plane analysis must find 4π² and the out-of-plane analysis π².

    A = 1000 and J = Ip keep the torsional eigenvalue G·J·A/(P·Ip) = 1000 far above the flexural ones.
    """
    # Pinned about both axes at both ends, twist (RY for a column along Y) restrained
    support_base = (True, True, True, False, True, False)
    support_top = (True, False, True, False, True, False)

    model = _build_column(1, E=1.0, G=1.0, nu=0.0, A=1e3, Iy=1.0, Iz=4.0, J=5.0, L=1.0,
                          support_base=support_base, support_top=support_top, P=-1.0)

    results_3d = model.analyze_buckling(num_modes=2)
    assert results_3d.load_multipliers[0] == pytest.approx(np.pi**2, rel=0.005)
    assert results_3d.classify_modes('XY')[0] == 'out-of-plane'
    assert results_3d.out_of_plane_share('XY')[0] == pytest.approx(1.0, abs=1e-9)

    results_in = model.analyze_buckling(num_modes=1, plane='XY')
    assert results_in.load_multipliers[0] == pytest.approx(4*np.pi**2, rel=0.005)
    assert results_in.classify_modes()[0] == 'in-plane'

    results_out = model.analyze_buckling(num_modes=1, plane='XY', out_of_plane=True)
    assert results_out.load_multipliers[0] == pytest.approx(np.pi**2, rel=0.005)
    assert results_out.classify_modes()[0] == 'out-of-plane'
    assert results_out.effective_length('Column', mode=0, plane='y') == pytest.approx(1.0, rel=0.005)


# ---------------------------------------------------------------------------
# Test 6: Lateral-torsional buckling against closed-form solutions (Iw = 0)
# ---------------------------------------------------------------------------
def test_ltb_simply_supported_uniform_moment():
    """Fork-supported beam under uniform moment: M_cr = (π/L)·sqrt(E·Iy·G·J) = π."""
    model = _uniform_moment_beam()

    results = model.analyze_buckling(num_modes=1, elements_per_member=8)

    assert results.load_multipliers[0] == pytest.approx(np.pi, rel=0.01)

    # The buckled shape is lateral displacement and twist only
    assert results.classify_modes('XY') == ['out-of-plane']


def test_ltb_simply_supported_central_point_load():
    """Fork-supported beam with a central point load at the shear centre: P_cr·L²/sqrt(EI·GJ) = 16.94."""
    model = _beam()
    model.def_support('N0', *FORK_I)
    model.def_support('N1', *FORK_J)
    model.add_member_pt_load('Beam', 'FY', -1.0, 0.5)

    results = model.analyze_buckling(num_modes=1, elements_per_member=16)

    assert results.load_multipliers[0] == pytest.approx(GAMMA_CENTRAL_POINT_LOAD, rel=0.005)


def test_ltb_simply_supported_uniform_load():
    """Fork-supported beam with a uniform load at the shear centre: q_cr·L³/sqrt(EI·GJ) = 28.3."""
    model = _beam()
    model.def_support('N0', *FORK_I)
    model.def_support('N1', *FORK_J)
    model.add_member_dist_load('Beam', 'FY', -1.0, -1.0)

    results = model.analyze_buckling(num_modes=1, elements_per_member=16)

    assert results.load_multipliers[0] == pytest.approx(GAMMA_UNIFORM_LOAD, rel=0.01)


def test_ltb_cantilever_tip_load():
    """Cantilever with a tip load at the shear centre: P_cr·L²/sqrt(EI·GJ) = 4.013."""
    model = _beam()
    model.def_support('N0', True, True, True, True, True, True)
    model.add_node_load('N1', 'FY', -1.0)

    results = model.analyze_buckling(num_modes=1, elements_per_member=16)

    assert results.load_multipliers[0] == pytest.approx(GAMMA_CANTILEVER_TIP_LOAD, rel=0.005)


def test_ltb_clamped_ends_uniform_moment():
    """
    Beam with lateral rotation and twist restrained at both ends (in-plane rotation free) under
    uniform moment. With Iw = 0 the lateral buckling length is L/2: M_cr = (2π/L)·sqrt(E·Iy·G·J).
    """
    model = _uniform_moment_beam(support_i=(True, True, True, True, True, False),
                                 support_j=(False, True, True, True, True, False))

    results = model.analyze_buckling(num_modes=1, elements_per_member=16)

    assert results.load_multipliers[0] == pytest.approx(2*np.pi, rel=0.01)


def test_axial_only_geometric_stiffness_misses_lateral_torsional_buckling():
    """
    Without the moment terms the cantilever carries no axial force, so the geometric stiffness
    matrix is empty and no buckling mode exists. This is the failure mode the moment terms fix.
    """
    model = _beam()
    model.def_support('N0', True, True, True, True, True, True)
    model.add_node_load('N1', 'FY', -1.0)

    with pytest.raises(ValueError, match='geometric stiffness matrix is zero'):
        model.analyze_buckling(num_modes=1, include_moments=False)


# ---------------------------------------------------------------------------
# Test 7: Equivalent torsion constant approximation of warping
# ---------------------------------------------------------------------------
def test_ltb_equivalent_torsion_constant_matches_closed_form_with_warping():
    """
    For a fork-supported beam in uniform bending the equivalent torsion constant
    J_eff = J + π²·E·Iw/(G·L²) reproduces the closed-form critical moment with warping exactly:
    M_cr = (π/L)·sqrt(E·Iy·(G·J + π²·E·Iw/L²)).
    """
    Iw = 0.5
    model = _uniform_moment_beam(Iw=Iw)

    results = model.analyze_buckling(num_modes=1, elements_per_member=16,
                                     warping='equivalent_torsion')

    M_cr = np.pi*np.sqrt(1.0*(1.0 + np.pi**2*Iw))
    assert results.load_multipliers[0] == pytest.approx(M_cr, rel=0.005)
    assert results.warping == 'equivalent_torsion'

    # The user's model is untouched by the override
    assert all(member._J_eff is None for member in model.members.values())


def test_equivalent_torsion_requires_warping_constant():
    """The warping approximation needs Iw on every section."""
    model = _uniform_moment_beam()

    with pytest.raises(ValueError, match='Iw'):
        model.analyze_buckling(warping='equivalent_torsion')


# ---------------------------------------------------------------------------
# Test 8: Rotated member consistency
# ---------------------------------------------------------------------------
def test_rotated_member_gives_same_critical_load():
    """
    The second half of a fork-supported beam with a central point load is rotated 90° about its
    axis, with Iy and Iz swapped so that the physical beam is unchanged. The in-plane moment is
    then My in the rotated half and Mz in the other, and the critical load must not change.
    """
    rotated = FEModel3D()
    rotated.add_material('Mat', 1.0, 1.0, 0.0, 0.0)
    rotated.add_section('Sec', 1e3, 1.0, 10.0, 1.0)
    rotated.add_section('SecRot', 1e3, 10.0, 1.0, 1.0)
    rotated.add_node('N0', 0.0, 0.0, 0.0)
    rotated.add_node('N1', 0.5, 0.0, 0.0)
    rotated.add_node('N2', 1.0, 0.0, 0.0)
    rotated.add_member('B1', 'N0', 'N1', 'Mat', 'Sec')
    rotated.add_member('B2', 'N1', 'N2', 'Mat', 'SecRot', rotation=90)
    rotated.def_support('N0', *FORK_I)
    rotated.def_support('N2', *FORK_J)
    rotated.add_node_load('N1', 'FY', -1.0)

    straight = _beam()
    straight.def_support('N0', *FORK_I)
    straight.def_support('N1', *FORK_J)
    straight.add_member_pt_load('Beam', 'FY', -1.0, 0.5)

    lam_rotated = rotated.analyze_buckling(num_modes=1, elements_per_member=8).load_multipliers[0]
    lam_straight = straight.analyze_buckling(num_modes=1, elements_per_member=16).load_multipliers[0]

    assert lam_rotated == pytest.approx(lam_straight, rel=1e-6)


# ---------------------------------------------------------------------------
# Test 9: Planar frame — in-plane, out-of-plane and 3D analyses
# ---------------------------------------------------------------------------
def test_portal_frame_modes_are_union_of_in_plane_and_out_of_plane_modes():
    """
    For a planar frame loaded in its plane the in-plane and out-of-plane DOFs are uncoupled, so
    the modes of the 3D analysis are the union of the modes of the two restricted analyses.
    """
    model = _portal_frame()

    # The unbraced frame's first several modes are all out-of-plane, so the restricted analyses
    # together have to supply at least as many modes as the 3D analysis is asked for
    results_3d = model.analyze_buckling(num_modes=6)
    results_in = model.analyze_buckling(num_modes=3, plane='XY')
    results_out = model.analyze_buckling(num_modes=6, plane='XY', out_of_plane=True)

    union = np.sort(np.concatenate([results_in.load_multipliers, results_out.load_multipliers]))
    np.testing.assert_allclose(results_3d.load_multipliers, union[:6], rtol=1e-6)
    assert results_3d.classify_modes('XY') == ['out-of-plane']*6

    # Asking the 3D analysis for enough modes brings the in-plane modes in at their place
    results_3d_many = model.analyze_buckling(num_modes=40)
    in_plane_3d = results_3d_many.load_multipliers[
        np.array(results_3d_many.classify_modes('XY')) == 'in-plane']
    np.testing.assert_allclose(in_plane_3d[:3], results_in.load_multipliers, rtol=1e-6)

    assert results_in.classify_modes() == ['in-plane']*3
    assert results_out.classify_modes() == ['out-of-plane']*6
    assert results_in.plane == 'XY' and not results_in.out_of_plane
    assert results_out.plane == 'XY' and results_out.out_of_plane


def test_unbraced_portal_frame_governed_by_lateral_torsional_buckling():
    """
    An unbraced IPE200 portal frame spanning 5 m under 10 kN/m carries about 30 kNm in the beam,
    close to the beam's closed-form critical moment of about 26 kNm (Iw = 0, fork supports). The
    out-of-plane critical load factor must therefore be near 1.5 while the in-plane factor is in
    the hundreds, and the axial-force-only matrix must overestimate the out-of-plane factor.
    """
    model = _portal_frame()

    lam_in = model.analyze_buckling(num_modes=1, plane='XY').load_multipliers[0]
    lam_out = model.analyze_buckling(num_modes=1, plane='XY', out_of_plane=True).load_multipliers[0]
    lam_out_axial_only = model.analyze_buckling(num_modes=1, plane='XY', out_of_plane=True,
                                                include_moments=False).load_multipliers[0]
    lam_out_warping = model.analyze_buckling(num_modes=1, plane='XY', out_of_plane=True,
                                             warping='equivalent_torsion').load_multipliers[0]
    lam_out_braced = _portal_frame(brace_eaves=True).analyze_buckling(
        num_modes=1, plane='XY', out_of_plane=True).load_multipliers[0]

    assert 1.0 < lam_out < 2.5
    assert lam_in > 50.0
    assert lam_out_axial_only > 1.5*lam_out
    assert lam_out < lam_out_warping < 1.2*lam_out
    assert lam_out_braced > lam_out


# ---------------------------------------------------------------------------
# Results handling and argument validation
# ---------------------------------------------------------------------------
def test_mode_shapes_are_stored_as_load_combinations():
    """Mode shapes land on the user's model as 'Buckling Mode n' combinations tagged 'buckling'."""
    model = _portal_frame()

    results = model.analyze_buckling(num_modes=3)
    results = model.analyze_buckling(num_modes=3)   # A repeat must not accumulate combinations

    combos = [name for name, combo in model.load_combos.items()
              if combo.combo_tags and 'buckling' in combo.combo_tags]
    assert combos == ['Buckling Mode 1', 'Buckling Mode 2', 'Buckling Mode 3']

    # Each mode is scaled to a largest component of 1 and the user's nodes carry the displacements
    np.testing.assert_allclose(np.max(np.abs(results.mode_shapes), axis=0), 1.0)
    assert 'Buckling Mode 1' in model.nodes['B'].DZ
    assert len(results.dof_map) == results.mode_shapes.shape[0] == results.free_dof_indices.size

    # The analysis ran on a copy with the subdivision nodes; the user's model is not that copy
    assert results.analysis_model is not model
    assert len(results.analysis_model.nodes) > len(model.nodes)
    assert results.elements_per_member == 8
    assert results.mode_count == 3 and not results.truncated
    assert model.solution == 'Buckling'


def test_static_analysis_clears_buckling_mode_combinations():
    """Like modal 'Mode n' combinations, 'Buckling Mode n' combinations are results, not loads."""
    model = _portal_frame()
    model.analyze_buckling(num_modes=2)
    assert 'Buckling Mode 1' in model.load_combos

    model.analyze()
    assert not any(name.startswith('Buckling Mode') for name in model.load_combos)
    assert model.solution != 'Buckling'


def test_invalid_arguments_are_rejected():
    """Bad arguments fail before any work is done."""
    model = _uniform_moment_beam()

    with pytest.raises(ValueError, match='plane'):
        model.analyze_buckling(plane='XQ')

    with pytest.raises(ValueError, match='out_of_plane'):
        model.analyze_buckling(out_of_plane=True)

    with pytest.raises(ValueError, match='warping'):
        model.analyze_buckling(warping='exact')

    with pytest.raises(ValueError, match='num_modes'):
        model.analyze_buckling(num_modes=0)

    with pytest.raises(ValueError, match='elements_per_member'):
        model.analyze_buckling(elements_per_member=0)


def test_member_geometric_stiffness_moment_terms():
    """The moment part of the member geometric stiffness matrix is symmetric and has the derived entries."""
    from Pynite.Member3D import Member3D

    L = 2.0
    kg = Member3D._kg_moments(L, Myi=1.0, Mzi=2.0, Myj=3.0, Mzj=4.0, Mxj=5.0)

    np.testing.assert_allclose(kg, kg.T)

    # Lateral displacement - twist coupling from the end moments
    assert kg[2, 3] == pytest.approx(2.0/L)      # w_i - θx_i : Mzi/L
    assert kg[2, 9] == pytest.approx(4.0/L)      # w_i - θx_j : Mzj/L
    assert kg[1, 3] == pytest.approx(1.0/L)      # v_i - θx_i : Myi/L
    # Twist - bending rotation coupling from the shear forces
    assert kg[3, 4] == pytest.approx((2.0 + 4.0)/6)
    assert kg[3, 5] == pytest.approx(-(1.0 + 3.0)/6)
    # Torque coupling between the bending planes
    assert kg[1, 4] == pytest.approx(5.0/L)
    assert kg[4, 11] == pytest.approx(5.0/2)

    # No moments, no terms
    assert not np.any(Member3D._kg_moments(L, 0.0, 0.0, 0.0, 0.0, 0.0))


# ---------------------------------------------------------------------------
# Entry point for direct execution
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    pytest.main([__file__, '-q'])
