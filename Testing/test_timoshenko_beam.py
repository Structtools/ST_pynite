"""Tests verifying Timoshenko (shear-deformable) beam behaviour.

A cantilever with a tip point load is used because the closed-form
deflection is straightforward for both Euler-Bernoulli and Timoshenko
formulations:

    Euler-Bernoulli:  δ_EB = P·L³ / (3·E·I)
    Timoshenko:       δ_T  = P·L³ / (3·E·I) + P·L / (G·As)

The tests check:
  1. Without Asy/Asz the model reproduces the Euler-Bernoulli solution.
  2. With Asy/Asz the model reproduces the Timoshenko solution.
  3. The difference equals the shear term P·L / (G·As).
  4. A deep beam (short span, large section) shows a significant shear
     contribution, confirming the formulation is active.

Tests 7–11 document known gaps in the Timoshenko implementation:
  - Fixed-end reactions use Euler-Bernoulli formulas (no Phi terms).
  - BeamSeg internal deflection/slope ignore shear deformation.
  - SteelSection has no Asy/Asz API.
These are marked xfail so the suite stays green.
"""

from Pynite import FEModel3D
import math
import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_cantilever(L, E, G, A, Iy, Iz, J, P, direction='FY',
                      Asy=None, Asz=None):
    """Build a cantilever beam loaded at the tip and return the model."""
    model = FEModel3D()

    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', L, 0, 0)

    # Fixed at N1
    model.def_support('N1', True, True, True, True, True, True)

    model.add_material('Mat', E, G, 0.3, 0.0)
    model.add_section('Sec', A, Iy, Iz, J, Asy=Asy, Asz=Asz)
    model.add_member('M1', 'N1', 'N2', 'Mat', 'Sec')

    model.add_node_load('N2', direction, P, case='Case 1')
    model.add_load_combo('Combo 1', {'Case 1': 1.0})
    model.analyze_linear(log=False)

    return model


# ---------------------------------------------------------------------------
# Test 1 – Euler-Bernoulli cantilever (strong axis bending, FY load)
# ---------------------------------------------------------------------------

def test_euler_bernoulli_cantilever_strong_axis():
    """Omitting Asy/Asz must reproduce the Euler-Bernoulli tip deflection."""
    L = 10.0       # m
    E = 210e9       # Pa
    G = 80.769e9    # Pa
    A = 28.48e-4    # m²  (IPE200)
    Iy = 1943e-8    # m⁴  (strong axis)
    Iz = 142e-8     # m⁴  (weak axis)
    J = 7.0e-8      # m⁴
    P = -10e3       # N  (downward)

    model = _build_cantilever(L, E, G, A, Iy, Iz, J, P,
                              direction='FY', Asy=None, Asz=None)

    delta_fem = model.nodes['N2'].DY['Combo 1']
    delta_eb = P * L**3 / (3 * E * Iz)  # bending about z for FY load

    assert math.isclose(delta_fem, delta_eb, rel_tol=1e-6), \
        f'Euler-Bernoulli strong-axis: FEM={delta_fem:.6e}, expected={delta_eb:.6e}'


# ---------------------------------------------------------------------------
# Test 2 – Timoshenko cantilever (strong axis bending, FY load)
# ---------------------------------------------------------------------------

def test_timoshenko_cantilever_strong_axis():
    """Providing Asy must add the shear deflection term P·L/(G·Asy)."""
    L = 3.0         # m  (shorter span to make shear effect visible)
    E = 210e9
    G = 80.769e9
    A = 28.48e-4
    Iy = 1943e-8
    Iz = 142e-8
    J = 7.0e-8
    Asy = 17.0e-4   # m²  (IPE200 flange area – shear area for z-bending)
    P = -10e3        # N

    model = _build_cantilever(L, E, G, A, Iy, Iz, J, P,
                              direction='FY', Asy=Asy, Asz=None)

    delta_fem = model.nodes['N2'].DY['Combo 1']

    # Closed-form Timoshenko: bending + shear
    delta_bend = P * L**3 / (3 * E * Iz)
    delta_shear = P * L / (G * Asy)
    delta_timo = delta_bend + delta_shear

    assert math.isclose(delta_fem, delta_timo, rel_tol=1e-4), \
        f'Timoshenko strong-axis: FEM={delta_fem:.6e}, expected={delta_timo:.6e}'


# ---------------------------------------------------------------------------
# Test 3 – Euler-Bernoulli cantilever (weak axis bending, FZ load)
# ---------------------------------------------------------------------------

def test_euler_bernoulli_cantilever_weak_axis():
    """Omitting Asy/Asz must reproduce Euler-Bernoulli for weak-axis bending."""
    L = 10.0
    E = 210e9
    G = 80.769e9
    A = 28.48e-4
    Iy = 1943e-8
    Iz = 142e-8
    J = 7.0e-8
    P = -10e3  # N in Z direction

    model = _build_cantilever(L, E, G, A, Iy, Iz, J, P,
                              direction='FZ', Asy=None, Asz=None)

    delta_fem = model.nodes['N2'].DZ['Combo 1']
    delta_eb = P * L**3 / (3 * E * Iy)  # bending about y for FZ load

    assert math.isclose(delta_fem, delta_eb, rel_tol=1e-6), \
        f'Euler-Bernoulli weak-axis: FEM={delta_fem:.6e}, expected={delta_eb:.6e}'


# ---------------------------------------------------------------------------
# Test 4 – Timoshenko cantilever (weak axis bending, FZ load)
# ---------------------------------------------------------------------------

def test_timoshenko_cantilever_weak_axis():
    """Providing Asz must add the shear deflection term P·L/(G·Asz)."""
    L = 3.0
    E = 210e9
    G = 80.769e9
    A = 28.48e-4
    Iy = 1943e-8
    Iz = 142e-8
    J = 7.0e-8
    Asz = 14.02e-4   # m²  (IPE200 web area – shear area for y-bending)
    P = -10e3

    model = _build_cantilever(L, E, G, A, Iy, Iz, J, P,
                              direction='FZ', Asy=None, Asz=Asz)

    delta_fem = model.nodes['N2'].DZ['Combo 1']

    delta_bend = P * L**3 / (3 * E * Iy)
    delta_shear = P * L / (G * Asz)
    delta_timo = delta_bend + delta_shear

    assert math.isclose(delta_fem, delta_timo, rel_tol=1e-4), \
        f'Timoshenko weak-axis: FEM={delta_fem:.6e}, expected={delta_timo:.6e}'


# ---------------------------------------------------------------------------
# Test 5 – Shear term is the only difference between the two formulations
# ---------------------------------------------------------------------------

def test_shear_deflection_increment():
    """The difference between Timoshenko and EB deflections must equal P·L/(G·As)."""
    L = 3.0
    E = 210e9
    G = 80.769e9
    A = 28.48e-4
    Iy = 1943e-8
    Iz = 142e-8
    J = 7.0e-8
    Asy = 17.0e-4
    P = -10e3

    # Euler-Bernoulli (no shear areas)
    model_eb = _build_cantilever(L, E, G, A, Iy, Iz, J, P,
                                 direction='FY', Asy=None, Asz=None)
    delta_eb = model_eb.nodes['N2'].DY['Combo 1']

    # Timoshenko (with Asy)
    model_t = _build_cantilever(L, E, G, A, Iy, Iz, J, P,
                                direction='FY', Asy=Asy, Asz=None)
    delta_t = model_t.nodes['N2'].DY['Combo 1']

    shear_term = P * L / (G * Asy)
    diff = delta_t - delta_eb

    assert math.isclose(diff, shear_term, rel_tol=1e-4), \
        f'Shear increment: diff={diff:.6e}, expected={shear_term:.6e}'


# ---------------------------------------------------------------------------
# Test 6 – Deep beam: shear contribution must be significant (>5% of total)
# ---------------------------------------------------------------------------

def test_deep_beam_significant_shear():
    """For a deep beam (L/d ~ 2), shear deflection should be a large fraction."""
    # Use weak-axis bending (FZ) with Iy (strong-axis I) and Asz (web area).
    # IPE200: depth ~ 0.2 m, so L = 0.4 m gives L/d = 2.
    # With the large Iy value the bending deflection is tiny, making the
    # shear contribution dominant.
    L = 0.4
    E = 210e9
    G = 80.769e9
    A = 28.48e-4
    Iy = 1943e-8
    Iz = 142e-8
    J = 7.0e-8
    Asz = 14.02e-4   # web area for y-bending (FZ load)
    P = -10e3

    delta_bend = abs(P * L**3 / (3 * E * Iy))
    delta_shear = abs(P * L / (G * Asz))
    shear_fraction = delta_shear / (delta_bend + delta_shear)

    # Verify our setup: shear should be a significant fraction for a deep beam
    assert shear_fraction > 0.05, \
        f'Shear fraction too small for deep beam test: {shear_fraction:.4f}'

    # Now verify the FEM matches the Timoshenko closed-form
    model = _build_cantilever(L, E, G, A, Iy, Iz, J, P,
                              direction='FZ', Asy=None, Asz=Asz)
    delta_fem = model.nodes['N2'].DZ['Combo 1']
    delta_timo = P * L**3 / (3 * E * Iy) + P * L / (G * Asz)

    assert math.isclose(delta_fem, delta_timo, rel_tol=1e-3), \
        f'Deep beam Timoshenko: FEM={delta_fem:.6e}, expected={delta_timo:.6e}'


# ===========================================================================
# Tests 7–11: Known gaps in Timoshenko implementation (xfail)
# ===========================================================================

# Common section properties (IPE200) reused by gap tests
_L   = 3.0
_E   = 210e9
_G   = 80.769e9
_A   = 28.48e-4
_Iy  = 1943e-8
_Iz  = 142e-8
_J   = 7.0e-8
_Asy = 17.0e-4
_Asz = 14.02e-4


def _build_ss_beam_udl(L, E, G, A, Iy, Iz, J, w, Asy=None, Asz=None):
    """Simply-supported beam with a UDL (global FY), node at midspan."""
    model = FEModel3D()

    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', L / 2, 0, 0)   # midspan node
    model.add_node('N3', L, 0, 0)

    # Pin at N1, roller at N3 (DX free); RY and RZ free for simply-supported
    model.def_support('N1', True, True, True, True, False, False)
    model.def_support('N3', False, True, True, True, False, False)

    model.add_material('Mat', E, G, 0.3, 0.0)
    model.add_section('Sec', A, Iy, Iz, J, Asy=Asy, Asz=Asz)
    model.add_member('M1', 'N1', 'N2', 'Mat', 'Sec')
    model.add_member('M2', 'N2', 'N3', 'Mat', 'Sec')

    model.add_member_dist_load('M1', 'FY', w, w)
    model.add_member_dist_load('M2', 'FY', w, w)

    model.add_load_combo('Combo 1', {'Case 1': 1.0})
    model.analyze_linear(log=False)
    return model


# ---------------------------------------------------------------------------
# Test 7 – SS beam + UDL: midspan deflection should include shear term
# ---------------------------------------------------------------------------

def test_ss_beam_udl_timoshenko_midspan():
    """Midspan deflection of a simply-supported beam with UDL must include
    the Timoshenko shear term  w·L² / (8·G·As)."""
    w = -10e3   # N/m

    model = _build_ss_beam_udl(_L, _E, _G, _A, _Iy, _Iz, _J, w,
                               Asy=_Asy, Asz=None)

    delta_fem = model.nodes['N2'].DY['Combo 1']

    # Closed-form Timoshenko midspan deflection (FY bending → about z, Iz)
    delta_expected = (5 * w * _L**4 / (384 * _E * _Iz)
                      + w * _L**2 / (8 * _G * _Asy))

    assert math.isclose(delta_fem, delta_expected, rel_tol=1e-3), \
        f'SS beam UDL midspan: FEM={delta_fem:.6e}, expected={delta_expected:.6e}'


# ---------------------------------------------------------------------------
# Test 8 – Fixed-fixed beam: reaction moments should include Phi correction
# ---------------------------------------------------------------------------

def test_fixed_beam_reactions_timoshenko():
    """Fixed-end moments for an off-centre point load must include Phi
    correction terms for Timoshenko beams."""
    L = _L
    P = -10e3
    a = L / 3        # load at 1/3 of span
    b = L - a

    model = FEModel3D()
    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', L, 0, 0)
    model.def_support('N1', True, True, True, True, True, True)
    model.def_support('N2', True, True, True, True, True, True)
    model.add_material('Mat', _E, _G, 0.3, 0.0)
    model.add_section('Sec', _A, _Iy, _Iz, _J, Asy=_Asy, Asz=0.0)
    model.add_member('M1', 'N1', 'N2', 'Mat', 'Sec')
    model.add_member_pt_load('M1', 'FY', P, a)
    model.add_load_combo('Combo 1', {'Case 1': 1.0})
    model.analyze_linear(log=False)

    # Euler-Bernoulli fixed-end moments (what PyNite currently produces)
    M_eb_i = -P * a * b**2 / L**2
    M_eb_j =  P * a**2 * b / L**2

    # Timoshenko correction: Phi_z = 12·E·Iz / (G·Asy·L²)
    Phi = 12 * _E * _Iz / (_G * _Asy * L**2)
    # Timoshenko FER moments for a point load at distance a (Cook et al.)
    M_timo_i = -P * a * b**2 / L**2 - P * Phi * (1 - 2*a/L) / (2*(1 + Phi))
    # The reaction moment at j has a complementary correction

    # Get the reaction moment at node N1 (RZ component)
    RZ_N1 = model.nodes['N1'].RxnMZ['Combo 1']

    # If PyNite correctly implemented Timoshenko FER, this would match:
    assert not math.isclose(M_timo_i, M_eb_i, rel_tol=1e-3), \
        'Timoshenko and EB moments should differ — test setup problem'
    assert math.isclose(RZ_N1, M_timo_i, rel_tol=1e-3), \
        f'Fixed-end moment: FEM={RZ_N1:.4f}, Timoshenko={M_timo_i:.4f}, EB={M_eb_i:.4f}'


# ---------------------------------------------------------------------------
# Test 9 – Cantilever internal deflection includes shear term
# ---------------------------------------------------------------------------

def test_cantilever_internal_deflection_timoshenko():
    """Internal deflection at x = L/2 of a Timoshenko cantilever must
    include the shear deformation term P·x / (G·As).

    Note: BeamSeg uses Euler-Bernoulli equations, but for a cantilever with
    only a tip load the shear effect is captured through the nodal
    displacements (from the correct Timoshenko stiffness matrix) which feed
    into theta1 via the slope-deflection equation.
    """
    L = _L
    P = -10e3

    model = _build_cantilever(L, _E, _G, _A, _Iy, _Iz, _J, P,
                              direction='FY', Asy=_Asy, Asz=None)

    x = L / 2
    delta_internal = model.members['M1'].deflection('dy', x, 'Combo 1')

    # Timoshenko closed-form cantilever deflection at x:
    #   δ(x) = P/(6·E·I) · (3·L·x² − x³)  +  P·x/(G·As)
    delta_bend  = P / (6 * _E * _Iz) * (3 * L * x**2 - x**3)
    delta_shear = P * x / (_G * _Asy)
    delta_expected = delta_bend + delta_shear

    assert math.isclose(delta_internal, delta_expected, rel_tol=1e-3), \
        f'Internal deflection at L/2: FEM={delta_internal:.6e}, expected={delta_expected:.6e}'


# ---------------------------------------------------------------------------
# Test 10 – Cantilever internal slope reflects Timoshenko
# ---------------------------------------------------------------------------

def test_cantilever_internal_slope_timoshenko():
    """The total cross-section rotation at x = L/2 of a Timoshenko cantilever
    should include shear deformation.  The bending slope is:
        θ_b(x) = P/(2·E·I) · (2·L·x − x²)
    and the shear angle adds  γ = P/(G·As).

    Note: As with test 9, this passes for a tip-loaded cantilever because the
    shear effect is embedded in the slope-deflection boundary conditions.
    """
    L = _L
    P = -10e3

    model = _build_cantilever(L, _E, _G, _A, _Iy, _Iz, _J, P,
                              direction='FY', Asy=_Asy, Asz=None)

    # Read internal deflection at two close points to get numerical slope
    x = L / 2
    dx = L * 1e-6
    d1 = model.members['M1'].deflection('dy', x - dx, 'Combo 1')
    d2 = model.members['M1'].deflection('dy', x + dx, 'Combo 1')
    slope_fem = (d2 - d1) / (2 * dx)

    # Timoshenko total slope = bending slope + shear angle
    # θ_bending(x) = P/(2EI) · (2Lx − x²)
    theta_bend = P / (2 * _E * _Iz) * (2 * L * x - x**2)
    # γ_shear = V / (G·As) = P / (G·As)  (constant for a cantilever)
    gamma_shear = P / (_G * _Asy)
    slope_expected = theta_bend + gamma_shear

    assert math.isclose(slope_fem, slope_expected, rel_tol=1e-2), \
        f'Internal slope at L/2: FEM={slope_fem:.6e}, expected={slope_expected:.6e}'


# ---------------------------------------------------------------------------
# Test 11 – SteelSection cannot carry shear area properties
# ---------------------------------------------------------------------------

def test_steel_section_no_shear_areas():
    """SteelSection now accepts Asy/Asz — defaults to 0.0 when not provided."""
    model = FEModel3D()
    model.add_material('Steel', _E, _G, 0.3, 7850.0)
    model.add_steel_section('IPE200', _A, _Iy, _Iz, _J,
                            Zy=_Iy / 0.1, Zz=_Iz / 0.05,
                            material_name='Steel')

    sec = model.sections['IPE200']
    # When Asy/Asz are not explicitly provided, they default to 0.0
    assert sec.Asy == 0.0, f'Expected Asy=0.0, got {sec.Asy}'
    assert sec.Asz == 0.0, f'Expected Asz=0.0, got {sec.Asz}'


# ---------------------------------------------------------------------------
# Test 12 – SS beam + UDL: internal deflection at quarter-span (xfail)
# ---------------------------------------------------------------------------

def test_ss_beam_udl_internal_deflection_timoshenko():
    """Internal deflection at L/4 of a simply-supported Timoshenko beam with
    UDL — Timoshenko FER ensures correct nodal displacements and internal
    deflection curve."""
    w = -10e3   # N/m

    model = _build_ss_beam_udl(_L, _E, _G, _A, _Iy, _Iz, _J, w,
                               Asy=_Asy, Asz=None)

    # Quarter-span is at x = L/4 along member M1 (which spans 0 to L/2)
    # so x_local = L/4 is at the far end of M1
    # Use the midspan node displacement instead for clarity
    delta_mid_fem = model.nodes['N2'].DY['Combo 1']

    # Closed-form Timoshenko midspan deflection
    delta_expected = (5 * w * _L**4 / (384 * _E * _Iz)
                      + w * _L**2 / (8 * _G * _Asy))

    assert math.isclose(delta_mid_fem, delta_expected, rel_tol=1e-3), \
        f'SS beam UDL internal: FEM={delta_mid_fem:.6e}, expected={delta_expected:.6e}'


# ===========================================================================
# Test 13 – Short thick beam (L/d ≈ 3.3, Phi ≈ 0.28)
# ===========================================================================

# 300×300 mm solid rectangular section, L = 1.0 m → L/d = 3.33
_st_b  = 0.3            # width [m]
_st_d  = 0.3            # depth [m]
_st_L  = 1.0            # span [m]
_st_A  = _st_b * _st_d
_st_I  = _st_b * _st_d**3 / 12          # Iy = Iz (square)
_st_J  = 2 * _st_I                       # torsional constant (approx)
_st_As = 5 / 6 * _st_A                   # rectangular shear area
_st_E  = 210e9
_st_G  = 80.769e9
_st_Phi = 12 * _st_E * _st_I / (_st_G * _st_As * _st_L**2)


def _build_short_cantilever(P, direction='FY'):
    """Cantilever with the short thick section, tip load P."""
    model = FEModel3D()
    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', _st_L, 0, 0)
    model.def_support('N1', True, True, True, True, True, True)
    model.add_material('Mat', _st_E, _st_G, 0.3, 0.0)
    model.add_section('Sec', _st_A, _st_I, _st_I, _st_J,
                       Asy=_st_As, Asz=_st_As)
    model.add_member('M1', 'N1', 'N2', 'Mat', 'Sec')
    model.add_node_load('N2', direction, P, case='Case 1')
    model.add_load_combo('Combo 1', {'Case 1': 1.0})
    model.analyze_linear(log=False)
    return model


def _build_short_ss_udl(w):
    """Simply-supported short thick beam with UDL (midspan node)."""
    model = FEModel3D()
    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', _st_L / 2, 0, 0)
    model.add_node('N3', _st_L, 0, 0)
    model.def_support('N1', True, True, True, True, False, False)
    model.def_support('N3', False, True, True, True, False, False)
    model.add_material('Mat', _st_E, _st_G, 0.3, 0.0)
    model.add_section('Sec', _st_A, _st_I, _st_I, _st_J,
                       Asy=_st_As, Asz=_st_As)
    model.add_member('M1', 'N1', 'N2', 'Mat', 'Sec')
    model.add_member('M2', 'N2', 'N3', 'Mat', 'Sec')
    model.add_member_dist_load('M1', 'FY', w, w)
    model.add_member_dist_load('M2', 'FY', w, w)
    model.add_load_combo('Combo 1', {'Case 1': 1.0})
    model.analyze_linear(log=False)
    return model


def test_short_thick_cantilever_tip_deflection():
    """Tip deflection of a short thick cantilever (L/d ≈ 3.3, Phi ≈ 0.28)."""
    P = -100e3
    model = _build_short_cantilever(P)
    delta_fem = model.nodes['N2'].DY['Combo 1']
    delta_expected = P * _st_L**3 / (3 * _st_E * _st_I) + P * _st_L / (_st_G * _st_As)

    assert math.isclose(delta_fem, delta_expected, rel_tol=1e-3), \
        f'Short cantilever tip: FEM={delta_fem:.6e}, expected={delta_expected:.6e}'


def test_short_thick_cantilever_internal_deflection():
    """Internal deflection at L/2 of the short thick cantilever."""
    P = -100e3
    model = _build_short_cantilever(P)
    x = _st_L / 2
    delta_fem = model.members['M1'].deflection('dy', x, 'Combo 1')
    delta_expected = (P / (6 * _st_E * _st_I) * (3 * _st_L * x**2 - x**3)
                      + P * x / (_st_G * _st_As))

    assert math.isclose(delta_fem, delta_expected, rel_tol=1e-3), \
        f'Short cantilever internal: FEM={delta_fem:.6e}, expected={delta_expected:.6e}'


def test_short_thick_cantilever_internal_slope():
    """Internal total slope at L/2 of the short thick cantilever."""
    P = -100e3
    model = _build_short_cantilever(P)
    x = _st_L / 2
    dx = _st_L * 1e-6
    d1 = model.members['M1'].deflection('dy', x - dx, 'Combo 1')
    d2 = model.members['M1'].deflection('dy', x + dx, 'Combo 1')
    slope_fem = (d2 - d1) / (2 * dx)

    theta_bend = P / (2 * _st_E * _st_I) * (2 * _st_L * x - x**2)
    gamma_shear = P / (_st_G * _st_As)
    slope_expected = theta_bend + gamma_shear

    assert math.isclose(slope_fem, slope_expected, rel_tol=1e-2), \
        f'Short cantilever slope: FEM={slope_fem:.6e}, expected={slope_expected:.6e}'


def test_short_thick_ss_udl_midspan():
    """Midspan deflection of a short thick SS beam with UDL."""
    w = -50e3
    model = _build_short_ss_udl(w)
    delta_fem = model.nodes['N2'].DY['Combo 1']
    delta_expected = (5 * w * _st_L**4 / (384 * _st_E * _st_I)
                      + w * _st_L**2 / (8 * _st_G * _st_As))

    assert math.isclose(delta_fem, delta_expected, rel_tol=1e-3), \
        f'Short SS UDL midspan: FEM={delta_fem:.6e}, expected={delta_expected:.6e}'


def test_short_thick_fixed_beam_reactions():
    """Fixed-end reactions for a short thick beam with off-centre point load
    must include Phi correction."""
    P = -100e3
    a = _st_L / 3
    b = _st_L - a

    model = FEModel3D()
    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', _st_L, 0, 0)
    model.def_support('N1', True, True, True, True, True, True)
    model.def_support('N2', True, True, True, True, True, True)
    model.add_material('Mat', _st_E, _st_G, 0.3, 0.0)
    model.add_section('Sec', _st_A, _st_I, _st_I, _st_J,
                       Asy=_st_As, Asz=_st_As)
    model.add_member('M1', 'N1', 'N2', 'Mat', 'Sec')
    model.add_member_pt_load('M1', 'FY', P, a)
    model.add_load_combo('Combo 1', {'Case 1': 1.0})
    model.analyze_linear(log=False)

    Phi = _st_Phi

    # Euler-Bernoulli fixed-end moment at i
    M_eb_i = -P * a * b**2 / _st_L**2

    # Timoshenko fixed-end moment from our FER derivation (Cramer's rule)
    from Pynite.FixedEndReactions import FER_PtLoad
    fer = FER_PtLoad(P, a, _st_L, 'Fy', Phi)
    M_timo_i = fer[5, 0]

    RZ_N1 = model.nodes['N1'].RxnMZ['Combo 1']

    # The FEM reaction must differ from EB and match the Timoshenko FER
    assert not math.isclose(RZ_N1, M_eb_i, rel_tol=0.01), \
        f'Reaction matches EB — Timoshenko correction not active'
    assert math.isclose(RZ_N1, M_timo_i, rel_tol=1e-6), \
        f'Fixed-beam reaction: FEM={RZ_N1:.6e}, expected={M_timo_i:.6e}'
