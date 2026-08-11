"""Timoshenko fixed-end reactions, checked against an independent oracle.

Why this module exists
----------------------
The Timoshenko fixed-end reaction (FER) code carried a sign error on every
Phi term for a long time, producing a 15% error in the support moment of a
fixed-fixed beam with a point load at L/3, and a 125% error (with a sign
flip) for an applied couple. The existing test suite stayed green
throughout. Two properties conspired to hide it:

1. **Symmetric loads cannot detect it.** A midspan point load gives PL/8 and
   a full-span UDL gives wL^2/12 *for any value of Phi*. Every Timoshenko FER
   test in the suite used a symmetric load, so all of them passed against
   both a correct and an incorrect implementation. Asymmetric load cases are
   mandatory, and this module is built around them.

2. **Global equilibrium is satisfied either way.** The j-end moment is
   back-computed from statics, so sum(F) and sum(M) are exactly zero even
   when the load is distributed wrongly between the two ends.
   `analyze_linear(check_statics=True)` cannot see the error.

The oracle
----------
The exact 2-node Timoshenko element is *exact* when every load acts at a
node, because no shape-function interpolation of the load is involved. So a
refined mesh with a node placed at the load position gives the true answer
while never calling FixedEndReactions at all. For point loads and couples
this is exact even at one element per span; the mesh refinement below is
belt-and-braces. Distributed loads use tributary lumping, which converges
rather than being exact, hence the looser tolerance on those cases.

This makes the oracle genuinely independent. A test that computes its
expected value by calling FER_PtLoad would only assert that the assembler
and the FER function agree with each other, which is true for any value of
the Phi terms and therefore cannot fail.
"""

import math

import pytest

from Pynite import FEModel3D


# ---------------------------------------------------------------------------
# Section: 300x300 mm rectangle. Phi = 12*E*I/(G*As*L^2) = 0.27 / L^2.
# ---------------------------------------------------------------------------

_E = 30e6
_G = 12e6
_NU = 0.25
_B = _H = 0.3
_A = _B * _H
_I = _B * _H**3 / 12
_J = 1e-3
_AS = 5.0 / 6.0 * _A


def _phi(L):
    """Timoshenko shear parameter for the reference section at span L."""
    return 12 * _E * _I / (_G * _AS * L**2)


# Spans chosen to sweep Phi across the interesting range. The Phi = 1.0 case
# is deliberate: the previous implementation had a spurious pole there
# (det = L^4*(Phi - 1)/12) and raised ZeroDivisionError. A correct
# implementation is perfectly smooth through it.
_L_PHI_019 = 1.2                        # Phi = 0.1875
_L_PHI_027 = 1.0                        # Phi = 0.27
_L_PHI_100 = math.sqrt(0.27)            # Phi = 1.0 exactly
_L_PHI_133 = 0.45                       # Phi = 1.3333

_SPANS = [_L_PHI_019, _L_PHI_027, _L_PHI_133]


# ---------------------------------------------------------------------------
# Model builders
# ---------------------------------------------------------------------------

def _fixed_fixed(n, L, shear=True, extra_x=()):
    """Fixed-fixed beam meshed into n elements, with optional extra nodes."""
    model = FEModel3D()
    model.add_material('Mat', _E, _G, _NU, 0.0)
    As = _AS if shear else 0.0
    model.add_section('Sec', _A, _I, _I, _J, Asy=As, Asz=As)

    xs = sorted(set([i * L / n for i in range(n + 1)] + list(extra_x)))
    for i, x in enumerate(xs):
        model.add_node(f'N{i}', x, 0, 0)
    model.def_support('N0', *(True,) * 6)
    model.def_support(f'N{len(xs) - 1}', *(True,) * 6)
    for i in range(len(xs) - 1):
        model.add_member(f'M{i}', f'N{i}', f'N{i + 1}', 'Mat', 'Sec')
    return model, xs


def _solve(model, plane='y'):
    """Analyse and return the i-end (force, moment) reaction pair."""
    model.add_load_combo('Combo 1', {'Case 1': 1.0})
    model.analyze_linear(log=False, check_statics=False)
    n0 = model.nodes['N0']
    if plane == 'y':
        return n0.RxnFY['Combo 1'], n0.RxnMZ['Combo 1']
    return n0.RxnFZ['Combo 1'], n0.RxnMY['Combo 1']


def _fer_result(L, kind, value, a, shear=True, plane='y'):
    """Single element, so the load goes through the FER code path."""
    model, _ = _fixed_fixed(1, L, shear=shear)
    direction = {'point': 'Fy', 'couple': 'Mz'}[kind] if plane == 'y' \
        else {'point': 'Fz', 'couple': 'My'}[kind]
    model.add_member_pt_load('M0', direction, value, a)
    return _solve(model, plane)


def _oracle_result(L, kind, value, a, shear=True, plane='y', n=4):
    """Refined mesh with the load applied AT A NODE. Never touches FER."""
    model, xs = _fixed_fixed(n, L, shear=shear, extra_x=(a,))
    direction = {'point': 'FY', 'couple': 'MZ'}[kind] if plane == 'y' \
        else {'point': 'FZ', 'couple': 'MY'}[kind]
    model.add_node_load(f'N{xs.index(a)}', direction, value)
    return _solve(model, plane)


def _oracle_distributed(L, w1, w2, x1, x2, n=2000):
    """Tributary-lumped nodal oracle for a linearly varying load.

    Converges rather than being exact, so callers use a looser tolerance.
    """
    model, xs = _fixed_fixed(n, L)
    for i, x in enumerate(xs):
        lo = max(x1, (xs[i - 1] + x) / 2 if i else x)
        hi = min(x2, (x + xs[i + 1]) / 2 if i < len(xs) - 1 else x)
        if hi > lo:
            w_mid = w1 + (w2 - w1) * ((lo + hi) / 2 - x1) / (x2 - x1)
            model.add_node_load(f'N{i}', 'FY', w_mid * (hi - lo))
    return _solve(model)


# ---------------------------------------------------------------------------
# Closed forms, typed independently of the implementation.
# Reaction convention: RxnFY = -R_i, RxnMZ = M_0.
# ---------------------------------------------------------------------------

def _closed_form_point(P, a, L, Phi):
    b = L - a
    fy = -P * b * (Phi * L**2 + b * (L + 2 * a)) / (L**3 * (1 + Phi))
    mz = -P * a * b * (b + Phi * L / 2) / (L**2 * (1 + Phi))
    return fy, mz


def _closed_form_couple(M, a, L, Phi):
    b = L - a
    fy = 6 * M * a * b / (L**3 * (1 + Phi))
    mz = M * b * ((2 * a - b) - Phi * L) / (L**2 * (1 + Phi))
    return fy, mz


def _assert_close(got, expected, rel_tol, label):
    for g, e, name in zip(got, expected, ('force', 'moment')):
        assert math.isclose(g, e, rel_tol=rel_tol, abs_tol=1e-9), \
            f'{label} [{name}]: got {g!r}, expected {e!r}'


# ---------------------------------------------------------------------------
# Case 1 - asymmetric point loads. This is the case the sign error broke.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('L', _SPANS)
@pytest.mark.parametrize('frac', [0.25, 1 / 3, 0.75])
def test_point_load_asymmetric(L, frac):
    a = frac * L
    got = _fer_result(L, 'point', -100.0, a)
    expected = _oracle_result(L, 'point', -100.0, a)
    _assert_close(got, expected, 1e-9, f'PtLoad L={L} a/L={frac} Phi={_phi(L):.4f}')


# ---------------------------------------------------------------------------
# Case 2 - midspan point load. Guards the symmetric case, and pins the fact
# that it is Phi-independent (which is exactly why it cannot detect bugs).
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('L', _SPANS)
def test_point_load_midspan_is_phi_independent(L):
    P = -100.0
    got = _fer_result(L, 'point', P, L / 2)
    assert math.isclose(got[0], -P / 2, rel_tol=1e-12)
    assert math.isclose(got[1], -P * L / 8, rel_tol=1e-12)

    no_shear = _fer_result(L, 'point', P, L / 2, shear=False)
    _assert_close(got, no_shear, 1e-12, 'midspan must not depend on Phi')


# ---------------------------------------------------------------------------
# Case 3 - concentrated couples. The old code was 125% out here, with a
# sign flip, because it wrongly gave the couple a shear-impulse term.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('L', _SPANS)
@pytest.mark.parametrize('frac', [0.25, 0.5, 0.8])
def test_couple_asymmetric(L, frac):
    a = frac * L
    got = _fer_result(L, 'couple', 50.0, a)
    expected = _oracle_result(L, 'couple', 50.0, a)
    _assert_close(got, expected, 1e-9, f'Couple L={L} a/L={frac} Phi={_phi(L):.4f}')


# ---------------------------------------------------------------------------
# Cases 4 and 5 - partial and trapezoidal distributed loads.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('w1,w2,label', [(-30.0, -30.0, 'partial UDL'),
                                         (-10.0, -40.0, 'trapezoid'),
                                         (-40.0, -10.0, 'trapezoid reversed')])
def test_partial_distributed_load(w1, w2, label):
    L, x1, x2 = 1.0, 0.2, 0.8
    model, _ = _fixed_fixed(1, L)
    model.add_member_dist_load('M0', 'Fy', w1, w2, x1, x2)
    got = _solve(model)
    expected = _oracle_distributed(L, w1, w2, x1, x2)
    _assert_close(got, expected, 2e-4, f'{label} Phi={_phi(L):.4f}')


# ---------------------------------------------------------------------------
# Case 6 - full-span UDL is exactly Phi-independent.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('L', _SPANS)
def test_full_span_udl_is_phi_independent(L):
    w = -30.0
    model, _ = _fixed_fixed(1, L)
    model.add_member_dist_load('M0', 'Fy', w, w, 0.0, L)
    got = _solve(model)
    assert math.isclose(got[0], -w * L / 2, rel_tol=1e-12)
    assert math.isclose(got[1], -w * L**2 / 12, rel_tol=1e-12)


# ---------------------------------------------------------------------------
# Case 7 - Phi = 1.0 exactly. The old determinant L^4*(Phi-1)/12 vanished
# here and raised ZeroDivisionError. Reachable with round-number inputs:
# L/h ~ 1.77 for a rectangular section with nu = 0.3.
# ---------------------------------------------------------------------------

def test_phi_exactly_one_is_not_singular():
    L = _L_PHI_100
    assert math.isclose(_phi(L), 1.0, rel_tol=1e-12), 'span should give Phi = 1'

    a = L / 3
    got = _fer_result(L, 'point', -100.0, a)
    expected = _oracle_result(L, 'point', -100.0, a)
    _assert_close(got, expected, 1e-9, 'PtLoad at Phi = 1.0')

    for value in (0.999, 1.0, 1.001):
        assert all(math.isfinite(v) for v in
                   _closed_form_point(-100.0, a, L, value))


@pytest.mark.parametrize('Phi', [0.0, 0.5, 1.0, 5.0, 50.0, 500.0])
def test_fer_smooth_across_phi_sweep(Phi):
    """FER_PtLoad must stay finite and match the closed form for any Phi >= 0."""
    from Pynite.FixedEndReactions import FER_PtLoad

    L, P, a = 1.0, -100.0, 1.0 / 3
    fer = FER_PtLoad(P, a, L, 'Fy', Phi)
    assert all(math.isfinite(v) for v in fer.flatten())

    fy, mz = _closed_form_point(P, a, L, Phi)
    # FER[1] = -R_i = RxnFY equivalent, FER[5] = M_0
    assert math.isclose(fer[1, 0], fy, rel_tol=1e-12)
    assert math.isclose(fer[5, 0], mz, rel_tol=1e-12)


# ---------------------------------------------------------------------------
# Case 8 - Phi = 0 regression. Frozen Euler-Bernoulli values, so a reviewer
# can confirm nothing changed for users who supply no shear areas.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('frac', [0.25, 1 / 3, 0.5, 0.75])
def test_phi_zero_reduces_to_euler_bernoulli_point_load(frac):
    L, P = 1.0, -100.0
    a, b = frac * L, L - frac * L
    got = _fer_result(L, 'point', P, a, shear=False)
    expected = (-P * b**2 * (L + 2 * a) / L**3, -P * a * b**2 / L**2)
    _assert_close(got, expected, 1e-13, f'EB PtLoad a/L={frac}')


@pytest.mark.parametrize('frac', [0.25, 0.5, 0.8])
def test_phi_zero_reduces_to_euler_bernoulli_couple(frac):
    L, M = 1.0, 50.0
    a, b = frac * L, L - frac * L
    got = _fer_result(L, 'couple', M, a, shear=False)
    expected = (6 * M * a * b / L**3, M * b * (2 * a - b) / L**2)
    _assert_close(got, expected, 1e-13, f'EB Couple a/L={frac}')


def test_phi_zero_reduces_to_euler_bernoulli_udl():
    L, w = 1.0, -30.0
    model, _ = _fixed_fixed(1, L, shear=False)
    model.add_member_dist_load('M0', 'Fy', w, w, 0.0, L)
    got = _solve(model)
    _assert_close(got, (-w * L / 2, -w * L**2 / 12), 1e-13, 'EB full UDL')


# ---------------------------------------------------------------------------
# Case 9 - Fz / My direction mirrors, exercising the other two branches of
# the 12-DOF sign mapping in _solve_FER_bending.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('frac', [0.25, 1 / 3, 0.75])
def test_point_load_fz_direction(frac):
    L = 1.0
    a = frac * L
    got = _fer_result(L, 'point', -100.0, a, plane='z')
    expected = _oracle_result(L, 'point', -100.0, a, plane='z')
    _assert_close(got, expected, 1e-9, f'PtLoad Fz a/L={frac}')


@pytest.mark.parametrize('frac', [0.25, 0.8])
def test_couple_my_direction(frac):
    L = 1.0
    a = frac * L
    got = _fer_result(L, 'couple', 50.0, a, plane='z')
    expected = _oracle_result(L, 'couple', 50.0, a, plane='z')
    _assert_close(got, expected, 1e-9, f'Couple My a/L={frac}')


# ---------------------------------------------------------------------------
# Unit level - catches transcription errors in the 12-DOF mapping that a
# model-level test could mask.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('Phi', [0.0, 0.27, 1.0, 1.3333])
@pytest.mark.parametrize('frac', [0.25, 1 / 3, 0.5, 0.75])
def test_fer_ptload_matches_closed_form(Phi, frac):
    from Pynite.FixedEndReactions import FER_PtLoad

    L, P = 1.0, -100.0
    a = frac * L
    fer = FER_PtLoad(P, a, L, 'Fy', Phi)
    fy, mz = _closed_form_point(P, a, L, Phi)

    assert math.isclose(fer[1, 0], fy, rel_tol=1e-12)
    assert math.isclose(fer[5, 0], mz, rel_tol=1e-12)
    # Equilibrium: the two end shears must carry the whole load.
    assert math.isclose(fer[1, 0] + fer[7, 0], -P, rel_tol=1e-12)


@pytest.mark.parametrize('Phi', [0.0, 0.27, 1.0, 1.3333])
@pytest.mark.parametrize('frac', [0.25, 0.5, 0.8])
def test_fer_moment_matches_closed_form(Phi, frac):
    from Pynite.FixedEndReactions import FER_Moment

    L, M = 1.0, 50.0
    a = frac * L
    fer = FER_Moment(M, a, L, 'Mz', Phi)
    fy, mz = _closed_form_couple(M, a, L, Phi)

    assert math.isclose(fer[1, 0], fy, rel_tol=1e-12)
    assert math.isclose(fer[5, 0], mz, rel_tol=1e-12)
    # A couple applies no net transverse load.
    assert math.isclose(fer[1, 0] + fer[7, 0], 0.0, abs_tol=1e-9)
