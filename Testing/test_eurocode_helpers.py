"""
Tests for the EN 1993-1-1 helpers in `Pynite.EurocodeHelpers`.

`get_critical_moment` is checked against the closed forms it reduces to (uniform moment with fork
supports, with and without warping restraint), against Pynite's own lateral-torsional buckling
analysis with a load applied off the shear centre, and for independence of the axis naming of
the section. `get_buckling_length` is checked to pass its arguments through.
"""

import warnings

import numpy as np
import pytest
from Pynite import FEModel3D
from Pynite.EurocodeHelpers import get_buckling_length, get_critical_moment

# IPE300 (N, m), Eurocode axis naming: Iy strong, Iz weak
E, G = 210e9, 81e9
IPE300 = dict(A=53.8e-4, Iy=8356e-8, Iz=603.8e-8, J=20.12e-8, Iw=125.9e-9)


def _beam(L=6.0, swap_axes=False, **section):
    """A beam along X with the IPE300 section, or with its Iy and Iz swapped (Pynite naming)."""
    props = dict(IPE300, **section)
    if swap_axes:
        props['Iy'], props['Iz'] = props['Iz'], props['Iy']
    model = FEModel3D()
    model.add_material('Steel', E, G, 0.3, 7850)
    model.add_section('Sec', **props)
    model.add_node('N0', 0.0, 0.0, 0.0)
    model.add_node('N1', L, 0.0, 0.0)
    model.add_member('Beam', 'N0', 'N1', 'Steel', 'Sec')
    return model


def _uniform_moment_mcr(L, k=1.0, kw=1.0):
    """M_cr = π²EIz/(kL)²·sqrt((k/kw)²·Iw/Iz + (kL)²GJ/(π²EIz)) of the IPE300."""
    Iz, Iw, J = IPE300['Iz'], IPE300['Iw'], IPE300['J']
    return np.pi**2*E*Iz/(k*L)**2*np.sqrt((k/kw)**2*Iw/Iz + (k*L)**2*G*J/(np.pi**2*E*Iz))


def test_critical_moment_uniform_moment_matches_closed_form():
    """C1 = 1, k = kw = 1 is the fork-supported beam in uniform moment, (π/L)·sqrt(EIz·(GJ + π²EIw/L²))."""
    L = 6.0
    model = _beam(L)
    Iz, Iw, J = IPE300['Iz'], IPE300['Iw'], IPE300['J']
    expected = np.pi/L*np.sqrt(E*Iz*(G*J + np.pi**2*E*Iw/L**2))
    assert get_critical_moment(model, 'Beam') == pytest.approx(expected, rel=1e-12)


def test_critical_moment_uses_warping_length_factor_and_given_length():
    """`kw` restrains warping and `L` replaces the member length as the length between restraints."""
    model = _beam(6.0)
    assert get_critical_moment(model, 'Beam', kw=0.5) == pytest.approx(_uniform_moment_mcr(6.0, kw=0.5), rel=1e-12)
    assert get_critical_moment(model, 'Beam', k=0.5, kw=0.5) == pytest.approx(_uniform_moment_mcr(6.0, k=0.5, kw=0.5), rel=1e-12)
    assert get_critical_moment(model, 'Beam', L=3.0) == pytest.approx(_uniform_moment_mcr(3.0), rel=1e-12)
    assert get_critical_moment(model, 'Beam', kw=0.5) > get_critical_moment(model, 'Beam')


def test_critical_moment_does_not_depend_on_axis_naming():
    """The minor axis is the one with the smaller second moment of area, whichever name it has."""
    eurocode = get_critical_moment(_beam(), 'Beam', C1=1.13, C2=0.45, z_g=0.15)
    pynite = get_critical_moment(_beam(swap_axes=True), 'Beam', C1=1.13, C2=0.45, z_g=0.15)
    assert pynite == pytest.approx(eurocode, rel=1e-12)


@pytest.mark.parametrize('z_g', [0.1, -0.1])
def test_critical_moment_load_height_matches_buckling_analysis(z_g):
    """
    A fork-supported beam with a central point load on the top (z_g > 0) or bottom (z_g < 0)
    flange, Iw = 0 and E·Iz = G·J = L = 1: the three-factor formula with C1 = 1.365 and
    C2 = 0.553 against Pynite's buckling analysis, to the 5 % accuracy of the formula.
    """
    model = FEModel3D()
    model.add_material('Mat', 1.0, 1.0, 0.0, 0.0)
    model.add_section('Sec', 1e3, 1.0, 10.0, 1.0, Iw=0.0)
    model.add_node('N0', 0.0, 0.0, 0.0)
    model.add_node('N1', 1.0, 0.0, 0.0)
    model.add_member('Beam', 'N0', 'N1', 'Mat', 'Sec')
    model.def_support('N0', True, True, True, True, False, False)
    model.def_support('N1', False, True, True, True, False, False)
    model.add_member_pt_load('Beam', 'FY', -1.0, 0.5, load_height=z_g)

    results = model.analyze_buckling(num_modes=1, elements_per_member=16)
    M_cr_fe = results.load_multipliers[0]/4   # M = P·L/4

    M_cr = get_critical_moment(model, 'Beam', C1=1.365, C2=0.553, z_g=z_g)
    assert M_cr == pytest.approx(M_cr_fe, rel=0.05)

    # The top flange lowers M_cr, the bottom flange raises it
    M_cr_0 = get_critical_moment(model, 'Beam', C1=1.365)
    assert (M_cr < M_cr_0) == (z_g > 0)


def test_critical_moment_warns_about_c2_without_load_height():
    """The load position term is C2·z_g, so C2 alone does nothing and is most likely a mistake."""
    model = _beam()
    with pytest.warns(UserWarning, match='z_g'):
        M_cr = get_critical_moment(model, 'Beam', C2=0.45)
    assert M_cr == pytest.approx(get_critical_moment(model, 'Beam'))

    with warnings.catch_warnings():
        warnings.simplefilter('error')
        get_critical_moment(model, 'Beam', C2=0.45, z_g=0.15)


def test_critical_moment_rejects_invalid_input():
    """A missing warping constant and non-positive lengths are errors."""
    model = _beam(Iw=None)
    with pytest.raises(AttributeError):
        get_critical_moment(model, 'Beam')

    model = _beam()
    for kwargs in ({'k': 0.0}, {'kw': -1.0}, {'L': 0.0}):
        with pytest.raises(ValueError):
            get_critical_moment(model, 'Beam', **kwargs)


def test_get_buckling_length_passes_arguments_through():
    """`get_buckling_length` is `BucklingResults.effective_length` by member name, including `x`."""
    model = FEModel3D()
    model.add_material('Mat', 1.0, 1.0, 0.0, 0.0)
    model.add_section('Sec', 1e3, 1.0, 1.0, 2.0)
    model.add_node('Base', 0.0, 0.0, 0.0)
    model.add_node('Top', 0.0, 1.0, 0.0)
    model.add_member('Column', 'Base', 'Top', 'Mat', 'Sec')
    model.def_support('Base', True, True, True, False, True, False)
    model.def_support('Top', True, False, True, False, True, False)
    model.add_node_load('Top', 'FY', -1.0)
    model.add_member_dist_load('Column', 'FY', -2.0, -2.0)

    results = model.analyze_buckling(num_modes=1, plane='XY')

    for x in (None, 1.0):
        assert get_buckling_length(model, 'Column', results, plane='z', x=x) == \
            results.effective_length('Column', plane='z', x=x)

    with pytest.raises(NameError):
        get_buckling_length(model, 'Missing', results)
