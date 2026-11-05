"""
Closed-form checks of the independent Rayleigh-Ritz lateral-torsional buckling
reference in ``Testing/polyframe/ltb_reference.py``, which is used as the third
opinion where PyNite and PolyFrame disagree.

1. Fork-supported beam, uniform moment, with warping: exact.
2. Fork-supported beam, UDL and central point load, Iw = 0: Timoshenko & Gere's
   28.3 and 16.94 coefficients (as the PyNite tests use).
3. Fork-supported beam, UDL on the top flange and the bottom flange, with warping:
   the three-factor formula (C1 = 1.132, C2 = 0.459) to within its own accuracy.
4. Torsional buckling of a column with twist held at one end only. Exact Vlasov
   solution: with warping free, the critical mode is a uniform twist rate with no
   warping energy, N_cr,T = G J / i0^2 (as if Iw = 0); with warping also fixed at the
   held end, L_T = 2 L; with twist held at both ends, L_T = L. The solver is told
   only the boundary conditions and must find all three.
5. Beam-column (uniform moment plus compression, fork supports): the exact
   interaction (M / M_cr0)^2 = (1 - N / N_z)(1 - N / N_T).

Units: SI.
"""

import numpy as np
import pytest

from Testing.polyframe.ltb_reference import (
    MemberLTB,
    critical_load_factor,
    m_cr_three_factor,
    m_cr_uniform_fork,
    n_cr_torsional,
)
from Testing.polyframe.sections import SECTIONS, STEEL_E, STEEL_G


def _ipe300(L=6.0, **kw):
    s = SECTIONS["IPE300"]
    return MemberLTB(
        L=L, E=STEEL_E, G=STEEL_G, Iz=s.Iz, J=s.It, Iw=s.Iw, i0_sq=(s.Iy + s.Iz) / s.A, **kw
    )


def test_uniform_moment_fork_supports_with_warping():
    member = _ipe300()
    factor = critical_load_factor(member, M=lambda x: np.ones_like(x))
    assert factor == pytest.approx(m_cr_uniform_fork(member), rel=1e-4)


@pytest.mark.parametrize(
    "load, coefficient",
    [
        ("udl", 28.3),
        ("point", 16.94),
    ],
)
def test_timoshenko_gere_coefficients_without_warping(load, coefficient):
    s = SECTIONS["IPE300"]
    member = MemberLTB(L=6.0, E=STEEL_E, G=STEEL_G, Iz=s.Iz, J=s.It, Iw=0.0, i0_sq=(s.Iy + s.Iz) / s.A)
    L = member.L
    root = np.sqrt(member.E * member.Iz * member.G * member.J)
    if load == "udl":
        M = lambda x: 1.0 * x * (L - x) / 2.0  # q = 1 N/m
        expected = coefficient * root / L**3  # q_cr
    else:
        M = lambda x: 1.0 * np.minimum(x, L - x) / 2.0  # P = 1 N at midspan
        expected = coefficient * root / L**2  # P_cr
    factor = critical_load_factor(member, M=M, n_v=16, n_phi=16)
    assert factor == pytest.approx(expected, rel=0.01)


@pytest.mark.parametrize("sign", [+1.0, -1.0])
def test_udl_load_height_matches_three_factor_formula(sign):
    """Top flange (destabilising) and bottom flange (stabilising), to the formula's accuracy."""
    member = _ipe300()
    L = member.L
    z_g = sign * 0.15
    q = 1.0  # N/m
    M = lambda x: q * x * (L - x) / 2.0
    factor = critical_load_factor(member, M=M, q=lambda x: q * np.ones_like(x), a=z_g)
    m_cr_ritz = factor * q * L**2 / 8.0
    m_cr_formula = m_cr_three_factor(member, C1=1.132, C2=0.459, z_g=z_g)
    assert m_cr_ritz == pytest.approx(m_cr_formula, rel=0.04)
    # and the direction is right
    m_cr_centre = critical_load_factor(member, M=M) * q * L**2 / 8.0
    assert (m_cr_ritz < m_cr_centre) == (sign > 0)


def test_torsional_buckling_end_conditions_are_found_from_the_warping_energy():
    s = SECTIONS["IPE300"]
    L = 4.0
    common = dict(L=L, E=STEEL_E, G=STEEL_G, Iz=1e3 * s.Iz, J=s.It, Iw=s.Iw, i0_sq=(s.Iy + s.Iz) / s.A)
    # Iz is inflated so that the flexural mode sits far above the torsional one
    one_end_free_warping = MemberLTB(**common, phi_held_at=(0.0,))
    one_end_fixed_warping = MemberLTB(**common, phi_held_at=(0.0,), warping_fixed_at=(0.0,))
    both_ends = MemberLTB(**common, phi_held_at=(0.0, 1.0))

    gj_over_i0 = s.It * STEEL_G / common["i0_sq"]
    assert critical_load_factor(one_end_free_warping, N=1.0, n_phi=16) == pytest.approx(gj_over_i0, rel=2e-3)
    assert critical_load_factor(one_end_fixed_warping, N=1.0, n_phi=16) == pytest.approx(
        n_cr_torsional(one_end_fixed_warping, 2.0 * L), rel=2e-3
    )
    assert critical_load_factor(both_ends, N=1.0, n_phi=16) == pytest.approx(
        n_cr_torsional(both_ends, L), rel=1e-3
    )


def test_beam_column_interaction_fork_supports():
    member = _ipe300(L=4.0)
    N_z = np.pi**2 * member.E * member.Iz / member.L**2
    N_T = n_cr_torsional(member, member.L)
    M0 = m_cr_uniform_fork(member)
    N, M = 100e3, 10e3  # N, N m, scaled together by the factor
    factor = critical_load_factor(member, M=lambda x: M * np.ones_like(x), N=N)
    lhs = (factor * M / M0) ** 2
    rhs = (1 - factor * N / N_z) * (1 - factor * N / N_T)
    assert lhs == pytest.approx(rhs, rel=1e-3)
