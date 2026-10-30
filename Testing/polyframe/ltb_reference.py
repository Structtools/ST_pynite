"""Independent lateral-torsional buckling reference: Rayleigh-Ritz on Vlasov theory.

A third opinion for the cases where PyNite and PolyFrame disagree. It shares no code
with either: one prismatic, doubly symmetric member; lateral displacement ``v(x)``
and twist ``phi(x)`` expanded in polynomial series that satisfy the essential
boundary conditions; strain energy with the **true warping term** ``E Iw phi''^2``
(no equivalent torsion constant, no element, no mesh); second-order work of the
in-plane moment ``M(x)``, the axial force ``N`` and a distributed load ``q`` acting at
height ``a`` above the shear centre. The critical load factor is the smallest positive
eigenvalue of the resulting generalised eigenproblem.

Because the warping energy is integrated directly, the solver does not assume any
warping length: a member whose twist is held at one end only simply has that
condition in its basis. That is what makes it an independent check of the warping
length PyNite's equivalent-torsion approximation needs for such a member. The exact
answer for pure torsional buckling is instructive: with the twist held at one end and
warping free everywhere, the critical mode is a uniform twist rate, which has no
warping energy, so ``N_cr,T = G J / i0^2`` exactly as if ``Iw = 0``; the textbook
``L_T = 2 L`` belongs to the case where warping is *fixed* at the held end
(``warping_fixed_at``). Under bending the warping contribution of a member held at
one end lies in between, which only a solver like this one can say.

Potential energy (Trahair, "Flexural-Torsional Buckling of Structures"; Timoshenko &
Gere ch. 6), with ``M`` positive sagging and ``q`` positive downward::

    Pi = 1/2 int (E Iz v''^2 + G J phi'^2 + E Iw phi''^2) dx          strain energy
       + lambda * [ int M v'' phi dx                                   moment
                  - 1/2 int N (v'^2 + i0^2 phi'^2) dx                  axial force
                  - 1/2 int q a phi^2 dx ]                             load height

so ``K a = lambda Kg a`` with ``K`` the strain-energy matrix and ``Kg`` the (indefinite)
second-order work matrix. The sign of the moment coupling only mirrors the eigenvector;
the signs of the ``N`` and ``q a`` terms are the physical ones (compression and a load
above the shear centre destabilise) and are checked against closed forms in
``Testing/test_ltb_reference.py``.

Boundary conditions are stated per quantity: ``v`` held at both ends (every case here;
PolyFrame holds the lateral displacement at every node), and ``phi`` held at both ends
(fork supports), at one end, or at the ends and interior points (a continuous beam's
supports). Warping is free wherever it is not held by a basis constraint, which is the
fork-support assumption of both tools.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from numpy.polynomial import Polynomial, legendre
from scipy.linalg import eigh

Function = Callable[[np.ndarray], np.ndarray]


@dataclass(frozen=True)
class MemberLTB:
    """A prismatic doubly symmetric member along x from 0 to L. SI units."""

    L: float
    E: float
    G: float
    Iz: float  # second moment of area about the weak axis (lateral bending)
    J: float  # St Venant torsion constant
    Iw: float  # warping constant
    i0_sq: float  # polar radius of gyration squared, (Iy + Iz) / A
    v_held_at: tuple[float, ...] = field(default=(0.0, 1.0))  # x / L
    phi_held_at: tuple[float, ...] = field(default=(0.0, 1.0))  # x / L
    warping_fixed_at: tuple[float, ...] = field(default=())  # x / L, phi' = 0 as well


def _basis(held_at: tuple[float, ...], n: int, clamped_at: tuple[float, ...] = ()) -> list[Polynomial]:
    """``n`` polynomials in xi = 2x/L - 1 that vanish where the quantity is held.

    ``clamped_at`` positions get a double root, so the first derivative vanishes
    there too (warping fixed). Shifted Legendre polynomials times the product of the
    constraint factors keep the Gram matrix well conditioned up to n of about 20.
    """
    constraint = Polynomial([1.0])
    for position in held_at:
        xi0 = 2.0 * position - 1.0
        constraint = constraint * Polynomial([-xi0, 1.0])
    for position in clamped_at:
        xi0 = 2.0 * position - 1.0
        constraint = constraint * Polynomial([-xi0, 1.0])
    out = []
    for k in range(n):
        coefficients = np.zeros(k + 1)
        coefficients[k] = 1.0
        legendre_k = Polynomial(legendre.leg2poly(coefficients))
        out.append(constraint * legendre_k)
    return out


def critical_load_factor(
    member: MemberLTB,
    M: Function | None = None,
    N: float = 0.0,
    q: Function | None = None,
    a: float = 0.0,
    n_v: int = 12,
    n_phi: int = 12,
    n_gauss: int = 96,
) -> float:
    """Smallest positive load factor at which the member buckles laterally.

    :param M: in-plane bending moment along the member, positive sagging, as a
        function of x in metres (N m). ``None`` for no moment.
    :param N: axial force, compression positive (N).
    :param q: distributed transverse load, positive downward, function of x (N/m).
        Only its second-order work through the load height is used; its moment is
        part of ``M``.
    :param a: height of the load above the shear centre (m); positive is destabilising
        for a downward load (top flange of a beam under gravity).
    """
    L = member.L
    xi, w = legendre.leggauss(n_gauss)
    x = 0.5 * L * (xi + 1.0)
    w = 0.5 * L * w  # weights on [0, L]
    scale = 2.0 / L  # d xi / d x

    def evaluate(basis: list[Polynomial], order: int) -> np.ndarray:
        rows = [p.deriv(order)(xi) * scale**order if order else p(xi) for p in basis]
        return np.array(rows)

    bv = _basis(member.v_held_at, n_v)
    bp = _basis(member.phi_held_at, n_phi, member.warping_fixed_at)
    v1, v2 = evaluate(bv, 1), evaluate(bv, 2)
    p0, p1, p2 = evaluate(bp, 0), evaluate(bp, 1), evaluate(bp, 2)

    EIz, GJ, EIw = member.E * member.Iz, member.G * member.J, member.E * member.Iw
    K = np.zeros((n_v + n_phi, n_v + n_phi))
    K[:n_v, :n_v] = (v2 * w) @ v2.T * EIz
    K[n_v:, n_v:] = (p1 * w) @ p1.T * GJ + (p2 * w) @ p2.T * EIw

    Kg = np.zeros_like(K)
    if M is not None:
        m = M(x)
        coupling = (v2 * (w * m)) @ p0.T  # int M v'' phi
        Kg[:n_v, n_v:] = -coupling
        Kg[n_v:, :n_v] = -coupling.T
    if N:
        Kg[:n_v, :n_v] += (v1 * w) @ v1.T * N
        Kg[n_v:, n_v:] += (p1 * w) @ p1.T * N * member.i0_sq
    if q is not None and a:
        Kg[n_v:, n_v:] += (p0 * (w * q(x) * a)) @ p0.T

    # Kg a = mu K a; lambda = 1 / mu, smallest positive lambda = largest positive mu
    mu = eigh(Kg, K, eigvals_only=True)
    positive = mu[mu > 1e-12]
    if positive.size == 0:
        return float("inf")
    return float(1.0 / positive.max())


# --- closed forms for the validation tests ------------------------------------


def m_cr_uniform_fork(member: MemberLTB) -> float:
    """Fork-supported beam in uniform bending: (pi/L) sqrt(E Iz (G J + pi^2 E Iw / L^2))."""
    L = member.L
    return (np.pi / L) * np.sqrt(
        member.E * member.Iz * (member.G * member.J + np.pi**2 * member.E * member.Iw / L**2)
    )


def m_cr_three_factor(member: MemberLTB, C1: float, C2: float, z_g: float) -> float:
    """ENV 1993-1-1 Annex F / NCCI SN003 three-factor formula, k = k_w = 1."""
    L = member.L
    EIz = member.E * member.Iz
    root = np.sqrt(
        member.Iw / member.Iz
        + L**2 * member.G * member.J / (np.pi**2 * EIz)
        + (C2 * z_g) ** 2
    )
    return C1 * np.pi**2 * EIz / L**2 * (root - C2 * z_g)


def n_cr_torsional(member: MemberLTB, L_T: float) -> float:
    """Torsional buckling of a column, twist held at two points L_T apart (or 2 L for one)."""
    return (member.G * member.J + np.pi**2 * member.E * member.Iw / L_T**2) / member.i0_sq
