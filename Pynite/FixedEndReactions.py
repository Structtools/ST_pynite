# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 20:58:03 2017

@author: D. Craig Brinck, SE
"""
# %%
from __future__ import annotations # Allows more recent type hints features
from typing import TYPE_CHECKING

from numpy import zeros

if TYPE_CHECKING:
    from typing import Literal
    from numpy import float64
    from numpy.typing import NDArray


def _solve_FER_bending(I1: float, I2: float, L: float, Phi: float, total_load: float, total_moment_about_i: float, Direction: str) -> 'NDArray[float64]':
    """Solve the 2x2 Timoshenko fixed-end beam system for bending FER.

    This solves the compatibility equations for a fixed-fixed Timoshenko beam
    to find the exact end reactions. When Phi=0 the results reduce to Euler-Bernoulli.

    Parameters
    ----------
    I1 : float
        integral of load(s) * b(s)^2 / 2 ds  (b = L - s)
    I2 : float
        integral_0^L (L - x) * Q(x) dx  -  (Phi*L^2/12) * S

        where Q(x) is the moment at x of the applied loads to the left of x,
        and S is the *shear impulse* of the applied loading:
            S = Q(L)  for transverse loads (point loads, distributed loads)
            S = 0     for applied couples

        The distinction matters. A concentrated couple makes M(x) jump but
        leaves V(x) continuous, so it contributes nothing to the shear
        compatibility integral and therefore carries no Phi term at all.
    L : float
        Member length
    Phi : float
        Timoshenko shear deformation parameter (12*E*I / (G*As*L^2)), 0 for EB
    total_load : float
        Total applied transverse load (integral of w ds, or P for point load)
    total_moment_about_i : float
        Total moment of the applied loads about the i-end
    Direction : str
        'Fy' or 'Fz' (or 'Mz'/'My' for moment loading sign mapping)
    """
    FER = zeros((12, 1))

    # The compatibility conditions are written in the (L - x)-weighted form,
    # which flips the sign of every Phi term relative to the x-weighted
    # derivation. Hence (2 - Phi) and (1 + Phi), not (2 + Phi) and (Phi - 1).
    #
    # Determinant of the 2x2 system:
    # | L       L^2/2            | | M_0 |   | I1 |
    # | L^2/2   L^3*(2-Phi)/12   | | R_i | = | I2 |
    #
    # det = -L^4*(1 + Phi)/12 never vanishes for Phi >= 0. Every correct
    # Timoshenko quantity carries (1 + Phi); a (Phi - 1) denominator here
    # would introduce a spurious pole at Phi = 1 (L/h ~ 1.77 for a
    # rectangular section with nu = 0.3, i.e. entirely reachable).
    det = -L**4 * (1 + Phi) / 12

    # Cramer's rule
    M_0 = (I1 * L**3 * (2 - Phi) / 12 - L**2 / 2 * I2) / det
    R_i = (L * I2 - L**2 / 2 * I1) / det

    R_j = total_load - R_i
    # Internal moment at j: M_at_j = M_0 + R_i*L - (moment of applied loads about j)
    # where moment_about_j = total_load*L - total_moment_about_i
    M_at_j = M_0 + R_i * L - total_load * L + total_moment_about_i
    # For the FER vector, the j-end moment entry (Fy direction) is -M_at_j
    M_j = -M_at_j

    # Map to 12-DOF FER vector.
    # Fy-plane (y forces, z moments): standard sign convention
    # Fz-plane (z forces, y moments): moment signs are negated
    # My direction: force signs are negated (vs Mz), moment signs same
    if Direction == 'Fy':
        FER[1, 0] = -R_i
        FER[5, 0] = M_0
        FER[7, 0] = -R_j
        FER[11, 0] = M_j
    elif Direction == 'Fz':
        FER[2, 0] = -R_i
        FER[4, 0] = -M_0
        FER[8, 0] = -R_j
        FER[10, 0] = -M_j
    elif Direction == 'Mz':
        FER[1, 0] = -R_i
        FER[5, 0] = M_0
        FER[7, 0] = -R_j
        FER[11, 0] = M_j
    elif Direction == 'My':
        FER[2, 0] = R_i
        FER[4, 0] = M_0
        FER[8, 0] = R_j
        FER[10, 0] = M_j

    return FER


# %%
def FER_PtLoad(P: float, x: float, L: float, Direction: Literal["Fy", "Fz"], Phi: float = 0.0) -> NDArray[float64]:
    """
    Returns the fixed end reaction vector for a point load

    Parameters:
    -----------
    P : float
        The magnitude of the point load
    x : float
        The location of the point load relative to the start of the member
    L : float
        The length of the member
    Direction : Literal["Fy", "Fz"]
        The direction of the point load. Must be one of the following:
            "Fy" = Force on the member's local y-axis
            "Fz" = Force on the member's local z-axis
    Phi : float
        Timoshenko shear deformation parameter (12*E*I/(G*As*L^2)). Default 0.0 (Euler-Bernoulli).
    """
    b = L - x

    # Integrals for point load P at position x: b(s) = L - s evaluated at s = x.
    # Shear impulse S = P*b, giving the -(Phi*L^2/12)*S term below.
    I1 = P * b**2 / 2
    I2 = P * (b**3 / 6 - Phi * L**2 * b / 12)
    total_load = P
    total_moment_about_i = P * x

    return _solve_FER_bending(I1, I2, L, Phi, total_load, total_moment_about_i, Direction)


def FER_Moment(M: float, x: float, L: float, Direction: Literal["My", "Mz"], Phi: float = 0.0) -> NDArray[float64]:
    """
    Returns the fixed end reaction vector for a concentrated moment

    Parameters
    ----------
    M : float
        The magnitude of the moment
    x : float
        The location of the moment relative to the start of the member
    L : float
        The length of the member
    Direction : Literal["My", "Mz"]
        The direction of the moment. Must be one of the following:
            "My" = Moment applied about the local y-axis
            "Mz" = Moment applied about the local z-axis
    Phi : float
        Timoshenko shear deformation parameter (12*E*I/(G*As*L^2)). Default 0.0 (Euler-Bernoulli).
    """
    b = L - x

    # For an applied moment M at position x on a fixed-fixed Timoshenko beam.
    #
    # There is NO Phi term here. A concentrated couple makes the bending
    # moment M(x) jump but leaves the shear V(x) continuous, so its shear
    # impulse S is zero and it drops out of the shear compatibility integral.
    # The Phi dependence of the result enters solely through the (1 + Phi)
    # and (2 - Phi) factors in _solve_FER_bending.
    I1 = -M * b
    I2 = -M * b**2 / 2
    total_load = 0.0
    total_moment_about_i = M

    return _solve_FER_bending(I1, I2, L, Phi, total_load, total_moment_about_i, Direction)


def FER_LinLoad(w1: float, w2: float, x1: float, x2: float, L: float, Direction: Literal["Fy", "Fz"], Phi: float = 0.0) -> NDArray[float64]:
    """
    Returns the fixed end reaction vector for a linear distributed load

    Parameters
    ----------
    w1 : float
        The load magnitude at the start location
    w2 : float
        The load magnitude at the end location
    x1 : float
        The start location of the distributed load
    x2 : float
        The end location of the distributed load
    L : float
        The length of the member
    Direction : Literal["Fy", "Fz"]
        The direction of the distributed load. Must be one of the following:
            "Fy" = Force on the member's local y-axis
            "Fz" = Force on the member's local z-axis
    Phi : float
        Timoshenko shear deformation parameter (12*E*I/(G*As*L^2)). Default 0.0 (Euler-Bernoulli).
    """
    # Closed-form integrals for linear load w(s) = w1 + (w2-w1)*(s-x1)/dx
    # over [x1, x2], with b(s) = L - s.
    # Parametrize: u = s - x1, b1 = L - x1, c = (w2-w1)/dx, dx = x2 - x1
    dx = x2 - x1
    b1 = L - x1
    c = (w2 - w1) / dx if dx != 0 else 0.0

    # J_n = integral_0^dx (w1 + c*u) * (b1 - u)^n du
    # J1 (n=1): b1*dx*(w1+w2)/2 - dx^2*(w1+2*w2)/6
    J1 = b1 * dx * (w1 + w2) / 2 - dx**2 * (w1 + 2*w2) / 6

    # J2 (n=2): w1*b1^2*dx - w1*b1*dx^2 + w1*dx^3/3 + c*b1^2*dx^2/2 - 2*c*b1*dx^3/3 + c*dx^4/4
    J2 = (w1*b1**2*dx - w1*b1*dx**2 + w1*dx**3/3
          + c*b1**2*dx**2/2 - 2*c*b1*dx**3/3 + c*dx**4/4)

    # J3 (n=3): w1*b1^3*dx - 3*w1*b1^2*dx^2/2 + w1*b1*dx^3 - w1*dx^4/4
    #           + c*b1^3*dx^2/2 - c*b1^2*dx^3 + 3*c*b1*dx^4/4 - c*dx^5/5
    J3 = (w1*b1**3*dx - 3*w1*b1**2*dx**2/2 + w1*b1*dx**3 - w1*dx**4/4
          + c*b1**3*dx**2/2 - c*b1**2*dx**3 + 3*c*b1*dx**4/4 - c*dx**5/5)

    I1 = J2 / 2
    # Shear impulse S = J1 (the total load moment arm integral), so the
    # Timoshenko term is subtracted, matching FER_PtLoad.
    I2 = J3 / 6 - Phi * L**2 / 12 * J1

    # Total load and moment about i-end
    total_load = dx * (w1 + w2) / 2
    # integral of w(s)*s ds from x1 to x2:
    total_moment_about_i = x1 * dx * (w1 + w2) / 2 + dx**2 * (w1 + 2*w2) / 6

    return _solve_FER_bending(I1, I2, L, Phi, total_load, total_moment_about_i, Direction)


# Returns the fixed end reaction vector for an axial point load
def FER_AxialPtLoad(P: float, x: float, L: float) -> NDArray[float64]:
    """
    Returns the fixed end reaction vector for an axial point load

    Parameters
    ----------
    P : float
        The magnitude of the axial point load
    x : float
        The location of the axial point load relative to the start of the member
    L : float
        The length of the member
    """

    # Create the fixed end reaction vector
    FER = zeros((12, 1))

    # Populate the fixed end reaction vector
    FER[0, 0] = -P*(L-x)/L
    FER[6, 0] = -P*x/L

    return FER


# Returns the fixed end reaction vector for a distributed axial load
def FER_AxialLinLoad(p1: float, p2: float, x1: float, x2: float, L: float) -> NDArray[float64]:
    """
    Returns the fixed end reaction vector for a distributed axial load

    Parameters
    ----------
    p1 : float
        The axial load magnitude at the start location
    p2 : float
        The axial load magnitude at the end location
    x1 : float
        The start location of the distributed axial load
    x2 : float
        The end location of the distributed axial load
    L : float
        The length of the member
    """

    # Create the fixed end reaction vector
    FER = zeros((12, 1))

    # Populate the fixed end reaction vector
    FER[0, 0] = 1/(6*L)*(x1-x2)*(3*L*p1+3*L*p2-2*p1*x1-p1*x2-p2*x1-2*p2*x2)
    FER[6, 0] = 1/(6*L)*(x1-x2)*(2*p1*x1+p1*x2+p2*x1+2*p2*x2)

    return FER


def FER_Torque(T: float, x: float, L: float) -> NDArray[float64]:
    """
    Returns the fixed end reaction vector for a concentrated torque

    Parameters
    ----------
    T : float
        The magnitude of the torque
    x : float
        The location of the torque relative to the start of the member
    L : float
        The length of the member
    """

    # Create the fixed end reaction vector
    FER = zeros((12, 1))

    # Populate the fixed end reaction vector
    FER[3, 0] = -T*(L - x)/L
    FER[9, 0] = -T*x/L

    return FER
