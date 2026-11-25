"""
EurocodeHelpers.py — EN 1993-1-1 integration utilities for PyNite buckling results.

Functions
---------
get_buckling_length
    Extract L_cr for a specific member from BucklingResults.
    For use in EN 1993-1-1 column buckling checks (N_b,Rd).

get_critical_moment
    Compute M_cr for lateral-torsional buckling using the 3-factor formula
    from EN 1993-1-1 Annex F (NCCI SN003).

    Note: PyNite uses 6-DOF beam elements (no warping DOF). The buckling
    analysis (``FEModel3D.analyze_buckling``) finds lateral-torsional modes
    without warping, or with warping approximated by an equivalent torsion
    constant; this closed-form C1 formula is the complement for a member
    check. Supply the warping constant I_w via ``add_section(..., Iw=...)``.
"""

from __future__ import annotations

import warnings

import numpy as np

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from Pynite.FEModel3D import FEModel3D
    from Pynite.Analysis import BucklingResults


def get_buckling_length(model: FEModel3D, member_name: str,
                        buckling_results: BucklingResults,
                        mode: int = 0, plane: str = 'y',
                        x: float | None = None) -> float:
    """Return the effective buckling length L_cr for a member.

    This is a convenience wrapper around
    :meth:`BucklingResults.effective_length`.  The result is suitable for use
    in EN 1993-1-1 clause 6.3.1 (column buckling) to compute the
    non-dimensional slenderness λ̄ and reduction factor χ.

    Parameters
    ----------
    model : FEModel3D
        The analysed model.
    member_name : str
        Name of the physical member in the model.
    buckling_results : BucklingResults
        Results from :meth:`FEModel3D.analyze_buckling`.
    mode : int, optional
        Buckling mode index (0 = critical/lowest). Default 0.
    plane : str, optional
        ``'y'`` → bending about local y-axis (uses section.Iy),
        ``'z'`` → bending about local z-axis (uses section.Iz).
        Default ``'y'``.
    x : float, optional
        Location along the member at which to take the axial force. Default
        ``None``: the largest compression along the member.

    Returns
    -------
    float
        L_cr in the same length units as the model geometry.
        Returns ``inf`` when the member is not in compression.
    """
    if member_name not in model.members:
        raise NameError(f"Member '{member_name}' does not exist in the model")
    return buckling_results.effective_length(member_name, mode=mode, plane=plane, x=x)


def get_critical_moment(model: FEModel3D, member_name: str,
                        combo_name: str = 'Combo 1',
                        C1: float = 1.0,
                        C2: float = 0.0,
                        k:  float = 1.0,
                        kw: float = 1.0,
                        z_g: float = 0.0,
                        L: float | None = None) -> float:
    """Compute M_cr for lateral-torsional buckling (LTB) of a doubly symmetric member.

    Uses the three-factor formula (ENV 1993-1-1 Annex F / NCCI SN003)::

        M_cr = C1 · π²EIz/(kL)² · ( sqrt( (k/kw)²·Iw/Iz + (kL)²·G·It/(π²EIz) + (C2·zg)² ) − C2·zg )

    where Iz is the second moment of area about the section's minor axis. The minor axis is taken
    as the one with the smaller second moment of area, so the result does not depend on whether
    the section follows Pynite's or the Eurocode's naming of the axes.

    This is the closed-form approach for 6-DOF FEM models that do **not**
    include warping (7th) DOFs.  The warping constant Iw must be supplied on
    the section (``add_section(..., Iw=...)``) before calling this function.

    Parameters
    ----------
    model : FEModel3D
        The analysed model.
    member_name : str
        Name of the physical member.
    combo_name : str, optional
        Not used. Kept for backward compatibility.
    C1 : float, optional
        Moment distribution factor (accounts for non-uniform bending).
        Typical values:

        * 1.000 — uniform moment (most conservative)
        * 1.132 — uniformly distributed load, simply supported
        * 1.285 — midspan point load, simply supported

        Default 1.0.
    C2 : float, optional
        Load position factor, used with `z_g`. Typical values: 0.459 for a
        uniformly distributed load and 0.553 for a midspan point load on a
        simply supported span. Default 0.0.
    k : float, optional
        Effective length factor for buckling (end rotation restraint).
        k = 1.0 → both ends free to rotate (default).
        k = 0.5 → both ends fully restrained.
    kw : float, optional
        Effective length factor for warping (end warping restraint).
        kw = 1.0 → both ends free to warp (default).
        kw = 0.5 → both ends fully restrained against warping.
    z_g : float, optional
        Distance from the shear centre to the point where the load is applied.
        Positive when the load acts towards the shear centre (a gravity load
        on the top flange, destabilising), negative when it acts away from it
        (stabilising), as for ``load_height`` in the member loads. Default 0.0.
    L : float, optional
        Length between the lateral restraints of the member. Default ``None``:
        the length of the physical member.

    Returns
    -------
    float
        M_cr in force · length units consistent with the model.

    Raises
    ------
    AttributeError
        If the member's section has no warping constant ``Iw``.
    ValueError
        If `k`, `kw` or `L` is not positive.

    Warns
    -----
    UserWarning
        If `C2` is given without `z_g`: the load position term is ``C2·z_g``,
        so `C2` alone has no effect.
    """
    member = model.members[member_name]

    E  = member.material.E
    G  = member.material.G
    Iz = min(member.section.Iy, member.section.Iz)  # minor axis (bending axis for LTB)
    It = member.section.J    # St. Venant torsion constant
    L  = member.L() if L is None else L  # length between lateral restraints

    Iw = getattr(member.section, 'Iw', None)
    if Iw is None:
        raise AttributeError(
            f"Section for member '{member_name}' has no 'Iw' (warping constant). "
            "Pass ``Iw=<value>`` to ``add_section`` before calling "
            "get_critical_moment()."
        )

    if not (k > 0 and kw > 0 and L > 0):
        raise ValueError(f'`k`, `kw` and `L` must be positive. k={k}, kw={kw}, L={L} were given.')

    if C2 != 0.0 and z_g == 0.0:
        warnings.warn(
            '`C2` has no effect without `z_g`: the load position term of M_cr is C2·z_g. Give '
            'the distance from the shear centre to the load with `z_g`.',
            UserWarning, stacklevel=2,
        )

    kL   = k * L
    pi2  = np.pi**2
    C2zg = C2 * z_g

    # Three-factor formula (C3 is zero for a doubly symmetric section)
    M_cr = C1 * (pi2 * E * Iz / kL**2) * (np.sqrt(
        (k / kw)**2 * Iw / Iz + kL**2 * G * It / (pi2 * E * Iz) + C2zg**2
    ) - C2zg)

    return M_cr
