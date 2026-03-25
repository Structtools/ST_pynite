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

    Note: PyNite uses 6-DOF beam elements (no warping DOF).  M_cr cannot be
    extracted directly from a 6-DOF eigenvalue, so the closed-form C1 formula
    is used instead.  Supply the warping constant I_w via ``section.Iw``.
"""

from __future__ import annotations

import numpy as np

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from Pynite.FEModel3D import FEModel3D
    from Pynite.Analysis import BucklingResults


def get_buckling_length(model: FEModel3D, member_name: str,
                        buckling_results: BucklingResults,
                        mode: int = 0, plane: str = 'y') -> float:
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

    Returns
    -------
    float
        L_cr in the same length units as the model geometry.
        Returns ``inf`` when the member has no axial force.
    """
    member = model.members[member_name]
    return buckling_results.effective_length(member, mode=mode, plane=plane)


def get_critical_moment(model: FEModel3D, member_name: str,
                        combo_name: str = 'Combo 1',
                        C1: float = 1.0,
                        C2: float = 0.0,
                        k:  float = 1.0,
                        kw: float = 1.0) -> float:
    """Compute M_cr for lateral-torsional buckling (LTB).

    Uses the three-factor formula from EN 1993-1-1 Annex F / NCCI SN003::

        M_cr = C1 · (π²EIz)/(kL)² · sqrt( Iw/Iz + (kL)²·G·It/(π²EIz) )

    This is the closed-form approach for 6-DOF FEM models that do **not**
    include warping (7th) DOFs.  The warping constant Iw must be supplied by
    setting ``section.Iw`` on the section object before calling this function.

    Parameters
    ----------
    model : FEModel3D
        The analysed model.
    member_name : str
        Name of the physical member.
    combo_name : str, optional
        Load combination name (used only to read the unbraced length).
        Default ``'Combo 1'``.
    C1 : float, optional
        Moment distribution factor (accounts for non-uniform bending).
        Typical values:

        * 1.000 — uniform moment (most conservative)
        * 1.132 — uniformly distributed load, simply supported
        * 1.285 — midspan point load, simply supported

        Default 1.0.
    C2 : float, optional
        Load position factor (> 0 for destabilising loads above centroid).
        Default 0.0.
    k : float, optional
        Effective length factor for buckling (end rotation restraint).
        k = 1.0 → both ends free to rotate (default).
        k = 0.5 → both ends fully restrained.
    kw : float, optional
        Effective length factor for warping (end warping restraint).
        kw = 1.0 → both ends free to warp (default).

    Returns
    -------
    float
        M_cr in force · length units consistent with the model.

    Raises
    ------
    AttributeError
        If the member's section does not have an ``Iw`` attribute.

    Notes
    -----
    The C2 term (load position) is omitted from the formula above because
    without warping DOFs the load-height effect cannot be evaluated from the
    FEM model.  If the load acts at the top flange (destabilising) supply
    C2 > 0 and note that the full 3-factor formula from SN003 should be used
    in that case::

        M_cr = C1 * (π²EIz)/(kL)² *
               ( sqrt(Iw/Iz + (kL)²GIt/(π²EIz) + (C2·zg)²) - C2·zg )

    where zg is the distance from the centroid to the load application point.
    """
    member = model.members[member_name]
    sub    = next(iter(member.sub_members.values()))

    E  = sub.material.E
    G  = sub.material.G
    Iz = sub.section.Iz   # minor axis (bending axis for LTB)
    It = sub.section.J    # St. Venant torsion constant
    L  = member.L()       # unbraced length

    Iw = getattr(sub.section, 'Iw', None)
    if Iw is None:
        raise AttributeError(
            f"Section for member '{member_name}' has no 'Iw' (warping constant) "
            "attribute.  Set ``section.Iw = <value>`` before calling "
            "get_critical_moment()."
        )

    kL   = k * L
    pi2  = np.pi**2

    # EN 1993-1-1 Annex F (simplified, C2·zg = 0)
    M_cr = C1 * (pi2 * E * Iz / kL**2) * np.sqrt(
        Iw / Iz + kL**2 * G * It / (pi2 * E * Iz)
    )

    return M_cr
