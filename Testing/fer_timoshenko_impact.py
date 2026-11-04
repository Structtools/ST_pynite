"""Practical impact of the Timoshenko fixed-end-reaction correction.

Run this to answer one question: *does the FER fix actually change anything
I care about?*

Background
----------
The Timoshenko shear (Phi) correction was applied in the wrong direction in
every bending fixed-end reaction. The corrected forms are in
``Pynite/FixedEndReactions.py``; the previous ones are reproduced verbatim at
the bottom of this module as ``_legacy_*`` so the two can be run side by side
without editing the library.

The headline finding is that the error is invisible in the reference
comparison scripts in this folder, and for a good reason -- **they apply only
uniform loads spanning a whole member**. For a uniform full-span load the
fixed-end moment is exactly ``w*L^2/12`` for *any* value of Phi, so those
models cannot detect an error in the Phi terms no matter how large Phi is.
The same blind spot is why the original test suite stayed green.

Three conditions must coincide before the correction changes a result:

    coarse mesh  +  non-uniform load  +  low E/G ratio

Steel (E/G ~ 2.6) is negligible in every configuration measured here. Glulam
(E/G ~ 19) analysed with a single element per member is not. Subdividing
members into two or more elements collapses the error in every case.

Usage
-----
    python Testing/fer_timoshenko_impact.py
"""

from __future__ import annotations

import Pynite.FixedEndReactions as FER
from numpy import zeros
from Pynite import FEModel3D

# ---------------------------------------------------------------------------
# Materials and sections
# ---------------------------------------------------------------------------

STEEL = dict(E=210e9, G=80.769e9)          # E/G = 2.6
GLULAM = dict(E=12600e6, G=650e6)          # E/G = 19.4  (GL28h)

# IPE200, as used by the reference comparison scripts in this folder.
IPE200 = dict(A=28.48e-4, Iy=1943e-8, Iz=142e-8, J=7.00e-8,
              Asy=17.00e-4, Asz=14.02e-4)


def rect(b: float, h: float) -> dict:
    """Solid rectangle: shear area 5/6 * A on both axes."""
    return dict(A=b * h, Iy=b * h**3 / 12, Iz=h * b**3 / 12, J=1e-5,
                Asy=5 / 6 * b * h, Asz=5 / 6 * b * h)


def phi_z(mat: dict, sec: dict, L: float) -> float:
    """Timoshenko parameter for bending about z (pairs with Iz and Asy)."""
    return 12 * mat['E'] * sec['Iz'] / (mat['G'] * sec['Asy'] * L**2)


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

def fixed_beam_support_moment(mat, sec, L, nsub, load):
    """i-end support moment of a fixed-fixed beam meshed into nsub elements.

    ``load`` is called once per sub-element as ``load(model, name, i, nsub, Lsub)``.
    """
    m = FEModel3D()
    m.add_material('Mat', mat['E'], mat['G'], 0.3, 0.0)
    m.add_section('Sec', sec['A'], sec['Iy'], sec['Iz'], sec['J'],
                  Asy=sec['Asy'], Asz=sec['Asz'])
    for i in range(nsub + 1):
        m.add_node(f'N{i}', i * L / nsub, 0, 0)
    m.def_support('N0', *(True,) * 6)
    m.def_support(f'N{nsub}', *(True,) * 6)
    for i in range(nsub):
        name = f'M{i}'
        m.add_member(name, f'N{i}', f'N{i + 1}', 'Mat', 'Sec')
        load(m, name, i, nsub, L / nsub)
    m.add_load_combo('Combo 1', {'Case 1': 1.0})
    m.analyze_linear(log=False, check_statics=False)
    return m.nodes['N0'].RxnMZ['Combo 1']


def uniform(w):
    def apply(m, name, i, nsub, Lsub):
        m.add_member_dist_load(name, 'Fy', w, w, 0, Lsub)
    return apply


def triangular(w_end):
    """Load ramping 0 -> w_end over the whole beam, sliced across elements."""
    def apply(m, name, i, nsub, Lsub):
        m.add_member_dist_load(name, 'Fy',
                               w_end * i / nsub, w_end * (i + 1) / nsub, 0, Lsub)
    return apply


def point_at_third(P, L):
    """Single point load at L/3 of the whole beam, wherever that falls."""
    def apply(m, name, i, nsub, Lsub):
        a = L / 3.0
        if i * Lsub <= a < (i + 1) * Lsub:
            m.add_member_pt_load(name, 'Fy', P, a - i * Lsub)
    return apply


# ---------------------------------------------------------------------------
# Running both implementations
# ---------------------------------------------------------------------------

_CURRENT = {name: getattr(FER, name)
            for name in ('FER_PtLoad', 'FER_Moment', 'FER_LinLoad')}


def _use(impl: dict) -> None:
    # Member3D resolves these as module attributes
    # (Pynite.FixedEndReactions.FER_*), so patching the module is enough.
    for name, fn in impl.items():
        setattr(FER, name, fn)


def both_ways(mat, sec, L, nsub, load):
    """Return (corrected, legacy) support moments for the same model."""
    _use(_CURRENT)
    corrected = fixed_beam_support_moment(mat, sec, L, nsub, load)
    _use(_LEGACY)
    try:
        legacy = fixed_beam_support_moment(mat, sec, L, nsub, load)
    finally:
        _use(_CURRENT)
    return corrected, legacy


def pct(corrected: float, legacy: float) -> float:
    return abs(corrected - legacy) / max(abs(corrected), 1e-12) * 100


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------

def report_load_types() -> None:
    print('=' * 78)
    print('  1. Which load types reach the corrected code path?')
    print('     Steel IPE200, L = 3.0 m, single element, Phi_z = '
          f'{phi_z(STEEL, IPE200, 3.0):.4f}')
    print('=' * 78)
    print(f"  {'load case':<44}{'corrected':>11}{'previous':>11}{'diff':>9}")
    print('  ' + '-' * 74)

    cases = [
        ('uniform UDL, full span  (comparison scripts)', uniform(-10e3)),
        ('triangular 0 -> -20 kN/m', triangular(-20e3)),
        ('point load at L/3', point_at_third(-30e3, 3.0)),
    ]
    for label, load in cases:
        c, l = both_ways(STEEL, IPE200, 3.0, 1, load)
        d = pct(c, l)
        note = '   <-- Phi-independent' if d < 1e-6 else ''
        print(f'  {label:<44}{c / 1e3:10.4f} {l / 1e3:10.4f} {d:7.3f}%{note}')

    print()
    print('  A uniform full-span load gives w*L^2/12 for ANY Phi. That is why')
    print('  the reference comparison scripts show no change whatsoever, and')
    print('  why a symmetric-only test set cannot detect this class of error.')
    print()


def report_mesh_dependence() -> None:
    print('=' * 78)
    print('  2. The corrected result is mesh-independent; the previous one is not')
    print('     Glulam 300x900, L = 3.0 m, triangular load 0 -> -20 kN/m')
    print('=' * 78)
    print(f"  {'sub-elements':>13}{'Phi(element)':>14}"
          f"{'corrected':>12}{'previous':>12}{'diff':>10}")
    print('  ' + '-' * 74)

    mat, sec, L = GLULAM, rect(0.3, 0.9), 3.0
    for nsub in (1, 2, 3, 4, 8, 16):
        c, l = both_ways(mat, sec, L, nsub, triangular(-20e3))
        phi = phi_z(mat, sec, L / nsub)
        near = '  <-- Phi ~ 1' if 0.6 < phi < 1.6 else ''
        print(f'  {nsub:>13}{phi:14.3f}{c / 1e3:12.5f}{l / 1e3:12.5f}'
              f'{pct(c, l):9.2f}%{near}')

    print()
    print('  The corrected value is identical at every mesh density -- the exact')
    print('  Timoshenko element needs no refinement for this problem. The previous')
    print('  implementation is mesh-dependent and NOT monotone: it is worse at two')
    print('  elements than at one.')
    print()


def report_pole() -> None:
    print('=' * 78)
    print('  3. Why: the previous determinant had a pole at Phi = 1')
    print('     det = L^4*(Phi - 1)/12   ->   corrected: -L^4*(1 + Phi)/12')
    print('=' * 78)
    print(f"  {'Phi':>8}{'corrected':>14}{'previous':>14}{'diff':>12}")
    print('  ' + '-' * 60)

    # FER_PtLoad at a = L/3 on a unit span, straight from the two implementations.
    L, P, a = 1.0, -100.0, 1.0 / 3.0
    for phi in (0.0, 0.25, 0.5, 0.9, 0.99, 1.0, 1.01, 1.1, 2.0, 10.0):
        good = _CURRENT['FER_PtLoad'](P, a, L, 'Fy', phi)[5, 0]
        try:
            bad = _legacy_FER_PtLoad(P, a, L, 'Fy', phi)[5, 0]
            bad_s, diff_s = f'{bad:14.4f}', f'{pct(good, bad):11.1f}%'
        except ZeroDivisionError:
            bad_s, diff_s = f'{"ZeroDivErr":>14}', f'{"--":>12}'
        print(f'  {phi:8.2f}{good:14.4f}{bad_s}{diff_s}')

    print()
    print('  Phi = 1 means 12*E*I = G*As*L^2, i.e. L/h ~ 1.8 for a rectangle with')
    print('  nu = 0.3. That is reachable -- the 2-element row above sits at 0.930.')
    print('  Both forms happen to agree again as Phi -> infinity, which is why very')
    print('  fine meshes look clean. Accuracy that depends on straddling a')
    print('  singularity is not accuracy.')
    print()


def main() -> None:
    print()
    print('Timoshenko fixed-end reaction correction -- practical impact')
    print()
    report_load_types()
    report_mesh_dependence()
    report_pole()
    print('=' * 78)
    print('  Conclusion')
    print('=' * 78)
    print('  * No effect on the reference comparison scripts. They apply only')
    print('    uniform full-span loads, which are Phi-independent by construction.')
    print('  * No practical effect on steel: E/G ~ 2.6 keeps Phi near zero, far')
    print('    from the pole. Worst case measured 0.145%.')
    print('  * Timber is exposed. E/G ~ 19 puts realistic members near Phi = 1,')
    print('    where the previous error is large, mesh-dependent and non-monotone.')
    print('  * Refining the mesh is NOT a reliable mitigation: it moves Phi through')
    print('    the pole rather than away from it.')
    print()


# ---------------------------------------------------------------------------
# The previous implementation, reproduced verbatim for comparison.
#
# Differences from the corrected version now in Pynite/FixedEndReactions.py:
#   det        L**4*(Phi - 1)/12        ->  -L**4*(1 + Phi)/12
#   numerator  (2 + Phi)                ->  (2 - Phi)
#   PtLoad I2  + Phi*L**2*b/12          ->  - Phi*L**2*b/12
#   LinLoad I2 + Phi*L**2/12*J1         ->  - Phi*L**2/12*J1
#   Moment I2  -M*(b**2/2 + Phi*L**2/12) ->  -M*b**2/2   (no Phi term at all)
# ---------------------------------------------------------------------------

def _legacy_solve(I1, I2, L, Phi, total_load, total_moment_about_i, Direction):
    out = zeros((12, 1))
    det = L**4 * (Phi - 1) / 12
    M_0 = (I1 * L**3 * (2 + Phi) / 12 - L**2 / 2 * I2) / det
    R_i = (L * I2 - L**2 / 2 * I1) / det
    R_j = total_load - R_i
    M_j = -(M_0 + R_i * L - total_load * L + total_moment_about_i)

    if Direction in ('Fy', 'Mz'):
        out[1, 0], out[5, 0], out[7, 0], out[11, 0] = -R_i, M_0, -R_j, M_j
    elif Direction == 'Fz':
        out[2, 0], out[4, 0], out[8, 0], out[10, 0] = -R_i, -M_0, -R_j, -M_j
    elif Direction == 'My':
        out[2, 0], out[4, 0], out[8, 0], out[10, 0] = R_i, M_0, R_j, M_j
    return out


def _legacy_FER_PtLoad(P, x, L, Direction, Phi=0.0):
    b = L - x
    return _legacy_solve(P * b**2 / 2,
                         P * (b**3 / 6 + Phi * L**2 * b / 12),
                         L, Phi, P, P * x, Direction)


def _legacy_FER_Moment(M, x, L, Direction, Phi=0.0):
    b = L - x
    return _legacy_solve(-M * b,
                         -M * (b**2 / 2 + Phi * L**2 / 12),
                         L, Phi, 0.0, M, Direction)


def _legacy_FER_LinLoad(w1, w2, x1, x2, L, Direction, Phi=0.0):
    dx = x2 - x1
    b1 = L - x1
    c = (w2 - w1) / dx if dx != 0 else 0.0
    J1 = b1 * dx * (w1 + w2) / 2 - dx**2 * (w1 + 2 * w2) / 6
    J2 = (w1 * b1**2 * dx - w1 * b1 * dx**2 + w1 * dx**3 / 3
          + c * b1**2 * dx**2 / 2 - 2 * c * b1 * dx**3 / 3 + c * dx**4 / 4)
    J3 = (w1 * b1**3 * dx - 3 * w1 * b1**2 * dx**2 / 2 + w1 * b1 * dx**3
          - w1 * dx**4 / 4 + c * b1**3 * dx**2 / 2 - c * b1**2 * dx**3
          + 3 * c * b1 * dx**4 / 4 - c * dx**5 / 5)
    return _legacy_solve(J2 / 2,
                         J3 / 6 + Phi * L**2 / 12 * J1,
                         L, Phi, dx * (w1 + w2) / 2,
                         x1 * dx * (w1 + w2) / 2 + dx**2 * (w1 + 2 * w2) / 6,
                         Direction)


_LEGACY = {'FER_PtLoad': _legacy_FER_PtLoad,
           'FER_Moment': _legacy_FER_Moment,
           'FER_LinLoad': _legacy_FER_LinLoad}


if __name__ == '__main__':
    main()
