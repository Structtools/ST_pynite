"""Verification of modal analysis against closed-form solutions and structural identities.

Units throughout are kN, m and s. Two consequences of that choice drive most of this file:

* Mass comes out in kN·s²/m, which is a tonne.
* `Material.rho` is a *weight* density, not a mass density. `add_member_self_weight` builds the
  self-weight load as `rho*A` and reads it as a force per unit length, and the mass matrix then
  divides that by gravity. So steel is 78.5 kN/m³ here rather than 7.85 t/m³, and `gravity` has to
  be passed as 9.81. Leaving `gravity` at its default of 1.0 overstates the mass by a factor of g
  and understates every frequency by a factor of sqrt(g) -- see `test_gravity_default_is_a_trap`.

The closed-form beam frequency is

    f_n = lambda_n**2/(2*pi*L**2) * sqrt(E*I/mu)

with `mu` the mass per unit length. The simply supported case is exact, with
`lambda_n**2 = n**2*pi**2`, and is the primary oracle. The other support conditions use the
constants below, and the tests confirm them independently: the finite element model does not know
what they are, so agreement to a few hundredths of a percent is a check on both sides.
"""

import numpy as np
import pytest

from Pynite import FEModel3D
from Pynite.ModalResults import ModalModelError

# Material and section properties. Steel in kN and m, with rho as a weight density.
E = 210e6          # kN/m2
G = 81e6           # kN/m2
NU = 0.3
RHO = 78.5         # kN/m3 -- weight density
GRAVITY = 9.81     # m/s2
A = 0.01           # m2
I = 8.333e-6       # m4
J = 1.67e-6        # m4
L = 6.0            # m

MU = RHO*A/GRAVITY  # tonne/m

# Frequency coefficients `lambda_n**2` for a uniform prismatic beam
LAMBDA_SQ = {
    'simply supported': (np.pi**2, 4*np.pi**2, 9*np.pi**2),
    'cantilever': (3.5160, 22.034),
    'fixed-pinned': (15.418, 49.965),
    'fixed-fixed': (22.373, 61.673),
    'free-free': (22.373, 61.673),
}


def beam_frequency(lambda_sq, length=L, EI=E*I, mu=MU):
    """Returns the closed-form natural frequency for a uniform prismatic beam."""

    return lambda_sq/(2*np.pi*length**2)*np.sqrt(EI/mu)


def beam_model(support='simply supported', n_user_nodes=2, self_weight=True, length=L,
               brace_out_of_plane=False):
    """Builds a single-member beam with the requested support condition.

    Only the two end nodes are created unless `n_user_nodes` asks for more, because the modal solver
    does its own subdivision. Every node is left free out of plane and the tests declare
    `plane='XY'`, so the solver restrains the interior nodes it creates as well.

    That leaves the model singular under a *static* solve, because nothing restrains the
    out-of-plane and torsional DOFs. Set `brace_out_of_plane` for the tests that solve statically.
    """

    model = FEModel3D()

    for i in range(n_user_nodes):
        model.add_node(f'N{i}', length*i/(n_user_nodes - 1), 0, 0)

    last = f'N{n_user_nodes - 1}'

    # Restrain the axial direction at one end for every case, so that the axial modes are real
    # rather than rigid-body translation
    if support == 'simply supported':
        model.def_support('N0', True, True, False, False, False, False)
        model.def_support(last, False, True, False, False, False, False)
    elif support == 'cantilever':
        model.def_support('N0', True, True, True, True, True, True)
    elif support == 'fixed-pinned':
        model.def_support('N0', True, True, True, True, True, True)
        model.def_support(last, False, True, False, False, False, False)
    elif support == 'fixed-fixed':
        model.def_support('N0', True, True, True, True, True, True)
        model.def_support(last, True, True, True, True, True, True)
    elif support == 'free-free':
        pass
    else:
        raise ValueError(f'Unknown support condition: {support}')

    # Hold the out-of-plane and torsional DOFs so that a static solve is well posed. This runs after
    # the support conditions above because `def_support` writes all six DOFs at once.
    if brace_out_of_plane:
        for node in model.nodes.values():
            node.support_DZ = True
            node.support_RX = True
            node.support_RY = True

    model.add_material('Steel', E, G, NU, RHO)
    model.add_section('Section', A, I, I, J)
    model.add_member('M1', 'N0', last, 'Steel', 'Section')

    if self_weight:
        model.add_member_self_weight('FY', -1.0, 'SW')
        model.add_load_combo('Mass', {'SW': 1.0})
    else:
        model.add_load_combo('Mass', {})

    return model


def modal(model, num_modes=3, **kwargs):
    """Runs a planar modal analysis with the conventions this file uses throughout."""

    kwargs.setdefault('plane', 'XY')

    return model.analyze_modal(num_modes=num_modes, mass_combo_name='Mass', gravity=GRAVITY,
                               **kwargs)


# ---------------------------------------------------------------------------------------------
# Closed-form oracles
# ---------------------------------------------------------------------------------------------

def test_simply_supported_beam_matches_exact_solution():
    """The simply supported case is exact, so it is the tightest oracle available."""

    results = modal(beam_model('simply supported'), num_modes=3)

    for mode, lambda_sq in enumerate(LAMBDA_SQ['simply supported']):

        expected = beam_frequency(lambda_sq)
        actual = results.frequencies[mode]

        assert actual == pytest.approx(expected, rel=2e-3), (
            f'Mode {mode + 1}: expected {expected:.4f} Hz, got {actual:.4f} Hz'
        )


@pytest.mark.parametrize('support', ['cantilever', 'fixed-pinned', 'fixed-fixed'])
def test_other_support_conditions_match_published_coefficients(support):
    """Confirms the tabulated frequency coefficients for the remaining support conditions.

    The model is built from geometry, stiffness and mass alone and has no knowledge of the
    coefficients, so this is an independent check of the constants in `LAMBDA_SQ` as well as of the
    solver.
    """

    results = modal(beam_model(support), num_modes=4)

    for mode, lambda_sq in enumerate(LAMBDA_SQ[support]):

        expected = beam_frequency(lambda_sq)

        # Find the computed bending mode closest to the target, because the axial modes of the
        # member are interleaved with the bending ones
        actual = min(results.frequencies, key=lambda f: abs(f - expected))

        assert actual == pytest.approx(expected, rel=5e-3), (
            f'{support} mode {mode + 1}: expected {expected:.4f} Hz, got {actual:.4f} Hz'
        )


def test_sdof_spring_mass_oscillator():
    """A mass on a spring, the simplest oracle there is: f = sqrt(k/m)/(2*pi)."""

    stiffness = 5000.0   # kN/m
    mass = 2.0           # tonne

    model = FEModel3D()

    # The spring needs a length so that its axis is defined. Placing it along Y and freeing only DY
    # at the loose end leaves a single degree of freedom acting along the spring.
    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', 0, -1.0, 0)
    model.def_support('N1', True, True, True, True, True, True)
    model.def_support('N2', True, False, True, True, True, True)
    model.add_spring('Spring', 'N1', 'N2', stiffness)
    model.add_node_mass('N2', mass)
    model.add_load_combo('Mass', {})

    results = model.analyze_modal(num_modes=1, mass_combo_name='Mass', gravity=GRAVITY,
                                  elements_per_member=1)

    expected = np.sqrt(stiffness/mass)/(2*np.pi)

    assert results.frequencies[0] == pytest.approx(expected, rel=1e-6)


def test_point_mass_at_midspan_matches_hand_sdof():
    """A dominant point mass at midspan reduces to a hand SDOF with k = 48EI/L**3."""

    mass = 1.0  # tonne, large enough that the beam's own mass barely matters

    model = beam_model('simply supported', self_weight=False)
    model.add_node('Mid', L/2, 0, 0)
    model.add_node_mass('Mid', mass)

    results = modal(model, num_modes=1)

    expected = np.sqrt(48*E*I/L**3/mass)/(2*np.pi)

    assert results.frequencies[0] == pytest.approx(expected, rel=1e-3)


def test_tip_mass_on_cantilever_matches_hand_sdof():
    """A dominant tip mass on a cantilever reduces to a hand SDOF with k = 3EI/L**3."""

    mass = 2.0  # tonne

    model = beam_model('cantilever', self_weight=False)
    model.add_node_mass('N1', mass)

    results = modal(model, num_modes=1)

    expected = np.sqrt(3*E*I/L**3/mass)/(2*np.pi)

    assert results.frequencies[0] == pytest.approx(expected, rel=1e-3)


def test_dunkerley_estimate_for_beam_with_point_mass():
    """A beam carrying a point mass falls between its two limiting cases, close to Dunkerley's.

    Dunkerley's approximation adds flexibilities, `1/f**2 = 1/f_beam**2 + 1/f_mass**2`, and
    slightly underestimates the true frequency. This pins the combined case against both limits
    and against that estimate.
    """

    mass = 0.5  # tonne, comparable to the beam's own 0.48 t

    beam_only = modal(beam_model('simply supported'), num_modes=1).frequencies[0]

    model = beam_model('simply supported')
    model.add_node('Mid', L/2, 0, 0)
    model.add_node_mass('Mid', mass)
    combined = modal(model, num_modes=1).frequencies[0]

    mass_only = np.sqrt(48*E*I/L**3/mass)/(2*np.pi)

    dunkerley = 1/np.sqrt(1/beam_only**2 + 1/mass_only**2)

    # Adding mass can only lower the frequency, and never below the mass acting alone
    assert combined < beam_only
    assert combined > dunkerley*0.98

    # Dunkerley underestimates, so the true answer sits just above it
    assert dunkerley <= combined*1.001


def test_axial_modes_of_a_fixed_free_bar():
    """Longitudinal modes of a bar fixed at one end: f_n = (2n-1)/(4L)*sqrt(EA/mu)."""

    # A stubby bar so the axial modes sit among the lower ones and are easy to identify
    length = 2.0
    model = beam_model('cantilever', length=length)

    # Axial displacement is interpolated linearly rather than cubically, so the axial modes need a
    # finer mesh than the bending modes to reach the same accuracy
    results = modal(model, num_modes=12, elements_per_member=16)

    # Mode 1 converges quickly; mode 2 is one harmonic further up and lags behind it
    for n, tolerance in ((1, 2e-3), (2, 1e-2)):

        expected = (2*n - 1)/(4*length)*np.sqrt(E*A/MU)
        actual = min(results.frequencies, key=lambda f: abs(f - expected))

        assert actual == pytest.approx(expected, rel=tolerance), (
            f'Axial mode {n}: expected {expected:.4f} Hz, got {actual:.4f} Hz'
        )


def test_portal_frame_sway_matches_hand_sdof():
    """A portal frame's sway mode against a hand SDOF idealisation.

    The beam is made stiff and heavy and the columns light, so the frame approaches the textbook
    case of a rigid mass on two fixed-base columns, `k = 2*12EI/h**3`.
    """

    height = 3.0
    span = 6.0
    mass = 10.0  # tonne, lumped at the two eaves

    model = FEModel3D()
    model.add_node('B1', 0, 0, 0)
    model.add_node('B2', span, 0, 0)
    model.add_node('E1', 0, height, 0)
    model.add_node('E2', span, height, 0)
    model.def_support('B1', True, True, True, True, True, True)
    model.def_support('B2', True, True, True, True, True, True)

    model.add_material('Steel', E, G, NU, 0.0)
    model.add_section('Column', A, I, I, J)

    # A beam stiff enough to act as the rigid link the hand idealisation assumes
    model.add_section('Beam', A*100, I*1000, I*1000, J*1000)

    model.add_member('C1', 'B1', 'E1', 'Steel', 'Column')
    model.add_member('C2', 'B2', 'E2', 'Steel', 'Column')
    model.add_member('Beam', 'E1', 'E2', 'Steel', 'Beam')

    model.add_node_mass('E1', mass/2)
    model.add_node_mass('E2', mass/2)
    model.add_load_combo('Mass', {})

    results = modal(model, num_modes=1)

    expected = np.sqrt(2*12*E*I/height**3/mass)/(2*np.pi)

    # A 5% band: the columns are not weightless in the model and the beam is stiff but not rigid
    assert results.frequencies[0] == pytest.approx(expected, rel=0.05)


# ---------------------------------------------------------------------------------------------
# Structural identities
# ---------------------------------------------------------------------------------------------

def test_element_consistent_mass_matrix_coefficients():
    """The element mass matrix must reproduce the classical rho*A*L/420 coefficients by hand.

    Bending coefficients 156, 22L, 54, -13L, 4L**2, -3L**2 and axial coefficients 140, 70, per
    Cook, Malkus, Plesha & Witt, *Concepts and Applications of Finite Element Analysis*.
    """

    from Pynite.Analysis import _prepare_model

    model = beam_model('simply supported')
    _prepare_model(model)

    member = list(list(model.members.values())[0].sub_members.values())[0]
    length = member.L()

    total_mass = RHO*A*length/GRAVITY
    scale = total_mass/420

    expected = scale*np.array([
        [140, 0,    0,    0,          0,        0,        70,  0,    0,    0,          0,        0       ],
        [0,   156,  0,    0,          0,        22*length, 0,  54,   0,    0,          0,       -13*length],
        [0,   0,    156,  0,         -22*length, 0,       0,   0,    54,   0,          13*length, 0      ],
        [0,   0,    0,    140*J/A,    0,        0,        0,   0,    0,    70*J/A,     0,        0       ],
        [0,   0,   -22*length, 0,     4*length**2, 0,     0,   0,   -13*length, 0,    -3*length**2, 0    ],
        [0,   22*length, 0, 0,        0,        4*length**2, 0, 13*length, 0, 0,       0,       -3*length**2],
        [70,  0,    0,    0,          0,        0,        140, 0,    0,    0,          0,        0       ],
        [0,   54,   0,    0,          0,        13*length, 0,  156,  0,    0,          0,       -22*length],
        [0,   0,    54,   0,         -13*length, 0,       0,   0,    156,  0,          22*length, 0      ],
        [0,   0,    0,    70*J/A,     0,        0,        0,   0,    0,    140*J/A,    0,        0       ],
        [0,   0,    13*length, 0,    -3*length**2, 0,     0,   0,    22*length, 0,     4*length**2, 0    ],
        [0,  -13*length, 0, 0,        0,       -3*length**2, 0, -22*length, 0, 0,      0,        4*length**2],
    ])

    actual = member.consistent_m('Mass', GRAVITY)

    assert np.allclose(actual, expected, rtol=1e-12, atol=1e-15)


def test_shape_function_matrix_partitions_unity():
    """A rigid unit translation must move every interior point by exactly one.

    This is what makes the load-derived consistent mass matrix reproduce the total mass exactly, so
    it is worth pinning separately from the totals below.
    """

    from Pynite.Analysis import _prepare_model

    model = beam_model('simply supported')
    _prepare_model(model)
    member = list(list(model.members.values())[0].sub_members.values())[0]

    for xi in (0.0, 0.1, 0.5, 0.75, 1.0):

        N = member._shape_matrix(xi)

        for direction in range(3):

            # A unit rigid translation of both end nodes in one direction
            q = np.zeros(12)
            q[direction] = 1.0
            q[direction + 6] = 1.0

            assert (N @ q)[direction] == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize('elements_per_member', [1, 2, 5, 8])
def test_total_assembled_mass_is_exact_for_self_mass(elements_per_member):
    """Total assembled mass must equal the model's mass, at any level of subdivision.

    The total is `r.T @ M @ r` for a rigid unit translation, not the sum of the diagonal: a
    consistent mass matrix shares mass between coupled DOFs, so its diagonal sums to only 312/420
    of the total in the transverse directions.
    """

    results = modal(beam_model('simply supported'), num_modes=2,
                    elements_per_member=elements_per_member)

    expected = MU*L

    for direction in 'XYZ':
        assert results.total_mass[direction] == pytest.approx(expected, rel=1e-9)

    # The per-node shares must add back up to the total
    assert sum(results.mass_per_node.values()) == pytest.approx(expected, rel=1e-9)


def test_total_assembled_mass_is_exact_for_load_derived_mass():
    """A distributed load converted to mass must assemble to exactly that mass.

    Load-derived mass is the path the serviceability mass combination uses, so its total matters
    more than any other. It is integrated over the loaded length rather than lumped at an estimated
    centroid, which is what makes the total exact rather than approximate.
    """

    model = beam_model('simply supported', self_weight=False)

    # A load chosen to carry exactly the same mass as the beam's own self-weight
    w = -MU*GRAVITY
    model.add_member_dist_load('M1', 'FY', w, w, case='Q')
    model.load_combos['Mass'].factors['Q'] = 1.0

    results = modal(model, num_modes=2)

    for direction in 'XYZ':
        assert results.total_mass[direction] == pytest.approx(MU*L, rel=1e-9)


def test_total_assembled_mass_is_exact_for_combined_sources():
    """Self-mass, load-derived mass and explicit nodal mass must add up without double counting."""

    model = beam_model('simply supported')

    w = -MU*GRAVITY
    model.add_member_dist_load('M1', 'FY', w, w, case='Q')
    model.load_combos['Mass'].factors['Q'] = 1.0

    model.add_node('Mid', L/2, 0, 0)
    model.add_node_mass('Mid', 0.25)

    results = modal(model, num_modes=2)

    # Self-weight plus an equal load-derived mass plus the nodal mass
    expected = MU*L + MU*L + 0.25

    assert results.total_mass['Y'] == pytest.approx(expected, rel=1e-9)


def test_load_derived_mass_is_partly_lost_when_lumped():
    """The lumped formulation still totals correctly, but brackets the frequency from below.

    Consistent mass converges on the exact frequency from above and lumped mass from below, so a
    coarse mesh straddling the closed-form answer is evidence that both paths are right.
    """

    exact = beam_frequency(np.pi**2)

    def f1(formulation):
        model = beam_model('simply supported', self_weight=False)
        w = -MU*GRAVITY
        model.add_member_dist_load('M1', 'FY', w, w, case='Q')
        model.load_combos['Mass'].factors['Q'] = 1.0

        return modal(model, num_modes=1, elements_per_member=2,
                     mass_formulation=formulation).frequencies[0]

    consistent = f1('consistent')
    lumped = f1('lumped')

    assert lumped < exact < consistent
    assert consistent == pytest.approx(exact, rel=0.01)
    assert lumped == pytest.approx(exact, rel=0.01)


def test_load_derived_mass_is_never_negative():
    """Every DOF must carry non-negative mass, whichever formulation is used.

    A negative diagonal makes the mass matrix indefinite, which surfaces as an eigensolver
    convergence failure rather than as the modelling error it is.
    """

    from Pynite.Analysis import _prepare_model

    for formulation in ('consistent', 'lumped'):
        for elements in (1, 2, 4):

            model = beam_model('simply supported', self_weight=False)
            w = -MU*GRAVITY
            model.add_member_dist_load('M1', 'FY', w, w, case='Q')
            model.load_combos['Mass'].factors['Q'] = 1.0

            mesh = model._modal_mesh_copy(elements, 'XY')
            _prepare_model(mesh)
            M = mesh.M('Mass', 'Y', GRAVITY, sparse=False, mass_formulation=formulation)

            assert np.all(np.diag(M) >= 0), (
                f'{formulation} mass with {elements} element(s) produced a negative diagonal'
            )


def test_load_derived_mass_is_sign_insensitive():
    """An upward load carries mass just as a downward one does."""

    def total(sign):
        model = beam_model('simply supported', self_weight=False)
        w = sign*MU*GRAVITY
        model.add_member_dist_load('M1', 'FY', w, w, case='Q')
        model.load_combos['Mass'].factors['Q'] = 1.0

        return modal(model, num_modes=1).total_mass['Y']

    assert total(+1) == pytest.approx(total(-1), rel=1e-12)


def test_local_direction_loads_convert_to_mass():
    """Loads given in member-local directions must convert without error.

    The member here runs along global X, so its local y-axis is the global Y-axis and a local `Fy`
    load has to produce the same mass as a global `FY` one.
    """

    def total(direction):
        model = beam_model('simply supported', self_weight=False)
        w = -MU*GRAVITY
        model.add_member_dist_load('M1', direction, w, w, case='Q')
        model.load_combos['Mass'].factors['Q'] = 1.0

        return modal(model, num_modes=1).total_mass['Y']

    assert total('Fy') == pytest.approx(total('FY'), rel=1e-12)


def test_mode_shapes_are_mass_orthonormal():
    """Mass normalization and orthogonality: `phi.T @ M11 @ phi` must be the identity."""

    results = modal(beam_model('simply supported'), num_modes=4)

    M = results.M.toarray()
    free = results.free_dof_indices
    M11 = M[np.ix_(free, free)]

    product = results.mode_shapes.T @ M11 @ results.mode_shapes

    assert np.allclose(product, np.eye(results.mode_count), atol=1e-8)


def test_rayleigh_quotient_returns_the_eigenvalues():
    """With mass-normalized mode shapes, `phi.T @ K11 @ phi` must return `omega**2`."""

    from Pynite.Analysis import _partition, _partition_D

    model = beam_model('simply supported')
    results = modal(model, num_modes=4)

    # Rebuild the free-DOF stiffness of the mesh the solver actually used
    mesh = model._modal_mesh_copy(8, 'XY')
    from Pynite.Analysis import _prepare_model, _set_force_timoshenko
    _prepare_model(mesh, 4)
    _set_force_timoshenko(mesh, True)
    D1, D2i, _ = _partition_D(mesh)
    K11 = _partition(mesh, mesh.Ke('Mass', sparse=True).tocsr(), D1, D2i)[0].toarray()

    quotient = np.diag(results.mode_shapes.T @ K11 @ results.mode_shapes)

    assert np.allclose(quotient, results.eigenvalues, rtol=1e-8)
    assert np.allclose(np.sqrt(quotient)/(2*np.pi), results.frequencies, rtol=1e-8)


# ---------------------------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------------------------

def test_frequencies_are_invariant_to_node_numbering():
    """Renumbering the model must not move the frequencies."""

    forward = modal(beam_model('simply supported', n_user_nodes=9), num_modes=3).frequencies

    # The same beam with its nodes defined in the opposite order
    reverse_model = FEModel3D()
    for i in reversed(range(9)):
        reverse_model.add_node(f'N{i}', L*i/8, 0, 0)
    reverse_model.def_support('N0', True, True, False, False, False, False)
    reverse_model.def_support('N8', False, True, False, False, False, False)
    reverse_model.add_material('Steel', E, G, NU, RHO)
    reverse_model.add_section('Section', A, I, I, J)
    reverse_model.add_member('M1', 'N0', 'N8', 'Steel', 'Section')
    reverse_model.add_member_self_weight('FY', -1.0, 'SW')
    reverse_model.add_load_combo('Mass', {'SW': 1.0})

    reverse = modal(reverse_model, num_modes=3).frequencies

    assert np.allclose(forward, reverse, rtol=1e-9)


def test_mode_shape_signs_are_invariant_to_node_numbering():
    """Mode shape signs must not flip when the model is renumbered.

    Without a deterministic convention the eigensolver's sign is arbitrary, which flips animations
    between runs and makes differential tests fail for no real reason. The antisymmetric modes are
    the ones that matter here: they have two equal and opposite peaks, so their sign is decided
    entirely by the tie-break, and a tie-break on DOF index would not survive renumbering.
    """

    forward_model = beam_model('simply supported', n_user_nodes=9)
    modal(forward_model, num_modes=3)

    reverse_model = FEModel3D()
    for i in reversed(range(9)):
        reverse_model.add_node(f'N{i}', L*i/8, 0, 0)
    reverse_model.def_support('N0', True, True, False, False, False, False)
    reverse_model.def_support('N8', False, True, False, False, False, False)
    reverse_model.add_material('Steel', E, G, NU, RHO)
    reverse_model.add_section('Section', A, I, I, J)
    reverse_model.add_member('M1', 'N0', 'N8', 'Steel', 'Section')
    reverse_model.add_member_self_weight('FY', -1.0, 'SW')
    reverse_model.add_load_combo('Mass', {'SW': 1.0})
    modal(reverse_model, num_modes=3)

    for mode in range(1, 4):

        a = np.array([forward_model.nodes[f'N{i}'].DY[f'Mode {mode}'] for i in range(9)])
        b = np.array([reverse_model.nodes[f'N{i}'].DY[f'Mode {mode}'] for i in range(9)])

        # The shapes must be the same to within numerical noise, signs included
        assert np.allclose(a, b, atol=1e-8), f'Mode {mode} signs or shape differ after renumbering'


def test_repeated_runs_are_identical():
    """Two runs of the same model must agree bit for bit, frequencies and mode shapes alike."""

    first = modal(beam_model('simply supported'), num_modes=3)
    second = modal(beam_model('simply supported'), num_modes=3)

    assert np.array_equal(first.frequencies, second.frequencies)
    assert np.array_equal(first.mode_shapes, second.mode_shapes)


# ---------------------------------------------------------------------------------------------
# Convergence study -- this is where the `elements_per_member` default comes from
# ---------------------------------------------------------------------------------------------

# Measured error in percent for a uniform simply supported beam, by elements per member. A member
# is one element between its end nodes, so without subdivision even mode 3 is out by a factor of
# four. The default of 8 keeps modes 1 to 3 within 0.15%.
CONVERGENCE = {
    1: (10.99, 27.16, 305.29),
    2: (0.39, 10.99, 23.99),
    3: (0.08, 1.18, 10.99),
    4: (0.03, 0.39, 1.83),
    6: (0.01, 0.08, 0.39),
    8: (0.00, 0.03, 0.13),
    12: (0.00, 0.01, 0.03),
    16: (0.00, 0.00, 0.01),
}


@pytest.mark.parametrize('elements_per_member', sorted(CONVERGENCE))
def test_convergence_study(elements_per_member):
    """Pins the measured convergence of modes 1 to 3 against elements per member.

    The table this checks is the evidence for the `elements_per_member` default and belongs in the
    documentation of any analysis that relies on it. Convergence is monotonic and from above, as it
    must be for a consistent mass matrix.
    """

    expected_errors = CONVERGENCE[elements_per_member]

    results = modal(beam_model('simply supported'), num_modes=3,
                    elements_per_member=elements_per_member)

    for mode, lambda_sq in enumerate(LAMBDA_SQ['simply supported']):

        exact = beam_frequency(lambda_sq)
        error = 100*(results.frequencies[mode] - exact)/exact

        # Convergence is from above for a consistent mass matrix
        assert error >= -0.01, f'Mode {mode + 1} converged from below: {error:+.3f}%'

        assert error == pytest.approx(expected_errors[mode], abs=0.02), (
            f'Mode {mode + 1} with {elements_per_member} element(s): expected '
            f'{expected_errors[mode]:+.2f}%, measured {error:+.3f}%'
        )


def test_default_subdivision_keeps_the_first_three_modes_accurate():
    """The shipped default has to be good enough to rely on without thinking about it."""

    results = modal(beam_model('simply supported'), num_modes=3)

    assert results.diagnostics.elements_per_member == 8

    for mode, lambda_sq in enumerate(LAMBDA_SQ['simply supported']):

        exact = beam_frequency(lambda_sq)
        error = abs(100*(results.frequencies[mode] - exact)/exact)

        assert error < 0.15, f'Mode {mode + 1} is out by {error:.3f}% at the default subdivision'


# ---------------------------------------------------------------------------------------------
# Units, robustness and isolation
# ---------------------------------------------------------------------------------------------

def test_gravity_default_is_a_trap():
    """Documents the consequence of leaving `gravity` at its default of 1.0.

    Mass is derived from loads by dividing by gravity, so the default overstates the mass by a
    factor of g and understates every frequency by sqrt(g). This is the single most likely
    silent-wrong-answer failure mode in the feature, and it is a trap rather than a bug: the fork
    is unit-agnostic and cannot know what g is. The caller is responsible for always passing it.
    """

    correct = modal(beam_model('simply supported'), num_modes=1).frequencies[0]

    wrong = beam_model('simply supported').analyze_modal(
        num_modes=1, mass_combo_name='Mass', plane='XY'
    ).frequencies[0]

    assert wrong == pytest.approx(correct/np.sqrt(GRAVITY), rel=1e-6)


def test_rho_is_a_weight_density():
    """Pins the meaning of `rho`, because getting it wrong is a factor-of-g error.

    `add_member_self_weight` reads `rho*A` as a force per unit length, so `rho` is a weight density
    and the assembled mass is `rho*A*L/g`. A caller who supplies a mass density instead lands a
    factor of g out.
    """

    results = modal(beam_model('simply supported'), num_modes=1)

    assert results.total_mass['Y'] == pytest.approx(RHO*A*L/GRAVITY, rel=1e-9)


def test_massless_model_is_rejected_with_a_useful_message():
    """A model with no mass is a modelling error, not a solver failure."""

    model = beam_model('simply supported', self_weight=False)

    with pytest.raises(ModalModelError, match='no mass'):
        modal(model, num_modes=1)


def test_fully_supported_model_is_rejected():
    """A model with no free DOFs has nothing to vibrate."""

    model = beam_model('simply supported')
    for name in list(model.nodes):
        model.def_support(name, True, True, True, True, True, True)

    with pytest.raises(ModalModelError, match='supported'):
        modal(model, num_modes=1, elements_per_member=1)


def test_more_modes_requested_than_exist_is_flagged_not_crashed():
    """Asking for more modes than the model has must degrade gracefully and say so."""

    results = modal(beam_model('simply supported'), num_modes=200, elements_per_member=1)

    assert results.diagnostics.requested_modes == 200
    assert results.diagnostics.converged_modes < 200
    assert results.diagnostics.truncated is True
    assert results.mode_count == results.diagnostics.converged_modes


def test_rigid_body_modes_are_filtered():
    """A free-floating model's rigid-body modes must be filtered, not reported as frequencies."""

    model = beam_model('free-free')

    # A single element cannot represent the first elastic mode of a free-free beam, so the default
    # subdivision is what makes this comparison meaningful
    results = modal(model, num_modes=8, check_stability=False)

    assert results.diagnostics.filtered_modes, 'No rigid-body modes were filtered'
    assert np.all(results.frequencies > 1.0)

    # The first elastic mode of a free-free beam shares the fixed-fixed coefficient
    expected = beam_frequency(LAMBDA_SQ['free-free'][0])
    assert results.frequencies[0] == pytest.approx(expected, rel=5e-3)


def test_tension_only_members_require_an_explicit_state():
    """A load-dependent stiffness has no single linear state, so it must not be guessed."""

    model = beam_model('simply supported')
    model.add_node('Mid', L/2, 0, 0)
    model.add_node('Anchor', L/2, -2.0, 0)
    model.def_support('Anchor', True, True, True, True, True, True)
    model.add_member('Tie', 'Anchor', 'Mid', 'Steel', 'Section', tension_only=True)

    with pytest.raises(ModalModelError, match='tension-only'):
        modal(model, num_modes=1)

    # An unsupported declaration is refused rather than quietly reinterpreted
    with pytest.raises(ModalModelError, match='all_active'):
        modal(model, num_modes=1, linear_state='state from combination Q')

    # The supported declaration goes through
    results = modal(model, num_modes=1, linear_state='all_active')
    assert results.frequencies[0] > 0


def test_invalid_arguments_are_rejected():
    """Bad arguments must fail immediately rather than part way through an assembly."""

    model = beam_model('simply supported')

    with pytest.raises(ValueError, match='num_modes'):
        modal(model, num_modes=0)

    with pytest.raises(ValueError, match='elements_per_member'):
        modal(model, num_modes=1, elements_per_member=0)

    with pytest.raises(ValueError, match='mass_formulation'):
        modal(model, num_modes=1, mass_formulation='consistant')

    with pytest.raises(ValueError, match='plane'):
        model.analyze_modal(num_modes=1, mass_combo_name='Mass', gravity=GRAVITY, plane='XW')


def test_negative_nodal_mass_is_rejected():
    """A negative mass is never physical and would make the mass matrix indefinite."""

    model = beam_model('simply supported')

    with pytest.raises(ValueError, match='negative'):
        model.add_node_mass('N0', -1.0)

    with pytest.raises(ValueError, match='negative'):
        model.add_node_mass('N0', 1.0, IY=-1.0)


def test_modal_analysis_leaves_the_users_model_alone():
    """Subdivision is an analysis-time refinement and must not be visible afterwards.

    The temporary nodes and the out-of-plane restraints both live on an internal copy, so the
    caller's model keeps the geometry and supports it was given, and the static results computed
    from it stay valid.
    """

    model = beam_model('simply supported', brace_out_of_plane=True)
    model.add_member_dist_load('M1', 'FY', -10.0, -10.0, case='Q')
    model.add_load_combo('Static', {'Q': 1.0})
    model.add_node('Mid', L/2, 0, 0)
    model.nodes['Mid'].support_DZ = True
    model.nodes['Mid'].support_RX = True
    model.nodes['Mid'].support_RY = True

    model.analyze_linear()
    static_before = {name: node.DY['Static'] for name, node in model.nodes.items()}
    nodes_before = set(model.nodes)
    plane_free_before = model.nodes['Mid'].support_DZ

    modal(model, num_modes=3)

    # Geometry and supports are untouched
    assert set(model.nodes) == nodes_before
    assert not any(name.startswith('_modal_') for name in model.nodes)
    assert model.nodes['Mid'].support_DZ == plane_free_before

    # The static results survive the modal run
    for name, node in model.nodes.items():
        assert node.DY['Static'] == static_before[name]

    # And a re-run of the static analysis reproduces them exactly
    model.analyze_linear()
    for name, node in model.nodes.items():
        assert node.DY['Static'] == static_before[name]


def test_static_results_are_unaffected_by_the_modal_feature():
    """A model that never asks for modes must solve statically exactly as it always did."""

    def deflection():
        model = beam_model('simply supported', brace_out_of_plane=True)
        model.add_member_dist_load('M1', 'FY', -10.0, -10.0, case='Q')
        model.add_load_combo('Static', {'Q': 1.0})
        model.add_node('Mid', L/2, 0, 0)
        model.nodes['Mid'].support_DZ = True
        model.nodes['Mid'].support_RX = True
        model.nodes['Mid'].support_RY = True
        model.analyze_linear()

        return model.nodes['Mid'].DY['Static']

    # The closed-form midspan deflection of a uniformly loaded simply supported beam
    expected = -5*10.0*L**4/(384*E*I)

    assert deflection() == pytest.approx(expected, rel=1e-3)

    # And it is reproducible
    assert deflection() == deflection()


def test_mode_shapes_are_reachable_from_the_result_and_the_model():
    """Mode shapes must be available both on the nodes and through the structured result."""

    model = beam_model('simply supported', n_user_nodes=3)
    results = modal(model, num_modes=2)

    # Through the model, under the `Mode n` load combinations
    assert 'Mode 1' in model.load_combos
    assert 'Mode 2' in model.load_combos
    midspan = model.nodes['N1'].DY['Mode 1']
    assert abs(midspan) > 0

    # Through the result, including the subdivision nodes the model never sees
    shape = results.mode_shape(1)
    assert ('N1', 'DY') in shape
    assert shape[('N1', 'DY')] == pytest.approx(midspan, rel=1e-9)
    assert any(name.startswith('_modal_') for name, _ in results.dof_map)

    with pytest.raises(IndexError):
        results.mode_shape(3)


def test_modal_combos_are_tagged_so_they_can_be_excluded():
    """The `Mode n` combos must carry a tag, so envelopes and reports can filter them out.

    They are real entries in `load_combos` holding mode shapes rather than displacements, so any
    code that walks every combination will pick them up unless it filters on the tag.
    """

    model = beam_model('simply supported', brace_out_of_plane=True)
    modal(model, num_modes=3)

    modal_combos = [name for name, combo in model.load_combos.items()
                    if combo.combo_tags is not None and 'modal' in combo.combo_tags]

    assert sorted(modal_combos) == ['Mode 1', 'Mode 2', 'Mode 3']

    # And they are cleared by the next analysis rather than accumulating
    model.analyze_linear()

    assert not [name for name, combo in model.load_combos.items()
                if combo.combo_tags is not None and 'modal' in combo.combo_tags]


def test_diagnostics_report_what_the_solver_did():
    """The diagnostics have to be complete enough to document an analysis from."""

    results = modal(beam_model('simply supported'), num_modes=3, mass_formulation='consistent',
                    elements_per_member=6)

    diagnostics = results.diagnostics

    assert diagnostics.requested_modes == 3
    assert diagnostics.converged_modes == 3
    assert diagnostics.truncated is False
    assert diagnostics.mass_formulation == 'consistent'
    assert diagnostics.elements_per_member == 6
    assert diagnostics.solver in ('sparse-shift-invert', 'dense')
    assert 'consistent mass' in diagnostics.summary()


def test_planar_analysis_excludes_out_of_plane_modes():
    """Declaring a plane must keep lateral and torsional modes out of the results.

    A cantilever with `Iy == Iz` has identical boundary conditions in both bending planes, so every
    bending mode is a degenerate pair. In three dimensions both halves of each pair come back and
    half the requested modes are spent on the lateral one; in plane, the same request buys four
    distinct in-plane modes.
    """

    planar = modal(beam_model('cantilever'), num_modes=4).frequencies

    spatial = beam_model('cantilever').analyze_modal(
        num_modes=4, mass_combo_name='Mass', gravity=GRAVITY
    ).frequencies

    # In plane, the four lowest modes are four distinct modes
    assert len(set(np.round(planar, 6))) == 4

    # In three dimensions they come in degenerate pairs, so the same four slots hold fewer
    assert len(set(np.round(spatial, 6))) < 4


def test_out_of_plane_modes_interleave_with_the_in_plane_ones():
    """Lateral modes are not merely extra, they land *between* the in-plane modes.

    This is what makes them a correctness problem rather than a nuisance: mode 2 of a
    three-dimensional run of a laterally braced beam is a lateral mode, so anything that reads "the
    second mode" gets an answer about the wrong plane. The interior nodes that subdivision creates
    are the ones left free here, which is why the caller cannot fix this from outside.
    """

    braced = beam_model('simply supported', brace_out_of_plane=True)
    planar = modal(braced, num_modes=4).frequencies

    spatial = beam_model('simply supported', brace_out_of_plane=True).analyze_modal(
        num_modes=4, mass_combo_name='Mass', gravity=GRAVITY
    ).frequencies

    # Both runs agree on the fundamental, which is an in-plane mode either way
    assert spatial[0] == pytest.approx(planar[0], rel=1e-6)

    # But the second mode of the 3D run is a lateral mode that the planar run never reports
    assert spatial[1] < planar[1]
    assert not np.any(np.isclose(planar, spatial[1], rtol=1e-6))


# ---------------------------------------------------------------------------------------------
# Sign conventions for self-weight
# ---------------------------------------------------------------------------------------------

def test_density_sign_is_ignored_when_assembling_mass():
    """A negative weight density must carry exactly the same mass as a positive one.

    There are two ways to express the direction of self-weight: a negative density with a positive
    load factor, or a positive density with a negative load factor. Both reach the same load, and
    both must reach the same mass, because mass is not a signed quantity. Downstream code relies on
    this to keep one self-weight convention across static and modal analysis, so it is a contract
    rather than an implementation detail -- a refactor that concluded "density is positive, drop the
    `abs`" would produce negative element masses and an indefinite mass matrix.
    """

    def run(rho, self_weight_factor):
        model = beam_model('simply supported', self_weight=False)
        model.materials['Steel'].rho = rho
        model.add_member_self_weight('FY', self_weight_factor, 'SW')
        model.load_combos['Mass'].factors['SW'] = 1.0

        return modal(model, num_modes=3)

    negative_density = run(-RHO, +1.0)
    positive_density = run(+RHO, -1.0)

    assert negative_density.total_mass['Y'] == pytest.approx(MU*L, rel=1e-12)
    assert np.array_equal(negative_density.frequencies, positive_density.frequencies)
    assert np.array_equal(negative_density.mode_shapes, positive_density.mode_shapes)


def test_all_four_sign_combinations_carry_the_same_mass():
    """Neither sign matters to the mass, only the product's magnitude."""

    totals = set()

    for rho in (-RHO, +RHO):
        for self_weight_factor in (-1.0, +1.0):

            model = beam_model('simply supported', self_weight=False)
            model.materials['Steel'].rho = rho
            model.add_member_self_weight('FY', self_weight_factor, 'SW')
            model.load_combos['Mass'].factors['SW'] = 1.0

            totals.add(round(modal(model, num_modes=1).total_mass['Y'], 12))

    assert totals == {round(MU*L, 12)}


def test_self_weight_factor_scales_the_mass():
    """A self-weight factor raised to account for connections must raise the mass too.

    `add_member_self_weight(direction, factor)` applies `factor*rho*A` as the load, so a factor of
    1.15 means the member weighs 15% more. That extra weight has to appear dynamically as well as
    statically, or a model tuned for connection weight would report frequencies that are too high.
    """

    for self_weight_factor in (1.0, 1.15, 2.0):

        model = beam_model('simply supported', self_weight=False)
        model.add_member_self_weight('FY', -self_weight_factor, 'SW')
        model.load_combos['Mass'].factors['SW'] = 1.0

        results = modal(model, num_modes=1)

        assert results.total_mass['Y'] == pytest.approx(MU*L*self_weight_factor, rel=1e-9), (
            f'A self-weight factor of {self_weight_factor} did not scale the mass'
        )


def test_opposing_self_weight_cases_add_rather_than_cancel():
    """Two self-weight contributions of opposing sign must add, not cancel.

    Mass has to be accumulated as magnitudes. Summing signed contributions and taking the absolute
    value only at the end lets a pair of opposing self-weight cases -- from mixed density sign
    conventions across materials, or a negative combination factor -- silently zero the self-mass.
    """

    model = beam_model('simply supported', self_weight=False)
    model.add_member_self_weight('FY', -1.0, 'D1')
    model.add_member_self_weight('FY', +1.0, 'D2')
    model.load_combos['Mass'].factors['D1'] = 1.0
    model.load_combos['Mass'].factors['D2'] = 1.0

    results = modal(model, num_modes=1)

    # Two self-weight loads means twice the mass, not zero
    assert results.total_mass['Y'] == pytest.approx(2*MU*L, rel=1e-9)


# ---------------------------------------------------------------------------------------------
# Participation
# ---------------------------------------------------------------------------------------------

def test_participating_mass_is_the_free_dof_mass():
    """Participating mass must be the mass on the free DOFs, below the assembled total.

    Mass held on a restrained DOF cannot move in any mode, so it is unreachable no matter how many
    modes are computed. Reporting it separately is what stops a consumer from dividing effective
    modal mass by the total and concluding that participation never completes.
    """

    results = modal(beam_model('simply supported'), num_modes=3)

    # In plane there is mobilizable mass, but less than the total, because the supports hold some
    for direction in 'XY':
        assert 0 < results.participating_mass[direction] < results.total_mass[direction]

    # Out of plane there is none at all: a planar analysis restrains every out-of-plane DOF, so no
    # mode can mobilize mass in that direction even though the mass is still assembled
    assert results.participating_mass['Z'] == 0.0
    assert results.total_mass['Z'] > 0.0


def test_participating_mass_equals_the_hand_calculation():
    """The exposed value must equal the `r.T @ M11 @ r` a consumer would compute by hand."""

    results = modal(beam_model('simply supported'), num_modes=3)

    M11 = results.M.toarray()[np.ix_(results.free_dof_indices, results.free_dof_indices)]

    for offset, direction in enumerate('XYZ'):

        r = np.zeros(results.M.shape[0])
        r[offset::6] = 1.0
        r11 = r[results.free_dof_indices]

        assert results.participating_mass[direction] == pytest.approx(float(r11 @ (M11 @ r11)),
                                                                     rel=1e-12)


def test_effective_mass_sums_to_the_participating_mass():
    """Summed over a complete set of modes, effective modal mass equals the participating mass.

    This is the completeness identity that participation post-processing rests on: the mode shapes
    span the free DOFs. A small model is used so every mode can be solved for, which makes the
    identity exact rather than approximate.
    """

    results = modal(beam_model('simply supported'), num_modes=400, elements_per_member=2)

    assert results.mode_count == len(results.free_dof_indices)

    for direction in 'XY':

        effective = results.effective_mass(direction)

        assert effective.sum() == pytest.approx(results.participating_mass[direction], rel=1e-9)

        # And the normalized form accounts for everything
        assert results.mass_participation(direction).sum() == pytest.approx(1.0, rel=1e-9)

    # A direction with no mobilizable mass reports zero participation rather than dividing by zero
    assert results.participating_mass['Z'] == 0.0
    assert np.all(results.effective_mass('Z') == 0.0)
    assert np.all(results.mass_participation('Z') == 0.0)


def test_participation_api_matches_the_hand_calculation():
    """The convenience methods must agree with the matrix algebra they replace."""

    results = modal(beam_model('simply supported'), num_modes=4)

    M11 = results.M.toarray()[np.ix_(results.free_dof_indices, results.free_dof_indices)]
    r = np.zeros(results.M.shape[0])
    r[1::6] = 1.0
    r11 = r[results.free_dof_indices]

    gamma = results.mode_shapes.T @ (M11 @ r11)

    assert np.allclose(results.participation_factors('Y'), gamma, rtol=1e-9)
    assert np.allclose(results.effective_mass('Y'), gamma**2, rtol=1e-9)
    assert np.allclose(results.mass_participation('Y'),
                       gamma**2/results.participating_mass['Y'], rtol=1e-9)


def test_mass_participation_identifies_the_governing_mode():
    """The governing mode is the lowest one with real participation, not necessarily mode 1.

    A simply supported beam's antisymmetric modes participate in nothing, which is exactly why "the
    frequency of the beam" has to be chosen on effective mass rather than by taking the first mode.
    """

    results = modal(beam_model('simply supported'), num_modes=6)

    participation = results.mass_participation('Y')

    # Mode 1 dominates, approaching the continuum value of 8/pi**2 of the total mass
    assert participation[0] > 0.8

    # Mode 2 is antisymmetric and participates in nothing
    assert participation[1] < 1e-6

    # Picking the governing mode the way a consumer would
    governing = int(np.argmax(participation > 0.05))
    assert governing == 0
    assert results.frequencies[governing] == pytest.approx(beam_frequency(np.pi**2), rel=2e-3)


def test_participation_rejects_a_bad_direction():
    """A mistyped direction must fail loudly rather than return zeros."""

    results = modal(beam_model('simply supported'), num_modes=1)

    with pytest.raises(ValueError, match='direction'):
        results.effective_mass('vertical')


def test_mass_per_node_is_the_same_in_every_direction():
    """Per-node mass shares must not depend on which direction they are taken in.

    A full translational row of the consistent mass matrix sums to half the element mass whether it
    is the axial row (140 + 70) or a transverse one (156 + 54), so the shares coincide. This is why
    the reported shares need no direction of their own.
    """

    results = modal(beam_model('simply supported', n_user_nodes=5), num_modes=2)

    M = results.M.toarray()

    for offset in range(3):

        r = np.zeros(M.shape[0])
        r[offset::6] = 1.0
        Mr = M @ r

        by_direction = {name: Mr[node_id*6 + offset]
                        for name, node_id in [(n, i) for i, n in enumerate(results.mass_per_node)]}

        # Compare against the reported shares, which are taken in Y
        for name, share in results.mass_per_node.items():
            assert by_direction[name] == pytest.approx(share, rel=1e-9)


# ---------------------------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------------------------

def test_diagnostics_record_the_run_provenance():
    """Every argument that changes the answer must be readable back off the result.

    Consumers need this for two things: a cache key that cannot drift out of step with the request
    it describes, and documenting the analysis assumption, which has to state the mass basis,
    gravity, the analysis plane and the linear-state declaration.
    """

    model = beam_model('simply supported')
    results = model.analyze_modal(num_modes=3, mass_combo_name='Mass', mass_direction='Y',
                                  gravity=GRAVITY, plane='XY', mass_formulation='lumped',
                                  elements_per_member=6)

    d = results.diagnostics

    assert d.mass_combo_name == 'Mass'
    assert d.mass_direction == 'Y'
    assert d.gravity == GRAVITY
    assert d.mass_formulation == 'lumped'
    assert d.elements_per_member == 6
    assert d.plane == 'XY'
    assert d.linear_state is None

    # Enough of it makes the one-line summary for a report to stand alone
    summary = d.summary()
    assert 'Mass' in summary
    assert str(GRAVITY) in summary
    assert 'XY plane' in summary
    assert 'lumped' in summary


def test_diagnostics_record_the_linear_state_declaration():
    """A declared linear state has to be recorded, because it is an analysis assumption."""

    model = beam_model('simply supported')
    model.add_node('Mid', L/2, 0, 0)
    model.add_node('Anchor', L/2, -2.0, 0)
    model.def_support('Anchor', True, True, True, True, True, True)
    model.add_member('Tie', 'Anchor', 'Mid', 'Steel', 'Section', tension_only=True)

    results = modal(model, num_modes=1, linear_state='all_active')

    assert results.diagnostics.linear_state == 'all_active'
    assert "linear state 'all_active'" in results.diagnostics.summary()


def test_diagnostics_report_a_three_dimensional_run():
    """A run with no plane declared must say so rather than leaving the field ambiguous.

    A cantilever is used because it is fully restrained at its base and so stable in three
    dimensions. The simply supported fixture is not: with nothing holding it out of plane it is free
    to translate and twist as a rigid body, and a three-dimensional run of it is refused as a
    mechanism.
    """

    results = beam_model('cantilever').analyze_modal(
        num_modes=2, mass_combo_name='Mass', gravity=GRAVITY
    )

    assert results.diagnostics.plane is None
    assert '3D' in results.diagnostics.summary()


def test_diagnostics_report_whether_shear_deformation_was_active():
    """Forcing Timoshenko on does nothing unless a section defines a shear area.

    `add_section` defaults both shear areas to zero, and the stiffness matrix skips the shear
    correction when they are, so a model built without them is solved as Euler-Bernoulli despite the
    forcing. Reporting which one was used keeps that visible, and flags in advance that supplying
    shear areas will move every modal frequency.
    """

    without = modal(beam_model('simply supported'), num_modes=1)

    assert without.diagnostics.shear_deformation is False
    assert 'Euler-Bernoulli' in without.diagnostics.summary()

    # The same beam with shear areas defined
    model = beam_model('simply supported')
    model.sections['Section'].Asy = 0.006
    model.sections['Section'].Asz = 0.004
    with_shear = modal(model, num_modes=1)

    assert with_shear.diagnostics.shear_deformation is True
    assert 'Timoshenko' in with_shear.diagnostics.summary()

    # Shear flexibility can only lower the frequency
    assert with_shear.frequencies[0] < without.frequencies[0]


# ---------------------------------------------------------------------------------------------
# Isolation, without subdivision
# ---------------------------------------------------------------------------------------------

def test_static_results_survive_modal_analysis_without_subdivision():
    """Isolation must not depend on the arguments given.

    Preparing a model for analysis clears every stored nodal displacement, so a modal run that
    solved on the caller's own model would erase their static results without raising anything. The
    `elements_per_member=1, plane=None` combination is the one that has nothing to subdivide and
    nothing to restrain, and so the one where an "only copy when needed" optimisation would skip the
    copy and do the damage.
    """

    model = beam_model('simply supported', brace_out_of_plane=True)
    model.add_member_dist_load('M1', 'FY', -10.0, -10.0, case='Q')
    model.add_load_combo('Static', {'Q': 1.0})
    model.add_node('Mid', L/2, 0, 0)
    model.nodes['Mid'].support_DZ = True
    model.nodes['Mid'].support_RX = True
    model.nodes['Mid'].support_RY = True

    model.analyze_linear()
    static_before = {name: node.DY['Static'] for name, node in model.nodes.items()}

    assert any(value != 0 for value in static_before.values()), 'The fixture computed no deflection'

    # No subdivision and no plane restraint: nothing for the analysis to change
    model.analyze_modal(num_modes=2, mass_combo_name='Mass', gravity=GRAVITY,
                        elements_per_member=1)

    for name, node in model.nodes.items():
        assert 'Static' in node.DY, f"Static results were erased from node {name}"
        assert node.DY['Static'] == static_before[name]


@pytest.mark.parametrize('kwargs', [
    {'elements_per_member': 1},
    {'elements_per_member': 1, 'plane': 'XY'},
    {'elements_per_member': 4},
    {'elements_per_member': 4, 'plane': 'XY'},
])
def test_model_geometry_survives_every_argument_combination(kwargs):
    """The caller's geometry and supports must come back untouched however modal was called."""

    model = beam_model('simply supported')
    nodes_before = set(model.nodes)
    supports_before = {name: (node.support_DX, node.support_DY, node.support_DZ,
                              node.support_RX, node.support_RY, node.support_RZ)
                       for name, node in model.nodes.items()}

    model.analyze_modal(num_modes=2, mass_combo_name='Mass', gravity=GRAVITY, **kwargs)

    assert set(model.nodes) == nodes_before
    assert not any(name.startswith('_modal_') for name in model.nodes)

    for name, node in model.nodes.items():
        assert (node.support_DX, node.support_DY, node.support_DZ,
                node.support_RX, node.support_RY, node.support_RZ) == supports_before[name]
