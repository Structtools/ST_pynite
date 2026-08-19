"""Structured results and error types for modal analysis.

`FEModel3D.analyze_modal` returns a `ModalResults` instance. It carries everything a caller needs
to post-process the eigenpairs itself: the frequencies, the mass-normalized and sign-fixed mode
shapes together with the DOF map that interprets them, the assembled mass matrix, the mass totals
for sanity-checking, and a record of what the solver actually did.

Nothing in this module knows anything about design codes. Participation factors, effective modal
mass, mode classification and acceptance criteria are all built on top of these results by the
caller.
"""

from __future__ import annotations  # Allows more recent type hints features
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Dict, List, Tuple

import numpy as np

if TYPE_CHECKING:
    from numpy import float64
    from numpy.typing import NDArray


class ModalModelError(Exception):
    """Raised when a model has no well-defined set of modes.

    This signals a modelling problem the caller has to resolve -- no mass, a mechanism, or a
    non-linear state that has no single linear interpretation -- as distinct from a numerical
    failure inside the eigensolver, which raises `ModalSolverError`.
    """


class ModalSolverError(Exception):
    """Raised when the eigensolver fails on a model that is otherwise well posed."""


@dataclass
class ModalDiagnostics:
    """A record of what the eigensolver was asked for and what it actually delivered.

    Every argument that changes the answer is echoed here, so that a result is self-describing.
    That matters for two things in particular: caching, where a key derived from the result cannot
    drift out of step with the request, and documenting an analysis assumption, where the mass basis,
    gravity, analysis plane and linear-state declaration all have to be stated.

    :param requested_modes: The number of modes the caller asked for.
    :param converged_modes: The number of modes returned. Never silently padded: if this is less
                            than `requested_modes` then `truncated` is `True` and the caller is
                            expected to surface that.
    :param solver: Which path was taken, `'sparse-shift-invert'` or `'dense'`.
    :param sigma: The shift used by the sparse solver, or `None` for the dense path.
    :param mass_combo_name: The load combination the mass was derived from.
    :param mass_direction: The direction whose load component was converted to mass.
    :param gravity: The acceleration used to convert loads to mass. This one is worth reading back:
                    it defaults to 1.0, which understates every frequency by a factor of
                    `sqrt(g)`.
    :param mass_formulation: How load-derived mass was distributed, `'consistent'` or `'lumped'`.
    :param elements_per_member: The subdivision level used for the run.
    :param plane: The plane the analysis was restricted to, or `None` for a full 3D analysis.
    :param linear_state: The linear-state declaration the run was given, or `None` if the model
                         needed none.
    :param shear_deformation: Whether shear deformation was actually active. Modal analysis forces
                              the Timoshenko formulation on, but a section with no shear area
                              defined degenerates to Euler-Bernoulli, so forcing it does not by
                              itself mean the model got it.
    :param stabilized_dof_count: How many DOFs had no mass and were given a negligible
                                 stabilization mass so the matrix pencil stayed solvable.
    :param filtered_modes: `(index, reason)` pairs for every raw eigenpair that was discarded.
    """

    requested_modes: int
    converged_modes: int
    solver: str
    sigma: float | None
    mass_combo_name: str
    mass_direction: str
    gravity: float
    mass_formulation: str
    elements_per_member: int
    plane: str | None
    linear_state: str | None
    shear_deformation: bool
    stabilized_dof_count: int
    filtered_modes: List[Tuple[int, str]] = field(default_factory=list)

    @property
    def truncated(self) -> bool:
        """Whether fewer modes were returned than were requested."""

        return self.converged_modes < self.requested_modes

    def summary(self) -> str:
        """Returns a one-line human-readable summary of the solve."""

        text = (f'{self.converged_modes} of {self.requested_modes} modes, {self.solver} solver, '
                f"mass from '{self.mass_combo_name}' ({self.mass_direction}) at g={self.gravity}, "
                f'{self.mass_formulation} mass, {self.elements_per_member} element(s) per member, '
                f'{"Timoshenko" if self.shear_deformation else "Euler-Bernoulli"}')

        text += f', {self.plane} plane' if self.plane else ', 3D'

        if self.linear_state:
            text += f", linear state '{self.linear_state}'"

        if self.filtered_modes:
            text += f', {len(self.filtered_modes)} mode(s) filtered'

        if self.stabilized_dof_count:
            text += f', {self.stabilized_dof_count} massless DOF(s) stabilized'

        return text


@dataclass
class ModalResults:
    """The result of a modal analysis.

    :param frequencies: Natural frequencies in Hz, ascending.
    :param omega: Circular natural frequencies in rad/s, ascending.
    :param eigenvalues: The eigenvalues of the pencil, equal to `omega**2`.
    :param mode_shapes: Mode shapes with one column per mode and one row per free DOF. Mass
                        normalized, so `mode_shapes.T @ M11 @ mode_shapes` is the identity, and
                        sign-fixed so that repeated runs and renumbered models agree.
    :param free_dof_indices: The global DOF index each row of `mode_shapes` corresponds to.
    :param dof_map: The `(node_name, dof)` pair each row of `mode_shapes` corresponds to, where
                    `dof` is one of `'DX'`, `'DY'`, `'DZ'`, `'RX'`, `'RY'`, `'RZ'`.
    :param M: The assembled global mass matrix, over all DOFs rather than just the free ones.
    :param total_mass: Total assembled mass, keyed by direction (`'X'`, `'Y'`, `'Z'`). The three
                       values are equal by construction -- a full translational row of the
                       consistent mass matrix sums to half the element mass in every direction,
                       axial and transverse alike -- and are reported separately so that a
                       direction-dependent mass definition would show up rather than being hidden.
    :param participating_mass: The mass on the *free* degrees of freedom, keyed by direction. This
                               is the mass that modes can actually mobilize, and therefore the
                               denominator for effective modal mass. Mass sitting on a restrained
                               DOF can never participate in any mode, so this is below
                               `total_mass`, and using `total_mass` instead makes participation
                               look permanently incomplete.
    :param mass_per_node: Each node's share of the total mass, summing to the total. This is the
                          node's row sum of the translational mass block, not its diagonal term,
                          because a consistent mass matrix spreads mass across DOFs. The shares are
                          the same whichever direction they are taken in, for the same reason the
                          directional totals agree.
    :param mesh_nodes: The coordinates of every node used in the solve, keyed by name. When the
                       run subdivided members, this includes the temporary nodes that only existed
                       during the analysis, so `dof_map` stays interpretable afterwards.
    :param diagnostics: What the solver was asked for and what it delivered.
    """

    frequencies: NDArray[float64]
    omega: NDArray[float64]
    eigenvalues: NDArray[float64]
    mode_shapes: NDArray[float64]
    free_dof_indices: NDArray
    dof_map: List[Tuple[str, str]]
    M: object
    total_mass: Dict[str, float]
    participating_mass: Dict[str, float]
    mass_per_node: Dict[str, float]
    mesh_nodes: Dict[str, Tuple[float, float, float]]
    diagnostics: ModalDiagnostics

    @property
    def mode_count(self) -> int:
        """The number of modes in this result."""

        return int(self.frequencies.size)

    @property
    def periods(self) -> NDArray[float64]:
        """The natural periods in seconds, `1/f`."""

        return 1.0/self.frequencies

    def __len__(self) -> int:
        return self.mode_count

    def _influence_vector(self, direction: str) -> NDArray[float64]:
        """Returns the rigid-body influence vector for a direction, restricted to the free DOFs.

        :param direction: `'X'`, `'Y'`, or `'Z'`.
        :type direction: str
        :raises ValueError: Occurs when the direction is not one of the three global axes.
        :return: The influence vector over the free degrees of freedom
        :rtype: NDArray[float64]
        """

        if direction not in ('X', 'Y', 'Z'):
            raise ValueError(f"`direction` must be 'X', 'Y', or 'Z'. '{direction}' was given.")

        # A unit rigid-body translation of the whole structure in this direction
        r = np.zeros(self.M.shape[0])
        r['XYZ'.index(direction)::6] = 1.0

        return r[self.free_dof_indices]

    def _free_mass_matrix(self):
        """Returns the mass matrix partitioned down to the free DOFs.

        :return: The free-DOF mass matrix, in whatever sparse or dense form `M` was given in
        """

        free = self.free_dof_indices

        # Slicing twice keeps a sparse matrix sparse, where fancy-indexing both axes at once would
        # not
        return self.M[free][:, free]

    def participation_factors(self, direction: str) -> NDArray[float64]:
        """Returns each mode's participation factor for a rigid-body translation in `direction`.

        `Gamma_i = phi_i.T @ M @ r / (phi_i.T @ M @ phi_i)`, and the generalized mass in the
        denominator is one because the mode shapes are mass-normalized.

        The sign of a participation factor follows the mode shape's sign convention, so it is
        deterministic but carries no physical meaning on its own. Square it, or take its magnitude,
        before comparing modes.

        :param direction: The global direction of interest, `'X'`, `'Y'`, or `'Z'`.
        :type direction: str
        :return: One participation factor per mode
        :rtype: NDArray[float64]
        """

        r = self._influence_vector(direction)
        M11 = self._free_mass_matrix()

        return np.asarray(self.mode_shapes.T @ (M11 @ r)).ravel()

    def effective_mass(self, direction: str) -> NDArray[float64]:
        """Returns each mode's effective modal mass for a rigid-body translation in `direction`.

        `m_eff,i = (phi_i.T @ M @ r)**2 / (phi_i.T @ M @ phi_i)`, which is the participation factor
        squared for mass-normalized mode shapes.

        Divide by `participating_mass[direction]` rather than by `total_mass[direction]` to get the
        fraction of mobilizable mass a mode carries. Summed over a complete set of modes this
        reaches the participating mass exactly.

        :param direction: The global direction of interest, `'X'`, `'Y'`, or `'Z'`.
        :type direction: str
        :return: One effective modal mass per mode
        :rtype: NDArray[float64]
        """

        return self.participation_factors(direction)**2

    def mass_participation(self, direction: str) -> NDArray[float64]:
        """Returns each mode's effective modal mass as a fraction of the participating mass.

        This is the quantity to read to answer "which mode actually matters in this direction?".
        The lowest mode with a significant fraction is the governing one, which is not always
        mode 1: the antisymmetric modes of a symmetric structure participate in nothing.

        :param direction: The global direction of interest, `'X'`, `'Y'`, or `'Z'`.
        :type direction: str
        :return: One fraction per mode
        :rtype: NDArray[float64]
        """

        participating = self.participating_mass[direction]

        # A structure with no mobilizable mass in this direction has no participation to report
        if participating == 0.0:
            return np.zeros(self.mode_count)

        return self.effective_mass(direction)/participating

    def mode_shape(self, mode: int) -> Dict[Tuple[str, str], float]:
        """Returns one mode shape keyed by `(node_name, dof)`.

        :param mode: The 1-based mode number.
        :type mode: int
        :raises IndexError: Occurs when the mode number is outside the range that was solved for.
        :return: The mode shape, keyed by node name and degree of freedom
        :rtype: Dict[Tuple[str, str], float]
        """

        if not 1 <= mode <= self.mode_count:
            raise IndexError(
                f'Mode {mode} is not available. This result holds modes 1 through {self.mode_count}.'
            )

        column = self.mode_shapes[:, mode - 1]

        return {key: float(value) for key, value in zip(self.dof_map, column)}

    def __repr__(self) -> str:
        if self.mode_count == 0:
            return 'ModalResults(no modes)'

        return (f'ModalResults({self.mode_count} modes, '
                f'f1={self.frequencies[0]:.4g} Hz, '
                f'f{self.mode_count}={self.frequencies[-1]:.4g} Hz)')


def fix_mode_shape_signs(phi: NDArray[float64],
                         dof_coords: List[Tuple[float, float, float, int]]) -> NDArray[float64]:
    """Applies a deterministic sign convention to a set of mode shapes, in place.

    An eigenvector is arbitrary in sign, so two runs of the same model -- or two models that differ
    only in the order their nodes were defined -- can return mode shapes that point opposite ways.
    That flips animations between runs and makes differential tests fail for no real reason.

    Each mode is flipped so that its largest-magnitude component is positive. Ties are broken on
    the *geometry* of the DOF rather than its index, because the index depends on the numbering
    that this convention exists to be independent of. An antisymmetric mode of a symmetric
    structure has two equal and opposite peaks and is decided entirely by the tie-break, so this
    distinction matters in practice rather than only in principle.

    :param phi: Mode shapes, one column per mode. Modified in place.
    :type phi: NDArray[float64]
    :param dof_coords: The `(X, Y, Z, dof)` key of each row of `phi`, used to break ties.
    :type dof_coords: List[Tuple[float, float, float, int]]
    :return: The sign-fixed mode shapes
    :rtype: NDArray[float64]
    """

    for j in range(phi.shape[1]):

        column = phi[:, j]
        magnitudes = np.abs(column)
        peak = magnitudes.max()

        # A mode of all zeros has no sign to fix
        if peak == 0.0:
            continue

        # Collect every component that shares the peak magnitude to within a tolerance, so that a
        # symmetric structure's equal-and-opposite peaks are recognized as the tie they are
        candidates = np.flatnonzero(magnitudes >= peak*(1 - 1e-9))

        # Choose between tied peaks by position, which does not depend on node numbering
        reference = min(candidates, key=lambda i: dof_coords[i])

        if column[reference] < 0:
            phi[:, j] = -column

    return phi
