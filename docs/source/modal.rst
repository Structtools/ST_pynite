==============
Modal Analysis
==============

Modal analysis computes a structure's natural frequencies and mode shapes. Pynite provides
the ``FEModel3D.analyze_modal()`` method which assembles the global stiffness and mass
matrices, solves the generalized eigenvalue problem, and returns a structured result. Mode
shapes are also stored as displacement results in modal load combinations.

Frame elements only. Plate and shell elements contribute no mass and no modes.

.. note::

   A task-oriented guide with worked examples, the unit contract, the convergence study behind
   the subdivision default, and post-processing recipes for participation factors and effective
   modal mass lives in ``docs/MODAL_ANALYSIS.md``.

How it works
============

- Subdivides every physical member into ``elements_per_member`` elements for the duration of
  the analysis. This happens on an internal copy, so the caller's model is never modified.
- Assembles the global stiffness matrix ``[K]`` and mass matrix ``[M]`` using the
  specified ``mass_combo_name`` and ``mass_direction``.
- Partitions out supported degrees of freedom before solving to avoid singularities.
- Solves the generalized eigenproblem ``[K]{φ} = λ[M]{φ}`` where ``λ = ω²``. Frequencies
  are returned in Hz as ``f = ω / (2π)``.
- Discards eigenpairs that are rigid-body artifacts or that live only on massless degrees of
  freedom, recording each one in the result's diagnostics.
- Mass-normalizes the mode shapes and applies a deterministic sign convention, so that
  results are reproducible and invariant to node and member numbering.
- Mode shapes (eigenvectors) are expanded back into the model's full DOF set and
  stored in load combinations named ``Mode 1``, ``Mode 2``, etc. Those combinations are
  tagged ``'modal'`` so that envelope and report code can filter them out.

Usage
=====

``analyze_modal`` signature
---------------------------

``FEModel3D.analyze_modal(num_modes: int = 12, mass_combo_name: str = 'Combo 1', mass_direction: str = 'Y', gravity: float = 1.0, log=False, check_stability=True, mass_formulation: str = 'consistent', elements_per_member: int = 8, plane: str | None = None, linear_state: str | None = None)``

Parameters
----------

- ``num_modes``: Number of modes to calculate (default: ``12``).
- ``mass_combo_name``: Name of the load combination used to convert loads to masses.
  Defaults to ``'Combo 1'``. If no load combos exist, a default is created automatically.
- ``mass_direction``: Direction to use when measuring loads for conversion to mass
  (``'X'``, ``'Y'``, or ``'Z'``, default ``'Y'``). This selects which *component* of each
  load is measured. The resulting mass is applied to all three translational directions, so
  sway modes are captured.
- ``gravity``: Acceleration used when converting loads to mass (default: ``1.0``). See the
  warning below: the default is almost never correct.
- ``log``: If ``True``, prints progress messages to the console.
- ``check_stability``: If ``True``, checks the stiffness matrix for unsupported DOFs
  and will raise an exception if instabilities are found.
- ``mass_formulation``: How load-derived mass is distributed, ``'consistent'`` (default) or
  ``'lumped'``. Consistent mass populates the rotational DOFs and converges from above;
  lumped mass is cheaper and converges from below. Self-weight is always consistent.
- ``elements_per_member``: Elements each physical member is subdivided into for the analysis
  (default: ``8``). One element per member leaves the higher modes badly wrong.
- ``plane``: Restricts the analysis to a plane by restraining the out-of-plane degrees of
  freedom at every node: ``'XY'``, ``'XZ'``, ``'YZ'``, or ``None`` (default) for a full
  three-dimensional analysis. Declare this for a planar model — subdivision creates interior
  nodes the caller cannot restrain, and left free they admit lateral and torsional modes.
- ``linear_state``: Required declaration when the model contains tension-only or
  compression-only elements, whose stiffness depends on the load they carry. Pass
  ``'all_active'`` to compute the modes with every such element engaged.

Return value
------------

A ``ModalResults`` object exposing ``frequencies`` (Hz), ``omega``, ``periods``,
``eigenvalues``, the mass-normalized ``mode_shapes`` with their ``free_dof_indices`` and
``dof_map``, the assembled mass matrix ``M``, ``total_mass``, ``participating_mass`` and
``mass_per_node`` for sanity-checking, ``mesh_nodes``, and ``diagnostics``.

It also provides ``participation_factors(direction)``, ``effective_mass(direction)`` and
``mass_participation(direction)``, which answer "which mode actually matters in this direction?"
without the caller having to partition the mass matrix themselves.

.. important::

   ``mass_participation`` divides by ``participating_mass``, not ``total_mass``. Mass held on a
   restrained degree of freedom cannot move in any mode, so no number of modes will recover it and
   effective modal mass accumulates towards the participating mass. Dividing by the total instead
   makes participation look permanently incomplete. For a planar analysis, ``participating_mass``
   out of plane is exactly zero while ``total_mass`` in that direction is unchanged.

Units
=====

Pynite is unit-agnostic, so consistency is the caller's responsibility. Two conventions cause
most incorrect modal results.

.. warning::

   ``gravity`` defaults to ``1.0``. Mass is derived from loads by dividing by gravity, so
   leaving the default in place overstates the mass by a factor of *g* and understates every
   frequency by a factor of √*g* — roughly 3.13 times too low. Always pass it explicitly.

.. warning::

   ``Material.rho`` is a **weight** density, not a mass density. ``add_member_self_weight()``
   builds the self-weight load as ``rho * A`` and reads it as a force per unit length, and the
   mass matrix then divides by ``gravity``. With base units of kN, m and s, steel is therefore
   ``78.5`` (kN/m³) rather than ``7.85`` (t/m³), and the assembled mass comes out in kN·s²/m,
   which is a tonne. Masses passed to ``add_node_mass()`` are in those mass units.

.. note::

   The **sign** of ``rho`` is ignored when mass is assembled; only its magnitude is used. A negative
   weight density is therefore a supported way to express the direction of self-weight, and
   ``rho = -78.5`` with ``factor = 1.0`` gives exactly the same frequencies as ``rho = +78.5`` with
   ``factor = -1.0``. Mass is not a signed quantity, so neither convention can produce a negative
   element mass.

Mass sources
============

Three sources add up, and none double-counts another:

- **Self-mass** from material density. Requires both a self-weight load and that load's case
  in the mass combination — density alone produces no mass. The factor given to
  ``add_member_self_weight`` scales the mass as well as the load, so raising it to account for
  connections raises the dynamic mass too.
- **Load-derived mass** from the mass combination, converted as ``load / gravity``. Mass is
  sign-insensitive, and distributed loads are integrated over their loaded length rather than
  lumped at an estimated centroid, so the assembled total is exact.
- **Explicit nodal point masses** via ``FEModel3D.add_node_mass()``. These are real masses in
  mass units, so they are neither divided by gravity nor scaled by a load factor.

Requirements and notes
======================

- SciPy is required for modal analysis. The method uses SciPy's sparse eigen-solver
  internally; if SciPy is not installed the method will raise an exception.
- The solver falls back to a dense solve for models too small for the sparse iterative
  solver, and reports partial convergence rather than padding the result. Check
  ``results.diagnostics.truncated``.
- Ensure materials have density with self-weight in the mass combination, that the mass
  combination contains loads that can be converted to mass, or that nodal point masses have
  been assigned. A ``ModalModelError`` is raised if no mass terms are found.
- ``total_mass`` is ``rᵀMr`` for a rigid unit translation, not the sum of the matrix diagonal.
  A consistent mass matrix shares mass between coupled DOFs, so its diagonal sums to only
  312/420 of the true total in the transverse directions.
- The analysis always runs on an internal copy of the model, whatever arguments are given, so the
  caller's geometry, supports and static results survive a modal run untouched.
- Modal analysis forces the Timoshenko formulation on, but the shear correction is skipped for a
  section whose shear area is zero, which is what ``add_section`` defaults to. A model built without
  shear areas is solved as Euler–Bernoulli despite the forcing;
  ``results.diagnostics.shear_deformation`` reports which one was actually used.

Errors
======

``ModalModelError`` signals a modelling problem that has to be resolved: no mass, no free
degrees of freedom, a mechanism, an undeclared tension-only or compression-only state, or a
model holding pushover results. ``ModalSolverError`` signals a numerical failure on a model
that is otherwise well posed. Both are importable from ``Pynite.ModalResults``.

Limitations
===========

- **No axial-load effect on frequency.** The geometric stiffness matrix is not included, so a
  compressed member's true frequency is lower than reported. Results for heavily loaded
  columns are non-conservative.
- **Elastic, undamped eigenvalues.** This is not a response calculation.
- **Frame elements only.** Plates and shells contribute no mass and no modes.

Example
-------

Create a simple model, add mass-producing loads or define material densities, and run modal
analysis::

   from Pynite import FEModel3D

   model = FEModel3D()
   model.add_node('N1', 0, 0, 0)
   model.add_node('N2', 6, 0, 0)
   model.def_support('N1', True, True, True, True, True, False)
   model.def_support('N2', False, True, True, True, True, False)

   model.add_material('Steel', 210e6, 81e6, 0.3, 78.5)   # rho in kN/m3
   model.add_section('Section', 0.00285, 1.42e-6, 1.94e-5, 6.98e-8)
   model.add_member('B1', 'N1', 'N2', 'Steel', 'Section')

   # Self-weight must exist as a load and be named in the mass combination
   model.add_member_self_weight('FY', -1.0, 'D')
   model.add_load_combo('Mass', {'D': 1.0})

   results = model.analyze_modal(num_modes=6, mass_combo_name='Mass', gravity=9.81,
                                 plane='XY', log=True)

   for i, f in enumerate(results.frequencies, start=1):
       print(f'Mode {i}: {f:.3f} Hz')

Mode displacement results are stored in load combinations named ``Mode 1``, ``Mode 2``, etc.,
and remain available on the model's nodes.
