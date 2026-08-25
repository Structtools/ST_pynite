# Modal analysis

Natural frequencies and mode shapes for frame models. This guide covers how to drive
`FEModel3D.analyze_modal()`, what it hands back, and the four things that most often make a modal
result wrong.

Plate and shell elements are not supported. Mass is assembled from members, springs and nodes only.

---

## Read this first: the unit contract

Pynite is unit-agnostic, which puts the burden of consistency on you. Two conventions cause almost
every wrong answer in modal analysis.

**`gravity` defaults to `1.0`, and that default is almost never what you want.** Mass is derived
from loads by dividing by gravity, so leaving the default in place overstates the mass by a factor
of *g* and understates every frequency by a factor of √*g* — about 3.13× too low. Always pass it:

```python
model.analyze_modal(num_modes=6, mass_combo_name='Mass', gravity=9.81)
```

**`Material.rho` is a weight density, not a mass density.** `add_member_self_weight()` builds the
self-weight load as `rho * A` and reads the result as a force per unit length; the mass matrix then
divides by `gravity`. So with base units of kN, m and s:

| Material | `rho` to use | Not |
|---|---|---|
| Steel | `78.5` kN/m³ | ~~7.85~~ |
| Concrete | `25.0` kN/m³ | ~~2.5~~ |
| Timber (C24) | `4.2` kN/m³ | ~~0.42~~ |

Mass then comes out in kN·s²/m, which is a **tonne**. Nodal masses you supply directly with
`add_node_mass()` are in those same mass units — tonnes, not kN.

A factor-of-1000 slip in density moves every frequency by √1000 ≈ 31.6×, so it is worth a sanity
check on every run. `results.total_mass` and `results.mass_per_node` exist for exactly that.

### The sign of `rho` is ignored

Only the magnitude of the density is used when mass is assembled — in all three mass paths. A
**negative weight density is a supported way to express the direction of self-weight**, which is
convenient in a Z-up model:

```python
# These two are equivalent, in both the static load and the modal mass
model.add_material('Steel', E, G, nu, rho=-78.5)   # sign in the density
model.add_member_self_weight('FZ', 1.0, 'D')

model.add_material('Steel', E, G, nu, rho=+78.5)   # sign in the factor
model.add_member_self_weight('FZ', -1.0, 'D')
```

Mass is not a signed quantity, so neither convention can produce a negative element mass or an
indefinite mass matrix. This is a contract, not an implementation detail —
`test_density_sign_is_ignored_when_assembling_mass` pins it.

---

## Quick start

```python
from Pynite import FEModel3D

model = FEModel3D()
model.add_node('N1', 0, 0, 0)
model.add_node('N2', 6, 0, 0)
model.def_support('N1', True, True, True, True, True, False)
model.def_support('N2', False, True, True, True, True, False)

model.add_material('Steel', E=210e6, G=81e6, nu=0.3, rho=78.5)   # kN/m3, weight density
model.add_section('IPE200', A=0.00285, Iy=1.42e-6, Iz=1.94e-5, J=6.98e-8)
model.add_member('B1', 'N1', 'N2', 'Steel', 'IPE200')

# Self-weight has to exist as a load and be named in the mass combination
model.add_member_self_weight('FY', -1.0, 'D')
model.add_load_combo('Mass', {'D': 1.0})

results = model.analyze_modal(num_modes=6, mass_combo_name='Mass', gravity=9.81, plane='XY')

for i, f in enumerate(results.frequencies, start=1):
    print(f'Mode {i}: {f:.3f} Hz  (T = {1/f:.4f} s)')

print(results.diagnostics.summary())
print(f'Total mass: {results.total_mass["Y"]:.3f} t')
```

---

## Where mass comes from

Three sources add up. Each can be used alone or together, and none of them double-counts another.

### 1. Self-mass, from material density

Requires **both** a self-weight load and that load's case in the mass combination. Density alone
produces no mass:

```python
model.add_member_self_weight('FY', -1.0, 'D')     # creates the load
model.add_load_combo('Mass', {'D': 1.0})          # brings it into the mass combination
```

Self-mass always uses the consistent mass matrix, regardless of `mass_formulation`.

The factor scales the mass as well as the load, so raising it to account for connections raises the
dynamic mass too:

```python
model.add_member_self_weight('FY', -1.15, 'D')     # 15% for connections -> 15% more mass
```

Both load factors compound: the one given to `add_member_self_weight` and the one the mass
combination applies to its case.

### 2. Load-derived mass, from a mass combination

Any non-self-weight load in the mass combination is converted to mass as `load / gravity`. This is
the path a serviceability mass combination such as *G + ψ₂Q* uses:

```python
model.add_member_dist_load('B1', 'FY', -2.5, -2.5, case='Q')
model.add_load_combo('Mass', {'D': 1.0, 'Q': 0.3})    # G + 0.3Q
```

Mass is sign-insensitive — an upward load carries mass just as a downward one does — and is applied
to **all three** translational directions, so sway modes of frames are captured. `mass_direction`
selects only which *component* of each load is measured, not which DOFs receive it.

**Load components in other directions are dropped, deliberately.** With `mass_direction='Z'` in a
Z-up model, a purely horizontal `FX` load contributes no mass at all. That is correct for mass
derived from gravity, which is the intended use: a horizontal load is wind or notional force, not
weight, and converting it to mass would invent mass that does not exist. If you have mass that is
genuinely not expressed as a gravity load, apply it with `add_node_mass()` instead of relying on the
conversion.

Distributed loads are integrated over their loaded length rather than lumped at an estimated
centroid, so the assembled total equals the applied mass exactly.

### 3. Explicit nodal point masses

For equipment, tanks, and any non-structural mass you have already worked out. This is a real mass
in mass units, so it is **not** divided by gravity and **not** scaled by a load factor:

```python
model.add_node_mass('N5', 2.4)                              # 2.4 tonnes
model.add_node_mass('N5', 0.8, IX=0.15, IY=0.4, IZ=0.4)     # with rotational inertia
```

Calls accumulate, so the node above ends up with 3.2 t. Negative masses are rejected.

### Checking the mass you actually assembled

```python
print(results.total_mass)        # {'X': 12.84, 'Y': 12.84, 'Z': 12.84}
print(results.mass_per_node)     # {'N1': 0.41, 'N2': 0.83, ...} sums to the total
```

`total_mass` is `rᵀMr` for a rigid unit translation, not the sum of the diagonal. A consistent mass
matrix shares mass between coupled DOFs, so its diagonal sums to only 312/420 of the true total in
the transverse directions — do not sanity-check against `M.diagonal()`.

---

## What you get back

`analyze_modal()` returns a `ModalResults`. Mode shapes are *also* stored on the nodes under load
combinations named `Mode 1`, `Mode 2`, … as before.

```python
results.frequencies      # Hz, ascending
results.omega            # rad/s
results.periods          # seconds
results.eigenvalues      # omega**2
results.mode_count

results.mode_shapes      # (n_free_dof, n_modes), mass-normalized and sign-fixed
results.free_dof_indices # global DOF index for each row of mode_shapes
results.dof_map          # [('N1', 'DY'), ...] one entry per row of mode_shapes
results.mode_shape(1)    # {('N1', 'DY'): 0.0123, ...} for mode 1

results.M                  # the assembled global mass matrix (sparse)
results.total_mass         # {'X': ..., 'Y': ..., 'Z': ...} everything assembled
results.participating_mass # {'X': ..., 'Y': ..., 'Z': ...} only what modes can mobilize
results.mass_per_node      # {node_name: mass}
results.mesh_nodes         # {node_name: (X, Y, Z)} including subdivision nodes
results.diagnostics

results.participation_factors('Z')   # per-mode gamma
results.effective_mass('Z')          # per-mode effective modal mass
results.mass_participation('Z')      # per-mode fraction of the participating mass
```

### Which mode actually matters

"The frequency of the beam" is the lowest mode with significant effective mass in the direction you
care about — **not necessarily mode 1**. A simply supported beam's antisymmetric modes have
effective mass of essentially zero, so mode 2 is invisible to a uniform vertical excitation.

```python
import numpy as np

participation = results.mass_participation('Y')          # fraction per mode
cumulative = np.cumsum(participation)

governing = int(np.argmax(participation > 0.05))         # first mode that carries real mass
print(f'Governing: mode {governing + 1} at {results.frequencies[governing]:.2f} Hz '
      f'({100*participation[governing]:.0f}% of participating mass)')

if cumulative[-1] < 0.90:
    print(f'Only {100*cumulative[-1]:.0f}% of mass captured — request more modes')
```

`participation_factors(direction)` and `effective_mass(direction)` give the unnormalized quantities
if you need them. Mode shapes are mass-normalized (`φᵀMφ = I`), so the generalized mass in the
denominator of the usual formulae is one:

- `Γᵢ = φᵢᵀ M r` — `participation_factors()`
- `m_eff,i = Γᵢ²` — `effective_mass()`
- `m_eff,i / participating_mass` — `mass_participation()`

### Why `participating_mass` and not `total_mass`

`mass_participation()` divides by `participating_mass`, and that distinction is the one thing here
worth reading twice. Mass sitting on a **restrained** DOF can never move in any mode, so no number
of modes will recover it. Effective modal mass therefore accumulates towards
`participating_mass[direction]`, not `total_mass[direction]`.

Divide by the total instead and participation looks permanently incomplete — for the beam in the
quick start above, capped around 84% no matter how many modes you compute. A planar analysis makes
this stark: `participating_mass` out of plane is exactly `0.0`, because every out-of-plane DOF is
restrained, while `total_mass` in that direction is unchanged.

Both are exposed so you never have to rebuild `M11` by index-slicing to find out.

---

## Accuracy: element subdivision

A physical member is a single element between its end nodes, and one element per member gets the
higher modes badly wrong. `analyze_modal` therefore subdivides every member for the duration of the
analysis, defaulting to **8 elements per member**.

Measured error against the exact solution for a uniform simply supported beam:

| Elements per member | f₁ | f₂ | f₃ |
|---|---|---|---|
| 1 | +10.99 % | +27.16 % | +305.29 % |
| 2 | +0.39 % | +10.99 % | +23.99 % |
| 4 | +0.03 % | +0.39 % | +1.83 % |
| **8 (default)** | **+0.00 %** | **+0.03 %** | **+0.13 %** |
| 16 | +0.00 % | +0.00 % | +0.01 % |

Convergence is monotonic and from above. The table is committed as
`Testing/test_modal_verification.py::test_convergence_study`, so it can be cited in analysis
documentation.

```python
results = model.analyze_modal(..., elements_per_member=16)   # more accuracy
results = model.analyze_modal(..., elements_per_member=1)    # no subdivision
```

**Your model is never modified, whatever arguments you pass.** The analysis always runs on an
internal copy — even with `elements_per_member=1` and no `plane`, where there is nothing to subdivide
or restrain — so your geometry, supports and static results all survive a modal run untouched. The
guarantee is unconditional rather than something that holds for most argument combinations, because
preparing a model for analysis clears every stored nodal displacement, and a copy skipped as an
optimisation would erase your static results with no error raised.

The subdivision nodes appear in `results.dof_map` and `results.mesh_nodes` — named
`_modal_<member>_<n>` — but never in `model.nodes`.

---

## Planar models: use `plane`

Subdivision creates interior nodes that you never see and therefore cannot restrain. Left free out
of plane, they admit lateral and torsional modes into the results, interleaved with the in-plane
ones. For a doubly symmetric section every frequency comes back *twice*.

If your model is planar, say so:

```python
results = model.analyze_modal(..., plane='XY')   # or 'XZ', 'YZ'
```

This restrains the out-of-plane and torsional DOFs at every node, including the interior ones, on
the internal copy. Omit it for a genuine three-dimensional analysis.

---

## Consistent vs lumped mass

```python
model.analyze_modal(..., mass_formulation='consistent')   # default
model.analyze_modal(..., mass_formulation='lumped')
```

Consistent mass populates the rotational DOFs and converges on the exact frequency **from above**.
Lumped mass puts load-derived mass on the translational DOFs only, is cheaper, and converges **from
below**. Running both brackets the exact answer, which is a cheap and rather effective check on a
coarse mesh. `mass_formulation` affects load-derived mass only; self-mass is always consistent.

---

## Determinism

Results are reproducible, and that is tested rather than assumed:

- Frequencies are invariant to node and member numbering.
- Mode shapes are mass-normalized and sign-fixed — the largest-magnitude component is positive,
  with ties broken on position rather than DOF index so that the convention survives renumbering.
- The eigensolver's starting vector is fixed, so repeated runs agree bit for bit.

Without this, animations flip direction between runs and differential tests fail spuriously.

---

## Diagnostics

`results.diagnostics` records what the solver was asked for and what it actually did. It is the
material for documenting an analysis assumption:

```python
d = results.diagnostics

# What was asked for, and what came back
d.requested_modes
d.converged_modes
d.truncated               # True if fewer converged than requested
d.solver                  # 'sparse-shift-invert' or 'dense'
d.sigma                   # the shift used
d.filtered_modes          # [(index, reason), ...] eigenpairs discarded as artifacts
d.stabilized_dof_count    # massless DOFs given a negligible mass to stay solvable

# Every argument that changes the answer, echoed back
d.mass_combo_name
d.mass_direction
d.gravity                 # read this one back; the 1.0 default is the classic silent error
d.mass_formulation
d.elements_per_member
d.plane
d.linear_state
d.shear_deformation       # whether forced Timoshenko actually did anything

d.summary()               # one-line human-readable version
```

The provenance fields exist so a result is self-describing, which matters in two places. **Cache
keys** can be derived from the result rather than re-threaded from the request, so they cannot drift
out of step with what they describe. **Analysis documentation** needs the mass basis, gravity,
analysis plane and linear-state declaration stated, and `summary()` renders all of it in one line:

```
6 of 6 modes, sparse-shift-invert solver, mass from 'Mass' (Y) at g=9.81, consistent mass,
8 element(s) per member, Euler-Bernoulli, XY plane
```

### `shear_deformation`: forced Timoshenko can be a no-op

Modal analysis forces the Timoshenko formulation on for every member, but the shear correction is
skipped when a section's shear area is zero — and `add_section` defaults both shear areas to zero.
**A model built without shear areas is solved as Euler–Bernoulli despite the forcing.**

`diagnostics.shear_deformation` reports which one you actually got. It is worth checking, because the
day you supply `Asy`/`Asz` for some unrelated reason, every modal frequency will change with no
change in the modal code:

```python
model.add_section('IPE200', A, Iy, Iz, J, Asy=0.0014, Asz=0.0019)   # now Timoshenko for real
```

Fewer modes than requested are reported honestly rather than padded, so **check `truncated`**. Modes
that are rigid-body artifacts, or that live only on massless DOFs, are filtered out with the reason
recorded in `filtered_modes`.

---

## Errors

Two exception types separate a modelling problem from a numerical one:

```python
from Pynite.ModalResults import ModalModelError, ModalSolverError
```

`ModalModelError` — the model has no well-defined set of modes, and you need to change it:

- no mass anywhere, or no mass on any free DOF
- every DOF is supported, so nothing can vibrate
- every eigenpair was a rigid-body artifact: the model is a mechanism
- the model contains tension-only or compression-only elements (see below)
- the model holds pushover results, so its state is past yield

`ModalSolverError` — the model is well posed but the eigensolver failed. Usually a conditioning
problem: check for wildly mismatched stiffnesses or masses.

`ValueError` is raised immediately for bad arguments, before any assembly work.

### Tension-only and compression-only elements

Their stiffness depends on the load they carry, so the model has no single linear state to take
modes of. Rather than pick one silently, `analyze_modal` requires a declaration:

```python
model.analyze_modal(..., linear_state='all_active')
```

`'all_active'` computes the modes with every such element engaged, and is the only value supported
at present. Taking the modes of a state derived from a particular load combination is not
implemented.

---

## Limitations

1. **No axial-load effect on frequency.** The geometric stiffness matrix is not included, so a
   compressed member's real frequency is *lower* than reported. Modal results for heavily loaded
   columns are non-conservative. Second-order elastic (P-Delta) models are accepted, but the modes
   are still computed from the elastic stiffness alone.
2. **Elastic, undamped eigenvalues.** This is not a response calculation. Acceptance against any
   comfort or serviceability criterion needs damping and a separate check.
3. **Frame elements only.** Plate and shell elements contribute no mass and no modes.
4. **A single beam is not a floor.** Real floors act as orthotropic plates with continuity, screed,
   partitions and non-structural mass. A beam-strip idealisation can err in *either* direction — it
   is not conservative.
5. **`plane` is not automatic.** Omit it on a planar model and lateral modes will be mixed into the
   results.

---

## Migrating from the previous behaviour

`analyze_modal()` keeps its original parameters and still populates `model.frequencies` and the
`Mode n` load combinations, so existing code keeps working. What changed:

- It now **returns a `ModalResults`** instead of `None`. (Its docstring previously promised a list
  of frequencies but returned nothing.)
- Members are **subdivided 8× by default**, so frequencies will differ from — and be considerably
  more accurate than — earlier results. Pass `elements_per_member=1` for the old meshing.
- Load-derived mass was previously lumped using an incorrect position term, which produced negative
  nodal masses and frequency errors around 14 %. It is now integrated correctly. **Any stored
  frequency that came from a mass combination should be recomputed.**
- Loads applied in member-local directions (`'Fx'`, `'Fy'`, `'Fz'`) raised a `TypeError` when
  converted to mass. They now work.
- The factor given to `add_member_self_weight()` was ignored when assembling mass, so a member whose
  self-weight had been raised to account for connections carried the unfactored mass. It now scales
  the mass as it always scaled the load. **Models using a self-weight factor other than ±1 will
  report lower frequencies than before, correctly.**
- Self-mass contributions are accumulated as magnitudes rather than signed values, so two
  self-weight cases of opposing sign now add instead of cancelling.
- `Mode n` combinations are tagged `'modal'`. Filter on that tag to keep them out of envelopes and
  reports:

  ```python
  design_combos = [name for name, combo in model.load_combos.items()
                   if combo.combo_tags is None or 'modal' not in combo.combo_tags]
  ```

---

## Verification

`Testing/test_modal_verification.py` holds the evidence base: closed-form beam oracles for four
support conditions, SDOF spring–mass and point-mass checks, a portal-frame sway comparison, axial
modes, the element mass matrix against hand-calculated coefficients, exact mass totals for each
source, orthonormality and the Rayleigh quotient, the effective-mass completeness identity,
renumbering invariance, the convergence study above, and the model-isolation guarantees.
