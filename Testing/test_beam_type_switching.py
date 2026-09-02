"""Tests for the beam_type switching mechanism (Bernoulli-Euler vs Timoshenko).

Covers:
  1. Validation of beam_type parameter
  2. beam_type='bernoulli' suppresses shear deformation even with Asy/Asz > 0
  3. _use_timoshenko property logic
  4. Eigenvalue analysis (buckling + modal) forces Timoshenko regardless of beam_type
  5. PhysMember propagation of beam_type and _force_timoshenko to sub-members
"""

import math
import pytest
import numpy as np
from Pynite import FEModel3D

# ---------------------------------------------------------------------------
# Section constants: 300×300 mm solid rectangular, L = 1 m (L/d ≈ 3.3)
# Shear effects are ~28% of total deflection — easily detectable.
# ---------------------------------------------------------------------------
_b  = 0.3                           # width [m]
_d  = 0.3                           # depth [m]
_L  = 1.0                           # span [m]
_A  = _b * _d
_I  = _b * _d**3 / 12               # Iy = Iz (square)
_J  = 2 * _I                        # torsional constant (approx)
_As = 5 / 6 * _A                    # rectangular shear area
_E  = 210e9
_G  = 80.769e9
_P  = -100e3                         # tip load [N]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_cantilever(beam_type='timoshenko', direction='FY'):
    """Cantilever with 300×300 square section, tip load, given beam_type."""
    model = FEModel3D()
    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', _L, 0, 0)
    model.def_support('N1', True, True, True, True, True, True)
    model.add_material('Mat', _E, _G, 0.3, 0.0)
    model.add_section('Sec', _A, _I, _I, _J, Asy=_As, Asz=_As)
    model.add_member('M1', 'N1', 'N2', 'Mat', 'Sec', beam_type=beam_type)
    model.add_node_load('N2', direction, _P, case='Case 1')
    model.add_load_combo('Combo 1', {'Case 1': 1.0})
    model.analyze_linear(log=False)
    return model


def _build_buckling_column(beam_type='timoshenko'):
    """Pinned-pinned column along global Y for buckling analysis."""
    n_elem = 10
    L = 1.0
    E = 1.0
    G = 0.4          # G = E/(2(1+nu)), nu=0.25
    A = 100.0         # large A to push torsional mode above flexural
    Iy = 1.0
    Iz = 1.0
    J = 2.0
    As = 5 / 6 * A    # shear area

    model = FEModel3D()
    h = L / n_elem

    model.add_material('Mat', E, G, 0.25, 1.0)
    model.add_section('Sec', A, Iy, Iz, J, Asy=As, Asz=As)

    for i in range(n_elem + 1):
        model.add_node(f'N{i+1}', 0.0, i * h, 0.0)

    # Pinned base
    model.def_support('N1', True, True, True, False, True, False)
    # Pinned roller top
    model.def_support(f'N{n_elem+1}', True, False, True, False, True, False)

    model.add_member('Column', 'N1', f'N{n_elem+1}', 'Mat', 'Sec',
                     beam_type=beam_type)

    model.add_node_load(f'N{n_elem+1}', 'FY', -1.0)

    return model


def _build_modal_beam(beam_type='timoshenko'):
    """Cantilever beam for modal analysis with self-weight as mass."""
    model = FEModel3D()

    L = 5.0
    E = 200e9
    G = 80e9
    rho = 7800.0
    A = 0.01
    I = 8.33e-6
    J = 1.67e-5
    As = 5 / 6 * A

    # 5 intermediate nodes
    for i in range(6):
        model.add_node(f'N{i+1}', L * i / 5, 0, 0)
        model.def_support(f'N{i+1}', True, False, True, True, True, False)

    model.def_support('N1', True, True, True, True, True, True)

    model.add_material('Steel', E, G, 0.3, rho)
    model.add_section('Sec', A, I, I, J, Asy=As, Asz=As)
    model.add_member('M1', 'N1', 'N6', 'Steel', 'Sec', beam_type=beam_type)

    model.add_member_self_weight('FY', -1.0, 'Mass')
    model.add_load_combo('MassCombo', {'Mass': 1.0})

    return model


# ===========================================================================
# Group 1: beam_type validation
# ===========================================================================

class TestBeamTypeValidation:

    def test_invalid_beam_type_raises(self):
        """Passing an invalid beam_type to add_member must raise ValueError."""
        model = FEModel3D()
        model.add_node('N1', 0, 0, 0)
        model.add_node('N2', 1, 0, 0)
        model.add_material('Mat', 210e9, 80e9, 0.3, 0.0)
        model.add_section('Sec', 0.01, 1e-4, 1e-4, 1e-5, 0.0, 0.0)

        with pytest.raises(ValueError, match="beam_type must be"):
            model.add_member('M1', 'N1', 'N2', 'Mat', 'Sec', beam_type='foo')

    def test_valid_beam_types_accepted(self):
        """Both 'timoshenko' and 'bernoulli' must be accepted without error."""
        model = FEModel3D()
        model.add_node('N1', 0, 0, 0)
        model.add_node('N2', 1, 0, 0)
        model.add_node('N3', 2, 0, 0)
        model.add_material('Mat', 210e9, 80e9, 0.3, 0.0)
        model.add_section('Sec', 0.01, 1e-4, 1e-4, 1e-5, 0.0, 0.0)

        model.add_member('M1', 'N1', 'N2', 'Mat', 'Sec', beam_type='timoshenko')
        model.add_member('M2', 'N2', 'N3', 'Mat', 'Sec', beam_type='bernoulli')

        assert model.members['M1'].beam_type == 'timoshenko'
        assert model.members['M2'].beam_type == 'bernoulli'


# ===========================================================================
# Group 2: beam_type='bernoulli' suppresses shear deformation
# ===========================================================================

class TestBernoulliSuppressesShear:

    def test_bernoulli_suppresses_shear_strong_axis(self):
        """beam_type='bernoulli' with non-zero Asy must give Euler-Bernoulli deflection."""
        model = _build_cantilever(beam_type='bernoulli', direction='FY')

        delta_fem = model.nodes['N2'].DY['Combo 1']
        delta_eb = _P * _L**3 / (3 * _E * _I)  # pure Euler-Bernoulli

        assert math.isclose(delta_fem, delta_eb, rel_tol=1e-6), \
            f'Expected EB: {delta_eb:.6e}, got: {delta_fem:.6e}'

    def test_bernoulli_suppresses_shear_weak_axis(self):
        """beam_type='bernoulli' with non-zero Asz must give Euler-Bernoulli for FZ."""
        model = _build_cantilever(beam_type='bernoulli', direction='FZ')

        delta_fem = model.nodes['N2'].DZ['Combo 1']
        delta_eb = _P * _L**3 / (3 * _E * _I)  # Iy = Iz for square section

        assert math.isclose(delta_fem, delta_eb, rel_tol=1e-6), \
            f'Expected EB: {delta_eb:.6e}, got: {delta_fem:.6e}'

    def test_bernoulli_vs_timoshenko_difference_equals_shear_term(self):
        """Difference between timoshenko and bernoulli deflections must equal P·L/(G·As)."""
        model_b = _build_cantilever(beam_type='bernoulli', direction='FY')
        model_t = _build_cantilever(beam_type='timoshenko', direction='FY')

        delta_b = model_b.nodes['N2'].DY['Combo 1']
        delta_t = model_t.nodes['N2'].DY['Combo 1']

        shear_term = _P * _L / (_G * _As)
        diff = delta_t - delta_b

        assert math.isclose(diff, shear_term, rel_tol=1e-4), \
            f'Shear increment: diff={diff:.6e}, expected={shear_term:.6e}'


# ===========================================================================
# Group 3: _use_timoshenko property
# ===========================================================================

class TestUseTimoshenkoProperty:

    def _make_member(self, beam_type):
        """Create a simple model and return its member."""
        model = FEModel3D()
        model.add_node('N1', 0, 0, 0)
        model.add_node('N2', 1, 0, 0)
        model.add_material('Mat', 210e9, 80e9, 0.3, 0.0)
        model.add_section('Sec', 0.01, 1e-4, 1e-4, 1e-5, 0.0, 0.0)
        model.add_member('M1', 'N1', 'N2', 'Mat', 'Sec', beam_type=beam_type)
        return model.members['M1']

    def test_use_timoshenko_reflects_beam_type(self):
        """_use_timoshenko must reflect the beam_type setting."""
        member_t = self._make_member('timoshenko')
        member_b = self._make_member('bernoulli')

        assert member_t._use_timoshenko is True
        assert member_b._use_timoshenko is False

    def test_force_flag_overrides_bernoulli(self):
        """_force_timoshenko=True must override beam_type='bernoulli'."""
        member = self._make_member('bernoulli')
        assert member._use_timoshenko is False

        member._force_timoshenko = True
        assert member._use_timoshenko is True

        member._force_timoshenko = False
        assert member._use_timoshenko is False


# ===========================================================================
# Group 4: Eigenvalue analysis forces Timoshenko
# ===========================================================================

class TestEigenvalueForceTimoshenko:

    def test_buckling_forces_timoshenko(self):
        """Buckling with beam_type='bernoulli' must match beam_type='timoshenko'
        (both are forced to Timoshenko during eigenvalue analysis)."""
        model_b = _build_buckling_column(beam_type='bernoulli')
        model_t = _build_buckling_column(beam_type='timoshenko')

        results_b = model_b.analyze_buckling(num_modes=1)
        results_t = model_t.analyze_buckling(num_modes=1)

        lam_b = results_b.load_multipliers[0]
        lam_t = results_t.load_multipliers[0]

        assert math.isclose(lam_b, lam_t, rel_tol=1e-6), \
            f'Buckling: bernoulli λ={lam_b:.6f}, timoshenko λ={lam_t:.6f}'

    def test_buckling_resets_force_flag(self):
        """After buckling analysis, _force_timoshenko must be False on all members."""
        model = _build_buckling_column(beam_type='bernoulli')
        model.analyze_buckling(num_modes=1)

        phys = model.members['Column']
        assert phys._force_timoshenko is False, \
            'PhysMember _force_timoshenko not reset after buckling'

        for sub in phys.sub_members.values():
            assert sub._force_timoshenko is False, \
                f'Sub-member {sub.name} _force_timoshenko not reset after buckling'

    def test_modal_forces_timoshenko(self):
        """Modal analysis with beam_type='bernoulli' must give same frequencies
        as beam_type='timoshenko' (both forced to Timoshenko)."""
        model_b = _build_modal_beam(beam_type='bernoulli')
        model_t = _build_modal_beam(beam_type='timoshenko')

        model_b.analyze_modal(num_modes=2, mass_combo_name='MassCombo',
                              mass_direction='Y', gravity=9.81, log=False)
        model_t.analyze_modal(num_modes=2, mass_combo_name='MassCombo',
                              mass_direction='Y', gravity=9.81, log=False)

        for f_b, f_t in zip(model_b.frequencies, model_t.frequencies):
            assert math.isclose(f_b, f_t, rel_tol=1e-6), \
                f'Modal: bernoulli f={f_b:.6f} Hz, timoshenko f={f_t:.6f} Hz'

    def test_modal_resets_force_flag(self):
        """After modal analysis, _force_timoshenko must be False on all members."""
        model = _build_modal_beam(beam_type='bernoulli')
        model.analyze_modal(num_modes=2, mass_combo_name='MassCombo',
                            mass_direction='Y', gravity=9.81, log=False)

        phys = model.members['M1']
        assert phys._force_timoshenko is False, \
            'PhysMember _force_timoshenko not reset after modal'

        for sub in phys.sub_members.values():
            assert sub._force_timoshenko is False, \
                f'Sub-member {sub.name} _force_timoshenko not reset after modal'


# ===========================================================================
# Group 5: PhysMember propagation
# ===========================================================================

class TestPhysMemberPropagation:

    def test_sub_members_inherit_beam_type(self):
        """Sub-members created during discretization must inherit beam_type."""
        model = FEModel3D()
        model.add_node('N1', 0, 0, 0)
        model.add_node('N2', 5, 0, 0)   # intermediate node
        model.add_node('N3', 10, 0, 0)
        model.def_support('N1', True, True, True, True, True, True)
        model.add_material('Mat', 210e9, 80e9, 0.3, 0.0)
        model.add_section('Sec', 0.01, 1e-4, 1e-4, 1e-5, Asy=0.008, Asz=0.008)

        model.add_member('M1', 'N1', 'N3', 'Mat', 'Sec', beam_type='bernoulli')
        model.add_node_load('N3', 'FY', -1000)
        model.add_load_combo('Combo 1', {'Case 1': 1.0})
        model.analyze_linear(log=False)

        phys = model.members['M1']
        assert len(phys.sub_members) == 2, \
            f'Expected 2 sub-members, got {len(phys.sub_members)}'

        for name, sub in phys.sub_members.items():
            assert sub.beam_type == 'bernoulli', \
                f'Sub-member {name} has beam_type={sub.beam_type}, expected bernoulli'

    def test_sub_members_inherit_force_timoshenko(self):
        """Sub-members must inherit _force_timoshenko from parent during discretization."""
        model = FEModel3D()
        model.add_node('N1', 0, 0, 0)
        model.add_node('N2', 5, 0, 0)   # intermediate node
        model.add_node('N3', 10, 0, 0)
        model.def_support('N1', True, True, True, True, True, True)
        model.add_material('Mat', 210e9, 80e9, 0.3, 0.0)
        model.add_section('Sec', 0.01, 1e-4, 1e-4, 1e-5, Asy=0.008, Asz=0.008)

        model.add_member('M1', 'N1', 'N3', 'Mat', 'Sec', beam_type='bernoulli')

        phys = model.members['M1']
        phys._force_timoshenko = True
        phys.descritize()

        assert len(phys.sub_members) == 2, \
            f'Expected 2 sub-members, got {len(phys.sub_members)}'

        for name, sub in phys.sub_members.items():
            assert sub._force_timoshenko is True, \
                f'Sub-member {name} has _force_timoshenko=False, expected True'
