"""
FEM Design comparison of the in-plane linear buckling analysis.

The scripts ``hinged_frame_comparison.py``, ``pitched_frame_comparison.py`` and
``braced_pitched_frame_comparison.py`` in this directory carry FEM Design 2024
critical parameters for their scenarios (an IPE300 portal with a corner hinge, an
IPE200 gable frame, and the same gable with a diagonal brace, all with fixed bases
and out-of-plane DOFs restrained). The scripts print a comparison table; this test
turns it into a gate.

Each scenario's first modes must lie within ±5 % of FEM Design's. The scripts
expose their scenario presets through ``_apply_scenario()`` and ``build_model()``,
so the same model the script prints is the one tested here.

Units: SI (N, m).
"""

import contextlib
import importlib
import io

import numpy as np
import pytest

MARGIN_PCT = 5.0

SCRIPTS = {
    "Testing.hinged_frame_comparison": "FEM_DESIGN_REF",
    "Testing.pitched_frame_comparison": "_FEM_DESIGN_REF",
    "Testing.braced_pitched_frame_comparison": "_FEM_DESIGN_REF",
}


def _scenarios():
    for module_name, attribute in SCRIPTS.items():
        module = importlib.import_module(module_name)
        for scenario in sorted(getattr(module, attribute)):
            yield pytest.param(module_name, attribute, scenario,
                               id=f"{module_name.split('.')[-1]}-scenario{scenario}")


@pytest.mark.parametrize("module_name,attribute,scenario", list(_scenarios()))
def test_fem_design_critical_parameters(module_name, attribute, scenario):
    module = importlib.import_module(module_name)
    reference = np.asarray(getattr(module, attribute)[scenario], dtype=float)

    module.SCENARIO = scenario
    module._apply_scenario()
    with contextlib.redirect_stdout(io.StringIO()):
        model = module.build_model()
        model = model[0] if isinstance(model, tuple) else model
        results = model.analyze_buckling(
            combo_name=module.COMBO_NAME, num_modes=max(len(reference), 4)
        )

    ours = np.asarray(results.load_multipliers[: len(reference)])
    assert ours.size == reference.size, (
        f"{module_name} scenario {scenario}: FEM Design lists {reference.size} modes, "
        f"PyNite found {ours.size}"
    )
    diff_pct = (ours - reference) / reference * 100.0
    worst = int(np.argmax(np.abs(diff_pct)))
    assert np.all(np.abs(diff_pct) <= MARGIN_PCT), (
        f"{module_name} scenario {scenario}: mode {worst + 1} PyNite {ours[worst]:.3f} vs "
        f"FEM Design {reference[worst]:.3f} ({diff_pct[worst]:+.2f} %); all modes: "
        + ", ".join(f"{d:+.2f} %" for d in diff_pct)
    )
