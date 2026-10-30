"""
PolyFrame comparison of the linear buckling analysis.

For every case in ``Testing/polyframe/cases.py`` that has a captured record in
``Testing/polyframe/corpus/``, PyNite's critical load factors are compared with
PolyFrame's general method (EN 1993-1-1 6.3.4):

* in plane:  ``analyze_buckling(plane='XY')``                   vs ``alpha_cr``
* out of plane: ``analyze_buckling(plane='XY', out_of_plane=True)`` vs ``alpha_cr,op``

The margin is ±5 % on the factor. A case whose catalogue entry declares a
``KnownDeviation`` for a key is reported as an expected failure while PyNite lies
**above** PolyFrame by between the margin and the deviation's bound (the
unconservative direction, which is where every recorded deviation sits), and fails
beyond the bound or if PyNite drops below PolyFrame by more than the margin. So a
characterised disagreement cannot grow unnoticed and cannot hide a change of sign.

pytest never starts PolyFrame: the records are captured by hand with
``python -m Testing.polyframe refresh`` on a licensed Windows machine, and a test
below checks that nothing under ``Testing/`` imports the bridge runner.

Units: SI (N, m). Steel S355, E = 210 GPa, G = E / 2.6.
"""

import ast
from pathlib import Path

import pytest

from Testing.polyframe import pynite_side
from Testing.polyframe.cases import FrameCase, all_cases
from Testing.polyframe.corpus import read_record
from Testing.polyframe.translate import model_hash, to_polyframe_model

MARGIN_PCT = 5.0

# PolyFrame's out-of-plane factor is compared with PyNite's first out-of-plane mode
# computed with the equivalent-torsion approximation of the warping stiffness on the
# twist-restraint spacing of each member (``pynite_side.warping_length_factor``).
# The no-warping value and the member-length value are reported in the summary.
VARIANT = "restraints"


def _captured_cases() -> list[FrameCase]:
    cases = []
    for case in all_cases():
        try:
            read_record(case.case_id)
        except FileNotFoundError:
            continue
        cases.append(case)
    return cases


_CASES = _captured_cases()
_IDS = [c.case_id for c in _CASES]

_RESULTS: dict[str, pynite_side.PyniteBuckling] = {}


def _pynite(case: FrameCase) -> pynite_side.PyniteBuckling:
    if case.case_id not in _RESULTS:
        _RESULTS[case.case_id] = pynite_side.run(case)
    return _RESULTS[case.case_id]


def _record(case: FrameCase) -> dict:
    record = read_record(case.case_id)
    assert record["provenance"]["model_sha256"] == model_hash(to_polyframe_model(case)), (
        f"{case.case_id}: the corpus record was captured from a different model; "
        "run `python -m Testing.polyframe refresh --case " + case.case_id + "`"
    )
    return record


def _compare(case: FrameCase, key: str, ours: float | None, theirs: float | None) -> None:
    assert theirs is not None, f"{case.case_id}: PolyFrame's {key} was not captured"
    assert ours is not None, f"{case.case_id}: PyNite found no {key}"
    diff_pct = (ours - theirs) / theirs * 100.0
    message = (
        f"{case.case_id} {key}: PyNite {ours:.3f} vs PolyFrame {theirs:.3f} "
        f"({diff_pct:+.1f} %, margin ±{MARGIN_PCT:.0f} %)"
    )
    deviation = case.deviation(key)
    if abs(diff_pct) <= MARGIN_PCT:
        return
    # One-sided: a known deviation only excuses PyNite being the higher (unconservative)
    # tool, up to its bound. PyNite falling below PolyFrame beyond the margin is new.
    if deviation is not None and 0 < diff_pct <= deviation.bound_pct:
        pytest.xfail(
            f"{message}; known deviation, PyNite high by up to {deviation.bound_pct:.0f} %: "
            f"{deviation.reason}"
        )
    pytest.fail(message)


@pytest.mark.skipif(not _CASES, reason="no PolyFrame corpus records captured")
@pytest.mark.parametrize("case", [c for c in _CASES if c.compare_in_plane],
                         ids=[c.case_id for c in _CASES if c.compare_in_plane])
def test_in_plane_critical_load_factor_matches_polyframe(case):
    """``alpha_cr`` of the first in-plane mode within ±5 % of PolyFrame's."""
    theirs = _record(case)["polyframe"]["alpha_cr_ip"]
    ours = _pynite(case).alpha_cr_ip
    _compare(case, "alpha_cr_ip", ours, theirs)


@pytest.mark.skipif(not _CASES, reason="no PolyFrame corpus records captured")
@pytest.mark.parametrize("case", [c for c in _CASES if c.compare_out_of_plane],
                         ids=[c.case_id for c in _CASES if c.compare_out_of_plane])
def test_out_of_plane_critical_load_factor_matches_polyframe(case):
    """``alpha_cr,op`` of the first out-of-plane mode within ±5 % of PolyFrame's."""
    theirs = _record(case)["polyframe"]["alpha_cr_op"]
    ours = _pynite(case).alpha_cr_op(VARIANT)
    _compare(case, "alpha_cr_op", ours, theirs)


# PolyFrame's "out-of-plane buckling mode: E3, ULS1" line is not the member of the
# critical mode: on the UDL portals it names a column while the factor it prints
# equals PyNite's beam mode to 1 %. It is recorded in the corpus and shown in the
# summary document, not gated.


def test_every_case_renders_both_models():
    """Every catalogue entry builds a PolyFrame model and a PyNite model."""
    for case in all_cases():
        model = to_polyframe_model(case)
        assert model["Elements"] and model["Nodes"] and model["LoadCases"]
        pynite_model = pynite_side.build_model(case)
        assert len(pynite_model.members) == len(case.elements)


def test_tests_never_start_polyframe():
    """No test module imports the bridge runner: capture is a manual step."""
    testing = Path(__file__).resolve().parent
    offenders = []
    for path in testing.glob("test_*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                names = [node.module] + [f"{node.module}.{a.name}" for a in node.names]
            if any(n.endswith("polyframe.runner") or n.endswith("polyframe.capture") for n in names):
                offenders.append(path.name)
    assert not offenders, f"these tests import the PolyFrame bridge: {offenders}"
