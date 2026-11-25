"""Turn one bridge outcome and its PDF into a corpus record."""

import datetime as dt
from pathlib import Path
from typing import Any

from .cases import FrameCase
from .corpus import hand_read
from .pdf_report import ReportValues, read_report
from .runner import CaseOutcome
from .translate import model_hash, to_polyframe_model


class CaptureError(RuntimeError):
    """A case produced no usable record."""


def _api_quantities(result: dict[str, Any] | None) -> dict[str, Any]:
    """Per-element axial force and moment extremes of ULS1, in kN and kNm."""
    out: dict[str, Any] = {"axial_kn": {}, "moment_knm": {}, "uls_utilization": None}
    if not result:
        return out
    out["uls_utilization"] = result.get("ULSUtilization")
    for combination in result.get("LoadCombinationResults") or []:
        if combination.get("Name") != "ULS1":
            continue
        for element in combination.get("FirstOrderSectionalForces") or []:
            rows = element.get("SectionalForces") or []
            number = str(element.get("Number"))
            n = [r["AxialForce"] for r in rows if r.get("AxialForce") is not None]
            m = [
                r["StrongAxisBendingMoment"]
                for r in rows
                if r.get("StrongAxisBendingMoment") is not None
            ]
            if n:
                # compression negative; the value of largest magnitude
                out["axial_kn"][number] = round(min(n, key=lambda v: -abs(v)) / 1e3, 4)
            if m:
                out["moment_knm"][number] = round(max(m, key=abs) / 1e3, 4)
    return out


def _merge_hand_read(case_id: str, values: ReportValues, polyframe: dict[str, Any]) -> None:
    hand = hand_read().get(case_id, {})
    for key in ("alpha_cr_ip", "alpha_cr_op"):
        if polyframe.get(key) is None and key in hand:
            polyframe[key] = hand[key]
            polyframe[f"{key}_source"] = "hand"
        elif key in hand and polyframe.get(key) is not None:
            if abs(hand[key] - polyframe[key]) > 0.011 * max(1.0, abs(hand[key])):
                values.messages.append(
                    f"hand-read {key} = {hand[key]} disagrees with the parsed "
                    f"{polyframe[key]} ({polyframe[f'{key}_source']})"
                )


def record_from_report(
    case: FrameCase,
    pdf_path: Path,
    *,
    result: dict[str, Any] | None,
    connector_version: str,
    captured_on: str | None = None,
    ocr: bool = True,
) -> dict[str, Any]:
    values = read_report(pdf_path, ocr=ocr)
    polyframe: dict[str, Any] = {
        "alpha_cr_ip": values.alpha_cr_ip,
        "alpha_cr_ip_source": (
            "text" if values.alpha_cr_ip_text is not None
            else "ocr" if values.alpha_cr_ip_ocr is not None else None
        ),
        "alpha_cr_op": values.alpha_cr_op,
        "alpha_cr_op_source": "ocr" if values.alpha_cr_op is not None else None,
        "alpha_cr_op_formula_ocr": values.alpha_cr_op_formula,
        "alpha_cr_op_plot_ocr": values.alpha_cr_op_plots,
        "alpha_ult_k": values.alpha_ult_k,
        "governing_element_op": values.governing_element_op,
        "util_stability": values.util_stability,
        "util_cross_section": values.util_cross_section,
        **_api_quantities(result),
    }
    _merge_hand_read(case.case_id, values, polyframe)

    if case.compare_out_of_plane and polyframe["alpha_cr_op"] is None:
        values.messages.append("no alpha_cr,op could be read; add it to hand_read.json")
    if case.compare_in_plane and polyframe["alpha_cr_ip"] is None:
        values.messages.append("no in-plane alpha_cr could be read; add it to hand_read.json")

    return {
        "case_id": case.case_id,
        "title": case.title,
        "provenance": {
            "kind": "polyframe_run",
            "method": "GeneralMethod",
            "application_version": values.application_version,
            "connector_version": connector_version,
            "captured_on": captured_on or dt.date.today().isoformat(),
            "model_sha256": model_hash(to_polyframe_model(case)),
            "report": pdf_path.name,
            "ocr": values.ocr_used,
        },
        "polyframe": polyframe,
        "messages": values.messages,
    }


def record_from_outcome(
    case: FrameCase, outcome: CaseOutcome, *, connector_version: str, ocr: bool = True
) -> dict[str, Any]:
    if not outcome.pdf_path:
        raise CaptureError(
            f"{case.case_id}: no PDF came back ({outcome.error or 'SavePDF returned false'}); "
            "the critical load factors exist only in the report"
        )
    if not outcome.calculated:
        # The PDF alone carries everything compared; note the missing forces.
        pass
    record = record_from_report(
        case, Path(outcome.pdf_path), result=outcome.result,
        connector_version=connector_version, ocr=ocr,
    )
    if outcome.error:
        record["messages"].append(f"bridge: {outcome.error}")
    return record
