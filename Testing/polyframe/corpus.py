"""The captured PolyFrame records: one JSON file per case in ``corpus/``.

A record holds what PolyFrame said and where it came from, never a PyNite value:

    {
      "case_id": ...,
      "provenance": {"application_version", "connector_version", "captured_on",
                     "model_sha256", "method": "GeneralMethod"},
      "polyframe": {"alpha_cr_ip", "alpha_cr_ip_source": "text" | "ocr" | "hand",
                    "alpha_cr_op", "alpha_cr_op_source": "ocr" | "hand",
                    "alpha_ult_k", "governing_element_op",
                    "util_stability", "util_cross_section", "uls_utilization",
                    "axial_kn": {"1": ..}, "moment_knm": {"1": ..}},
      "messages": [...]
    }

``hand_read.json`` beside the records holds values typed in from a report by a
person when OCR could not read them; the capture merges them with
``alpha_cr_*_source = "hand"`` so that the origin stays visible.
"""

import json
from pathlib import Path
from typing import Any

CORPUS_DIR = Path(__file__).resolve().parent / "corpus"
REPORTS_DIR = CORPUS_DIR / "reports"
HAND_READ = CORPUS_DIR / "hand_read.json"


def record_path(case_id: str) -> Path:
    return CORPUS_DIR / f"{case_id}.json"


def read_record(case_id: str) -> dict[str, Any]:
    path = record_path(case_id)
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def write_record(record: dict[str, Any]) -> Path:
    CORPUS_DIR.mkdir(parents=True, exist_ok=True)
    path = record_path(record["case_id"])
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def hand_read() -> dict[str, dict[str, float]]:
    if not HAND_READ.exists():
        return {}
    return json.loads(HAND_READ.read_text(encoding="utf-8"))


def captured_case_ids() -> list[str]:
    if not CORPUS_DIR.exists():
        return []
    return sorted(p.stem for p in CORPUS_DIR.glob("*.json") if p.name != HAND_READ.name)
