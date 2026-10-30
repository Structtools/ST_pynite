"""Render a ``FrameCase`` as a PolyFrame data model (JSON-serialisable dict).

The bridge deserialises the dict straight into ``PolyFrameDataModel`` with
case-insensitive property names and a string-enum converter, so property names are
PascalCase and enums are their C# names. Units are SI base: metres and newtons.

Conventions measured on PolyFrame 1.1.11 (structural_backend,
``docs/POLYFRAME_DISCOVERY.md``): the plane is X right, Y up; a positive
``GlobalY`` load acts **downward**, a positive ``GlobalX`` load acts along +X; the
supports are ``GlobalXDeflectionSupported`` / ``GlobalYDeflectionSupported`` /
``GlobalZRotationSupported``. Out-of-plane, the general method restrains the lateral
displacement at every node and prints that as "Stability Supports FZ" in the
report; it is not an input.
"""

import hashlib
import json
from typing import Any

from .cases import FrameCase, LineLoad, PointLoad

POLYFRAME_MODEL_VERSION = 8

_END_TOLERANCE_M = 1e-9

_RESTRAINTS: dict[str, tuple[bool, bool, bool]] = {
    #            DX     DY     RZ
    "pinned": (True, True, False),
    "fixed": (True, True, True),
    "roller-y": (False, True, False),
    "roller-x": (True, False, False),
}

_DIRECTION = {"down": "GlobalY", "X": "GlobalX"}

_LOAD_GROUP = {"G": "PermanentLoad", "W": "WindLoad"}
_DURATION = {"G": "Permanent", "W": "ShortTerm"}


def kn_to_n(value: float) -> float:
    return value * 1e3


def _node(case: FrameCase, node_id: int) -> dict[str, Any]:
    node = case.node(node_id)
    out: dict[str, Any] = {
        "Number": node.id,
        "XCoordinate": node.x,
        "YCoordinate": node.y,
    }
    if node.support:
        dx, dy, rz = _RESTRAINTS[node.support]
        out.update(
            GlobalXDeflectionSupported=dx,
            GlobalYDeflectionSupported=dy,
            GlobalZRotationSupported=rz,
        )
    return out


def _section_numbers(case: FrameCase) -> dict[str, int]:
    names = list(dict.fromkeys(e.section for e in case.elements))
    return {name: i for i, name in enumerate(names, start=1)}


def _point_load(case: FrameCase, load: PointLoad) -> dict[str, Any]:
    out: dict[str, Any] = {
        "ElementsString": str(load.element),
        "LoadAttackPoint": load.attack,
        "LoadDirection": _DIRECTION[load.direction],
        "Force": kn_to_n(load.kn),
    }
    length = case.element_length(load.element)
    if load.at <= _END_TOLERANCE_M:
        out["ElementStartLoad"] = True
    elif abs(load.at - length) <= _END_TOLERANCE_M:
        out["ElementEndLoad"] = True
    else:
        out["Location"] = load.at
    return out


def _line_load(case: FrameCase, load: LineLoad) -> dict[str, Any]:
    return {
        "ElementsString": str(load.element),
        "LoadAttackPoint": load.attack,
        "StartLocation": 0.0,
        "EndLocation": case.element_length(load.element),
        "ForceStart": kn_to_n(load.kn_per_m),
        "ForceEnd": kn_to_n(load.kn_per_m),
        "LoadDirection": _DIRECTION[load.direction],
        "IsSelfWheight": False,
    }


def to_polyframe_model(case: FrameCase) -> dict[str, Any]:
    sections = _section_numbers(case)
    load_cases: list[dict[str, Any]] = []
    for name in case.load_cases():
        load_case: dict[str, Any] = {
            "Name": name,
            "LoadGroup": _LOAD_GROUP[name],
            "LoadDuration": _DURATION[name],
            "IncludeSelfWeight": False,
        }
        point = [_point_load(case, l) for l in case.point_loads if l.load_case == name]
        line = [_line_load(case, l) for l in case.line_loads if l.load_case == name]
        if point:
            load_case["PointLoads"] = point
        if line:
            load_case["DistributedLoads"] = line
        load_cases.append(load_case)

    names = case.load_cases()
    return {
        "Version": POLYFRAME_MODEL_VERSION,
        "DesignSettings": {
            "DesignCode": "EurocodeDKNA",
            "ConsequenceClass": "CC2",
            "BucklingDesignMethod": "GeneralMethod",
        },
        "DocumentationSettings": {"DocumentTitle": case.title},
        "Nodes": [_node(case, n.id) for n in case.nodes],
        "CrossSections": [
            {
                "Number": number,
                "CrossSectionName": name,
                "CrossSectionType": "Steel",
                "MaterialName": "S355",
            }
            for name, number in sections.items()
        ],
        "Elements": [
            {
                "Number": e.number,
                "NodeNumbers": [e.start, e.end],
                "ElementType": "SteelElement",
                "SteelElementProperties": {
                    "SectionNumber": sections[e.section],
                    "SectionName": e.section,
                },
            }
            for e in case.elements
        ],
        "LoadGroups": [
            {"Name": _LOAD_GROUP[name], "Type": _LOAD_GROUP[name]} for name in names
        ],
        "LoadCases": load_cases,
        "ULSLoadCombinations": [
            {
                "Name": "ULS1",
                "LoadCasesAndFactors": [
                    {"LoadCase": name, "Factor": 1.0} for name in names
                ],
                "LoadDuration": "ShortTerm" if "W" in names else "Permanent",
            }
        ],
        "SLSLoadCombinations": [],
        "ALSLoadCombinations": [],
        "FireLoadCombinations": [],
    }


def model_hash(model: dict[str, Any]) -> str:
    """A stable hash of the model a record was captured from."""
    text = json.dumps(model, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()
