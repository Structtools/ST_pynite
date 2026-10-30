"""Third opinion on the out-of-plane cases: the governing member, isolated, in the
independent Rayleigh-Ritz solver (``ltb_reference``).

For a case, the member PyNite's first out-of-plane mode lives in is cut out of the
frame with its first-order forces from a PyNite static solve (which agree with
PolyFrame's to 0.2 %), and given the end conditions both tools assume: lateral
displacement held at both ends (every node is a lateral restraint), twist held at an
end where a non-collinear member meets it (frame corner, stub), free at a support
alone (a column base), warping free. A distributed load on the member keeps its load
height. The solver integrates the true warping energy, so it needs no warping
length; its factor is compared with PolyFrame's and with PyNite's three warping
variants.

Isolating a member assumes the rest of the frame acts as rigid forks at its ends,
which the stiffness ratios support here (a column's twist is held by the beam's
weak-axis bending, 50 times stiffer than the column's torsion; a rafter's twist by
the column's weak-axis bending, 300 times stiffer). It is an approximation both
tools share, not a bias towards either.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from . import pynite_side
from .cases import FrameCase, all_cases, get_case
from .corpus import read_record
from .ltb_reference import MemberLTB, critical_load_factor
from .sections import SECTIONS, STEEL_E, STEEL_G

SAMPLES = 61


@dataclass
class ReferenceCheck:
    case_id: str
    element: int
    polyframe: float | None
    pynite_restraints: float | None
    pynite_member_length: float | None
    pynite_no_warping: float | None
    ritz: float
    twist_held: tuple[bool, bool]
    note: str = ""


def _in_plane_moment(member, x: np.ndarray, combo: str) -> np.ndarray:
    """The member's in-plane bending moment along it, whichever local axis carries it."""
    my = np.array([member.moment("My", xi, combo) for xi in x])
    mz = np.array([member.moment("Mz", xi, combo) for xi in x])
    return mz if np.abs(mz).max() >= np.abs(my).max() else my


def check_case(case: FrameCase, element: int | None = None) -> ReferenceCheck:
    try:
        record = read_record(case.case_id)
        polyframe = record["polyframe"].get("alpha_cr_op")
    except FileNotFoundError:
        polyframe = None

    ours = pynite_side.run(case)
    number = element if element is not None else ours.op_element_restraints
    e = case.element(number)
    section = SECTIONS[e.section]
    L = case.element_length(number)

    # First-order forces from a static solve of the frame
    static = pynite_side.build_model(case)
    static.analyze_linear(check_statics=False)
    member = static.members[f"E{number}"]
    x = np.linspace(0.0, L, SAMPLES)
    axial = np.array([member.axial(xi, "ULS1") for xi in x])
    compression = float(np.mean(axial))  # PyNite's member.axial(): compression positive
    moment = _in_plane_moment(member, x, "ULS1")

    a_node, b_node = case.node(e.start), case.node(e.end)
    horizontal = abs(a_node.y - b_node.y) < 1e-9
    q_value, a_height = 0.0, 0.0
    if horizontal:
        # Sagging positive: a downward load bends the member with the top in compression
        mid = moment[SAMPLES // 2]
        if abs(mid) > 1e-9:
            moment = moment * np.sign(mid)
        for load in case.line_loads:
            if load.element == number:
                q_value += load.kn_per_m * 1e3
                a_height = pynite_side._load_height(case, number, load.attack, load.kn_per_m)

    fork = pynite_side.is_collinear(case)
    held_start = pynite_side._twist_restrained(case, number, e.start, fork)
    held_end = pynite_side._twist_restrained(case, number, e.end, fork)
    phi_held = tuple(p for p, held in ((0.0, held_start), (1.0, held_end)) if held)

    ltb = MemberLTB(
        L=L, E=STEEL_E, G=STEEL_G, Iz=section.Iz, J=section.It, Iw=section.Iw,
        i0_sq=(section.Iy + section.Iz) / section.A, phi_held_at=phi_held,
    )
    M = lambda s: np.interp(s, x, moment)  # noqa: E731
    q = (lambda s: q_value * np.ones_like(s)) if q_value else None
    ritz = critical_load_factor(ltb, M=M, N=compression, q=q, a=a_height, n_v=14, n_phi=14)

    return ReferenceCheck(
        case_id=case.case_id, element=number, polyframe=polyframe,
        pynite_restraints=ours.alpha_cr_op_restraints,
        pynite_member_length=ours.alpha_cr_op_eqtorsion,
        pynite_no_warping=ours.alpha_cr_op_ignore,
        ritz=ritz, twist_held=(held_start, held_end),
        note=f"N = {compression / 1e3:.1f} kN, |M|max = {np.abs(moment).max() / 1e3:.1f} kNm"
        + (f", q = {q_value / 1e3:.0f} kN/m at a = {a_height:+.2f} m" if q_value else ""),
    )


def check_two_span() -> ReferenceCheck:
    """The two-span beam needs the whole 6 m member: v held at three points."""
    case = get_case("beam-ipe300-two-span")
    polyframe = read_record(case.case_id)["polyframe"]["alpha_cr_op"]
    ours = pynite_side.run(case)
    s = SECTIONS["IPE300"]
    L, l, q = 6.0, 3.0, 10e3
    R = 3 * q * l / 8  # end reaction of a two-equal-span beam
    M = lambda x: np.where(x <= l, R * x - q * x**2 / 2, R * (L - x) - q * (L - x) ** 2 / 2)  # noqa: E731
    ltb = MemberLTB(L=L, E=STEEL_E, G=STEEL_G, Iz=s.Iz, J=s.It, Iw=s.Iw, i0_sq=(s.Iy + s.Iz) / s.A,
                    v_held_at=(0.0, 0.5, 1.0), phi_held_at=(0.0, 1.0))
    ritz = critical_load_factor(ltb, M=M, n_v=18, n_phi=18)
    return ReferenceCheck(
        case_id=case.case_id, element=1, polyframe=polyframe,
        pynite_restraints=ours.alpha_cr_op_restraints, pynite_member_length=ours.alpha_cr_op_eqtorsion,
        pynite_no_warping=ours.alpha_cr_op_ignore, ritz=ritz, twist_held=(True, True),
        note="whole 6 m beam, lateral displacement held at the three supports, twist at the ends",
    )


DEFAULT_CASES = (
    "beam-ipe300-udl",
    "beam-ipe300-udl-top-flange",
    "beam-ipe300-midspan-point",
    "column-ipe300-pinned-axial-moment",
    "column-ipe300-cantilever",
    "portal-ipe300-pinned-100-wind",
    "portal-ipe300-fixed-100-wind",
    "portal-ipe300-pinned-300-wind",
    "portal-he200a-pinned-100-wind",
    "portal-ipe300-pinned-udl-wind",
    "portal-ipe300-fixed-udl-wind",
    "portal-he200a-fixed-udl-wind",
    "portal-ipe300-pinned-udl-top-flange",
    "gable-ipe300-pinned",
)


def run_all(case_ids=DEFAULT_CASES) -> list[ReferenceCheck]:
    out = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        for cid in case_ids:
            out.append(check_case(get_case(cid)))
        out.append(check_two_span())
    return out


def _pct(a: float | None, ref: float) -> str:
    return "-" if a is None else f"{(a - ref) / ref * 100:+.1f} %"


def format_table(checks: list[ReferenceCheck]) -> str:
    lines = [
        "| Case | Member (twist held start/end) | Ritz reference | PolyFrame | PyNite, restraint spacing | PyNite, member length | PyNite, no warping | Forces |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for c in checks:
        held = "/".join("yes" if h else "no" for h in c.twist_held)
        lines.append(
            f"| `{c.case_id}` | E{c.element} ({held}) | **{c.ritz:.3f}** "
            f"| {c.polyframe if c.polyframe is not None else '-'} ({_pct(c.polyframe, c.ritz)}) "
            f"| {c.pynite_restraints:.3f} ({_pct(c.pynite_restraints, c.ritz)}) "
            f"| {c.pynite_member_length:.3f} ({_pct(c.pynite_member_length, c.ritz)}) "
            f"| {c.pynite_no_warping:.3f} ({_pct(c.pynite_no_warping, c.ritz)}) | {c.note} |"
        )
    return "\n".join(lines)


if __name__ == "__main__":
    print(format_table(run_all()))
