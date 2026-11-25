"""Render a ``FrameCase`` as a PyNite model and run both buckling analyses.

Geometry
--------
The frame lies in PyNite's global XY plane (X right, Y up), the same plane and axes
as PolyFrame. Members are rotated by 90 degrees so that the Eurocode strong axis
``Iy`` bends in the plane of the frame, as the existing frame tests do.

Supports
--------
In plane, a support preset restrains DX / DY / RZ as PolyFrame does. Out of plane,
PolyFrame's general method restrains the lateral displacement (its "FZ") at every
node of the frame and nothing else - no twist, no out-of-plane rotation, at fixed
bases included (the report's "Stability Supports" table has an X in the FZ column
of every node and the MX / MY columns empty). The PyNite model does the same: DZ
restrained at every node, RX and RY free everywhere.

A single member (a column, a beam, a continuous beam) is a twist mechanism under
those supports, so for a case whose elements are all collinear the supported nodes
also restrain the rotation about the member axis - fork supports, the usual
assumption behind a lateral restraint, and the only one that gives a simply
supported beam a finite ``M_cr``. ``build_model`` takes ``fork`` to force either
choice, and ``out_of_plane_rotations`` to restrain RX and RY at the supports
altogether.

Loads
-----
PolyFrame has element loads only; a point load at an element end is applied here
as a nodal load so that it does not depend on the member subdivision. A vertical
load with a top- or bottom-flange attack point carries the load height ``z_g`` of
EN 1993-1-1 Annex F (positive towards the shear centre): ``+h/2`` for a gravity load
on the top flange, ``-h/2`` on the bottom flange.

Warping
-------
PyNite's ``warping='equivalent_torsion'`` replaces ``J`` by ``J + pi^2 E Iw /
(G L^2)`` with ``L`` the member length, which is exact for a member whose twist is
restrained at both ends (fork supports). PolyFrame's model restrains no twist at a
support, so a column standing on a base is twist-free there and the distance
between its twist restraints is ``2 L``; the same holds for each span of a
continuous beam at an intermediate support. Comparing both tools showed exactly
that: beam-governed cases match the equivalent-torsion value, column-governed cases
match it once the column's warping length is ``2 L``, and a two-span beam once each
span's is. ``warping_length_factor`` encodes the rule - a member end is
twist-restrained when another, non-collinear member meets it there - and the
``restraints`` variant below passes it to PyNite with ``def_warping_length`` (the
library's own per-member warping length; ``math.inf`` when neither end is held).
The rule was found by matching PolyFrame and has a closed-form core (torsional
buckling with twist held at one end, ``L_T = 2 L``, tested in
``Testing/test_buckling.py``); its use under bending is checked against the
independent Rayleigh-Ritz reference in ``ltb_reference.py``.

Analysis
--------
``analyze_buckling`` is run four times on fresh models: ``plane='XY'`` for the
in-plane factor, and ``plane='XY', out_of_plane=True`` with ``warping='ignore'``,
with ``warping='equivalent_torsion'`` on the member length, and with the
equivalent torsion constant on the twist-restraint spacing. The last is the one
compared with PolyFrame.
"""

import warnings
from dataclasses import dataclass, field

import numpy as np

from .cases import FrameCase
from .sections import SECTIONS, STEEL_E, STEEL_G, STEEL_NU

ELEMENTS_PER_MEMBER = 8

_IN_PLANE: dict[str, tuple[bool, bool, bool]] = {
    #            DX     DY     RZ
    "pinned": (True, True, False),
    "fixed": (True, True, True),
    "roller-y": (False, True, False),
    "roller-x": (True, False, False),
}


def _node_name(node_id: int) -> str:
    return f"N{node_id}"


def _member_name(number: int) -> str:
    return f"E{number}"


def _load_height(case: FrameCase, element: int, attack: str, kn: float) -> float:
    if attack == "Middle":
        return 0.0
    h = SECTIONS[case.element(element).section].h
    towards_centre = (attack == "Top") == (kn > 0)
    return 0.5 * h if towards_centre else -0.5 * h


def is_collinear(case: FrameCase) -> bool:
    """All elements on one straight line (a column or a beam, possibly in pieces)."""
    first = case.element(1)
    a, b = case.node(first.start), case.node(first.end)
    for element in case.elements:
        for node_id in (element.start, element.end):
            c = case.node(node_id)
            cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
            if abs(cross) > 1e-9:
                return False
    return True


def _twist_dofs(case: FrameCase, node_id: int) -> tuple[bool, bool]:
    """(RX, RY) to restrain for a fork support at a node: the axis of its members."""
    rx = ry = False
    for element in case.elements:
        if node_id not in (element.start, element.end):
            continue
        a, b = case.node(element.start), case.node(element.end)
        if abs(a.y - b.y) < 1e-9:
            rx = True  # horizontal member: twist is RX
        elif abs(a.x - b.x) < 1e-9:
            ry = True  # vertical member: twist is RY
        else:
            rx = ry = True
    return rx, ry


def _collinear(case: FrameCase, a: int, b: int) -> bool:
    ea, eb = case.element(a), case.element(b)
    p, q = case.node(ea.start), case.node(ea.end)
    r, s = case.node(eb.start), case.node(eb.end)
    dx1, dy1 = q.x - p.x, q.y - p.y
    dx2, dy2 = s.x - r.x, s.y - r.y
    return abs(dx1 * dy2 - dy1 * dx2) < 1e-9 * max(1.0, abs(dx1 * dx2 + dy1 * dy2))


def _twist_restrained(case: FrameCase, element: int, node_id: int, fork: bool) -> bool:
    """Is the twist of ``element`` held at ``node_id``?

    By another member meeting it at an angle (a frame corner, a stub), or by a fork
    support when ``fork`` is on. A support alone holds no twist in PolyFrame's model.
    """
    node = case.node(node_id)
    if fork and node.support:
        return True
    for other in case.elements:
        if other.number == element or node_id not in (other.start, other.end):
            continue
        if not _collinear(case, element, other.number):
            return True
    return False


def warping_length_factor(case: FrameCase, element: int, fork: bool | None = None) -> float | None:
    """Twist-restraint spacing over member length: 1, 2, or ``None`` for unrestrained.

    1 when both ends are twist-restrained, 2 when one is (the member warps as the
    half of a member twice as long), ``None`` when neither is, in which case the
    warping stiffness has no restraint to work against and is dropped.
    """
    if fork is None:
        fork = is_collinear(case)
    e = case.element(element)
    held = [_twist_restrained(case, element, n, fork) for n in (e.start, e.end)]
    if all(held):
        return 1.0
    if any(held):
        return 2.0
    return None


def build_model(case: FrameCase, out_of_plane_rotations: bool = False,
                fork: bool | None = None, warping_by_restraints: bool = False):
    """The PyNite ``FEModel3D`` of a case, ready for ``analyze_buckling``.

    :param out_of_plane_rotations: restrain RX and RY at the supported nodes as
        well. ``False`` (default) mirrors PolyFrame's stability supports.
    :param fork: restrain the twist about the member axis at supported nodes.
        ``None`` (default) does so exactly when the case is a single straight
        member, where it is needed for stability.
    :param warping_by_restraints: declare each member's warping length to PyNite
        (``def_warping_length``) as ``warping_length_factor`` times its length, so
        that ``warping='equivalent_torsion'`` works with the twist-restraint spacing
        instead of the member length.
    """
    import math

    from Pynite import FEModel3D

    if fork is None:
        fork = is_collinear(case)

    model = FEModel3D()
    model.add_material("S355", STEEL_E, STEEL_G, STEEL_NU, 0.0)
    for name in dict.fromkeys(e.section for e in case.elements):
        s = SECTIONS[name]
        model.add_section(name, s.A, s.Iy, s.Iz, s.It, Asy=s.Avy, Asz=s.Avz, Iw=s.Iw)

    for node in case.nodes:
        model.add_node(_node_name(node.id), node.x, node.y, 0.0)
        if node.support:
            dx, dy, rz = _IN_PLANE[node.support]
            rot = out_of_plane_rotations
            rx, ry = _twist_dofs(case, node.id) if fork else (False, False)
            model.def_support(_node_name(node.id), dx, dy, True, rot or rx, rot or ry, rz)
        else:
            model.def_support(_node_name(node.id), False, False, True, False, False, False)

    for element in case.elements:
        model.add_member(
            _member_name(element.number),
            _node_name(element.start),
            _node_name(element.end),
            "S355",
            element.section,
            rotation=90,
        )
        if warping_by_restraints:
            factor = warping_length_factor(case, element.number, fork)
            length = math.inf if factor is None else factor * case.element_length(element.number)
            model.def_warping_length(_member_name(element.number), length)

    for load in case.point_loads:
        direction = "FY" if load.direction == "down" else "FX"
        value = -load.kn * 1e3 if load.direction == "down" else load.kn * 1e3
        length = case.element_length(load.element)
        element = case.element(load.element)
        z_g = _load_height(case, load.element, load.attack, load.kn)
        if load.at <= 1e-9 and z_g == 0.0:
            model.add_node_load(_node_name(element.start), direction, value, case=load.load_case)
        elif abs(load.at - length) <= 1e-9 and z_g == 0.0:
            model.add_node_load(_node_name(element.end), direction, value, case=load.load_case)
        else:
            model.add_member_pt_load(
                _member_name(load.element), direction, value, load.at,
                case=load.load_case, load_height=z_g,
            )

    for load in case.line_loads:
        direction = "FY" if load.direction == "down" else "FX"
        value = -load.kn_per_m * 1e3 if load.direction == "down" else load.kn_per_m * 1e3
        z_g = _load_height(case, load.element, load.attack, load.kn_per_m)
        model.add_member_dist_load(
            _member_name(load.element), direction, value, value,
            case=load.load_case, load_height=z_g,
        )

    model.add_load_combo("ULS1", {name: 1.0 for name in case.load_cases()})
    return model


@dataclass
class PyniteBuckling:
    """Both factors of one case, plus what identifies the governing member."""

    alpha_cr_ip: float | None
    alpha_cr_op_ignore: float | None
    alpha_cr_op_eqtorsion: float | None
    alpha_cr_op_restraints: float | None
    op_element_ignore: int | None
    op_element_eqtorsion: int | None
    op_element_restraints: int | None
    axial_kn: dict[int, float] = field(default_factory=dict)
    in_plane_modes: list[float] = field(default_factory=list)
    out_of_plane_modes_ignore: list[float] = field(default_factory=list)
    out_of_plane_modes_eqtorsion: list[float] = field(default_factory=list)
    out_of_plane_modes_restraints: list[float] = field(default_factory=list)

    VARIANTS = ("ignore", "equivalent_torsion", "restraints")

    def alpha_cr_op(self, variant: str) -> float | None:
        return {
            "ignore": self.alpha_cr_op_ignore,
            "equivalent_torsion": self.alpha_cr_op_eqtorsion,
            "restraints": self.alpha_cr_op_restraints,
        }[variant]

    def op_element(self, variant: str) -> int | None:
        return {
            "ignore": self.op_element_ignore,
            "equivalent_torsion": self.op_element_eqtorsion,
            "restraints": self.op_element_restraints,
        }[variant]


def _governing_member(results, dof: int = 2) -> int | None:
    """The physical member carrying the largest out-of-plane displacement of mode 1.

    ``dof = 2`` is DZ. Rows of the mode shape belong to nodes of the analysis
    model; each member's sub-members name the nodes they span.
    """
    if results.mode_count == 0:
        return None
    model = results.analysis_model
    mode = results.mode_shapes[:, 0]
    rows: dict[str, float] = {}
    for (node, name), value in zip(results.dof_map, mode):
        if name == "DZ":
            rows[node] = rows.get(node, 0.0) + float(value) ** 2
    best, best_value = None, -1.0
    for member_name, member in model.members.items():
        nodes = set()
        for sub in member.sub_members.values():
            nodes.add(sub.i_node.name)
            nodes.add(sub.j_node.name)
        value = sum(rows.get(n, 0.0) for n in nodes)
        if value > best_value:
            best, best_value = member_name, value
    return int(best[1:]) if best else None


def run(case: FrameCase, num_modes: int = 3, elements_per_member: int = ELEMENTS_PER_MEMBER,
        out_of_plane_rotations: bool = False, fork: bool | None = None) -> PyniteBuckling:
    """In-plane and out-of-plane critical load factors of a case."""
    ip_modes: list[float] = []
    alpha_ip = None
    if case.compare_in_plane:
        model = build_model(case, out_of_plane_rotations, fork)
        try:
            res = model.analyze_buckling(
                combo_name="ULS1", num_modes=num_modes, plane="XY",
                elements_per_member=elements_per_member,
            )
            ip_modes = [float(v) for v in res.load_multipliers]
            alpha_ip = ip_modes[0] if ip_modes else None
        except (ValueError, RuntimeError):
            alpha_ip = None

    out: dict[str, tuple[float | None, int | None, list[float]]] = {}
    for variant in PyniteBuckling.VARIANTS:
        by_restraints = variant == "restraints"
        model = build_model(case, out_of_plane_rotations, fork, warping_by_restraints=by_restraints)
        with warnings.catch_warnings():
            # The member-length variant is computed on purpose, to report what the
            # default gives; PyNite's warning about it is the point, not noise here.
            if variant == "equivalent_torsion":
                warnings.simplefilter("ignore", UserWarning)
            res = model.analyze_buckling(
                combo_name="ULS1", num_modes=num_modes, plane="XY", out_of_plane=True,
                elements_per_member=elements_per_member,
                warping="ignore" if variant == "ignore" else "equivalent_torsion",
            )
        modes = [float(v) for v in res.load_multipliers]
        out[variant] = (modes[0] if modes else None, _governing_member(res), modes)

    # First-order axial force per element (compression negative), from the last run
    axial: dict[int, float] = {}
    static = build_model(case, out_of_plane_rotations, fork)
    static.analyze_linear(check_statics=False)
    for element in case.elements:
        member = static.members[_member_name(element.number)]
        axial[element.number] = float(np.mean([
            member.axial(x, "ULS1") for x in (0.0, member.L())
        ])) / 1e3

    return PyniteBuckling(
        alpha_cr_ip=alpha_ip,
        alpha_cr_op_ignore=out["ignore"][0],
        alpha_cr_op_eqtorsion=out["equivalent_torsion"][0],
        alpha_cr_op_restraints=out["restraints"][0],
        op_element_ignore=out["ignore"][1],
        op_element_eqtorsion=out["equivalent_torsion"][1],
        op_element_restraints=out["restraints"][1],
        axial_kn=axial,
        in_plane_modes=ip_modes,
        out_of_plane_modes_ignore=out["ignore"][2],
        out_of_plane_modes_eqtorsion=out["equivalent_torsion"][2],
        out_of_plane_modes_restraints=out["restraints"][2],
    )
