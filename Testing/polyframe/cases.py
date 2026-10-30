"""The PolyFrame buckling case catalogue.

One entry is one planar steel frame, described once and rendered into a PolyFrame
model (``translate.to_polyframe_model``) and a PyNite model
(``pynite_side.build_model``). Every case is S355, DK national annex, PolyFrame's
general method (EN 1993-1-1 6.3.4), with one ULS combination ``1.0 G + 1.0 W`` so
that a critical load factor means the same on both sides.

The cases are grouped so that a disagreement points at one thing:

* **columns** - one vertical element. Pure axial load isolates the weak-axis
  flexural factor; a horizontal load adds in-plane bending and therefore the
  interaction between flexural and lateral-torsional buckling.
* **beams** - one horizontal element with no axial force, so the out-of-plane
  factor is ``M_cr / M_Ed`` alone. The load attack point (top / middle / bottom
  flange) probes the load-height term of the geometric stiffness.

PolyFrame's general method refuses a structure whose elements are all collinear
("Structural design is not available because of warnings" - measured on a column,
a cantilever, a simply supported beam, a two-span beam and a column drawn as two
elements), so every column and beam carries an **unloaded 1 m stub** at a right
angle to make it a frame. The stub carries no force. Its tip is a node, and the
general method restrains the lateral displacement of every node, so the stub also
acts as a twist restraint (a fork) at the member end it hangs from - in both tools,
which model the same geometry. For the same reason a free end cannot exist out of
plane in PolyFrame's general method: every node is laterally held.
* **portals** - the frames of the earlier PolyFrame harness (fixed and pinned
  bases, two sections, with and without wind), a symmetric one without any
  bending, a gable and a two-bay frame.

Plane: x to the right, y up (PolyFrame's X and Y). Loads are stated in physical
words - ``"down"`` acts downward, ``"X"`` acts along +x - and each renderer owns
its own sign convention. Loads are light enough that the in-plane factor is above
10 wherever possible, because PolyFrame then prints it on the text layer of the
report; below 10 it is read from an image by OCR.
"""

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Literal

Support = Literal["pinned", "fixed", "roller-y", "roller-x"]
Direction = Literal["down", "X"]
Attack = Literal["Top", "Middle", "Bottom"]

HEIGHT_M = 4.0
SPAN_M = 6.0


@dataclass(frozen=True)
class Node:
    id: int
    x: float
    y: float
    support: Support | None = None


@dataclass(frozen=True)
class Element:
    number: int
    start: int
    end: int
    section: str = "IPE300"


@dataclass(frozen=True)
class PointLoad:
    element: int
    at: float
    kn: float
    direction: Direction = "down"
    load_case: str = "G"
    attack: Attack = "Middle"


@dataclass(frozen=True)
class LineLoad:
    element: int
    kn_per_m: float
    direction: Direction = "down"
    load_case: str = "G"
    attack: Attack = "Middle"


@dataclass(frozen=True)
class KnownDeviation:
    """A characterised disagreement that still gates.

    ``key`` is ``"alpha_cr_ip"`` or ``"alpha_cr_op"``. One-sided: the test marks the
    case as an expected failure while PyNite is *above* PolyFrame by between the 5 %
    margin and ``bound_pct`` (the unconservative direction every recorded deviation
    has), and fails outright beyond the bound or when PyNite falls below PolyFrame by
    more than the margin. A bound is a watch on a disagreement whose cause is stated,
    not a widened margin; ``bound_pct`` of 100 means the case is a watch only.
    """

    key: str
    bound_pct: float
    reason: str


@dataclass(frozen=True)
class FrameCase:
    case_id: str
    title: str
    nodes: tuple[Node, ...]
    elements: tuple[Element, ...]
    point_loads: tuple[PointLoad, ...] = ()
    line_loads: tuple[LineLoad, ...] = ()
    compare_in_plane: bool = True
    compare_out_of_plane: bool = True
    known_deviations: tuple[KnownDeviation, ...] = ()
    notes: str = ""
    group: str = "portal"

    def node(self, node_id: int) -> Node:
        for node in self.nodes:
            if node.id == node_id:
                return node
        raise KeyError(f"{self.case_id}: no node {node_id}")

    def element(self, number: int) -> Element:
        for element in self.elements:
            if element.number == number:
                return element
        raise KeyError(f"{self.case_id}: no element {number}")

    def element_length(self, number: int) -> float:
        element = self.element(number)
        a, b = self.node(element.start), self.node(element.end)
        return ((b.x - a.x) ** 2 + (b.y - a.y) ** 2) ** 0.5

    def load_cases(self) -> tuple[str, ...]:
        names = {load.load_case for load in (*self.point_loads, *self.line_loads)}
        return tuple(name for name in ("G", "W") if name in names)

    def deviation(self, key: str) -> KnownDeviation | None:
        for deviation in self.known_deviations:
            if deviation.key == key:
                return deviation
        return None

    def __post_init__(self) -> None:
        ids = [n.id for n in self.nodes]
        if sorted(ids) != list(range(1, len(ids) + 1)):
            raise ValueError(f"{self.case_id}: nodes must be numbered 1..n")
        numbers = [e.number for e in self.elements]
        if sorted(numbers) != list(range(1, len(numbers) + 1)):
            raise ValueError(f"{self.case_id}: elements must be numbered 1..n")
        if not any(n.support for n in self.nodes):
            raise ValueError(f"{self.case_id}: no supports")
        for load in self.point_loads:
            if load.at < 0 or load.at > self.element_length(load.element) + 1e-9:
                raise ValueError(f"{self.case_id}: point load off element {load.element}")
        for load in (*self.point_loads, *self.line_loads):
            if load.attack != "Middle" and load.direction != "down":
                raise ValueError(
                    f"{self.case_id}: a top/bottom flange attack point is only "
                    "defined for a vertical load"
                )
        for load in self.line_loads:
            element = self.element(load.element)
            a, b = self.node(element.start), self.node(element.end)
            if load.direction == "down" and abs(a.y - b.y) > 1e-9:
                raise ValueError(
                    f"{self.case_id}: a 'down' line load needs a horizontal element"
                )
            if load.direction == "X" and abs(a.x - b.x) > 1e-9:
                raise ValueError(
                    f"{self.case_id}: an 'X' line load needs a vertical element"
                )


# --- building blocks -------------------------------------------------------


STUB_M = 1.0


def _column_nodes(base: Support, top: Support | None, height: float = HEIGHT_M):
    """Base, top, and the tip of a horizontal stub at the top."""
    return (
        Node(1, 0.0, 0.0, base),
        Node(2, 0.0, height, top),
        Node(3, STUB_M, height),
    )


def _column_elements(section: str = "IPE300"):
    return (Element(1, 1, 2, section), Element(2, 2, 3, section))


def _beam_nodes(span: float = SPAN_M):
    """Two supports and the tips of a vertical stub above each."""
    return (
        Node(1, 0.0, 0.0, "pinned"),
        Node(2, span, 0.0, "roller-y"),
        Node(3, 0.0, STUB_M),
        Node(4, span, STUB_M),
    )


def _beam_elements(section: str = "IPE300"):
    return (Element(1, 1, 2, section), Element(2, 1, 3, section), Element(3, 2, 4, section))


def _portal_nodes(base: Support):
    return (
        Node(1, 0.0, 0.0, base),
        Node(2, 0.0, HEIGHT_M),
        Node(3, SPAN_M, HEIGHT_M),
        Node(4, SPAN_M, 0.0, base),
    )


def _portal_elements(section: str = "IPE300"):
    """Left column 1 (base to eaves), beam 2, right column 3 (base to eaves)."""
    return (
        Element(1, 1, 2, section),
        Element(2, 2, 3, section),
        Element(3, 4, 3, section),
    )


def _eaves_loads(kn: float, wind_kn: float = 5.0):
    loads = [
        PointLoad(1, HEIGHT_M, kn, "down", "G"),
        PointLoad(3, HEIGHT_M, kn, "down", "G"),
    ]
    if wind_kn:
        loads.append(PointLoad(1, HEIGHT_M, wind_kn, "X", "W"))
    return tuple(loads)


_FIXED_BASE_REASON = (
    "PolyFrame's out-of-plane factor on a FIXED-base portal is 6-8 % below PyNite's "
    "(measured 2026-10-02: PyNite +7.3 %, +7.6 %, +6.5 % on the three fixed-base frames; "
    "every pinned-base frame agrees within 2.1 %). Fixing the bases changes nothing out "
    "of plane in either model - PolyFrame holds no twist at a fixed base and PyNite "
    "mirrors that - so the difference sits in the in-plane moment distribution, and "
    "specifically in the FIXING MOMENT AT A TWIST-FREE NODE: on the pinned frames, where "
    "the column base moment is zero, the independent Rayleigh-Ritz Vlasov solver "
    "(ltb_reference.py) reproduces PolyFrame to 0.1-0.6 % and PyNite to 0.8-2.7 %; on "
    "the fixed-base column it gives 6.46 against PyNite 6.31 and PolyFrame 5.88, and on "
    "the rafters 4.52 / 6.87 against PyNite 4.08 / 6.39 and PolyFrame 3.79 / 6.00. The "
    "second-order work of a moment applied at an end that is free to twist depends on "
    "how that moment is modelled (PyNite: nodal moments do no second-order work; the "
    "reference: the internal moment field works everywhere; PolyFrame: unknown), and a "
    "fixed base that holds no twist is PolyFrame's model, not a real detail - a real "
    "fixed base holds the twist and the question disappears. DECISION 2026-10-02: "
    "accepted as safe-side on both sides. A real fixed base holds the twist, and PyNite "
    "with the twist held at the bases gives 15.1 / 4.50 / 6.92 against the twist-free "
    "6.31 / 4.08 / 6.39, so the twist-free model both tools use is conservative for a "
    "real fixed base; both tools also sit below the rigid-fork member reference on the "
    "rafters. The remaining 6-8 % is PolyFrame being the more conservative of two safe "
    "values. The bound stays just above the measured difference as a watch."
)

_TOP_FLANGE_FRAME_REASON = (
    "PolyFrame's load-height penalty inside a frame is not the one it applies to a "
    "simply supported beam. On `beam-ipe300-udl-*` both tools agree: a top-flange load "
    "lowers alpha_cr,op by 22 % (PolyFrame) / 23 % (PyNite) and a bottom-flange load "
    "raises it by 27 % / 30 %. On this portal PolyFrame drops from 3.99 to 1.44 (-64 %) "
    "for the top flange while a bottom-flange load raises it only to 4.43 (+11 %, "
    "probe of 2026-10-02); PyNite gives 4.03 -> 2.64 (-34 %). ADJUDICATED 2026-10-02 by "
    "the independent Rayleigh-Ritz Vlasov solver (ltb_reference.py) on the rafter with "
    "its first-order moments, 10 kN/m at +0.15 m and fork ends: 2.655, within 0.5 % of "
    "PyNite's 2.64; its top-flange / shear-centre ratio is 0.63 against PyNite's 0.66 "
    "and PolyFrame's 0.36. PolyFrame is the outlier, consistent with a C2-type "
    "load-height term evaluated with the frame beam's end moments. This is a finding "
    "about PolyFrame; the bound only keeps the case under watch."
)

_LTB_NOTE = (
    "Out-of-plane only: there is no axial force, so alpha_cr,op = M_cr / M_Ed and "
    "the in-plane factor is not defined."
)


def _catalogue() -> list[FrameCase]:
    cases: list[FrameCase] = []

    # -- columns -----------------------------------------------------------
    cases += [
        FrameCase(
            case_id="column-ipe300-pinned-axial",
            title="IPE300 column 4 m, pinned base, roller top, 100 kN axial",
            group="column",
            nodes=_column_nodes("pinned", "roller-x"),
            elements=_column_elements(),
            point_loads=(PointLoad(1, HEIGHT_M, 100.0),),
            notes=(
                "Pure compression. In plane: Euler strong axis, alpha_cr = N_cr,y / N = "
                "108. Out of plane: Euler weak axis, 782 kN / 100 kN = 7.82, with no "
                "bending to couple with. Anchors both tools' flexural terms."
            ),
        ),
        FrameCase(
            case_id="column-he200a-pinned-axial",
            title="HE200A column 4 m, pinned base, roller top, 100 kN axial",
            group="column",
            nodes=_column_nodes("pinned", "roller-x"),
            elements=_column_elements("HE200A"),
            point_loads=(PointLoad(1, HEIGHT_M, 100.0),),
            notes="The stocky-weak-axis twin of the case above (alpha_cr,op about 17.3).",
        ),
        FrameCase(
            case_id="column-ipe300-pinned-axial-moment",
            title="IPE300 column 4 m, pinned base, roller top, 100 kN axial + 5 kN at mid-height",
            group="column",
            nodes=_column_nodes("pinned", "roller-x"),
            elements=_column_elements(),
            point_loads=(
                PointLoad(1, HEIGHT_M, 100.0),
                PointLoad(1, HEIGHT_M / 2, 5.0, "X", "W"),
            ),
            notes=(
                "Compression plus a 5 kNm mid-height moment: the out-of-plane factor "
                "couples weak-axis flexural and lateral-torsional buckling."
            ),
        ),
        FrameCase(
            case_id="column-ipe300-cantilever",
            title="IPE300 cantilever column 4 m, fixed base, 50 kN axial + 5 kN horizontal at the head",
            group="column",
            nodes=_column_nodes("fixed", None),
            elements=_column_elements(),
            point_loads=(
                PointLoad(1, HEIGHT_M, 50.0),
                PointLoad(1, HEIGHT_M, 5.0, "X", "W"),
            ),
            notes=(
                "A free end in plane: L_cr = 2 L, alpha_cr = 54. Out of plane the "
                "general method holds the head laterally like every node, so the "
                "column is pinned-pinned there (alpha_cr,op about 15.6, less for the "
                "bending): the earlier harness's finding that PolyFrame takes L for a "
                "cantilever's lateral-torsional length is built into the method."
            ),
        ),
    ]

    # -- beams -------------------------------------------------------------
    cases += [
        FrameCase(
            case_id="beam-ipe300-udl",
            title="IPE300 beam 6 m, simply supported, 5 kN/m at the shear centre",
            group="beam",
            nodes=_beam_nodes(),
            elements=_beam_elements(),
            line_loads=(LineLoad(1, 5.0),),
            compare_in_plane=False,
            notes=_LTB_NOTE + (
                " M_Ed = 22.5 kNm. With warping M_cr (C1 = 1.13) is about 102 kNm "
                "(alpha 4.5); without warping about 85 kNm (alpha 3.8), so this case "
                "tells whether PolyFrame includes the warping stiffness."
            ),
        ),
        FrameCase(
            case_id="beam-ipe300-udl-top-flange",
            title="IPE300 beam 6 m, simply supported, 5 kN/m on the top flange",
            group="beam",
            nodes=_beam_nodes(),
            elements=_beam_elements(),
            line_loads=(LineLoad(1, 5.0, attack="Top"),),
            compare_in_plane=False,
            notes=_LTB_NOTE + " Destabilising load height z_g = +h/2 = 0.15 m.",
        ),
        FrameCase(
            case_id="beam-ipe300-udl-bottom-flange",
            title="IPE300 beam 6 m, simply supported, 5 kN/m on the bottom flange",
            group="beam",
            nodes=_beam_nodes(),
            elements=_beam_elements(),
            line_loads=(LineLoad(1, 5.0, attack="Bottom"),),
            compare_in_plane=False,
            notes=_LTB_NOTE + " Stabilising load height z_g = -h/2 = -0.15 m.",
        ),
        FrameCase(
            case_id="beam-ipe300-midspan-point",
            title="IPE300 beam 6 m, simply supported, 20 kN at midspan",
            group="beam",
            nodes=_beam_nodes(),
            elements=_beam_elements(),
            point_loads=(PointLoad(1, SPAN_M / 2, 20.0),),
            compare_in_plane=False,
            notes=_LTB_NOTE + " M_Ed = 30 kNm, moment-distribution factor C1 about 1.35.",
        ),
        FrameCase(
            case_id="beam-he200a-udl",
            title="HE200A beam 6 m, simply supported, 10 kN/m at the shear centre",
            group="beam",
            nodes=_beam_nodes(),
            elements=_beam_elements("HE200A"),
            line_loads=(LineLoad(1, 10.0),),
            compare_in_plane=False,
            notes=_LTB_NOTE + " A wide-flange section with a large warping constant.",
        ),
        FrameCase(
            case_id="beam-ipe200-udl-4m",
            title="IPE200 beam 4 m, simply supported, 5 kN/m at the shear centre",
            group="beam",
            nodes=_beam_nodes(4.0),
            elements=_beam_elements("IPE200"),
            line_loads=(LineLoad(1, 5.0),),
            compare_in_plane=False,
            notes=_LTB_NOTE + " A short, light section; M_Ed = 10 kNm.",
        ),
        FrameCase(
            case_id="beam-ipe300-two-span",
            title="IPE300 beam 2 x 3 m over three supports, 10 kN/m",
            group="beam",
            nodes=(
                Node(1, 0.0, 0.0, "pinned"),
                Node(2, 3.0, 0.0, "roller-y"),
                Node(3, 6.0, 0.0, "roller-y"),
                Node(4, 0.0, STUB_M),
                Node(5, 6.0, STUB_M),
            ),
            elements=(Element(1, 1, 2), Element(2, 2, 3), Element(3, 1, 4), Element(4, 3, 5)),
            line_loads=(LineLoad(1, 10.0), LineLoad(2, 10.0)),
            compare_in_plane=False,
            notes=_LTB_NOTE + (
                " Hogging over the middle support. Both tools restrain the lateral "
                "displacement at every node, so each span is a 3 m segment."
            ),
        ),
    ]

    # -- portals -----------------------------------------------------------
    cases += [
        FrameCase(
            case_id="portal-ipe300-pinned-100-wind",
            title="IPE300 portal 6 x 4 m, pinned bases, 100 kN on each column + 5 kN wind",
            nodes=_portal_nodes("pinned"),
            elements=_portal_elements(),
            point_loads=_eaves_loads(100.0),
            notes=(
                "The frame of the earlier harness: PolyFrame printed alpha_cr = 17.2 "
                "(in plane) and alpha_cr,op = 5.63 for the right column."
            ),
        ),
        FrameCase(
            case_id="portal-ipe300-fixed-100-wind",
            title="IPE300 portal 6 x 4 m, fixed bases, 100 kN on each column + 5 kN wind",
            nodes=_portal_nodes("fixed"),
            elements=_portal_elements(),
            point_loads=_eaves_loads(100.0),
            known_deviations=(KnownDeviation("alpha_cr_op", 10.0, _FIXED_BASE_REASON),),
            notes="Fixed-base twin: PolyFrame printed 69.2 and 5.88.",
        ),
        FrameCase(
            case_id="portal-ipe300-pinned-300-wind",
            title="IPE300 portal 6 x 4 m, pinned bases, 300 kN on each column + 5 kN wind",
            nodes=_portal_nodes("pinned"),
            elements=_portal_elements(),
            point_loads=_eaves_loads(300.0),
            notes=(
                "alpha_cr = 5.72 in plane, below 10, so PolyFrame analyses the frame to "
                "second order and prints the factor inside an image (OCR). Out of "
                "plane 2.40."
            ),
        ),
        FrameCase(
            case_id="portal-he200a-pinned-100-wind",
            title="HE200A portal 6 x 4 m, pinned bases, 100 kN on each column + 5 kN wind",
            nodes=_portal_nodes("pinned"),
            elements=_portal_elements("HE200A"),
            point_loads=_eaves_loads(100.0),
            notes="PolyFrame printed 7.65 in plane and 10.42 out of plane.",
        ),
        FrameCase(
            case_id="portal-ipe300-pinned-100-symmetric",
            title="IPE300 portal 6 x 4 m, pinned bases, 100 kN on each column, no wind",
            nodes=_portal_nodes("pinned"),
            elements=_portal_elements(),
            point_loads=_eaves_loads(100.0, wind_kn=0.0),
            notes=(
                "No bending anywhere: the out-of-plane factor is the weak-axis Euler "
                "load of a column, 782 kN / 100 kN = 7.82, in both tools if they agree "
                "on the column's end conditions."
            ),
        ),
        FrameCase(
            case_id="portal-ipe300-pinned-udl-wind",
            title="IPE300 portal 6 x 4 m, pinned bases, 10 kN/m on the beam + 5 kN wind",
            nodes=_portal_nodes("pinned"),
            elements=_portal_elements(),
            line_loads=(LineLoad(2, 10.0),),
            point_loads=(PointLoad(1, HEIGHT_M, 5.0, "X", "W"),),
            notes="Beam in bending with light column compression: beam LTB can govern.",
        ),
        FrameCase(
            case_id="portal-ipe300-fixed-udl-wind",
            title="IPE300 portal 6 x 4 m, fixed bases, 10 kN/m on the beam + 5 kN wind",
            nodes=_portal_nodes("fixed"),
            elements=_portal_elements(),
            line_loads=(LineLoad(2, 10.0),),
            point_loads=(PointLoad(1, HEIGHT_M, 5.0, "X", "W"),),
            known_deviations=(KnownDeviation("alpha_cr_op", 10.0, _FIXED_BASE_REASON),),
        ),
        FrameCase(
            case_id="portal-ipe300-pinned-udl-top-flange",
            title="IPE300 portal 6 x 4 m, pinned bases, 10 kN/m on the beam's top flange + 5 kN wind",
            nodes=_portal_nodes("pinned"),
            elements=_portal_elements(),
            line_loads=(LineLoad(2, 10.0, attack="Top"),),
            point_loads=(PointLoad(1, HEIGHT_M, 5.0, "X", "W"),),
            known_deviations=(KnownDeviation("alpha_cr_op", 100.0, _TOP_FLANGE_FRAME_REASON),),
            notes="The load-height term inside a frame.",
        ),
        FrameCase(
            case_id="portal-he200a-fixed-udl-wind",
            title="HE200A portal 6 x 4 m, fixed bases, 10 kN/m on the beam + 5 kN wind",
            nodes=_portal_nodes("fixed"),
            elements=_portal_elements("HE200A"),
            line_loads=(LineLoad(2, 10.0),),
            point_loads=(PointLoad(1, HEIGHT_M, 5.0, "X", "W"),),
            known_deviations=(KnownDeviation("alpha_cr_op", 10.0, _FIXED_BASE_REASON),),
        ),
        FrameCase(
            case_id="gable-ipe300-pinned",
            title="IPE300 gable portal, 6 m span, 4 m eaves, 1 m rise, 20 kN at the apex + 5 kN wind",
            nodes=(
                Node(1, 0.0, 0.0, "pinned"),
                Node(2, 0.0, HEIGHT_M),
                Node(3, SPAN_M / 2, HEIGHT_M + 1.0),
                Node(4, SPAN_M, HEIGHT_M),
                Node(5, SPAN_M, 0.0, "pinned"),
            ),
            elements=(
                Element(1, 1, 2),
                Element(2, 2, 3),
                Element(3, 3, 4),
                Element(4, 5, 4),
            ),
            point_loads=(
                PointLoad(2, (9.0 + 1.0) ** 0.5, 20.0),
                PointLoad(1, HEIGHT_M, 5.0, "X", "W"),
            ),
            notes="Inclined rafters of 3.162 m carrying compression and bending.",
        ),
        FrameCase(
            case_id="two-bay-ipe300-pinned-udl-wind",
            title="IPE300 two-bay portal 2 x 6 m, pinned bases, 10 kN/m on both beams + 10 kN wind",
            nodes=(
                Node(1, 0.0, 0.0, "pinned"),
                Node(2, 0.0, HEIGHT_M),
                Node(3, SPAN_M, HEIGHT_M),
                Node(4, SPAN_M, 0.0, "pinned"),
                Node(5, 2 * SPAN_M, HEIGHT_M),
                Node(6, 2 * SPAN_M, 0.0, "pinned"),
            ),
            elements=(
                Element(1, 1, 2),
                Element(2, 2, 3),
                Element(3, 4, 3),
                Element(4, 3, 5),
                Element(5, 6, 5),
            ),
            line_loads=(LineLoad(2, 10.0), LineLoad(4, 10.0)),
            point_loads=(PointLoad(1, HEIGHT_M, 10.0, "X", "W"),),
            notes="The middle column carries twice the compression of the outer ones.",
        ),
    ]
    return cases


@lru_cache(maxsize=1)
def all_cases() -> tuple[FrameCase, ...]:
    cases = _catalogue()
    ids = [c.case_id for c in cases]
    duplicates = sorted({i for i in ids if ids.count(i) > 1})
    if duplicates:
        raise ValueError(f"duplicate case ids: {', '.join(duplicates)}")
    return tuple(cases)


def get_case(case_id: str) -> FrameCase:
    for case in all_cases():
        if case.case_id == case_id:
            return case
    raise KeyError(f"no case '{case_id}'")
