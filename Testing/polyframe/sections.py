"""Steel sections and material shared by the PolyFrame and PyNite renderings.

SI base units: m, m^2, m^4, m^6, Pa. ``Iy`` is the strong axis and ``Iz`` the weak
axis in Eurocode notation. PolyFrame resolves the same profiles by name from its
own library, so the PolyFrame model only carries the name; these numbers are the
PyNite side. They are the tabulated values of the European catalogue (ArcelorMittal
sales programme), the same ones PolyFrame prints on its report (IPE300:
``Iy = 83,561,092 mm^4``, ``A = 5,381 mm^2``).

Shear areas follow EN 1993-1-1 6.2.6(3): ``Avz`` is the web shear area for shear
in the plane of the web (strong-axis bending) and ``Avy`` is the flange area for
shear parallel to the flanges (weak-axis bending). PyNite's Timoshenko element takes
them as ``Asz`` and ``Asy``; PolyFrame includes shear deformation as well.
"""

from dataclasses import dataclass

STEEL_E = 210e9
STEEL_NU = 0.3
STEEL_G = STEEL_E / (2.0 * (1.0 + STEEL_NU))
STEEL_FY_S355 = 355e6


@dataclass(frozen=True)
class SteelSection:
    name: str
    h: float
    b: float
    tw: float
    tf: float
    r: float
    A: float
    Iy: float
    Iz: float
    It: float
    Iw: float

    @property
    def Avz(self) -> float:
        """Web shear area, EN 1993-1-1 6.2.6(3) for rolled I sections."""
        return self.A - 2.0 * self.b * self.tf + (self.tw + 2.0 * self.r) * self.tf

    @property
    def Avy(self) -> float:
        """Flange shear area, EN 1993-1-1 6.2.6(3)."""
        return 2.0 * self.b * self.tf


SECTIONS: dict[str, SteelSection] = {
    "IPE200": SteelSection(
        "IPE200", h=0.200, b=0.100, tw=0.0056, tf=0.0085, r=0.012,
        A=28.48e-4, Iy=1943e-8, Iz=142.4e-8, It=6.98e-8, Iw=12.99e-9,
    ),
    "IPE300": SteelSection(
        "IPE300", h=0.300, b=0.150, tw=0.0071, tf=0.0107, r=0.015,
        A=53.81e-4, Iy=8356e-8, Iz=603.8e-8, It=20.12e-8, Iw=125.9e-9,
    ),
    "HE200A": SteelSection(
        "HE200A", h=0.190, b=0.200, tw=0.0065, tf=0.0100, r=0.018,
        A=53.83e-4, Iy=3692e-8, Iz=1336e-8, It=20.98e-8, Iw=108.0e-9,
    ),
}
