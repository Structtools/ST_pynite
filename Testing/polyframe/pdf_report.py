"""Read the critical load factors off PolyFrame's PDF report.

PolyFrame's API returns forces, reactions and aggregate utilisations and nothing
else. The two numbers this harness compares live only in the report:

* **in-plane** ``alpha_cr`` - on the text layer when it is above 10 ("The in-plane
  critical load factor 17.2 is larger than 10 ..."), otherwise inside a rendered
  formula image (``alpha_cr = 5.72``) because PolyFrame then runs a second-order
  analysis with an imperfection;
* **out-of-plane** ``alpha_cr,op`` - always inside images: the formula
  (``alpha_cr,op = 5.63``) in the stability check and the captions of the two
  buckling-shape plots in the appendix ("Critical load factor: alpha_cr = 5.63").

The images are read with RapidOCR after rendering the page region at 300 and
500 dpi with a white margin, which reads the LaTeX-style glyphs reliably; every
value is cross-checked between renderings and between the formula and the plot
captions, and a disagreement is recorded as a message rather than resolved
silently. The text layer also gives the governing out-of-plane element
("E3, ULS1") and the whole-percentage utilisations.

OCR is optional at import time: ``read_report`` works without it and simply
returns ``None`` for the image-only values, which the capture then flags.
"""

import io
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

_WATERMARK = "Student"

_VERSION = re.compile(r"^Version\s+(\d+\.\d+\.\d+)\b")
_ALPHA_CR_TEXT = re.compile(
    r"in-plane critical load factor\s+(\d+(?:\.\d+)?)\s+is larger than", re.IGNORECASE
)
_PERCENT = re.compile(r"^(\d+(?:\.\d+)?)\s*%\s*[≤<]\s*100\s*%")
_GOVERNING = re.compile(r"^E(\d+)\s*,")
_PLOT_FACTOR = re.compile(r"factor[^0-9]{0,12}(\d+(?:\.\d+)?)", re.IGNORECASE)

_OOP_HEADING = "Lateral and lateral torsional buckling"


@dataclass
class ReportValues:
    application_version: str | None = None
    alpha_cr_ip_text: float | None = None
    alpha_cr_ip_ocr: float | None = None
    alpha_cr_op_formula: float | None = None
    alpha_cr_op_plots: list[float] = field(default_factory=list)
    alpha_ult_k: list[float] = field(default_factory=list)
    governing_element_op: int | None = None
    util_stability: float | None = None
    util_cross_section: float | None = None
    messages: list[str] = field(default_factory=list)
    ocr_used: bool = False

    @property
    def alpha_cr_ip(self) -> float | None:
        return self.alpha_cr_ip_text if self.alpha_cr_ip_text is not None else self.alpha_cr_ip_ocr

    @property
    def alpha_cr_op(self) -> float | None:
        if self.alpha_cr_op_plots:
            return self.alpha_cr_op_plots[0]
        return self.alpha_cr_op_formula


class PdfToolingMissingError(RuntimeError):
    """PyMuPDF is not installed."""


def clean_lines(raw: list[str]) -> list[str]:
    return [s for line in raw if (s := line.strip()) and s not in _WATERMARK]


def _open(path: Path):
    try:
        import pymupdf  # noqa: PLC0415
    except ImportError as exc:
        raise PdfToolingMissingError(
            "PyMuPDF is needed to read PolyFrame's report: "
            "`uv sync --group polyframe` (or `pip install pymupdf`)."
        ) from exc
    return pymupdf, pymupdf.open(path)


# --- text layer --------------------------------------------------------------


def parse_text(lines: list[str], values: ReportValues) -> None:
    for line in lines:
        if values.application_version is None:
            m = _VERSION.match(line)
            if m:
                values.application_version = m.group(1)
        m = _ALPHA_CR_TEXT.search(line)
        if m:
            values.alpha_cr_ip_text = float(m.group(1))

    for i, line in enumerate(lines):
        if line.startswith(_OOP_HEADING):
            for candidate in lines[i + 1 : i + 4]:
                m = _GOVERNING.match(candidate)
                if m:
                    values.governing_element_op = int(m.group(1))
                    break

    for previous, line in zip(lines, lines[1:]):
        m = _PERCENT.match(line)
        if not m:
            continue
        value = float(m.group(1)) / 100.0
        if previous.lower().startswith("stability check"):
            values.util_stability = value
        elif previous.startswith("CS"):
            values.util_cross_section = value


# --- images ------------------------------------------------------------------


class _Ocr:
    """RapidOCR, loaded once and only when an image has to be read."""

    _engine = None

    @classmethod
    def available(cls) -> bool:
        try:
            import rapidocr_onnxruntime  # noqa: F401, PLC0415
            import PIL  # noqa: F401, PLC0415
        except ImportError:
            return False
        return True

    @classmethod
    def tokens(cls, png: bytes) -> list[str]:
        import numpy as np  # noqa: PLC0415
        from PIL import Image, ImageOps  # noqa: PLC0415
        from rapidocr_onnxruntime import RapidOCR  # noqa: PLC0415

        if cls._engine is None:
            cls._engine = RapidOCR()
        image = Image.open(io.BytesIO(png)).convert("RGB")
        image = ImageOps.expand(image, border=40, fill=(255, 255, 255))
        result, _ = cls._engine(np.array(image))
        return [row[1] for row in (result or [])]


def _normalise(tokens: list[str]) -> str:
    text = " ".join(tokens)
    # The Greek alpha renders as any of these in OCR output.
    for glyph in ("α", "Q", "&", "O", "a", "Ω", "o"):
        text = text.replace(f"{glyph}cr", "cr").replace(f"{glyph}ult", "ult")
    return text.replace(",", ".").replace(" .", ".").replace(". ", ".")


_DECIMAL = r"(\d+\.\d+)"
# ``label = number``, ``= number label`` or ``number label`` (OCR returns the tokens
# in either order and sometimes drops the equals sign).
_FORMULA = {
    "alpha_cr_op": re.compile(
        rf"(?:cr\.op\s*=\s*{_DECIMAL}|=?\s*{_DECIMAL}\s+\S*cr\.op)"
    ),
    "alpha_ult_k": re.compile(
        rf"(?:ult\.k\s*=\s*{_DECIMAL}|=?\s*{_DECIMAL}\s+\S*ult\.k)"
    ),
    "alpha_cr_ip": re.compile(
        rf"(?:(?<![a-z.])cr\s*=\s*{_DECIMAL}|=?\s*{_DECIMAL}\s+\S*cr(?![.,a-z]))"
    ),
}


def _read_formula(text: str) -> tuple[str, float] | None:
    """(kind, value) of the single ``label = number`` formula an image shows.

    The slenderness images (``lambda = sqrt(alpha_ult,k / alpha_cr)``) carry both
    labels and the reduction-factor image (``chi_op = min(chi, chi_LT)``) carries
    ``min``; both are skipped, their number is not a factor. A neighbouring image's
    label caught by the render margin cannot pass either, because it never sits
    next to the number.
    """
    if ("ult" in text and "cr" in text) or "min" in text:
        return None
    for kind, pattern in _FORMULA.items():
        m = pattern.search(text)
        if m:
            value = next(g for g in m.groups() if g)
            return kind, float(value)
    return None


def _vote(readings: list[float]) -> float | None:
    if not readings:
        return None
    value, count = Counter(readings).most_common(1)[0]
    return value


def _page_images(page):
    for info in page.get_image_info():
        w, h = info["width"], info["height"]
        if w >= 800 and h >= 800:  # the PolyFrame logo
            continue
        if w < 60:
            continue
        yield info


def parse_images(pymupdf, document, values: ReportValues) -> None:
    if not _Ocr.available():
        values.messages.append("OCR unavailable (rapidocr_onnxruntime / Pillow missing)")
        return
    values.ocr_used = True

    formula_votes: dict[str, list[float]] = {"alpha_cr_op": [], "alpha_ult_k": [], "alpha_cr_ip": []}
    plot_votes: list[float] = []

    for page in document:
        text = page.get_text()
        stability_page = "Stability Check" in text
        plots_page = "Stability Plots" in text
        if not (stability_page or plots_page):
            continue
        # Only the stability block of the results page carries the factors: stop
        # at the cross-section heading so the per-check formulas are not OCR'd.
        limit = None
        if stability_page:
            hits = page.search_for("Cross-Section CS")
            if hits:
                limit = min(h.y0 for h in hits)
        for info in _page_images(page):
            bbox = pymupdf.Rect(info["bbox"])
            if limit is not None and bbox.y0 >= limit:
                continue
            readings: dict[str, list[float]] = {}
            captions: list[float] = []
            for dpi in (300, 500):
                clip = pymupdf.Rect(bbox.x0 - 8, bbox.y0 - 8, bbox.x1 + 8, bbox.y1 + 8)
                pix = page.get_pixmap(clip=clip, dpi=dpi, alpha=False)
                tokens = _Ocr.tokens(pix.tobytes("png"))
                if not tokens:
                    continue
                joined = " ".join(tokens)
                if plots_page:
                    m = _PLOT_FACTOR.search(joined.replace(",", "."))
                    if m:
                        captions.append(float(m.group(1)))
                    continue
                formula = _read_formula(_normalise(tokens))
                if formula is not None:
                    kind, number = formula
                    readings.setdefault(kind, []).append(number)
            for kind, numbers in readings.items():
                voted = _vote(numbers)
                if voted is not None:
                    formula_votes[kind].append(voted)
                if len(set(numbers)) > 1:
                    values.messages.append(f"OCR renderings disagree for {kind}: {numbers}")
            if captions:
                plot_votes.append(_vote(captions))

    values.alpha_cr_op_formula = _vote(formula_votes["alpha_cr_op"])
    values.alpha_ult_k = formula_votes["alpha_ult_k"]
    values.alpha_cr_ip_ocr = _vote(formula_votes["alpha_cr_ip"])
    values.alpha_cr_op_plots = plot_votes

    if plot_votes and len(set(plot_votes)) > 1:
        values.messages.append(f"plot captions disagree on alpha_cr,op: {plot_votes}")
    if (
        values.alpha_cr_op_formula is not None
        and plot_votes
        and abs(values.alpha_cr_op_formula - plot_votes[0]) > 0.011
    ):
        values.messages.append(
            f"formula ({values.alpha_cr_op_formula}) and plot caption ({plot_votes[0]}) "
            "disagree on alpha_cr,op"
        )
    if (
        values.alpha_cr_ip_text is not None
        and values.alpha_cr_ip_ocr is not None
        and abs(values.alpha_cr_ip_text - values.alpha_cr_ip_ocr) > 0.011 * max(1.0, values.alpha_cr_ip_text)
    ):
        values.messages.append(
            f"text ({values.alpha_cr_ip_text}) and OCR ({values.alpha_cr_ip_ocr}) "
            "disagree on alpha_cr"
        )


# --- entry point -------------------------------------------------------------


def read_report(path: Path, ocr: bool = True) -> ReportValues:
    pymupdf, document = _open(Path(path))
    values = ReportValues()
    with document:
        lines: list[str] = []
        for page in document:
            lines.extend(page.get_text().splitlines())
        parse_text(clean_lines(lines), values)
        if ocr:
            parse_images(pymupdf, document, values)
    return values
