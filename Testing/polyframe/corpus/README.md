# PolyFrame corpus

One JSON record per catalogue case, captured from PolyFrame 1.1.11 (connector 1.1.8.0)
on 2026-10-02 with `python -m Testing.polyframe refresh`; the PDF report each record
was read from goes in `reports/`, which is gitignored (student licence, below) and kept
only on the capturing machine. Records carry the SHA-256 of the PolyFrame model they
were captured from; `Testing/test_polyframe_buckling.py` fails when the catalogue entry
no longer hashes to it.

**Licence.** These reports were produced under a PolyStruc **student licence** (every page
carries a "Student License" footer). PolyStruc's terms restrict student licences to
educational use. The numbers are what they are, but before any out-of-plane value from
this corpus is cited in project documentation, and before this repository is made public,
recapture under a commercial licence (`refresh`, about 25 minutes) and re-run the gate;
the hash check makes any changed number visible. No further student-licence captures.

**Where the numbers come from.** `alpha_cr_ip` is PolyFrame's text layer when above 10 and
OCR of its formula image otherwise; `alpha_cr_op` is OCR of the formula image,
cross-checked against the two plot captions in the appendix (`alpha_cr_op_formula_ocr`,
`alpha_cr_op_plot_ocr`). `*_source` says which. A value typed in by hand from the report
goes in `hand_read.json` and is marked `hand`.

**What is gated.** The in-plane and out-of-plane factors, both at ±5 %. Known deviations
are one-sided: they excuse PyNite being *above* PolyFrame up to a stated bound, never
below. See `../README.md`.
