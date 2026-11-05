"""PolyFrame buckling comparison harness.

Compares the in-plane and out-of-plane critical load factors of PyNite's
``FEModel3D.analyze_buckling`` with PolyFrame's general method (EN 1993-1-1
6.3.4), which reports an in-plane ``alpha_cr`` and an out-of-plane
``alpha_cr,op`` per frame.

The package has two halves that never run in the same process:

* capture (``python -m Testing.polyframe refresh``) drives PolyFrame through the
  C# bridge in ``bridge/`` on a licensed Windows machine, reads the PDF report
  and writes one JSON record per case to ``corpus/``;
* comparison (``pytest Testing/test_polyframe_buckling.py``) reads the corpus and
  runs PyNite. It never starts PolyFrame.

See ``README.md`` in this directory.
"""
