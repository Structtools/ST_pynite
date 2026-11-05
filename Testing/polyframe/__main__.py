"""CLI for the PolyFrame buckling corpus.

    uv run --group polyframe python -m Testing.polyframe list
    uv run --group polyframe python -m Testing.polyframe refresh [--case <id> ...]
    uv run --group polyframe python -m Testing.polyframe reparse [--case <id> ...]
    uv run --group polyframe python -m Testing.polyframe pynite  [--case <id> ...]
    uv run --group polyframe python -m Testing.polyframe summary
    uv run --group polyframe python -m Testing.polyframe reference [--case <id> ...]

``refresh`` is the only command that starts PolyFrame (Windows, licence signed in
once by hand). ``reparse`` re-reads the saved reports with the current parser.
``pynite`` prints PyNite's factors. ``summary`` writes
``docs/POLYFRAME_BUCKLING_COMPARISON.md`` from the corpus and PyNite.
"""

import argparse
import sys
import tempfile
from pathlib import Path

from . import pynite_side
from .capture import CaptureError, record_from_outcome, record_from_report
from .cases import all_cases, get_case
from .corpus import REPORTS_DIR, read_record, write_record
from .runner import CaseJob, PolyFrameUnavailableError, check_available, run_bridge
from .summary import write_summary
from .translate import model_hash, to_polyframe_model


def _selected(args: argparse.Namespace):
    if args.case:
        return [get_case(c) for c in args.case]
    return list(all_cases())


def _cmd_list(args: argparse.Namespace) -> int:
    for case in all_cases():
        try:
            record = read_record(case.case_id)
            stale = record["provenance"]["model_sha256"] != model_hash(to_polyframe_model(case))
            state = "STALE" if stale else f"captured {record['provenance']['captured_on']}"
            pf = record["polyframe"]
            state += f"  ip={pf.get('alpha_cr_ip')} op={pf.get('alpha_cr_op')}"
        except FileNotFoundError:
            state = "not captured"
        print(f"{case.case_id:<40} {state}")
    return 0


def _cmd_refresh(args: argparse.Namespace) -> int:
    cases = _selected(args)
    try:
        check_available()
    except PolyFrameUnavailableError as exc:
        print(f"cannot capture: {exc}", file=sys.stderr)
        return 2
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    jobs = [
        CaseJob(case.case_id, to_polyframe_model(case), REPORTS_DIR / f"{case.case_id}.pdf")
        for case in cases
    ]
    work_dir = Path(tempfile.mkdtemp(prefix="polyframe-"))
    print(f"running {len(jobs)} case(s) through PolyFrame (job files in {work_dir})")
    outcome = run_bridge(jobs, work_dir)
    failures = 0
    for case in cases:
        try:
            record = record_from_outcome(
                case, outcome.outcome_for(case.case_id), connector_version=outcome.connector_version
            )
        except CaptureError as exc:
            print(f"FAILED {exc}")
            failures += 1
            continue
        path = write_record(record)
        pf = record["polyframe"]
        flag = "  !! " + "; ".join(record["messages"]) if record["messages"] else ""
        print(f"{case.case_id:<40} ip={pf['alpha_cr_ip']} op={pf['alpha_cr_op']} -> {path.name}{flag}")
    return 1 if failures else 0


def _cmd_reparse(args: argparse.Namespace) -> int:
    failures = 0
    for case in _selected(args):
        pdf = REPORTS_DIR / f"{case.case_id}.pdf"
        if not pdf.exists():
            print(f"{case.case_id}: no report")
            failures += 1
            continue
        try:
            old = read_record(case.case_id)
            connector = old["provenance"]["connector_version"]
            captured_on = old["provenance"]["captured_on"]
            result = None
        except FileNotFoundError:
            connector, captured_on, result = "unknown", None, None
        record = record_from_report(
            case, pdf, result=result, connector_version=connector, captured_on=captured_on
        )
        if result is None:
            try:
                for key in ("axial_kn", "moment_knm", "uls_utilization"):
                    record["polyframe"][key] = old["polyframe"].get(key)
            except NameError:
                pass
        write_record(record)
        pf = record["polyframe"]
        flag = "  !! " + "; ".join(record["messages"]) if record["messages"] else ""
        print(f"{case.case_id:<40} ip={pf['alpha_cr_ip']} op={pf['alpha_cr_op']}{flag}")
    return 1 if failures else 0


def _cmd_pynite(args: argparse.Namespace) -> int:
    for case in _selected(args):
        r = pynite_side.run(case)
        ip = f"{r.alpha_cr_ip:.3f}" if r.alpha_cr_ip is not None else "-"
        print(
            f"{case.case_id:<40} ip={ip:>8}  op(ignore)={r.alpha_cr_op_ignore:.3f} [E{r.op_element_ignore}]"
            f"  op(eq.torsion, L)={r.alpha_cr_op_eqtorsion:.3f} [E{r.op_element_eqtorsion}]"
            f"  op(eq.torsion, restraints)={r.alpha_cr_op_restraints:.3f} [E{r.op_element_restraints}]"
        )
    return 0


def _cmd_summary(args: argparse.Namespace) -> int:
    path = write_summary()
    print(f"wrote {path}")
    return 0


def _cmd_reference(args: argparse.Namespace) -> int:
    from .reference_checks import DEFAULT_CASES, format_table, run_all

    ids = tuple(args.case) if args.case else DEFAULT_CASES
    print(format_table(run_all(ids)))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m Testing.polyframe")
    sub = parser.add_subparsers(dest="command", required=True)
    for name, fn in (
        ("list", _cmd_list),
        ("refresh", _cmd_refresh),
        ("reparse", _cmd_reparse),
        ("pynite", _cmd_pynite),
        ("summary", _cmd_summary),
        ("reference", _cmd_reference),
    ):
        p = sub.add_parser(name)
        p.add_argument("--case", action="append", help="case id (repeatable)")
        p.set_defaults(fn=fn)
    args = parser.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
