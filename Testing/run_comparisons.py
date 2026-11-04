"""Run the reference-comparison scripts in this folder and collect their output.

These scripts are not pytest tests -- each one builds a model, runs an analysis and
prints a report meant to be compared against another program (PolyFrame, FEM Design,
...). This runner executes all of them in separate subprocesses, prints a status table
and writes a single markdown summary file containing every script's output so the whole
run can be handed to Claude (or anyone else) for analysis.

Usage (from the repository root)::

    python Testing/run_comparisons.py                    # run all, write the summary
    python Testing/run_comparisons.py --list             # list what would run
    python Testing/run_comparisons.py -k hinged          # only matching scripts
    python Testing/run_comparisons.py --quiet            # don't echo output to console
    python Testing/run_comparisons.py --timeout 120      # per-script time limit
    python Testing/run_comparisons.py --summary out.md   # write the summary elsewhere
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

TESTING_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTING_DIR.parent

DEFAULT_TIMEOUT = 600  # Seconds per script.
DEFAULT_SUMMARY = TESTING_DIR / 'comparison_summary.md'

# The comparison scripts, in the order they should be run. Each entry is
# (file name, one line description used in the summary).
COMPARISON_SCRIPTS: list[tuple[str, str]] = [
    ('braced_pitched_frame_comparison.py',
     'Pitched-roof portal frame with a diagonal brace (A -> C)'),
    ('buckling_comparison.py',
     'Two-bay single-storey portal frame: column and beam buckling scenarios'),
    ('hinged_frame_comparison.py',
     'Single-bay portal frame with a moment release at the top-left corner'),
    ('pitched_frame_comparison.py',
     'Single-bay pitched-roof (gable) portal frame'),
]


def resolve_scripts(pattern: str | None = None) -> list[tuple[Path, str]]:
    """Return the (path, description) pairs to run, skipping anything missing."""

    scripts: list[tuple[Path, str]] = []

    for name, description in COMPARISON_SCRIPTS:

        if pattern and pattern.lower() not in name.lower():
            continue

        path = TESTING_DIR / name

        if not path.exists():
            print(f'Warning: {name} not found in {TESTING_DIR} - skipping')
            continue

        scripts.append((path, description))

    return scripts


def child_environment() -> dict[str, str]:
    """Build the environment for a child process."""

    env = os.environ.copy()

    # Import `Pynite` from this repository rather than an installed copy.
    python_path = [str(REPO_ROOT)]
    if env.get('PYTHONPATH'):
        python_path.append(env['PYTHONPATH'])
    env['PYTHONPATH'] = os.pathsep.join(python_path)

    # Keep the captured output decodable regardless of the console code page, and keep
    # any stray matplotlib call from opening a window that blocks the run.
    env['PYTHONIOENCODING'] = 'utf-8'
    env.setdefault('MPLBACKEND', 'Agg')

    return env


def run_script(script: Path, timeout: int) -> tuple[str, float, str]:
    """Run one comparison script in a subprocess.

    :return: A ``(status, elapsed_seconds, output)`` tuple where status is one of
        ``'passed'``, ``'failed'`` or ``'timed out'``.
    """

    command = [sys.executable, str(script)]
    start = time.perf_counter()

    try:
        completed = subprocess.run(command, cwd=str(REPO_ROOT), env=child_environment(),
                                   capture_output=True, text=True, encoding='utf-8',
                                   errors='replace', timeout=timeout)
    except subprocess.TimeoutExpired as error:
        elapsed = time.perf_counter() - start
        output = (error.stdout or '') + (error.stderr or '')
        return 'timed out', elapsed, output

    elapsed = time.perf_counter() - start
    output = completed.stdout + completed.stderr
    status = 'passed' if completed.returncode == 0 else 'failed'

    return status, elapsed, output


def git_revision() -> str:
    """Return a short description of the checked out revision, or 'unknown'."""

    try:
        completed = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'],
                                   cwd=str(REPO_ROOT), capture_output=True, text=True,
                                   timeout=15)
        branch = subprocess.run(['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
                                cwd=str(REPO_ROOT), capture_output=True, text=True,
                                timeout=15)
    except Exception:
        return 'unknown'

    if completed.returncode != 0:
        return 'unknown'

    revision = completed.stdout.strip()
    if branch.returncode == 0 and branch.stdout.strip():
        revision += f' ({branch.stdout.strip()})'

    return revision


def indent(text: str, prefix: str = '    | ') -> str:
    """Indent a block of captured output so it reads as belonging to one script."""

    lines = text.rstrip().splitlines()
    return '\n'.join(prefix + line for line in lines)


def write_summary(path: Path, results: list[dict], total_elapsed: float) -> None:
    """Write a markdown summary of the run, output included, ready to hand off."""

    passed = sum(1 for result in results if result['status'] == 'passed')
    failed = [result for result in results if result['status'] != 'passed']

    lines: list[str] = []
    lines.append('# Pynite reference-comparison run')
    lines.append('')
    lines.append(f'- Run at: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    lines.append(f'- Repository: `{REPO_ROOT}`')
    lines.append(f'- Revision: {git_revision()}')
    lines.append(f'- Python: {sys.version.split()[0]}')
    lines.append(f'- Result: {passed} passed, {len(failed)} failed of {len(results)} '
                 f'script(s) in {total_elapsed:.1f}s')
    lines.append('')
    lines.append('These scripts print analysis results intended for comparison against '
                 'an external reference (PolyFrame, FEM Design, hand calculations). '
                 'The full console output of each run is reproduced below.')
    lines.append('')

    # Status table
    lines.append('## Status')
    lines.append('')
    lines.append('| Script | Status | Time | What it models |')
    lines.append('| --- | --- | --- | --- |')
    for result in results:
        lines.append(f'| `{result["name"]}` | {result["status"]} '
                     f'| {result["elapsed"]:.1f}s | {result["description"]} |')
    lines.append('')

    if failed:
        lines.append('Failed or timed out: '
                     + ', '.join(f'`{result["name"]}`' for result in failed))
        lines.append('')

    # Full output per script
    lines.append('## Output')
    lines.append('')
    for result in results:
        lines.append(f'### {result["name"]}')
        lines.append('')
        lines.append(f'{result["description"]}. Status: **{result["status"]}** '
                     f'({result["elapsed"]:.1f}s).')
        lines.append('')
        lines.append('```text')
        lines.append(result['output'].rstrip() or '(no output)')
        lines.append('```')
        lines.append('')

    path.write_text('\n'.join(lines), encoding='utf-8')


def main(argv: list[str] | None = None) -> int:

    parser = argparse.ArgumentParser(
        description='Run the reference-comparison scripts in the Testing folder.')
    parser.add_argument('-k', '--filter', dest='pattern', metavar='TEXT',
                        help='only run scripts whose file name contains TEXT '
                             '(case insensitive)')
    parser.add_argument('--list', action='store_true',
                        help='list the scripts that would run, then exit')
    parser.add_argument('--quiet', action='store_true',
                        help='only echo the output of scripts that fail (the summary '
                             'file always contains everything)')
    parser.add_argument('--timeout', type=int, default=DEFAULT_TIMEOUT, metavar='SECONDS',
                        help=f'time limit per script (default {DEFAULT_TIMEOUT})')
    parser.add_argument('--stop-on-error', action='store_true',
                        help='stop at the first script that fails')
    parser.add_argument('--summary', type=Path, default=DEFAULT_SUMMARY, metavar='PATH',
                        help=f'where to write the markdown summary '
                             f'(default {DEFAULT_SUMMARY.name} in this folder)')
    parser.add_argument('--no-summary', action='store_true',
                        help='skip writing the markdown summary file')

    args = parser.parse_args(argv)

    scripts = resolve_scripts(args.pattern)

    if not scripts:
        print('No comparison scripts to run'
              + (f' matching {args.pattern!r}' if args.pattern else ''))
        return 1

    if args.list:
        print(f'{len(scripts)} comparison script(s) in {TESTING_DIR}:')
        for path, description in scripts:
            print(f'  {path.name:<40} {description}')
        return 0

    print(f'Running {len(scripts)} comparison script(s) from {TESTING_DIR}\n')

    results: list[dict] = []
    total_start = time.perf_counter()

    for number, (path, description) in enumerate(scripts, start=1):

        print(f'[{number}/{len(scripts)}] {path.name} ... ', end='', flush=True)
        status, elapsed, output = run_script(path, args.timeout)
        print(f'{status} ({elapsed:.1f}s)')

        if output.strip() and (status != 'passed' or not args.quiet):
            print(indent(output))
            print()

        results.append({'name': path.name, 'description': description,
                        'status': status, 'elapsed': elapsed, 'output': output})

        if status != 'passed' and args.stop_on_error:
            print('Stopping at the first failure (--stop-on-error)')
            break

    total_elapsed = time.perf_counter() - total_start

    # Status table
    name_width = max(len(result['name']) for result in results)
    print('=' * (name_width + 22))
    print('Summary')
    print('=' * (name_width + 22))

    for result in results:
        print(f'{result["name"]:<{name_width}}  {result["status"]:<9}  '
              f'{result["elapsed"]:>6.1f}s')

    passed = sum(1 for result in results if result['status'] == 'passed')
    failed = [result['name'] for result in results if result['status'] != 'passed']

    print('=' * (name_width + 22))
    print(f'{passed} passed, {len(failed)} failed of {len(results)} run '
          f'in {total_elapsed:.1f}s')

    if failed:
        print('\nFailed scripts:')
        for name in failed:
            print(f'  {name}')

    if not args.no_summary:
        summary_path = args.summary if args.summary.is_absolute() \
            else (Path.cwd() / args.summary)
        write_summary(summary_path, results, total_elapsed)
        print(f'\nSummary written to {summary_path}')
        print('Hand that file to Claude for a separate analysis of the results.')

    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
