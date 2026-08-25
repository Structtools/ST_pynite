"""Run every other example script in this folder.

Each example is executed in its own subprocess so a crash, a `sys.exit()` call or a
leaked VTK/PyVista window in one example cannot affect the others. A summary table is
printed at the end and the exit code is non-zero if any example failed.

By default the examples run *headless*: the renderers are forced into off-screen mode
and matplotlib uses the non-interactive `Agg` backend, so nothing blocks waiting for a
window to be closed. The analysis and post-processing code still runs in full.

Usage (from anywhere)::

    python Examples/run_all_examples.py                  # headless, all examples
    python Examples/run_all_examples.py --list           # just list what would run
    python Examples/run_all_examples.py -k simple_beam   # only matching examples
    python Examples/run_all_examples.py --gui            # show the render windows
    python Examples/run_all_examples.py --quiet          # only show output of failures
    python Examples/run_all_examples.py --timeout 120    # per-example time limit

`--gui` will pause on every `render_model()` call until you close the window.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

EXAMPLES_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXAMPLES_DIR.parent
THIS_FILE = Path(__file__).resolve()

DEFAULT_TIMEOUT = 600  # Seconds per example. Some meshed models take a few minutes.

# Flag used when this file re-invokes itself to run a single example in a child process.
CHILD_FLAG = '--run-single'


# %% Child process: run one example
def _stub_renderer_module(module_name: str) -> None:
    """Register a do-nothing stand-in for a renderer module in `sys.modules`.

    This is only used as a fallback when the real module can't be imported (typically
    because `vtk` or `pyvista` isn't installed). It lets the rest of the example run so
    we still get analysis coverage instead of an import error.
    """

    import types

    class _StubRenderer:
        """Accepts any attribute assignment and swallows any method call."""

        def __init__(self, *args, **kwargs) -> None:
            pass

        def __getattr__(self, name: str):
            return lambda *args, **kwargs: None

        def __setattr__(self, name: str, value) -> None:
            pass

    module = types.ModuleType(module_name)
    module.Renderer = _StubRenderer
    module.render_model = lambda *args, **kwargs: None
    sys.modules[module_name] = module


def _make_renderers_headless() -> None:
    """Force Pynite's renderers to draw off-screen instead of opening a window."""

    # With the `Agg` backend every `plt.show()` in the examples warns that the canvas is
    # non-interactive. That's expected here, so keep it out of the captured output.
    import warnings

    warnings.filterwarnings('ignore', message='FigureCanvasAgg is non-interactive')

    # PyVista's global off-screen switch. `Pynite.Rendering.Renderer.__init__` reads this
    # when it builds its plotter, so it has to be set before any renderer is created.
    try:
        import pyvista as pv

        pv.OFF_SCREEN = True
    except Exception:
        pass

    # Pynite.Rendering (PyVista based). `render_model()` already supports off-screen
    # rendering, so we just force that argument on.
    try:
        from Pynite import Rendering
    except Exception:
        _stub_renderer_module('Pynite.Rendering')
    else:
        _original_render = Rendering.Renderer.render_model

        def _headless_render(self, reset_camera: bool = True, off_screen: bool = False):
            return _original_render(self, reset_camera=reset_camera, off_screen=True)

        Rendering.Renderer.render_model = _headless_render

        _original_screenshot = Rendering.Renderer.screenshot

        def _headless_screenshot(self, filepath: str = './Pynite_Image.png',
                                 interact: bool = True, reset_camera: bool = False):
            return _original_screenshot(self, filepath=filepath, interact=False,
                                        reset_camera=reset_camera)

        Rendering.Renderer.screenshot = _headless_screenshot

    # Pynite.Visualization (VTK based). Its `render_model()` blocks in the VTK interactor
    # unless `interact=False`, and the render window itself is switched off-screen.
    try:
        from Pynite import Visualization
    except Exception:
        _stub_renderer_module('Pynite.Visualization')
    else:
        _original_vtk_render = Visualization.Renderer.render_model

        def _headless_vtk_render(self, interact: bool = True, reset_camera: bool = True):
            try:
                self.window.SetOffScreenRendering(1)
            except Exception:
                pass
            return _original_vtk_render(self, interact=False, reset_camera=reset_camera)

        Visualization.Renderer.render_model = _headless_vtk_render

        _original_vtk_screenshot = Visualization.Renderer.screenshot

        def _headless_vtk_screenshot(self, filepath: str = 'console', interact: bool = True,
                                     reset_camera: bool = True):
            return _original_vtk_screenshot(self, filepath=filepath, interact=False,
                                            reset_camera=reset_camera)

        Visualization.Renderer.screenshot = _headless_vtk_screenshot


def run_single_example(script: Path, headless: bool) -> None:
    """Execute one example script in this process as if it were `__main__`."""

    if headless:
        _make_renderers_headless()

    source = script.read_text(encoding='utf-8')
    code = compile(source, str(script), 'exec')

    # Examples are written to be run from their own folder, and some of them read or
    # write files relative to the working directory.
    os.chdir(script.parent)

    globals_dict = {
        '__name__': '__main__',
        '__file__': str(script),
        '__builtins__': __builtins__,
    }
    exec(code, globals_dict)


# %% Parent process: discover and run every example
def discover_examples(pattern: str | None = None) -> list[Path]:
    """Return the example scripts to run, sorted by name, excluding this runner."""

    scripts = [path for path in sorted(EXAMPLES_DIR.glob('*.py'))
               if path.resolve() != THIS_FILE and not path.name.startswith('_')]

    if pattern:
        needle = pattern.lower()
        scripts = [path for path in scripts if needle in path.name.lower()]

    return scripts


def child_environment(headless: bool) -> dict[str, str]:
    """Build the environment for a child process."""

    env = os.environ.copy()

    # Make sure the example imports the `Pynite` package from this repository rather than
    # an installed copy.
    python_path = [str(REPO_ROOT)]
    if env.get('PYTHONPATH'):
        python_path.append(env['PYTHONPATH'])
    env['PYTHONPATH'] = os.pathsep.join(python_path)

    # Keep the captured output decodable regardless of the console code page.
    env['PYTHONIOENCODING'] = 'utf-8'

    if headless:
        # Matplotlib is imported lazily by `Member3D.plot_*()`, so setting the backend
        # through the environment is enough to keep `plt.show()` from blocking.
        env['MPLBACKEND'] = 'Agg'
        env['PYVISTA_OFF_SCREEN'] = 'true'

    return env


def run_example(script: Path, headless: bool, timeout: int) -> tuple[str, float, str]:
    """Run one example in a subprocess.

    :return: A ``(status, elapsed_seconds, output)`` tuple where status is one of
        ``'passed'``, ``'failed'`` or ``'timed out'``.
    """

    command = [sys.executable, str(THIS_FILE), CHILD_FLAG, str(script)]
    if not headless:
        command.append('--gui')

    start = time.perf_counter()

    try:
        completed = subprocess.run(command, cwd=str(EXAMPLES_DIR),
                                   env=child_environment(headless),
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


def indent(text: str, prefix: str = '    | ') -> str:
    """Indent a block of captured output so it reads as belonging to one example."""

    lines = text.rstrip().splitlines()
    return '\n'.join(prefix + line for line in lines)


def main(argv: list[str] | None = None) -> int:

    parser = argparse.ArgumentParser(
        description='Run every other example script in the Examples folder.')
    parser.add_argument('-k', '--filter', dest='pattern', metavar='TEXT',
                        help='only run examples whose file name contains TEXT '
                             '(case insensitive)')
    parser.add_argument('--list', action='store_true',
                        help='list the examples that would run, then exit')
    parser.add_argument('--gui', action='store_true',
                        help='run the examples unmodified, showing the render windows '
                             '(each one blocks until you close it)')
    parser.add_argument('--quiet', action='store_true',
                        help='only print the output of examples that fail')
    parser.add_argument('--timeout', type=int, default=DEFAULT_TIMEOUT, metavar='SECONDS',
                        help=f'time limit per example (default {DEFAULT_TIMEOUT})')
    parser.add_argument('--stop-on-error', action='store_true',
                        help='stop at the first example that fails')
    parser.add_argument(CHILD_FLAG, dest='single', metavar='SCRIPT',
                        help=argparse.SUPPRESS)  # Internal: run one example in-process

    args = parser.parse_args(argv)

    # Child mode: run the one script we were handed and let any exception propagate so
    # the parent sees a non-zero exit code and a traceback.
    if args.single:
        run_single_example(Path(args.single).resolve(), headless=not args.gui)
        return 0

    scripts = discover_examples(args.pattern)

    if not scripts:
        print('No examples found' + (f' matching {args.pattern!r}' if args.pattern else ''))
        return 1

    if args.list:
        print(f'{len(scripts)} example(s) in {EXAMPLES_DIR}:')
        for script in scripts:
            print(f'  {script.name}')
        return 0

    headless = not args.gui
    mode = 'headless' if headless else 'GUI'
    print(f'Running {len(scripts)} example(s) from {EXAMPLES_DIR} in {mode} mode\n')

    results: list[tuple[str, str, float]] = []
    total_start = time.perf_counter()

    for number, script in enumerate(scripts, start=1):

        print(f'[{number}/{len(scripts)}] {script.name} ... ', end='', flush=True)
        status, elapsed, output = run_example(script, headless, args.timeout)
        print(f'{status} ({elapsed:.1f}s)')

        if output.strip() and (status != 'passed' or not args.quiet):
            print(indent(output))

        results.append((script.name, status, elapsed))

        if status != 'passed' and args.stop_on_error:
            print('\nStopping at the first failure (--stop-on-error)')
            break

    total_elapsed = time.perf_counter() - total_start

    # Summary
    name_width = max(len(name) for name, _, _ in results)
    print('\n' + '=' * (name_width + 22))
    print('Summary')
    print('=' * (name_width + 22))

    for name, status, elapsed in results:
        print(f'{name:<{name_width}}  {status:<9}  {elapsed:>6.1f}s')

    passed = sum(1 for _, status, _ in results if status == 'passed')
    failed = [name for name, status, _ in results if status != 'passed']

    print('=' * (name_width + 22))
    print(f'{passed} passed, {len(failed)} failed of {len(results)} run '
          f'in {total_elapsed:.1f}s')

    if failed:
        print('\nFailed examples:')
        for name in failed:
            print(f'  {name}')

    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
