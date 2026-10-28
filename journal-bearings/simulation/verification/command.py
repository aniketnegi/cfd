"""Command-line reporting for verification postprocessing."""
import argparse
import contextlib
import io
from pathlib import Path
import sys
import time
import traceback


def require(condition, message):
    if not condition:
        raise ValueError(message)


def parser(description):
    result = argparse.ArgumentParser(description=description)
    result.add_argument('directory', type=Path, help='existing verification run directory')
    detail = result.add_mutually_exclusive_group()
    detail.add_argument('-v', '--verbose', action='store_true', help='show diagnostics and tracebacks')
    detail.add_argument('-q', '--quiet', action='store_true', help='show only errors')
    return result


def execute(action, args, *, label):
    started = time.monotonic()
    log = io.StringIO()
    try:
        with contextlib.redirect_stdout(sys.stderr if args.verbose else log):
            result = action()
    except KeyboardInterrupt:
        print(f'{label}: interrupted', file=sys.stderr)
        return 130
    except Exception as error:
        if args.verbose:
            traceback.print_exc()
        else:
            print(f'{label}: {error}', file=sys.stderr)
        return 1
    if not args.quiet:
        print(f'ok {label} [{time.monotonic() - started:.1f}s]', file=sys.stderr)
        if isinstance(result, Path):
            print(f'output {result.resolve()}', file=sys.stderr)
    return 0
