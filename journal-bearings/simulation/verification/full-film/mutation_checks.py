"""Check that deliberate shared-core physics defects fail the conical tests.

An unmodified control and each mutation run in an isolated source tree. Logs
use a fresh output directory; production files and saved cases are untouched.
This is a narrow sensitivity check, not proof against every defect.
"""
import argparse
import sys
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from command import execute

ROOT = Path(__file__).resolve().parents[3]
SOLVER = Path("simulation/reynolds/src/full_film.jl")
TESTS = Path("simulation/tests/conical_checks.jl")
SOURCES = (
    SOLVER, TESTS,
    Path("simulation/reynolds/src/conical.jl"),
    Path("simulation/reynolds/src/cylindrical.jl"),
    Path("simulation/reynolds/src/pressure_vtk.jl"),
    Path("bearing.jl"), Path("mesh/2d/grid.jl"),
    Path("mesh/2d/quality.jl"), Path("mesh/2d/metrics.jl"), Path("mesh/vtk.jl"),
)
MUTATIONS = {
    "missing_cone_distance_factor": (
        "width = axes.z_widths[j] / ca", "width = axes.z_widths[j]"),
    "reversed_wedge_source": (
        "b[owner] -= adv[i,j]; b[neighbour] += adv[i,j]",
        "b[owner] += adv[i,j]; b[neighbour] -= adv[i,j]"),
}


def run_mutations(output):
    if output is None:
        parent = ROOT / "out/reviews"
        parent.mkdir(parents=True, exist_ok=True)
        output = Path(tempfile.mkdtemp(prefix="conical-mutations.", dir=parent))
    else:
        output.mkdir(parents=True, exist_ok=False)
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       JULIA_NUM_THREADS="1")
    results = []
    with tempfile.TemporaryDirectory() as directory:
        sandbox = Path(directory)
        for path in SOURCES:
            target = sandbox / path
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / path, target)
        source = (sandbox / SOLVER).read_text()

        def run_tests(name):
            run = subprocess.run(
                ["julia", f"--project={ROOT}", "--threads=1", "--gcthreads=1",
                 "--check-bounds=yes", str(sandbox / TESTS)],
                cwd=sandbox, env=environment, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, timeout=120)
            (output / f"{name}.log").write_text(run.stdout)
            return run

        control = run_tests("control")
        if control.returncode != 0:
            raise RuntimeError(f"Unmodified control failed; inspect {output / 'control.log'}")
        for name, (before, after) in MUTATIONS.items():
            if source.count(before) != 1:
                raise RuntimeError(f"Mutation {name} must match exactly once in {SOLVER}")
            (sandbox / SOLVER).write_text(source.replace(before, after))
            run = run_tests(f"mutation-{name}")
            detected = run.returncode != 0 and "Test Failed" in run.stdout
            results.append(dict(mutation=name, exit_code=run.returncode, detected=detected))
            print(f"{name}: detected={detected}", flush=True)
            if not detected:
                raise RuntimeError(f"Mutation {name} was not detected; inspect {output}")
    hashes = {str(path): hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
              for path in (*SOURCES, Path(__file__).resolve().relative_to(ROOT),
                           Path("simulation/verification/command.py"))}
    (output / "mutations.json").write_text(json.dumps({
        "source_sha256": hashes, "control_exit_code": control.returncode,
        "results": results}, indent=2) + "\n")
    print(f"Saved mutation checks to {output}")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", type=Path,
                        help="fresh output directory (default: unique directory under out/reviews)")
    detail = parser.add_mutually_exclusive_group()
    detail.add_argument("-v", "--verbose", action="store_true", help="show diagnostics and tracebacks")
    detail.add_argument("-q", "--quiet", action="store_true", help="show only errors")
    args = parser.parse_args()
    if args.output is not None and (args.output.exists() or args.output.is_symlink()):
        parser.error(f"output already exists: {args.output}")
    return execute(lambda: run_mutations(args.output), args, label="mutation checks")


if __name__ == "__main__":
    raise SystemExit(main())
