#!/usr/bin/env python3
"""Resource-aware queue of isolated single-threaded Julia mesh studies."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time
import tomllib

SOURCE_FILES = [
    "bearing.jl", "mesh/2d/conical/grid.jl", "mesh/2d/conical/quality.jl",
    "mesh/2d/conical/metrics.jl", "mesh/2d/conical/surface_feed.jl",
    "mesh/2d/conical/inlet/controlled.jl",
    "mesh/2d/conical/inlet/seed.jl", "mesh/2d/conical/inlet/placement.jl",
    "mesh/examples/compare_inlet_iterations.jl", "mesh/examples/inlet_long_run.jl",
]
CASES = ["reference-fine", "reference", "diameter-8", "combined-odd",
         "held-out", "reference-medium", "squarer-fine"]
GIB = 1024**3


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    temporary.replace(path)


def available_memory():
    info = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
    available = int(info["MemAvailable"].split()[0]) * 1024
    limit, current = Path("/sys/fs/cgroup/memory.max"), Path("/sys/fs/cgroup/memory.current")
    if limit.exists() and current.exists() and limit.read_text().strip() != "max":
        available = min(available, max(0, int(limit.read_text()) - int(current.read_text())))
    return available


def resident_memory(pid):
    try:
        return int(next(line.split()[1] for line in Path(f"/proc/{pid}/status").read_text().splitlines()
                        if line.startswith("VmRSS:"))) * 1024
    except (OSError, StopIteration):
        return 0


def read_toml(path):
    try:
        return tomllib.loads(path.read_text())
    except (OSError, tomllib.TOMLDecodeError):
        return {}


def terminate_owned(process):
    if process.poll() is None:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass


def read_progress(path):
    latest = read_toml(path)
    keys = ["phase", "iteration", "accepted_steps", "energy", "scaled_gradient_norm",
            "snapshot_iteration", "quality_iteration", "quality_pass", "termination_reason"]
    return {key: latest[key] for key in keys if key in latest}


def run_queue(args):
    root = Path(__file__).resolve().parents[2]
    directory = args.directory.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    archive = directory / "source"
    hashes = {}
    files = SOURCE_FILES + ["mesh/examples/queue_inlet_long_runs.py"]
    files += [name for name in ["Project.toml", "Manifest.toml"] if (root / name).exists()]
    for name in files:
        source, target = root / name, archive / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        hashes[name] = hashlib.sha256(target.read_bytes()).hexdigest()
        assert hashes[name] == hashlib.sha256(source.read_bytes()).hexdigest()
    estimates = {"reference-medium": 1.5, "reference-fine": args.job_gib,
                 "squarer-fine": args.job_gib}
    job_reservations = {name: min(args.job_gib, estimates.get(name, 1.)) for name in args.cases}
    atomic_json(directory / "inputs.json", {
        "created_utc": utc_now(), "root": str(root), "max_workers": args.workers,
        "max_iterations": args.iterations, "reserve_gib": args.reserve_gib,
        "job_reservation_gib": args.job_gib, "job_reservations_gib": job_reservations, "debug": args.debug,
        "source_sha256": hashes, "cases": args.cases,
        "note": "Frozen sources; resource controls affect scheduling, never optimizer stopping",
    })
    (directory / "queue.pid").write_text(str(os.getpid()) + "\n")
    jobs = [{"case": name, "state": "queued", "pid": None} for name in args.cases]
    active = {}
    stopping = False

    def request_stop(signum, frame):
        nonlocal stopping
        stopping = True

    signal.signal(signal.SIGTERM, request_stop)
    signal.signal(signal.SIGINT, request_stop)
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               JULIA_NUM_THREADS="1", DEBUG="1" if args.debug else "0")
    resources = (directory / "resources.tsv").open("w", buffering=1)
    resources.write("utc\tload1\tavailable_gib\trunning\tqueued\n")
    try:
        while True:
            for case, (process, stdout, stderr) in list(active.items()):
                code = process.poll()
                if code is not None:
                    stdout.close(); stderr.close()
                    job = next(job for job in jobs if job["case"] == case)
                    result = read_toml(directory / case / "result.toml")
                    state = "cancelled" if stopping else "finished" if code == 0 and result else "failed"
                    job.update(state=state, exit_code=code, finished_utc=utc_now(), result=result,
                               progress=read_progress(directory / case / "latest.toml"))
                    del active[case]
            if stopping:
                for process, _, _ in active.values():
                    terminate_owned(process)
                for job in jobs:
                    if job["state"] == "queued":
                        job["state"] = "cancelled"
            load = os.getloadavg()[0]
            available = available_memory()
            pending_reservation = sum(max(0, job_reservations[case] * GIB - resident_memory(p.pid))
                                      for case, (p, _, _) in active.items())
            next_job = next((job for job in jobs if job["state"] == "queued"), None)
            next_reservation = job_reservations[next_job["case"]] if next_job else args.job_gib
            cpu_limit = max(1, len(os.sched_getaffinity(0)) - 2)
            can_start = (not stopping and not (directory / "STOP").exists()
                         and len(active) < args.workers and load < cpu_limit
                         and available - pending_reservation >= (args.reserve_gib + next_reservation) * GIB)
            if can_start and next_job:
                case = next_job["case"]
                stdout = (directory / f"{case}.stdout.log").open("w")
                stderr = (directory / f"{case}.stderr.log").open("w")
                command = ["nice", "-n", "10"]
                if shutil.which("ionice"):
                    command += ["ionice", "-c", "3"]
                command += ["julia", f"--project={archive}", "--threads=1", "--gcthreads=1",
                            "--heap-size-hint=2G", "--check-bounds=yes", "--depwarn=error",
                            str(archive / "mesh/examples/inlet_long_run.jl"), case,
                            str(directory / case), str(args.iterations)]
                if args.debug:
                    command += ["--debug"]
                process = subprocess.Popen(command, cwd=archive, env=env, stdout=stdout,
                                           stderr=stderr, start_new_session=True)
                active[case] = process, stdout, stderr
                next_job.update(state="running", pid=process.pid, started_utc=utc_now(), command=command)
            for job in jobs:
                if job["state"] == "running":
                    job["progress"] = read_progress(directory / job["case"] / "latest.toml")
                    job["rss_gib"] = resident_memory(job["pid"]) / GIB
            remaining = sum(job["state"] == "queued" for job in jobs)
            status = {"updated_utc": utc_now(), "queue_pid": os.getpid(), "running": len(active),
                      "queued": remaining, "paused": (directory / "STOP").exists(),
                      "available_gib": available / GIB, "load1": load, "jobs": jobs}
            atomic_json(directory / "status.json", status)
            resources.write(f"{utc_now()}\t{load}\t{available/GIB}\t{len(active)}\t{remaining}\n")
            report = ["# Live L-BFGS long-run study", "", f"Updated: {status['updated_utc']}", "",
                      "Historical stability_steps=0 mode; process completion is not mesh acceptance.", "",
                      "|Case|Queue state|Iteration|Stop reason|Quality (last sampled)|", "|---|---|---:|---|---|"]
            for job in jobs:
                progress = job.get("progress", {})
                result = job.get("result", {})
                if job["state"] in ("failed", "cancelled"):
                    result = {"termination_reason": f"process_exit_{job.get('exit_code', 'cancelled')}"}
                report.append(f"|{job['case']}|{job['state']}|{progress.get('iteration','—')}|"
                              f"{result.get('termination_reason',progress.get('termination_reason','—'))}|"
                              f"{result.get('quality_pass',progress.get('quality_pass','—'))}|")
            report += ["", "Read status.json, CASE/latest.toml, CASE/progress.tsv and CASE/quality.tsv while running.",
                       "Snapshots: CASE/latest-state.jls + CASE/seed.jls (original anchors); not a terminal-mesh restart API.",
                       "Touch STOP to pause new launches; remove it to resume. SIGTERM queue PID cancels only this study's process groups."]
            temporary = directory / "report.md.tmp"
            temporary.write_text("\n".join(report) + "\n")
            temporary.replace(directory / "report.md")
            if not active and not remaining:
                break
            time.sleep(5)
    finally:
        for process, stdout, stderr in active.values():
            terminate_owned(process)
            stdout.close()
            stderr.close()
        resources.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="fresh study directory")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--iterations", type=int, default=20000)
    parser.add_argument("--reserve-gib", type=float, default=3.)
    parser.add_argument("--job-gib", type=float, default=2.25)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--cases", nargs="+", choices=CASES, default=CASES)
    args = parser.parse_args()
    if len(set(args.cases)) != len(args.cases):
        parser.error("case names must be unique")
    if args.workers < 1 or args.iterations < 1 or args.reserve_gib < 0 or args.job_gib <= 0:
        parser.error("invalid resource/budget setting")
    run_queue(args)
