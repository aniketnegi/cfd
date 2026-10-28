"""Exercise verification CLI contracts without running refinement sweeps.

Run with: python simulation/verification/interface_checks.py
"""
import json
import os
import select
import time
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
VERIFY = ROOT / 'simulation/verification'
RUNNERS = [
    'full-film/run.jl', 'full-film/field_checks.jl',
    'full-film/reference_profiles.jl', 'full-film/geometry_review.jl',
    'full-film/conical_verification.jl',
    'swift-stieber/run_cylindrical.jl', 'swift-stieber/run_conical.jl',
    'swift-stieber/compare_fd.jl', 'swift-stieber/fd_checks.jl',
    'swift-stieber/fv_flux_checks.jl', 'swift-stieber/long_bearing.jl',
    'swift-stieber/conical_mms.jl', 'swift-stieber/conical_checks.jl',
    'jfo/slider/jfo_slider_verification.jl',
    'jfo/conical/conical_jfo_verification.jl',
]
ENV = dict(os.environ, JULIA_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           PYTHONDONTWRITEBYTECODE='1')


class Interfaces(unittest.TestCase):
    def command(self, command, cwd):
        return subprocess.run(command, cwd=cwd, env=ENV, capture_output=True,
                              text=True, timeout=90)

    def test_julia_entry_points(self):
        with tempfile.TemporaryDirectory() as directory:
            for runner in RUNNERS:
                path = VERIFY / runner
                command = ['julia', f'--project={ROOT}', str(path)]
                with self.subTest(runner=runner):
                    help_run = self.command(command + ['--help'], directory)
                    self.assertEqual(help_run.returncode, 0, help_run.stderr)
                    self.assertIn('Usage:', help_run.stdout)
                    self.assertEqual(help_run.stderr, '')
                    bad = self.command(command + ['--unknown'], directory)
                    self.assertEqual(bad.returncode, 2, bad.stderr)
                    self.assertEqual(bad.stdout, '')
                    self.assertNotIn('Stacktrace', bad.stderr)
                    included = self.command(['julia', f'--project={ROOT}', '-e',
                                             f'include({json.dumps(str(path))})'], directory)
                    self.assertEqual(included.returncode, 0, included.stderr)
                    self.assertEqual(included.stdout, '')
                    self.assertEqual(list(Path(directory).iterdir()), [])

    def test_shared_options_and_failure_reporting(self):
        code = '''
include(joinpath(ARGS[1], "cli.jl"))
using Test
mktempdir() do dir
    path = joinpath(dir, "fresh")
    @test BearingCLI.study_options(["--full", path]; full=true).full
    @test BearingCLI.study_options([path, "--full"]; full=true).full
    @test_throws ArgumentError BearingCLI.study_options([path, "extra"])
    @test_throws ArgumentError BearingCLI.study_options([path, "-q", "-v"])
    @test_throws ArgumentError BearingCLI.study_options([dir])
    @test_throws ArgumentError BearingCLI.study_options([path, "--full"])
    @test_throws ArgumentError BearingCLI.study_options([path, "--full", "--full"]; full=true)
    symlink(joinpath(dir, "absent"), path)
    @test_throws ArgumentError BearingCLI.study_options([path])
    @test BearingCLI.report(() -> error("fixture failure"); quiet=true) == 1
    @test BearingCLI.report(() -> throw(InterruptException()); quiet=true) == 130
end
'''
        result = self.command(['julia', f'--project={ROOT}', '-e', code, str(ROOT)], ROOT)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('fixture failure', result.stderr)
        self.assertNotIn('Stacktrace', result.stderr)

    def test_terminal_progress(self):
        master, slave = os.openpty()
        code = '''
include(joinpath(ARGS[1], "cli.jl"))
exit(BearingCLI.report(; label="busy fixture") do
    start = time_ns()
    while (time_ns() - start) < 3e9
        yield()  # Fitted placement yields at its owned observer snapshots.
    end
end)
'''
        process = subprocess.Popen(['julia', f'--project={ROOT}', '-e', code, str(ROOT)],
                                   stdout=subprocess.PIPE, stderr=slave, env=ENV)
        os.close(slave)
        transcript = b''
        live_dot = False
        deadline = time.monotonic() + 30
        try:
            while time.monotonic() < deadline:
                ready, _, _ = select.select([master], [], [], 0.1)
                if ready:
                    try:
                        chunk = os.read(master, 4096)
                    except OSError:
                        break
                    transcript += chunk
                    if b'.' in transcript and b'[OK]' not in transcript:
                        live_dot = True
                if process.poll() is not None and not ready:
                    break
            self.assertEqual(process.wait(timeout=5), 0, transcript)
            self.assertTrue(live_dot, transcript)
            self.assertIn(b'[OK]', transcript)
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
            process.stdout.close()
            os.close(master)

    def test_mesh_streams_and_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            command = ['julia', f'--project={ROOT}', str(ROOT / 'mesh/2d/generate.jl'),
                       '--geometry', 'unported', '--journal-radius-mm', '50',
                       '--length-mm', '100', '--clearance-mm', '.05',
                       '--cone-angle-deg', '10', '--eccentricity-ratio', '.6',
                       '--n-theta', '16', '--n-axial', '4']
            for mode in ['', '-q', '-v']:
                path = run / f'mesh{mode}.vtk'
                result = self.command(command + ['--output', str(path)] + ([mode] if mode else []), run)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout, '')
                if mode == '-q':
                    self.assertEqual(result.stderr, '')
                elif mode == '-v':
                    self.assertIn('Surface quality', result.stderr)
                else:
                    self.assertIn('output ', result.stderr)
                    self.assertNotIn('Surface quality', result.stderr)
                self.assertEqual(path.read_bytes(), (run / 'mesh.vtk').read_bytes())
            blocked = run / 'blocked'
            blocked.write_text('not a directory')
            failed = self.command(command + ['--output', str(blocked / 'mesh.vtk'), '-q'], run)
            self.assertEqual(failed.returncode, 1, failed.stderr)
            self.assertEqual(failed.stdout, '')
            self.assertTrue(failed.stderr)
            self.assertNotIn('Stacktrace', failed.stderr)
            existing = self.command(command + ['--output', str(run / 'mesh.vtk')], run)
            self.assertEqual(existing.returncode, 2, existing.stderr)

    def test_python_validation_survives_optimization(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            (run / 'inputs.toml').write_text('[source_sha256]\n')
            for suite in ['full-film', 'swift-stieber']:
                path = VERIFY / suite / 'analyse.py'
                # An empty/unsupported run must fail even with assertions disabled.
                if suite == 'full-film':
                    (run / 'results.tsv').write_text('model\tL_D\tepsilon\n')
                else:
                    (run / 'results.tsv').unlink()
                for optimize in [[], ['-O']]:
                    result = self.command([sys.executable, *optimize, str(path), '-q', directory], ROOT)
                    self.assertEqual(result.returncode, 1, result.stderr)
                    self.assertEqual(result.stdout, '')
                    self.assertNotIn('Traceback', result.stderr)
                for name in ['summary.json', 'analysis.json']:
                    self.assertFalse((run / name).exists())
            # A real failed tolerance check, not just a missing file.
            (run / 'fd-mms.tsv').write_text('positive_pressure_residual\n1.0\n')
            result = self.command([sys.executable, '-O', str(VERIFY / 'swift-stieber/analyse.py'),
                                   directory], ROOT)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertIn('exceeds tolerance', result.stderr)

    def test_python_help_errors_and_plot_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            for suite in ['full-film', 'swift-stieber']:
                for name in ['analyse.py', 'plot.py']:
                    command = [sys.executable, str(VERIFY / suite / name)]
                    help_run = self.command(command + ['--help'], ROOT)
                    self.assertEqual(help_run.returncode, 0, help_run.stderr)
                    missing = self.command(command + [directory + '/missing'], ROOT)
                    self.assertEqual(missing.returncode, 2, missing.stderr)
                    conflict = self.command(command + [directory, '-q', '-v'], ROOT)
                    self.assertEqual(conflict.returncode, 2, conflict.stderr)
            for name in ['results.tsv', 'mms.tsv']:
                (Path(directory) / name).write_text('nt\n32\n')
            result = self.command([sys.executable, str(VERIFY / 'swift-stieber/plot.py'), directory], ROOT)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertIn('select --table', result.stderr)
            self.assertFalse((Path(directory) / 'verification.png').exists())


if __name__ == '__main__':
    unittest.main()
