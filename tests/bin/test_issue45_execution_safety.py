"""Failure and restoration contracts for the opt-in Issue45 execution helper."""

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

HELPERS = Path(__file__).resolve().parents[2] / "benchmarks/issue45"
sys.path.insert(0, str(HELPERS))
import run_benchmarks as benchmark_runner  # noqa: E402
import run_generation as runner  # noqa: E402


class GenerationRestorationTests(unittest.TestCase):
    def prepare(self, directory):
        keys = {
            "generate.nf": "staged_wrapper_sha256",
            "samples.csv": "samplesheet_sha256",
            "nextflow.config": "nextflow_config_sha256",
        }
        manifest = {}
        (directory / "bin").mkdir()
        helper = directory / "bin/synthetic_helper.py"
        helper.write_text("synthetic fixture")
        manifest["production_bin_sha256"] = {
            helper.name: hashlib.sha256(helper.read_bytes()).hexdigest()
        }
        for name, key in keys.items():
            data = name.encode()
            (directory / name).write_bytes(data)
            manifest[key] = hashlib.sha256(data).hexdigest()
        (directory / "production_lineage_manifest.json").write_text(json.dumps(manifest))
        return {
            "Id": "a" * 64,
            "Image": "sha256:" + "b" * 64,
            "State": {"Running": True, "Paused": False},
            "Config": {"Env": ["TOKEN=do-not-record-this"], "Image": "trial"},
            "HostConfig": {
                "PortBindings": {"1919/tcp": [{"HostIp": "127.0.0.1", "HostPort": "1919"}]}
            },
        }

    def invoke(self, directory, before, memory_values, stop_error=False):
        child = Mock()
        child.poll.return_value = 1
        child.returncode = 1
        calls = []

        def execute(arguments, **kwargs):
            calls.append(arguments)
            if stop_error and arguments[1] == "stop":
                raise subprocess.CalledProcessError(1, arguments)
            return subprocess.CompletedProcess(arguments, 0)

        with (
            patch.object(
                sys,
                "argv",
                ["run_generation.py", "--run-dir", str(directory), "--expected-trial-id", "a" * 12],
            ),
            patch.object(runner, "inspect", return_value=before),
            patch.object(runner, "serving", return_value=True),
            patch.object(runner, "memory", side_effect=memory_values),
            patch.object(runner.signal, "signal"),
            patch.object(runner.subprocess, "check_output", return_value=""),
            patch.object(runner.subprocess, "run", side_effect=execute),
            patch.object(runner.subprocess, "Popen", return_value=child),
            patch.object(runner.shutil, "disk_usage", return_value=Mock(free=4_000_000_000_000)),
        ):
            with self.assertRaises((ValueError, RuntimeError, subprocess.CalledProcessError)):
                runner.main()
        return calls

    def test_nextflow_failure_restores_exact_container_and_does_not_record_secrets(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            before = self.prepare(directory)
            calls = self.invoke(directory, before, [(120 * 1024**3, 0)] * 3)
            record_text = (directory / "resource-preparation.private.json").read_text()
            record = json.loads(record_text)
            self.assertTrue(record["restored"])
            self.assertEqual(1, record["nextflow_exit_code"])
            self.assertIn(["docker", "start", before["Id"]], calls)
            self.assertNotIn("do-not-record-this", record_text)

    def test_insufficient_headroom_after_stop_still_restores_trial(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            before = self.prepare(directory)
            calls = self.invoke(directory, before, [(3 * 1024**3, 0)] * 3)
            record = json.loads((directory / "resource-preparation.private.json").read_text())
            self.assertTrue(record["restored"])
            self.assertIn("Less than 110 GiB", record["failure"])
            self.assertEqual(["docker", "start", before["Id"]], calls[-1])

    def test_failed_stop_command_still_attempts_original_running_state(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            before = self.prepare(directory)
            calls = self.invoke(directory, before, [(3 * 1024**3, 0)] * 2, stop_error=True)
            self.assertEqual(["docker", "start", before["Id"]], calls[-1])
            record = json.loads((directory / "resource-preparation.private.json").read_text())
            self.assertTrue(record["restored"])

    def test_changed_launch_file_is_rejected_before_stopping_trial(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            before = self.prepare(directory)
            (directory / "nextflow.config").write_text("changed config")
            calls = self.invoke(directory, before, [])
            self.assertEqual([], calls)
            self.assertFalse((directory / "resource-preparation.private.json").exists())

    def test_changed_frozen_production_script_is_rejected_before_stopping_trial(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            before = self.prepare(directory)
            (directory / "bin/synthetic_helper.py").write_text("changed")
            calls = self.invoke(directory, before, [])
            self.assertEqual([], calls)
            self.assertFalse((directory / "resource-preparation.private.json").exists())


class LaunchHeadroomTests(unittest.TestCase):
    """A resumed suite must not be refused for space its own finished work occupies."""

    def guard(self, root, free):
        guard = benchmark_runner.ResourceGuard(root, None)
        guard.record = {}
        with patch.object(benchmark_runner.shutil, "disk_usage", return_value=Mock(free=free)):
            guard.check_launch_headroom(benchmark_runner.LAUNCH_MEMORY_GATE_BYTES)
        return guard.record

    def test_space_already_written_by_this_run_counts_toward_the_suite_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "E0").mkdir()
            (root / "E0/workspace.bin").write_bytes(b"x" * 4096)
            free = benchmark_runner.LAUNCH_STORAGE_GATE_BYTES - 1024
            with self.assertRaisesRegex(ValueError, "headroom gate"):
                self.guard(root, free - 8192)
            record = self.guard(root, free)
            self.assertEqual(free, record["storage_free_at_launch_bytes"])
            self.assertGreaterEqual(record["storage_budget_bytes"], free + 4096)

    def test_hard_linked_files_are_counted_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = root / "a.bin"
            original.write_bytes(b"y" * 2048)
            (root / "b.bin").hardlink_to(original)
            self.assertEqual(2048, benchmark_runner.occupied(root))

    def test_insufficient_memory_fails_the_gate_however_much_storage_exists(self):
        with tempfile.TemporaryDirectory() as tmp:
            guard = benchmark_runner.ResourceGuard(Path(tmp), None)
            guard.record = {}
            with (
                patch.object(
                    benchmark_runner.shutil, "disk_usage", return_value=Mock(free=9 * 10**12)
                ),
                self.assertRaisesRegex(ValueError, "headroom gate"),
            ):
                guard.check_launch_headroom(3 * 1024**3)


class ExecutingHelperTests(unittest.TestCase):
    """The lineage must cover every helper that can change a measurement, and nothing else."""

    def test_every_helper_is_classified_as_executing_or_reporting_only(self):
        present = {path.name for path in HELPERS.glob("*.py")}
        classified = set(benchmark_runner.EXECUTING_HELPERS) | set(
            benchmark_runner.REPORTING_ONLY_HELPERS
        )
        self.assertEqual(
            set(),
            present - classified,
            "a new helper must be declared as executing or reporting-only, so the lineage "
            "either covers it or documents why it does not",
        )
        self.assertEqual(set(), classified - present)
        self.assertEqual(
            set(),
            set(benchmark_runner.EXECUTING_HELPERS) & set(benchmark_runner.REPORTING_ONLY_HELPERS),
        )

    def test_modules_the_controller_imports_are_all_recorded_as_executing(self):
        imported = {
            Path(module.__file__).name
            for module in sys.modules.values()
            if getattr(module, "__file__", None)
            and Path(module.__file__).parent == HELPERS
            and Path(module.__file__).name not in set(benchmark_runner.REPORTING_ONLY_HELPERS)
        }
        self.assertTrue(imported)
        self.assertEqual(set(), imported - set(benchmark_runner.EXECUTING_HELPERS))


class LineageRecordTests(unittest.TestCase):
    def test_recorded_resource_policy_names_only_real_module_constants(self):
        import ast

        source = (HELPERS / "run_benchmarks.py").read_text()
        referenced = {
            node.attr
            for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Attribute)
            and isinstance(node.value, ast.Name)
            and node.value.id == "execute_experiments"
        }
        self.assertTrue(referenced)
        for name in sorted(referenced):
            with self.subTest(name=name):
                self.assertTrue(hasattr(benchmark_runner.execute_experiments, name))


class MeasurementTests(unittest.TestCase):
    def test_failed_command_records_exit_code_and_refuses_existing_measurement(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "measurement.json"
            args = [
                sys.executable,
                str(HELPERS / "measure_process.py"),
                str(output),
                sys.executable,
                "-c",
                "raise SystemExit(7)",
            ]
            first = subprocess.run(args, capture_output=True, check=False)
            self.assertEqual(7, first.returncode)
            original = output.read_bytes()
            self.assertEqual(7, json.loads(original)["exit_code"])
            second = subprocess.run(args, capture_output=True, check=False)
            self.assertNotEqual(0, second.returncode)
            self.assertEqual(original, output.read_bytes())


class BenchmarkRestorationTests(unittest.TestCase):
    def test_benchmark_failure_and_launch_headroom_failure_both_restore_trial(self):
        for available in (3 * 1024**3, 120 * 1024**3):
            with self.subTest(available=available), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                before = GenerationRestorationTests().prepare(root)
                calls = []

                def execute(arguments, **kwargs):
                    calls.append(arguments)
                    return subprocess.CompletedProcess(arguments, 0)

                with (
                    patch.object(benchmark_runner, "inspect", return_value=before),
                    patch.object(benchmark_runner, "serving", return_value=True),
                    patch.object(benchmark_runner, "memory", return_value=(available, 0)),
                    patch.object(benchmark_runner.subprocess, "check_output", return_value=""),
                    patch.object(benchmark_runner.subprocess, "run", side_effect=execute),
                    patch.object(benchmark_runner.threading, "Thread"),
                    patch.object(
                        benchmark_runner.shutil,
                        "disk_usage",
                        return_value=Mock(free=4_000_000_000_000),
                    ),
                ):
                    with self.assertRaises((RuntimeError, ValueError)):
                        with benchmark_runner.ResourceGuard(root, "a" * 12):
                            raise RuntimeError("synthetic benchmark failure")
                record_text = (root / "resource-preparation.private.json").read_text()
                self.assertTrue(json.loads(record_text)["restored"])
                self.assertEqual(["docker", "start", before["Id"]], calls[-1])
                self.assertNotIn("do-not-record-this", record_text)

    def test_no_trial_pause_measures_without_stopping_any_user_workload(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with (
                patch.object(benchmark_runner, "memory", return_value=(120 * 1024**3, 0)),
                patch.object(benchmark_runner, "inspect") as examine,
                patch.object(benchmark_runner.subprocess, "check_output", return_value=""),
                patch.object(benchmark_runner.subprocess, "run") as execute,
                patch.object(benchmark_runner.threading, "Thread"),
                patch.object(
                    benchmark_runner.shutil, "disk_usage", return_value=Mock(free=4_000_000_000_000)
                ),
            ):
                with benchmark_runner.ResourceGuard(root, None):
                    pass
            examine.assert_not_called()
            self.assertNotIn(
                "stop", [call[1] for call in execute.call_args_list if len(call[1]) > 1]
            )
            record = json.loads((root / "resource-preparation.private.json").read_text())
            self.assertEqual([], record["paused_workloads"])
            self.assertNotIn("stop_command", record)

    def test_no_trial_pause_still_enforces_the_launch_headroom_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with (
                patch.object(benchmark_runner, "memory", return_value=(3 * 1024**3, 0)),
                patch.object(benchmark_runner.subprocess, "check_output", return_value=""),
                patch.object(benchmark_runner.subprocess, "run"),
                patch.object(benchmark_runner.threading, "Thread"),
                patch.object(
                    benchmark_runner.shutil, "disk_usage", return_value=Mock(free=4_000_000_000_000)
                ),
            ):
                with self.assertRaisesRegex(ValueError, "headroom gate"):
                    with benchmark_runner.ResourceGuard(root, None):
                        self.fail("Must not enter")

    def test_active_request_prevents_any_trial_stop(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            before = GenerationRestorationTests().prepare(root)
            with (
                patch.object(benchmark_runner, "inspect", return_value=before),
                patch.object(
                    benchmark_runner.subprocess, "check_output", return_value="active connection"
                ),
                patch.object(benchmark_runner.subprocess, "run") as execute,
            ):
                with self.assertRaisesRegex(ValueError, "busy"):
                    with benchmark_runner.ResourceGuard(root, "a" * 12):
                        self.fail("Must not enter")
            execute.assert_not_called()


if __name__ == "__main__":
    unittest.main()
