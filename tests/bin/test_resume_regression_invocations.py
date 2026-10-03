"""Each resume assertion must read that invocation's trace, preserving history."""

import importlib.util
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/check_issue69_resume.py"
SPEC = importlib.util.spec_from_file_location("resume_regression", SCRIPT)
resume_regression = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(resume_regression)


class ResumeInvocationTests(unittest.TestCase):
    def test_changed_input_check_cannot_read_the_previous_cached_trace(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            statuses = iter(["CACHED", "COMPLETED"])
            log_names = []

            def nextflow(command, **kwargs):
                trace = Path(command[command.index("-with-trace") + 1])
                # Reproduce Nextflow's refusal to overwrite an existing trace.
                if not trace.exists():
                    trace.write_text(
                        "name\tstatus\nADZUKI_SNP_PIPELINE:FASTP (sample:group)\t"
                        + next(statuses)
                        + "\n"
                    )
                log_names.append(kwargs["stdout"].name)
                return subprocess.CompletedProcess(command, 0)

            with patch.object(resume_regression.subprocess, "run", side_effect=nextflow):
                first = resume_regression.run_pipeline(
                    root, "interrupted", root / "samples.csv", resume=True
                )[2]
                second = resume_regression.run_pipeline(
                    root, "interrupted", root / "samples.csv", resume=True
                )[2]
            self.assertNotEqual(first, second)
            self.assertIn("\tCACHED\n", first.read_text())
            self.assertIn("\tCOMPLETED\n", second.read_text())
            self.assertNotEqual(log_names[0], log_names[1])
            self.assertTrue(all(Path(name).is_file() for name in log_names))


if __name__ == "__main__":
    unittest.main()
