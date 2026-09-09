"""Regression guard for Issue #51's host-user Docker execution contract.

Not a `bin/` CLI script test -- what this pins is `nextflow.config`'s
`docker.runOptions`, which is plain text (Groovy DSL) with no CLI surface
to exercise. `tests/modules/docker_host_user_ownership.nf.test` covers the
runtime half of the contract (a container task's output really is owned by
the launching user) but only for the profile nf-test runs under,
`test,docker`. The `docker_amd64` half cannot be covered that way: this
project's CI runs on linux/amd64, where `--platform linux/amd64` changes
nothing observable, so an emulation-path regression would pass a runtime
test silently.

The specific regression these tests exist for is the one Issue #51 called
out: `docker.runOptions` is a plain config key, so a profile that assigns
it replaces the value outright rather than appending to it. `-u` set once
outside `profiles` is therefore dropped the moment `docker_amd64` assigns
`--platform linux/amd64` to the same key, which would leave the Apple
Silicon path writing root-owned files again while the amd64 path stayed
fixed. Nextflow 26.04 offers no `docker.platform` setting to carry the
platform separately (assigning one parses and is then ignored, emitting no
`--platform` at all), and the config syntax rejects a shared variable
declaration ("Variable declarations cannot be mixed with config
statements"), so the option string is deliberately written twice. These
tests are what keeps the two copies in agreement.

Run with: python3 -m unittest discover -s tests/bin -v
"""

from __future__ import annotations

import re
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# The exact option string, `$(id -u)`/`$(id -g)` included: the IDs must be
# resolved by the task launcher's shell at run time, so that the owner is
# whoever actually started the run. A literal UID baked in here would be
# wrong on every other host.
HOST_USER_OPTION = "-u $(id -u):$(id -g)"

PLATFORM_OPTION = "--platform linux/amd64"

RUN_OPTIONS_ASSIGNMENT = re.compile(r"docker\.runOptions\s*=\s*'([^']*)'")


def _config_text() -> str:
    return (REPO_ROOT / "nextflow.config").read_text(encoding="utf-8")


def _extract_profile_block(config_text: str, profile: str) -> str:
    """Return the `<profile> { ... }` block's raw text, braces included."""
    marker = f"\n    {profile} {{"
    start = config_text.index(marker) + 1
    depth = 0
    for index in range(start, len(config_text)):
        if config_text[index] == "{":
            depth += 1
        elif config_text[index] == "}":
            depth -= 1
            if depth == 0:
                return config_text[start : index + 1]
    raise AssertionError(f"unbalanced braces in profile block: {profile}")


class DockerHostUserOptionTests(unittest.TestCase):
    """Every `docker.runOptions` value must carry the host-user option."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.config_text = _config_text()

    def test_run_options_is_assigned_outside_any_profile(self) -> None:
        # The base assignment has to precede `profiles {` so that it
        # applies to every Docker run, including profile combinations
        # this repository does not document.
        profiles_start = self.config_text.index("\nprofiles {")
        base_assignments = RUN_OPTIONS_ASSIGNMENT.findall(
            self.config_text[:profiles_start]
        )
        self.assertEqual(
            base_assignments,
            [HOST_USER_OPTION],
            "nextflow.config no longer sets the host-user docker.runOptions "
            "outside profiles",
        )

    def test_every_run_options_assignment_keeps_the_host_user_option(self) -> None:
        # This is the guard against the last-assignment-wins trap: any
        # new or edited assignment that forgets `-u` silently reverts
        # root ownership for whichever profile selection reaches it.
        assignments = RUN_OPTIONS_ASSIGNMENT.findall(self.config_text)
        self.assertGreaterEqual(len(assignments), 2)
        for value in assignments:
            with self.subTest(run_options=value):
                self.assertIn(HOST_USER_OPTION, value)

    def test_docker_amd64_keeps_both_the_host_user_and_platform_options(self) -> None:
        block = _extract_profile_block(self.config_text, "docker_amd64")
        assignments = RUN_OPTIONS_ASSIGNMENT.findall(block)
        self.assertEqual(len(assignments), 1)
        self.assertIn(HOST_USER_OPTION, assignments[0])
        self.assertIn(PLATFORM_OPTION, assignments[0])

    def test_platform_option_is_still_carried_by_run_options(self) -> None:
        # Nextflow 26.04 accepts a `docker.platform` setting and then
        # ignores it -- no `--platform` reaches the emitted `docker run`
        # -- so moving the platform out of `runOptions` would disable
        # amd64 emulation without any error. Both the fact that
        # `--platform` appears only in `runOptions` and that no
        # `docker.platform` key exists are pinned here.
        settings = [
            line for line in self.config_text.splitlines()
            if not line.strip().startswith("//")
        ]
        self.assertEqual([line for line in settings if "docker.platform" in line], [])
        platform_lines = [line for line in settings if PLATFORM_OPTION in line]
        self.assertEqual(len(platform_lines), 1)
        self.assertIn("docker.runOptions", platform_lines[0])

    def test_docker_profile_does_not_reassign_run_options(self) -> None:
        # The plain `docker` profile inherits the base assignment. If it
        # started assigning the key too, `-profile docker_amd64,docker`
        # would drop `--platform` again -- the mirror image of the bug
        # this Issue fixed.
        block = _extract_profile_block(self.config_text, "docker")
        self.assertNotIn("docker.runOptions", block)


if __name__ == "__main__":
    unittest.main()
