"""Exercise the actual embedded farm scripts without GitHub or Docker."""
import base64
import gzip
import io
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import textwrap
import unittest

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / ".github/workflows/campaign-farm.yml"
MACRO = "/remoll/seed @@SEED@@\n/remoll/filename @@OUTPUT@@\n" + "".join(
    f"/remoll/SD/enable {det}\n" for det in (32, 33, 34, 35)) + "/run/beamOn 10\n"


def script(marker):
    # Read one literal run block; no YAML dependency is needed for these tests.
    lines = WORKFLOW.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == marker)
    start = next(i for i in range(start + 1, len(lines)) if lines[i].strip() == "run: |") + 1
    end = start
    while end < len(lines) and (not lines[end].strip() or lines[end].startswith("          ")):
        end += 1
    return textwrap.dedent("\n".join(lines[start:end])) + "\n"


def bundle(items):
    data = io.BytesIO()
    with tarfile.open(fileobj=data, mode="w:gz") as archive:
        for name, content in items:
            member = tarfile.TarInfo(name)
            payload = content.encode()
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))
    return base64.b64encode(data.getvalue()).decode()


class FarmScripts(unittest.TestCase):
    def materialize(self, macro=MACRO, **updates):
        with tempfile.TemporaryDirectory() as directory:
            env = dict(os.environ, MODE="macro", INDEX="0", SEED_BASE="710000", LABEL="check",
                       MACRO=base64.b64encode(gzip.compress(macro.encode())).decode() + "\n",
                       BUNDLE="", MEMBER="", GITHUB_OUTPUT=directory + "/outputs")
            env.update(updates)
            result = subprocess.run(["bash", "-e", "-c", script("- name: Materialize this job's macro")],
                                    cwd=directory, env=env, text=True, capture_output=True)
            path = Path(directory) / "job/run.mac"
            output = Path(directory) / "outputs"
            self.outputs = output.read_text() if output.exists() else ""
            return result.returncode, path.read_text() if path.exists() else result.stderr

    def test_distinct_worker_seeds(self):
        for index in (0, 1):
            code, text = self.materialize(INDEX=str(index))
            self.assertEqual(code, 0, text)
            self.assertIn(f"/remoll/seed {710000 + index}\n", text)
            self.assertIn(f"/work/job/check_w{index:02d}.root", text)

    def test_fixed_missing_comment_only_and_duplicate_seed(self):
        for bad in (MACRO.replace("@@SEED@@", "42"),
                    MACRO.replace("/remoll/seed @@SEED@@", "# @@SEED@@"),
                    MACRO + "/remoll/seed 42\n"):
            self.assertNotEqual(self.materialize(bad)[0], 0)

    def test_last_filter_setting_before_run(self):
        for commands, expected in (("/remoll/evgen/transportNeymanEnable YES\n", 1),
                                   ("/remoll/evgen/transportNeymanEnable true\n/remoll/evgen/transportNeymanEnable false\n", 0)):
            code, text = self.materialize(commands + MACRO)
            self.assertEqual(code, 0, text)
            self.assertIn(f"filtered={expected}\n", self.outputs)

    def test_missing_output_or_multiple_runs(self):
        for bad in (MACRO.replace("@@OUTPUT@@", "fixed.root"),
                    MACRO.replace("/run/beamOn 10", "/run/beamOn 0"),
                    MACRO + "/run/beamOn 10\n"):
            self.assertNotEqual(self.materialize(bad)[0], 0)

    def test_bundle_with_plan_and_trailing_newline(self):
        code, text = self.materialize(MODE="bundle", MEMBER="check.mac",
                                     BUNDLE=bundle([("check.mac", MACRO), ("check.neyman.tsv", "plan")]) + "\n")
        self.assertEqual(code, 0, text)
        self.assertIn("/work/job/check.root", text)

    def test_unsafe_or_duplicate_bundle_names(self):
        for items in ([("../check.mac", MACRO)], [("check.mac", MACRO), ("check.mac", MACRO)]):
            self.assertNotEqual(self.materialize(MODE="bundle", MEMBER="check.mac", BUNDLE=bundle(items))[0], 0)

    def test_seed_after_run_is_rejected(self):
        text = MACRO.replace("/remoll/seed @@SEED@@\n", "") + "/remoll/seed @@SEED@@\n"
        self.assertNotEqual(self.materialize(text)[0], 0)

    def test_bundled_output_cannot_fake_success(self):
        payload = bundle([("check.mac", "/thisCommandDoesNotExist\n" + MACRO),
                          ("check.root", "this is not ROOT output")])
        self.assertNotEqual(self.materialize(MODE="bundle", MEMBER="check.mac", BUNDLE=payload)[0], 0)

    def test_run_clears_old_output_and_rejects_command_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / "job").mkdir()
            (path / "job/check.root").write_text("stale output")
            (path / "docker").write_text("#!/bin/bash\nif [[ \"$1\" == run ]]; then echo 'COMMAND NOT FOUND'; fi\n")
            (path / "docker").chmod(0o755)
            env = dict(os.environ, PATH=directory + os.pathsep + os.environ["PATH"],
                       GITHUB_WORKSPACE=directory, IMAGE="example@sha256:" + "a" * 64,
                       JOB_NAME="check", THREADS="2")
            result = subprocess.run(["bash", "-e", "-c", script("- name: Run remoll")],
                                    cwd=directory, env=env, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((path / "job/check.root").exists())

    def test_one_image_digest_for_all_workers(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            digest = "sha256:" + "a" * 64
            (path / "docker").write_text("#!/bin/bash\nprintf 'Digest: " + digest +
                                        "\\n  Digest: sha256:child\\n'\n")
            (path / "docker").chmod(0o755)
            env = dict(os.environ, PATH=directory + os.pathsep + os.environ["PATH"],
                       REPOSITORY="gahljust/remoll", IMAGE_TAG="gahljust-working", GITHUB_OUTPUT=directory + "/output")
            result = subprocess.run(["bash", "-e", "-c", script("- name: Pin the campaign image")],
                                    cwd=directory, env=env, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((path / "output").read_text(), "ref=ghcr.io/gahljust/remoll/remoll-campaign@" + digest + "\n")
        self.assertIn("IMAGE: ${{ needs.plan.outputs.image }}", WORKFLOW.read_text())
        self.assertNotIn("IMAGE_TAG", script("- name: Run remoll"))

    def test_invalid_seed_and_label(self):
        for value in ("0", "-2", "2147483648", "1+2"):
            self.assertNotEqual(self.materialize(SEED_BASE=value)[0], 0)
        self.assertNotEqual(self.materialize(LABEL="$(touch unwanted)")[0], 0)

    def test_no_input_interpolation_in_run_shell(self):
        self.assertNotIn("${{", script("- name: Run remoll"))
        self.assertNotIn("${{", script("- name: Validate the saved ROOT output"))

    def test_docker_context_isolation(self):
        self.assertIn(".git", (ROOT / ".dockerignore").read_text().splitlines())
        self.assertNotIn(".git", (ROOT / "Dockerfile.campaign.dockerignore").read_text().splitlines())
        self.assertIn("persist-credentials: false", (ROOT / ".github/workflows/campaign-image.yml").read_text())
        self.assertIn("rm -rf /jlab/remoll/.git", (ROOT / "Dockerfile.campaign").read_text())


if __name__ == "__main__":
    unittest.main()
