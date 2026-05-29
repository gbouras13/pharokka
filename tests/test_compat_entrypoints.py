"""
Smoke tests for backward-compatible entry-point scripts.

These tests verify that each legacy script name (pharokka.py,
pharokka_proteins.py, etc.) is importable, prints a deprecation warning to
stderr, and exits cleanly when called with --help / -h.

No full pipeline is run — we only confirm the entry points are wired up
correctly and forward to the right subcommand.
"""

import subprocess
import sys

import pytest


def _run(cmd: list[str]) -> subprocess.CompletedProcess:
    """Run a command and return the result (stdout + stderr captured)."""
    return subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def assert_deprecation_warning(result: subprocess.CompletedProcess, old_name: str, new_cmd: str) -> None:
    """Check that the deprecation warning was emitted on stderr."""
    assert old_name in result.stderr, (
        f"Expected deprecation warning mentioning '{old_name}' in stderr.\n"
        f"stderr was:\n{result.stderr}"
    )
    assert new_cmd in result.stderr, (
        f"Expected new command '{new_cmd}' mentioned in stderr.\n"
        f"stderr was:\n{result.stderr}"
    )


def assert_help_shown(result: subprocess.CompletedProcess) -> None:
    """Check that --help produced usage output (not a traceback)."""
    combined = result.stdout + result.stderr
    assert "usage:" in combined.lower(), (
        f"Expected 'usage:' in output but got:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
    assert result.returncode == 0, (
        f"Expected exit code 0 but got {result.returncode}.\n"
        f"stdout: {result.stdout}\nstderr: {result.stderr}"
    )


# ---------------------------------------------------------------------------
# pharokka.py  →  pharokka run
# ---------------------------------------------------------------------------

def test_pharokka_py_help():
    result = _run(["pharokka.py", "--help"])
    assert_deprecation_warning(result, "pharokka.py", "pharokka run")
    assert_help_shown(result)


def test_pharokka_py_importable():
    """pharokka_scripts module must be importable (entry-point wiring)."""
    result = _run([sys.executable, "-c", "from pharokka_scripts import pharokka_py"])
    assert result.returncode == 0, f"Import failed:\n{result.stderr}"


# ---------------------------------------------------------------------------
# pharokka_proteins.py  →  pharokka proteins
# ---------------------------------------------------------------------------

def test_pharokka_proteins_py_help():
    result = _run(["pharokka_proteins.py", "--help"])
    assert_deprecation_warning(result, "pharokka_proteins.py", "pharokka proteins")
    assert_help_shown(result)


# ---------------------------------------------------------------------------
# install_databases.py  →  pharokka install
# ---------------------------------------------------------------------------

def test_install_databases_py_help():
    result = _run(["install_databases.py", "--help"])
    assert_deprecation_warning(result, "install_databases.py", "pharokka install")
    assert_help_shown(result)


# ---------------------------------------------------------------------------
# pharokka_plotter.py  →  pharokka plot
# ---------------------------------------------------------------------------

def test_pharokka_plotter_py_help():
    result = _run(["pharokka_plotter.py", "--help"])
    assert_deprecation_warning(result, "pharokka_plotter.py", "pharokka plot")
    assert_help_shown(result)


# ---------------------------------------------------------------------------
# pharokka_multiplotter.py  →  pharokka multiplot
# ---------------------------------------------------------------------------

def test_pharokka_multiplotter_py_help():
    result = _run(["pharokka_multiplotter.py", "--help"])
    assert_deprecation_warning(result, "pharokka_multiplotter.py", "pharokka multiplot")
    assert_help_shown(result)


# ---------------------------------------------------------------------------
# create_custom_hmm.py  →  pharokka create-hmm
# ---------------------------------------------------------------------------

def test_create_custom_hmm_py_help():
    result = _run(["create_custom_hmm.py", "--help"])
    assert_deprecation_warning(result, "create_custom_hmm.py", "pharokka create-hmm")
    assert_help_shown(result)


# ---------------------------------------------------------------------------
# cli.py invocation-name detection (symlink / alias path)
# ---------------------------------------------------------------------------

def test_cli_detects_legacy_name_via_argv(tmp_path):
    """
    Simulate the cli.py basename-detection path: call pharokka.cli.main()
    with sys.argv[0] set to 'pharokka.py' and '--help' as the argument.
    The dispatcher should inject 'run', print the deprecation warning, and
    show pharokka-run help.
    """
    script = tmp_path / "smoke.py"
    script.write_text(
        "import sys\n"
        "sys.argv = ['pharokka.py', '--help']\n"
        "from pharokka.cli import main\n"
        "main()\n"
    )
    result = _run([sys.executable, str(script)])
    assert "pharokka.py" in result.stderr
    assert "pharokka run" in result.stderr
    combined = result.stdout + result.stderr
    assert "usage:" in combined.lower()
    assert result.returncode == 0
