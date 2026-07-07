"""
pytest configuration shared across all test modules.

Installs the loguru → sys.exit(1) sink so that unit tests which call
pharokka validation functions directly (e.g. validate_fasta) get the same
SystemExit behaviour as a real pharokka run.  Without this, logger.error()
would only log to stderr and the assertRaises(SystemExit) assertions in
test_input_commands.py would never fire.
"""

import sys

import pytest
from loguru import logger


@pytest.fixture(autouse=True, scope="session")
def loguru_exit_on_error():
    """Add a loguru sink that calls sys.exit(1) on ERROR for the whole test run.

    Scoped to ``session`` rather than the default per-function so the sink
    is registered once and torn down once across all tests — the sink has
    no per-test state, so re-installing for every test is pure overhead.
    """
    sink_id = logger.add(lambda _: sys.exit(1), level="ERROR")
    yield
    logger.remove(sink_id)
