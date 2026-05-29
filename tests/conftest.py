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


@pytest.fixture(autouse=True)
def loguru_exit_on_error():
    """Add a loguru sink that calls sys.exit(1) on ERROR, then clean it up."""
    sink_id = logger.add(lambda _: sys.exit(1), level="ERROR")
    yield
    logger.remove(sink_id)
