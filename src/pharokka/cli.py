"""
pharokka CLI dispatcher.

Provides the `pharokka` entry point with subcommands:

  pharokka run        - annotate phage genome(s)
  pharokka proteins   - annotate proteins only
  pharokka install    - install/download databases
  pharokka plot       - plot a single phage genome
  pharokka multiplot  - plot multiple phages from a genbank file
  pharokka create-hmm - build a custom HMM database from MSAs
"""

import os
import sys

# ---------------------------------------------------------------------------
# Old script name → subcommand mapping (for backward compatibility)
# ---------------------------------------------------------------------------
# If pharokka is invoked via a symlink or alias that uses a legacy script
# name (e.g. 'pharokka.py'), we inject the correct subcommand so that the
# call is transparently forwarded to the right handler.
_LEGACY_NAMES: dict = {
    "pharokka.py":              "run",
    "pharokka_proteins.py":     "proteins",
    "install_databases.py":     "install",
    "pharokka_plotter.py":      "plot",
    "pharokka_multiplotter.py": "multiplot",
    "create_custom_hmm.py":     "create-hmm",
}

USAGE = """\
Usage: pharokka <command> [options]

Commands:
  run         Annotate phage genome(s)
  proteins    Annotate proteins only
  install     Install/download databases
  plot        Plot a single phage genome
  multiplot   Plot multiple phages from a genbank file
  create-hmm  Create a custom HMM database from MSAs

Run `pharokka <command> --help` for help on a specific command.
"""


def main():
    # Make every logger.error() call exit with status 1.
    # Matches the behaviour of the old pharokka.py entry point
    # (https://github.com/gbouras13/plassembler/pull/69).
    from loguru import logger as _logger
    _logger.add(lambda _: sys.exit(1), level="ERROR")

    # Backward compat: detect legacy invocation names and inject subcommand.
    # This covers users who have a symlink  pharokka.py → pharokka  or an
    # alias, rather than the installed compat script in scripts/.
    _basename = os.path.basename(sys.argv[0])
    if _basename in _LEGACY_NAMES:
        _sub = _LEGACY_NAMES[_basename]
        print(
            f"[pharokka] DeprecationWarning: '{_basename}' is deprecated and will be "
            f"removed in a future release. Use 'pharokka {_sub}' instead.",
            file=sys.stderr,
        )
        sys.argv.insert(1, _sub)

    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print(USAGE)
        sys.exit(0)

    if sys.argv[1] in ("-V", "--version", "version"):
        from pharokka.version import __version__
        print(f"pharokka v{__version__}")
        sys.exit(0)

    subcommand = sys.argv[1]

    # Remove the subcommand token so each module's own argparse parser
    # receives only its own arguments (sys.argv[0] stays as the program name).
    sys.argv.pop(1)

    if subcommand == "run":
        from pharokka.run import main as _main
        _main()

    elif subcommand == "proteins":
        from pharokka.proteins import main as _main
        _main()

    elif subcommand == "install":
        from pharokka.install import main as _main
        _main()

    elif subcommand == "plot":
        from pharokka.plot_entry import main as _main
        _main()

    elif subcommand == "multiplot":
        from pharokka.multiplot_entry import main as _main
        _main()

    elif subcommand in ("create-hmm", "create_hmm"):
        from pharokka.create_custom_hmm import main as _main
        _main()

    else:
        print(f"pharokka: unknown command '{subcommand}'\n", file=sys.stderr)
        print(USAGE, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
