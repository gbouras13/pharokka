"""
pharokka_scripts — backward-compatibility console_scripts shims.

Why this module exists (and why it is NOT inside the pharokka package)
-----------------------------------------------------------------------
pip generates console_script wrappers that look like:

    #!/path/to/python
    from <module> import <func>
    <func>()

When the wrapper is named ``pharokka.py`` (to match the old CLI), Python
adds ``bin/`` to ``sys.path[0]`` before executing it.  Because the script
file is *named* ``pharokka.py``, a top-level ``from pharokka.X import ...``
would find the script itself instead of the installed package, producing:

    ModuleNotFoundError: No module named 'pharokka.X'; 'pharokka' is not a package

Putting the shim functions here (module name ``pharokka_scripts``) avoids
the collision entirely — there is no ``pharokka_scripts.py`` in ``bin/``,
so Python resolves the import from site-packages without issue.

Each shim then removes ``bin/`` from ``sys.path`` *before* importing from
``pharokka.*``, so subsequent imports are also safe.
"""

import os
import sys


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _fix_sys_path() -> None:
    """Remove the entry-point script's directory from sys.path.

    When pip installs e.g. ``bin/pharokka.py``, running that script adds
    ``bin/`` as ``sys.path[0]``.  Any ``.py`` file in ``bin/`` with the
    same stem as a package will shadow that package.  Removing ``bin/``
    before importing from ``pharokka.*`` restores correct resolution.
    """
    if not sys.path:
        return
    # sys.argv[0] is the entry-point script that was executed
    script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    sys.path = [p for p in sys.path if os.path.abspath(p) != script_dir]


def _warn(old: str, new: str) -> None:
    print(
        f"[pharokka] DeprecationWarning: '{old}' is deprecated and will be "
        f"removed in a future release. Use '{new}' instead.",
        file=sys.stderr,
    )


# ---------------------------------------------------------------------------
# Shim entry points (registered in pyproject.toml [project.scripts])
# ---------------------------------------------------------------------------

def pharokka_main() -> None:
    """pharokka  →  pharokka.cli:main (via sys.path-safe shim).

    The main ``pharokka`` entry point is also routed through this module so
    that ``bin/pharokka.py`` (the legacy backward-compat script) does not
    shadow the ``pharokka`` package when ``bin/`` is on ``sys.path[0]``.
    """
    _fix_sys_path()
    from pharokka.cli import main
    main()


def pharokka_py() -> None:
    """pharokka.py  →  pharokka run"""
    _fix_sys_path()
    _warn("pharokka.py", "pharokka run")
    from pharokka.run import main
    main()


def pharokka_proteins_py() -> None:
    """pharokka_proteins.py  →  pharokka proteins"""
    _fix_sys_path()
    _warn("pharokka_proteins.py", "pharokka proteins")
    from pharokka.proteins import main
    main()


def install_databases_py() -> None:
    """install_databases.py  →  pharokka install"""
    _fix_sys_path()
    _warn("install_databases.py", "pharokka install")
    from pharokka.install import main
    main()


def pharokka_plotter_py() -> None:
    """pharokka_plotter.py  →  pharokka plot"""
    _fix_sys_path()
    _warn("pharokka_plotter.py", "pharokka plot")
    from pharokka.plot_entry import main
    main()


def pharokka_multiplotter_py() -> None:
    """pharokka_multiplotter.py  →  pharokka multiplot"""
    _fix_sys_path()
    _warn("pharokka_multiplotter.py", "pharokka multiplot")
    from pharokka.multiplot_entry import main
    main()


def create_custom_hmm_py() -> None:
    """create_custom_hmm.py  →  pharokka create-hmm"""
    _fix_sys_path()
    _warn("create_custom_hmm.py", "pharokka create-hmm")
    from pharokka.create_custom_hmm import main
    main()
