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
    """Remove any sys.path entry that contains a ``pharokka.py`` file.

    When pip installs the ``pharokka.py`` backward-compat console_script,
    it writes ``bin/pharokka.py``.  Python then adds ``bin/`` to
    ``sys.path[0]`` whenever any script in that directory is executed,
    causing ``import pharokka`` to find the file ``bin/pharokka.py``
    instead of the installed ``pharokka`` package.

    Rather than guessing the offending directory from ``sys.argv[0]``
    (which is unreliable on HPC systems where argv[0] may be a relative
    path, a symlink, or a wrapper script), we scan every entry in sys.path
    and drop any directory that contains a ``pharokka.py`` file.

    Short-circuits immediately if ``sys.path[0]`` doesn't contain
    ``pharokka.py`` — the common case when called from the main
    ``pharokka`` entry point, which has no such file in its script
    directory.
    """
    # Fast path: nothing to do if no entry on sys.path shadows pharokka.
    if not any(
        p and os.path.isfile(os.path.join(p, "pharokka.py"))
        for p in sys.path
    ):
        return
    sys.path = [
        p for p in sys.path
        if not (p and os.path.isfile(os.path.join(p, "pharokka.py")))
    ]


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


def _legacy_shim(old: str, new: str, import_main):
    """Common boilerplate for the legacy script shims.

    Fixes sys.path, prints the deprecation notice, registers the
    loguru→sys.exit error sink (these shims bypass ``cli.main()`` which
    normally registers it), then calls the subcommand's ``main()``.
    """
    _fix_sys_path()
    _warn(old, new)
    from pharokka.util import register_error_sink
    register_error_sink()
    import_main()()


def pharokka_py() -> None:
    """pharokka.py  →  pharokka run"""
    def _imp():
        from pharokka.run import main
        return main
    _legacy_shim("pharokka.py", "pharokka run", _imp)


def pharokka_proteins_py() -> None:
    """pharokka_proteins.py  →  pharokka proteins"""
    def _imp():
        from pharokka.proteins import main
        return main
    _legacy_shim("pharokka_proteins.py", "pharokka proteins", _imp)


def install_databases_py() -> None:
    """install_databases.py  →  pharokka install"""
    def _imp():
        from pharokka.install import main
        return main
    _legacy_shim("install_databases.py", "pharokka install", _imp)


def pharokka_plotter_py() -> None:
    """pharokka_plotter.py  →  pharokka plot"""
    def _imp():
        from pharokka.plot_entry import main
        return main
    _legacy_shim("pharokka_plotter.py", "pharokka plot", _imp)


def pharokka_multiplotter_py() -> None:
    """pharokka_multiplotter.py  →  pharokka multiplot"""
    def _imp():
        from pharokka.multiplot_entry import main
        return main
    _legacy_shim("pharokka_multiplotter.py", "pharokka multiplot", _imp)


def create_custom_hmm_py() -> None:
    """create_custom_hmm.py  →  pharokka create-hmm"""
    def _imp():
        from pharokka.create_custom_hmm import main
        return main
    _legacy_shim("create_custom_hmm.py", "pharokka create-hmm", _imp)
