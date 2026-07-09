# Changelog

## v1.10.0 — CLI Redesign & Polars Refactor

### New subcommand-based CLI

`pharokka` v1.10.0 introduces a unified `pharokka` entry point with subcommands, replacing the collection of separate `.py` scripts:

```
pharokka <command> [options]

Commands:
  run         Annotate phage genome(s)
  proteins    Annotate proteins only
  install     Install/download databases
  plot        Plot a single phage genome
  multiplot   Plot multiple phages from a genbank file
  create-hmm  Create a custom HMM database from MSAs
```

All arguments are identical to the old scripts — only the invocation syntax changes.

#### Quick migration reference

| Old command (v1.9.x and earlier) | New command (v1.10.0+) |
|---|---|
| `pharokka.py -i ...` | `pharokka run -i ...` |
| `pharokka_proteins.py -i ...` | `pharokka proteins -i ...` |
| `install_databases.py -o ...` | `pharokka install -o ...` |
| `pharokka_plotter.py -i ...` | `pharokka plot -i ...` |
| `pharokka_multiplotter.py -g ...` | `pharokka multiplot -g ...` |
| `create_custom_hmm.py -i ...` | `pharokka create-hmm -i ...` |

### Backward compatibility — old script names still work

**You do not need to update your existing pipelines immediately.** All legacy script names are still installed alongside the new `pharokka` command. They forward every call transparently to the correct subcommand, with the only difference being a one-line deprecation notice on stderr:

```
[pharokka] DeprecationWarning: 'pharokka.py' is deprecated and will be
removed in a future release. Use 'pharokka run' instead.
```

This means existing workflows, Snakemake rules, Nextflow processes, shell scripts, and Galaxy wrappers that call `pharokka.py` (or any of the other old script names) will continue to run without modification.

The deprecated script names will be removed in a future major release.

### Other changes in v1.10.0

* **Polars rewrite**: the internal data processing pipeline has been rewritten from pandas to [polars](https://pola.rs) for improved performance and memory efficiency.
* **`src/` layout**: the package now uses a proper `src/` layout for cleaner installation and import isolation.
* **Numeric contig header bug fix**: tmRNA and CRISPR features were silently not counted for contigs with purely numeric headers (e.g. those produced by Unicycler). This is now fixed.
* **Meta mode coordinate bug fix**: `top_hits_vfdb.tsv` and `top_hits_card.tsv` reported inverted coordinates (start > stop) for negative-strand CDS hits in meta mode. This was a pandas reference side-effect — `locus_df` captured coordinates before the GFF positional swap, so negative-strand starts and stops were reversed. Both modes now correctly report GFF-order coordinates (start ≤ stop, strand indicates direction).
* **tRNA anticodon position bug fix**: the `/anticodon=(pos:X..Y,…)` qualifier in `pharokka.gff`, `pharokka.gbk` and `pharokka.tbl` was attached to the wrong tRNA on any contig with more than one tRNA. tRNAscan-SE's `.sec` output is sorted by tRNA ID; the GFF is sorted by genomic position; the v1.9.1 implementation matched the two by row index, so positions ended up scrambled. Now keyed by the tRNAscan ID (e.g. `MW460250_1.trna3`) for unambiguous matching.
* **No-tRNA crash fix**: phages with no tRNAs (or CRISPRs/tmRNAs) produce an empty `trnascan_out.gff`. Reading a 0-byte file with polars raised `OSError: Exec format error (os error 8)` on Linux (and `NoDataError` on other platforms), crashing post-processing. Empty or missing feature GFFs are now handled by a single `read_feature_gff` helper that returns an empty table.
* **Numerical precision**: float handling has been cleaned up throughout for more consistent output.
* **Parallel external tools**: tRNAscan-SE, MinCED, Aragorn and mash now run via a thread pool rather than sequential batch-and-wait loops, speeding up runs.
* **Test portability**: test database path now resolves via `$PHAROKKA_DB` environment variable (defaults to `tests/test_data/database`), removing hardcoded developer paths.
* **[pixi](https://pixi.sh) support**: `pharokka` can now be installed and developed reproducibly with [pixi](https://pixi.sh). `pixi.lock` pins the full dependency stack (Python libraries + external tools) across `linux-64`, `linux-aarch64`, `osx-64` and `osx-arm64`, and `pixi shell` sets up `pharokka` with all of its dependencies from source. See the [installation guide](install.md#source).
* **Developer tooling**: continuous integration migrated to pixi, with `ruff` formatting/linting and golden-output regression tests added.

---

## v1.9.0 (12 January 2026)

* Adds `pyrodigal-rv` as a gene predictor option for RNA phages (`-g pyrodigal-rv`)
* Fixes bug with incorrect translation table in `-g prodigal` + meta mode
* Adds `--reverse_mmseqs2` flag (PHROG MMseqs2 profile database as target rather than query; only recommended for enormous datasets)
* Adds `--sensitivity` CLI option to control MMseqs2 profile search sensitivity

## v1.8.0 (14 September 2025)

**Database reinstall required** — the MMseqs2 PHROG profile database format changed for compatibility with MMseqs2 v14+.

* Rebuilt MMseqs2 PHROG profile database
* Updated Phanotate to v1.6.7 (Python 3.13 support)
* Updated INPHARED Mash sketch database to 9 Aug 2025 release

For full release notes, see [HISTORY.md](https://github.com/gbouras13/pharokka/blob/main/HISTORY.md).
