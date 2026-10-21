# Changelog

## v1.11.0 — ncRNA annotation with Infernal and Rfam

### ncRNA annotation, on by default

`pharokka run` now annotates non-coding RNAs by scanning the genome against
[Rfam](https://rfam.org) 15.1 (4,227 covariance models) with
[Infernal](http://eddylab.org/infernal/) `cmscan`.  This picks up structured
RNAs that pharokka previously could not see at all — riboswitches, ribozymes,
regulatory sRNAs, group I/II introns and leader elements.

It runs **by default**, because it is inexpensive for a phage genome: roughly
5 seconds for a typical 40 kb phage and 14 seconds for a 140 kb phage on 8
threads. `--threads` scales it well (4–5x on 8 cores) even for a single genome.

To turn it off:

```bash
pharokka run -i phage.fasta -o output -d database --skip_rfam
```

**In meta mode (`-m`) it is skipped by default**, since runtime scales with
assembly size at roughly 2 min/Mbp on 8 threads. Use `--meta_rfam` to run it
in meta mode anyway.

New outputs:

* `{prefix}_ncrna.tsv` — one row per ncRNA, with Rfam accession, family, type,
  coordinates, bit score and E-value.
* `{prefix}_cmscan.tblout` — the raw Infernal output.
* `ncRNA` features in the `.gff` and `.gbk`, and an `ncRNAs` row per contig in
  `{prefix}_cds_functions.tsv`.

By default, Rfam tRNA (RF00005) and tmRNA (RF00023) hits are discarded, since
tRNAscan-SE and ARAGORN already annotate these and are more sensitive on phage
sequence.  Pass `--rfam_keep_trna` to keep them.

**Rfam does not replace tRNAscan-SE, ARAGORN or MinCED** — it is purely
additive.

### vrRNA guide array detection (VIPR systems), on by default

`pharokka run` now calls **vrRNA guide arrays** de novo, with no database and no new dependency.

vrRNAs are the guides of [VIPR](https://www.biorxiv.org/content/10.64898/2026.10.15.759431v1), an
RNA-guided DNA silencing system encoded by phages (Lopez *et al.* 2026): a *vipr* effector and a
*vap* accessory gene next to a tandem array of 6-16 `GGYNN` units, `GGY` invariant and the `NN`
pairs guide-defining, closed by a structured 3' tail. No annotation tool saw this feature class
before. Because the motif is the model, detection needs nothing downloaded — the implementation is
pure Python and numpy in `src/pharokka/vrrna.py`.

New outputs:

* `{prefix}_vrrna.tsv` — one row per array, with coordinates, unit count, conformity, empirical
  p-value, tier, cluster membership, both tail measures, the guide dinucleotides and the tail
  sequence.
* `ncRNA` features with `product=vrRNA guide array` in the `.gff`, `.tbl` and `.gbk`, carrying all
  of the above in the `note`.
* a `vrRNA_arrays` row per contig in `{prefix}_cds_functions.tsv`, always written.

To turn it off: `--skip_vrrna`.

#### The null is the filter

`GGYNN` carries two obligate Gs, so raw calls accumulate with GC content rather than with
biology — on a 70% GC actinobacteriophage a naive caller fires constantly and every one of those
calls is noise. Each array is therefore scored against a null built by shuffling its own contig,
and the user-facing knob is the resulting p-value, **not** a minimum number of units:

```bash
pharokka run -i phage.fasta -o output -d database --vrrna_max_p 0.01
```

The default is `0.05`. `0.01` is the right setting for bulk or high-GC work — it is where the
remaining GC artefacts drop out, for roughly 1.9 points of recall. `--vrrna_shuffles` (default
100) sets how many shuffles, and the smallest attainable p-value is `1 / (shuffles + 1)`.
`--vrrna_null_model dinuc` offers an Altschul-Erikson dinucleotide-preserving shuffle as a
control; it is not the default because, measured, it is the *weaker* null.

#### Three things are reported, not gated on

* **the 3' tail**, because it is type-specific: informative when it matches and silent when it
  does not. `tail_identity` scores against the type II consensus, `tail_agreement` against the
  other arrays on the same contig.
* **clustering** (`cluster_id`, `cluster_size`), because three arrays within 5 kb is a much better
  VIPR locus than one. It is reported as architecture and never multiplied into a joint p-value —
  all of a contig's arrays are scored against the same null draw, so three arrays at p = 0.01 are
  three correlated observations, not p = 1e-6.
* **the tier**: `strict` calls come from the genome-wide pass at >= 6 units; `relaxed` calls are
  >= 4 units within 2 kb of a strict call, scored against a null drawn from that neighbourhood.

#### Accuracy and cost

Measured against the 585 published arrays of the Lopez *et al.* supplement: recall 584/585 =
99.8%, exact boundaries 84.1%, 0.282 discordant calls per 100 kb — against a published 342/348 =
98.3% at a higher false-call rate.

The scan is ~11 ms for a median phage genome (measured over 2,364 genomes, 173 Mb). The null
re-scans the contig 100 times and costs seconds, but only runs on a contig that produced at least
one raw call. On SUSP1 (90 kb) the whole feature took 4.1 s against PHANOTATE's 47 s in the same
run.

### Plotting

`pharokka plot` and `pharokka multiplot` now draw ncRNAs. Rfam hits (orange) and vrRNA guide
arrays (teal) get a ring of their own just inside the CDS ring, with one legend entry each and
labels thinned the way the tRNA/tmRNA/CRISPR labels are. They are not given the black "extras"
treatment because, unlike tRNAs, ncRNAs routinely sit *inside* a CDS — all three Rfam sRNAs of
NC_004617.1 and three of the four vrRNA arrays of NC_028808.2 are wholly contained in one — so at
the CDS radius they read as part of the gene they overlap. The ring is only drawn when the contig
has ncRNAs, so plots of genomes without any are unchanged.

Two plotting bugs fixed along the way:

* `pharokka multiplot` asked pyCirclize for CDS features only (its default), so **no** tRNA,
  tmRNA or CRISPR was ever plotted by the multiplotter, despite the code that filters for them.
* CRISPR label thinning looped over the tmRNA count while indexing into the CRISPR lists, which
  dropped CRISPR labels on a genome with fewer tmRNAs than CRISPRs and raised `IndexError` on one
  with more. The multiplotter also never drew CRISPR features themselves, listing `tmRNA` twice
  where `repeat_region` was meant.

### Requirements — action needed when upgrading

* **Infernal >= 1.1.4 must be installed** (`conda install -c bioconda infernal`).
  It is only checked when Rfam annotation will actually run.
* **The v1.11.0 database is required**, which adds the pressed Rfam covariance
  models. Re-run `pharokka install` to update.

Because ncRNA annotation is on by default, running v1.11.0 against a v1.10.x or
older database will fail with an explanatory error. Either update the database
or pass `--skip_rfam`, which restores the previous behaviour exactly.

### Other changes

* Fixed the database tarball filename being hardcoded to `v1.8.0` rather than
  derived from the database version.
* The PHROGs database version marker file is now derived from the database
  version instead of being hardcoded.

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
