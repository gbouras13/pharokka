"""De novo detection of vrRNA guide arrays (VIPR systems).

VIPR is an RNA-guided DNA silencing system encoded by phages and other mobile
elements (Lopez *et al.* 2026, bioRxiv 10.64898/2026.10.15.759431): a Vipr
effector plus an accessory Vap protein, sitting next to an array of vrRNA
guides.  The array is a run of 6-16 tandem ``GGYNN`` units - ``GGY`` is the
invariant scaffold, ``NN`` is guide-defining - followed by a structured 3'
tail.  No annotation tool currently calls them, and they need no database:
the motif is the model.

The caller here is a port of the standalone implementation that was
benchmarked against the 585 published arrays of the Lopez *et al.* supplement
(recall 584/585 = 99.8%, exact boundaries 84.1%, 0.282 discordant calls per
100 kb, against a published 342/348 = 98.3% at a higher false-call rate).
The algorithm functions below - ``periodicity_windows``, ``call_in_window``,
``resolve``, ``scan``, ``null_max_units`` and the two shufflers - are that
version unchanged, because every one of those numbers was measured on it.
Changing them silently invalidates the benchmark.

The pipeline, which is also the order the functions appear in:

    periodicity_windows -> call_in_window -> shuffled null (p-value)
      -> relaxed second tier -> tail extraction -> polars DataFrame
      -> add_locus_tags() -> GFF3 / .tbl / .gbk

Four things are deliberate:

**The null is the filter, not a unit count.**  ``GGYNN`` carries two obligate
G's, so the naive false-positive rate climbs steeply with GC content - a 70%
GC actinobacteriophage will hand you arrays all day.  Each contig is therefore
scored against shuffles of itself, and the knob the user turns is
``--vrrna_max_p``, not a minimum number of units.  The default 0.05 with 100
shuffles is calibrated-to-conservative at high GC; 0.01 roughly halves the
false-call rate for ~1.9 points of recall and is where the remaining GC
artefacts drop out.

**The tail is reported, never gated on.**  The 3' tail is type-specific.  It
is strong evidence when it matches and says nothing at all when it does not,
and gating on it would discard the one locus in the benchmark that is
independently confirmed.  So ``tail_identity`` (against the type II consensus)
and ``tail_agreement`` (against the other arrays in the same contig) are
columns, not filters.

**Clustering is architecture, not a compounded probability.**  Three arrays
within 5 kb is a stronger VIPR locus than one, and the coordinates are free.
But all of a contig's arrays are scored against the same null draw, so their
p-values are correlated - three arrays at p = 0.01 are emphatically not
p = 1e-6.  ``cluster_size`` is reported as a description of the locus and
nothing is multiplied.

**The second tier lowers the unit floor, not the statistical bar.**  Once an
array has cleared the null, a shorter neighbour within ~2 kb of it is likely
real - the SUSP1 locus packs four arrays into 1.6 kb.  But a 4-unit array can
never clear a *genome-wide* null: 96% of shuffles of the SUSP1 genome contain a
4-unit run somewhere, which is the right answer to the wrong question.  Relaxed
calls are therefore scored against a null drawn from their own neighbourhood,
at the same ``--vrrna_max_p``.  Measured, that is what separates the two cases:
on Medea1 the three 4-5 unit neighbours of the real array score p = 0.29-0.53
locally and are dropped, while in the SUSP1 neighbourhood 85% of shuffles
produce nothing at all, so a genuine 5-unit satellite would be called at
p = 0.05.  The relaxed tier is also strictly additive - it can never change,
suppress or re-cut a strict call, so the strict output is identical to the
caller the benchmark was measured on.
"""

import os
import random

import numpy as np
import polars as pl
from loguru import logger

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Consensus vrRNA tail, derived from the SUSP1 vrRNA tails (NC_028808.2).
# Type II specific - see the module docstring on why it is never a gate.
DEFAULT_TAIL = "GGGTAGTTAGTCTCAGCTACATTAAAATGAGAC"

# Periodicity-prefilter density.  0.45 rather than the 0.55 of the first
# implementation: a 6-8 unit array is only 30-40 nt and cannot fill a 40 nt
# window at the higher density.  Worth 1 point of recall on native arrays and
# roughly double on solo ones, at a *lower* discordant rate.
DEFAULT_DENSITY = 0.45

# Bases two calls must share before one suppresses the other, in resolve().
MIN_OVERLAP = 5

# Strict, genome-wide pass.  6 rather than the published 9: every array of 9
# units or more was already found, and the whole recall deficit sat at 6-8.
DEFAULT_MIN_UNITS = 6
DEFAULT_MAX_UNITS = 40
DEFAULT_MIN_CONFORMITY = 0.80

# Relaxed, second-tier pass: only inside RELAXED_FLANK_BP of an array that has
# already cleared the null, or of a Vipr protein hit.
RELAXED_MIN_UNITS = 4
RELAXED_FLANK_BP = 2000

# Arrays this close together on a contig are called one cluster.
CLUSTER_WINDOW_BP = 5000

# Defaults for the per-contig shuffled null.
DEFAULT_SHUFFLES = 100
DEFAULT_MAX_P = 0.05
NULL_SEED = 42

VRRNA_TSV_COLUMNS = [
    "contig",
    "locus_tag",
    "start",
    "stop",
    "strand",
    "units",
    "conformity",
    "score",
    "pvalue",
    "tier",
    "cluster_id",
    "cluster_size",
    "tail_identity",
    "tail_agreement",
    "gc",
    "guides",
    "tail",
]

_SCHEMA = {
    "contig": pl.Utf8,
    "locus_tag": pl.Utf8,
    "start": pl.Int64,
    "stop": pl.Int64,
    "strand": pl.Utf8,
    "units": pl.Int64,
    "conformity": pl.Float64,
    "score": pl.Float64,
    "pvalue": pl.Float64,
    "tier": pl.Utf8,
    "cluster_id": pl.Int64,
    "cluster_size": pl.Int64,
    "tail_identity": pl.Float64,
    "tail_agreement": pl.Float64,
    "gc": pl.Float64,
    "guides": pl.Utf8,
    "tail": pl.Utf8,
}


def revcomp(seq):
    return seq.translate(COMP)[::-1]


# ---------------------------------------------------------------- units


def unit_ok(u):
    """All three fixed positions of a GGYNN unit match."""
    return len(u) == 5 and u[0] == "G" and u[1] == "G" and u[2] in "CT"


def fixed_matches(u):
    if len(u) != 5:
        return 0
    return (u[0] == "G") + (u[1] == "G") + (u[2] in "CT")


def unambiguous_nn(u):
    return len(u) == 5 and u[3] in "ACGT" and u[4] in "ACGT"


# ---------------------------------------------------- periodicity prefilter


def periodicity_windows(
    seq, period=5, window=40, density=DEFAULT_DENSITY, merge=300, pad=200
):
    """Candidate regions with period-5 self-similarity.

    Builds a match-at-lag track (``seq[i] == seq[i + period]``), convolves it
    over a window, and keeps runs above a density threshold.  Returns merged,
    padded ``(start, end)`` intervals.

    This is a ~5x speedup over scanning everything and it *lowers* the
    discordant-call rate, so it earns its place twice.  The density floor is
    what short arrays fail: a 6-8 unit array is 30-40 nt, so it cannot fill the
    40 nt window on its own and is diluted by its flanks.
    """
    n = len(seq)
    if n <= period + window:
        return [(0, n)]
    arr = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
    match = (arr[:-period] == arr[period:]).astype(np.float32)
    if len(match) < window:
        return [(0, n)]
    kern = np.ones(window, dtype=np.float32)
    dens = np.convolve(match, kern, mode="valid") / window
    hits = np.flatnonzero(dens >= density)
    if hits.size == 0:
        return []
    # group consecutive hit positions into runs, then merge and pad
    splits = np.flatnonzero(np.diff(hits) > merge)
    groups = np.split(hits, splits + 1)
    out = []
    for g in groups:
        s = max(0, int(g[0]) - pad)
        e = min(n, int(g[-1]) + window + pad)
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


# ------------------------------------------------------------- array calling


def unit_tracks(seq):
    """Per-position ``(fixed_matches, unambiguous_nn)`` for every 5-mer start.

    Exactly ``fixed_matches(seq[p:p + 5])`` and ``unambiguous_nn(seq[p:p + 5])``
    for p in ``range(len(seq) - 4)``, computed as two numpy passes instead of
    one Python slice and two function calls per position.

    This is a speedup, not a change: the per-unit predicates above remain the
    definition and ``tests/test_vrrna.py`` checks the two agree.  It matters
    because the shuffled null re-scans a whole genome 100 times over, and the
    string slicing was 80% of that cost - a 90 kb genome went from 11 s to
    well under 3 s when these two arrays stopped being rebuilt a character at
    a time.
    """
    arr = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
    n = arr.size
    if n < 5:
        return np.zeros(0, dtype=np.int8), np.zeros(0, dtype=bool)
    a, c, g, t = (ord(x) for x in "ACGT")
    fixed = (
        (arr[0 : n - 4] == g).astype(np.int8)
        + (arr[1 : n - 3] == g)
        + ((arr[2 : n - 2] == c) | (arr[2 : n - 2] == t))
    )
    acgt = (arr == a) | (arr == c) | (arr == g) | (arr == t)
    nn_ok = acgt[3 : n - 1] & acgt[4:n]
    return fixed, nn_ok


def call_in_window(
    seq, w_start, w_end, strand, min_units, max_units, min_conformity, tracks=None
):
    """Scan one window, all five register offsets, extend to the longest
    conforming block, then trim.

    Coordinates are 0-based half-open on ``seq`` as given, i.e. on the reverse
    complement when ``strand`` is ``-``; ``scan()`` flips them back.

    ``tracks`` is the ``unit_tracks(seq)`` pair.  It is computed here when not
    supplied; callers that scan the same sequence more than once - every one
    of them, since there are five register offsets - pass it in.
    """
    fixed_tr, nn_tr = unit_tracks(seq) if tracks is None else tracks
    ok_tr = fixed_tr == 3
    calls = []
    for offset in range(5):
        base = w_start + offset
        idx = range(base, w_end - 5, 5)
        n_idx = len(idx)
        if n_idx < min_units:
            continue
        sl = slice(base, base + 5 * n_idx, 5)
        fixed = fixed_tr[sl].tolist()
        oks = ok_tr[sl].tolist()
        nn_ok = nn_tr[sl]
        # Only a conforming unit can open an array, and in a 40% GC genome
        # barely 1% of units conform, so walking every index to find them is
        # 99% wasted work.  Jumping start to start over flatnonzero() is the
        # same traversal with the no-op steps removed - the body below is
        # reached for exactly the same i, in the same order.
        starts = np.flatnonzero(ok_tr[sl]).tolist()
        ptr = 0
        while ptr < len(starts):
            i = starts[ptr]
            # Extend to the LAST end whose whole block still clears the
            # conformity floor.  Breaking at the first prefix below it instead
            # cannot bridge a non-conforming unit early in a short array, even
            # when the block as a whole passes: GGCGC GCGGC GGCGC GGTCA TACCC
            # GGTAG GGCCG is 0.810 overall but 0.667 after two units.  Block
            # conformity is worth 7 points of boundary accuracy and 4 more
            # arrays over the prefix rule.
            hit = tot = 0
            end = i
            for j in range(i, min(i + max_units, n_idx)):
                hit += fixed[j]
                tot += 3
                if oks[j] and hit / tot >= min_conformity:
                    end = j
            # trim trailing non-conforming units: array must END on a good unit
            while end > i and not oks[end]:
                end -= 1
            n_units = end - i + 1
            if n_units >= min_units:
                conf = sum(fixed[i : end + 1]) / (3 * n_units)
                if conf >= min_conformity and nn_ok[i : end + 1].all():
                    block_oks = oks[i : end + 1]
                    block_seq = seq[idx[i] : idx[end] + 5]
                    nn = [block_seq[5 * k + 3 : 5 * k + 5] for k in range(n_units)]
                    p3 = {x[0] for x in nn}
                    p4 = {x[1] for x in nn}
                    run = best_run = 0
                    for is_ok in block_oks:
                        run = run + 1 if is_ok else 0
                        best_run = max(best_run, run)
                    if len(p3) >= 2 and len(p4) >= 2 and best_run >= 2:
                        s, e = idx[i], idx[end] + 5
                        calls.append(
                            {
                                "strand": strand,
                                "offset": offset,
                                "start": s,
                                "end": e,
                                "units": n_units,
                                "conformity": round(conf, 3),
                                "score": round(n_units * best_run * conf, 1),
                                "guide": "".join(nn),
                                "tail": seq[e : e + 60],
                            }
                        )
                resume = end + 1
            else:
                resume = i + 1
            while ptr < len(starts) and starts[ptr] < resume:
                ptr += 1
    return calls


def tail_identity(tail, consensus, min_overlap=25):
    """Best fractional identity between tail and consensus over any offset."""
    if len(tail) < min_overlap or not consensus:
        return 0.0
    best = 0.0
    for off in range(-len(consensus) + min_overlap, len(tail) - min_overlap + 1):
        a_s, b_s = max(0, off), max(0, -off)
        ln = min(len(tail) - a_s, len(consensus) - b_s)
        if ln < min_overlap:
            continue
        m = sum(1 for k in range(ln) if tail[a_s + k] == consensus[b_s + k])
        best = max(best, m / ln)
    return round(best, 3)


def resolve(calls, min_overlap=MIN_OVERLAP, seed=None):
    """More units wins, then higher conformity, then smaller span.

    Suppression needs ``min_overlap`` bases in common, not merely any overlap:
    adjacent tandem arrays in different registers abut with exactly 1 bp of
    shared span, so an any-overlap rule deletes one of every such pair.  The
    SUSP1 locus itself has two.

    ``seed`` pre-populates the kept list, so those calls survive unconditionally
    and everything in ``calls`` is judged against them.  That is how the
    relaxed tier is kept strictly additive to the strict one.
    """
    calls = sorted(
        calls, key=lambda c: (-c["units"], -c["conformity"], c["end"] - c["start"])
    )
    kept = list(seed) if seed else []
    for c in calls:
        if any(
            c["strand"] == k["strand"]
            and min(c["end"], k["end"]) - max(c["start"], k["start"]) >= min_overlap
            for k in kept
        ):
            continue
        kept.append(c)
    return kept


def scan(
    seq,
    min_units=DEFAULT_MIN_UNITS,
    max_units=DEFAULT_MAX_UNITS,
    min_conformity=DEFAULT_MIN_CONFORMITY,
    consensus=DEFAULT_TAIL,
    use_prefilter=True,
    density=DEFAULT_DENSITY,
):
    """Call arrays on both strands of one contig.

    Returns a list of dicts sorted by forward-strand start, each carrying
    ``fstart``/``fend`` as 0-based half-open forward-strand coordinates.
    """
    length = len(seq)
    out = []
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        windows = (
            periodicity_windows(s, density=density) if use_prefilter else [(0, len(s))]
        )
        calls = []
        tracks = unit_tracks(s)
        for w_start, w_end in windows:
            calls += call_in_window(
                s,
                w_start,
                w_end,
                strand,
                min_units,
                max_units,
                min_conformity,
                tracks=tracks,
            )
        for c in calls:
            c["tail_id"] = tail_identity(c["tail"], consensus)
            if strand == "-":
                c["fstart"], c["fend"] = length - c["end"], length - c["start"]
            else:
                c["fstart"], c["fend"] = c["start"], c["end"]
        out += calls
    return sorted(resolve(out), key=lambda c: c["fstart"])


# --------------------------------------------------------------- null models


def _mono_shuffler(seq, rng):
    """Mononucleotide shuffle: preserves base composition only."""
    chars = list(seq)
    while True:
        rng.shuffle(chars)
        yield "".join(chars)


def _last_edges_reach_root(lastedge, verts, root):
    """Do the chosen last-edges form a spanning in-tree rooted at ``root``?"""
    for v in verts:
        seen = {v}
        u = v
        while u != root:
            u = lastedge.get(u)
            if u is None or u in seen:
                return False
            seen.add(u)
    return True


def dinuc_shuffle(seq, rng):
    """One Altschul-Erikson dinucleotide-preserving shuffle.

    Altschul & Erikson 1985, Mol. Biol. Evol. 2:526.  A sequence is an Eulerian
    path through the graph whose vertices are bases and whose edges are the
    dinucleotides, so any other Eulerian path over the same edge multiset has
    *identical* mononucleotide and dinucleotide counts.  A uniform random such
    path is drawn by choosing, for every vertex except the sequence's last
    base, one outgoing edge to traverse last; those edges must form a spanning
    in-tree rooted at the last base (otherwise the walk strands itself in a
    cycle), and the rest of each vertex's edges are permuted freely.

    First and last base are preserved, which is what makes the edge multiset of
    the result equal to that of the input.
    """
    n = len(seq)
    if n < 3:
        return seq
    succ = {}
    for i in range(n - 1):
        a = seq[i]
        if a in succ:
            succ[a].append(seq[i + 1])
        else:
            succ[a] = [seq[i + 1]]
    last = seq[-1]
    need = [v for v in succ if v != last]
    while True:
        lastedge = {v: succ[v][rng.randrange(len(succ[v]))] for v in need}
        if _last_edges_reach_root(lastedge, need, last):
            break
    order = {}
    for v, edges in succ.items():
        e = list(edges)
        if v in lastedge:
            e.remove(lastedge[v])  # same-character edges are interchangeable
            rng.shuffle(e)
            e.append(lastedge[v])
        else:
            rng.shuffle(e)
        order[v] = e
    cur = seq[0]
    out = [cur]
    ptr = dict.fromkeys(succ, 0)
    push = out.append
    for _ in range(n - 1):
        k = ptr[cur]
        ptr[cur] = k + 1
        cur = order[cur][k]
        push(cur)
    return "".join(out)


def _dinuc_shuffler(seq, rng):
    while True:
        yield dinuc_shuffle(seq, rng)


NULL_MODELS = {"mono": _mono_shuffler, "dinuc": _dinuc_shuffler}
DEFAULT_NULL_MODEL = "mono"


def null_max_units(
    seq,
    n_shuffles,
    min_units,
    max_units=DEFAULT_MAX_UNITS,
    min_conformity=DEFAULT_MIN_CONFORMITY,
    seed=NULL_SEED,
    model=DEFAULT_NULL_MODEL,
):
    """Distribution of the best array size under shuffling.

    GGYNN carries two obligate G's, so the false-positive rate rises steeply
    with GC content.  Shuffling preserves base composition, which makes this an
    empirical, per-contig control for that bias.

    ``model`` selects what the shuffle holds fixed:

    ``mono``
        base composition only (the default)
    ``dinuc``
        base *and* dinucleotide composition (Altschul-Erikson)

    Measured, ``dinuc`` is the *weaker* null: GC-rich phage DNA is
    GG-*depleted* relative to its own base composition (f(GG)/f(G)^2 ~ 0.84
    above 60% GC), so holding dinucleotides fixed makes ~15% fewer GGY
    trinucleotides than a mononucleotide shuffle and lowers the threshold by
    ~0.5 units.  It is also anti-conservative at alpha <= 0.01, which is the
    regime that matters.  ``mono`` is the default for those reasons, and
    ``dinuc`` is offered as the control it is.

    The null is drawn without the prefilter - it has to see everything a
    shuffle can produce, and shuffled sequence has no periodicity for the
    prefilter to find.
    """
    rng = random.Random(seed)
    try:
        gen = NULL_MODELS[model](seq, rng)
    except KeyError:
        raise ValueError(
            f"unknown null model {model!r}; choose from {sorted(NULL_MODELS)}"
        )
    out = []
    for _ in range(n_shuffles):
        sh = next(gen)
        cs = []
        for strand, s in (("+", sh), ("-", revcomp(sh))):
            cs += call_in_window(
                s,
                0,
                len(s),
                strand,
                min_units,
                max_units,
                min_conformity,
                tracks=unit_tracks(s),
            )
        out.append(max((c["units"] for c in cs), default=0))
    return out


def empirical_p(units, nulls):
    """One-sided empirical p with the usual +1 correction.

    ``nulls`` is the per-shuffle maximum array size over whatever sequence was
    shuffled, so this is the fraction of shuffles that produced an array at
    least this long *anywhere in that sequence*.  For the strict tier that
    sequence is the whole contig, which is why the arrays of one contig have
    correlated p-values and must not be multiplied together.
    """
    n = len(nulls)
    if n == 0:
        return float("nan")
    ge = sum(1 for m in nulls if m >= units)
    return round((ge + 1) / (n + 1), 4)


# ------------------------------------------------------- the two-tier caller


def _relaxed_regions(length, anchors, flank=RELAXED_FLANK_BP):
    """Merged ``flank``-padded intervals around ``anchors``.

    ``anchors`` are ``(start, end)`` pairs in forward-strand 0-based half-open
    coordinates - accepted arrays, or Vipr/Vap protein hits.
    """
    spans = sorted((max(0, s - flank), min(length, e + flank)) for s, e in anchors)
    merged = []
    for s, e in spans:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def _local_null(
    seq,
    region,
    n_shuffles,
    min_units,
    max_units=DEFAULT_MAX_UNITS,
    min_conformity=DEFAULT_MIN_CONFORMITY,
    seed=NULL_SEED,
    model=DEFAULT_NULL_MODEL,
):
    """The shuffled null for one relaxed region, drawn from that region only.

    The strict tier asks "could a shuffle of this contig produce an array this
    long anywhere?", which is the right question for a genome-wide search and
    the wrong one for a 4 kb window next to a confirmed array - 96% of
    shuffles of the SUSP1 genome contain a 4-unit run somewhere, so no short
    array could ever pass.  Scoping the shuffle to the window asks the
    question the relaxed tier is actually asking.

    The region includes the accepted array that anchored it, whose GGYNN
    content inflates the G count of the thing being shuffled.  That biases the
    null upwards, i.e. towards rejecting the relaxed call, which is the safe
    direction.
    """
    start, end = region
    return null_max_units(
        seq[start:end],
        n_shuffles,
        min_units,
        max_units,
        min_conformity,
        seed=seed,
        model=model,
    )


def _call_regions(seq, regions, min_units, max_units, min_conformity, consensus):
    """Run both strands over a set of forward-strand regions, no prefilter."""
    length = len(seq)
    rc = revcomp(seq)
    tracks = {"+": unit_tracks(seq), "-": unit_tracks(rc)}
    out = []
    for f_start, f_end in regions:
        for strand, s, w_start, w_end in (
            ("+", seq, f_start, f_end),
            ("-", rc, length - f_end, length - f_start),
        ):
            for c in call_in_window(
                s,
                w_start,
                w_end,
                strand,
                min_units,
                max_units,
                min_conformity,
                tracks=tracks[strand],
            ):
                c["tail_id"] = tail_identity(c["tail"], consensus)
                if strand == "-":
                    c["fstart"], c["fend"] = length - c["end"], length - c["start"]
                else:
                    c["fstart"], c["fend"] = c["start"], c["end"]
                out.append(c)
    return out


def scan_contig(
    seq,
    min_units=DEFAULT_MIN_UNITS,
    max_units=DEFAULT_MAX_UNITS,
    min_conformity=DEFAULT_MIN_CONFORMITY,
    consensus=DEFAULT_TAIL,
    use_prefilter=True,
    density=DEFAULT_DENSITY,
    n_shuffles=DEFAULT_SHUFFLES,
    max_p=DEFAULT_MAX_P,
    null_model=DEFAULT_NULL_MODEL,
    seed=NULL_SEED,
    relaxed_min_units=RELAXED_MIN_UNITS,
    relaxed_flank=RELAXED_FLANK_BP,
    protein_hits=None,
):
    """Call, score and accept vrRNA arrays on one contig.

    Four steps:

    1. strict genome-wide pass at ``min_units``;
    2. a whole-contig shuffled null, which gives every call an empirical p,
       and arrays with ``p > max_p`` are dropped;
    3. a relaxed pass at ``relaxed_min_units`` restricted to within
       ``relaxed_flank`` of a surviving array or of a ``protein_hits`` entry,
       each call scored against a null drawn from its own region and held to
       the same ``max_p`` - see the module docstring on why the null has to be
       local here;
    4. tail agreement and clustering over whatever survived.

    ``protein_hits`` is the hook for a Vipr/Vap HMM: a list of ``(start, stop)``
    1-based inclusive coordinates on this contig.  A locus with a Vipr protein
    but no array that clears the null genome-wide is exactly the case the
    relaxed tier exists for.  Nothing calls it yet; the HMM is not built.

    :return: list of call dicts, sorted by start, with ``pvalue`` and ``tier``
    """
    length = len(seq)
    strict = scan(
        seq, min_units, max_units, min_conformity, consensus, use_prefilter, density
    )
    for c in strict:
        c["tier"] = "strict"

    anchors = []
    if protein_hits:
        # 1-based inclusive in, 0-based half-open out
        anchors += [(max(0, int(s) - 1), min(length, int(e))) for s, e in protein_hits]

    nulls = []
    if n_shuffles > 0 and strict:
        nulls = null_max_units(
            seq,
            n_shuffles,
            min_units,
            max_units,
            min_conformity,
            seed=seed,
            model=null_model,
        )

    for c in strict:
        c["pvalue"] = empirical_p(c["units"], nulls) if nulls else float("nan")
    accepted = [c for c in strict if not nulls or c["pvalue"] <= max_p]

    anchors += [(c["fstart"], c["fend"]) for c in accepted]
    if anchors and relaxed_min_units < min_units:
        regions = _relaxed_regions(length, anchors, relaxed_flank)
        relaxed = []
        for region in regions:
            local = (
                _local_null(
                    seq,
                    region,
                    n_shuffles,
                    relaxed_min_units,
                    max_units,
                    min_conformity,
                    seed=seed,
                    model=null_model,
                )
                if n_shuffles > 0
                else []
            )
            for c in _call_regions(
                seq, [region], relaxed_min_units, max_units, min_conformity, consensus
            ):
                c["tier"] = "relaxed"
                c["pvalue"] = empirical_p(c["units"], local) if local else float("nan")
                if not local or c["pvalue"] <= max_p:
                    relaxed.append(c)
        # Seeded with the accepted calls, so a relaxed re-discovery of one is
        # suppressed rather than duplicating or re-cutting it: the strict tier
        # comes out of scan_contig() exactly as scan() produced it.
        accepted = resolve(relaxed, seed=accepted)

    accepted.sort(key=lambda c: (c["fstart"], c["fend"]))
    _annotate_tail_agreement(accepted)
    _annotate_clusters(accepted)
    return accepted


def _annotate_tail_agreement(calls):
    """Best tail identity of each array to any *other* array on the contig.

    The tail is type-specific, so a shared tail between two arrays in one
    genome is good evidence they belong to one VIPR system even when neither
    matches the type II consensus.  With a single array there is nothing to
    compare against and the field is null.
    """
    for c in calls:
        others = [o["tail"] for o in calls if o is not c]
        c["tail_agreement"] = (
            max((tail_identity(c["tail"], o) for o in others), default=None)
            if others
            else None
        )


def _annotate_clusters(calls, window=CLUSTER_WINDOW_BP):
    """Single-linkage grouping of arrays within ``window`` bp on one contig.

    Reported, never multiplied: every array on a contig is scored against the
    same null draw, so a three-array cluster at p = 0.01 is three correlated
    observations of one locus and not p = 1e-6.  What it is good for is
    architecture - a cluster looks like a VIPR locus, a lone marginal array
    looks like the null's own false-positive rate.
    """
    cluster_id = 0
    members = []
    last_end = None
    for c in calls:
        if last_end is None or c["fstart"] - last_end > window:
            cluster_id += 1
            members.append([])
        members[-1].append(c)
        c["cluster_id"] = cluster_id
        last_end = max(last_end or 0, c["fend"])
    for group in members:
        for c in group:
            c["cluster_size"] = len(group)


# ----------------------------------------------------------------- dataframe


def empty_vrrna_df():
    """An empty dataframe with the vrRNA schema, so downstream code is uniform."""
    return pl.DataFrame(
        {name: pl.Series([], dtype=dtype) for name, dtype in _SCHEMA.items()}
    )


def _gc(seq):
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def calls_to_records(contig, seq, calls):
    """Turn one contig's calls into TSV rows.

    Coordinates become 1-based inclusive here, which is what the rest of
    pharokka - and GFF3, the .tbl and GenBank - all use.  The caller itself
    works in 0-based half-open coordinates, matching the published array
    table.
    """
    gc = round(_gc(seq), 4)
    records = []
    for c in calls:
        records.append(
            {
                "contig": contig,
                "locus_tag": "",
                "start": c["fstart"] + 1,
                "stop": c["fend"],
                "strand": c["strand"],
                "units": c["units"],
                "conformity": float(c["conformity"]),
                "score": float(c["score"]),
                "pvalue": float(c["pvalue"]),
                "tier": c["tier"],
                "cluster_id": c["cluster_id"],
                "cluster_size": c["cluster_size"],
                "tail_identity": float(c["tail_id"]),
                "tail_agreement": c["tail_agreement"],
                "gc": gc,
                "guides": c["guide"],
                "tail": c["tail"],
            }
        )
    return records


def read_fasta(path):
    """Minimal FASTA reader: yields ``(id, uppercased sequence)``."""
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield name.split()[0], "".join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line.upper())
    if name:
        yield name.split()[0], "".join(seq)


def scan_fasta(
    fasta_path,
    max_p=DEFAULT_MAX_P,
    n_shuffles=DEFAULT_SHUFFLES,
    null_model=DEFAULT_NULL_MODEL,
    min_units=DEFAULT_MIN_UNITS,
    min_conformity=DEFAULT_MIN_CONFORMITY,
    seed=NULL_SEED,
    protein_hits=None,
):
    """Scan every contig of a FASTA and return the vrRNA dataframe.

    :param protein_hits: optional ``{contig: [(start, stop), ...]}`` of Vipr /
        Vap protein hits, 1-based inclusive, seeding the relaxed tier.  The
        Vipr HMM that would fill this does not exist yet.
    :return: polars dataframe with VRRNA_TSV_COLUMNS, ``locus_tag`` left empty
    """
    protein_hits = protein_hits or {}
    records = []
    for contig, seq in read_fasta(fasta_path):
        if not seq:
            continue
        calls = scan_contig(
            seq,
            min_units=min_units,
            min_conformity=min_conformity,
            n_shuffles=n_shuffles,
            max_p=max_p,
            null_model=null_model,
            seed=seed,
            protein_hits=protein_hits.get(contig),
        )
        records += calls_to_records(contig, seq, calls)

    if not records:
        return empty_vrrna_df()

    return pl.DataFrame(records, schema=_SCHEMA).sort(["contig", "start"])


def add_locus_tags(df, locustag, contig_count):
    """Assigns pharokka locus tags to vrRNA rows.

    Follows the tRNA and ncRNA convention in post_processing.create_gff():
    per-contig numbering when there are multiple contigs, otherwise a single
    run of numbers prefixed with the locustag.
    """
    if df.height == 0:
        return df

    if contig_count > 1:
        return df.with_columns(
            (
                pl.col("contig")
                + pl.lit("_vrRNA_")
                + (pl.col("contig").cum_count().over("contig")).cast(pl.Utf8)
            ).alias("locus_tag")
        )

    return df.with_columns(
        (
            pl.lit(f"{locustag}_vrRNA_") + (pl.int_range(1, pl.len() + 1)).cast(pl.Utf8)
        ).alias("locus_tag")
    )


def vrrna_note(row):
    """The free-text ``note`` for one array, for GFF3, the .tbl and GenBank.

    No ``;``, ``=`` or ``,``: all three are structural in GFF3 attributes, and
    a note containing one silently corrupts the record.
    """
    bits = [
        f"VIPR vrRNA guide array ({row['tier']} call)",
        f"{row['units']} GGYNN units",
        f"fixed-position conformity {row['conformity']:.3f}",
        f"empirical p {row['pvalue']:.4f} vs shuffled contig",
    ]
    if row["cluster_size"] and row["cluster_size"] > 1:
        bits.append(
            f"member of vrRNA cluster {row['cluster_id']} on this contig - "
            f"{row['cluster_size']} arrays within {CLUSTER_WINDOW_BP} bp "
            "(correlated evidence as they share one null draw)"
        )
    bits.append(f"3' tail identity to type II consensus {row['tail_identity']:.3f}")
    if row["tail_agreement"] is not None:
        bits.append(
            f"best tail identity to another array on this contig "
            f"{row['tail_agreement']:.3f}"
        )
    bits.append(f"guide dinucleotides {row['guides']}")
    return " | ".join(bits)


def write_vrrna_tsv(df, out_dir, prefix):
    """Writes the ``{prefix}_vrrna.tsv`` output."""
    out_path = os.path.join(out_dir, f"{prefix}_vrrna.tsv")
    df.select(VRRNA_TSV_COLUMNS).write_csv(
        out_path, separator="\t", quote_style="never"
    )
    logger.info(f"{df.height} vrRNA array(s) written to {out_path}")
