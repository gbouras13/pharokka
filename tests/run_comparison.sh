#!/usr/bin/env bash
# run_comparison.sh
# Runs every output-producing test case from test_overall.py with both:
#   - bioconda 1.9.1  (pharokka.py)
#   - refactored dev  (pharokka run)
# then compares key output files ignoring timestamps.
#
# Usage (from repo root): bash tests/run_comparison.sh
# Outputs go to /tmp/pharokka_compare/

set -uo pipefail

DB="${PHAROKKA_DB:-tests/test_data/database}"
T=8
BASE="/tmp/pharokka_compare"
D="tests/test_data"
CUSTOM_HMM="${D}/custom_db/microvirus.h3m"

# Set PHAROKKA_NEW_PATH to the bin/ of your refactored pharokka conda env
# Set PHAROKKA_REF_PATH to the bin/ of the pharokka 1.9.1 bioconda env
# e.g.:
#   export PHAROKKA_NEW_PATH="$HOME/miniforge3/envs/pharokka_env/bin"
#   export PHAROKKA_REF_PATH="$HOME/miniforge3/envs/pharokka_1.9.1/bin"
NEW_PATH="${PHAROKKA_NEW_PATH:?Set PHAROKKA_NEW_PATH to the refactored pharokka env bin/ path}"
REF_PATH="${PHAROKKA_REF_PATH:?Set PHAROKKA_REF_PATH to the pharokka 1.9.1 env bin/ path}"

mkdir -p "$BASE"

PASS=0
FAIL=0
ERRORS=()
TP=4   # threads for proteins tests

run_case() {
    local name="$1"; shift
    local args="$*"

    local out_new="$BASE/${name}_new"
    local out_ref="$BASE/${name}_ref"

    printf "\n══ %-45s ══\n" "$name"

    # new (refactored): pharokka run ...
    printf "  [%s] → new (pharokka run) starting...\n" "$(date '+%H:%M:%S')"
    if eval "PATH='$NEW_PATH:$PATH' pharokka run $args -d '$DB' -o '$out_new' -t $T -f" \
            > "$BASE/${name}_new.log" 2>&1; then
        printf "  [%s] → new DONE ✓\n" "$(date '+%H:%M:%S')"
    else
        printf "  [%s] → new FAILED ✗ (see %s)\n" "$(date '+%H:%M:%S')" "$BASE/${name}_new.log"
        ERRORS+=("$name: new run FAILED — see $BASE/${name}_new.log")
        FAIL=$((FAIL+1))
        return
    fi

    # ref (bioconda 1.9.1): pharokka.py ...
    printf "  [%s] → ref (pharokka.py 1.9.1) starting...\n" "$(date '+%H:%M:%S')"
    if eval "PATH='$REF_PATH:$PATH' pharokka.py $args -d '$DB' -o '$out_ref' -t $T -f" \
            > "$BASE/${name}_ref.log" 2>&1; then
        printf "  [%s] → ref DONE ✓\n" "$(date '+%H:%M:%S')"
    else
        printf "  [%s] → ref FAILED ✗ (see %s)\n" "$(date '+%H:%M:%S')" "$BASE/${name}_ref.log"
        ERRORS+=("$name: ref run FAILED — see $BASE/${name}_ref.log")
        FAIL=$((FAIL+1))
        return
    fi

    # compare
    printf "  [%s] → diff starting...\n" "$(date '+%H:%M:%S')"
    if PATH="$NEW_PATH:$PATH" python tests/compare_outputs.py "$out_new" "$out_ref"; then
        printf "  [%s] → diff MATCH ✓\n" "$(date '+%H:%M:%S')"
        PASS=$((PASS+1))
    else
        printf "  [%s] → diff MISMATCH ✗\n" "$(date '+%H:%M:%S')"
        ERRORS+=("$name: OUTPUT MISMATCH")
        FAIL=$((FAIL+1))
    fi
}

run_case_proteins() {
    local name="$1"; shift
    local args="$*"

    local out_new="$BASE/${name}_new"
    local out_ref="$BASE/${name}_ref"

    printf "\n══ %-45s ══\n" "$name"

    # new (refactored): pharokka proteins ...
    printf "  [%s] → new (pharokka proteins) starting...\n" "$(date '+%H:%M:%S')"
    if eval "PATH='$NEW_PATH:$PATH' pharokka proteins $args -d '$DB' -o '$out_new' -t $TP -f" \
            > "$BASE/${name}_new.log" 2>&1; then
        printf "  [%s] → new DONE ✓\n" "$(date '+%H:%M:%S')"
    else
        printf "  [%s] → new FAILED ✗ (see %s)\n" "$(date '+%H:%M:%S')" "$BASE/${name}_new.log"
        ERRORS+=("$name: new run FAILED — see $BASE/${name}_new.log")
        FAIL=$((FAIL+1))
        return
    fi

    # ref (bioconda 1.9.1): pharokka_proteins.py ...
    printf "  [%s] → ref (pharokka_proteins.py 1.9.1) starting...\n" "$(date '+%H:%M:%S')"
    if eval "PATH='$REF_PATH:$PATH' pharokka_proteins.py $args -d '$DB' -o '$out_ref' -t $TP -f" \
            > "$BASE/${name}_ref.log" 2>&1; then
        printf "  [%s] → ref DONE ✓\n" "$(date '+%H:%M:%S')"
    else
        printf "  [%s] → ref FAILED ✗ (see %s)\n" "$(date '+%H:%M:%S')" "$BASE/${name}_ref.log"
        ERRORS+=("$name: ref run FAILED — see $BASE/${name}_ref.log")
        FAIL=$((FAIL+1))
        return
    fi

    # compare
    printf "  [%s] → diff starting...\n" "$(date '+%H:%M:%S')"
    if PATH="$NEW_PATH:$PATH" python tests/compare_outputs.py "$out_new" "$out_ref"; then
        printf "  [%s] → diff MATCH ✓\n" "$(date '+%H:%M:%S')"
        PASS=$((PASS+1))
    else
        printf "  [%s] → diff MISMATCH ✗\n" "$(date '+%H:%M:%S')"
        ERRORS+=("$name: OUTPUT MISMATCH")
        FAIL=$((FAIL+1))
    fi
}

# ── all test cases from test_overall.py (output-producing only) ─────────────

run_case "overall" \
    "-i ${D}/overall/Standard_examples/SAOMS1.fasta -l PHARTEST"

run_case "overall_reverse_mmseqs_sensitivity" \
    "-i ${D}/overall/Standard_examples/SAOMS1.fasta --reverse_mmseqs2 -g pyrodigal-rv --sensitivity 0.5 -l PHARTEST"

run_case "overall_trna_anticodon" \
    "-i ${D}/overall/bug_examples/AJ251789_trna_anticodon.fa --fast -l PHARTEST"

run_case "overall_pyhmmer_alphabet" \
    "-i ${D}/overall/bug_examples/pyhmmer_alphabet_issue357.fasta --fast -l PHARTEST"

run_case "overall_vfdb_issue_410" \
    "-i ${D}/overall/bug_examples/vfdb_issue_410.fasta --fast -l PHARTEST"

run_case "overall_mash_distance" \
    "-i ${D}/overall/Standard_examples/SAOMS1.fasta --mash_distance 0.05 -l PHARTEST"

run_case "overall_crispr" \
    "-i ${D}/overall/CRISPR_example/Biggiephage_A_fullcontig_CasΦ1.fasta -l PHARTEST"

run_case "overall_crispr_minced_args" \
    "-i ${D}/overall/CRISPR_example/Biggiephage_A_fullcontig_CasΦ1.fasta -g prodigal --minced_args 'minNR 2 -minRL 21' -l PHARTEST"

run_case "overall_vfdb" \
    "-i ${D}/overall/VFDB_example/NC_004617.fasta --skip_extra_annotations -l PHARTEST"

run_case "overall_amr" \
    "-i ${D}/overall/AMR_example/NC_007458.fasta --skip_mash -l PHARTEST"

run_case "overall_tmrna" \
    "-i ${D}/overall/tmRNA_example/NC_051700.fasta -l PHARTEST"

run_case "overall_numeric_header_prodigal" \
    "-i ${D}/overall/Standard_examples/SAOMS1_numeric_header.fasta -g prodigal --fast -m"

run_case "overall_numeric_header_phanotate" \
    "-i ${D}/overall/Standard_examples/SAOMS1_numeric_header.fasta --fast -l PHARTEST"

run_case "overall_numeric_header_prodigal_gv" \
    "-i ${D}/overall/Standard_examples/SAOMS1_numeric_header.fasta -m --fast"

run_case "meta" \
    "-i ${D}/overall/Meta_example/combined_meta.fasta -m"

run_case "meta_unicycler_header_prodigal_gv" \
    "-i ${D}/overall/Meta_example/combined_meta_unicycler_headers.fasta -g prodigal-gv -m"

run_case "pyrodigal_gv" \
    "-i ${D}/overall/Meta_example/combined_meta.fasta -m -g prodigal-gv"

run_case "prodigal_rv" \
    "-i ${D}/overall/rv_example/Tymoviridae.fna -g pyrodigal-rv -l PHARTEST"

run_case "meta_dnaapler_all_bug" \
    "-i ${D}/overall/Meta_example/combined_meta.fasta -m -s --dnaapler --meta_hmm"

run_case "overall_locus" \
    "-i ${D}/overall/Standard_examples/SAOMS1.fasta -l SAOMS1 -p SAOMS1"

run_case "custom" \
    "-i ${D}/overall/custom_examples/MH649026.fasta --custom_hmm ${CUSTOM_HMM} -l PHARTEST"

run_case "custom_meta" \
    "-i ${D}/overall/custom_examples/hundred_microviruses.fasta --custom_hmm ${CUSTOM_HMM} -m"

run_case "overall_prodigal" \
    "-i ${D}/overall/Standard_examples/SAOMS1.fasta -g prodigal -l PHARTEST"

run_case "overall_stop_recode" \
    "-i ${D}/overall/stop_recoding/table_4/SRR1747055_scaffold_7.fa -g prodigal -c 4 -l PHARTEST"

run_case "overall_dnaapler" \
    "-i ${D}/overall/Standard_examples/SAOMS1.fasta --dnaapler -l PHARTEST"

run_case "overall_fast" \
    "-i ${D}/overall/Standard_examples/SAOMS1.fasta --fast -l PHARTEST"

run_case "overall_mmseqs" \
    "-i ${D}/overall/Standard_examples/SAOMS1.fasta --mmseqs2_only -l PHARTEST"

run_case "meta_no_cds_contig" \
    "-i ${D}/overall/Meta_example/fake_meta.fa -m"

run_case "terminase" \
    "-i ${D}/overall/Standard_examples/SAOMS1.fasta --terminase --terminase_start 340 --terminase_strand neg -l PHARTEST"

run_case "overall_genbank" \
    "-i ${D}/overall/genbank_examples/SAOMS1.gbk --genbank -l PHARTEST"

run_case "overall_genbank_meta" \
    "-i ${D}/overall/genbank_examples/hundred_microviruses.gbk --genbank -m --meta_hmm"

# ── all test cases from test_proteins.py (output-producing only) ────────────

run_case_proteins "proteins" \
    "-i ${D}/proteins/phanotate.faa"

run_case_proteins "proteins_with_vfdb_card" \
    "-i ${D}/proteins/vfdb_card.faa"

run_case_proteins "proteins_hmm_only" \
    "-i ${D}/proteins/phanotate.faa --hmm_only"

run_case_proteins "proteins_mmseqs_only" \
    "-i ${D}/proteins/phanotate.faa --mmseqs2_only"

run_case_proteins "proteins_reverse_mmseqs_sensitivity" \
    "-i ${D}/proteins/phanotate.faa --reverse_mmseqs2 --sensitivity 0.5"

# ── summary ──────────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════════════"
printf "  Results: %d passed, %d failed\n" "$PASS" "$FAIL"
echo "══════════════════════════════════════════════════════"
if [ ${#ERRORS[@]} -gt 0 ]; then
    echo "  Failures:"
    for e in "${ERRORS[@]}"; do
        echo "    ✗ $e"
    done
    exit 1
else
    echo "  All cases match."
    exit 0
fi
