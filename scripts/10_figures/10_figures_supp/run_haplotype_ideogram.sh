#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# run_haplotype_ideogram.sh
#
# Generates publication-ready and poster-ready butterfly haplotype ideograms
# with HiFi Single, HiFi Dual, and ONT Dual assemblies compared against the
# T2T reference.
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(cd "$SCRIPT_DIR/../.." && pwd)"

PYTHON="/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/ds/bin/python3"
if [[ ! -x "$PYTHON" ]]; then
    PYTHON="python3"
fi

PLOT_SCRIPT="$SCRIPT_DIR/plot_ideogram_haplotype_avgcov_v6.py"

# Data & Reference
CEN="$BASE/data/t2t/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff"
TELO_P="$BASE/data/t2t/bTaeGut7.T2T.fasta_terminal_telomeres.bed"

# HiFi Single inputs
S_DIR="$BASE/results/single"
S_TSV="$S_DIR/03_coverage_single/chain_cov.target_cov.tsv"
S_COLL_CHAIN="$S_DIR/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.target.collinear.chain"
S_NC_CHAIN="$S_DIR/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.target.non-collinear.chain"
S_UNLOC_CHAIN="$S_DIR/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.unloc.rbest.chain"
S_PAIRS="$S_DIR/02_chainpipeline_single/t2t.vs.single.best_chrom_pairs.tsv"
S_BED="$S_DIR/05_hapmers_single/single.switch_blocks.final.bed"
S_TELO="$S_DIR/04_telomeres_single/single_telomere_presence.tsv"
S_COV_SUM="$S_DIR/03_coverage_single/single.coverage_summary.tsv"
S_GAPS="$S_DIR/06_gaps_single/single_combined.renamed.sorted.reoriented.annotated.gaps.bed"

# HiFi Dual inputs
D_DIR="$BASE/results/dual"
D_TSV="$D_DIR/03_coverage_dual/chain_cov.target_cov.tsv"
D_COLL_CHAIN="$D_DIR/02_chainpipeline_dual/t2t.vs.dual.T2T.vs.ASM.target.collinear.chain"
D_NC_CHAIN="$D_DIR/02_chainpipeline_dual/t2t.vs.dual.T2T.vs.ASM.target.non-collinear.chain"
D_UNLOC_CHAIN="$D_DIR/02_chainpipeline_dual/t2t.vs.dual.T2T.vs.ASM.unloc.rbest.chain"
D_PAIRS="$D_DIR/02_chainpipeline_dual/t2t.vs.dual.best_chrom_pairs.tsv"
D_BED="$D_DIR/05_hapmers_dual/dual.switch_blocks.final.bed"
D_TELO="$D_DIR/04_telomeres_dual/dual_telomere_presence.tsv"
D_COV_SUM="$D_DIR/03_coverage_dual/dual.coverage_summary.tsv"
D_GAPS="$D_DIR/06_gaps_dual/dual_combined.renamed.sorted.reoriented.annotated.gaps.bed"

# ONT Dual inputs
O_DIR="$BASE/results/ont"
O_COLL_CHAIN="$O_DIR/02_chainpipeline_ont/t2t.vs.ont.T2T.vs.ASM.target.collinear.chain"
O_NC_CHAIN="$O_DIR/02_chainpipeline_ont/t2t.vs.ont.T2T.vs.ASM.target.non-collinear.chain"
O_PAIRS="$O_DIR/02_chainpipeline_ont/t2t.vs.ont.best_chrom_pairs.tsv"
O_BED="$O_DIR/05_hapmers_ont/ont.switch_blocks.final.bed"
O_TELO="$O_DIR/04_telomeres_ont/ont_telomere_presence.tsv"
O_COV_SUM="$O_DIR/03_coverage_ont/ont.coverage_summary.tsv"
O_GAPS="$O_DIR/06_gaps_ont/asm3_ONT_combined.sorted.reoriented.annotated.gaps.bed"

# Output directory
OUT_DIR="$BASE/results/figures"
mkdir -p "$OUT_DIR"

echo "=== Generating Paper Butterfly Ideogram (PNG, 300 DPI) ==="
"$PYTHON" "$PLOT_SCRIPT" \
    --single-tsv "$S_TSV" \
    --dual-tsv "$D_TSV" \
    --single-chain "$S_COLL_CHAIN" \
    --single-nc-chain "$S_NC_CHAIN" \
    --single-unlocs-chain "$S_UNLOC_CHAIN" \
    --single-pairs "$S_PAIRS" \
    --single-bed "$S_BED" \
    --single-telomere-presence "$S_TELO" \
    --single-coverage-summary "$S_COV_SUM" \
    --single-annotated-gaps "$S_GAPS" \
    --dual-chain "$D_COLL_CHAIN" \
    --dual-nc-chain "$D_NC_CHAIN" \
    --dual-unlocs-chain "$D_UNLOC_CHAIN" \
    --dual-pairs "$D_PAIRS" \
    --dual-bed "$D_BED" \
    --dual-telomere-presence "$D_TELO" \
    --dual-coverage-summary "$D_COV_SUM" \
    --dual-annotated-gaps "$D_GAPS" \
    --ont-dual-chain "$O_COLL_CHAIN" \
    --ont-dual-nc-chain "$O_NC_CHAIN" \
    --ont-dual-pairs "$O_PAIRS" \
    --ont-dual-bed "$O_BED" \
    --ont-telomere-presence "$O_TELO" \
    --ont-coverage-summary "$O_COV_SUM" \
    --ont-annotated-gaps "$O_GAPS" \
    --centromeres "$CEN" \
    --telo-p-bed "$TELO_P" \
    --layout butterfly \
    --style paper \
    --format png \
    --dpi 300 \
    --output "$OUT_DIR/bTaeGut7_haplotype_ideogram_paper.png"

echo "=== Generating Paper Butterfly Ideogram (PDF) ==="
"$PYTHON" "$PLOT_SCRIPT" \
    --single-tsv "$S_TSV" \
    --dual-tsv "$D_TSV" \
    --single-chain "$S_COLL_CHAIN" \
    --single-nc-chain "$S_NC_CHAIN" \
    --single-unlocs-chain "$S_UNLOC_CHAIN" \
    --single-pairs "$S_PAIRS" \
    --single-bed "$S_BED" \
    --single-telomere-presence "$S_TELO" \
    --single-coverage-summary "$S_COV_SUM" \
    --single-annotated-gaps "$S_GAPS" \
    --dual-chain "$D_COLL_CHAIN" \
    --dual-nc-chain "$D_NC_CHAIN" \
    --dual-unlocs-chain "$D_UNLOC_CHAIN" \
    --dual-pairs "$D_PAIRS" \
    --dual-bed "$D_BED" \
    --dual-telomere-presence "$D_TELO" \
    --dual-coverage-summary "$D_COV_SUM" \
    --dual-annotated-gaps "$D_GAPS" \
    --ont-dual-chain "$O_COLL_CHAIN" \
    --ont-dual-nc-chain "$O_NC_CHAIN" \
    --ont-dual-pairs "$O_PAIRS" \
    --ont-dual-bed "$O_BED" \
    --ont-telomere-presence "$O_TELO" \
    --ont-coverage-summary "$O_COV_SUM" \
    --ont-annotated-gaps "$O_GAPS" \
    --centromeres "$CEN" \
    --telo-p-bed "$TELO_P" \
    --layout butterfly \
    --style paper \
    --format pdf \
    --output "$OUT_DIR/bTaeGut7_haplotype_ideogram_paper.pdf"

echo "=== Finished! Figures located in $OUT_DIR ==="
