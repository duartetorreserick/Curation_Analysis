#!/bin/bash
# =============================================================================
# run_synteny_ideogram_macro.sh
#
# Generates the Supplementary Macrochromosome Synteny Ideogram comparing T2T
# with HiFi Single, HiFi Dual, and ONT Dual assemblies in butterfly layout.
#
# Output:
#   results/figures/10_figures_supp/bTaeGut7_synteny_ideogram_macro_<TIMESTAMP>.png
#   results/figures/10_figures_supp/bTaeGut7_synteny_ideogram_macro_<TIMESTAMP>.pdf
# =============================================================================

set -euo pipefail

BASE="/lustre/fs5/vgl/scratch/eduarte/curations/birds/bTaeGut7/redo_curation_analysis"
PYTHON="/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/ds/bin/python3"
SCRIPT="${BASE}/scripts/10_figures/10_figures_supp/plot_synteny_ideogram_macro.py"
OUT_DIR="${BASE}/results/figures/10_figures_supp"

mkdir -p "$OUT_DIR"

TIMESTAMP="${1:-$(date +%Y%m%d_%H%M%S)}"
OUT_PNG="${OUT_DIR}/bTaeGut7_synteny_ideogram_macro_${TIMESTAMP}.png"
OUT_PDF="${OUT_DIR}/bTaeGut7_synteny_ideogram_macro_${TIMESTAMP}.pdf"

echo "=== Running Macrochromosome Synteny Ideogram ==="
echo "Timestamp: $TIMESTAMP"
echo "Output PNG: $OUT_PNG"
echo "Output PDF: $OUT_PDF"

SINGLE_CHAIN_DIR="$BASE/results/single/02_chainpipeline_single/02_chainpipeline_single_strictLinearGap"
DUAL_CHAIN_DIR="$BASE/results/dual/02_chainpipeline_dual/02_chainpipeline_dual_strictLinearGap"
ONT_CHAIN_DIR="$BASE/results/ont/02_chainpipeline_ont/02_chainpipeline_ont_strictLinearGap"

$PYTHON "$SCRIPT" \
    --t2t-tsv "$BASE/results/single/03_coverage_single/03_coverage_single_strictLinearGap/collinear.single.target_cov.tsv" \
    --single-tsv "$BASE/results/single/03_coverage_single/03_coverage_single_strictLinearGap/collinear.single.query_cov.tsv" \
    --dual-tsv "$BASE/results/dual/03_coverage_dual/03_coverage_dual_strictLinearGap/collinear.dual.query_cov.tsv" \
    --ont-tsv "$BASE/results/ont/03_coverage_ont/03_coverage_ont_strictLinearGap/collinear.ont.query_cov.tsv" \
    --single-pairs "$SINGLE_CHAIN_DIR/t2t.vs.single_strictGap.best_chrom_pairs.tsv" \
    --dual-pairs "$DUAL_CHAIN_DIR/t2t.vs.dual_strictGap.best_chrom_pairs.tsv" \
    --ont-pairs "$ONT_CHAIN_DIR/t2t.vs.ont_strictGap.best_chrom_pairs.tsv" \
    --single-chain "$SINGLE_CHAIN_DIR/t2t.vs.single_strictGap.T2T.vs.ASM.target.collinear.chain" \
    --single-nc-chain "$SINGLE_CHAIN_DIR/t2t.vs.single_strictGap.T2T.vs.ASM.target.non-collinear.chain" \
    --dual-chain "$DUAL_CHAIN_DIR/t2t.vs.dual_strictGap.T2T.vs.ASM.target.collinear.chain" \
    --dual-nc-chain "$DUAL_CHAIN_DIR/t2t.vs.dual_strictGap.T2T.vs.ASM.target.non-collinear.chain" \
    --ont-chain "$ONT_CHAIN_DIR/t2t.vs.ont_strictGap.T2T.vs.ASM.target.collinear.chain" \
    --ont-nc-chain "$ONT_CHAIN_DIR/t2t.vs.ont_strictGap.T2T.vs.ASM.target.non-collinear.chain" \
    --single-rec-chain "$BASE/results/single/02_chainpipeline_single/02_chainpipeline_uncovered_single/asm_recovered_uncovered.chain" \
    --dual-rec-chain "$BASE/results/dual/02_chainpipeline_dual/02_chainpipeline_uncovered_dual/asm_recovered_uncovered.chain" \
    --ont-rec-chain "$BASE/results/ont/02_chainpipeline_ont/02_chainpipeline_uncovered_ont/asm_recovered_uncovered.chain" \
    --single-gaps "$BASE/results/single/06_gaps_single/single_combined.renamed.sorted.reoriented.annotated.gaps.bed" \
    --dual-gaps "$BASE/results/dual/06_gaps_dual/dual_combined.renamed.sorted.reoriented.annotated.gaps.bed" \
    --ont-gaps "$BASE/results/ont/06_gaps_ont/asm3_ONT_combined.sorted.reoriented.annotated.gaps.bed" \
    --single-bed "$BASE/results/single/05_hapmers_single/single.switch_blocks.final.bed" \
    --dual-bed "$BASE/results/dual/05_hapmers_dual/dual.switch_blocks.final.bed" \
    --ont-bed "$BASE/results/ont/05_hapmers_ont/ont.switch_blocks.final.bed" \
    --single-telomeres "$BASE/results/single/04_telomeres_single/single_telomere_presence.tsv" \
    --dual-telomeres "$BASE/results/dual/04_telomeres_dual/dual_telomere_presence.tsv" \
    --ont-telomeres "$BASE/results/ont/04_telomeres_ont/ont_telomere_presence.tsv" \
    --centromeres "$BASE/data/t2t/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff" \
    --telo-p-bed "$BASE/data/t2t/bTaeGut7.T2T.fasta_terminal_telomeres.bed" \
    --single-cov-sum "$BASE/results/single/03_coverage_single/03_coverage_single_strictLinearGap/single.coverage_summary.tsv" \
    --dual-cov-sum "$BASE/results/dual/03_coverage_dual/03_coverage_dual_strictLinearGap/dual.coverage_summary.tsv" \
    --ont-cov-sum "$BASE/results/ont/03_coverage_ont/03_coverage_ont_strictLinearGap/ont.coverage_summary.tsv" \
    --output-png "$OUT_PNG" \
    --output-pdf "$OUT_PDF" \
    --dpi 300

echo "=== Generation Complete! ==="
