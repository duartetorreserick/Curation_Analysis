#!/bin/bash
# =============================================================================
# run_synteny_ideogram_stacked.sh
#
# Generates the Supplementary Macrochromosome Stacked Multi-Alignment
# Synteny Ideogram comparing a single shared T2T reference chromosome at the top
# against the three assemblies (HiFi Single, HiFi Dual, ONT Dual) stacked below it.
#
# Output:
#   results/figures/10_figures_supp/bTaeGut7_synteny_ideogram_stacked_<TIMESTAMP>.png
#   results/figures/10_figures_supp/bTaeGut7_synteny_ideogram_stacked_<TIMESTAMP>.pdf
# =============================================================================

set -euo pipefail

BASE="/lustre/fs5/vgl/scratch/eduarte/curations/birds/bTaeGut7/redo_curation_analysis"
PYTHON="/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/ds/bin/python3"
SCRIPT="${BASE}/scripts/10_figures/10_figures_supp/plot_synteny_ideogram_stacked.py"
OUT_DIR="${BASE}/results/figures/10_figures_supp"

mkdir -p "$OUT_DIR"

TIMESTAMP="${1:-$(date +%Y%m%d_%H%M%S)}"
OUT_PNG="${OUT_DIR}/bTaeGut7_synteny_ideogram_stacked_${TIMESTAMP}.png"
OUT_PDF="${OUT_DIR}/bTaeGut7_synteny_ideogram_stacked_${TIMESTAMP}.pdf"

echo "=== Running Stacked Multi-Alignment Synteny Ideogram ==="
echo "Timestamp: $TIMESTAMP"
echo "Output PNG: $OUT_PNG"
echo "Output PDF: $OUT_PDF"

$PYTHON "$SCRIPT" \
    --t2t-tsv "$BASE/results/single/03_coverage_single/chain_cov.target_cov.tsv" \
    --single-tsv "$BASE/results/single/03_coverage_single/chain_cov.query_cov.tsv" \
    --dual-tsv "$BASE/results/dual/03_coverage_dual/chain_cov.query_cov.tsv" \
    --ont-tsv "$BASE/results/ont/03_coverage_ont/chain_cov.query_cov.tsv" \
    --single-pairs "$BASE/results/single/02_chainpipeline_single/t2t.vs.single.best_chrom_pairs.tsv" \
    --dual-pairs "$BASE/results/dual/02_chainpipeline_dual/t2t.vs.dual.best_chrom_pairs.tsv" \
    --ont-pairs "$BASE/results/ont/02_chainpipeline_ont/t2t.vs.ont.best_chrom_pairs.tsv" \
    --single-chain "$BASE/results/single/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.target.collinear.chain" \
    --single-nc-chain "$BASE/results/single/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.target.non-collinear.chain" \
    --dual-chain "$BASE/results/dual/02_chainpipeline_dual/t2t.vs.dual.T2T.vs.ASM.target.collinear.chain" \
    --dual-nc-chain "$BASE/results/dual/02_chainpipeline_dual/t2t.vs.dual.T2T.vs.ASM.target.non-collinear.chain" \
    --ont-chain "$BASE/results/ont/02_chainpipeline_ont/t2t.vs.ont.T2T.vs.ASM.target.collinear.chain" \
    --ont-nc-chain "$BASE/results/ont/02_chainpipeline_ont/t2t.vs.ont.T2T.vs.ASM.target.non-collinear.chain" \
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
    --output-png "$OUT_PNG" \
    --output-pdf "$OUT_PDF" \
    --dpi 300

echo "=== Generation Complete! ==="
