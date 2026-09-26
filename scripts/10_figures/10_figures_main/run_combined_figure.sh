#!/usr/bin/env bash
# run_combined_figure.sh
#
# Generates the multi-panel combined figure (Panels a–g) and category statistics
# for the bTaeGut7 assembly QC (Single HiFi, Dual HiFi, and ONT Dual) with
# synteny ideograms (Panels d, e, f).
#
# Panels:
#   Top row:
#     a: Per-category collinear coverage distribution (Macro, Micro, Dot split-violins)
#     b: Telomere completeness lollipop
#     c: Average genome coverage bars (Single, Dual, ONT)
#   Ideogram rows:
#     d: Macrochromosomes (Chr 4, 5, W maternal; Chr 4, 5 paternal)
#     e: Microchromosomes (Chr 11, 14, 18 butterfly)
#     f: Dot chromosomes (Chr 32, 34, 35 butterfly)
#   Bottom row:
#     g: Gene-model error rates (Frameshifts, Premature stop codons, PCGs)
#   Legend:
#     Comprehensive publication legend band at the bottom.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Conda environment / python executable
if command -v conda &>/dev/null && conda env list | grep -q "^ds "; then
    PYTHON="conda run -n ds python3"
elif [ -x "/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/ds/bin/python3" ]; then
    PYTHON="/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/ds/bin/python3"
else
    PYTHON="python3"
fi

# T2T Reference Annotations
CEN="$BASE/data/t2t/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff"
PBED="$BASE/data/t2t/bTaeGut7.T2T.fasta_terminal_telomeres.bed"
ANNOT="$BASE/data/annotations/annotation_HiFi_ONT.tsv"

# Chain and coverage directories
SINGLE_CHAIN_DIR="$BASE/results/single/02_chainpipeline_single/02_chainpipeline_single_strictLinearGap"
DUAL_CHAIN_DIR="$BASE/results/dual/02_chainpipeline_dual/02_chainpipeline_dual_strictLinearGap"
ONT_CHAIN_DIR="$BASE/results/ont/02_chainpipeline_ont/02_chainpipeline_ont_strictLinearGap"

# Output Paths
OUTDIR="$BASE/results/figures/10_figures_main"
STATS_DIR="$OUTDIR/stats"
mkdir -p "$OUTDIR" "$STATS_DIR"

DATE=$(date +%Y%m%d_%H%M%S)
OUT="$OUTDIR/bTaeGut7_combined_figure_synteny_${DATE}"

# Stats TSV (defaults to pre-generated wide stats file; can be overridden via $1)
STATS_FILE="${1:-$STATS_DIR/bTaeGut7_stats_by_category_wide.tsv}"

if [ ! -f "$STATS_FILE" ]; then
    # Fallback to the latest available *_wide.tsv if default not found
    LATEST_STATS=$(ls -t "$STATS_DIR"/*wide.tsv 2>/dev/null | head -n 1 || true)
    if [ -n "$LATEST_STATS" ] && [ -f "$LATEST_STATS" ]; then
        STATS_FILE="$LATEST_STATS"
    else
        echo "Warning: No stats TSV found at $STATS_FILE."
        echo "You can generate it by running: $BASE/scripts/10_figures/10_figures_supp/run_summarize_stats_by_category.sh"
    fi
fi

# =============================================================================
# Plot combined multi-panel figure (Panels a–g)
# =============================================================================
echo "========================================================================"
echo "Plotting combined figure (Panels a–g)"
echo "Using Stats TSV: $STATS_FILE"
echo "========================================================================"
$PYTHON "$SCRIPT_DIR/plot_combined_figure.py" \
  --t2t-tsv             "$BASE/results/single/03_coverage_single/03_coverage_single_strictLinearGap/collinear.single.target_cov.tsv" \
  --single-tsv          "$BASE/results/single/03_coverage_single/03_coverage_single_strictLinearGap/collinear.single.query_cov.tsv" \
  --dual-tsv            "$BASE/results/dual/03_coverage_dual/03_coverage_dual_strictLinearGap/collinear.dual.query_cov.tsv" \
  --ont-tsv             "$BASE/results/ont/03_coverage_ont/03_coverage_ont_strictLinearGap/collinear.ont.query_cov.tsv" \
  --single-pairs        "$SINGLE_CHAIN_DIR/t2t.vs.single_strictGap.best_chrom_pairs.tsv" \
  --dual-pairs          "$DUAL_CHAIN_DIR/t2t.vs.dual_strictGap.best_chrom_pairs.tsv" \
  --ont-pairs           "$ONT_CHAIN_DIR/t2t.vs.ont_strictGap.best_chrom_pairs.tsv" \
  --single-chain        "$SINGLE_CHAIN_DIR/t2t.vs.single_strictGap.T2T.vs.ASM.target.collinear.chain" \
  --single-nc-chain     "$SINGLE_CHAIN_DIR/t2t.vs.single_strictGap.T2T.vs.ASM.target.non-collinear.chain" \
  --dual-chain          "$DUAL_CHAIN_DIR/t2t.vs.dual_strictGap.T2T.vs.ASM.target.collinear.chain" \
  --dual-nc-chain       "$DUAL_CHAIN_DIR/t2t.vs.dual_strictGap.T2T.vs.ASM.target.non-collinear.chain" \
  --ont-chain           "$ONT_CHAIN_DIR/t2t.vs.ont_strictGap.T2T.vs.ASM.target.collinear.chain" \
  --ont-nc-chain        "$ONT_CHAIN_DIR/t2t.vs.ont_strictGap.T2T.vs.ASM.target.non-collinear.chain" \
  --single-rec-chain    "$BASE/results/single/02_chainpipeline_single/02_chainpipeline_uncovered_single/asm_recovered_uncovered.chain" \
  --dual-rec-chain      "$BASE/results/dual/02_chainpipeline_dual/02_chainpipeline_uncovered_dual/asm_recovered_uncovered.chain" \
  --ont-rec-chain       "$BASE/results/ont/02_chainpipeline_ont/02_chainpipeline_uncovered_ont/asm_recovered_uncovered.chain" \
  --single-gaps         "$BASE/results/single/06_gaps_single/single_combined.renamed.sorted.reoriented.annotated.gaps.bed" \
  --dual-gaps           "$BASE/results/dual/06_gaps_dual/dual_combined.renamed.sorted.reoriented.annotated.gaps.bed" \
  --ont-gaps            "$BASE/results/ont/06_gaps_ont/asm3_ONT_combined.sorted.reoriented.annotated.gaps.bed" \
  --single-bed          "$BASE/results/single/05_hapmers_single/single.switch_blocks.final.bed" \
  --dual-bed            "$BASE/results/dual/05_hapmers_dual/dual.switch_blocks.final.bed" \
  --ont-bed             "$BASE/results/ont/05_hapmers_ont/ont.switch_blocks.final.bed" \
  --single-telomeres    "$BASE/results/single/04_telomeres_single/single_telomere_presence.tsv" \
  --dual-telomeres      "$BASE/results/dual/04_telomeres_dual/dual_telomere_presence.tsv" \
  --ont-telomeres       "$BASE/results/ont/04_telomeres_ont/ont_telomere_presence.tsv" \
  --single-cov-sum      "$BASE/results/single/03_coverage_single/03_coverage_single_strictLinearGap/single.coverage_summary.tsv" \
  --dual-cov-sum        "$BASE/results/dual/03_coverage_dual/03_coverage_dual_strictLinearGap/dual.coverage_summary.tsv" \
  --ont-cov-sum         "$BASE/results/ont/03_coverage_ont/03_coverage_ont_strictLinearGap/ont.coverage_summary.tsv" \
  --telo-p-bed          "$PBED" \
  --centromeres         "$CEN" \
  --stats-tsv           "$STATS_FILE" \
  --annotation-tsv      "$ANNOT" \
  --output              "$OUT" \
  --dpi                 300

echo "========================================================================"
echo "Done!"
echo "Generated files:"
echo "  Figure (PNG):     ${OUT}.png"
echo "  Figure (PDF):     ${OUT}.pdf"
echo "========================================================================"
