#!/usr/bin/env bash
# run_combined_figure_boxplot_v2_updated_categories_panelg.sh
#
# Generates the multi-panel combined figure (Panels a–g) and category statistics
# for the bTaeGut7 assembly QC (Single HiFi, Dual HiFi, and ONT Dual).
#
# Panels:
#   Top row:
#     a: Per-category boxplot / stats (Macro, Micro, Dot)
#     b: Telomere completeness lollipop
#     c: Average genome coverage bars (Single, Dual, ONT)
#   Ideogram rows:
#     d: Macrochromosomes (Chr 4, W, Z butterfly)
#     e: Microchromosomes (Chr 10, Chr 15 butterfly)
#     f: Dot chromosomes (Chr 31, Chr 32 butterfly)
#   Bottom row:
#     g: Gene-model error rates dot plot (Frameshifts, Premature stop codons, PCGs)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(cd "$SCRIPT_DIR/../.." && pwd)"

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
PBED="$BASE/results/single/04_telomeres_single/p.terminal.telomeres.gap.bed"
ANNOT="$BASE/data/annotations/annotation_HiFi_ONT.tsv"

# Output Paths
OUTDIR="$BASE/results/figures"
STATS_DIR="$BASE/results/stats"
mkdir -p "$OUTDIR" "$STATS_DIR"

STATS_PREFIX="$STATS_DIR/bTaeGut7_stats_by_category_updated_categories"
OUT="$OUTDIR/bTaeGut7_combined_boxplot_v2_updated_categories_panelg"

# =============================================================================
# Step 1 — Regenerate per-category stats TSV with updated classification
# =============================================================================
echo "========================================================================"
echo "Step 1: Regenerating stats TSV with updated categories"
echo "========================================================================"
$PYTHON "$SCRIPT_DIR/summarize_stats_by_category_v2.py" \
  --single-tsv          "$BASE/results/single/03_coverage_single/chain_cov.target_cov.tsv" \
  --dual-tsv            "$BASE/results/dual/03_coverage_dual/chain_cov.target_cov.tsv" \
  --single-chain        "$BASE/results/single/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.target.collinear.chain" \
  --dual-chain          "$BASE/results/dual/02_chainpipeline_dual/t2t.vs.dual.T2T.vs.ASM.target.collinear.chain" \
  --single-nc-chain     "$BASE/results/single/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.target.non-collinear.chain" \
  --dual-nc-chain       "$BASE/results/dual/02_chainpipeline_dual/t2t.vs.dual.T2T.vs.ASM.target.non-collinear.chain" \
  --single-unlocs-chain "$BASE/results/single/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.unloc.rbest.chain" \
  --dual-unlocs-chain   "$BASE/results/dual/02_chainpipeline_dual/t2t.vs.dual.T2T.vs.ASM.unloc.rbest.chain" \
  --ont-dual-chain      "$BASE/results/ont/02_chainpipeline_ont/t2t.vs.ont.T2T.vs.ASM.target.collinear.chain" \
  --ont-dual-nc-chain   "$BASE/results/ont/02_chainpipeline_ont/t2t.vs.ont.T2T.vs.ASM.target.non-collinear.chain" \
  --single-pairs        "$BASE/results/single/02_chainpipeline_single/t2t.vs.single.best_chrom_pairs.tsv" \
  --dual-pairs          "$BASE/results/dual/02_chainpipeline_dual/t2t.vs.dual.best_chrom_pairs.tsv" \
  --single-bed          "$BASE/results/single/05_hapmers_single/single.switch_blocks.final.bed" \
  --dual-bed            "$BASE/results/dual/05_hapmers_dual/dual.switch_blocks.final.bed" \
  --single-telomere-presence "$BASE/results/single/04_telomeres_single/single_telomere_presence.tsv" \
  --dual-telomere-presence   "$BASE/results/dual/04_telomeres_dual/dual_telomere_presence.tsv" \
  --ont-telomere-presence    "$BASE/results/ont/04_telomeres_ont/ont_telomere_presence.tsv" \
  --output              "$STATS_PREFIX"

# =============================================================================
# Step 2 — Plot combined multi-panel figure
# =============================================================================
echo "========================================================================"
echo "Step 2: Plotting combined figure (Panels a–g)"
echo "========================================================================"
$PYTHON "$SCRIPT_DIR/plot_combined_figure_boxplot_v2_updated_categories_panelg.py" \
  --single-tsv          "$BASE/results/single/03_coverage_single/chain_cov.target_cov.tsv" \
  --dual-tsv            "$BASE/results/dual/03_coverage_dual/chain_cov.target_cov.tsv" \
  --single-chain        "$BASE/results/single/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.target.collinear.chain" \
  --dual-chain          "$BASE/results/dual/02_chainpipeline_dual/t2t.vs.dual.T2T.vs.ASM.target.collinear.chain" \
  --single-nc-chain     "$BASE/results/single/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.target.non-collinear.chain" \
  --dual-nc-chain       "$BASE/results/dual/02_chainpipeline_dual/t2t.vs.dual.T2T.vs.ASM.target.non-collinear.chain" \
  --single-unlocs-chain "$BASE/results/single/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.unloc.rbest.chain" \
  --dual-unlocs-chain   "$BASE/results/dual/02_chainpipeline_dual/t2t.vs.dual.T2T.vs.ASM.unloc.rbest.chain" \
  --single-pairs        "$BASE/results/single/02_chainpipeline_single/t2t.vs.single.best_chrom_pairs.tsv" \
  --dual-pairs          "$BASE/results/dual/02_chainpipeline_dual/t2t.vs.dual.best_chrom_pairs.tsv" \
  --single-bed          "$BASE/results/single/05_hapmers_single/single.switch_blocks.final.bed" \
  --dual-bed            "$BASE/results/dual/05_hapmers_dual/dual.switch_blocks.final.bed" \
  --ont-dual-chain      "$BASE/results/ont/02_chainpipeline_ont/t2t.vs.ont.T2T.vs.ASM.target.collinear.chain" \
  --ont-dual-nc-chain   "$BASE/results/ont/02_chainpipeline_ont/t2t.vs.ont.T2T.vs.ASM.target.non-collinear.chain" \
  --single-telomere-presence "$BASE/results/single/04_telomeres_single/single_telomere_presence.tsv" \
  --dual-telomere-presence   "$BASE/results/dual/04_telomeres_dual/dual_telomere_presence.tsv" \
  --ont-telomere-presence    "$BASE/results/ont/04_telomeres_ont/ont_telomere_presence.tsv" \
  --single-coverage-summary  "$BASE/results/single/03_coverage_single/single.coverage_summary.tsv" \
  --dual-coverage-summary    "$BASE/results/dual/03_coverage_dual/dual.coverage_summary.tsv" \
  --ont-coverage-summary     "$BASE/results/ont/03_coverage_ont/ont.coverage_summary.tsv" \
  --telo-p-bed          "$PBED" \
  --centromeres         "$CEN" \
  --stats-tsv           "${STATS_PREFIX}_wide.tsv" \
  --annotation-tsv      "$ANNOT" \
  --output              "$OUT" \
  --no-timestamp \
  --style               paper \
  --format              png \
  --simplify

echo "========================================================================"
echo "Done!"
echo "Generated files:"
echo "  Stats TSV (wide): ${STATS_PREFIX}_wide.tsv"
echo "  Stats TSV (long): ${STATS_PREFIX}_long.tsv"
echo "  Figure (PNG):     ${OUT}.png"
echo "  Figure (PDF):     ${OUT}.pdf"
echo "  Figure (SVG):     ${OUT}.svg"
echo "========================================================================"
