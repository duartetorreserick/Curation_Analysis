#!/usr/bin/env bash
# =============================================================================
# run_summarize_stats_by_category.sh
#
# Independent script to compute and summarize per-category genome assembly
# QC statistics for Single HiFi, Dual HiFi, and ONT Dual assemblies.
#
# Metrics computed:
#   - Collinear coverage %
#   - Non-collinear coverage %
#   - Hap-mer switch error %
#   - Telomere completeness %
#
# Outputs:
#   results/figures/10_figures_main/stats/bTaeGut7_stats_by_category_wide.tsv
#   results/figures/10_figures_main/stats/bTaeGut7_stats_by_category_long.tsv
# =============================================================================

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

# Chain and coverage directories
SINGLE_CHAIN_DIR="$BASE/results/single/02_chainpipeline_single/02_chainpipeline_single_strictLinearGap"
DUAL_CHAIN_DIR="$BASE/results/dual/02_chainpipeline_dual/02_chainpipeline_dual_strictLinearGap"
ONT_CHAIN_DIR="$BASE/results/ont/02_chainpipeline_ont/02_chainpipeline_ont_strictLinearGap"

# Output Paths
OUTDIR="$BASE/results/figures/10_figures_main/stats"
mkdir -p "$OUTDIR"

# Optional prefix passed as first argument (defaults to standard stats prefix)
OUTPUT_PREFIX="${1:-$OUTDIR/bTaeGut7_stats_by_category}"

echo "========================================================================"
echo "Computing per-category assembly QC statistics"
echo "Output prefix: $OUTPUT_PREFIX"
echo "========================================================================"

$PYTHON "$BASE/scripts/10_figures/10_figures_supp/summarize_stats_by_category.py" \
  --single-tsv          "$BASE/results/single/03_coverage_single/03_coverage_single_strictLinearGap/collinear.single.target_cov.tsv" \
  --dual-tsv            "$BASE/results/dual/03_coverage_dual/03_coverage_dual_strictLinearGap/collinear.dual.target_cov.tsv" \
  --single-chain        "$SINGLE_CHAIN_DIR/t2t.vs.single_strictGap.T2T.vs.ASM.target.collinear.chain" \
  --dual-chain          "$DUAL_CHAIN_DIR/t2t.vs.dual_strictGap.T2T.vs.ASM.target.collinear.chain" \
  --single-nc-chain     "$SINGLE_CHAIN_DIR/t2t.vs.single_strictGap.T2T.vs.ASM.target.non-collinear.chain" \
  --dual-nc-chain       "$DUAL_CHAIN_DIR/t2t.vs.dual_strictGap.T2T.vs.ASM.target.non-collinear.chain" \
  --single-unlocs-chain "$SINGLE_CHAIN_DIR/t2t.vs.single_strictGap.T2T.vs.ASM.unloc.rbest.chain" \
  --dual-unlocs-chain   "$DUAL_CHAIN_DIR/t2t.vs.dual_strictGap.T2T.vs.ASM.unloc.rbest.chain" \
  --ont-dual-chain      "$ONT_CHAIN_DIR/t2t.vs.ont_strictGap.T2T.vs.ASM.target.collinear.chain" \
  --ont-dual-nc-chain   "$ONT_CHAIN_DIR/t2t.vs.ont_strictGap.T2T.vs.ASM.target.non-collinear.chain" \
  --single-pairs        "$SINGLE_CHAIN_DIR/t2t.vs.single_strictGap.best_chrom_pairs.tsv" \
  --dual-pairs          "$DUAL_CHAIN_DIR/t2t.vs.dual_strictGap.best_chrom_pairs.tsv" \
  --single-bed          "$BASE/results/single/05_hapmers_single/single.switch_blocks.final.reoriented.bed" \
  --dual-bed            "$BASE/results/dual/05_hapmers_dual/dual.switch_blocks.final.reoriented.bed" \
  --single-telomere-presence "$BASE/results/single/04_telomeres_single/single_telomere_presence.tsv" \
  --dual-telomere-presence   "$BASE/results/dual/04_telomeres_dual/dual_telomere_presence.tsv" \
  --ont-telomere-presence    "$BASE/results/ont/04_telomeres_ont/ont_telomere_presence.tsv" \
  --output              "$OUTPUT_PREFIX"

echo "========================================================================"
echo "Done! Generated files:"
echo "  Stats TSV (wide): ${OUTPUT_PREFIX}_wide.tsv"
echo "  Stats TSV (long): ${OUTPUT_PREFIX}_long.tsv"
echo "========================================================================"
