#!/usr/bin/env bash
# run_individual_synteny.sh
#
# Generates high-resolution individual chromosome synteny QC figures for all chromosomes
# comparing HiFi Single, HiFi Dual, and ONT Dual assemblies against the T2T reference.
#
# Usage:
#   bash run_individual_synteny.sh [CHROMOSOMES]
#
# Examples:
#   bash run_individual_synteny.sh                   # Plots all 41 chromosomes
#   bash run_individual_synteny.sh "chr3,chr11,chr30" # Plots selected chromosomes
#

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

CHROMS="${1:-all}"

# T2T Reference Annotations
CEN="$BASE/data/t2t/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff"
PBED="$BASE/data/t2t/bTaeGut7.T2T.fasta_terminal_telomeres.bed"

# Chain and coverage directories
SINGLE_CHAIN_DIR="$BASE/results/single/02_chainpipeline_single/02_chainpipeline_single_strictLinearGap"
DUAL_CHAIN_DIR="$BASE/results/dual/02_chainpipeline_dual/02_chainpipeline_dual_strictLinearGap"
ONT_CHAIN_DIR="$BASE/results/ont/02_chainpipeline_ont/02_chainpipeline_ont_strictLinearGap"

# Output Directory
OUTDIR="$BASE/results/figures/10_individual_synteny"
mkdir -p "$OUTDIR"

echo "========================================================================"
echo "Generating Individual Synteny Figures"
echo "Target Chromosomes: $CHROMS"
echo "Output Directory:   $OUTDIR"
echo "========================================================================"

$PYTHON "$SCRIPT_DIR/plot_individual_synteny.py" \
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
  --single-bed          "$BASE/results/single/05_hapmers_single/single.switch_blocks.final.reoriented.bed" \
  --dual-bed            "$BASE/results/dual/05_hapmers_dual/dual.switch_blocks.final.reoriented.bed" \
  --ont-bed             "$BASE/results/ont/05_hapmers_ont/ont.switch_blocks.final.reoriented.bed" \
  --single-telomeres    "$BASE/results/single/04_telomeres_single/single_telomere_presence.tsv" \
  --dual-telomeres      "$BASE/results/dual/04_telomeres_dual/dual_telomere_presence.tsv" \
  --ont-telomeres       "$BASE/results/ont/04_telomeres_ont/ont_telomere_presence.tsv" \
  --telo-p-bed          "$PBED" \
  --centromeres         "$CEN" \
  --chromosomes         "$CHROMS" \
  --output-dir          "$OUTDIR" \
  --dpi                 300

echo "========================================================================"
echo "Done! Figures saved in: $OUTDIR"
echo "========================================================================"
