#!/usr/bin/env bash
# =============================================================================
# run_summarize_switch_block_lengths.sh
#
# Computes the total sum of switch block lengths in maternal and paternal
# scaffolds per assembly and per category.
# Saves outputs in results/inspection/07_inspection_hapmers/.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

PY_BIN="/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/ds/bin/python3"
SCRIPT="${SCRIPT_DIR}/summarize_switch_block_lengths.py"
OUT_DIR="${BASE_DIR}/results/inspection/07_inspection_hapmers"

mkdir -p "${OUT_DIR}"

echo "========================================================================"
echo "Computing Switch Block Lengths by Category & Haplotype"
echo "Repository Base: ${BASE_DIR}"
echo "Output Directory: ${OUT_DIR}"
echo "========================================================================"

"${PY_BIN}" "${SCRIPT}" \
    --base-dir "${BASE_DIR}" \
    --out-dir "${OUT_DIR}"

echo "========================================================================"
echo "Summary complete!"
echo "Generated files:"
echo "  - Details: ${OUT_DIR}/switch_blocks_details_with_haplotype_category.tsv"
echo "  - Summary: ${OUT_DIR}/switch_blocks_summary_by_category_and_haplotype.tsv"
echo "========================================================================"
