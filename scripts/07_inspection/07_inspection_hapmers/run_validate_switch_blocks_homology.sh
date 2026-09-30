#!/usr/bin/env bash
# =============================================================================
# run_validate_switch_blocks_homology.sh
# 
# Runs the switch blocks homology validation script for all assemblies
# and saves outputs in results/inspection/07_inspection_hapmers/.
# =============================================================================

set -euo pipefail

# Directory paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

PY_BIN="/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/ds/bin/python3"
VALIDATION_SCRIPT="${SCRIPT_DIR}/validate_switch_blocks_homology.py"
OUT_DIR="${BASE_DIR}/results/inspection/07_inspection_hapmers"

mkdir -p "${OUT_DIR}"

echo "========================================================================"
echo "Running Switch Blocks Homology Validation across Assemblies"
echo "Repository Base: ${BASE_DIR}"
echo "Output Directory: ${OUT_DIR}"
echo "========================================================================"

"${PY_BIN}" "${VALIDATION_SCRIPT}" \
    --base-dir "${BASE_DIR}" \
    --out-dir "${OUT_DIR}"

echo "========================================================================"
echo "Validation complete!"
echo "Generated files:"
echo "  - Single:   ${OUT_DIR}/single_switch_blocks_homology_validation.tsv"
echo "  - Dual:     ${OUT_DIR}/dual_switch_blocks_homology_validation.tsv"
echo "  - Combined: ${OUT_DIR}/switch_blocks_homology_validation_combined.tsv"
echo "========================================================================"
