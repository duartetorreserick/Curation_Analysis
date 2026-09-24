#!/bin/bash
set -euo pipefail

# ==============================================================================
# Run query-centric 1x netting and recover chains for uncovered assembly regions
# (ONT strictGap dataset)
# ==============================================================================

WORKDIR="/lustre/fs5/vgl/scratch/eduarte/curations/birds/bTaeGut7/redo_curation_analysis"
cd "${WORKDIR}"

ENV_PYTHON="/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/ds/bin/python3"
UCSC_BIN="/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/ucsc-tools/bin"
SCRIPT="scripts/02_chainpipeline/recover_query_uncovered_chains.py"

CHAIN_DIR="results/ont/02_chainpipeline_ont/02_chainpipeline_ont_strictLinearGap"
OUTDIR="results/ont/02_chainpipeline_ont/02_chainpipeline_uncovered_ont"

mkdir -p "${OUTDIR}"

echo "=== Running Query-Centric Recovery for Uncovered Assembly Regions (ONT) ==="
echo "Output Directory: ${OUTDIR}"

"${ENV_PYTHON}" "${SCRIPT}" \
  --sorted-1to1-chain "${CHAIN_DIR}/t2t.vs.ont_strictGap.T2T.vs.ASM.sorted.1to1.chain" \
  --collinear-chain "${CHAIN_DIR}/t2t.vs.ont_strictGap.T2T.vs.ASM.target.collinear.chain" \
  --noncollinear-chain "${CHAIN_DIR}/t2t.vs.ont_strictGap.T2T.vs.ASM.target.non-collinear.chain" \
  --asm-sizes "${CHAIN_DIR}/asm3_ONT_combined.sorted.reoriented.sizes" \
  --t2t-sizes "${CHAIN_DIR}/bTaeGut7.T2T.sizes" \
  --outdir "${OUTDIR}" \
  --min-uncovered-size 20000 \
  --ucsc-bin "${UCSC_BIN}"

echo "=== Run Complete ==="
