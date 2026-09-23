#!/bin/bash
# =============================================================================
# RunPipeline.sh  —  bTaeGut7 harmonized chain pipeline
#
# T2T is assumed to already be sanitized (spaces replaced with underscores).
#
# Stage order:
#   1   Index T2T + raw ASM  (parallel, immediate)
#   2   Initial FastGA       (raw ASM vs T2T)
#   3   Chromosome pairs     (Hungarian or provided)
#   4a  Sort ASM        ─┐  (parallel after pairs)
#   4b  Misorientation  ─┘
#   5   Reorient             (needs 4a + 4b)
#   6   Re-index reoriented  (for Final FastGA + chain tools)
#   7   Final FastGA         (reoriented ASM vs T2T → PSL)
#   8   Chain / Net / Synteny on BOTH sides (target + swapped query)
#       — separa unlocs: los SUPERs siguen al filtro 1-to-1, los unlocs a 8b
#   8b  Unlocs: rbest de los dos lados en su propio espacio (no compiten
#       contra su cromosoma padre en chainNet)
#   9   Reciprocal intersection + collinear / non-collinear split
#       — los unlocs nunca son colineales; solo los que caen en territorio
#         no reclamado por los SUPERs entran a non-collinear
#
# Usage:
#   bash RunPipeline.sh \
#     --prefix  PREFIX \
#     --outdir  /path/to/workdir \
#     --t2t     bTaeGut7.T2T.sanitized.fa \
#     --asm     assembly.fa \
#     [--pairs  pairs.tsv]    (skip Hungarian if provided)
# =============================================================================

set -euo pipefail

# ── Tool paths ────────────────────────────────────────────────────────────────
SCRIPTS="/lustre/fs5/vgl/scratch/eduarte/curations/birds/bTaeGut7/redo_curation_analysis/scripts"
CONDA_SH="/lustre/fs5/vgl/scratch/eduarte/miniconda3/etc/profile.d/conda.sh"
ENV_UCSC="/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/ucsc-tools"
ENV_DS="/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/ds"
FASTGA_PATH="/lustre/fs5/vgl/scratch/eduarte/softwares/fastGA/FASTGA-1.5"
GFASTATS="/lustre/fs5/vgl/store/eduarte/programs/gfastats/build/bin/gfastats"
SAMTOOLS="/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/vgp/bin/samtools"
SORT_SCRIPT="$SCRIPTS/02_chainpipeline/sort_fasta_pipeline.sh"
HUNGARIAN_PY="$SCRIPTS/02_chainpipeline/hungarian_pairs.py"
FILTER_PY="$SCRIPTS/02_chainpipeline/filter_collinear.py"
LINEAR_GAP="/lustre/fs5/vgl/scratch/eduarte/curations/birds/bTaeGut7/redo_curation_analysis/results/inspection/strict.linearGap"

# ── SLURM defaults ────────────────────────────────────────────────────────────
PARTITION="vgl_b"
ACCOUNT="vgl_condo_bank"
THREADS=32

# ── Parse arguments ───────────────────────────────────────────────────────────
PREFIX="" OUTDIR="" T2T_FA="" ASM_FA="" USER_PAIRS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --prefix) PREFIX="$2";     shift 2 ;;
        --outdir) OUTDIR="$2";     shift 2 ;;
        --t2t)    T2T_FA="$2";     shift 2 ;;
        --asm)    ASM_FA="$2";     shift 2 ;;
        --pairs)  USER_PAIRS="$2"; shift 2 ;;
        *) echo "ERROR: Unknown argument $1"; exit 1 ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
[[ -z "$PREFIX" ]] && { echo "ERROR: --prefix is required"; exit 1; }
[[ -z "$OUTDIR" ]] && { echo "ERROR: --outdir is required"; exit 1; }
[[ -z "$T2T_FA" ]] && { echo "ERROR: --t2t is required"; exit 1; }
[[ -z "$ASM_FA" ]] && { echo "ERROR: --asm is required"; exit 1; }
[[ ! -f "$T2T_FA" ]] && { echo "ERROR: T2T not found: $T2T_FA"; exit 1; }
[[ ! -f "$ASM_FA" ]] && { echo "ERROR: ASM not found: $ASM_FA"; exit 1; }
[[ ! -f "$SORT_SCRIPT"  ]] && { echo "ERROR: not found: $SORT_SCRIPT"; exit 1; }
[[ ! -x "$GFASTATS"     ]] && { echo "ERROR: not executable: $GFASTATS"; exit 1; }
[[ ! -x "$SAMTOOLS"     ]] && { echo "ERROR: not executable: $SAMTOOLS"; exit 1; }
[[ ! -f "$FILTER_PY"    ]] && { echo "ERROR: not found: $FILTER_PY"; exit 1; }
[[ ! -f "$HUNGARIAN_PY" && -z "$USER_PAIRS" ]] && \
    { echo "ERROR: not found: $HUNGARIAN_PY"; exit 1; }
[[ -n "$USER_PAIRS" && ! -f "$USER_PAIRS" ]] && \
    { echo "ERROR: Pairs file not found: $USER_PAIRS"; exit 1; }

# Resolve all inputs to absolute paths — SLURM jobs must not depend on cwd
T2T_FA=$(realpath "$T2T_FA")
ASM_FA=$(realpath "$ASM_FA")
[[ -n "$USER_PAIRS" ]] && USER_PAIRS=$(realpath "$USER_PAIRS")
TIMESTAMP=$(date +%Y%m%d-%H%M%S)
OUTDIR="${OUTDIR}_${TIMESTAMP}"
mkdir -p "$OUTDIR"
WORKDIR=$(realpath "$OUTDIR")
LOGDIR="${WORKDIR}/logs"
mkdir -p "$LOGDIR"

# ── Derived names ─────────────────────────────────────────────────────────────
strip_ext() { basename "$1" | sed 's/\.\(fa\|fna\|fasta\)$//'; }

T2T_BASE=$(strip_ext "$T2T_FA")
ASM_BASE=$(strip_ext "$ASM_FA")
ASM_SORTED_FA="${WORKDIR}/${ASM_BASE}.sorted.fa"
ASM_SORTED_BASE="${ASM_BASE}.sorted"
ASM_REOR_FA="${WORKDIR}/${ASM_BASE}.sorted.reoriented.fa"
ASM_REOR_BASE="${ASM_BASE}.sorted.reoriented"

INIT_NAME="${PREFIX}.T2T.vs.ASM.initial"
PAIRS="${WORKDIR}/${PREFIX}.best_chrom_pairs.tsv"
SAK="${WORKDIR}/${PREFIX}.sak.tsv"

echo "============================================="
echo " Pipeline: $PREFIX"
echo " Workdir:  $WORKDIR"
echo " T2T:      $T2T_FA"
echo " ASM:      $ASM_FA"
echo " Pairs:    ${USER_PAIRS:-(Hungarian from initial FastGA)}"
echo "============================================="

# ── SLURM submit helper ───────────────────────────────────────────────────────
submit() {
    local job_name=$1 log=$2 time=$3 mem=$4 dep=$5
    shift 5
    local dep_flag=""
    [[ -n "$dep" ]] && dep_flag="--dependency=afterok:${dep}"
    sbatch --parsable \
        --job-name="$job_name" \
        --output="${LOGDIR}/${log}" \
        --partition="$PARTITION" \
        --account="$ACCOUNT" \
        --time="$time" \
        --mem="$mem" \
        --cpus-per-task="$THREADS" \
        $dep_flag \
        "$@"
}

# =============================================================================
# STAGE 1: Index T2T and raw ASM
#   faToTwoBit → .2bit   (axtChain)
#   twoBitInfo → .sizes  (chainNet)
#   FAtoGDB + GIXmake    (FastGA)
# Both start immediately in parallel.
# =============================================================================
JID_IDX_T2T=$(submit "Index_T2T" "Index_T2T_%j.log" "4:00:00" "50G" "" \
    --wrap="
set -euo pipefail
export PATH='/usr/bin:/bin:/lustre/fs5/vgl/scratch/eduarte/miniconda3/bin:${FASTGA_PATH}:\$PATH'
set +u; source '${CONDA_SH}' && conda activate '${ENV_UCSC}'; set -u
cd '${WORKDIR}'
echo '--- Indexing T2T: ${T2T_BASE} ---'
faToTwoBit '${T2T_FA}' '${T2T_BASE}'.2bit
twoBitInfo '${T2T_BASE}'.2bit '${T2T_BASE}'.sizes
FAtoGDB    '${T2T_FA}' '${T2T_BASE}'
GIXmake    '${T2T_BASE}'
echo 'Done: ${T2T_BASE}.2bit / .sizes / .gix'
")

JID_IDX_ASM=$(submit "Index_ASM" "Index_ASM_%j.log" "4:00:00" "50G" "" \
    --wrap="
set -euo pipefail
export PATH='/usr/bin:/bin:/lustre/fs5/vgl/scratch/eduarte/miniconda3/bin:${FASTGA_PATH}:\$PATH'
set +u; source '${CONDA_SH}' && conda activate '${ENV_UCSC}'; set -u
cd '${WORKDIR}'
echo '--- Indexing raw ASM: ${ASM_BASE} ---'
faToTwoBit '${ASM_FA}' '${ASM_BASE}'.2bit
twoBitInfo '${ASM_BASE}'.2bit '${ASM_BASE}'.sizes
FAtoGDB    '${ASM_FA}' '${ASM_BASE}'
GIXmake    '${ASM_BASE}'
echo 'Done: ${ASM_BASE}.2bit / .sizes / .gix'
")

echo "Stage 1   Index T2T: $JID_IDX_T2T"
echo "Stage 1   Index ASM: $JID_IDX_ASM"

# =============================================================================
# STAGE 2: Initial FastGA alignment (raw ASM vs sanitized T2T)
# =============================================================================
JID_INIT=$(submit "FastGA_Init" "FastGA_Init_%j.log" "24:00:00" "200G" \
    "${JID_IDX_T2T}:${JID_IDX_ASM}" \
    --wrap="
set -euo pipefail
export PATH='/usr/bin:/bin:/lustre/fs5/vgl/scratch/eduarte/miniconda3/bin:${FASTGA_PATH}:\$PATH'
cd '${WORKDIR}'
echo '--- FastGA initial: ${INIT_NAME} ---'
FastGA -T${THREADS} -k '${ASM_BASE}'.gix '${T2T_BASE}'.gix -1:'${INIT_NAME}'.1aln
ALNtoPAF -T${THREADS} '${INIT_NAME}'.1aln > '${INIT_NAME}'.paf
echo 'Done: ${INIT_NAME}.paf'
")
echo "Stage 2   FastGA:    $JID_INIT"

# =============================================================================
# STAGE 3: Chromosome pairs
#   --pairs provided → copy as-is
#   --pairs absent   → Hungarian algorithm (scipy) on initial PAF
# =============================================================================
if [[ -n "$USER_PAIRS" ]]; then
    JID_PAIRS=$(submit "Pairs" "Pairs_%j.log" "0:10:00" "4G" "$JID_INIT" \
        --wrap="
set -euo pipefail
cd '${WORKDIR}'
cp '${USER_PAIRS}' '${PAIRS}'
echo \"Pairs (provided): ${PAIRS} (\$(wc -l < '${PAIRS}') pairs)\"
")
else
    JID_PAIRS=$(submit "Pairs" "Pairs_%j.log" "1:00:00" "32G" "$JID_INIT" \
        --wrap="
set -euo pipefail
cd '${WORKDIR}'
echo '--- Hungarian chromosome pairs ---'
'${ENV_DS}/bin/python3' '${HUNGARIAN_PY}' '${INIT_NAME}'.paf '${PAIRS}'
echo \"Pairs (Hungarian): ${PAIRS} (\$(wc -l < '${PAIRS}') pairs)\"
")
fi
echo "Stage 3   Pairs:     $JID_PAIRS"

# =============================================================================
# STAGE 4a: Sort ASM to T2T chromosome order  (parallel with Stage 4b)
# =============================================================================
JID_SORT=$(submit "Sort_ASM" "Sort_ASM_%j.log" "2:00:00" "16G" "$JID_PAIRS" \
    --wrap="
set -euo pipefail
export PATH='/usr/bin:/bin:/lustre/fs5/vgl/scratch/eduarte/miniconda3/bin:\$PATH'
set +u; source '${CONDA_SH}'; set -u
cd '${WORKDIR}'
echo '--- Sorting ASM to T2T order ---'
bash '${SORT_SCRIPT}' '${PAIRS}' '${T2T_FA}' '${ASM_FA}' '${ASM_SORTED_FA}' '${ASM_FA}.fai'
echo \"Done: ${ASM_SORTED_FA}  (\$(grep -c '^>' '${ASM_SORTED_FA}') sequences)\"
")
echo "Stage 4a  Sort ASM:  $JID_SORT"

# =============================================================================
# STAGE 4b: Detect misoriented sequences  (parallel with Stage 4a)
# Computes +/- strand alignment fraction per ASM sequence from initial PAF.
# Sequences with >50% minus-strand bases → SAK file (RVCP entries).
# =============================================================================
JID_MISO=$(submit "Misoriented" "Misoriented_%j.log" "1:00:00" "8G" "$JID_PAIRS" \
    --wrap="
set -euo pipefail
export PATH='/usr/bin:/bin:/lustre/fs5/vgl/scratch/eduarte/miniconda3/bin:\$PATH'
cd '${WORKDIR}'
echo '--- Detecting misoriented sequences ---'
awk -v min_pct=50 -v min_aln=50000 '
NR==FNR {
    key=\$1; sub(/_[^_]+$/,\"\",key); t2t_to_asm[key]=\$2; next
}
{
    asm_name=\$1; t2t_name=\$6; strand=\$5; bases=\$10   # PAF: \$5=strand, \$10=matches
    if (bases < min_aln) next
    t2t_key=t2t_name; sub(/_[^_]+$/,\"\",t2t_key)
    if (t2t_to_asm[t2t_key]==asm_name) {
        seen[asm_name]=1
        if (strand==\"+\") plus[asm_name]+=bases
        if (strand==\"-\") minus[asm_name]+=bases
    }
}
END {
    for (s in seen) {
        total=plus[s]+minus[s]
        if (total>0 && (minus[s]/total*100)>=min_pct)
            print \"RVCP\t\"s
    }
}' '${PAIRS}' '${INIT_NAME}'.paf | sort -k2,2 > '${SAK}'
echo \"Done: ${SAK} (\$(wc -l < '${SAK}') sequences to reorient)\"
cat '${SAK}'
")
echo "Stage 4b  Misori:    $JID_MISO"

# =============================================================================
# STAGE 5: Reorient sorted ASM using SAK file
# Waits for both Sort (ASM_SORTED_FA) and Misorientation (SAK).
# =============================================================================
JID_REOR=$(submit "Reorient" "Reorient_%j.log" "4:00:00" "50G" \
    "${JID_SORT}:${JID_MISO}" \
    --wrap="
set -euo pipefail
export PATH='/usr/bin:/bin:/lustre/fs5/vgl/scratch/eduarte/miniconda3/bin:\$PATH'
cd '${WORKDIR}'
echo '--- Reorienting sorted ASM ---'
if [[ -s '${SAK}' ]]; then
    '${GFASTATS}' '${ASM_SORTED_FA}' -k '${SAK}' -o '${ASM_REOR_FA}'
    echo 'Reorientation applied from ${SAK}'
else
    cp '${ASM_SORTED_FA}' '${ASM_REOR_FA}'
    echo 'No misoriented sequences — copied as-is'
fi
'${SAMTOOLS}' faidx '${ASM_REOR_FA}'
echo \"Done: ${ASM_REOR_FA}\"
")
echo "Stage 5   Reorient:  $JID_REOR"

# =============================================================================
# STAGE 6: Re-index reoriented ASM
# Produces the indices used by Final FastGA and all chain tools.
# =============================================================================
JID_REIDX=$(submit "Reindex_ASM" "Reindex_ASM_%j.log" "4:00:00" "50G" "$JID_REOR" \
    --wrap="
set -euo pipefail
export PATH='/usr/bin:/bin:/lustre/fs5/vgl/scratch/eduarte/miniconda3/bin:${FASTGA_PATH}:\$PATH'
set +u; source '${CONDA_SH}' && conda activate '${ENV_UCSC}'; set -u
cd '${WORKDIR}'
echo '--- Re-indexing reoriented ASM: ${ASM_REOR_BASE} ---'
faToTwoBit '${ASM_REOR_FA}' '${ASM_REOR_BASE}'.2bit
twoBitInfo '${ASM_REOR_BASE}'.2bit '${ASM_REOR_BASE}'.sizes
FAtoGDB    '${ASM_REOR_FA}' '${ASM_REOR_BASE}'
GIXmake    '${ASM_REOR_BASE}'
echo 'Done: ${ASM_REOR_BASE}.2bit / .sizes / .gix'
")
echo "Stage 6   Reindex:   $JID_REIDX"

# =============================================================================
# STAGE 7: Final FastGA alignment (reoriented ASM vs T2T)
# Produces PSL needed by axtChain in Stage 8.
# Depends on Stage 6 (reoriented indices) + Stage 1 (T2T index already done).
# =============================================================================
NAME="${PREFIX}.T2T.vs.ASM"

JID_FASTGA=$(submit "FastGA_Final" "FastGA_Final_%j.log" "24:00:00" "200G" \
    "${JID_REIDX}:${JID_IDX_T2T}" \
    --wrap="
set -euo pipefail
export PATH='/usr/bin:/bin:/lustre/fs5/vgl/scratch/eduarte/miniconda3/bin:${FASTGA_PATH}:\$PATH'
cd '${WORKDIR}'
echo '--- FastGA final: ${NAME} ---'
FastGA -T${THREADS} -k '${ASM_REOR_BASE}'.gix '${T2T_BASE}'.gix -1:'${NAME}'.1aln
ALNtoPAF -T${THREADS} '${NAME}'.1aln > '${NAME}'.paf
ALNtoPSL -T${THREADS} '${NAME}'.1aln > '${NAME}'.psl
echo 'Done: ${NAME}.psl'
")
echo "Stage 7   FastGA Final: $JID_FASTGA"

# =============================================================================
# STAGE 8: Chain / Net / Synteny — built independently from BOTH sides
#
#   axtChain + chainSort + pair filter  →  sorted.1to1.chain
#
#   Target side:  net on T2T axis           → target.syn.net.chain
#   Query side:   chainSwap (ASM as target) → query.syn.net.chain
#
# The swap is required: netChainSubset needs the net and the chain to share a
# coordinate space, and query.net lives in ASM-as-target space.
# Chain ids survive chainSwap, so Stage 9 can intersect the two sets by id.
# =============================================================================
JID_CHAIN=$(submit "Chain" "Chain_%j.log" "24:00:00" "200G" \
    "${JID_FASTGA}:${JID_PAIRS}" \
    --wrap="
set -euo pipefail
export PATH='${ENV_UCSC}/bin:/usr/bin:/bin:/lustre/fs5/vgl/scratch/eduarte/miniconda3/bin:\$PATH'
cd '${WORKDIR}'

echo '--- axtChain ---'
axtChain -linearGap='${LINEAR_GAP}' -psl '${NAME}'.psl \
  '${T2T_BASE}'.2bit '${ASM_REOR_BASE}'.2bit '${NAME}'.chain
chainSort '${NAME}'.chain '${NAME}'.sorted.chain

echo '--- Separar unlocs (no colocados) de los SUPERs ---'
awk '/^chain/ {keep = (\$8 !~ /_unloc_/)} keep' \
  '${NAME}'.sorted.chain > '${NAME}'.placed.chain
awk 'NR==FNR {t2t_key=\$1; sub(/_[^_]+$/,\"\",t2t_key); pairs[t2t_key,\$2]=1; next}
     /^chain/ {t2t_key=\$3; sub(/_[^_]+$/,\"\",t2t_key); q=\$8
               if (q ~ /_unloc_/) {parent=q; sub(/_unloc_[0-9]+/,\"\",parent)
                                   keep=((t2t_key SUBSEP parent) in pairs)}
               else keep=0} keep' \
  '${PAIRS}' '${NAME}'.sorted.chain > '${NAME}'.unloc.chain
echo \"  colocados: \$(grep -c '^chain' '${NAME}'.placed.chain) | unloc: \$(grep -c '^chain' '${NAME}'.unloc.chain || true)\"

echo '--- 1-to-1 chromosome pair filter (solo colocados) ---'
awk 'NR==FNR {t2t_key=\$1; sub(/_[^_]+$/,\"\",t2t_key); pairs[t2t_key,\$2]=1; next}
     /^chain/ {t2t_key=\$3; sub(/_[^_]+$/,\"\",t2t_key);
               keep=(t2t_key SUBSEP \$8 in pairs)} keep' \
  '${PAIRS}' '${NAME}'.placed.chain > '${NAME}'.sorted.1to1.chain
echo \"  1-to-1 chains: \$(grep -c '^chain' '${NAME}'.sorted.1to1.chain)\"

echo '--- Target-side net (T2T as target) ---'
chainNet '${NAME}'.sorted.1to1.chain \
  '${T2T_BASE}'.sizes '${ASM_REOR_BASE}'.sizes \
  '${NAME}'.target.net /dev/null
netSyntenic '${NAME}'.target.net '${NAME}'.target.syn.net
netChainSubset '${NAME}'.target.syn.net '${NAME}'.sorted.1to1.chain stdout \
  | chainStitchId stdin '${NAME}'.target.syn.net.chain
echo \"  target syn chains: \$(grep -c '^chain' '${NAME}'.target.syn.net.chain)\"

echo '--- Cascada rbest: el net del query opera sobre la salida del target ---'
# Reciprocidad estructural: una cadena sobrevive solo si el net la eligio en el
# eje T2T y DESPUES en el eje del ensamblaje, sobre el conjunto ya reducido.
chainStitchId '${NAME}'.target.syn.net.chain stdout \
  | chainSwap stdin stdout | chainSort stdin stdout > '${NAME}'.q.tBest.chain
chainPreNet '${NAME}'.q.tBest.chain '${ASM_REOR_BASE}'.sizes '${T2T_BASE}'.sizes stdout \
  | chainNet -minSpace=1 -minScore=0 stdin '${ASM_REOR_BASE}'.sizes '${T2T_BASE}'.sizes stdout /dev/null \
  | netSyntenic stdin stdout > '${NAME}'.q.rbest.net
netChainSubset '${NAME}'.q.rbest.net '${NAME}'.q.tBest.chain stdout \
  | chainStitchId stdin stdout > '${NAME}'.q.rbest.chain
chainSwap '${NAME}'.q.rbest.chain stdout | chainSort stdin stdout > '${NAME}'.t.rbest.chain
chainNet -minSpace=1 -minScore=0 '${NAME}'.t.rbest.chain \
  '${T2T_BASE}'.sizes '${ASM_REOR_BASE}'.sizes '${NAME}'.t.rbest.net /dev/null
netChainSubset '${NAME}'.t.rbest.net '${NAME}'.t.rbest.chain '${NAME}'.target.rbest.chain
echo \"  rbest chains:      \$(grep -c '^chain' '${NAME}'.target.rbest.chain)\"

echo 'Done: ${NAME}.target.rbest.chain'
")
echo "Stage 8   Chain/Net:    $JID_CHAIN"

# =============================================================================
# STAGE 8b: Unlocs — cobertura única de los dos lados, en su propio espacio
#
# Los unlocs se netean SOLOS (su chain file no contiene SUPERs), así que nunca
# compiten contra su cromosoma padre en chainNet. Receta rbest de UCSC, igual
# que RunPipelineCombined_Unlocs.sh:
#   7a stitch/swap/sort  7b preNet/net/syntenic  7c subset/stitch
#   7d swap/sort         7e net (eje T2T)        7f subset
# Salida: unloc.rbest.chain — se clasifica en Stage 9, nunca como colineal.
# =============================================================================
JID_UNLOC=$(submit "Unloc_Rbest" "Unloc_Rbest_%j.log" "8:00:00" "100G" "$JID_CHAIN" \
    --wrap="
set -euo pipefail
export PATH='${ENV_UCSC}/bin:/usr/bin:/bin:/lustre/fs5/vgl/scratch/eduarte/miniconda3/bin:\$PATH'
cd '${WORKDIR}'

if [[ ! -s '${NAME}'.unloc.chain ]]; then
    : > '${NAME}'.unloc.rbest.chain
    echo 'Sin unlocs en el ensamblaje — nada que hacer'
    exit 0
fi

T2T_SIZES='${T2T_BASE}.sizes'
ASM_SIZES='${ASM_REOR_BASE}.sizes'
P=unloc_rbest_tmp

echo '--- 7a: chainStitchId | chainSwap | chainSort ---'
chainStitchId '${NAME}'.unloc.chain stdout | chainSwap stdin stdout | chainSort stdin stdout \
    > \"\${P}.query.tBest.chain\"

echo '--- 7b: chainPreNet | chainNet | netSyntenic ---'
chainPreNet \"\${P}.query.tBest.chain\" \"\$ASM_SIZES\" \"\$T2T_SIZES\" stdout \
    | chainNet -minSpace=1 -minScore=0 stdin \"\$ASM_SIZES\" \"\$T2T_SIZES\" stdout /dev/null \
    | netSyntenic stdin stdout \
    > \"\${P}.query.rbest.net\"

echo '--- 7c: netChainSubset | chainStitchId ---'
netChainSubset \"\${P}.query.rbest.net\" \"\${P}.query.tBest.chain\" stdout \
    | chainStitchId stdin stdout \
    > \"\${P}.query.rbest.chain\"

echo '--- 7d: chainSwap | chainSort ---'
chainSwap \"\${P}.query.rbest.chain\" stdout | chainSort stdin stdout \
    > \"\${P}.target.rbest.chain\"

echo '--- 7e: chainNet (eje T2T) ---'
chainNet -minSpace=1 -minScore=0 \"\${P}.target.rbest.chain\" \"\$T2T_SIZES\" \"\$ASM_SIZES\" \
    \"\${P}.target.rbest.net\" /dev/null

echo '--- 7f: netChainSubset ---'
netChainSubset \"\${P}.target.rbest.net\" \"\${P}.target.rbest.chain\" '${NAME}'.unloc.rbest.chain

rm -f unloc_rbest_tmp.*
echo \"  unloc rbest chains: \$(grep -c '^chain' '${NAME}'.unloc.rbest.chain || true)\"
")
echo "Stage 8b  Unloc rbest:  $JID_UNLOC"

# =============================================================================
# STAGE 9: Reciprocal intersection + collinearity split
#
#   intersect chain ids  target.syn.net.chain ∩ query.syn.net.chain
#     → 1x coverage on both sides; everything else is noise (paralogs,
#       secondary alignments) and is dropped
#
#   then split the survivors by strand and target ordering:
#     qStrand '-'                 → non-collinear (inversion)
#     '+' but tStart < prev tEnd  → non-collinear (transposition)
#     '+' and in order            → collinear
#
#   Los unlocs entran por --unloc-chain y NUNCA son colineales. Solo pasan a
#   non-collinear los que caen en territorio de T2T que los SUPERs no reclaman;
#   los que caen sobre región ya cubierta salen a unloc.claimed.chain, porque
#   son duplicación / haplotipo alternativo, no reordenamiento.
# =============================================================================
JID_FILTER=$(submit "Collinear" "Collinear_%j.log" "2:00:00" "32G" \
    "${JID_CHAIN}:${JID_UNLOC}" \
    --wrap="
set -euo pipefail
cd '${WORKDIR}'
echo '--- Collinearity split ---'
'${ENV_DS}/bin/python3' '${FILTER_PY}' \
  '${NAME}'.target.rbest.chain \
  '${NAME}'.target.collinear.chain \
  '${NAME}'.target.non-collinear.chain \
  --unloc-chain '${NAME}'.unloc.rbest.chain \
  --unloc-claimed-out '${NAME}'.unloc.claimed.chain \
  --max-claim-overlap 0.5
echo 'Done: ${NAME}.target.{collinear,non-collinear}.chain + ${NAME}.unloc.claimed.chain'
")
echo "Stage 9   Collinear:    $JID_FILTER"

echo ""
echo "============================================="
echo " All jobs submitted."
echo " Final outputs:"
echo "   ${NAME}.target.collinear.chain"
echo "   ${NAME}.target.non-collinear.chain   (incluye unlocs no reclamados)"
echo "   ${NAME}.unloc.claimed.chain          (unlocs sobre región ya cubierta)"
echo "============================================="
echo " Dependency graph:"
echo "  1   Index T2T   ($JID_IDX_T2T)  ┐ parallel"
echo "  1   Index ASM   ($JID_IDX_ASM)  ┘"
echo "  2   FastGA      ($JID_INIT)      ← 1"
echo "  3   Pairs       ($JID_PAIRS)     ← 2"
echo "  4a  Sort        ($JID_SORT)      ← 3  ┐ parallel"
echo "  4b  Misori      ($JID_MISO)      ← 3  ┘"
echo "  5   Reorient    ($JID_REOR)      ← 4a + 4b"
echo "  6   Reindex     ($JID_REIDX)     ← 5"
echo "  7   FastGA Final($JID_FASTGA)    ← 6 + 1"
echo "  8   Chain/Net   ($JID_CHAIN)     ← 7 + 3"
echo "  8b  Unloc rbest ($JID_UNLOC)     ← 8"
echo "  9   Collinear   ($JID_FILTER)    ← 8 + 8b"
echo ""
echo " Monitor: squeue -j ${JID_IDX_T2T},${JID_IDX_ASM},${JID_INIT},${JID_PAIRS},${JID_SORT},${JID_MISO},${JID_REOR},${JID_REIDX},${JID_FASTGA},${JID_CHAIN},${JID_UNLOC},${JID_FILTER}"
