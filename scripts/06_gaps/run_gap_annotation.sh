#!/usr/bin/env bash
# run_gap_annotation.sh
# Annotates physical gaps (from gfastats) into ASSEMBLY and CURATION
# for DUAL, SINGLE, and ONT assemblies using first-principles coordinate matching.

set -euo pipefail

BASE=/lustre/fs5/vgl/scratch/eduarte/curations/birds/bTaeGut7/redo_curation_analysis
SCRIPT=$BASE/scripts/06_gaps/annotate_assembly_and_curation_gaps.py

echo "=== [1/3] Annotating Gaps: DUAL ==="
python3 "$SCRIPT" \
    --mode dual \
    --gaps-bed "$BASE/results/dual/06_gaps_dual/dual_combined.renamed.gaps.bed" \
    --fai "$BASE/results/dual/02_chainpipeline_dual/dual_combined.renamed.sorted.reoriented.fa.fai" \
    --scaffold-map "$BASE/data/dual/scaffold_hap_to_combined_super.tsv" \
    --agp "$BASE/data/dual/corrected.agp" \
    --output-annotated "$BASE/results/dual/06_gaps_dual/dual_combined.renamed.sorted.reoriented.annotated.gaps.bed" \
    --output-curation-bed "$BASE/results/dual/06_gaps_dual/dual_combined.renamed.sorted.reoriented.curation_gaps.bed" \
    --output-summary "$BASE/results/dual/06_gaps_dual/dual_gap_annotation_summary.txt"

echo "=== [2/3] Annotating Gaps: SINGLE ==="
python3 "$SCRIPT" \
    --mode single \
    --gaps-bed "$BASE/results/single/06_gaps_single/single_combined.renamed.sorted.reoriented.gaps.bed" \
    --fai "$BASE/results/single/02_chainpipeline_single/single_combined.renamed.sorted.reoriented.fa.fai" \
    --scaffold-map "$BASE/data/single/scaffold_hap_to_combined_super.tsv" \
    --hap1-agp "$BASE/data/single/Hap1.corrected.agp" \
    --hap2-agp "$BASE/data/single/Hap2.corrected.agp" \
    --output-annotated "$BASE/results/single/06_gaps_single/single_combined.renamed.sorted.reoriented.annotated.gaps.bed" \
    --output-curation-bed "$BASE/results/single/06_gaps_single/single_combined.renamed.sorted.reoriented.curation_gaps.bed" \
    --output-summary "$BASE/results/single/06_gaps_single/single_gap_annotation_summary.txt"

echo "=== [3/3] Annotating Gaps: ONT ==="
python3 "$SCRIPT" \
    --mode ont \
    --gaps-bed "$BASE/results/ont/06_gaps_ont/ont_combined.renamed.gaps.bed" \
    --fai "$BASE/results/ont/02_chainpipeline_ont/asm3_ONT_combined.sorted.reoriented.fa.fai" \
    --scaffold-map "$BASE/data/ont/scaffold_hap_to_combined_super.tsv" \
    --agp "$BASE/data/ont/corrected.agp" \
    --output-annotated "$BASE/results/ont/06_gaps_ont/asm3_ONT_combined.sorted.reoriented.annotated.gaps.bed" \
    --output-curation-bed "$BASE/results/ont/06_gaps_ont/asm3_ONT_combined.sorted.reoriented.curation_gaps.bed" \
    --output-summary "$BASE/results/ont/06_gaps_ont/ont_gap_annotation_summary.txt"

echo "=== All gap annotations complete ==="
