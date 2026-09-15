#!/bin/bash
# Construye single_pat.renamed.fasta / single_mat.renamed.fasta:
#   - los SUPER_* salen de single_pat / single_mat (tal cual)
#   - los unlocs salen del curado (H1 o H2) al que pertenece su SUPER padre,
#     identificado por longitud (assign_unlocs.awk)
#   - los Scaffold_* se descartan
#   - todo se renombra a Pat_SUPER_#[_unloc_#].H#  /  Mat_...
set -euo pipefail

BASE=/lustre/fs5/vgl/scratch/eduarte/curations/birds/bTaeGut7/redo_curation_analysis
D=$BASE/data/single
OUT=${1:-$BASE/results/single/01_formatting_assemblies_single}
SAM=/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/vgp/bin/samtools

H1=$D/bTaeGut7.H1.hap1.cur.20260324.fasta
H2=$D/bTaeGut7.H2.hap2.cur.20260324.fasta
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
mkdir -p "$OUT"

# ── 1) Plan: hap / fuente / nombre_original / nombre_nuevo ───────────────────
awk -F'\t' -v OFS='\t' -f "$BASE/scripts/assign_unlocs.awk" \
    "$H1.fai" "$H2.fai" "$D/single_mat.fasta.fai" "$D/single_pat.fasta.fai" \
    > "$OUT/hap_assignment.tsv"

# ── 2) Extraer y renombrar por (hap, fuente), luego reordenar ────────────────
for hap in pat mat; do
    single=$D/single_${hap}.fasta
    : > "$TMP/$hap.all.fa"

    for src in SINGLE HAP1cur HAP2cur; do
        case $src in
            SINGLE)  fa=$single ;;
            HAP1cur) fa=$H1 ;;
            HAP2cur) fa=$H2 ;;
        esac
        awk -F'\t' -v h="$hap" -v s="$src" '$1==h && $2==s {print $3}' \
            "$OUT/hap_assignment.tsv" > "$TMP/$hap.$src.regions"
        [[ -s "$TMP/$hap.$src.regions" ]] || continue

        # extraer, y reescribir el header con el nombre nuevo de esa fuente
        $SAM faidx -n 60 -r "$TMP/$hap.$src.regions" "$fa" \
        | awk -F'\t' -v h="$hap" -v s="$src" -v map="$OUT/hap_assignment.tsv" '
            BEGIN { while ((getline l < map) > 0) {
                        split(l, f, "\t")
                        if (f[1]==h && f[2]==s) new[f[3]] = f[4]
                    } }
            /^>/ { n = substr($0,2); sub(/[ \t].*/, "", n)
                   if (!(n in new)) { print "ERROR: header inesperado: " n > "/dev/stderr"; exit 1 }
                   print ">" new[n]; next }
            { print }
          ' >> "$TMP/$hap.all.fa"
    done

    # orden final = orden del plan (cada SUPER seguido de sus unlocs)
    $SAM faidx "$TMP/$hap.all.fa"
    awk -F'\t' -v h="$hap" '$1==h {print $4}' "$OUT/hap_assignment.tsv" > "$TMP/$hap.order"
    $SAM faidx -n 60 -r "$TMP/$hap.order" "$TMP/$hap.all.fa" > "$OUT/single_${hap}.renamed.fasta"
    $SAM faidx "$OUT/single_${hap}.renamed.fasta"

    printf '%s: %s secuencias, %s bp\n' "single_${hap}.renamed.fasta" \
        "$(grep -c '^>' "$OUT/single_${hap}.renamed.fasta")" \
        "$(awk -F'\t' '{t+=$2} END{print t}' "$OUT/single_${hap}.renamed.fasta.fai")"
done
