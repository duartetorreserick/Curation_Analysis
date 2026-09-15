#!/bin/bash
# DUAL: construye dual_pat.renamed.fasta / dual_mat.renamed.fasta
#   - los SUPER_* salen de dual_pat / dual_mat (tal cual)
#   - los unlocs salen del curado (H1 o H2) que indica el sufijo .H del padre
#   - los Scaffold_* se descartan
#   - todo se renombra a Mat_<nombre> / Pat_<nombre>, conservando el .H original
set -euo pipefail

BASE=/lustre/fs5/vgl/scratch/eduarte/curations/birds/bTaeGut7/redo_curation_analysis
D=$BASE/data/dual
OUT=${1:-$BASE/results/dual/01_formatting_assemblies_dual}
SAM=/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/vgp/bin/samtools
AWK=$BASE/scripts/01_formatting_assemblies/01_formatting_assemblies_dual/assign_unlocs.awk

# Localizados por patron: la extension de los curados varia (.fasta / .fast)
H1=$(ls "$D"/*hap1.cur* 2>/dev/null | grep -v '\.fai$' | head -1)
H2=$(ls "$D"/*hap2.cur* 2>/dev/null | grep -v '\.fai$' | head -1)
MAT=$D/mat/dual_mat.fasta
PAT=$D/pat/dual_pat.fasta
[[ -n $H1 && -n $H2 ]] || { echo "ERROR: no encuentro hap1.cur / hap2.cur en $D"; exit 1; }
echo "  H1: $(basename "$H1")"
echo "  H2: $(basename "$H2")"

for f in "$H1" "$H2" "$MAT" "$PAT"; do
    [[ -f $f ]] || { echo "ERROR: no existe $f"; exit 1; }
done
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
mkdir -p "$OUT"

# ── 0) Indexar lo que falte (data/dual no trae .fai) ────────────────────────
for f in "$H1" "$H2" "$MAT" "$PAT"; do
    [[ -f "$f.fai" ]] || { echo "--- indexando $(basename "$f") ---"; $SAM faidx "$f"; }
done

# ── 1) Plan: hap / fuente / nombre_original / nombre_nuevo ──────────────────
awk -F'\t' -v OFS='\t' -f "$AWK" \
    "$H1.fai" "$H2.fai" "$MAT.fai" "$PAT.fai" > "$OUT/hap_assignment.tsv"

# ── 2) Extraer y renombrar por (hap, fuente), luego reordenar ───────────────
for hap in pat mat; do
    [[ $hap == mat ]] && dual=$MAT || dual=$PAT
    : > "$TMP/$hap.all.fa"

    for src in DUAL HAP1cur HAP2cur; do
        case $src in
            DUAL)    fa=$dual ;;
            HAP1cur) fa=$H1 ;;
            HAP2cur) fa=$H2 ;;
        esac
        awk -F'\t' -v h="$hap" -v s="$src" '$1==h && $2==s {print $3}' \
            "$OUT/hap_assignment.tsv" > "$TMP/$hap.$src.regions"
        [[ -s "$TMP/$hap.$src.regions" ]] || continue

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

    $SAM faidx "$TMP/$hap.all.fa"
    awk -F'\t' -v h="$hap" '$1==h {print $4}' "$OUT/hap_assignment.tsv" > "$TMP/$hap.order"
    $SAM faidx -n 60 -r "$TMP/$hap.order" "$TMP/$hap.all.fa" > "$OUT/dual_${hap}.renamed.fasta"
    $SAM faidx "$OUT/dual_${hap}.renamed.fasta"

    printf '%s: %s secuencias, %s bp\n' "dual_${hap}.renamed.fasta" \
        "$(grep -c '^>' "$OUT/dual_${hap}.renamed.fasta")" \
        "$(awk -F'\t' '{t+=$2} END{print t}' "$OUT/dual_${hap}.renamed.fasta.fai")"
done
