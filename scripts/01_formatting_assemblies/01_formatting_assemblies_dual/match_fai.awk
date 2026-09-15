# DUAL: lista, por scaffold de dual_mat / dual_pat, de que curado viene y que
# unlocs arrastra.
#
# A diferencia de la version single, aqui NO se cruza por longitud: los nombres
# ya traen el sufijo .H1/.H2, asi que el curado de origen se lee del nombre.
# La longitud se usa solo como comprobacion.
#
# Uso:
#   awk -F'\t' -v OFS='\t' -v MAT=mat.tsv -v PAT=pat.tsv -f match_fai.awk \
#       H1cur.fai H2cur.fai dual_mat.fasta.fai dual_pat.fasta.fai

function ver(f,   b) {
    b = f; sub(/.*\//, "", b)
    if (b ~ /hap1\.cur/) return "HAP1cur"
    if (b ~ /hap2\.cur/) return "HAP2cur"
    return "?"
}
function src_of(name,   h) {
    h = name; sub(/.*\./, "", h)
    return (h == "H1") ? "HAP1cur" : (h == "H2") ? "HAP2cur" : ""
}
function parent_of(name,   p) { p = name; sub(/_unloc_[0-9]+/, "", p); return p }
function idx_of(name,   n)    { n = name; sub(/.*_unloc_/, "", n); sub(/\..*/, "", n); return n + 0 }

# ---- 1) curados: longitudes y unlocs por padre ----
FILENAME !~ /dual_(mat|pat)/ {
    if ($1 ~ /^Scaffold_/) next
    v = ver(FILENAME)
    len[v, $1] = $2
    if ($1 ~ /_unloc_[0-9]+/) {
        p = parent_of($1); n = idx_of($1)
        ulname[v, p, n] = $1 ":" $2
        if (n > ulmax[v, p]) ulmax[v, p] = n
    }
    next
}

# ---- 2) dual_mat / dual_pat ----
FNR == 1 { out = (FILENAME ~ /dual_mat/) ? MAT : PAT
           print "#scaffold" OFS "length" OFS "curado" OFS "len_ok" OFS "n_unlocs" OFS "unlocs" > out }
{
    if ($1 ~ /^Scaffold_/) next
    out = (FILENAME ~ /dual_mat/) ? MAT : PAT
    src = src_of($1)
    ok  = ((src, $1) in len) ? (len[src, $1] == $2 ? "si" : "NO") : "ausente"

    u = ""; k = 0
    for (n = 1; n <= ulmax[src, $1]; n++)
        if ((src, $1, n) in ulname) { u = (u == "" ? "" : u ",") ulname[src, $1, n]; k++ }
    print $1, $2, src, ok, k, (u == "" ? "-" : u) > out
}
