# Cruza los .fai por LONGITUD (los nombres SUPER_* se repiten entre mat y pat).
# Uso:  awk -f match_fai.awk -v MAT=mat.tsv -v PAT=pat.tsv  REF1.fai [REF2.fai ...] single_mat.fai single_pat.fai
function ver(f,   b) {                      # etiqueta de version segun el nombre del archivo
    b = f; sub(/.*\//, "", b)
    if (b ~ /hap1\.cur/) return "HAP1cur"
    if (b ~ /hap2\.cur/) return "HAP2cur"
    if (b ~ /hap1/)      return "HAP1"
    if (b ~ /hap2/)      return "HAP2"
    sub(/\.f(ast)?a\.fai$/, "", b); return b
}

# ---- 1) archivos de referencia: indice longitud->nombre + unlocs por padre ----
FILENAME !~ /single_(mat|pat)/ {
    v = ver(FILENAME)
    if (!(v in vseen)) { vseen[v] = 1; vers[++nv] = v }
    byLen[v, $2] = $1
    if ($1 ~ /_unloc_[0-9]+$/) {
        p = $1; sub(/_unloc_[0-9]+$/, "", p)
        ul[v, p] = (ul[v, p] == "" ? "" : ul[v, p] ",") $1 ":" $2
        nul[v, p]++
    }
    next
}

# ---- 2) single_mat / single_pat: buscar el scaffold equivalente ----
FNR == 1 && !hdr[FILENAME ~ /single_mat/ ? "M" : "P"]++ {
    h = "#single_scaffold" OFS "length"
    for (i = 1; i <= nv; i++) h = h OFS vers[i] "_scaffold" OFS vers[i] "_unlocs"
    print h > (FILENAME ~ /single_mat/ ? MAT : PAT)
}
{
    out  = (FILENAME ~ /single_mat/) ? MAT : PAT
    line = $1 OFS $2
    for (i = 1; i <= nv; i++) {
        v = vers[i]
        m = ((v, $2) in byLen) ? byLen[v, $2] : "-"
        u = (m != "-" && ((v, m) in ul))     ? ul[v, m]  : "-"
        line = line OFS m OFS u
    }
    print line > out
}
