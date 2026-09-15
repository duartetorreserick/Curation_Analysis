# Asigna cada unloc de los ensamblajes curados al haplotipo (mat/pat) de su SUPER padre.
# El padre se identifica por LONGITUD, porque los nombres SUPER_* se repiten entre mat y pat.
# Los Scaffold_* se descartan.
#
# Uso:
#   awk -F'\t' -v OFS='\t' -f assign_unlocs.awk \
#       H1cur.fai H2cur.fai single_mat.fai single_pat.fai
#
# Salida (stdout), en el orden final deseado del FASTA:
#   hap <TAB> fuente <TAB> nombre_original <TAB> nombre_nuevo
#   fuente = SINGLE | HAP1cur | HAP2cur   (de que FASTA extraer la secuencia)
#   nombre_nuevo = Pat_SUPER_3.H1 / Mat_SUPER_3_unloc_1.H2  (.H# = curado de origen)

function ver(f,   b) {
    b = f; sub(/.*\//, "", b)
    if (b ~ /hap1\.cur/) return "HAP1cur"
    if (b ~ /hap2\.cur/) return "HAP2cur"
    return "?"
}
function rename(hap, name, src) {
    return (hap == "pat" ? "Pat_" : "Mat_") name "." (src == "HAP1cur" ? "H1" : "H2")
}

# ---- 1) curados: indice longitud->nombre, y unlocs indexados por su numero ----
FILENAME !~ /single_(mat|pat)/ {
    v = ver(FILENAME)
    if ($1 ~ /^Scaffold_/) next
    byLen[v, $2] = $1
    if (match($1, /_unloc_[0-9]+$/)) {
        p = substr($1, 1, RSTART - 1)
        n = substr($1, RSTART + 7) + 0
        ulname[v, p, n] = $1
        if (n > ulmax[v, p]) ulmax[v, p] = n
    }
    next
}

# ---- 2) single_mat / single_pat: emitir el SUPER y luego sus unlocs ----
{
    if ($1 ~ /^Scaffold_/) { skipped++; next }
    hap = (FILENAME ~ /single_mat/) ? "mat" : "pat"

    src = ""
    if (("HAP1cur", $2) in byLen) src = "HAP1cur"
    if (("HAP2cur", $2) in byLen) src = (src == "") ? "HAP2cur" : "AMBIGUO"
    if (src == "" || src == "AMBIGUO") {
        printf("ERROR: %s %s (%s bp) -> %s\n", hap, $1, $2,
               (src == "") ? "sin match por longitud" : "match en AMBOS curados") > "/dev/stderr"
        errs++
        next
    }

    print hap, "SINGLE", $1, rename(hap, $1, src)      # el SUPER, tal cual del single
    parent = byLen[src, $2]
    for (n = 1; n <= ulmax[src, parent]; n++)
        if ((src, parent, n) in ulname) {
            u = ulname[src, parent, n]
            print hap, src, u, rename(hap, u, src)     # sus unlocs, del curado que le toca
            nunloc[hap]++
        }
    nsuper[hap]++
}

END {
    printf("mat: %d SUPER + %d unlocs | pat: %d SUPER + %d unlocs | Scaffold_* descartados: %d | errores: %d\n",
           nsuper["mat"], nunloc["mat"], nsuper["pat"], nunloc["pat"], skipped+0, errs+0) > "/dev/stderr"
    if (errs) exit 1
}
