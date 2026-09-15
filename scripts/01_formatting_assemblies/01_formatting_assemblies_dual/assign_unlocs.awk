# DUAL: asigna a cada secuencia su haplotipo (Mat_/Pat_) y su curado de origen (.H1/.H2).
#
#   Mat_/Pat_  <- del archivo de entrada (dual_mat.fasta / dual_pat.fasta)
#   .H1/.H2    <- DETERMINADO por longitud contra los curados, no leido del nombre.
#                 Cualquier sufijo .H que ya venga en el header se descarta y se
#                 reemplaza por el que dicta la longitud; si no coinciden se avisa.
#
# Los unlocs heredan el haplotipo de su SUPER padre y salen del mismo curado.
# Los Scaffold_* se descartan.
#
# Uso:
#   awk -F'\t' -v OFS='\t' -f assign_unlocs.awk \
#       H1cur.fai H2cur.fai dual_mat.fasta.fai dual_pat.fasta.fai
#
# Salida (stdout), en el orden final del FASTA:
#   hap <TAB> fuente <TAB> nombre_original <TAB> nombre_nuevo
#   fuente = DUAL | HAP1cur | HAP2cur   (de que FASTA extraer la secuencia)

function ver(f,   b) {
    b = f; sub(/.*\//, "", b)
    if (b ~ /hap1\.cur/) return "HAP1cur"
    if (b ~ /hap2\.cur/) return "HAP2cur"
    return "?"
}
function stripH(name,   s) { s = name; sub(/\.H[0-9]+$/, "", s); return s }
function tag(src)          { return (src == "HAP1cur") ? ".H1" : ".H2" }
# SUPER_Z_unloc_20.H1 -> SUPER_Z.H1  (sin anclar: el .H va DESPUES del indice)
function parent_of(name,   p) { p = name; sub(/_unloc_[0-9]+/, "", p); return p }
function idx_of(name,   n)    { n = name; sub(/.*_unloc_/, "", n); sub(/\..*/, "", n); return n + 0 }

# ---- 1) curados: indice longitud->nombre, y unlocs por padre ----
FILENAME !~ /dual_(mat|pat)/ {
    if ($1 ~ /^Scaffold_/) next
    v = ver(FILENAME)
    if (($2) in seenlen && vlen[$2] != v) ambig[$2] = 1     # misma longitud en ambos curados
    seenlen[$2] = 1; vlen[$2] = v
    byLen[v, $2] = $1
    if ($1 ~ /_unloc_[0-9]+/) {
        p = parent_of($1); n = idx_of($1)
        ulname[v, p, n] = $1
        if (n > ulmax[v, p]) ulmax[v, p] = n
    }
    next
}

# ---- 2) dual_mat / dual_pat ----
{
    if ($1 ~ /^Scaffold_/) { skipped++; next }
    hap    = (FILENAME ~ /dual_mat/) ? "mat" : "pat"
    prefix = (hap == "mat") ? "Mat_" : "Pat_"

    # --- determinar el curado ---
    # El sufijo del nombre es la clave: en este ensamblaje dos secuencias
    # distintas comparten longitud (SUPER_35.H1 y SUPER_33.H2, 2050511 bp), asi
    # que la longitud sola es ambigua. La longitud se usa para VERIFICAR.
    src = ""
    if ($1 ~ /\.H1$/) src = "HAP1cur"
    else if ($1 ~ /\.H2$/) src = "HAP2cur"

    if (src == "") {                       # sin sufijo: caer a la longitud
        if (("HAP1cur", $2) in byLen) src = "HAP1cur"
        if (("HAP2cur", $2) in byLen) src = (src == "") ? "HAP2cur" : "AMBIGUO"
        if (src == "" || src == "AMBIGUO") {
            printf("ERROR: %s %s (%s bp): sin sufijo .H y %s\n", hap, $1, $2,
                   (src == "") ? "sin match por longitud" : "longitud en AMBOS curados") > "/dev/stderr"
            errs++; next
        }
    } else if (!((src, $2) in byLen)) {     # sufijo presente: verificar longitud
        printf("ERROR: %s %s (%s bp): el sufijo dice %s pero no hay secuencia de esa longitud ahi\n",
               hap, $1, $2, tag(src)) > "/dev/stderr"
        errs++; next
    }

    base = stripH($1)
    print hap, "DUAL", $1, prefix base tag(src)            # el SUPER

    parent = byLen[src, $2]                                # su nombre dentro del curado
    for (n = 1; n <= ulmax[src, parent]; n++)
        if ((src, parent, n) in ulname) {
            u = ulname[src, parent, n]
            print hap, src, u, prefix stripH(u) tag(src)   # sus unlocs
            nunloc[hap]++
        }
    nsuper[hap]++
}

END {
    printf("mat: %d SUPER + %d unlocs | pat: %d SUPER + %d unlocs | Scaffold_* descartados: %d | avisos: %d | errores: %d\n",
           nsuper["mat"], nunloc["mat"], nsuper["pat"], nunloc["pat"], skipped+0, warns+0, errs+0) > "/dev/stderr"
    if (errs) exit 1
}
