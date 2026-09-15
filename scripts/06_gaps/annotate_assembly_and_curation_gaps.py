#!/usr/bin/env python3
"""
annotate_assembly_and_curation_gaps.py

Classifies physical gaps (obtained via gfastats on the final FASTA) into:
  - CURATION: gaps introduced during manual curation / contig joining (lines 'U' in AGP)
  - ASSEMBLY: original gaps from the initial assembler contigs

Transforms curation AGP gap coordinates from old/pre-reorientation coordinates
to final assembly coordinates using the scaffold mapping, sequence lengths (from .fai),
and the reorientation list (rvcp.sak).

Outputs standard format:
  scaffold   start   end   label (ASSEMBLY or CURATION)
"""

import os
import sys
import argparse
from collections import defaultdict


def load_fai_lengths(fai_path):
    """Load sequence lengths from .fai index."""
    lengths = {}
    if not fai_path or not os.path.exists(fai_path):
        return lengths
    with open(fai_path) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                lengths[parts[0]] = int(parts[1])
    return lengths


def load_scaffold_map(map_paths):
    """
    Parses scaffold mapping TSV file(s).
    Supports formats:
      1) combined_SUPER, mat_pat, haplotype, scaffold, SUPER_internal, SUPER_renamed
      2) scaffold, target_name
      3) old_name -> new_name
    Returns dict: {scaffold_old: scaffold_new}
    """
    scaf_map = {}
    if not map_paths:
        return scaf_map

    for p in map_paths:
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                # Check format 1 (Dual / Single table)
                if parts[0] == 'combined_SUPER':
                    continue
                if len(parts) >= 4 and parts[3].startswith('Scaffold_'):
                    comb = parts[0]
                    matpat = parts[1]
                    scaf = parts[3]
                    prefix = 'Mat_' if matpat.lower() == 'mat' else 'Pat_'
                    base = comb.replace('.mat', '').replace('.pat', '')
                    scaf_map[scaf] = prefix + base
                elif len(parts) >= 2:
                    scaf_map[parts[0]] = parts[1]
                    scaf_map[parts[1]] = parts[1]  # self-map
    return scaf_map


def load_rvcp_supers(rvcp_path):
    """
    Load set of supers that were reoriented (reverse-complemented) from rvcp.sak.
    Lines typically: RVCP SUPER_2 or just SUPER_2.
    """
    rvcp = set()
    if not rvcp_path or not os.path.exists(rvcp_path):
        return rvcp
    with open(rvcp_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) == 1:
                rvcp.add(parts[0])
            elif len(parts) >= 2:
                rvcp.add(parts[1])
    return rvcp


def parse_agp_curation_gaps(agp_paths, scaf_map, rvcp_supers, fai_lens):
    """
    Extracts 'U' lines from AGP file(s) and converts them to 0-based BED coordinates
    in the final assembly space.
    Returns dict: {scaffold: list of (start, end)}
    """
    cur_gaps = defaultdict(list)
    if not agp_paths:
        return cur_gaps

    for agp_path in agp_paths:
        if not os.path.exists(agp_path):
            print(f"Warning: AGP file not found: {agp_path}", file=sys.stderr)
            continue

        with open(agp_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                p = line.split('\t')
                if len(p) < 5:
                    continue

                # Check if it's a gap record (component_type == 'U')
                if p[4] != 'U':
                    continue

                scaf = p[0]
                beg = int(p[1])  # 1-based
                end = int(p[2])  # 1-based

                # Determine final scaffold name
                final_scaf = scaf_map.get(scaf, scaf)

                # Check if this scaffold was reoriented
                # Match against base super name (e.g. SUPER_2)
                super_base = final_scaf.replace('Mat_', '').replace('Pat_', '').split('.')[0]
                # In single/dual, typically H2 was reoriented with rvcp
                is_rv = False
                if super_base in rvcp_supers:
                    if '.H2' in final_scaf or '.pat' in final_scaf or '.mat' in final_scaf or not any('.H' in s for s in fai_lens):
                        is_rv = True

                if is_rv and final_scaf in fai_lens:
                    L = fai_lens[final_scaf]
                    # Invert coordinates: 0-based start = L - end, end = L - (beg - 1)
                    s_new = L - end
                    e_new = L - (beg - 1)
                    if s_new >= 0 and e_new <= L:
                        cur_gaps[final_scaf].append((min(s_new, e_new), max(s_new, e_new)))
                else:
                    s_new = beg - 1
                    e_new = end
                    cur_gaps[final_scaf].append((s_new, e_new))

    return cur_gaps


def load_curation_bed(bed_paths):
    """Load pre-projected curation gaps BED directly if provided."""
    cur_gaps = defaultdict(list)
    if not bed_paths:
        return cur_gaps
    for bp in bed_paths:
        if not os.path.exists(bp):
            continue
        with open(bp) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                p = line.split('\t')
                chrom = p[0]
                s = int(p[1])
                e = int(p[2])
                cur_gaps[chrom].append((s, e))
    return cur_gaps


def main():
    parser = argparse.ArgumentParser(description="Annotate gfastats physical gaps into ASSEMBLY and CURATION.")
    parser.add_argument("--gaps-bed", required=True, help="Input BED with all physical gaps from gfastats")
    parser.add_argument("--agp", nargs='*', help="One or more curation AGP files")
    parser.add_argument("--curation-bed", nargs='*', help="Alternative pre-projected curation gaps BED")
    parser.add_argument("--scaffold-map", nargs='*', help="TSV mapping AGP scaffold names to final names")
    parser.add_argument("--rvcp-sak", help="File listing reoriented supers (rvcp.sak)")
    parser.add_argument("--fai", help="FASTA index (.fai) for sequence lengths")
    parser.add_argument("--output-annotated", required=True, help="Output annotated gaps BED (scaffold, start, end, label)")
    parser.add_argument("--output-summary", required=True, help="Output summary text file")
    parser.add_argument("--asm-name", default="Assembly", help="Assembly label (Single, Dual, ONT)")
    parser.add_argument("--tolerance", type=int, default=10, help="Coordinate matching tolerance in bp (default: 10)")

    args = parser.parse_args()

    # 1. Load sequence lengths
    fai_lens = load_fai_lengths(args.fai)

    # 2. Load scaffold mapping
    scaf_map = load_scaffold_map(args.scaffold_map)

    # 3. Load reorientation list
    rvcp_supers = load_rvcp_supers(args.rvcp_sak)

    # 4. Extract curation gaps
    cur_gaps = defaultdict(list)
    if args.curation_bed:
        cur_bed_gaps = load_curation_bed(args.curation_bed)
        for k, v in cur_bed_gaps.items():
            cur_gaps[k].extend(v)

    if args.agp:
        agp_gaps = parse_agp_curation_gaps(args.agp, scaf_map, rvcp_supers, fai_lens)
        for k, v in agp_gaps.items():
            cur_gaps[k].extend(v)

    total_cur_targets = sum(len(v) for v in cur_gaps.values())
    print(f"[{args.asm_name}] Loaded {total_cur_targets} candidate curation gaps from AGP/BED.", file=sys.stderr)

    # 5. Read physical gaps from gfastats
    annotated_rows = []
    n_curation = 0
    n_assembly = 0
    tol = args.tolerance

    with open(args.gaps_bed) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            scaf = parts[0]
            s = int(parts[1])
            e = int(parts[2])

            is_curation = False
            for cs, ce in cur_gaps.get(scaf, []):
                # Overlap check or close proximity check
                if (max(s, cs) < min(e, ce)) or (abs(s - cs) <= tol and abs(e - ce) <= tol):
                    is_curation = True
                    break

            label = "CURATION" if is_curation else "ASSEMBLY"
            if is_curation:
                n_curation += 1
            else:
                n_assembly += 1

            annotated_rows.append((scaf, s, e, label))

    # 6. Write output annotated gaps file
    os.makedirs(os.path.dirname(os.path.abspath(args.output_annotated)), exist_ok=True)
    with open(args.output_annotated, 'w') as fh:
        for scaf, s, e, label in annotated_rows:
            fh.write(f"{scaf}\t{s}\t{e}\t{label}\n")

    print(f"[{args.asm_name}] Wrote {len(annotated_rows)} annotated gaps to {args.output_annotated}", file=sys.stderr)

    # 7. Write output summary
    os.makedirs(os.path.dirname(os.path.abspath(args.output_summary)), exist_ok=True)
    n_total = len(annotated_rows)
    pct_cur = (n_curation / n_total * 100.0) if n_total > 0 else 0.0
    pct_asm = (n_assembly / n_total * 100.0) if n_total > 0 else 0.0

    summary_text = (
        f"===========================================================\n"
        f" Gap Annotation Summary: {args.asm_name}\n"
        f"===========================================================\n"
        f"Total Physical Gaps (gfastats) : {n_total}\n"
        f"  - Curation Gaps (joins)      : {n_curation} ({pct_cur:.1f}%)\n"
        f"  - Assembly Gaps (contigs)    : {n_assembly} ({pct_asm:.1f}%)\n"
        f"Candidate Curation Joins (AGP) : {total_cur_targets}\n"
        f"===========================================================\n"
    )

    with open(args.output_summary, 'w') as fh:
        fh.write(summary_text)

    print(summary_text, file=sys.stderr)


if __name__ == "__main__":
    main()
