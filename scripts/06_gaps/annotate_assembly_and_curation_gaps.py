#!/usr/bin/env python3
"""
annotate_assembly_and_curation_gaps.py

Classifies physical gaps (obtained via gfastats -bg on the final reoriented FASTA)
into:
  - CURATION: joins introduced during manual curation (lines 'U' in AGP)
  - ASSEMBLY: original gaps from the initial assembler contigs

Uses first-principles orientation detection:
For each chromosome with length L (from .fai), curation joins from the AGP
are tested against physical gaps in both:
  - Forward (+) orientation: [beg - 1, end]
  - Inverted (-) orientation: [L - end, L - (beg - 1)]

This eliminates fragile dependencies on intermediate rvcp.sak files.
"""

import os
import sys
import argparse
from collections import defaultdict


def load_fai_lengths(fai_path):
    """Load sequence lengths from .fai index."""
    lengths = {}
    if not fai_path or not os.path.exists(fai_path):
        raise FileNotFoundError(f"FAI file not found: {fai_path}")
    with open(fai_path) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                lengths[parts[0]] = int(parts[1])
    return lengths


def load_scaffold_maps(map_path, mode="dual"):
    """
    Loads scaffold mapping from scaffold_hap_to_combined_super.tsv.
    Format:
      combined_SUPER  mat_pat  haplotype  scaffold  SUPER_internal  SUPER_renamed

    Returns:
      For single mode: (h1_map, h2_map) where key is scaffold (e.g. Scaffold_1)
      For dual/ont mode: (scaf_map, None)
    """
    if not map_path or not os.path.exists(map_path):
        raise FileNotFoundError(f"Scaffold map file not found: {map_path}")

    h1_map = {}
    h2_map = {}
    unified_map = {}

    with open(map_path) as fh:
        header = None
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if header is None:
                header = parts
                continue
            combined = parts[0]
            hap = parts[2]
            scaf = parts[3]

            if hap == 'H1':
                h1_map[scaf] = combined
            elif hap == 'H2':
                h2_map[scaf] = combined

            unified_map[scaf] = combined

    if mode.lower() == "single":
        return h1_map, h2_map
    return unified_map, None


def parse_agp_curation_gaps(agp_file, scaf_map, fai_lens):
    """
    Parses an AGP file and extracts non-unloc curation gaps (component_type == 'U').
    Excludes gaps adjacent to Unloc contigs and gaps beyond sequence length in FAI.
    Returns list of dicts:
      [{'chrom': final_chrom, 'beg': beg_0, 'end': end_0, 'scaf': scaf, 'agp': agp_file}, ...]
    """
    cur_gaps = []
    if not agp_file or not os.path.exists(agp_file):
        raise FileNotFoundError(f"AGP file not found: {agp_file}")

    with open(agp_file) as fh:
        lines = [l.strip().split('\t') for l in fh if l.strip() and not l.startswith('#')]

    for i, p in enumerate(lines):
        if len(p) >= 5 and p[4] == 'U':
            scaf = p[0]
            # Check surrounding contigs for 'Unloc'
            prev_c = lines[i - 1][5] if i > 0 and len(lines[i - 1]) > 5 else ''
            next_c = lines[i + 1][5] if i + 1 < len(lines) and len(lines[i + 1]) > 5 else ''
            if 'unloc' in prev_c.lower() or 'unloc' in next_c.lower():
                continue

            if scaf in scaf_map:
                final_chrom = scaf_map[scaf]
                beg_0 = int(p[1]) - 1
                end_0 = int(p[2])
                # Only keep if within chromosome length in FAI
                if final_chrom in fai_lens and end_0 <= fai_lens[final_chrom]:
                    cur_gaps.append({
                        'chrom': final_chrom,
                        'beg': beg_0,
                        'end': end_0,
                        'scaf': scaf,
                        'agp': os.path.basename(agp_file)
                    })
    return cur_gaps


def load_physical_gaps(bed_path):
    """
    Loads physical gaps from gfastats -bg output.
    Returns:
      ordered_gaps: list of (chrom, start, end)
      gaps_by_chrom: dict of {chrom: list of (start, end)}
    """
    ordered_gaps = []
    gaps_by_chrom = defaultdict(list)
    if not bed_path or not os.path.exists(bed_path):
        raise FileNotFoundError(f"Gaps BED not found: {bed_path}")

    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            chrom = parts[0]
            s = int(parts[1])
            e = int(parts[2])
            ordered_gaps.append((chrom, s, e))
            gaps_by_chrom[chrom].append((s, e))
    return ordered_gaps, gaps_by_chrom


def main():
    parser = argparse.ArgumentParser(
        description="Annotate gfastats physical gaps into ASSEMBLY and CURATION using first-principles orientation."
    )
    parser.add_argument("--mode", choices=["dual", "single", "ont"], required=True,
                        help="Assembly mode: dual, single, or ont")
    parser.add_argument("--gaps-bed", required=True,
                        help="Physical gaps BED from gfastats on final reoriented FASTA")
    parser.add_argument("--fai", required=True,
                        help="FASTA index (.fai) for sequence lengths")
    parser.add_argument("--scaffold-map", required=True,
                        help="scaffold_hap_to_combined_super.tsv file")
    parser.add_argument("--agp", nargs="*",
                        help="AGP file(s) for dual or ont mode")
    parser.add_argument("--hap1-agp",
                        help="Hap1 AGP file (required for single mode)")
    parser.add_argument("--hap2-agp",
                        help="Hap2 AGP file (required for single mode)")
    parser.add_argument("--tolerance", type=int, default=50,
                        help="Midpoint matching tolerance in bp (default: 50)")
    parser.add_argument("--output-annotated", required=True,
                        help="Output BED: chrom, start, end, label (ASSEMBLY or CURATION)")
    parser.add_argument("--output-curation-bed",
                        help="Optional output BED containing only CURATION gaps")
    parser.add_argument("--output-summary", required=True,
                        help="Output summary text report")

    args = parser.parse_args()

    # 1. Load sequence lengths
    fai_lens = load_fai_lengths(args.fai)

    # 2. Load scaffold maps
    map1, map2 = load_scaffold_maps(args.scaffold_map, mode=args.mode)

    # 3. Parse AGP curation gaps
    candidate_gaps = []
    if args.mode == "single":
        if not args.hap1_agp or not args.hap2_agp:
            parser.error("--mode single requires both --hap1-agp and --hap2-agp")
        candidate_gaps.extend(parse_agp_curation_gaps(args.hap1_agp, map1, fai_lens))
        candidate_gaps.extend(parse_agp_curation_gaps(args.hap2_agp, map2, fai_lens))
    else:
        if not args.agp:
            parser.error(f"--mode {args.mode} requires --agp")
        for agp_f in args.agp:
            candidate_gaps.extend(parse_agp_curation_gaps(agp_f, map1, fai_lens))

    print(f"[{args.mode.upper()}] Parsed {len(candidate_gaps)} candidate curation joins from AGP.", file=sys.stderr)

    # 4. Load physical gaps from gfastats
    ordered_phys, phys_by_chrom = load_physical_gaps(args.gaps_bed)
    print(f"[{args.mode.upper()}] Loaded {len(ordered_phys)} physical gaps from gfastats.", file=sys.stderr)

    # 5. Match candidate curation gaps against physical gaps (testing + and - orientations)
    matched_phys_keys = set()  # set of (chrom, s, e) identified as CURATION
    chrom_orientations = defaultdict(set)
    unmatched_candidates = []

    for cgap in candidate_gaps:
        chrom = cgap['chrom']
        s_agp = cgap['beg']
        e_agp = cgap['end']
        scaf = cgap['scaf']
        L = fai_lens[chrom]

        # Forward (+) coordinates
        mid_dir = (s_agp + e_agp) / 2.0
        # Inverted (-) coordinates
        s_inv = L - e_agp
        e_inv = L - s_agp
        mid_inv = (s_inv + e_inv) / 2.0

        best_diff = float('inf')
        best_ori = None
        best_gap = None

        for ps, pe in phys_by_chrom.get(chrom, []):
            pmid = (ps + pe) / 2.0
            d_dir = abs(pmid - mid_dir)
            d_inv = abs(pmid - mid_inv)

            if d_dir < best_diff:
                best_diff = d_dir
                best_ori = '+'
                best_gap = (chrom, ps, pe)
            if d_inv < best_diff:
                best_diff = d_inv
                best_ori = '-'
                best_gap = (chrom, ps, pe)

        if best_diff <= args.tolerance and best_gap is not None:
            matched_phys_keys.add(best_gap)
            chrom_orientations[chrom].add(best_ori)
        else:
            unmatched_candidates.append(cgap)

    # 6. Annotate all physical gaps
    annotated_rows = []
    n_curation = 0
    n_assembly = 0
    curation_rows = []

    for chrom, s, e in ordered_phys:
        is_cur = (chrom, s, e) in matched_phys_keys
        label = "CURATION" if is_cur else "ASSEMBLY"
        if is_cur:
            n_curation += 1
            curation_rows.append((chrom, s, e, label))
        else:
            n_assembly += 1
        annotated_rows.append((chrom, s, e, label))

    # 7. Write outputs
    os.makedirs(os.path.dirname(os.path.abspath(args.output_annotated)), exist_ok=True)
    with open(args.output_annotated, 'w') as fh:
        for chrom, s, e, label in annotated_rows:
            fh.write(f"{chrom}\t{s}\t{e}\t{label}\n")
    print(f"[{args.mode.upper()}] Wrote {len(annotated_rows)} annotated gaps to {args.output_annotated}", file=sys.stderr)

    if args.output_curation_bed:
        os.makedirs(os.path.dirname(os.path.abspath(args.output_curation_bed)), exist_ok=True)
        with open(args.output_curation_bed, 'w') as fh:
            for chrom, s, e, label in curation_rows:
                fh.write(f"{chrom}\t{s}\t{e}\t{label}\n")
        print(f"[{args.mode.upper()}] Wrote {len(curation_rows)} curation gaps to {args.output_curation_bed}", file=sys.stderr)

    # 8. Generate Summary Report
    total_phys = len(ordered_phys)
    pct_cur = (n_curation / total_phys * 100.0) if total_phys > 0 else 0.0
    pct_asm = (n_assembly / total_phys * 100.0) if total_phys > 0 else 0.0

    lines = [
        "================================================================================",
        f" Gap Annotation Summary: {args.mode.upper()}",
        "================================================================================",
        f"Total Physical Gaps (gfastats -bg) : {total_phys}",
        f"  - Curation Gaps (joins)          : {n_curation} ({pct_cur:.1f}%)",
        f"  - Assembly Gaps (contigs)        : {n_assembly} ({pct_asm:.1f}%)",
        f"Candidate AGP Curation Joins       : {len(candidate_gaps)}",
        f"  - Matched in Physical Gaps       : {len(matched_phys_keys)}",
        f"  - Unmatched Candidates           : {len(unmatched_candidates)}",
        "--------------------------------------------------------------------------------",
        " Inferred Chromosome Orientations (from Curation Gaps):",
        "--------------------------------------------------------------------------------",
    ]

    for chrom in sorted(chrom_orientations.keys()):
        oris = "/".join(sorted(chrom_orientations[chrom]))
        ori_label = "FORWARD (+)" if oris == "+" else ("INVERTED (-)" if oris == "-" else f"MIXED ({oris})")
        lines.append(f"  {chrom:<25} : {ori_label}")

    if unmatched_candidates:
        lines.append("--------------------------------------------------------------------------------")
        lines.append(" Unmatched AGP Candidates (scaffolds without corresponding physical gaps):")
        lines.append("--------------------------------------------------------------------------------")
        for u in unmatched_candidates:
            lines.append(f"  {u['chrom']} ({u['scaf']}) [{u['agp']}] : {u['beg']}-{u['end']}")

    lines.append("================================================================================\n")
    summary_text = "\n".join(lines)

    os.makedirs(os.path.dirname(os.path.abspath(args.output_summary)), exist_ok=True)
    with open(args.output_summary, 'w') as fh:
        fh.write(summary_text)

    print(summary_text, file=sys.stderr)


if __name__ == "__main__":
    main()
