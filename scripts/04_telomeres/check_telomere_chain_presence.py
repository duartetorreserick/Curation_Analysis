#!/usr/bin/env python3
"""
check_telomere_chain_presence.py

Evaluates the presence and absence of telomeres for each chromosome across
collinear and non-collinear chain alignments, with Teloscope-based rescue for
ends where the chain did not reach the physical chromosome terminus.

Output columns:
  chromosome            T2T reference chromosome
  scaffold              Paired assembly scaffold from best_chrom_pairs.tsv
  collinear             Telomeres covered collinearly (pq, p, q, none)
  non-collinear         Telomeres covered non-collinearly (pq, p, q, none)
  teloscope_rescued     Whether Teloscope rescued any telomere (yes, no)
  type_of_tele_rescued  Which arm was rescued by Teloscope (p, q, pq, none)

Usage:
  python3 check_telomere_chain_presence.py \
      --t2t-bed data/t2t/bTaeGut7.T2T.fasta_terminal_telomeres.bed \
      --pairs-file results/single/02_chainpipeline_single/t2t.vs.single.best_chrom_pairs.tsv \
      --collinear-chain results/single/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.target.collinear.chain \
      --noncollinear-chain results/single/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.target.non-collinear.chain \
      --assembly-telo-bed results/single/04_telomeres_single/single_combined.renamed.fasta_terminal_telomeres.bed \
      --output-tsv results/single/04_telomeres_single/single_telomere_presence.tsv \
      --output-summary results/single/04_telomeres_single/single_telomere_summary.txt \
      --asm-name Single
"""

import os
import sys
import argparse
from collections import defaultdict


def parse_t2t_bed(bed_path):
    """
    Parses T2T terminal telomeres BED file.
    Returns:
      chrom_order: list of chromosome names preserving original order
      t2t_by_chrom: dict of {chrom: {arm: (start, end)}}
    """
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"T2T BED file not found: {bed_path}")

    chrom_order = []
    t2t_by_chrom = defaultdict(dict)

    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 5:
                continue
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            arm = cols[4]

            if chrom not in t2t_by_chrom:
                chrom_order.append(chrom)
            t2t_by_chrom[chrom][arm] = (start, end)

    return chrom_order, t2t_by_chrom


def parse_pairs_file(pairs_path):
    """
    Parses best_chrom_pairs.tsv file mapping T2T chromosome -> assembly scaffold.
    Returns dict: {t2t_chrom: assembly_scaffold}
    """
    pairs = {}
    if not pairs_path or not os.path.exists(pairs_path):
        return pairs

    with open(pairs_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) >= 2:
                pairs[cols[0]] = cols[1]
    return pairs


def parse_assembly_telo_bed(bed_path):
    """
    Parses assembly terminal telomeres BED file (Teloscope / FindTelomeres).
    Returns dict: {scaffold: {'p': True/False, 'q': True/False}}
    """
    telos = defaultdict(lambda: {'p': False, 'q': False})
    if not bed_path or not os.path.exists(bed_path):
        return telos

    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) >= 5:
                scaf = cols[0]
                arm = cols[4]
                if arm in ('p', 'q'):
                    telos[scaf][arm] = True
    return telos


def parse_chain_intervals(chain_path):
    """
    Parses target alignment intervals from a UCSC chain file.
    Returns dict: {t_chrom: list of (t_start, t_end)}
    """
    intervals = defaultdict(list)
    if not chain_path or not os.path.exists(chain_path):
        return intervals

    with open(chain_path) as fh:
        t_name = None
        pos_t = 0

        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                if not line:
                    t_name = None
                continue

            if line.startswith('chain'):
                parts = line.split()
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                t_name = parts[2]
                pos_t = int(parts[5])
                continue

            if t_name is None:
                continue

            parts = line.split()
            if not parts:
                continue

            size = int(parts[0])
            intervals[t_name].append((pos_t, pos_t + size))

            if len(parts) >= 2:
                dt = int(parts[1])
                pos_t += size + dt
            else:
                pos_t += size

    # Sort intervals for each target chromosome
    for chrom in intervals:
        intervals[chrom].sort(key=lambda x: x[0])

    return intervals


def check_overlap(chrom, start, end, chain_ivs):
    """
    Checks if interval [start, end) on chrom overlaps any chain interval.
    """
    for a, b in chain_ivs.get(chrom, []):
        if a < end and b > start:
            return True
    return False


def evaluate_assembly(chrom_order, t2t_by_chrom, pairs, coll_ivs, nc_ivs, asm_telos):
    """
    Evaluates each chromosome for collinear, non-collinear, and rescued telomeres.
    Returns list of row dicts.
    """
    rows = []

    for chrom in chrom_order:
        arms_dict = t2t_by_chrom[chrom]
        scaf = pairs.get(chrom, "unassigned")

        coll_arms = []
        nc_arms = []
        rescued_arms = []

        for arm in ['p', 'q']:
            if arm not in arms_dict:
                continue
            start, end = arms_dict[arm]
            c_ok = check_overlap(chrom, start, end, coll_ivs)
            nc_ok = check_overlap(chrom, start, end, nc_ivs)

            if c_ok:
                coll_arms.append(arm)
            elif nc_ok:
                nc_arms.append(arm)
            else:
                # Missing in chain alignments: check if paired scaffold has same arm in Teloscope
                if scaf != "unassigned" and asm_telos[scaf][arm]:
                    rescued_arms.append(arm)

        coll_str = ''.join(coll_arms) if coll_arms else 'none'
        nc_str = ''.join(nc_arms) if nc_arms else 'none'
        rescued_str = 'yes' if rescued_arms else 'no'
        type_rescued_str = ''.join(rescued_arms) if rescued_arms else 'none'

        rows.append({
            'chromosome': chrom,
            'scaffold': scaf,
            'collinear': coll_str,
            'non-collinear': nc_str,
            'teloscope_rescued': rescued_str,
            'type_of_tele_rescued': type_rescued_str,
            'coll_arms': coll_arms,
            'nc_arms': nc_arms,
            'rescued_arms': rescued_arms
        })

    return rows


def write_output_tsv(rows, out_path):
    """Writes results to TSV file with exact requested columns."""
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    headers = [
        "chromosome",
        "scaffold",
        "collinear",
        "non-collinear",
        "teloscope_rescued",
        "type_of_tele_rescued"
    ]
    with open(out_path, 'w') as fh:
        fh.write('\t'.join(headers) + '\n')
        for r in rows:
            fh.write(f"{r['chromosome']}\t{r['scaffold']}\t{r['collinear']}\t{r['non-collinear']}\t{r['teloscope_rescued']}\t{r['type_of_tele_rescued']}\n")


def write_summary(rows, out_path, asm_name="Assembly"):
    """Writes a human-readable summary report with statistics."""
    total_chroms = len(rows)
    total_telos = sum(len(r['coll_arms']) + len(r['nc_arms']) + len(r['rescued_arms']) for r in rows)

    coll_telo_count = sum(len(r['coll_arms']) for r in rows)
    nc_telo_count = sum(len(r['nc_arms']) for r in rows)
    rescued_telo_count = sum(len(r['rescued_arms']) for r in rows)
    expected_telos = total_chroms * 2  # 80 * 2 = 160

    # Count complete (2-telomere) chromosomes
    complete_before = sum(1 for r in rows if len(set(r['coll_arms'] + r['nc_arms'])) == 2)
    complete_after = sum(1 for r in rows if len(set(r['coll_arms'] + r['nc_arms'] + r['rescued_arms'])) == 2)

    lines = [
        f"===========================================================",
        f" Telomere Analysis & Teloscope Rescue Summary: {asm_name}",
        f"===========================================================",
        f"Total T2T Reference Chromosomes : {total_chroms}",
        f"Total Reference Telomere Ends   : {expected_telos} (80 p + 80 q)",
        f"-----------------------------------------------------------",
        f"Telomeres Covered in Chains:",
        f"  - Collinear                   : {coll_telo_count} ({coll_telo_count/expected_telos*100:.1f}%)",
        f"  - Non-collinear               : {nc_telo_count} ({nc_telo_count/expected_telos*100:.1f}%)",
        f"  - Subtotal in Chains          : {coll_telo_count + nc_telo_count} ({(coll_telo_count + nc_telo_count)/expected_telos*100:.1f}%)",
        f"-----------------------------------------------------------",
        f"Teloscope Rescue:",
        f"  - Telomeres Rescued           : {rescued_telo_count}",
        f"  - Total Telomeres Accounted   : {total_telos} ({total_telos/expected_telos*100:.1f}%)",
        f"  - Missing Telomeres           : {expected_telos - total_telos} ({(expected_telos - total_telos)/expected_telos*100:.1f}%)",
        f"-----------------------------------------------------------",
        f"Complete (2-Telomere) Chromosomes:",
        f"  - Before Teloscope Rescue     : {complete_before} / {total_chroms} ({complete_before/total_chroms*100:.1f}%)",
        f"  - AFTER Teloscope Rescue      : {complete_after} / {total_chroms} ({complete_after/total_chroms*100:.1f}%)",
        f"  - Net Gain                    : +{complete_after - complete_before} complete chromosomes",
        f"==========================================================="
    ]

    rescued_list = [
        f"  {r['chromosome']} ({r['scaffold']}): rescued arm(s) [{r['type_of_tele_rescued']}], collinear [{r['collinear']}], non-collinear [{r['non-collinear']}]"
        for r in rows if r['teloscope_rescued'] == 'yes'
    ]
    if rescued_list:
        lines.append(f"\nRescued Chromosomes Details ({len(rescued_list)}):")
        lines.extend(rescued_list)

    summary_text = '\n'.join(lines) + '\n'

    if out_path:
        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        with open(out_path, 'w') as fh:
            fh.write(summary_text)

    return summary_text


def main():
    parser = argparse.ArgumentParser(description="Check telomere presence/absence across chains with Teloscope rescue.")
    parser.add_argument("--t2t-bed", required=True, help="T2T reference terminal telomeres BED")
    parser.add_argument("--pairs-file", required=True, help="best_chrom_pairs.tsv mapping T2T to assembly scaffolds")
    parser.add_argument("--collinear-chain", required=True, help="Target collinear chain file")
    parser.add_argument("--noncollinear-chain", required=True, help="Target non-collinear chain file")
    parser.add_argument("--assembly-telo-bed", required=True, help="Assembly terminal telomeres BED from Teloscope")
    parser.add_argument("--output-tsv", required=True, help="Path for output TSV")
    parser.add_argument("--output-summary", help="Optional path for output summary TXT")
    parser.add_argument("--asm-name", default="Assembly", help="Assembly label (e.g. Single, Dual, ONT)")

    args = parser.parse_args()

    print(f"Loading T2T reference telomeres from {args.t2t_bed}...", file=sys.stderr)
    chrom_order, t2t_by_chrom = parse_t2t_bed(args.t2t_bed)

    print(f"Loading chromosome pairs from {args.pairs_file}...", file=sys.stderr)
    pairs = parse_pairs_file(args.pairs_file)

    print(f"Loading assembly telomeres from {args.assembly_telo_bed}...", file=sys.stderr)
    asm_telos = parse_assembly_telo_bed(args.assembly_telo_bed)

    print(f"Parsing collinear chain intervals from {args.collinear_chain}...", file=sys.stderr)
    coll_ivs = parse_chain_intervals(args.collinear_chain)

    print(f"Parsing non-collinear chain intervals from {args.noncollinear_chain}...", file=sys.stderr)
    nc_ivs = parse_chain_intervals(args.noncollinear_chain)

    print("Evaluating telomere presence and Teloscope rescue...", file=sys.stderr)
    rows = evaluate_assembly(chrom_order, t2t_by_chrom, pairs, coll_ivs, nc_ivs, asm_telos)

    write_output_tsv(rows, args.output_tsv)
    print(f"Wrote output TSV: {args.output_tsv}", file=sys.stderr)

    summary_text = write_summary(rows, args.output_summary, args.asm_name)
    print(summary_text)


if __name__ == "__main__":
    main()
