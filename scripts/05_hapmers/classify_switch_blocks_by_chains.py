#!/usr/bin/env python3
"""
classify_switch_blocks_by_chains.py

Classifies switch blocks (hapmers) into collinear, non-collinear, unaligned,
or mixed categories by intersecting them directly with UCSC alignment chains.

Labels (in output BED):
  - collinear:     >= 90% of the block falls within the collinear chain
  - non_collinear: >= 90% of the block falls within the non-collinear chain
  - unaligned:     >= 90% of the block falls outside both chains (gaps/divergent)
  - mixed:         block spans multiple states without any single state reaching 90%

The script reads tName (T2T reference) and qName (assembly scaffold) directly from
the chain headers, requiring no external pairing file.

Usage:
  python3 classify_switch_blocks_by_chains.py \
      --switch-bed results/single/05_hapmers_single/single.single_combined.100_20000.phased_block.renamed.switch_blocks.filtered.bed \
      --collinear-chain results/single/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.target.collinear.chain \
      --noncollinear-chain results/single/02_chainpipeline_single/t2t.vs.single.T2T.vs.ASM.target.non-collinear.chain \
      --output-bed results/single/05_hapmers_single/single_switch_blocks_annotated.bed \
      --output-tsv results/single/05_hapmers_single/single_switch_blocks_chain_breakdown.tsv \
      --asm-name Single
"""

import os
import sys
import argparse
from collections import defaultdict


def parse_chain_blocks(chain_path):
    """
    Parses alignment blocks from a UCSC chain file indexed by query sequence (qName).
    Returns dict:
      {
        qName: list of (qs, qe, tName, ts, te, tStrand, qStrand)
      }
    where qs, qe are 0-based forward-strand query coordinates.
    """
    query_blocks = defaultdict(list)
    if not chain_path or not os.path.exists(chain_path):
        return query_blocks

    with open(chain_path) as fh:
        hdr = None
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            if line.startswith('chain'):
                p = line.split()
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                hdr = {
                    'tName': p[2],
                    'tSize': int(p[3]),
                    'tStrand': p[4],
                    'tStart': int(p[5]),
                    'tEnd': int(p[6]),
                    'qName': p[7],
                    'qSize': int(p[8]),
                    'qStrand': p[9],
                    'qStart': int(p[10]),
                    'qEnd': int(p[11])
                }
                t_pos = hdr['tStart']
                q_pos = hdr['qStart']
                q_strand = hdr['qStrand']
                q_size = hdr['qSize']
                t_strand = hdr['tStrand']
                continue

            if hdr is None:
                continue

            parts = line.split()
            if not parts:
                continue

            size = int(parts[0])
            dt = int(parts[1]) if len(parts) > 1 else 0
            dq = int(parts[2]) if len(parts) > 2 else 0

            if q_strand == '-':
                qs = q_size - (q_pos + size)
                qe = q_size - q_pos
            else:
                qs = q_pos
                qe = q_pos + size

            ts = t_pos
            te = t_pos + size

            query_blocks[hdr['qName']].append((min(qs, qe), max(qs, qe), hdr['tName'], ts, te, t_strand, q_strand))

            t_pos += size + dt
            q_pos += size + dq

    for qn in query_blocks:
        query_blocks[qn].sort(key=lambda x: x[0])

    return query_blocks


def get_overlap_and_liftover(b_start, b_end, blocks):
    """
    Computes overlap between [b_start, b_end) and query blocks.
    Projects intersecting sub-intervals to target (T2T) coordinates.
    Returns:
      total_cov_bp: total base pairs overlapping
      t2t_targets: dict {tName: list of (p_start, p_end)}
    """
    total_cov_bp = 0
    t2t_targets = defaultdict(list)

    for qs, qe, tn, ts, te, t_strand, q_strand in blocks:
        ov_s = max(b_start, qs)
        ov_e = min(b_end, qe)
        if ov_s < ov_e:
            cov = ov_e - ov_s
            total_cov_bp += cov

            offset_s = ov_s - qs
            offset_e = ov_e - qs

            if q_strand == '-':
                p_s = ts + (qe - ov_e)
                p_e = ts + (qe - ov_s)
            else:
                p_s = ts + offset_s
                p_e = ts + offset_e

            t2t_targets[tn].append((min(p_s, p_e), max(p_s, p_e)))

    return total_cov_bp, t2t_targets


def summarize_target_coords(t2t_targets):
    """
    Summarizes projected target coordinates.
    Returns (tName, min_start, max_end) or ('none', 0, 0).
    """
    if not t2t_targets:
        return 'none', 0, 0

    best_target = None
    best_bp = -1
    for tn, ivs in t2t_targets.items():
        bp = sum(e - s for s, e in ivs)
        if bp > best_bp:
            best_bp = bp
            best_target = tn

    all_s = [s for s, e in t2t_targets[best_target]]
    all_e = [e for s, e in t2t_targets[best_target]]
    return best_target, min(all_s), max(all_e)


def determine_label(coll_pct, nc_pct, unaligned_pct):
    """
    Assigns one of the 4 agreed labels:
      - collinear:     >= 90% collinear
      - non_collinear: >= 90% non-collinear
      - unaligned:     >= 90% unaligned
      - mixed:         any other combination
    """
    if coll_pct >= 90.0:
        return "collinear"
    elif nc_pct >= 90.0:
        return "non_collinear"
    elif unaligned_pct >= 90.0:
        return "unaligned"
    else:
        return "mixed"


def build_composition_detail(coll_pct, nc_pct, unaligned_pct):
    """Returns a readable composition string."""
    parts = []
    if coll_pct > 0.05:
        parts.append(f"collinear({coll_pct:.1f}%)")
    if nc_pct > 0.05:
        parts.append(f"non_collinear({nc_pct:.1f}%)")
    if unaligned_pct > 0.05:
        parts.append(f"unaligned({unaligned_pct:.1f}%)")
    return ' + '.join(parts) if parts else "unaligned(100.0%)"


def main():
    parser = argparse.ArgumentParser(description="Classify switch blocks into collinear, non-collinear, unaligned, or mixed.")
    parser.add_argument("--switch-bed", required=True, help="Input filtered switch blocks BED")
    parser.add_argument("--collinear-chain", required=True, help="Collinear chain file")
    parser.add_argument("--noncollinear-chain", required=True, help="Non-collinear chain file")
    parser.add_argument("--output-bed", required=True, help="Output annotated BED (Option A: original coords + label)")
    parser.add_argument("--output-tsv", required=True, help="Output detailed breakdown TSV")
    parser.add_argument("--asm-name", default="Assembly", help="Assembly label (Single, Dual, ONT)")

    args = parser.parse_args()

    print(f"[{args.asm_name}] Loading collinear chain: {args.collinear_chain}...", file=sys.stderr)
    coll_blocks = parse_chain_blocks(args.collinear_chain)

    print(f"[{args.asm_name}] Loading non-collinear chain: {args.noncollinear_chain}...", file=sys.stderr)
    nc_blocks = parse_chain_blocks(args.noncollinear_chain)

    # Read switch blocks
    records = []
    if os.path.exists(args.switch_bed) and os.path.getsize(args.switch_bed) > 0:
        with open(args.switch_bed) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                scaf = parts[0]
                bs = int(parts[1])
                be = int(parts[2])
                extra = parts[3:] if len(parts) > 3 else []
                records.append((scaf, bs, be, extra, line))

    rows = []
    bed_entries = []

    for scaf, bs, be, extra, orig_line in records:
        blen = be - bs
        coll_bp, coll_targets = get_overlap_and_liftover(bs, be, coll_blocks.get(scaf, []))
        nc_bp, nc_targets = get_overlap_and_liftover(bs, be, nc_blocks.get(scaf, []))

        # Unaligned bp is whatever bases are not covered by either chain
        unaligned_bp = max(0, blen - (coll_bp + nc_bp))

        coll_pct = (coll_bp / blen) * 100.0 if blen > 0 else 0.0
        nc_pct = (nc_bp / blen) * 100.0 if blen > 0 else 0.0
        unaligned_pct = (unaligned_bp / blen) * 100.0 if blen > 0 else 0.0

        label = determine_label(coll_pct, nc_pct, unaligned_pct)
        detail = build_composition_detail(coll_pct, nc_pct, unaligned_pct)

        # Liftover T2T coordinates: prioritize chain with highest coverage
        combined_targets = defaultdict(list)
        if coll_bp >= nc_bp and coll_targets:
            tn, ts, te = summarize_target_coords(coll_targets)
        elif nc_targets:
            tn, ts, te = summarize_target_coords(nc_targets)
        elif coll_targets:
            tn, ts, te = summarize_target_coords(coll_targets)
        else:
            tn, ts, te = "unaligned", 0, 0

        hapmer_id = extra[0] if extra else "switch_block"

        rows.append({
            'scaffold': scaf,
            'start': bs,
            'end': be,
            'length': blen,
            'label': label,
            't2t_chrom': tn,
            't2t_start': ts,
            't2t_end': te,
            'collinear_bp': coll_bp,
            'collinear_pct': f"{coll_pct:.1f}%",
            'noncollinear_bp': nc_bp,
            'noncollinear_pct': f"{nc_pct:.1f}%",
            'unaligned_bp': unaligned_bp,
            'unaligned_pct': f"{unaligned_pct:.1f}%",
            'composition_detail': detail,
            'hapmer_info': '\t'.join(extra)
        })

        # BED format: scaffold, start, end, label, hapmer_id, length, extra...
        bed_entries.append(f"{scaf}\t{bs}\t{be}\t{label}\t{hapmer_id}\t{blen}")

    # Write output BED (Option A: original scaffold coordinates with label)
    os.makedirs(os.path.dirname(os.path.abspath(args.output_bed)), exist_ok=True)
    with open(args.output_bed, 'w') as fh:
        for entry in bed_entries:
            fh.write(entry + '\n')
    print(f"[{args.asm_name}] Wrote annotated BED ({len(bed_entries)} blocks): {args.output_bed}", file=sys.stderr)

    # Write output detailed TSV
    os.makedirs(os.path.dirname(os.path.abspath(args.output_tsv)), exist_ok=True)
    tsv_headers = [
        "scaffold", "start", "end", "length", "label",
        "t2t_chrom", "t2t_start", "t2t_end",
        "collinear_bp", "collinear_pct",
        "noncollinear_bp", "noncollinear_pct",
        "unaligned_bp", "unaligned_pct",
        "composition_detail", "hapmer_info"
    ]
    with open(args.output_tsv, 'w') as fh:
        fh.write('\t'.join(tsv_headers) + '\n')
        for r in rows:
            fh.write(f"{r['scaffold']}\t{r['start']}\t{r['end']}\t{r['length']}\t{r['label']}\t"
                     f"{r['t2t_chrom']}\t{r['t2t_start']}\t{r['t2t_end']}\t"
                     f"{r['collinear_bp']}\t{r['collinear_pct']}\t"
                     f"{r['noncollinear_bp']}\t{r['noncollinear_pct']}\t"
                     f"{r['unaligned_bp']}\t{r['unaligned_pct']}\t"
                     f"{r['composition_detail']}\t{r['hapmer_info']}\n")
    print(f"[{args.asm_name}] Wrote detailed TSV: {args.output_tsv}", file=sys.stderr)

    # Print summary breakdown to terminal
    n_total = len(rows)
    n_coll = sum(1 for r in rows if r['label'] == 'collinear')
    n_nc = sum(1 for r in rows if r['label'] == 'non_collinear')
    n_un = sum(1 for r in rows if r['label'] == 'unaligned')
    n_mix = sum(1 for r in rows if r['label'] == 'mixed')

    print(f"\n===========================================================", file=sys.stderr)
    print(f" Switch Blocks Chain Classification Summary: {args.asm_name}", file=sys.stderr)
    print(f"===========================================================", file=sys.stderr)
    print(f"Total Filtered Switch Blocks : {n_total}", file=sys.stderr)
    print(f"  - Collinear (>= 90%)       : {n_coll}" + (f" ({n_coll/n_total*100:.1f}%)" if n_total else ""), file=sys.stderr)
    print(f"  - Non-collinear (>= 90%)   : {n_nc}" + (f" ({n_nc/n_total*100:.1f}%)" if n_total else ""), file=sys.stderr)
    print(f"  - Mixed (spans multiple)   : {n_mix}" + (f" ({n_mix/n_total*100:.1f}%)" if n_total else ""), file=sys.stderr)
    print(f"  - Unaligned (>= 90%)       : {n_un}" + (f" ({n_un/n_total*100:.1f}%)" if n_total else ""), file=sys.stderr)
    print(f"===========================================================\n", file=sys.stderr)


if __name__ == "__main__":
    main()
