#!/usr/bin/env python3
"""
recover_query_uncovered_chains.py

Identifies regions of the assembly not covered by any chain in the main pipeline
(collinear or non-collinear), bridges gaps <= max_skip_gap (default: 20 kb),
computes query-centric 1x netting using UCSC tools, and extracts the chains
aligning those previously uncovered assembly regions.

Author: Antigravity
Date: 2026-09-24
"""

import os
import sys
import argparse
import subprocess
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract uncovered assembly regions from main chain pipeline and recover query-centric 1x chains."
    )
    parser.add_argument(
        "--sorted-1to1-chain", required=True,
        help="Path to t2t.vs.asm.sorted.1to1.chain (from Stage 8 before target chainNet)"
    )
    parser.add_argument(
        "--collinear-chain", required=True,
        help="Path to t2t.vs.asm.target.collinear.chain"
    )
    parser.add_argument(
        "--noncollinear-chain", required=True,
        help="Path to t2t.vs.asm.target.non-collinear.chain"
    )
    parser.add_argument(
        "--asm-sizes", required=True,
        help="Assembly chromosome sizes file (*.sizes)"
    )
    parser.add_argument(
        "--t2t-sizes", required=True,
        help="T2T reference chromosome sizes file (*.sizes)"
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Output directory for generated BED, chain, and summary files"
    )
    parser.add_argument(
        "--max-skip-gap", type=int, default=20000,
        help="Maximum query gap size (bp) inside/between alignments to skip/bridge (default: 20000 bp)"
    )
    parser.add_argument(
        "--min-uncovered-size", type=int, default=20000,
        help="Minimum size (bp) of an uncovered region to report and recover (default: 20000 bp)"
    )
    parser.add_argument(
        "--ucsc-bin", default="/lustre/fs5/vgl/scratch/eduarte/miniconda3/envs/ucsc-tools/bin",
        help="Directory containing UCSC tool binaries (chainSwap, chainSort, chainPreNet, chainNet, netChainSubset)"
    )
    return parser.parse_args()


def load_sizes(sizes_path):
    """Load sequence names and lengths from a UCSC .sizes file."""
    sizes = {}
    with open(sizes_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def extract_covered_query_blocks(chain_path):
    """
    Parses a chain file and returns all aligned query intervals:
    dict: q_name -> list of (q_start, q_end)
    """
    blocks = defaultdict(list)
    if not chain_path or not os.path.exists(chain_path):
        return blocks

    with open(chain_path) as f:
        in_chain = False
        q_name = None
        q_pos = 0

        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('chain'):
                parts = line.split()
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                q_name = parts[7]
                q_strand = parts[9]
                q_size = int(parts[8])
                q_raw_start = int(parts[10])
                q_raw_end = int(parts[11])
                if q_strand == '-':
                    q_pos = q_size - q_raw_end
                else:
                    q_pos = q_raw_start
                in_chain = True
            else:
                parts = line.split()
                size = int(parts[0])
                dq = int(parts[2]) if len(parts) > 2 else 0

                b_start = q_pos
                b_end = q_pos + size
                if b_end > b_start:
                    blocks[q_name].append((b_start, b_end))

                q_pos += size + dq
                if len(parts) == 1:
                    in_chain = False

    return blocks


def merge_and_bridge_intervals(intervals, max_skip_gap=20000):
    """
    Merges overlapping intervals and bridges gaps <= max_skip_gap.
    Returns sorted, disjoint list of [start, end].
    """
    if not intervals:
        return []

    sorted_ivs = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = []
    cur_s, cur_e = sorted_ivs[0]

    for s, e in sorted_ivs[1:]:
        if s <= cur_e + max_skip_gap:
            # Overlaps or gap is <= max_skip_gap: bridge across
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged


def compute_uncovered_intervals(asm_sizes, covered_blocks, max_skip_gap=20000, min_uncovered_size=20000):
    """
    Subtracts merged covered intervals from chromosome total lengths.
    Returns list of dicts with gap info.
    """
    uncovered = []

    for q_name, total_len in asm_sizes.items():
        q_ivs = covered_blocks.get(q_name, [])
        merged = merge_and_bridge_intervals(q_ivs, max_skip_gap=max_skip_gap)

        if not merged:
            # Whole sequence uncovered
            if total_len >= min_uncovered_size:
                uncovered.append({
                    'chrom': q_name,
                    'start': 0,
                    'end': total_len,
                    'size': total_len,
                    'type': 'whole_chromosome_uncovered'
                })
            continue

        # Check 5' end
        if merged[0][0] >= min_uncovered_size:
            uncovered.append({
                'chrom': q_name,
                'start': 0,
                'end': merged[0][0],
                'size': merged[0][0],
                'type': '5prime_subtelomeric_gap'
            })

        # Check internal gaps
        for i in range(len(merged) - 1):
            gap_s = merged[i][1]
            gap_e = merged[i + 1][0]
            gap_len = gap_e - gap_s
            if gap_len >= min_uncovered_size:
                uncovered.append({
                    'chrom': q_name,
                    'start': gap_s,
                    'end': gap_e,
                    'size': gap_len,
                    'type': 'internal_structural_gap'
                })

        # Check 3' end
        if (total_len - merged[-1][1]) >= min_uncovered_size:
            uncovered.append({
                'chrom': q_name,
                'start': merged[-1][1],
                'end': total_len,
                'size': total_len - merged[-1][1],
                'type': '3prime_subtelomeric_gap'
            })

    return uncovered


def run_cmd(cmd, env=None):
    """Run a shell command, printing and verifying exit code."""
    print(f"[CMD] {cmd}")
    res = subprocess.run(cmd, shell=True, env=env, check=True)
    return res


def run_query_1x_netting(sorted_1to1_chain, asm_sizes_path, t2t_sizes_path, outdir, ucsc_bin):
    """
    Uses UCSC tools to build a query-centric 1x net from sorted.1to1.chain:
    1. chainSwap (ASM becomes Target, T2T becomes Query)
    2. chainSort by target
    3. chainPreNet & chainNet (with ASM.sizes as target, T2T.sizes as query)
    4. netChainSubset -> asm_1x.chain (1x coverage from assembly side)
    5. chainSwap back -> asm_1x_t2t_target.chain (T2T target orientation)
    """
    env = os.environ.copy()
    env["PATH"] = f"{ucsc_bin}:{env.get('PATH', '')}"

    asm_as_target_chain = os.path.join(outdir, "asm_as_target.chain")
    asm_net = os.path.join(outdir, "asm_as_target.net")
    asm_1x_chain = os.path.join(outdir, "asm_query_1x_asm_target.chain")
    asm_1x_t2t_target_chain = os.path.join(outdir, "asm_query_1x_t2t_target.chain")

    print("\n=== Stage B: Query-side 1x Netting via UCSC Tools ===")

    # 1. Swap chain so ASM is target, and sort
    cmd1 = f"chainSwap '{sorted_1to1_chain}' stdout | chainSort stdin '{asm_as_target_chain}'"
    run_cmd(cmd1, env=env)

    # 2. chainPreNet & chainNet
    cmd2 = (
        f"chainPreNet '{asm_as_target_chain}' '{asm_sizes_path}' '{t2t_sizes_path}' stdout | "
        f"chainNet -minSpace=1 -minScore=0 stdin '{asm_sizes_path}' '{t2t_sizes_path}' '{asm_net}' /dev/null"
    )
    run_cmd(cmd2, env=env)

    # 3. netChainSubset to extract 1x chains
    cmd3 = f"netChainSubset '{asm_net}' '{asm_as_target_chain}' '{asm_1x_chain}'"
    run_cmd(cmd3, env=env)

    # 4. Swap back to standard T2T-as-target orientation
    cmd4 = f"chainSwap '{asm_1x_chain}' stdout | chainSort stdin '{asm_1x_t2t_target_chain}'"
    run_cmd(cmd4, env=env)

    print(f"  Query-centric 1x chains (ASM target): {asm_1x_chain}")
    print(f"  Query-centric 1x chains (T2T target): {asm_1x_t2t_target_chain}")
    return asm_1x_chain, asm_1x_t2t_target_chain


def parse_chains_with_blocks(chain_path, t2t_is_target=True):
    """
    Parses a chain file and stores each chain record with its alignment blocks.
    Returns: list of dicts:
      {
        'header': line,
        'tName': ..., 'tStart': ..., 'tEnd': ..., 'tStrand': ...,
        'qName': ..., 'qStart': ..., 'qEnd': ..., 'qStrand': ...,
        'score': ..., 'id': ...,
        'blocks': [(q_s, q_e, t_s, t_e), ...],
        'lines': [raw_block_lines]
      }
    """
    chains = []
    if not os.path.exists(chain_path):
        return chains

    with open(chain_path) as f:
        cur = None
        for line in f:
            raw = line
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('chain'):
                if cur:
                    chains.append(cur)
                p = line.split()
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                q_size = int(p[8])
                q_strand = p[9]
                q_start = int(p[10])
                q_end = int(p[11])
                q_pos = (q_size - q_end) if q_strand == '-' else q_start
                t_pos = int(p[5])

                cur = {
                    'header': line,
                    'score': int(p[1]),
                    'tName': p[2], 'tSize': int(p[3]), 'tStrand': p[4],
                    'tStart': int(p[5]), 'tEnd': int(p[6]),
                    'qName': p[7], 'qSize': q_size, 'qStrand': q_strand,
                    'qStart': q_start, 'qEnd': q_end,
                    'id': p[12] if len(p) > 12 else '0',
                    'q_pos': q_pos,
                    't_pos': t_pos,
                    'blocks': [],
                    'lines': []
                }
            else:
                if cur:
                    p = line.split()
                    cur['lines'].append(raw)
                    size = int(p[0])
                    dt = int(p[1]) if len(p) > 1 else 0
                    dq = int(p[2]) if len(p) > 2 else 0

                    bs_q = cur['q_pos']
                    be_q = cur['q_pos'] + size
                    bs_t = cur['t_pos']
                    be_t = cur['t_pos'] + size

                    cur['blocks'].append((bs_q, be_q, bs_t, be_t))
                    cur['q_pos'] += size + dq
                    cur['t_pos'] += size + dt

        if cur:
            chains.append(cur)

    return chains


def recover_uncovered_chains(uncovered_regions, asm_1x_t2t_target_chain, outdir):
    """
    Matches the uncovered regions against the query-side 1x net chains (in T2T-as-target orientation).
    Writes:
      1. asm_recovered_uncovered_chains.bed
      2. asm_recovered_uncovered.chain
      3. uncovered_and_recovery_summary.tsv
    """
    print("\n=== Stage C: Intersecting and Recovering Chains for Uncovered Regions ===")

    chains = parse_chains_with_blocks(asm_1x_t2t_target_chain, t2t_is_target=True)

    # Index chains by qName
    chains_by_q = defaultdict(list)
    for c in chains:
        chains_by_q[c['qName']].append(c)

    recovered_bed_records = []
    recovered_chain_ids = set()
    summary_records = []

    for uncov in uncovered_regions:
        q_chrom = uncov['chrom']
        u_s = uncov['start']
        u_e = uncov['end']
        u_size = uncov['size']
        u_type = uncov['type']

        overlapping_chains = []
        covered_bases_in_gap = 0

        for c in chains_by_q.get(q_chrom, []):
            # Check overlap between [u_s, u_e] and chain query bounds [qStart, qEnd]
            if not (c['qEnd'] <= u_s or c['qStart'] >= u_e):
                # Calculate aligned block overlap
                ov_bases = 0
                for bq_s, bq_e, bt_s, bt_e in c['blocks']:
                    os_ = max(bq_s, u_s)
                    oe_ = min(bq_e, u_e)
                    if oe_ > os_:
                        ov_bases += (oe_ - os_)
                        # Convert query overlap coordinates to target span
                        offset_s = os_ - bq_s
                        offset_e = oe_ - bq_s
                        rec_t_s = bt_s + offset_s
                        rec_t_e = bt_s + offset_e

                        recovered_bed_records.append({
                            'q_chrom': q_chrom,
                            'q_start': os_,
                            'q_end': oe_,
                            'gap_id': f"{q_chrom}:{u_s}-{u_e}",
                            't_chrom': c['tName'],
                            't_start': rec_t_s,
                            't_end': rec_t_e,
                            'strand': c['qStrand'],
                            'chain_id': c['id'],
                            'score': c['score']
                        })

                if ov_bases > 0:
                    covered_bases_in_gap += ov_bases
                    overlapping_chains.append((c, ov_bases))
                    recovered_chain_ids.add(c['id'])

        cov_pct = (covered_bases_in_gap / u_size * 100.0) if u_size > 0 else 0.0

        target_targets = sorted(list(set(
            f"{c[0]['tName']}:{c[0]['tStart']}-{c[0]['tEnd']}({c[0]['qStrand']})"
            for c in overlapping_chains
        )))
        target_summary_str = ";".join(target_targets) if target_targets else "None"

        summary_records.append({
            'chrom': q_chrom,
            'gap_start': u_s,
            'gap_end': u_e,
            'gap_size_bp': u_size,
            'gap_type': u_type,
            'recovered_aligned_bases': covered_bases_in_gap,
            'recovered_pct': round(cov_pct, 2),
            'num_recovering_chains': len(overlapping_chains),
            'target_mapping_regions': target_summary_str
        })

    # Write BED of uncovered regions
    uncov_bed_path = os.path.join(outdir, "asm_uncovered_main_pipeline.bed")
    with open(uncov_bed_path, 'w') as f:
        f.write("#chrom\tstart\tend\tgap_size_bp\tgap_type\n")
        for u in uncovered_regions:
            f.write(f"{u['chrom']}\t{u['start']}\t{u['end']}\t{u['size']}\t{u['type']}\n")
    print(f"  Wrote: {uncov_bed_path} ({len(uncovered_regions)} uncovered intervals)")

    # Write BED of recovered segments
    rec_bed_path = os.path.join(outdir, "asm_recovered_uncovered_chains.bed")
    with open(rec_bed_path, 'w') as f:
        f.write("#q_chrom\tq_start\tq_end\tgap_id\tt_chrom\tt_start\tt_end\tstrand\tchain_id\tscore\n")
        for r in recovered_bed_records:
            f.write(
                f"{r['q_chrom']}\t{r['q_start']}\t{r['q_end']}\t{r['gap_id']}\t"
                f"{r['t_chrom']}\t{r['t_start']}\t{r['t_end']}\t{r['strand']}\t"
                f"{r['chain_id']}\t{r['score']}\n"
            )
    print(f"  Wrote: {rec_bed_path} ({len(recovered_bed_records)} aligned segments recovered)")

    # Write filtered chain file containing only the chains that recover uncovered regions
    rec_chain_path = os.path.join(outdir, "asm_recovered_uncovered.chain")
    with open(rec_chain_path, 'w') as f:
        for c in chains:
            if c['id'] in recovered_chain_ids:
                f.write(c['header'] + '\n')
                for line in c['lines']:
                    f.write(line)
                f.write('\n')
    print(f"  Wrote: {rec_chain_path} ({len(recovered_chain_ids)} recovering chains)")

    # Write summary TSV
    summary_tsv_path = os.path.join(outdir, "uncovered_and_recovery_summary.tsv")
    with open(summary_tsv_path, 'w') as f:
        cols = [
            'chrom', 'gap_start', 'gap_end', 'gap_size_bp', 'gap_type',
            'recovered_aligned_bases', 'recovered_pct', 'num_recovering_chains',
            'target_mapping_regions'
        ]
        f.write("\t".join(cols) + "\n")
        for s in summary_records:
            f.write("\t".join(str(s[c]) for c in cols) + "\n")
    print(f"  Wrote: {summary_tsv_path}")

    # Print high-level console summary of notable structural gaps
    print("\n--- Key Recovered Structural Gaps (>= 100 kb) ---")
    header_fmt = "{:<20} {:>12} {:>12} {:>12} {:>12} {:>10} {:<35}"
    print(header_fmt.format("Chromosome", "Start", "End", "Size(bp)", "Recovered(bp)", "%Recov", "Target Region"))
    print("-" * 115)
    for s in summary_records:
        if s['gap_size_bp'] >= 100000:
            target_disp = s['target_mapping_regions'][:33] + ".." if len(s['target_mapping_regions']) > 35 else s['target_mapping_regions']
            print(header_fmt.format(
                s['chrom'],
                f"{s['gap_start']:,}",
                f"{s['gap_end']:,}",
                f"{s['gap_size_bp']:,}",
                f"{s['recovered_aligned_bases']:,}",
                f"{s['recovered_pct']}%",
                target_disp
            ))


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("=================================================================")
    print(" Recovering Query-Centric 1x Chains for Uncovered Assembly Gaps")
    print("=================================================================")
    print(f"Sorted 1-to-1 Chain : {args.sorted_1to1_chain}")
    print(f"Collinear Chain     : {args.collinear_chain}")
    print(f"Non-collinear Chain : {args.noncollinear_chain}")
    print(f"Assembly Sizes      : {args.asm_sizes}")
    print(f"T2T Sizes           : {args.t2t_sizes}")
    print(f"Max Skip Gap        : {args.max_skip_gap:,} bp")
    print(f"Min Uncovered Size  : {args.min_uncovered_size:,} bp")
    print(f"Output Directory    : {args.outdir}")
    print("=================================================================\n")

    # 1. Load sequence lengths
    asm_sizes = load_sizes(args.asm_sizes)

    # 2. Extract covered blocks from main pipeline chains
    print("=== Stage A: Extracting Main Pipeline Covered Spans ===")
    col_blocks = extract_covered_query_blocks(args.collinear_chain)
    nc_blocks  = extract_covered_query_blocks(args.noncollinear_chain)

    all_covered = defaultdict(list)
    for qn, blks in col_blocks.items():
        all_covered[qn].extend(blks)
    for qn, blks in nc_blocks.items():
        all_covered[qn].extend(blks)

    print(f"  Collinear query chromosomes    : {len(col_blocks)}")
    print(f"  Non-collinear query chromosomes: {len(nc_blocks)}")

    # 3. Compute uncovered intervals (bridging gaps <= max_skip_gap)
    uncovered = compute_uncovered_intervals(
        asm_sizes,
        all_covered,
        max_skip_gap=args.max_skip_gap,
        min_uncovered_size=args.min_uncovered_size
    )
    print(f"  Total uncovered regions >= {args.min_uncovered_size:,} bp: {len(uncovered)}")

    # 4. Run UCSC query-side 1x netting
    asm_1x_asm_target, asm_1x_t2t_target = run_query_1x_netting(
        args.sorted_1to1_chain,
        args.asm_sizes,
        args.t2t_sizes,
        args.outdir,
        args.ucsc_bin
    )

    # 5. Recover chains for uncovered regions
    recover_uncovered_chains(uncovered, asm_1x_t2t_target, args.outdir)

    print("\n=== Process Completed Successfully! ===")


if __name__ == "__main__":
    main()
