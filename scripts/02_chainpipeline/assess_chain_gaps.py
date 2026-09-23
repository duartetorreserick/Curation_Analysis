#!/usr/bin/env python3
"""
assess_chain_gaps.py — Analyze gap length distributions within UCSC chain files.

For each chain file (collinear, non-collinear), this script analyzes the internal
gaps defined on block lines (size  dt  dq):
  - Single-sided target gaps (dt > 0, dq == 0)
  - Single-sided query gaps (dt == 0, dq > 0)
  - Double-sided gaps (dt > 0, dq > 0)

Usage:
  python3 assess_chain_gaps.py [CHAIN_FILES ...]
  or
  python3 assess_chain_gaps.py --results-dir /path/to/results/
"""

import sys
import os
import glob
import argparse
import numpy as np


def parse_chain_gaps(chain_file):
    """Parse gap lengths from a chain file."""
    t_only_gaps = []
    q_only_gaps = []
    both_dt = []
    both_dq = []
    both_max = []
    
    chain_count = 0
    total_aligned_bases = 0
    chain_aligned_lengths = []
    
    with open(chain_file, 'r') as fh:
        cur_aligned = 0
        in_chain = False
        
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('chain'):
                if in_chain:
                    chain_aligned_lengths.append(cur_aligned)
                    total_aligned_bases += cur_aligned
                chain_count += 1
                cur_aligned = 0
                in_chain = True
                continue
            
            parts = line.split()
            if len(parts) == 3:
                size, dt, dq = int(parts[0]), int(parts[1]), int(parts[2])
                cur_aligned += size
                
                if dt > 0 and dq == 0:
                    t_only_gaps.append(dt)
                elif dt == 0 and dq > 0:
                    q_only_gaps.append(dq)
                elif dt > 0 and dq > 0:
                    both_dt.append(dt)
                    both_dq.append(dq)
                    both_max.append(max(dt, dq))
            elif len(parts) == 1:
                size = int(parts[0])
                cur_aligned += size
        
        if in_chain:
            chain_aligned_lengths.append(cur_aligned)
            total_aligned_bases += cur_aligned

    return {
        'chain_count': chain_count,
        'total_aligned_bases': total_aligned_bases,
        'chain_aligned_lengths': chain_aligned_lengths,
        't_only_gaps': np.array(t_only_gaps, dtype=np.int64),
        'q_only_gaps': np.array(q_only_gaps, dtype=np.int64),
        'single_gaps': np.array(t_only_gaps + q_only_gaps, dtype=np.int64),
        'both_dt': np.array(both_dt, dtype=np.int64),
        'both_dq': np.array(both_dq, dtype=np.int64),
        'both_max': np.array(both_max, dtype=np.int64),
    }


def compute_stats(arr):
    if len(arr) == 0:
        return None
    return {
        'count': len(arr),
        'min': int(np.min(arr)),
        'p25': float(np.percentile(arr, 25)),
        'median': float(np.median(arr)),
        'p75': float(np.percentile(arr, 75)),
        'p90': float(np.percentile(arr, 90)),
        'p95': float(np.percentile(arr, 95)),
        'p99': float(np.percentile(arr, 99)),
        'max': int(np.max(arr)),
    }


def get_bins(arr):
    if len(arr) == 0:
        return {}
    bins = [
        ("1-10 bp", np.sum((arr >= 1) & (arr <= 10))),
        ("11-50 bp", np.sum((arr >= 11) & (arr <= 50))),
        ("51-100 bp", np.sum((arr >= 51) & (arr <= 100))),
        ("101-500 bp", np.sum((arr >= 101) & (arr <= 500))),
        ("501-2000 bp", np.sum((arr >= 501) & (arr <= 2000))),
        ("2001-10000 bp", np.sum((arr >= 2001) & (arr <= 10000))),
        (">10000 bp", np.sum(arr > 10000)),
    ]
    return bins


def print_report(label, data):
    print("=" * 80)
    print(f" FILE: {label}")
    print(f" Total Chains: {data['chain_count']:,} | Total Aligned: {data['total_aligned_bases']:,} bp")
    print("=" * 80)
    
    categories = [
        ("Single-sided Gaps (All: Target + Query)", data['single_gaps']),
        ("  - Target-only gaps (dt > 0, dq = 0)", data['t_only_gaps']),
        ("  - Query-only gaps (dt = 0, dq > 0)", data['q_only_gaps']),
        ("Double-sided Gaps (Both dt > 0 and dq > 0) [max(dt, dq)]", data['both_max']),
    ]
    
    for title, arr in categories:
        stats = compute_stats(arr)
        print(f"\n--- {title} ---")
        if stats is None:
            print("  No gaps found in this category.")
            continue
        print(f"  Count:  {stats['count']:,}")
        print(f"  Min:    {stats['min']:,} bp")
        print(f"  25%:    {stats['p25']:.1f} bp")
        print(f"  Median: {stats['median']:.1f} bp")
        print(f"  75%:    {stats['p75']:.1f} bp")
        print(f"  90%:    {stats['p90']:.1f} bp")
        print(f"  95%:    {stats['p95']:.1f} bp")
        print(f"  99%:    {stats['p99']:.1f} bp")
        print(f"  Max:    {stats['max']:,} bp")
        
        bins = get_bins(arr)
        total = stats['count']
        print("\n  Size distribution breakdown:")
        for b_name, b_count in bins:
            pct = (b_count / total * 100) if total > 0 else 0
            bar = "#" * int(pct / 2)
            print(f"    {b_name:<15} : {b_count:>8,} ({pct:>5.1f}%)  {bar}")


def main():
    parser = argparse.ArgumentParser(description="Analyze gap lengths in chain files.")
    parser.add_argument("files", nargs="*", help="Chain files to analyze")
    parser.add_argument("--results-dir", help="Base results dir containing assembly subdirs")
    args = parser.parse_args()
    
    chain_files = list(args.files)
    if args.results_dir:
        pattern = os.path.join(args.results_dir, "**", "*collinear*.chain")
        found = glob.glob(pattern, recursive=True)
        # Sort to ensure consistent order
        chain_files.extend(sorted(found))
    
    if not chain_files:
        print("ERROR: No chain files specified or found.", file=sys.stderr)
        sys.exit(1)
        
    for cf in chain_files:
        if not os.path.isfile(cf):
            continue
        rel_name = os.path.basename(cf)
        parent = os.path.basename(os.path.dirname(os.path.dirname(cf)))
        label = f"{parent} / {os.path.basename(os.path.dirname(cf))} / {rel_name}"
        data = parse_chain_gaps(cf)
        print_report(label, data)
        print("\n")


if __name__ == "__main__":
    main()
