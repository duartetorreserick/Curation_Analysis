#!/usr/bin/env python3
"""
chain_coverage.py — Compute per-sequence and genome-wide coverage statistics
from one or more UCSC chain files.

For each chain file the script collects alignment blocks for target (T2T) and
query (curated assembly), merges overlapping intervals per sequence, and
reports:
  - per-sequence coverage tables  (<prefix>.target_cov.tsv, <prefix>.query_cov.tsv)
  - a human-readable summary       (<prefix>.summary.txt)

Usage
-----
  python3 chain_coverage.py \
      --chain  FILE [FILE ...] \
      [--target-sizes FILE] \
      [--query-sizes  FILE] \
      [--output-prefix STR]          (default: chain_cov)
      [--min-score INT]              (skip chains with score < INT; default: 0)

Chain file format (UCSC)
------------------------
  chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
  size  [dt  dq]
  ...
  size
  <blank line>

Both target and query coordinates are always in positive-strand space.
"""

import argparse
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# Interval merging
# ---------------------------------------------------------------------------

def merge_intervals(intervals):
    """Return merged, non-overlapping intervals from a list of (start, end) tuples."""
    if not intervals:
        return []
    intervals.sort()
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def covered_bases(intervals):
    """Total bases covered by a list of (start, end) intervals (half-open)."""
    return sum(e - s for s, e in merge_intervals(intervals))


# ---------------------------------------------------------------------------
# Chain parser
# ---------------------------------------------------------------------------

def parse_chain_file(path, min_score=0):
    """
    Yield alignment-block intervals for each chain.

    Yields dicts with keys:
      score, tName, tSize, tStrand, tStart, tEnd,
      qName, qSize, qStrand, qStart, qEnd, id,
      target_intervals  [(tStart, tEnd), ...],
      query_intervals   [(qStart, qEnd), ...]
    """
    with open(path) as fh:
        header = None
        blocks = []

        for raw in fh:
            line = raw.rstrip('\n')

            # Skip comment / empty lines that are not block separators
            if line.startswith('#'):
                continue

            if line.startswith('chain'):
                # Save previous chain if present
                if header is not None:
                    yield _build_chain(header, blocks)
                parts = line.split()
                header = {
                    'score':   int(parts[1]),
                    'tName':   parts[2],
                    'tSize':   int(parts[3]),
                    'tStrand': parts[4],
                    'tStart':  int(parts[5]),
                    'tEnd':    int(parts[6]),
                    'qName':   parts[7],
                    'qSize':   int(parts[8]),
                    'qStrand': parts[9],
                    'qStart':  int(parts[10]),
                    'qEnd':    int(parts[11]),
                    'id':      parts[12] if len(parts) > 12 else '.',
                    'score_val': int(parts[1]),
                }
                blocks = []

            elif line == '':
                # Block separator — save chain
                if header is not None:
                    yield _build_chain(header, blocks)
                    header = None
                    blocks = []

            else:
                # Alignment block data line: "size [dt dq]"
                parts = line.split()
                if parts:
                    blocks.append(tuple(int(x) for x in parts))

        # End of file — save last chain
        if header is not None:
            yield _build_chain(header, blocks)


def _build_chain(header, blocks):
    """Reconstruct per-position intervals for target and query from block data."""
    t_pos = header['tStart']
    q_pos = header['qStart']
    t_intervals = []
    q_intervals = []

    for i, block in enumerate(blocks):
        size = block[0]
        t_intervals.append((t_pos, t_pos + size))
        q_intervals.append((q_pos, q_pos + size))
        # Last block has no gap fields
        if len(block) == 3:
            dt, dq = block[1], block[2]
            t_pos += size + dt
            q_pos += size + dq
        else:
            t_pos += size
            q_pos += size

    chain = dict(header)
    chain['target_intervals'] = t_intervals
    chain['query_intervals']  = q_intervals
    return chain


# ---------------------------------------------------------------------------
# Coverage accumulation
# ---------------------------------------------------------------------------

def accumulate_coverage(chain_files, min_score=0):
    """
    Parse all chain files and accumulate raw intervals per sequence.

    Returns:
      target_ivs  {tName: [(start, end), ...]}
      target_sizes {tName: int}    (from chain headers — max seen size)
      query_ivs   {qName: [(start, end), ...]}
      query_sizes {qName: int}
    """
    target_ivs   = defaultdict(list)
    target_sizes = {}
    query_ivs    = defaultdict(list)
    query_sizes  = {}

    for path in chain_files:
        print(f"  Reading: {path}", file=sys.stderr)
        n = 0
        for chain in parse_chain_file(path, min_score):
            if chain['score_val'] < min_score:
                continue
            n += 1
            tn = chain['tName']
            qn = chain['qName']
            # Track largest size seen per sequence (should be constant, but safe)
            target_sizes[tn] = max(target_sizes.get(tn, 0), chain['tSize'])
            query_sizes[qn]  = max(query_sizes.get(qn,  0), chain['qSize'])
            target_ivs[tn].extend(chain['target_intervals'])
            query_ivs[qn].extend(chain['query_intervals'])
        print(f"    → {n} chains loaded", file=sys.stderr)

    return target_ivs, target_sizes, query_ivs, query_sizes


# ---------------------------------------------------------------------------
# Size file loader
# ---------------------------------------------------------------------------

def load_sizes(path):
    """Load a two-column TSV (name, size) into a dict."""
    sizes = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            sizes[parts[0]] = int(parts[1])
    return sizes


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def write_coverage_table(ivs, sizes, out_path, label):
    """
    Write per-sequence coverage table.

    Columns: name, size, covered_bp, uncovered_bp, pct_covered
    """
    rows = []
    for name, intervals in ivs.items():
        seq_size = sizes.get(name, 0)
        cov = covered_bases(intervals)
        uncov = seq_size - cov if seq_size else 0
        pct = cov / seq_size * 100 if seq_size else 0.0
        rows.append((name, seq_size, cov, uncov, pct))

    # Sequences present in sizes but absent from any chain
    for name, seq_size in sizes.items():
        if name not in ivs:
            rows.append((name, seq_size, 0, seq_size, 0.0))

    rows.sort(key=lambda r: (-r[1], r[0]))

    with open(out_path, 'w') as fh:
        fh.write(f"# {label} per-sequence coverage\n")
        fh.write("name\tsize\tcovered_bp\tuncovered_bp\tpct_covered\n")
        for name, size, cov, uncov, pct in rows:
            fh.write(f"{name}\t{size}\t{cov}\t{uncov}\t{pct:.4f}\n")

    return rows


def write_summary(target_rows, query_rows, target_total, query_total, out_path, chain_files):
    """Write a human-readable summary."""
    t_cov_total = sum(r[2] for r in target_rows)
    q_cov_total = sum(r[2] for r in query_rows)

    t_size_total = target_total if target_total else sum(r[1] for r in target_rows)
    q_size_total = query_total  if query_total  else sum(r[1] for r in query_rows)

    t_pct = t_cov_total / t_size_total * 100 if t_size_total else 0.0
    q_pct = q_cov_total / q_size_total * 100 if q_size_total else 0.0

    lines = [
        "=" * 60,
        "  Chain Coverage Summary",
        "=" * 60,
        "",
        "Chain files:",
    ]
    for f in chain_files:
        lines.append(f"  {f}")
    lines += [
        "",
        f"{'':>30}  {'sequences':>10}  {'size (bp)':>15}  {'covered (bp)':>15}  {'% covered':>10}",
        "-" * 90,
        f"{'TARGET (T2T)':>30}  {len(target_rows):>10}  {t_size_total:>15,}  {t_cov_total:>15,}  {t_pct:>9.2f}%",
        f"{'QUERY (curated assembly)':>30}  {len(query_rows):>10}  {q_size_total:>15,}  {q_cov_total:>15,}  {q_pct:>9.2f}%",
        "-" * 90,
        "",
        "Per-sequence details written to separate TSV files.",
    ]

    with open(out_path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

    # Also print to stdout
    print('\n'.join(lines))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Compute coverage statistics from UCSC chain files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--chain', nargs='+', required=True, metavar='FILE',
                        help='One or more chain files to analyse')
    parser.add_argument('--target-sizes', metavar='FILE',
                        help='Two-column TSV with target (T2T) sequence lengths')
    parser.add_argument('--query-sizes', metavar='FILE',
                        help='Two-column TSV with query (curated) sequence lengths')
    parser.add_argument('--output-prefix', default='chain_cov', metavar='STR',
                        help='Prefix for output files (default: chain_cov)')
    parser.add_argument('--min-score', type=int, default=0, metavar='INT',
                        help='Skip chains with score below this value (default: 0)')
    args = parser.parse_args()

    print("chain_coverage.py — Computing coverage from chain files", file=sys.stderr)
    print(f"  Output prefix : {args.output_prefix}", file=sys.stderr)
    print(f"  Min chain score: {args.min_score}", file=sys.stderr)

    # ---- Accumulate intervals ------------------------------------------------
    target_ivs, target_sizes_hdr, query_ivs, query_sizes_hdr = \
        accumulate_coverage(args.chain, min_score=args.min_score)

    # ---- Load external size files (override header-derived sizes) -----------
    if args.target_sizes:
        target_sizes = load_sizes(args.target_sizes)
        print(f"  Target sizes loaded from: {args.target_sizes} ({len(target_sizes)} seqs)",
              file=sys.stderr)
    else:
        target_sizes = target_sizes_hdr

    if args.query_sizes:
        query_sizes = load_sizes(args.query_sizes)
        print(f"  Query sizes loaded from: {args.query_sizes} ({len(query_sizes)} seqs)",
              file=sys.stderr)
    else:
        query_sizes = query_sizes_hdr

    # ---- Write per-sequence tables ------------------------------------------
    target_out = f"{args.output_prefix}.target_cov.tsv"
    query_out  = f"{args.output_prefix}.query_cov.tsv"

    print(f"\nWriting target coverage → {target_out}", file=sys.stderr)
    target_rows = write_coverage_table(target_ivs, target_sizes, target_out, "TARGET (T2T)")

    print(f"Writing query coverage  → {query_out}", file=sys.stderr)
    query_rows = write_coverage_table(query_ivs, query_sizes, query_out,
                                      "QUERY (curated assembly)")

    # ---- Summary ------------------------------------------------------------
    target_total = sum(target_sizes.values()) if target_sizes else 0
    query_total  = sum(query_sizes.values())  if query_sizes  else 0

    summary_out = f"{args.output_prefix}.summary.txt"
    print(f"Writing summary         → {summary_out}\n", file=sys.stderr)
    write_summary(target_rows, query_rows, target_total, query_total,
                  summary_out, args.chain)


if __name__ == '__main__':
    main()
