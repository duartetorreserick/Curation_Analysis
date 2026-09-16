#!/usr/bin/env python3
"""
summarize_stats_by_category_v2.py

Computes per-category genome assembly statistics for Single, Dual, and ONT Dual assemblies.

Categories
----------
  macro   — chromosomes 1–8, Z, W (excluding dot 16, 25)
  micro   — chromosomes 9–28 (excluding dot chromosomes 16, 25)
  dot     — dot chromosomes (16, 25) and chromosomes 29+
  genome  — all chromosomes combined

Metrics
-------
  collinear_pct       collinear chain coverage / total category length × 100
  noncollinear_pct    non-collinear chain coverage (including unlocs if provided) / total category length × 100
  switch_pct          hapmer switch-error bp (lifted to T2T) / total category length × 100
  telomere_pct        telomere ends detected / expected × 100  (2 per chrom, 2n=80 → 160)

Outputs
-------
  PREFIX_wide.tsv   one row per assembly × category
  PREFIX_long.tsv   long-form melted metrics for plotting
"""

import argparse
import os
import re
from collections import defaultdict
import pandas as pd


DOT_CHROMS = {'16', '25', '29', '30', '31', '32', '33', '34', '35', '36', '37'}
CATEGORY_ORDER = ['macro', 'micro', 'dot', 'genome']
ASSEMBLY_ORDER = ['Single', 'Dual', 'ONT_Dual']


def chrom_token(name):
    m = re.search(r'(?:chromosome_|SUPER_)([0-9]+[A-Za-z]*|[A-Za-z]+)', name)
    if m:
        return m.group(1)
    m = re.match(r'^chr([0-9]+[A-Za-z]*|[A-Za-z]+)(?:_(?:mat|pat))?$', name)
    if m:
        return m.group(1)
    return name


def classify(name):
    tok = chrom_token(name)
    if tok.upper() in ('Z', 'W'):
        return 'macro'
    if tok in DOT_CHROMS:
        return 'dot'
    m = re.match(r'^(\d+)', tok)
    if m:
        n = int(m.group(1))
        if n <= 8:
            return 'macro'
        elif n <= 28:
            return 'micro'
        else:
            return 'dot'
    return 'macro'


def _merge(intervals):
    if not intervals:
        return []
    ivs = sorted(intervals)
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def _bp(intervals):
    return sum(e - s for s, e in intervals)


def _parse_chain_target(path):
    if not path or not os.path.exists(path):
        return
    with open(path) as fh:
        hdr = None
        blks = []
        for raw in fh:
            line = raw.rstrip('\n')
            if not line or line.startswith('#'):
                if hdr is not None and not line:
                    yield _emit(hdr, blks)
                    hdr = None
                    blks = []
                continue
            if line.startswith('chain'):
                if hdr is not None:
                    yield _emit(hdr, blks)
                p = line.split()
                hdr = dict(tName=p[2], tSize=int(p[3]), tStart=int(p[5]))
                blks = []
            else:
                parts = line.split()
                if parts:
                    blks.append(tuple(int(x) for x in parts))
        if hdr is not None:
            yield _emit(hdr, blks)


def _emit(hdr, blks):
    pos = hdr['tStart']
    ivs = []
    for blk in blks:
        size = blk[0]
        ivs.append((pos, pos + size))
        pos += size + (blk[1] if len(blk) == 3 else 0)
    return hdr['tName'], hdr['tSize'], ivs


def compute_covered_bp(chain_paths, seq_sizes):
    if isinstance(chain_paths, str):
        chain_paths = [chain_paths]
    covered = defaultdict(list)
    for cp in chain_paths:
        if not cp or not os.path.exists(cp):
            continue
        for tname, _, ivs in _parse_chain_target(cp):
            if tname in seq_sizes:
                covered[tname].extend(ivs)
    return {nm: _bp(_merge(covered.get(nm, []))) for nm in seq_sizes}


def _load_pairs(path):
    pairs = {}
    with open(path) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= 2:
                pairs[p[0]] = p[1]
    return pairs


def _load_switch_bed(path):
    blocks = defaultdict(list)
    if not path or not os.path.exists(path):
        return dict(blocks)
    with open(path) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= 3:
                try:
                    blocks[p[0]].append((int(p[1]), int(p[2])))
                except ValueError:
                    pass
    return dict(blocks)


def _build_liftover(chain_path):
    liftover = defaultdict(lambda: defaultdict(list))

    def _proc(hdr, blks):
        t_pos = hdr['tStart']
        q_pos = hdr['qStart']
        q_strand = hdr['qStrand']
        q_size = hdr['qSize']
        for blk in blks:
            sz = blk[0]
            dt = blk[1] if len(blk) > 1 else 0
            dq = blk[2] if len(blk) > 2 else 0
            qs = (q_size - (q_pos + sz)) if q_strand == '-' else q_pos
            qe = (q_size - q_pos) if q_strand == '-' else q_pos + sz
            liftover[hdr['qName']][hdr['tName']].append((qs, qe, t_pos))
            t_pos += sz + dt
            q_pos += sz + dq

    with open(chain_path) as fh:
        hdr = None
        blks = []
        for raw in fh:
            line = raw.rstrip('\n')
            if not line or line.startswith('#'):
                if hdr is not None and not line:
                    _proc(hdr, blks)
                    hdr = None
                    blks = []
                continue
            if line.startswith('chain'):
                if hdr is not None:
                    _proc(hdr, blks)
                p = line.split()
                hdr = {
                    'tName': p[2],
                    'tSize': int(p[3]),
                    'tStrand': p[4],
                    'tStart': int(p[5]),
                    'qName': p[7],
                    'qSize': int(p[8]),
                    'qStrand': p[9],
                    'qStart': int(p[10]),
                }
                blks = []
            else:
                parts = line.split()
                if parts:
                    blks.append(tuple(int(x) for x in parts))
        if hdr is not None:
            _proc(hdr, blks)

    return {
        qn: {tn: sorted(bs, key=lambda x: x[0]) for tn, bs in td.items()}
        for qn, td in liftover.items()
    }


def compute_switch_bp(pairs_path, bed_path, chain_path, seq_sizes):
    pairs = _load_pairs(pairs_path)
    switches = _load_switch_bed(bed_path)
    liftover = _build_liftover(chain_path)

    result = {nm: 0 for nm in seq_sizes}
    for t2t_name, scaffold in pairs.items():
        raw = switches.get(scaffold, [])
        chain_blocks = liftover.get(scaffold, {}).get(t2t_name, [])
        if not raw or not chain_blocks:
            continue
        t2t_ivs = []
        for h_s, h_e in raw:
            for q_s, q_e, t_start in chain_blocks:
                ov_s = max(h_s, q_s)
                ov_e = min(h_e, q_e)
                if ov_s >= ov_e:
                    continue
                offset = ov_s - q_s
                t2t_ivs.append((t_start + offset, t_start + offset + (ov_e - ov_s)))
        if t2t_ivs:
            result[t2t_name] = _bp(_merge(t2t_ivs))
    return result


def load_telomere_presence_tsv(path):
    """
    Parses modern <asm>_telomere_presence.tsv:
    Columns: chromosome, scaffold, collinear, non-collinear, teloscope_rescued, type_of_tele_rescued
    Returns {chrom: set(arm)} containing all detected arms.
    """
    telo = defaultdict(set)
    if not path or not os.path.exists(path):
        return dict(telo)
    df = pd.read_csv(path, sep='\t')
    for _, row in df.iterrows():
        chrom = str(row['chromosome']).strip()
        coll = str(row.get('collinear', '')).strip().lower()
        nc = str(row.get('non-collinear', '')).strip().lower()
        rescued = str(row.get('teloscope_rescued', '')).strip().lower() == 'yes'
        rescued_type = str(row.get('type_of_tele_rescued', '')).strip().lower()

        for arm in ['p', 'q']:
            if arm in coll or arm in nc:
                telo[chrom].add(arm)
            elif rescued and arm in rescued_type:
                telo[chrom].add(arm)
    return dict(telo)


def load_legacy_hifi_telomeres(path):
    telo_s = defaultdict(set)
    telo_d = defaultdict(set)
    if not path or not os.path.exists(path):
        return dict(telo_s), dict(telo_d)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 5:
                continue
            _, chrom, arm, single_status, dual_status = cols[:5]
            if single_status.strip().lower() == 'ok':
                telo_s[chrom].add(arm)
            if dual_status.strip().lower() == 'ok':
                telo_d[chrom].add(arm)
    return dict(telo_s), dict(telo_d)


def load_legacy_ont_telomeres(path):
    telo = defaultdict(set)
    if not path or not os.path.exists(path):
        return dict(telo)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 4:
                continue
            chrom, arm, nc_status, coll_status = cols[:4]
            if coll_status.strip().lower() == 'ok' or nc_status.strip().lower() == 'ok':
                telo[chrom].add(arm)
    return dict(telo)


def aggregate_coverage(seq_sizes, cov_bp, nc_bp, sw_bp=None):
    agg = defaultdict(lambda: [0, 0, 0, 0])
    for chrom, size in seq_sizes.items():
        cat = classify(chrom)
        agg[cat][0] += size
        agg[cat][1] += cov_bp.get(chrom, 0)
        agg[cat][2] += nc_bp.get(chrom, 0)
        agg[cat][3] += sw_bp.get(chrom, 0) if sw_bp else 0

    genome = [sum(v[i] for v in agg.values()) for i in range(4)]
    agg['genome'] = genome
    return {cat: tuple(v) for cat, v in agg.items()}


def aggregate_telomeres(seq_sizes, telo_dict):
    agg = defaultdict(lambda: [0, 0])
    for chrom in seq_sizes:
        cat = classify(chrom)
        agg[cat][1] += 2  # 2 arms expected per T2T chromosome
        agg[cat][0] += len(telo_dict.get(chrom, set()))

    genome = [sum(v[i] for v in agg.values()) for i in range(2)]
    agg['genome'] = genome
    return {cat: tuple(v) for cat, v in agg.items()}


def main():
    ap = argparse.ArgumentParser(
        description='Per-category genome assembly QC statistics (v2).',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument('--single-tsv', required=True, metavar='FILE')
    ap.add_argument('--dual-tsv', required=True, metavar='FILE')
    ap.add_argument('--single-chain', required=True, metavar='FILE')
    ap.add_argument('--dual-chain', required=True, metavar='FILE')
    ap.add_argument('--single-nc-chain', required=True, metavar='FILE')
    ap.add_argument('--dual-nc-chain', required=True, metavar='FILE')
    ap.add_argument('--single-unlocs-chain', required=False, metavar='FILE', default=None)
    ap.add_argument('--dual-unlocs-chain', required=False, metavar='FILE', default=None)
    ap.add_argument('--ont-dual-chain', required=True, metavar='FILE')
    ap.add_argument('--ont-dual-nc-chain', required=True, metavar='FILE')
    ap.add_argument('--single-pairs', required=True, metavar='FILE')
    ap.add_argument('--dual-pairs', required=True, metavar='FILE')
    ap.add_argument('--single-bed', required=True, metavar='FILE')
    ap.add_argument('--dual-bed', required=True, metavar='FILE')
    ap.add_argument('--single-telomere-presence', required=False, metavar='FILE', default=None)
    ap.add_argument('--dual-telomere-presence', required=False, metavar='FILE', default=None)
    ap.add_argument('--ont-telomere-presence', required=False, metavar='FILE', default=None)
    ap.add_argument('--telo-combined-report', required=False, metavar='FILE', default=None)
    ap.add_argument('--ont-telo-combined-report', required=False, metavar='FILE', default=None)
    ap.add_argument('--single-fai', required=False, metavar='FILE', default=None)
    ap.add_argument('--dual-fai', required=False, metavar='FILE', default=None)
    ap.add_argument('--output', required=True, metavar='PREFIX',
                    help='Output prefix — writes PREFIX_wide.tsv and PREFIX_long.tsv')
    args = ap.parse_args()

    print('Loading chromosome sizes...')
    def _load_tsv(path):
        df = pd.read_csv(path, sep='\t', comment='#')
        df.columns = [c.strip() for c in df.columns]
        return df.set_index('name')['size'].to_dict()

    sizes_s = _load_tsv(args.single_tsv)
    sizes_d = _load_tsv(args.dual_tsv)
    seq_sizes = {**sizes_d, **sizes_s}

    print('Computing collinear coverage...')
    cov_s = compute_covered_bp(args.single_chain, seq_sizes)
    cov_d = compute_covered_bp(args.dual_chain, seq_sizes)
    cov_o = compute_covered_bp(args.ont_dual_chain, seq_sizes)

    print('Computing non-collinear coverage...')
    s_nc_chains = [args.single_nc_chain]
    if args.single_unlocs_chain:
        s_nc_chains.append(args.single_unlocs_chain)
    d_nc_chains = [args.dual_nc_chain]
    if args.dual_unlocs_chain:
        d_nc_chains.append(args.dual_unlocs_chain)

    nc_s = compute_covered_bp(s_nc_chains, seq_sizes)
    nc_d = compute_covered_bp(d_nc_chains, seq_sizes)
    nc_o = compute_covered_bp(args.ont_dual_nc_chain, seq_sizes)

    print('Computing switch-error bp (Single)...')
    sw_s = compute_switch_bp(args.single_pairs, args.single_bed, args.single_chain, seq_sizes)
    print('Computing switch-error bp (Dual)...')
    sw_d = compute_switch_bp(args.dual_pairs, args.dual_bed, args.dual_chain, seq_sizes)
    sw_o = {nm: 0 for nm in seq_sizes}

    agg_s = aggregate_coverage(seq_sizes, cov_s, nc_s, sw_s)
    agg_d = aggregate_coverage(seq_sizes, cov_d, nc_d, sw_d)
    agg_o = aggregate_coverage(seq_sizes, cov_o, nc_o, sw_o)

    print('Loading telomeres...')
    telo_s_dict = {}
    telo_d_dict = {}
    telo_o_dict = {}

    if args.single_telomere_presence and os.path.exists(args.single_telomere_presence):
        telo_s_dict = load_telomere_presence_tsv(args.single_telomere_presence)
    elif args.telo_combined_report:
        telo_s_dict, _ = load_legacy_hifi_telomeres(args.telo_combined_report)

    if args.dual_telomere_presence and os.path.exists(args.dual_telomere_presence):
        telo_d_dict = load_telomere_presence_tsv(args.dual_telomere_presence)
    elif args.telo_combined_report:
        _, telo_d_dict = load_legacy_hifi_telomeres(args.telo_combined_report)

    if args.ont_telomere_presence and os.path.exists(args.ont_telomere_presence):
        telo_o_dict = load_telomere_presence_tsv(args.ont_telomere_presence)
    elif args.ont_telo_combined_report:
        telo_o_dict = load_legacy_ont_telomeres(args.ont_telo_combined_report)

    telo_agg_s = aggregate_telomeres(seq_sizes, telo_s_dict)
    telo_agg_d = aggregate_telomeres(seq_sizes, telo_d_dict)
    telo_agg_o = aggregate_telomeres(seq_sizes, telo_o_dict)

    n_chroms = defaultdict(int)
    for chrom in seq_sizes:
        n_chroms[classify(chrom)] += 1
    n_chroms['genome'] = sum(n_chroms[c] for c in CATEGORY_ORDER[:-1])

    wide_rows = []
    specs = [
        ('Single', agg_s, telo_agg_s, True),
        ('Dual', agg_d, telo_agg_d, True),
        ('ONT_Dual', agg_o, telo_agg_o, False),
    ]

    for asm, agg, telo_agg, has_switch in specs:
        for cat in CATEGORY_ORDER:
            if cat not in agg:
                continue
            total_bp, coll_bp, nc_bp, sw_bp = agg[cat]
            telo_present, telo_expected = telo_agg.get(cat, (0, 0))

            coll_pct = 100.0 * coll_bp / total_bp if total_bp else 0.0
            nc_pct = 100.0 * nc_bp / total_bp if total_bp else 0.0
            sw_pct = 100.0 * sw_bp / total_bp if (total_bp and has_switch) else 0.0
            telo_pct = 100.0 * telo_present / telo_expected if telo_expected else 0.0

            wide_rows.append({
                'assembly': asm,
                'category': cat,
                'n_chroms': n_chroms[cat],
                'total_bp': total_bp,
                'collinear_bp': coll_bp,
                'noncollinear_bp': nc_bp,
                'switch_bp': sw_bp if has_switch else 0,
                'collinear_pct': round(coll_pct, 4),
                'noncollinear_pct': round(nc_pct, 4),
                'switch_pct': round(sw_pct, 4),
                'n_telo_present': telo_present,
                'n_telo_expected': telo_expected,
                'telomere_pct': round(telo_pct, 4),
            })

    wide_df = pd.DataFrame(wide_rows)

    long_df = wide_df.melt(
        id_vars=['assembly', 'category'],
        value_vars=['collinear_pct', 'noncollinear_pct', 'switch_pct', 'telomere_pct'],
        var_name='metric',
        value_name='value',
    )
    long_df['metric'] = long_df['metric'].str.replace('_pct', '', regex=False)
    long_df['assembly'] = pd.Categorical(long_df['assembly'], categories=ASSEMBLY_ORDER, ordered=True)
    long_df['category'] = pd.Categorical(long_df['category'], categories=CATEGORY_ORDER, ordered=True)
    long_df = long_df.sort_values(['assembly', 'category', 'metric']).reset_index(drop=True)

    out_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(out_dir, exist_ok=True)

    wide_path = f'{args.output}_wide.tsv'
    long_path = f'{args.output}_long.tsv'
    wide_df.to_csv(wide_path, sep='\t', index=False, float_format='%.4f')
    long_df.to_csv(long_path, sep='\t', index=False, float_format='%.4f')

    print(f'\nSaved: {wide_path}')
    print(f'Saved: {long_path}')

    print('\n=== Coverage and telomere statistics by category ===\n')
    display_cols = ['assembly', 'category', 'n_chroms', 'total_bp',
                    'collinear_pct', 'noncollinear_pct', 'switch_pct', 'telomere_pct']
    print(wide_df[display_cols].to_string(index=False, float_format=lambda x: f'{x:.2f}'))


if __name__ == '__main__':
    main()
