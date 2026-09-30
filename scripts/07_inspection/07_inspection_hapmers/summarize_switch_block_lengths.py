#!/usr/bin/env python3
"""
summarize_switch_block_lengths.py

Computes the total sum of switch block lengths (in bp), count of switch blocks,
and count of opposite hap-mer markers grouped by:
  - Assembly (Single, Dual, ONT Dual)
  - Reference Haplotype (Maternal, Paternal, and Combined Total)
  - Chromosome Category (Macro, Micro, Dot, and Genome Total)

Uses the corresponding T2T reference chromosome from best_chrom_pairs.tsv
to establish the true reference haplotype and chromosome category.
Outputs TSV and summary reports to results/inspection/07_inspection_hapmers/.
"""

import os
import re
import sys
import argparse
from collections import defaultdict
import pandas as pd

DOT_CHROMS = {'16', '25', '29', '30', '31', '32', '33', '34', '35', '36', '37'}


def chrom_token(name):
    """Extract standard chromosome token (e.g. 1, 1A, 4, 30, Z, W) from name."""
    m = re.search(r'(?:chromosome_|SUPER_)([0-9]+[A-Za-z]*|[A-Za-z]+)', name)
    if m:
        return m.group(1)
    m = re.match(r'^chr([0-9]+[A-Za-z]*|[A-Za-z]+)(?:_(?:mat|pat))?$', name)
    if m:
        return m.group(1)
    return name


def classify_category(name):
    """Categorize into macro, micro, or dot chromosome based on T2T chrom name."""
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


def determine_reference_haplotype(t2t_chrom_name, scaffold_name):
    """Determines if the reference chromosome/scaffold is Maternal or Paternal."""
    name_check = t2t_chrom_name if t2t_chrom_name != 'Unassigned' else scaffold_name
    if 'Mat_' in name_check or 'mat' in name_check.lower():
        return 'Maternal'
    elif 'Pat_' in name_check or 'pat' in name_check.lower():
        return 'Paternal'
    return 'Unknown'


def load_best_chrom_pairs(pairs_path):
    """Loads best_chrom_pairs.tsv into a dict: scaffold -> assigned T2T chromosome."""
    asm_to_t2t = {}
    if os.path.exists(pairs_path):
        with open(pairs_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                p = line.split('\t')
                if len(p) >= 2:
                    t2t_c, scaf = p[0], p[1]
                    asm_to_t2t[scaf] = t2t_c
    return asm_to_t2t


def load_switch_blocks(switch_file):
    """Loads switch blocks from TSV or BED file."""
    records = []
    if not os.path.exists(switch_file):
        return records

    with open(switch_file) as fh:
        hdr = None
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            p = line.split('\t')
            if not hdr and ('scaffold' in p[0].lower() or 'chrom' in p[0].lower()):
                hdr = [x.lower() for x in p]
                continue

            if hdr and 'scaffold' in hdr and 'start' in hdr:
                row = dict(zip(hdr, p))
                scaf = row['scaffold']
                s = int(row['start'])
                e = int(row['end'])
                length = int(row.get('length', e - s))
                hap_info = row.get('hapmer_info', '')
                records.append({
                    'scaffold': scaf,
                    'start': s,
                    'end': e,
                    'length': length,
                    'hapmer_info': hap_info
                })
            elif len(p) >= 3:
                scaf = p[0]
                s = int(p[1])
                e = int(p[2])
                length = e - s
                hap_info = '\t'.join(p[3:]) if len(p) > 3 else ''
                records.append({
                    'scaffold': scaf,
                    'start': s,
                    'end': e,
                    'length': length,
                    'hapmer_info': hap_info
                })
    return records


def parse_hapmer_info(hap_str):
    """Parses Merqury hapmer info string."""
    if not hap_str:
        return 'Unknown', 0, 0, 0.0
    parts = hap_str.strip().split()
    db = parts[0] if len(parts) > 0 else 'Unknown'
    exp_cnt = int(parts[2]) if len(parts) > 2 and parts[2].isdigit() else 0
    opp_cnt = int(parts[3]) if len(parts) > 3 and parts[3].isdigit() else 0
    dens = float(parts[5]) if len(parts) > 5 else 0.0
    return db, exp_cnt, opp_cnt, dens


def main():
    parser = argparse.ArgumentParser(description="Summarize switch block lengths per haplotype and category.")
    parser.add_argument('--base-dir', default='/lustre/fs5/vgl/scratch/eduarte/curations/birds/bTaeGut7/redo_curation_analysis',
                        help='Base repository directory')
    parser.add_argument('--out-dir', default=None, help='Output directory')
    args = parser.parse_args()

    base = args.base_dir
    out_dir = args.out_dir or os.path.join(base, 'results/inspection/07_inspection_hapmers')
    os.makedirs(out_dir, exist_ok=True)

    configs = [
        {
            'name': 'Single',
            'switch_file': f'{base}/results/single/05_hapmers_single/single.switch_blocks.final.tsv',
            'pairs_file': f'{base}/results/single/02_chainpipeline_single/02_chainpipeline_single_strictLinearGap/t2t.vs.single_strictGap.best_chrom_pairs.tsv'
        },
        {
            'name': 'Dual',
            'switch_file': f'{base}/results/dual/05_hapmers_dual/dual.switch_blocks.final.tsv',
            'pairs_file': f'{base}/results/dual/02_chainpipeline_dual/02_chainpipeline_dual_strictLinearGap/t2t.vs.dual_strictGap.best_chrom_pairs.tsv'
        }
    ]

    detailed_rows = []

    for cfg in configs:
        if not os.path.exists(cfg['switch_file']):
            continue

        asm_name = cfg['name']
        asm_to_t2t = load_best_chrom_pairs(cfg['pairs_file'])
        blocks = load_switch_blocks(cfg['switch_file'])

        for b in blocks:
            scaf = b['scaffold']
            bs = b['start']
            be = b['end']
            blen = b['length']
            hap_str = b['hapmer_info']
            hap_db, exp_cnt, opp_cnt, dens = parse_hapmer_info(hap_str)

            assigned_t2t = asm_to_t2t.get(scaf, 'Unassigned')
            haplotype = determine_reference_haplotype(assigned_t2t, scaf)
            category = classify_category(assigned_t2t if assigned_t2t != 'Unassigned' else scaf)
            t2t_token = chrom_token(assigned_t2t if assigned_t2t != 'Unassigned' else scaf)

            detailed_rows.append({
                'assembly': asm_name,
                'haplotype': haplotype,
                'category': category,
                'chrom_token': t2t_token,
                'assigned_t2t_chrom': assigned_t2t,
                'scaffold': scaf,
                'start': bs,
                'end': be,
                'length_bp': blen,
                'hapmer_db': hap_db,
                'expected_hapmer_cnt': exp_cnt,
                'opposite_hapmer_cnt': opp_cnt,
                'opposite_density_pct': dens
            })

    df_details = pd.DataFrame(detailed_rows)
    details_tsv = os.path.join(out_dir, 'switch_blocks_details_with_haplotype_category.tsv')
    df_details.to_csv(details_tsv, sep='\t', index=False)
    print(f"Saved detailed switch blocks TSV: {details_tsv} ({len(df_details)} blocks)")

    # Aggregate by assembly, haplotype, and category
    summary_list = []
    assemblies = df_details['assembly'].unique()
    categories = ['macro', 'micro', 'dot']
    haplotypes = ['Maternal', 'Paternal']

    for asm in assemblies:
        df_a = df_details[df_details['assembly'] == asm]
        for hap in haplotypes:
            df_ah = df_a[df_a['haplotype'] == hap]
            for cat in categories:
                df_ahc = df_ah[df_ah['category'] == cat]
                tot_bp = df_ahc['length_bp'].sum() if len(df_ahc) > 0 else 0
                num_blocks = len(df_ahc)
                tot_opp = df_ahc['opposite_hapmer_cnt'].sum() if len(df_ahc) > 0 else 0
                summary_list.append({
                    'assembly': asm,
                    'haplotype': hap,
                    'category': cat,
                    'num_switch_blocks': num_blocks,
                    'total_switch_length_bp': tot_bp,
                    'total_switch_length_mb': round(tot_bp / 1_000_000.0, 4),
                    'total_opposite_hapmer_markers': tot_opp
                })
            # Subtotal for haplotype
            tot_hap_bp = df_ah['length_bp'].sum() if len(df_ah) > 0 else 0
            summary_list.append({
                'assembly': asm,
                'haplotype': hap,
                'category': 'ALL_CATEGORIES_SUBTOTAL',
                'num_switch_blocks': len(df_ah),
                'total_switch_length_bp': tot_hap_bp,
                'total_switch_length_mb': round(tot_hap_bp / 1_000_000.0, 4),
                'total_opposite_hapmer_markers': df_ah['opposite_hapmer_cnt'].sum() if len(df_ah) > 0 else 0
            })

        # Assembly Total
        tot_asm_bp = df_a['length_bp'].sum() if len(df_a) > 0 else 0
        summary_list.append({
            'assembly': asm,
            'haplotype': 'COMBINED_TOTAL',
            'category': 'GENOME_TOTAL',
            'num_switch_blocks': len(df_a),
            'total_switch_length_bp': tot_asm_bp,
            'total_switch_length_mb': round(tot_asm_bp / 1_000_000.0, 4),
            'total_opposite_hapmer_markers': df_a['opposite_hapmer_cnt'].sum() if len(df_a) > 0 else 0
        })

    df_summary = pd.DataFrame(summary_list)
    summary_tsv = os.path.join(out_dir, 'switch_blocks_summary_by_category_and_haplotype.tsv')
    df_summary.to_csv(summary_tsv, sep='\t', index=False)
    print(f"Saved summary TSV: {summary_tsv}")

    # Generate pretty wide table for display
    print("\n" + "=" * 110)
    print("TOTAL SUM OF SWITCH BLOCK LENGTHS (BP) BY ASSEMBLY, HAPLOTYPE & CATEGORY")
    print("=" * 110)
    pivot_bp = df_summary[df_summary['category'].isin(['macro', 'micro', 'dot'])].pivot_table(
        index=['assembly', 'category'],
        columns='haplotype',
        values='total_switch_length_bp',
        aggfunc='sum',
        fill_value=0
    )
    pivot_bp['Combined_Total_bp'] = pivot_bp.sum(axis=1)
    pivot_bp['Combined_Total_Mb'] = (pivot_bp['Combined_Total_bp'] / 1_000_000.0).round(4)
    print(pivot_bp.to_string())
    print("=" * 110)


if __name__ == '__main__':
    main()
