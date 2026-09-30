#!/usr/bin/env python3
"""
validate_switch_blocks_homology.py

Inspects phase switch error blocks detected in an assembly and evaluates their
alignment at exact base-pair resolution and alignment scores by parsing chain alignment blocks:
  1. The Assigned T2T Reference Chromosome (assigned_score, assigned_aligned_bp, assigned_cov_pct)
  2. The Homologous T2T Reference Chromosome (homolog_score, homolog_aligned_bp, homolog_cov_pct)
  3. The Global Best T2T Target across the entire genome (global_best_score, global_best_aligned_bp)

Outputs comprehensive TSV reports to results/07_inspection_hapmers/.
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
    """Categorize into macro, micro, or dot chromosome."""
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


def load_best_chrom_pairs(pairs_path):
    """
    Loads best_chrom_pairs.tsv and builds:
      - asm_to_t2t: scaffold -> assigned T2T chromosome
      - t2t_to_homolog: t2t_chrom -> homologous partner t2t_chrom (same token, opposite side)
    """
    asm_to_t2t = {}
    all_t2t_chroms = []
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
                    all_t2t_chroms.append(t2t_c)

    token_map = defaultdict(list)
    for c in set(all_t2t_chroms):
        tok = chrom_token(c)
        token_map[tok].append(c)

    t2t_to_homolog = {}
    for tok, chrom_list in token_map.items():
        mat_c = [c for c in chrom_list if 'Mat_' in c or 'mat' in c.lower()]
        pat_c = [c for c in chrom_list if 'Pat_' in c or 'pat' in c.lower()]
        if mat_c and pat_c:
            for m in mat_c:
                t2t_to_homolog[m] = pat_c[0]
            for p in pat_c:
                t2t_to_homolog[p] = mat_c[0]

    return asm_to_t2t, t2t_to_homolog


def load_switch_blocks(switch_file):
    """
    Loads switch blocks from TSV or BED file.
    """
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


def parse_chain_blocks_for_switch_blocks(chain_path, switch_blocks):
    """
    Parses a .chain file block by block (size dt dq).
    Calculates exact base pairs aligned and highest chain score for each switch block.
    """
    blocks_by_scaf = defaultdict(list)
    for idx, sb in enumerate(switch_blocks):
        blocks_by_scaf[sb['scaffold']].append((idx, sb['start'], sb['end']))

    # block_idx -> target_chrom -> {'aligned_bp': int, 'best_score': int, 't_min': int, 't_max': int, 'strand': str}
    block_target_stats = defaultdict(lambda: defaultdict(lambda: {
        'aligned_bp': 0, 'best_score': 0, 't_min': float('inf'), 't_max': float('-inf'),
        'q_min': float('inf'), 'q_max': float('-inf'), 'strand': '+'
    }))

    if not os.path.exists(chain_path):
        return block_target_stats

    with open(chain_path) as fh:
        cur_q_in_targets = False
        target_sb_list = []
        t_pos, q_pos = 0, 0
        q_size = 0
        q_strand = '+'
        t_name = ''
        chain_score = 0

        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue

            if line.startswith('chain'):
                p = line.split()
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                chain_score = int(p[1])
                t_name = p[2]
                q_name = p[7]
                q_size = int(p[8])
                q_strand = p[9]
                t_pos = int(p[5])
                q_pos = int(p[10])

                if q_name in blocks_by_scaf:
                    cur_q_in_targets = True
                    target_sb_list = blocks_by_scaf[q_name]
                else:
                    cur_q_in_targets = False
                continue

            if cur_q_in_targets:
                p = line.split()
                size = int(p[0])
                dt = int(p[1]) if len(p) > 1 else 0
                dq = int(p[2]) if len(p) > 2 else 0

                # Block coordinates on query
                if q_strand == '-':
                    qs_block = q_size - (q_pos + size)
                    qe_block = q_size - q_pos
                else:
                    qs_block = q_pos
                    qe_block = q_pos + size

                ts_block = t_pos
                te_block = t_pos + size

                # Check overlap with any switch block on this query
                for b_idx, bs, be in target_sb_list:
                    ov_s = max(qs_block, bs)
                    ov_e = min(qe_block, be)
                    if ov_s < ov_e:
                        ov_len = ov_e - ov_s
                        st = block_target_stats[b_idx][t_name]
                        st['aligned_bp'] += ov_len
                        st['best_score'] = max(st['best_score'], chain_score)
                        st['t_min'] = min(st['t_min'], ts_block)
                        st['t_max'] = max(st['t_max'], te_block)
                        st['q_min'] = min(st['q_min'], ov_s)
                        st['q_max'] = max(st['q_max'], ov_e)
                        st['strand'] = q_strand

                t_pos += size + dt
                q_pos += size + dq

    return block_target_stats


def evaluate_assembly(asm_name, switch_file, pairs_file, chain_file):
    """Evaluates homology status for all switch blocks in a single assembly."""
    asm_to_t2t, t2t_to_homolog = load_best_chrom_pairs(pairs_file)
    blocks = load_switch_blocks(switch_file)
    target_stats_by_block = parse_chain_blocks_for_switch_blocks(chain_file, blocks)

    results = []
    for idx, b in enumerate(blocks):
        scaf = b['scaffold']
        bs = b['start']
        be = b['end']
        blen = b['length']
        hap_str = b['hapmer_info']
        hap_db, exp_cnt, opp_cnt, dens = parse_hapmer_info(hap_str)

        assigned_t2t = asm_to_t2t.get(scaf, 'Unassigned')
        homolog_t2t = t2t_to_homolog.get(assigned_t2t, 'No_Homolog')
        cat = classify_category(assigned_t2t if assigned_t2t != 'Unassigned' else scaf)

        t_stats = target_stats_by_block.get(idx, {})

        # 1. Assigned Target Stats
        if assigned_t2t in t_stats:
            ass_bp = t_stats[assigned_t2t]['aligned_bp']
            ass_score = t_stats[assigned_t2t]['best_score']
            ass_pct = (ass_bp / blen) * 100.0
            st = t_stats[assigned_t2t]
            ass_coords = f"{assigned_t2t}:{st['t_min']:,}-{st['t_max']:,}({st['strand']})"
        else:
            ass_bp = 0
            ass_score = 0
            ass_pct = 0.0
            ass_coords = 'None'

        # 2. Homolog Target Stats
        if homolog_t2t in t_stats:
            hom_bp = t_stats[homolog_t2t]['aligned_bp']
            hom_score = t_stats[homolog_t2t]['best_score']
            hom_pct = (hom_bp / blen) * 100.0
            st = t_stats[homolog_t2t]
            hom_coords = f"{homolog_t2t}:{st['t_min']:,}-{st['t_max']:,}({st['strand']})"
        else:
            hom_bp = 0
            hom_score = 0
            hom_pct = 0.0
            hom_coords = 'None'

        # 3. Global Best Target Stats
        if t_stats:
            # Rank primarily by aligned base pairs (and break ties with chain score)
            glob_tName = max(t_stats.keys(), key=lambda k: (t_stats[k]['aligned_bp'], t_stats[k]['best_score']))
            glob_bp = t_stats[glob_tName]['aligned_bp']
            glob_score = t_stats[glob_tName]['best_score']
            glob_pct = (glob_bp / blen) * 100.0
            st = t_stats[glob_tName]
            glob_coords = f"{glob_tName}:{st['t_min']:,}-{st['t_max']:,}({st['strand']})"
        else:
            glob_tName = 'None'
            glob_bp = 0
            glob_score = 0
            glob_pct = 0.0
            glob_coords = 'None'

        # 4. Status Determination:
        # A switch block is a true homolog switch if:
        # - Homolog aligned bp > assigned aligned bp, OR
        # - Homolog score > assigned score with >= 90% coverage on homolog, OR
        # - Global best target is the homologous chromosome
        if (hom_bp > ass_bp or (hom_score > ass_score and hom_pct >= 90.0) or glob_tName == homolog_t2t) and hom_pct >= 5.0:
            status = 'TRUE_HOMOLOG_SWITCH'
        elif ass_bp >= hom_bp and ass_pct >= 5.0 and glob_tName == assigned_t2t:
            status = 'SAME_CHROMOSOME'
        elif glob_pct >= 5.0 and glob_tName not in (assigned_t2t, homolog_t2t):
            status = 'OFF_TARGET_ALIGNMENT'
        elif glob_bp > 0:
            status = 'REPETITIVE_ALIGNMENT'
        else:
            status = 'UNALIGNED'

        results.append({
            'assembly': asm_name,
            'category': cat,
            'scaffold': scaf,
            'start': bs,
            'end': be,
            'length_bp': blen,
            'hapmer_db': hap_db,
            'expected_hapmer_cnt': exp_cnt,
            'opposite_hapmer_cnt': opp_cnt,
            'opposite_density_pct': dens,
            'assigned_t2t_chrom': assigned_t2t,
            'assigned_score': ass_score,
            'assigned_aligned_bp': ass_bp,
            'assigned_cov_pct': round(ass_pct, 2),
            'assigned_target_coords': ass_coords,
            'homolog_t2t_chrom': homolog_t2t,
            'homolog_score': hom_score,
            'homolog_aligned_bp': hom_bp,
            'homolog_cov_pct': round(hom_pct, 2),
            'homolog_target_coords': hom_coords,
            'global_best_t2t_chrom': glob_tName,
            'global_best_score': glob_score,
            'global_best_aligned_bp': glob_bp,
            'global_best_cov_pct': round(glob_pct, 2),
            'global_best_coords': glob_coords,
            'validation_status': status
        })
    return results


def main():
    parser = argparse.ArgumentParser(description="Validate switch blocks against T2T homologous chromosomes.")
    parser.add_argument('--base-dir', default='/lustre/fs5/vgl/scratch/eduarte/curations/birds/bTaeGut7/redo_curation_analysis',
                        help='Base repository directory')
    parser.add_argument('--out-dir', default=None, help='Output directory for validation TSVs')
    args = parser.parse_args()

    base = args.base_dir
    out_dir = args.out_dir or os.path.join(base, 'results/inspection/07_inspection_hapmers')
    os.makedirs(out_dir, exist_ok=True)

    configs = [
        {
            'name': 'Single',
            'switch_file': f'{base}/results/single/05_hapmers_single/single.switch_blocks.final.tsv',
            'pairs_file': f'{base}/results/single/02_chainpipeline_single/02_chainpipeline_single_strictLinearGap/t2t.vs.single_strictGap.best_chrom_pairs.tsv',
            'chain_file': f'{base}/results/single/02_chainpipeline_single/02_chainpipeline_single_strictLinearGap/t2t.vs.single_strictGap.T2T.vs.ASM.chain'
        },
        {
            'name': 'Dual',
            'switch_file': f'{base}/results/dual/05_hapmers_dual/dual.switch_blocks.final.tsv',
            'pairs_file': f'{base}/results/dual/02_chainpipeline_dual/02_chainpipeline_dual_strictLinearGap/t2t.vs.dual_strictGap.best_chrom_pairs.tsv',
            'chain_file': f'{base}/results/dual/02_chainpipeline_dual/02_chainpipeline_dual_strictLinearGap/t2t.vs.dual_strictGap.T2T.vs.ASM.chain'
        }
    ]

    all_results = []
    for cfg in configs:
        if os.path.exists(cfg['switch_file']):
            print(f"Validating switch blocks for {cfg['name']}...")
            res = evaluate_assembly(cfg['name'], cfg['switch_file'], cfg['pairs_file'], cfg['chain_file'])
            all_results.extend(res)
            df_asm = pd.DataFrame(res)
            asm_out = os.path.join(out_dir, f"{cfg['name'].lower()}_switch_blocks_homology_validation.tsv")
            df_asm.to_csv(asm_out, sep='\t', index=False)
            print(f"  -> Saved: {asm_out} ({len(res)} switch blocks)")

    if all_results:
        df_all = pd.DataFrame(all_results)
        comb_out = os.path.join(out_dir, 'switch_blocks_homology_validation_combined.tsv')
        df_all.to_csv(comb_out, sep='\t', index=False)
        print(f"\nSaved combined validation TSV: {comb_out}")

        print("\n" + "=" * 95)
        print("SUMMARY: SWITCH BLOCKS HOMOLOGY VALIDATION BY CATEGORY & STATUS")
        print("=" * 95)
        summary = df_all.groupby(['assembly', 'category', 'validation_status']).size().unstack(fill_value=0)
        print(summary.to_string())
        print("=" * 95)


if __name__ == '__main__':
    main()
