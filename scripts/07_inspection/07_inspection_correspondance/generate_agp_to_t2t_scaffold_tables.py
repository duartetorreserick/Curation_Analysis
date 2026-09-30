#!/usr/bin/env python3
"""
generate_agp_to_t2t_scaffold_tables.py

Generates cross-assembly mapping tables connecting:
  1. AGP Painted Scaffolds (Scaffold_1, Scaffold_2, ...)
  2. AGP Component Scaffolds / Scaffolds H1/H2 (scaffold_1.H1, scaffold_35.H2, ...)
  3. Pretext Internal SUPER name (SUPER_internal)
  4. Standardized / Harmonized SUPER name (SUPER_renamed)
  5. Final Curated Chromosome name (combined_SUPER, e.g. Mat_SUPER_1.H1, Mat_SUPER_2)
  6. Corresponding T2T reference chromosome (e.g. chr1, chr2, chrZ, chrW)
  7. Exact T2T sequence header accession (e.g. Mat_NC_133026.1_chromosome_2)

Outputs:
  - Component-level mapping table (every individual component block in AGP):
    results/inspection/02_scaffold_agp_to_t2t/agp_scaffold_to_t2t_all_assemblies.tsv (.csv)
  - Painted scaffold summary table (consolidated per curated chromosome):
    results/inspection/02_scaffold_agp_to_t2t/agp_painted_scaffold_summary.tsv (.csv)

Usage:
  python3 scripts/07_inspection/generate_agp_to_t2t_scaffold_tables.py
"""

import os
import sys
import csv
import re
from collections import OrderedDict


def parse_t2t_name(t2t_seq):
    if not t2t_seq or t2t_seq == "N/A":
        return "N/A", "N/A"
    m = re.search(r'chromosome_(.+)$', t2t_seq)
    chrom = f"chr{m.group(1)}" if m else t2t_seq
    return chrom, t2t_seq


def main():
    base = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
    out_dir = os.path.join(base, "results/inspection/02_scaffold_agp_to_t2t")
    os.makedirs(out_dir, exist_ok=True)

    # 1. Dual best chrom pairs & scaffold map
    dual_best = {}
    dual_best_path = os.path.join(base, "results/dual/02_chainpipeline_dual/t2t.vs.dual.best_chrom_pairs.tsv")
    with open(dual_best_path) as f:
        for line in f:
            if line.strip():
                t2t, sup = line.strip().split('\t')
                dual_best[sup] = t2t
    # Dual Pat_SUPER_37.H2 biological reference pairing
    dual_best['Pat_SUPER_37.H2'] = 'Pat_CM109787.1_chromosome_37'

    dual_scaf_map = {}
    with open(os.path.join(base, "data/dual/scaffold_hap_to_combined_super.tsv")) as f:
        header = next(f).strip().split('\t')
        for line in f:
            if line.strip():
                row = dict(zip(header, line.strip().split('\t')))
                dual_scaf_map[row['scaffold']] = row

    # 2. Single best chrom pairs & scaffold map
    single_best = {}
    single_best_path = os.path.join(base, "results/single/02_chainpipeline_single/t2t.vs.single.best_chrom_pairs.tsv")
    with open(single_best_path) as f:
        for line in f:
            if line.strip():
                t2t, sup = line.strip().split('\t')
                single_best[sup] = t2t

    single_scaf_map = {}
    with open(os.path.join(base, "data/single/scaffold_hap_to_combined_super.tsv")) as f:
        header = next(f).strip().split('\t')
        for line in f:
            if line.strip():
                row = dict(zip(header, line.strip().split('\t')))
                single_scaf_map[(row['scaffold'], row['haplotype'])] = row

    # 3. ONT best chrom pairs & scaffold map
    ont_best = {}
    ont_best_path = os.path.join(base, "results/ont/02_chainpipeline_ont/t2t.vs.ont.best_chrom_pairs.tsv")
    with open(ont_best_path) as f:
        for line in f:
            if line.strip():
                t2t, sup = line.strip().split('\t')
                ont_best[sup] = t2t

    ont_scaf_map = {}
    with open(os.path.join(base, "data/ont/scaffold_hap_to_combined_super.tsv")) as f:
        header = next(f).strip().split('\t')
        for line in f:
            if line.strip():
                row = dict(zip(header, line.strip().split('\t')))
                ont_scaf_map[row['scaffold']] = row

    all_records = []
    painted_summaries = []

    # Process Dual
    dual_agp_path = os.path.join(base, "data/dual/corrected.agp")
    dual_painted_comps = OrderedDict()
    with open(dual_agp_path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            p = line.strip().split('\t')
            scaf = p[0]
            if any('Painted' in x for x in p) and p[4] == 'W':
                comp = p[5]
                dual_painted_comps.setdefault(scaf, []).append((comp, p[1], p[2], p[8]))

    for scaf, comps in dual_painted_comps.items():
        s_info = dual_scaf_map.get(scaf, {})
        final_super = s_info.get('combined_SUPER', 'N/A')
        s_int = s_info.get('SUPER_internal', 'N/A')
        s_ren = s_info.get('SUPER_renamed', 'N/A')
        t2t_raw = dual_best.get(final_super, 'N/A')
        t2t_chr, t2t_acc = parse_t2t_name(t2t_raw)
        hap = s_info.get('haplotype', 'N/A')
        mp = s_info.get('mat_pat', 'N/A')

        unique_comps = []
        for c in comps:
            if c[0] not in unique_comps:
                unique_comps.append(c[0])

        painted_summaries.append({
            'Assembly': 'Dual',
            'Haplotype': hap,
            'Parentage': mp,
            'Painted_Scaffold': scaf,
            'SUPER_internal': s_int,
            'SUPER_renamed': s_ren,
            'Final_Super_Name': final_super,
            'T2T_Chromosome': t2t_chr,
            'T2T_Accession': t2t_acc,
            'Scaffold_H1_H2_List': ", ".join(unique_comps),
            'Num_Constituent_Scaffolds': len(unique_comps)
        })

        for c, beg, end, ori in comps:
            all_records.append({
                'Assembly': 'Dual',
                'Haplotype': hap,
                'Parentage': mp,
                'Painted_Scaffold': scaf,
                'Scaffold_H1_H2': c,
                'Component_Start': beg,
                'Component_End': end,
                'Orientation': ori,
                'SUPER_internal': s_int,
                'SUPER_renamed': s_ren,
                'Final_Super_Name': final_super,
                'T2T_Chromosome': t2t_chr,
                'T2T_Accession': t2t_acc
            })

    # Process ONT
    ont_agp_path = os.path.join(base, "data/ont/corrected.agp")
    ont_painted_comps = OrderedDict()
    with open(ont_agp_path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            p = line.strip().split('\t')
            scaf = p[0]
            if any('Painted' in x for x in p) and p[4] == 'W':
                comp = p[5]
                ont_painted_comps.setdefault(scaf, []).append((comp, p[1], p[2], p[8]))

    for scaf, comps in ont_painted_comps.items():
        s_info = ont_scaf_map.get(scaf, {})
        final_super = s_info.get('combined_SUPER', 'N/A')
        s_int = s_info.get('SUPER_internal', 'N/A')
        s_ren = s_info.get('SUPER_renamed', 'N/A')
        t2t_raw = ont_best.get(final_super, 'N/A')
        t2t_chr, t2t_acc = parse_t2t_name(t2t_raw)
        hap = s_info.get('haplotype', 'N/A')
        mp = s_info.get('mat_pat', 'N/A')

        unique_comps = []
        for c in comps:
            if c[0] not in unique_comps:
                unique_comps.append(c[0])

        painted_summaries.append({
            'Assembly': 'ONT',
            'Haplotype': hap,
            'Parentage': mp,
            'Painted_Scaffold': scaf,
            'SUPER_internal': s_int,
            'SUPER_renamed': s_ren,
            'Final_Super_Name': final_super,
            'T2T_Chromosome': t2t_chr,
            'T2T_Accession': t2t_acc,
            'Scaffold_H1_H2_List': ", ".join(unique_comps),
            'Num_Constituent_Scaffolds': len(unique_comps)
        })

        for c, beg, end, ori in comps:
            all_records.append({
                'Assembly': 'ONT',
                'Haplotype': hap,
                'Parentage': mp,
                'Painted_Scaffold': scaf,
                'Scaffold_H1_H2': c,
                'Component_Start': beg,
                'Component_End': end,
                'Orientation': ori,
                'SUPER_internal': s_int,
                'SUPER_renamed': s_ren,
                'Final_Super_Name': final_super,
                'T2T_Chromosome': t2t_chr,
                'T2T_Accession': t2t_acc
            })

    # Process Single
    for h_tag, agp_file in [('H1', 'data/single/Hap1.corrected.agp'), ('H2', 'data/single/Hap2.corrected.agp')]:
        s_painted_comps = OrderedDict()
        with open(os.path.join(base, agp_file)) as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                p = line.strip().split('\t')
                scaf = p[0]
                if any('Painted' in x for x in p) and p[4] == 'W':
                    comp = p[5]
                    s_painted_comps.setdefault(scaf, []).append((comp, p[1], p[2], p[8]))

        for scaf, comps in s_painted_comps.items():
            s_info = single_scaf_map.get((scaf, h_tag), {})
            final_super = s_info.get('combined_SUPER', 'N/A')
            s_int = s_info.get('SUPER_internal', 'N/A')
            s_ren = s_info.get('SUPER_renamed', 'N/A')
            t2t_raw = single_best.get(final_super, 'N/A')
            t2t_chr, t2t_acc = parse_t2t_name(t2t_raw)
            mp = s_info.get('mat_pat', 'N/A')

            unique_comps = []
            for c in comps:
                if c[0] not in unique_comps:
                    unique_comps.append(c[0])

            painted_summaries.append({
                'Assembly': 'Single',
                'Haplotype': h_tag,
                'Parentage': mp,
                'Painted_Scaffold': scaf,
                'SUPER_internal': s_int,
                'SUPER_renamed': s_ren,
                'Final_Super_Name': final_super,
                'T2T_Chromosome': t2t_chr,
                'T2T_Accession': t2t_acc,
                'Scaffold_H1_H2_List': ", ".join(unique_comps),
                'Num_Constituent_Scaffolds': len(unique_comps)
            })

            for c, beg, end, ori in comps:
                all_records.append({
                    'Assembly': 'Single',
                    'Haplotype': h_tag,
                    'Parentage': mp,
                    'Painted_Scaffold': scaf,
                    'Scaffold_H1_H2': c,
                    'Component_Start': beg,
                    'Component_End': end,
                    'Orientation': ori,
                    'SUPER_internal': s_int,
                    'SUPER_renamed': s_ren,
                    'Final_Super_Name': final_super,
                    'T2T_Chromosome': t2t_chr,
                    'T2T_Accession': t2t_acc
                })

    rec_fields = [
        'Assembly', 'Haplotype', 'Parentage', 'Painted_Scaffold', 'Scaffold_H1_H2',
        'Component_Start', 'Component_End', 'Orientation', 'SUPER_internal',
        'SUPER_renamed', 'Final_Super_Name', 'T2T_Chromosome', 'T2T_Accession'
    ]
    tsv_rec_path = os.path.join(out_dir, "agp_scaffold_to_t2t_all_assemblies.tsv")
    csv_rec_path = os.path.join(out_dir, "agp_scaffold_to_t2t_all_assemblies.csv")

    with open(tsv_rec_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rec_fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(all_records)

    with open(csv_rec_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rec_fields)
        writer.writeheader()
        writer.writerows(all_records)

    sum_fields = [
        'Assembly', 'Haplotype', 'Parentage', 'Painted_Scaffold', 'SUPER_internal',
        'SUPER_renamed', 'Final_Super_Name', 'T2T_Chromosome', 'T2T_Accession',
        'Scaffold_H1_H2_List', 'Num_Constituent_Scaffolds'
    ]
    tsv_sum_path = os.path.join(out_dir, "agp_painted_scaffold_summary.tsv")
    csv_sum_path = os.path.join(out_dir, "agp_painted_scaffold_summary.csv")

    with open(tsv_sum_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=sum_fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(painted_summaries)

    with open(csv_sum_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=sum_fields)
        writer.writeheader()
        writer.writerows(painted_summaries)

    print(f"Successfully generated mapping tables in {out_dir}:")
    print(f"  - Detailed component table: {tsv_rec_path} ({len(all_records)} records)")
    print(f"  - Painted scaffold summary: {tsv_sum_path} ({len(painted_summaries)} records)")


if __name__ == '__main__':
    main()
