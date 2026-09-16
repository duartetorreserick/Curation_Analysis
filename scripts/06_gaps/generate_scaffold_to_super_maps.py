#!/usr/bin/env python3
"""
generate_scaffold_to_super_maps.py

Generates the authoritative scaffold-to-chromosome mapping table
(scaffold_hap_to_combined_super.tsv) for Dual, Single, or ONT assemblies.

Maps each scaffold from the curated AGP (Pretext output) to its:
  1. combined_SUPER: Final chromosome name in the chain/final FASTA
  2. mat_pat: Maternal or Paternal assignment
  3. haplotype: H1 or H2
  4. scaffold: AGP scaffold identifier (e.g. Scaffold_1)
  5. SUPER_internal: Internal SUPER name assigned during Pretext curation
  6. SUPER_renamed: Standard reference chromosome SUPER name

Usage Examples:
  # Explicit arguments for DUAL:
  python3 scripts/06_gaps/generate_scaffold_to_super_maps.py \
      --mode dual \
      --hap1-inter data/dual/hap1.inter_chr.tsv \
      --hap2-inter data/dual/hap2.inter_chr.tsv \
      --hap2-vs-hap1 data/dual/hap2.vs.hap1.tsv \
      --hap-assignment results/dual/01_formatting_assemblies_dual/hap_assignment.tsv \
      --fai results/dual/02_chainpipeline_dual/dual_combined.renamed.sorted.reoriented.fa.fai \
      --out data/dual/scaffold_hap_to_combined_super.tsv

  # Explicit arguments for SINGLE:
  python3 scripts/06_gaps/generate_scaffold_to_super_maps.py \
      --mode single \
      --hap1-inter data/single/hap1.inter_chr.tsv \
      --hap2-inter data/single/hap2.inter_chr.tsv \
      --hap2-vs-hap1 data/single/hap2.vs.hap1.tsv \
      --hap-assignment results/single/01_formatting_assemblies_single/hap_assignment.tsv \
      --fai results/single/02_chainpipeline_single/single_combined.renamed.sorted.reoriented.fa.fai \
      --out data/single/scaffold_hap_to_combined_super.tsv

  # Explicit arguments for ONT:
  python3 scripts/06_gaps/generate_scaffold_to_super_maps.py \
      --mode ont \
      --hap1-inter data/ont/hap1.inter_chr.tsv \
      --hap2-inter data/ont/hap2.inter_chr.tsv \
      --new-names data/ont/new_name_change.tsv \
      --fai results/ont/02_chainpipeline_ont/asm3_ONT_combined.sorted.reoriented.fa.fai \
      --out data/ont/scaffold_hap_to_combined_super.tsv
"""

import os
import sys
import argparse


def read_tsv_dict(filepath, key_col=0, val_col=1, skip_header=False):
    """Read a two-column mapping TSV into a dict."""
    mapping = {}
    if not filepath or not os.path.exists(filepath):
        raise FileNotFoundError(f"Missing required file: {filepath}")
    with open(filepath) as f:
        if skip_header:
            next(f)
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) > max(key_col, val_col):
                mapping[parts[key_col]] = parts[val_col]
    return mapping


def load_fai_names(fai_path):
    """Load ordered sequence names from a .fai file."""
    if not fai_path or not os.path.exists(fai_path):
        return []
    with open(fai_path) as f:
        return [line.split()[0] for line in f if line.strip()]


def write_tsv(filepath, rows):
    """Write rows to TSV file with standard 6-column header."""
    headers = ["combined_SUPER", "mat_pat", "haplotype", "scaffold", "SUPER_internal", "SUPER_renamed"]
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    with open(filepath, "w") as f:
        f.write("\t".join(headers) + "\n")
        for r in rows:
            line = "\t".join(str(r[h]) for h in headers)
            f.write(line + "\n")
    print(f"SUCCESS: Wrote {len(rows)} records to: {filepath}")


def run_dual(hap1_file, hap2_file, hap2_vs_hap1_file, assign_file, fai_file, out_file):
    print("=== Generating Dual scaffold mapping ===")
    h1_map = read_tsv_dict(hap1_file)
    h2_map = read_tsv_dict(hap2_file)
    h2_to_h1 = read_tsv_dict(hap2_vs_hap1_file)

    # Load hap_assignment: key = e.g. 'SUPER_1.H1', val = ('mat', 'Mat_SUPER_1.H1')
    hap_assign = {}
    with open(assign_file) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.split()
                hap_assign[parts[2]] = (parts[0], parts[3])

    rows = []
    # 1. Hap1
    for scaf, s_int in h1_map.items():
        key = f"{s_int}.H1"
        if key in hap_assign:
            mp, final_name = hap_assign[key]
            rows.append({
                "combined_SUPER": final_name,
                "mat_pat": mp,
                "haplotype": "H1",
                "scaffold": scaf,
                "SUPER_internal": s_int,
                "SUPER_renamed": s_int
            })
        else:
            print(f"Warning: Dual H1 {scaf} ({key}) not found in hap_assignment", file=sys.stderr)

    # 2. Hap2
    for scaf, s_int in h2_map.items():
        s_ren = h2_to_h1.get(s_int, s_int)
        if s_int == "SUPER_39" or s_ren == "Scaffold_193":
            s_ren = "SUPER_39"
        key = f"{s_ren}.H2"
        if key in hap_assign:
            mp, final_name = hap_assign[key]
            rows.append({
                "combined_SUPER": final_name,
                "mat_pat": mp,
                "haplotype": "H2",
                "scaffold": scaf,
                "SUPER_internal": s_int,
                "SUPER_renamed": s_ren
            })
        else:
            print(f"Warning: Dual H2 {scaf} ({key}) not found in hap_assignment", file=sys.stderr)

    # Validate against final FASTA .fai if provided
    if fai_file:
        fai_names = load_fai_names(fai_file)
        row_finals = [r["combined_SUPER"] for r in rows]
        if fai_names:
            diff = set(fai_names) - set(row_finals)
            if diff:
                print(f"ERROR: Missing in Dual mapping: {diff}", file=sys.stderr)
            else:
                print(f"Validation: Dual mapped all {len(rows)}/{len(fai_names)} sequences matching .fai perfectly.")

    write_tsv(out_file, rows)


def run_single(hap1_file, hap2_file, hap2_vs_hap1_file, assign_file, fai_file, out_file):
    print("=== Generating Single scaffold mapping ===")
    h1_map = read_tsv_dict(hap1_file)
    h2_map = read_tsv_dict(hap2_file)
    h2_to_h1 = read_tsv_dict(hap2_vs_hap1_file)

    def remap_h2(s_int):
        if "_unloc" in s_int:
            base, unloc = s_int.split("_unloc", 1)
            return h2_to_h1.get(base, base) + "_unloc" + unloc
        return h2_to_h1.get(s_int, s_int)

    hap_assign_h1 = {}
    hap_assign_h2 = {}
    with open(assign_file) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.split()
                mp, super_name, final_name = parts[0], parts[2], parts[3]
                if final_name.endswith('.H1'):
                    hap_assign_h1[super_name] = (mp, final_name)
                elif final_name.endswith('.H2'):
                    hap_assign_h2[super_name] = (mp, final_name)

    rows = []
    # 1. Hap1 (Entries ending in .H1)
    for scaf, s_int in h1_map.items():
        if s_int in hap_assign_h1:
            mp, fn = hap_assign_h1[s_int]
        else:
            print(f"Warning: Single H1 {scaf} ({s_int}) not found in hap_assignment", file=sys.stderr)
            continue
        rows.append({
            "combined_SUPER": fn,
            "mat_pat": mp,
            "haplotype": "H1",
            "scaffold": scaf,
            "SUPER_internal": s_int,
            "SUPER_renamed": s_int
        })

    # 2. Hap2 (Entries ending in .H2)
    for scaf, s_int in h2_map.items():
        s_ren = remap_h2(s_int)
        if s_ren in hap_assign_h2:
            mp, fn = hap_assign_h2[s_ren]
        else:
            print(f"Warning: Single H2 {scaf} ({s_int} -> {s_ren}) not found in hap_assignment", file=sys.stderr)
            continue
        rows.append({
            "combined_SUPER": fn,
            "mat_pat": mp,
            "haplotype": "H2",
            "scaffold": scaf,
            "SUPER_internal": s_int,
            "SUPER_renamed": s_ren
        })

    # Validate against final FASTA .fai if provided
    if fai_file:
        fai_names = load_fai_names(fai_file)
        row_finals = [r["combined_SUPER"] for r in rows]
        if fai_names:
            diff = set(fai_names) - set(row_finals)
            if diff:
                print(f"ERROR: Missing in Single mapping: {diff}", file=sys.stderr)
            else:
                print(f"Validation: Single mapped all {len(rows)}/{len(fai_names)} sequences matching .fai perfectly.")

    write_tsv(out_file, rows)


def run_ont(hap1_file, hap2_file, rename_file, fai_file, out_file):
    print("=== Generating ONT scaffold mapping ===")
    h1_map = read_tsv_dict(hap1_file)
    h2_map = read_tsv_dict(hap2_file)
    rename_map = read_tsv_dict(rename_file, skip_header=True)

    rows = []
    # 1. Hap1: 39 Maternal autosomes + SUPER_W (Maternal) + SUPER_Z (Paternal)
    for scaf, s_int in h1_map.items():
        s_ren = rename_map.get(s_int, s_int)
        if s_int == "SUPER_Z":
            mp = "pat"
            fn = "Pat_SUPER_Z"
        elif s_int == "SUPER_W":
            mp = "mat"
            fn = "Mat_SUPER_W"
        else:
            mp = "mat"
            fn = f"Mat_{s_ren}"

        rows.append({
            "combined_SUPER": fn,
            "mat_pat": mp,
            "haplotype": "H1",
            "scaffold": scaf,
            "SUPER_internal": s_int,
            "SUPER_renamed": s_ren
        })

    # 2. Hap2: 39 Paternal autosomes
    for scaf, s_int in h2_map.items():
        s_ren = rename_map.get(s_int, s_int)
        mp = "pat"
        fn = f"Pat_{s_ren}"
        rows.append({
            "combined_SUPER": fn,
            "mat_pat": mp,
            "haplotype": "H2",
            "scaffold": scaf,
            "SUPER_internal": s_int,
            "SUPER_renamed": s_ren
        })

    # Validate against final FASTA .fai if provided
    if fai_file:
        fai_names = load_fai_names(fai_file)
        row_finals = [r["combined_SUPER"] for r in rows]
        if fai_names:
            diff = set(fai_names) - set(row_finals)
            if diff:
                print(f"ERROR: Missing in ONT mapping: {diff}", file=sys.stderr)
            else:
                print(f"Validation: ONT mapped all {len(rows)}/{len(fai_names)} sequences matching .fai perfectly.")

    write_tsv(out_file, rows)


def main():
    parser = argparse.ArgumentParser(description="Generate scaffold-to-super mapping table with explicit file inputs.")
    parser.add_argument("--mode", required=True, choices=["dual", "single", "ont", "all"],
                        help="Assembly mode to process")
    parser.add_argument("--hap1-inter", help="Path to hap1.inter_chr.tsv")
    parser.add_argument("--hap2-inter", help="Path to hap2.inter_chr.tsv")
    parser.add_argument("--hap2-vs-hap1", help="Path to hap2.vs.hap1.tsv")
    parser.add_argument("--hap-assignment", help="Path to hap_assignment.tsv (for dual/single)")
    parser.add_argument("--new-names", help="Path to new_name_change.tsv (for ont)")
    parser.add_argument("--fai", help="Path to final assembly .fai index for validation")
    parser.add_argument("--out", help="Path to output mapping TSV file")
    parser.add_argument("--base-dir", default="/lustre/fs5/vgl/scratch/eduarte/curations/birds/bTaeGut7/redo_curation_analysis",
                        help="Base directory when using --mode all")

    args = parser.parse_args()

    if args.mode == "dual":
        if not (args.hap1_inter and args.hap2_inter and args.hap2_vs_hap1 and args.hap_assignment and args.out):
            parser.error("--mode dual requires: --hap1-inter, --hap2-inter, --hap2-vs-hap1, --hap-assignment, --out")
        run_dual(args.hap1_inter, args.hap2_inter, args.hap2_vs_hap1, args.hap_assignment, args.fai, args.out)

    elif args.mode == "single":
        if not (args.hap1_inter and args.hap2_inter and args.hap2_vs_hap1 and args.hap_assignment and args.out):
            parser.error("--mode single requires: --hap1-inter, --hap2-inter, --hap2-vs-hap1, --hap-assignment, --out")
        run_single(args.hap1_inter, args.hap2_inter, args.hap2_vs_hap1, args.hap_assignment, args.fai, args.out)

    elif args.mode == "ont":
        if not (args.hap1_inter and args.hap2_inter and args.new_names and args.out):
            parser.error("--mode ont requires: --hap1-inter, --hap2-inter, --new-names, --out")
        run_ont(args.hap1_inter, args.hap2_inter, args.new_names, args.fai, args.out)

    elif args.mode == "all":
        b = args.base_dir
        run_dual(
            os.path.join(b, "data/dual/hap1.inter_chr.tsv"),
            os.path.join(b, "data/dual/hap2.inter_chr.tsv"),
            os.path.join(b, "data/dual/hap2.vs.hap1.tsv"),
            os.path.join(b, "results/dual/01_formatting_assemblies_dual/hap_assignment.tsv"),
            os.path.join(b, "results/dual/02_chainpipeline_dual/dual_combined.renamed.sorted.reoriented.fa.fai"),
            os.path.join(b, "data/dual/scaffold_hap_to_combined_super.tsv")
        )
        run_single(
            os.path.join(b, "data/single/hap1.inter_chr.tsv"),
            os.path.join(b, "data/single/hap2.inter_chr.tsv"),
            os.path.join(b, "data/single/hap2.vs.hap1.tsv"),
            os.path.join(b, "results/single/01_formatting_assemblies_single/hap_assignment.tsv"),
            os.path.join(b, "results/single/02_chainpipeline_single/single_combined.renamed.sorted.reoriented.fa.fai"),
            os.path.join(b, "data/single/scaffold_hap_to_combined_super.tsv")
        )
        run_ont(
            os.path.join(b, "data/ont/hap1.inter_chr.tsv"),
            os.path.join(b, "data/ont/hap2.inter_chr.tsv"),
            os.path.join(b, "data/ont/new_name_change.tsv"),
            os.path.join(b, "results/ont/02_chainpipeline_ont/asm3_ONT_combined.sorted.reoriented.fa.fai"),
            os.path.join(b, "data/ont/scaffold_hap_to_combined_super.tsv")
        )


if __name__ == "__main__":
    main()
