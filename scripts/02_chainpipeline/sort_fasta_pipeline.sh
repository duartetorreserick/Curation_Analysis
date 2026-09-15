#!/bin/bash
# Pipeline to sort DUAL fasta to match T2T fasta order using the TSV mapping
# Usage: bash sort_fasta_pipeline.sh

set -euo pipefail

TSV=$1
T2T_FASTA=$2
DUAL_FASTA=$3
OUTPUT_FASTA=$4
REGIONS="ordered_regions.txt"
dual_fai=$5
# Step 1: Index the DUAL fasta
conda run -n vgp samtools faidx "$DUAL_FASTA"

# Step 2: Generate ordered regions list from T2T order + TSV mapping
python3 << PYEOF
import sys

tsv_file = "${TSV}"
t2t_fasta = "${T2T_FASTA}"
dual_fai = "${dual_fai}"

# Read TSV mapping: T2T_name -> SUPER_name
t2t_to_super = {}
with open(tsv_file) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) >= 2:
            t2t_to_super[parts[0]] = parts[1]

# Get T2T header order
t2t_order = []
with open(t2t_fasta) as f:
    for line in f:
        if line.startswith('>'):
            t2t_order.append(line[1:].split()[0])

# Get all DUAL sequence names from .fai
all_dual = []
with open(dual_fai) as f:
    for line in f:
        all_dual.append(line.split('\t')[0])

# Group unloc sequences by parent (including .mat/.pat suffix so the key
# matches the super_name values that come from the pairs TSV).
unloc_by_parent = {}
for name in all_dual:
    if '_unloc_' in name:
        base   = name.rsplit('_unloc_', 1)[0]          # e.g. 'SUPER_3'
        suffix = ('.' + name.rsplit('.', 1)[1]) if '.' in name else ''  # e.g. '.mat'
        parent = base + suffix                          # e.g. 'SUPER_3.mat'
        unloc_by_parent.setdefault(parent, []).append(name)
for parent in unloc_by_parent:
    # Strip suffix before int() so '2.mat' → 2
    unloc_by_parent[parent].sort(key=lambda x: int(x.rsplit('_', 1)[1].split('.')[0]))

# Build final order
written = set()
final_order = []
for t2t_name in t2t_order:
    if t2t_name in t2t_to_super:
        super_name = t2t_to_super[t2t_name]
        final_order.append(super_name)
        written.add(super_name)
        if super_name in unloc_by_parent:
            for unloc in unloc_by_parent[super_name]:
                final_order.append(unloc)
                written.add(unloc)
    else:
        print(f"WARNING: {t2t_name} not in TSV", file=sys.stderr)

# Add remaining sequences not in mapping
for name in all_dual:
    if name not in written:
        final_order.append(name)
        print(f"WARNING: {name} appended at end (not in mapping)", file=sys.stderr)

# Write ordered regions
with open("ordered_regions.txt", "w") as out:
    for name in final_order:
        out.write(name + "\n")

print(f"Total: {len(final_order)} sequences", file=sys.stderr)
PYEOF

# Step 3: Extract sequences in the desired order
conda run -n vgp samtools faidx "$DUAL_FASTA" $(cat "$REGIONS" | tr '\n' ' ') > "$OUTPUT_FASTA"

echo "Done. Sorted fasta written to $OUTPUT_FASTA"
