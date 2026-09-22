#!/usr/bin/env python3
"""
crop_macro_cytogenetic.py

Crops key benchmark regions from the generated macrochromosome synteny ideogram
for visual verification and embedding in artifacts.
"""

import os
from PIL import Image

SRC = "/lustre/fs5/vgl/scratch/eduarte/curations/birds/bTaeGut7/redo_curation_analysis/results/figures/10_figures_supp/bTaeGut7_synteny_ideogram_macro_20260922_140218.png"
ART_DIR = "/ru-auth/local/home/eduarte/.gemini/antigravity-cli/brain/98399b56-8147-461c-b26d-40b42a769012/figures"
os.makedirs(ART_DIR, exist_ok=True)

img = Image.open(SRC)
W, H = img.size
print(f"Source image dimensions: {W} x {H}")

# 1. Chr 1 crop (top of butterfly plot)
# X: 0 to W, Y: 0.03*H to 0.16*H
c1 = img.crop((0, int(0.02 * H), W, int(0.15 * H)))
c1_path = f"{ART_DIR}/crop_macro_chr1.png"
c1.save(c1_path)
print(f"Saved Chr 1 crop: {c1_path}")

# 2. Chr 3 & 4 crop (inversions, switch blocks, query insertions)
# Y: 0.22*H to 0.42*H
c34 = img.crop((0, int(0.20 * H), W, int(0.38 * H)))
c34_path = f"{ART_DIR}/crop_macro_chr3_4.png"
c34.save(c34_path)
print(f"Saved Chr 3-4 crop: {c34_path}")

# 3. Chr W | Z crop (bottom row: Chr Z unaligned region + unlinked scaffold u1)
# Y: 0.82*H to 0.96*H
cwz = img.crop((0, int(0.82 * H), W, int(0.95 * H)))
cwz_path = f"{ART_DIR}/crop_macro_chrzw.png"
cwz.save(cwz_path)
print(f"Saved Chr W/Z crop: {cwz_path}")

# 4. Legend crop
cleg = img.crop((int(0.05 * W), int(0.95 * H), int(0.95 * W), H))
cleg_path = f"{ART_DIR}/crop_macro_legend.png"
cleg.save(cleg_path)
print(f"Saved Legend crop: {cleg_path}")
