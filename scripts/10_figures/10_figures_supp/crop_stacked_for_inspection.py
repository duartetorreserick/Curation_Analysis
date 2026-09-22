#!/usr/bin/env python3
"""
crop_stacked_for_inspection.py
Creates crops of the generated stacked synteny ideogram for artifact inspection.
"""

import os
from PIL import Image

SRC = "/lustre/fs5/vgl/scratch/eduarte/curations/birds/bTaeGut7/redo_curation_analysis/results/figures/10_figures_supp/bTaeGut7_synteny_ideogram_stacked_20260922_133432.png"
ART_DIR = "/ru-auth/local/home/eduarte/.gemini/antigravity-cli/brain/98399b56-8147-461c-b26d-40b42a769012/figures"
os.makedirs(ART_DIR, exist_ok=True)

img = Image.open(SRC)
w, h = img.size
print(f"Image dimensions: {w} x {h}")

# Copy full image
dst_full = os.path.join(ART_DIR, "bTaeGut7_synteny_ideogram_stacked_20260922_133432.png")
img.save(dst_full)
print(f"Saved full image to {dst_full}")

# Crop 1: Top header & Chr 1 (top ~15% of image)
crop1 = img.crop((0, 0, w, int(h * 0.16)))
dst_crop1 = os.path.join(ART_DIR, "stacked_chr1_crop.png")
crop1.save(dst_crop1)
print(f"Saved Chr 1 crop to {dst_crop1}")

# Crop 2: Chr 3 & Chr 4 (around 22% to 40% of image)
crop2 = img.crop((0, int(h * 0.22), w, int(h * 0.40)))
dst_crop2 = os.path.join(ART_DIR, "stacked_chr3_4_crop.png")
crop2.save(dst_crop2)
print(f"Saved Chr 3-4 crop to {dst_crop2}")

# Crop 3: Bottom Chr W / Z & Legend (bottom ~16% of image)
crop3 = img.crop((0, int(h * 0.83), w, h))
dst_crop3 = os.path.join(ART_DIR, "stacked_chrz_legend_crop.png")
crop3.save(dst_crop3)
print(f"Saved Chr Z & legend crop to {dst_crop3}")
