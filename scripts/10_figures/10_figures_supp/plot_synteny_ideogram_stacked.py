#!/usr/bin/env python3
"""
plot_synteny_ideogram_stacked.py

Generates publication-quality supplementary butterfly ideograms with stacked multi-alignment
synteny maps comparing a single shared T2T reference chromosome at the top against
the three assemblies (HiFi Single, HiFi Dual, and ONT Dual) stacked below it.

Layout per chromosome row:
  - Bar 0 (Top): T2T Reference chromosome in solid slate gray (#475569)
      * Centromere constriction marker (#C084FC)
      * Terminal telomeres (#10B981)
  - Ribbon Space 0: Synteny ribbons connecting T2T Reference to HiFi Single
  - Bar 1: HiFi Single assembly scaffold in solid royal blue (#2563EB)
      * Native curation gaps (protruding high-contrast tick marks)
      * Native assembly gaps (charcoal ticks)
      * Native hap-mer switch blocks (red solid blocks)
      * Native telomeres (green dots)
  - Ribbon Space 1: Synteny ribbons connecting HiFi Single to HiFi Dual
  - Bar 2: HiFi Dual assembly scaffold in solid emerald green (#059669)
      * Native curation gaps, assembly gaps, hap-mers, telomeres
  - Ribbon Space 2: Synteny ribbons connecting HiFi Dual to ONT Dual
  - Bar 3: ONT Dual assembly scaffold in solid warm orange (#EA580C)
      * Native curation gaps, assembly gaps, hap-mers, telomeres
"""

import argparse
import os
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
from matplotlib.gridspec import GridSpec
import numpy as np


# =============================================================================
# Constants & Colors
# =============================================================================

MACRO_TOKENS = ['1', '1A', '2', '3', '4', '4A', '5', '6', '7', '8', 'ZW']

# Bar Colors (Solid saturated fills)
COL_T2T_BG     = '#475569'   # Solid slate gray
COL_T2T_BORDER = '#1E293B'

COL_S_BG       = '#2563EB'   # Solid royal blue
COL_S_BORDER   = '#1D4ED8'

COL_D_BG       = '#059669'   # Solid emerald green
COL_D_BORDER   = '#047857'

COL_O_BG       = '#EA580C'   # Solid warm orange
COL_O_BORDER   = '#C2410C'

# Ribbon Colors
# Tier 0 (T2T <-> Single): Sky blue
COL_RIB0_COLL   = '#38BDF8'
COL_RIB0_COLL_E = '#0284C7'
COL_RIB0_NC     = '#F59E0B'
COL_RIB0_NC_E   = '#B45309'

# Tier 1 (Single <-> Dual): Cyan/Teal
COL_RIB1_COLL   = '#2DD4BF'
COL_RIB1_COLL_E = '#0F766E'
COL_RIB1_NC     = '#F59E0B'
COL_RIB1_NC_E   = '#B45309'

# Tier 2 (Dual <-> ONT): Amber/Orange
COL_RIB2_COLL   = '#FBBF24'
COL_RIB2_COLL_E = '#D97706'
COL_RIB2_NC     = '#F43F5E'
COL_RIB2_NC_E   = '#BE123C'

# Feature Colors
COL_CUR_GAP    = '#FEF08A'   # High-contrast bright yellow/white on colored bars
COL_CUR_GAP_ED = '#DC2626'   # Red shadow/edge
COL_ASM_GAP    = '#0F172A'   # Deep charcoal
COL_SWITCH     = '#EF4444'   # Vibrant red block
COL_SWITCH_ED  = '#7F1D1D'
COL_TELO       = '#10B981'   # Bright emerald green
COL_CENTRO     = '#C084FC'   # Bright purple constriction


# =============================================================================
# Parsers & Data Loaders
# =============================================================================

def load_seq_sizes(tsv_path):
    """Loads target or query sequence sizes from TSV."""
    sizes = {}
    if not tsv_path or not os.path.exists(tsv_path):
        return sizes
    with open(tsv_path) as fh:
        for line in fh:
            if line.startswith('#'): continue
            p = line.strip().split()
            if len(p) >= 2 and p[0] != 'name':
                sizes[p[0]] = int(p[1])
    return sizes


def load_chrom_pairs(path):
    """Loads t2t_chrom -> asm_scaffold mapping."""
    pairs = {}
    if not path or not os.path.exists(path):
        return pairs
    with open(path) as fh:
        for line in fh:
            p = line.strip().split()
            if len(p) >= 2:
                pairs[p[0]] = p[1]
    return pairs


def load_centromeres(gff_path):
    """Loads centromere spans from centromere detector GFF."""
    centromeres = {}
    if not gff_path or not os.path.exists(gff_path):
        return centromeres
    with open(gff_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split('\t')
            if len(parts) < 5: continue
            m = re.match(r'^chr([0-9]+[A-Za-z]*|[A-Za-z]+)_(pat|mat)$', parts[0])
            if m:
                tok, side = m.group(1), m.group(2)
                cs, ce = int(parts[3]), int(parts[4])
                key = (tok.upper(), side.lower())
                if key not in centromeres:
                    centromeres[key] = (cs, ce)
                else:
                    prev_s, prev_e = centromeres[key]
                    centromeres[key] = (min(prev_s, cs), max(prev_e, ce))
    return centromeres


def load_p_arm_bed(bed_path):
    """Loads p-arm telomere positions to determine chromosome orientation."""
    p_starts = {}
    if not bed_path or not os.path.exists(bed_path):
        return p_starts
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            p = line.split('\t')
            if len(p) >= 5:
                if p[4].lower() == 'p' or p[3].lower() == 'p':
                    p_starts[p[0]] = int(p[1])
            elif len(p) >= 2:
                p_starts[p[0]] = int(p[1])
    return p_starts


def build_flip_set_from_centromeres_and_telomeres(seq_sizes, centromeres, p_starts):
    """
    Determines which chromosomes should be flipped so p-arm or centromere
    is consistently oriented toward the center of the butterfly plot.
    """
    flip_set = set()
    for chrom, sz in seq_sizes.items():
        if chrom in p_starts:
            if p_starts[chrom] > sz / 2:
                flip_set.add(chrom)
                continue
        m = re.search(r'chromosome_(\w+)', chrom)
        tok = m.group(1).upper() if m else ''
        side = 'pat' if ('Pat' in chrom or tok == 'Z') else 'mat'
        key = (tok, side)
        if key in centromeres:
            cs, ce = centromeres[key]
            if (cs + ce) / 2 > sz / 2:
                flip_set.add(chrom)
    return flip_set


def load_telomere_presence_tsv(path):
    """Loads telomere presence calls for an assembly."""
    telo = defaultdict(set)
    if not path or not os.path.exists(path):
        return dict(telo)
    with open(path) as f:
        header = None
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            p = line.split('\t')
            if header is None:
                header = [c.lower() for c in p]
                continue
            chrom = p[0]
            coll_val = p[2].lower() if len(p) > 2 else 'none'
            nc_val   = p[3].lower() if len(p) > 3 else 'none'
            for arm in ['p', 'q']:
                if arm in coll_val or arm in nc_val:
                    telo[chrom].add(arm)
    return dict(telo)


def load_t2t_terminal_telomeres(bed_path):
    """Loads T2T terminal telomeres."""
    telo = defaultdict(set)
    if not bed_path or not os.path.exists(bed_path):
        return dict(telo)
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            p = line.split('\t')
            if len(p) >= 5:
                chrom = p[0]
                arm = p[4].lower()
                if arm in ('p', 'q'):
                    telo[chrom].add(arm)
    return dict(telo)


def load_gaps_bed(bed_path):
    """
    Loads native assembly gaps from annotated gaps BED.
    Returns: {scaffold: [(start, end, 'CURATION' | 'ASSEMBLY')]}
    """
    gaps = defaultdict(list)
    if not bed_path or not os.path.exists(bed_path):
        return dict(gaps)
    with open(bed_path) as fh:
        for line in fh:
            p = line.strip().split()
            if len(p) >= 4:
                gaps[p[0]].append((int(p[1]), int(p[2]), p[3].upper()))
    return dict(gaps)


def load_switch_blocks(bed_path):
    """
    Loads native hap-mer switch error blocks from BED.
    Returns: {scaffold: [(start, end)]}
    """
    sw = defaultdict(list)
    if not bed_path or not os.path.exists(bed_path):
        return dict(sw)
    with open(bed_path) as fh:
        for line in fh:
            p = line.strip().split()
            if len(p) >= 3:
                sw[p[0]].append((int(p[1]), int(p[2])))
    return dict(sw)


def extract_synteny_ribbons(chain_path, max_gap=25000, min_block_size=5000, is_collinear=True):
    """
    Extracts contiguous syntenic ribbon segments connecting (tStart, tEnd) to (qStart, qEnd).
    Merges internal alignment blocks separated by indels <= max_gap.
    """
    ribbons = []
    if not chain_path or not os.path.exists(chain_path):
        return ribbons

    def _proc(hdr, blks):
        t_pos = hdr['tStart']
        q_pos = hdr['qStart']
        q_size = hdr['qSize']
        qStrand = hdr['qStrand']

        cur_ts = None; cur_te = None
        cur_qs = None; cur_qe = None

        for blk in blks:
            sz = blk[0]
            dt = blk[1] if len(blk) > 1 else 0
            dq = blk[2] if len(blk) > 2 else 0

            blk_ts = t_pos
            blk_te = t_pos + sz
            if qStrand == '-':
                blk_qs = q_size - (q_pos + sz)
                blk_qe = q_size - q_pos
            else:
                blk_qs = q_pos
                blk_qe = q_pos + sz

            if cur_ts is None:
                cur_ts, cur_te = blk_ts, blk_te
                cur_qs, cur_qe = blk_qs, blk_qe
            else:
                gap_t = blk_ts - cur_te
                gap_q = blk_qs - cur_qe if qStrand == '+' else cur_qs - blk_qe

                if gap_t <= max_gap and gap_q <= max_gap and gap_t >= 0 and gap_q >= 0:
                    cur_te = blk_te
                    if qStrand == '+':
                        cur_qe = blk_qe
                    else:
                        cur_qs = blk_qs
                else:
                    if (cur_te - cur_ts) >= min_block_size:
                        ribbons.append((
                            hdr['tName'], cur_ts, cur_te,
                            hdr['qName'],
                            min(cur_qs, cur_qe), max(cur_qs, cur_qe),
                            qStrand, is_collinear
                        ))
                    cur_ts, cur_te = blk_ts, blk_te
                    cur_qs, cur_qe = blk_qs, blk_qe

            t_pos += sz + dt
            q_pos += sz + dq

        if cur_ts is not None and (cur_te - cur_ts) >= min_block_size:
            ribbons.append((
                hdr['tName'], cur_ts, cur_te,
                hdr['qName'],
                min(cur_qs, cur_qe), max(cur_qs, cur_qe),
                qStrand, is_collinear
            ))

    with open(chain_path) as fh:
        hdr = None
        blks = []
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'): continue
            if line.startswith('chain'):
                if hdr and blks:
                    _proc(hdr, blks)
                p = line.split()
                hdr = {
                    'tName': p[2], 'tSize': int(p[3]), 'tStart': int(p[5]), 'tEnd': int(p[6]),
                    'qName': p[7], 'qSize': int(p[8]), 'qStrand': p[9],
                    'qStart': int(p[10]), 'qEnd': int(p[11])
                }
                blks = []
            else:
                blks.append([int(x) for x in line.split()])
        if hdr and blks:
            _proc(hdr, blks)

    return ribbons


def derive_pairwise_synteny(ribbons_A, ribbons_B, min_overlap=5000, max_gap=25000):
    """
    Derives direct pairwise synteny ribbons between Assembly A and Assembly B
    through their common coordinate space on T2T.
    Returns: list of (qA_start, qA_end, qB_start, qB_end, rel_strand, is_collinear)
    """
    raw_blocks = []
    for rA in ribbons_A:
        tsA, teA = rA[1], rA[2]
        qsA, qeA = rA[4], rA[5]
        strA = rA[6]
        colA = rA[7]
        lenA = max(1, teA - tsA)

        for rB in ribbons_B:
            tsB, teB = rB[1], rB[2]
            qsB, qeB = rB[4], rB[5]
            strB = rB[6]
            colB = rB[7]
            lenB = max(1, teB - tsB)

            # Overlap on T2T
            ov_s = max(tsA, tsB)
            ov_e = min(teA, teB)
            ov_len = ov_e - ov_s

            if ov_len >= min_overlap:
                f0_A = (ov_s - tsA) / lenA
                f1_A = (ov_e - tsA) / lenA
                if strA == '+':
                    qA0 = qsA + f0_A * (qeA - qsA)
                    qA1 = qsA + f1_A * (qeA - qsA)
                else:
                    qA0 = qeA - f1_A * (qeA - qsA)
                    qA1 = qeA - f0_A * (qeA - qsA)

                f0_B = (ov_s - tsB) / lenB
                f1_B = (ov_e - tsB) / lenB
                if strB == '+':
                    qB0 = qsB + f0_B * (qeB - qsB)
                    qB1 = qsB + f1_B * (qeB - qsB)
                else:
                    qB0 = qeB - f1_B * (qeB - qsB)
                    qB1 = qeB - f0_B * (qeB - qsB)

                rel_strand = '+' if strA == strB else '-'
                is_coll = (colA and colB and rel_strand == '+')

                raw_blocks.append((
                    int(min(qA0, qA1)), int(max(qA0, qA1)),
                    int(min(qB0, qB1)), int(max(qB0, qB1)),
                    rel_strand, is_coll
                ))

    if not raw_blocks:
        return []

    # Sort and merge contiguous blocks
    raw_blocks.sort(key=lambda x: (x[0], x[2]))
    merged = []
    cur = raw_blocks[0]

    for nxt in raw_blocks[1:]:
        gap_A = nxt[0] - cur[1]
        gap_B = nxt[2] - cur[3] if cur[4] == '+' else cur[2] - nxt[3]

        if (nxt[4] == cur[4] and nxt[5] == cur[5] and
            0 <= gap_A <= max_gap and 0 <= gap_B <= max_gap):
            cur = (
                cur[0], max(cur[1], nxt[1]),
                min(cur[2], nxt[2]), max(cur[3], nxt[3]),
                cur[4], cur[5]
            )
        else:
            if (cur[1] - cur[0]) >= min_overlap:
                merged.append(cur)
            cur = nxt

    if (cur[1] - cur[0]) >= min_overlap:
        merged.append(cur)

    return merged


# =============================================================================
# Drawing Helpers
# =============================================================================

def draw_synteny_ribbon(ax, x0_top, x1_top, y_top, x0_bot, x1_bot, y_bot,
                       color, edge_color, alpha=0.40, inverted=False):
    """Draws a smooth cubic Bezier synteny ribbon between adjacent bars."""
    if x1_top <= x0_top or x1_bot <= x0_bot:
        return
    dy = y_top - y_bot
    cy_top = y_top - dy * 0.45
    cy_bot = y_bot + dy * 0.45

    if not inverted:
        pts = [
            (x0_top, y_top),
            (x1_top, y_top),
            (x1_top, cy_top), (x1_bot, cy_bot), (x1_bot, y_bot),
            (x0_bot, y_bot),
            (x0_bot, cy_bot), (x0_top, cy_top), (x0_top, y_top)
        ]
    else:
        # Cross / twist for inversion
        pts = [
            (x0_top, y_top),
            (x1_top, y_top),
            (x1_top, cy_top), (x0_bot, cy_bot), (x0_bot, y_bot),
            (x1_bot, y_bot),
            (x1_bot, cy_bot), (x0_top, cy_top), (x0_top, y_top)
        ]
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4
    ]
    path = Path(pts, codes)
    patch = mpatches.PathPatch(path, facecolor=color, edgecolor=edge_color,
                               alpha=alpha, lw=0.45, zorder=4)
    ax.add_patch(patch)


def draw_ideogram_bar(ax, x0, x1, y_center, h_bar, fc, ec, lw=0.9, zorder=6):
    """Draws a solid filled chromosome ideogram bar with a clean outline."""
    rect = mpatches.Rectangle((x0, y_center - h_bar / 2), x1 - x0, h_bar,
                              facecolor=fc, edgecolor=ec, lw=lw, zorder=zorder)
    ax.add_patch(rect)


def draw_centromere_marker(ax, cs, ce, y_center, h_bar, color=COL_CENTRO):
    """Draws a centromere constriction marker on the bar."""
    cw = max(ce - cs, 120_000)
    cx = (cs + ce) / 2
    x0 = max(0, cx - cw / 2)
    x1 = cx + cw / 2
    rect = mpatches.Rectangle((x0, y_center - h_bar / 2), x1 - x0, h_bar,
                              facecolor=color, edgecolor='#581C87', lw=0.7, zorder=7)
    ax.add_patch(rect)


def draw_telomere_marker(ax, x_pos, y_center, h_bar, color=COL_TELO):
    """Draws a small telomere cap dot at the chromosome end."""
    ax.plot([x_pos], [y_center], marker='o', markersize=3.8,
            markerfacecolor=color, markeredgecolor='#064E3B', markeredgewidth=0.6,
            zorder=10)


def draw_curation_gap_tick(ax, x_pos, y_center, h_bar):
    """Draws a high-contrast vertical tick mark for curation joins."""
    y0 = y_center - h_bar * 0.65
    y1 = y_center + h_bar * 0.65
    # Red background line slightly wider
    ax.plot([x_pos, x_pos], [y0, y1], color=COL_CUR_GAP_ED, lw=1.8, zorder=8)
    # Bright yellow center line
    ax.plot([x_pos, x_pos], [y0, y1], color=COL_CUR_GAP, lw=1.0, zorder=9)


def draw_assembly_gap_tick(ax, x_pos, y_center, h_bar):
    """Draws a charcoal tick mark for native contig gaps."""
    y0 = y_center - h_bar * 0.50
    y1 = y_center + h_bar * 0.50
    ax.plot([x_pos, x_pos], [y0, y1], color=COL_ASM_GAP, lw=0.9, zorder=8)


def draw_switch_block(ax, x0, x1, y_center, h_bar):
    """Draws a red switch block on the assembly bar."""
    w = max(x1 - x0, 80_000)
    rect = mpatches.Rectangle((x0, y_center - h_bar / 2), w, h_bar,
                              facecolor=COL_SWITCH, edgecolor=COL_SWITCH_ED,
                              lw=0.7, zorder=7)
    ax.add_patch(rect)


# =============================================================================
# Main Butterfly Plotting Engine
# =============================================================================

def plot_butterfly_stacked(
    t2t_sizes, flip_set,
    s_query_sizes, s_pairs, s_col_ribbons, s_nc_ribbons, s_gaps, s_sw, s_telo,
    d_query_sizes, d_pairs, d_col_ribbons, d_nc_ribbons, d_gaps, d_sw, d_telo,
    o_query_sizes, o_pairs, o_col_ribbons, o_nc_ribbons, o_gaps, o_sw, o_telo,
    centromeres, t2t_telo, out_png, out_pdf, dpi=300
):
    """Builds the 4-tier stacked multi-alignment ideogram figure."""
    tokens = MACRO_TOKENS
    n_tokens = len(tokens)

    # 1. Global X scale
    max_global_size = max(t2t_sizes.values()) if t2t_sizes else 160_000_000
    for qs in (s_query_sizes, d_query_sizes, o_query_sizes):
        if qs:
            max_global_size = max(max_global_size, max(qs.values()))
    x_margin = max_global_size * 0.02

    # 2. Geometry layout
    # Per chromosome group:
    # Bar 0: T2T Reference
    # Ribbon 0: T2T <-> Single
    # Bar 1: HiFi Single
    # Ribbon 1: Single <-> Dual
    # Bar 2: HiFi Dual
    # Ribbon 2: Dual <-> ONT
    # Bar 3: ONT Dual
    H_BAR       = 0.22                              # Bar thickness
    H_RIBBON    = 0.42                              # Synteny space between adjacent bars
    ASM_STEP    = H_BAR + H_RIBBON                  # 0.64
    GROUP_H     = 3 * ASM_STEP + H_BAR              # 4 bars + 3 ribbon tracks = 2.14
    CHROM_GAP   = 0.65                              # Gap between chromosome groups
    ROW_H       = GROUP_H + CHROM_GAP               # 2.79

    total_plot_height = n_tokens * ROW_H + 1.2

    # 3. Figure & GridSpec
    fig_w, fig_h = 20.0, 26.0
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')

    # GridSpec: [ Maternal (Left) | Labels (Center) | Paternal (Right) ]
    gs = GridSpec(1, 3, figure=fig, width_ratios=[10.0, 1.85, 10.0],
                  wspace=0.03, left=0.035, right=0.965, top=0.95, bottom=0.05)

    ax_mat = fig.add_subplot(gs[0, 0])
    ax_lbl = fig.add_subplot(gs[0, 1])
    ax_pat = fig.add_subplot(gs[0, 2])

    y_top = total_plot_height
    y_bot = 0.0

    for ax in (ax_mat, ax_lbl, ax_pat):
        ax.set_ylim(y_bot, y_top)
        ax.set_yticks([])
        ax.set_xticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    # Maternal: 0 near center (right), max_size at left
    ax_mat.set_xlim(max_global_size + x_margin, -x_margin)
    # Paternal: 0 near center (left), max_size at right
    ax_pat.set_xlim(-x_margin, max_global_size + x_margin)
    # Labels
    ax_lbl.set_xlim(0, 1)

    # Pre-index ribbons by (tName, qName)
    def index_ribbons(rib_list):
        idx = defaultdict(list)
        for r in rib_list:
            idx[(r[0], r[3])].append(r)
        return idx

    s_col_idx = index_ribbons(s_col_ribbons)
    s_nc_idx  = index_ribbons(s_nc_ribbons)
    d_col_idx = index_ribbons(d_col_ribbons)
    d_nc_idx  = index_ribbons(d_nc_ribbons)
    o_col_idx = index_ribbons(o_col_ribbons)
    o_nc_idx  = index_ribbons(o_nc_ribbons)

    # 4. Render Chromosome Groups
    for row_idx, tok in enumerate(tokens):
        y_group_top = y_top - 0.6 - row_idx * ROW_H
        y_t2t = y_group_top - H_BAR / 2
        y_s   = y_t2t - ASM_STEP
        y_d   = y_s - ASM_STEP
        y_o   = y_d - ASM_STEP

        # Divider line between chromosome groups
        if row_idx > 0:
            y_div = y_group_top + CHROM_GAP / 2
            ax_mat.axhline(y_div, color='#E2E8F0', lw=0.8, ls=':')
            ax_pat.axhline(y_div, color='#E2E8F0', lw=0.8, ls=':')
            ax_lbl.axhline(y_div, color='#CBD5E1', lw=0.8, ls='--')

        # Chromosome Badge Header in Center Column
        lbl_text = 'Chr W | Z' if tok == 'ZW' else f'Chr {tok}'
        ax_lbl.text(0.5, y_group_top + 0.14, lbl_text, ha='center', va='bottom',
                    fontsize=12, fontweight='bold', color='#0F172A',
                    bbox=dict(boxstyle='round,pad=0.28', facecolor='#F8FAFC', edgecolor='#CBD5E1', lw=0.8))

        # Center Column labels for each tier
        ax_lbl.text(0.5, y_t2t, 'T2T Ref', ha='center', va='center',
                    fontsize=7.2, color='#475569', fontweight='bold')
        ax_lbl.text(0.5, (y_t2t + y_s) / 2, '↕', ha='center', va='center',
                    fontsize=8.0, color='#94A3B8')

        ax_lbl.text(0.5, y_s, 'HiFi Single', ha='center', va='center',
                    fontsize=7.2, color='#1D4ED8', fontweight='bold')
        ax_lbl.text(0.5, (y_s + y_d) / 2, '↕', ha='center', va='center',
                    fontsize=8.0, color='#94A3B8')

        ax_lbl.text(0.5, y_d, 'HiFi Dual', ha='center', va='center',
                    fontsize=7.2, color='#047857', fontweight='bold')
        ax_lbl.text(0.5, (y_d + y_o) / 2, '↕', ha='center', va='center',
                    fontsize=8.0, color='#94A3B8')

        ax_lbl.text(0.5, y_o, 'ONT Dual', ha='center', va='center',
                    fontsize=7.2, color='#C2410C', fontweight='bold')

        # Draw Maternal and Paternal
        for side, ax in [('mat', ax_mat), ('pat', ax_pat)]:
            # Resolve T2T chromosome name
            if tok == 'ZW':
                t2t_chrom = 'Mat_NC_133064.1_chromosome_W' if side == 'mat' else 'Mat_NC_133063.1_chromosome_Z'
            else:
                t2t_chrom = next(
                    (k for k in t2t_sizes.keys()
                     if f'chromosome_{tok}' in k and (side == 'pat' if ('Pat' in k or tok == 'Z') else side == 'mat')),
                    None
                )

            if not t2t_chrom:
                continue

            t2t_len = t2t_sizes.get(t2t_chrom, 0)
            if t2t_len == 0:
                continue

            flip = t2t_chrom in flip_set
            p_arm = 'p' if not flip else 'q'
            q_arm = 'q' if not flip else 'p'

            # Resolve query scaffolds & lengths
            q_s = s_pairs.get(t2t_chrom)
            len_s = s_query_sizes.get(q_s, t2t_len) if q_s else t2t_len

            q_d = d_pairs.get(t2t_chrom)
            len_d = d_query_sizes.get(q_d, t2t_len) if q_d else t2t_len

            q_o = o_pairs.get(t2t_chrom)
            len_o = o_query_sizes.get(q_o, t2t_len) if q_o else t2t_len

            # ---------------------------------------------------------
            # 4a. Bar 0: T2T Reference
            # ---------------------------------------------------------
            draw_ideogram_bar(ax, 0, t2t_len, y_t2t, H_BAR,
                              COL_T2T_BG, COL_T2T_BORDER, lw=0.8, zorder=6)

            # Centromere
            tok_cen = 'W' if (tok == 'ZW' and side == 'mat') else ('Z' if (tok == 'ZW' and side == 'pat') else tok)
            cen_span = centromeres.get((tok_cen.upper(), side.lower()))
            if cen_span:
                cs = (t2t_len - cen_span[1]) if flip else cen_span[0]
                ce = (t2t_len - cen_span[0]) if flip else cen_span[1]
                draw_centromere_marker(ax, min(cs, ce), max(cs, ce), y_t2t, H_BAR, COL_CENTRO)

            # Telomeres T2T
            t2t_arms = t2t_telo.get(t2t_chrom, set())
            if p_arm in t2t_arms: draw_telomere_marker(ax, 0, y_t2t, H_BAR, COL_TELO)
            if q_arm in t2t_arms: draw_telomere_marker(ax, t2t_len, y_t2t, H_BAR, COL_TELO)
            if not t2t_arms:
                draw_telomere_marker(ax, 0, y_t2t, H_BAR, COL_TELO)
                draw_telomere_marker(ax, t2t_len, y_t2t, H_BAR, COL_TELO)

            # ---------------------------------------------------------
            # 4b. Ribbon Space 0: T2T <-> HiFi Single
            # ---------------------------------------------------------
            rib0_col = s_col_idx.get((t2t_chrom, q_s), [])
            rib0_nc  = s_nc_idx.get((t2t_chrom, q_s), [])

            for r in rib0_col:
                _, ts, te, _, qs, qe, qstr, _ = r
                if flip:
                    t0, t1 = t2t_len - te, t2t_len - ts
                    q0, q1 = len_s - qe, len_s - qs
                else:
                    t0, t1 = ts, te
                    q0, q1 = qs, qe
                draw_synteny_ribbon(ax, min(t0, t1), max(t0, t1), y_t2t - H_BAR/2,
                                    min(q0, q1), max(q0, q1), y_s + H_BAR/2,
                                    COL_RIB0_COLL, COL_RIB0_COLL_E, alpha=0.40,
                                    inverted=(qstr == '-'))

            for r in rib0_nc:
                _, ts, te, _, qs, qe, qstr, _ = r
                if flip:
                    t0, t1 = t2t_len - te, t2t_len - ts
                    q0, q1 = len_s - qe, len_s - qs
                else:
                    t0, t1 = ts, te
                    q0, q1 = qs, qe
                draw_synteny_ribbon(ax, min(t0, t1), max(t0, t1), y_t2t - H_BAR/2,
                                    min(q0, q1), max(q0, q1), y_s + H_BAR/2,
                                    COL_RIB0_NC, COL_RIB0_NC_E, alpha=0.60,
                                    inverted=(qstr == '-'))

            # ---------------------------------------------------------
            # 4c. Bar 1: HiFi Single
            # ---------------------------------------------------------
            draw_ideogram_bar(ax, 0, len_s, y_s, H_BAR,
                              COL_S_BG, COL_S_BORDER, lw=0.8, zorder=6)

            if q_s:
                # Gaps
                for gs, ge, gtype in s_gaps.get(q_s, []):
                    gx = ((len_s - (gs + ge) / 2) if flip else ((gs + ge) / 2))
                    if 'CURATION' in gtype:
                        draw_curation_gap_tick(ax, gx, y_s, H_BAR)
                    else:
                        draw_assembly_gap_tick(ax, gx, y_s, H_BAR)
                # Switch blocks
                for sw_s, sw_e in s_sw.get(q_s, []):
                    p0 = (len_s - sw_e) if flip else sw_s
                    p1 = (len_s - sw_s) if flip else sw_e
                    draw_switch_block(ax, min(p0, p1), max(p0, p1), y_s, H_BAR)
                # Telomeres
                asm_arms = s_telo.get(t2t_chrom, set())
                if p_arm in asm_arms: draw_telomere_marker(ax, 0, y_s, H_BAR, COL_TELO)
                if q_arm in asm_arms: draw_telomere_marker(ax, len_s, y_s, H_BAR, COL_TELO)

            # ---------------------------------------------------------
            # 4d. Ribbon Space 1: HiFi Single <-> HiFi Dual
            # ---------------------------------------------------------
            all_s_rib = rib0_col + rib0_nc
            rib_d_col = d_col_idx.get((t2t_chrom, q_d), [])
            rib_d_nc  = d_nc_idx.get((t2t_chrom, q_d), [])
            all_d_rib = rib_d_col + rib_d_nc

            pairwise_sd = derive_pairwise_synteny(all_s_rib, all_d_rib)
            for r in pairwise_sd:
                qa0, qa1, qb0, qb1, rel_str, is_coll = r
                if flip:
                    s0, s1 = len_s - qa1, len_s - qa0
                    d0, d1 = len_d - qb1, len_d - qb0
                else:
                    s0, s1 = qa0, qa1
                    d0, d1 = qb0, qb1
                c_fill = COL_RIB1_COLL if is_coll else COL_RIB1_NC
                c_edge = COL_RIB1_COLL_E if is_coll else COL_RIB1_NC_E
                draw_synteny_ribbon(ax, min(s0, s1), max(s0, s1), y_s - H_BAR/2,
                                    min(d0, d1), max(d0, d1), y_d + H_BAR/2,
                                    c_fill, c_edge, alpha=0.40 if is_coll else 0.60,
                                    inverted=(rel_str == '-'))

            # ---------------------------------------------------------
            # 4e. Bar 2: HiFi Dual
            # ---------------------------------------------------------
            draw_ideogram_bar(ax, 0, len_d, y_d, H_BAR,
                              COL_D_BG, COL_D_BORDER, lw=0.8, zorder=6)

            if q_d:
                # Gaps
                for gs, ge, gtype in d_gaps.get(q_d, []):
                    gx = ((len_d - (gs + ge) / 2) if flip else ((gs + ge) / 2))
                    if 'CURATION' in gtype:
                        draw_curation_gap_tick(ax, gx, y_d, H_BAR)
                    else:
                        draw_assembly_gap_tick(ax, gx, y_d, H_BAR)
                # Switch blocks
                for sw_s, sw_e in d_sw.get(q_d, []):
                    p0 = (len_d - sw_e) if flip else sw_s
                    p1 = (len_d - sw_s) if flip else sw_e
                    draw_switch_block(ax, min(p0, p1), max(p0, p1), y_d, H_BAR)
                # Telomeres
                asm_arms = d_telo.get(t2t_chrom, set())
                if p_arm in asm_arms: draw_telomere_marker(ax, 0, y_d, H_BAR, COL_TELO)
                if q_arm in asm_arms: draw_telomere_marker(ax, len_d, y_d, H_BAR, COL_TELO)

            # ---------------------------------------------------------
            # 4f. Ribbon Space 2: HiFi Dual <-> ONT Dual
            # ---------------------------------------------------------
            rib_o_col = o_col_idx.get((t2t_chrom, q_o), [])
            rib_o_nc  = o_nc_idx.get((t2t_chrom, q_o), [])
            all_o_rib = rib_o_col + rib_o_nc

            pairwise_do = derive_pairwise_synteny(all_d_rib, all_o_rib)
            for r in pairwise_do:
                qa0, qa1, qb0, qb1, rel_str, is_coll = r
                if flip:
                    d0, d1 = len_d - qa1, len_d - qa0
                    o0, o1 = len_o - qb1, len_o - qb0
                else:
                    d0, d1 = qa0, qa1
                    o0, o1 = qb0, qb1
                c_fill = COL_RIB2_COLL if is_coll else COL_RIB2_NC
                c_edge = COL_RIB2_COLL_E if is_coll else COL_RIB2_NC_E
                draw_synteny_ribbon(ax, min(d0, d1), max(d0, d1), y_d - H_BAR/2,
                                    min(o0, o1), max(o0, o1), y_o + H_BAR/2,
                                    c_fill, c_edge, alpha=0.40 if is_coll else 0.60,
                                    inverted=(rel_str == '-'))

            # ---------------------------------------------------------
            # 4g. Bar 3: ONT Dual
            # ---------------------------------------------------------
            draw_ideogram_bar(ax, 0, len_o, y_o, H_BAR,
                              COL_O_BG, COL_O_BORDER, lw=0.8, zorder=6)

            if q_o:
                # Gaps
                for gs, ge, gtype in o_gaps.get(q_o, []):
                    gx = ((len_o - (gs + ge) / 2) if flip else ((gs + ge) / 2))
                    if 'CURATION' in gtype:
                        draw_curation_gap_tick(ax, gx, y_o, H_BAR)
                    else:
                        draw_assembly_gap_tick(ax, gx, y_o, H_BAR)
                # Switch blocks
                for sw_s, sw_e in o_sw.get(q_o, []):
                    p0 = (len_o - sw_e) if flip else sw_s
                    p1 = (len_o - sw_s) if flip else sw_e
                    draw_switch_block(ax, min(p0, p1), max(p0, p1), y_o, H_BAR)
                # Telomeres
                asm_arms = o_telo.get(t2t_chrom, set())
                if p_arm in asm_arms: draw_telomere_marker(ax, 0, y_o, H_BAR, COL_TELO)
                if q_arm in asm_arms: draw_telomere_marker(ax, len_o, y_o, H_BAR, COL_TELO)

    # 5. Header Titles
    ax_mat.set_title('Maternal Haplotype', fontsize=16, fontweight='bold',
                     pad=18, color='#0F172A')
    ax_lbl.set_title('Chromosome', fontsize=12, fontweight='bold',
                     pad=18, color='#475569')
    ax_pat.set_title('Paternal Haplotype', fontsize=16, fontweight='bold',
                     pad=18, color='#0F172A')

    # 6. Scale Bars at bottom
    scale_len = 20_000_000 # 20 Mb
    scale_y = 0.50

    # Maternal scale bar
    ax_mat.plot([0, scale_len], [scale_y, scale_y], color='#0F172A', lw=2.2)
    ax_mat.plot([0, 0], [scale_y - 0.12, scale_y + 0.12], color='#0F172A', lw=2.0)
    ax_mat.plot([scale_len, scale_len], [scale_y - 0.12, scale_y + 0.12], color='#0F172A', lw=2.0)
    ax_mat.text(scale_len / 2, scale_y + 0.20, '20 Mb', ha='center', va='bottom',
                fontsize=9.5, fontweight='bold', color='#0F172A')

    # Paternal scale bar
    ax_pat.plot([0, scale_len], [scale_y, scale_y], color='#0F172A', lw=2.2)
    ax_pat.plot([0, 0], [scale_y - 0.12, scale_y + 0.12], color='#0F172A', lw=2.0)
    ax_pat.plot([scale_len, scale_len], [scale_y - 0.12, scale_y + 0.12], color='#0F172A', lw=2.0)
    ax_pat.text(scale_len / 2, scale_y + 0.20, '20 Mb', ha='center', va='bottom',
                fontsize=9.5, fontweight='bold', color='#0F172A')

    # 7. Comprehensive Legend
    legend_elements = [
        # Ideogram Bars
        mpatches.Patch(facecolor=COL_T2T_BG, edgecolor=COL_T2T_BORDER, label='T2T Reference Bar'),
        mpatches.Patch(facecolor=COL_S_BG, edgecolor=COL_S_BORDER, label='HiFi Single Assembly'),
        mpatches.Patch(facecolor=COL_D_BG, edgecolor=COL_D_BORDER, label='HiFi Dual Assembly'),
        mpatches.Patch(facecolor=COL_O_BG, edgecolor=COL_O_BORDER, label='ONT Dual Assembly'),
        # Ribbons
        mpatches.Patch(facecolor=COL_RIB0_COLL, edgecolor=COL_RIB0_COLL_E, alpha=0.5, label='Collinear Synteny Ribbon'),
        mpatches.Patch(facecolor=COL_RIB0_NC, edgecolor=COL_RIB0_NC_E, alpha=0.7, label='Non-collinear / Inverted Ribbon'),
        # Native Features
        plt.Line2D([0], [0], color=COL_CUR_GAP, markeredgecolor=COL_CUR_GAP_ED, lw=2.2, label='Curation Gap (Join)'),
        plt.Line2D([0], [0], color=COL_ASM_GAP, lw=1.5, label='Assembly Gap'),
        mpatches.Patch(facecolor=COL_SWITCH, edgecolor=COL_SWITCH_ED, label='Hap-mer Switch Error'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_CENTRO, markeredgecolor='#581C87',
                   markersize=7, label='Centromere'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_TELO, markeredgecolor='#064E3B',
                   markersize=7, label='Telomere')
    ]

    fig.legend(handles=legend_elements, loc='lower center',
               bbox_to_anchor=(0.5, 0.015), ncol=6, frameon=True,
               facecolor='#F8FAFC', edgecolor='#E2E8F0', fontsize=8.8,
               handlelength=1.5, handleheight=0.9, columnspacing=1.6)

    # 8. Save
    print(f'Saving PNG: {out_png}')
    fig.savefig(out_png, dpi=dpi, facecolor='white', bbox_inches='tight')
    if out_pdf:
        print(f'Saving PDF: {out_pdf}')
        fig.savefig(out_pdf, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    print('Stacked synteny ideogram successfully generated!')


# =============================================================================
# CLI Main
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description='Stacked Multi-Alignment Macrochromosome Synteny Ideogram')
    parser.add_argument('--t2t-tsv', required=True)
    parser.add_argument('--single-tsv', required=True)
    parser.add_argument('--dual-tsv', required=True)
    parser.add_argument('--ont-tsv', required=True)

    parser.add_argument('--single-pairs', required=True)
    parser.add_argument('--dual-pairs', required=True)
    parser.add_argument('--ont-pairs', required=True)

    parser.add_argument('--single-chain', required=True)
    parser.add_argument('--single-nc-chain', required=True)
    parser.add_argument('--dual-chain', required=True)
    parser.add_argument('--dual-nc-chain', required=True)
    parser.add_argument('--ont-chain', required=True)
    parser.add_argument('--ont-nc-chain', required=True)

    parser.add_argument('--single-gaps', required=True)
    parser.add_argument('--dual-gaps', required=True)
    parser.add_argument('--ont-gaps', required=True)

    parser.add_argument('--single-bed', required=True)
    parser.add_argument('--dual-bed', required=True)
    parser.add_argument('--ont-bed', required=True)

    parser.add_argument('--single-telomeres', required=True)
    parser.add_argument('--dual-telomeres', required=True)
    parser.add_argument('--ont-telomeres', required=True)

    parser.add_argument('--centromeres', required=True)
    parser.add_argument('--telo-p-bed', required=True)

    parser.add_argument('--output-png', required=True)
    parser.add_argument('--output-pdf', default=None)
    parser.add_argument('--dpi', type=int, default=300)

    return parser.parse_args()


def main():
    args = parse_args()

    print('Loading sequence sizes...')
    t2t_sizes    = load_seq_sizes(args.t2t_tsv)
    s_query_sizes = load_seq_sizes(args.single_tsv)
    d_query_sizes = load_seq_sizes(args.dual_tsv)
    o_query_sizes = load_seq_sizes(args.ont_tsv)

    print('Loading chromosome pairings...')
    s_pairs = load_chrom_pairs(args.single_pairs)
    d_pairs = load_chrom_pairs(args.dual_pairs)
    o_pairs = load_chrom_pairs(args.ont_pairs)

    print('Loading centromeres and telomeres...')
    centromeres = load_centromeres(args.centromeres)
    p_starts    = load_p_arm_bed(args.telo_p_bed)
    flip_set    = build_flip_set_from_centromeres_and_telomeres(t2t_sizes, centromeres, p_starts)
    t2t_telo    = load_t2t_terminal_telomeres(args.telo_p_bed)

    print('Loading gaps, hap-mers, and assembly telomeres...')
    s_gaps = load_gaps_bed(args.single_gaps)
    d_gaps = load_gaps_bed(args.dual_gaps)
    o_gaps = load_gaps_bed(args.ont_gaps)

    s_sw = load_switch_blocks(args.single_bed)
    d_sw = load_switch_blocks(args.dual_bed)
    o_sw = load_switch_blocks(args.ont_bed)

    s_telo = load_telomere_presence_tsv(args.single_telomeres)
    d_telo = load_telomere_presence_tsv(args.dual_telomeres)
    o_telo = load_telomere_presence_tsv(args.ont_telomeres)

    print('Extracting synteny ribbons...')
    s_col_ribbons = extract_synteny_ribbons(args.single_chain, is_collinear=True)
    s_nc_ribbons  = extract_synteny_ribbons(args.single_nc_chain, is_collinear=False)

    d_col_ribbons = extract_synteny_ribbons(args.dual_chain, is_collinear=True)
    d_nc_ribbons  = extract_synteny_ribbons(args.dual_nc_chain, is_collinear=False)

    o_col_ribbons = extract_synteny_ribbons(args.ont_chain, is_collinear=True)
    o_nc_ribbons  = extract_synteny_ribbons(args.ont_nc_chain, is_collinear=False)

    print(f'  Single: {len(s_col_ribbons)} collinear, {len(s_nc_ribbons)} non-collinear ribbons')
    print(f'  Dual:   {len(d_col_ribbons)} collinear, {len(d_nc_ribbons)} non-collinear ribbons')
    print(f'  ONT:    {len(o_col_ribbons)} collinear, {len(o_nc_ribbons)} non-collinear ribbons')

    print('Building stacked multi-alignment synteny figure...')
    plot_butterfly_stacked(
        t2t_sizes=t2t_sizes,
        flip_set=flip_set,
        s_query_sizes=s_query_sizes,
        s_pairs=s_pairs,
        s_col_ribbons=s_col_ribbons,
        s_nc_ribbons=s_nc_ribbons,
        s_gaps=s_gaps,
        s_sw=s_sw,
        s_telo=s_telo,
        d_query_sizes=d_query_sizes,
        d_pairs=d_pairs,
        d_col_ribbons=d_col_ribbons,
        d_nc_ribbons=d_nc_ribbons,
        d_gaps=d_gaps,
        d_sw=d_sw,
        d_telo=d_telo,
        o_query_sizes=o_query_sizes,
        o_pairs=o_pairs,
        o_col_ribbons=o_col_ribbons,
        o_nc_ribbons=o_nc_ribbons,
        o_gaps=o_gaps,
        o_sw=o_sw,
        o_telo=o_telo,
        centromeres=centromeres,
        t2t_telo=t2t_telo,
        out_png=args.output_png,
        out_pdf=args.output_pdf,
        dpi=args.dpi
    )


if __name__ == '__main__':
    main()
