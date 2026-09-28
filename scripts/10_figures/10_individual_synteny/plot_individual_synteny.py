#!/usr/bin/env python3
"""
plot_individual_synteny.py

Generates ultra-detailed, publication-quality individual synteny QC figures for each
chromosome in the zebra finch (bTaeGut7) curation analysis.

Features per chromosome:
  - High-resolution butterfly synteny ideograms (Maternal on left, Paternal on right).
  - All 3 assembly technologies compared against T2T reference:
      1. HiFi Single-Assembly (Hic & Trio binned)
      2. HiFi Dual-Assembly (Hifiasm dual)
      3. ONT Dual-Assembly (Verkko/ONT dual)
  - Full annotation feature set:
      * Collinear ribbons (semi-transparent shaded)
      * Non-collinear synteny ribbons (amber)
      * Recovered / uncovered synteny ribbons (pink)
      * Centromeric constrictions (hourglass shape)
      * Telomere presence indicators (filled circle = collinear, open circle = non-collinear)
      * Assembly N-gaps (solid black bars) and curation join markers (crimson red lines)
      * Haplotype switch blocks (magenta highlights)
      * Unlocalized scaffolds (dashed borders with leader lines)
  - Dynamically scaled Mb coordinate axes and grid.
  - Informative header card and publication legend.
"""

import argparse
import os
import re
import math
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import matplotlib.ticker as mticker

# =============================================================================
# Constants & Colors
# =============================================================================

INSERTION_GAP_THRESHOLD = 20_000
MIN_NC_RIBBON_BP        = 10_000

COL_T2T          = '#E2E8F0'
COL_UNALIGNED    = '#FFFFFF'
COL_BORDER       = '#1E293B'
COL_CENTROMERE   = '#475569'
COL_ASM_GAP      = '#000000'
COL_CUR_GAP      = '#EF4444'
COL_SWITCH       = '#D946EF'
COL_UNLOC_EDGE   = '#64748B'
COL_RIB_COLL     = '#E2E8F0'
COL_RIB_NC       = '#F59E0B'
COL_RIB_REC      = '#FDA4AF'

COL_COV_S_DARK   = '#22D3EE'
COL_COV_S_LIGHT  = '#A5F3FC'
COL_COV_D_DARK   = '#818CF8'
COL_COV_D_LIGHT  = '#C7D2FE'
COL_COV_O_DARK   = '#4ADE80'
COL_COV_O_LIGHT  = '#BBF7D0'

MACRO_SET = {'1', '1A', '2', '3', '4', '4A', '5', '6', '7', '8', 'W', 'Z'}
MICRO_SET = {'9', '10', '11', '12', '13', '14', '15', '17', '18', '19', '20', '21', '22', '23', '24', '26', '27', '28'}
DOT_SET   = {'16', '25', '29', '30', '31', '32', '33', '34', '35', '36', '37'}

ALL_CHROM_ORDER = [
    '1', '1A', '2', '3', '4', '4A', '5', '6', '7', '8', 'W', 'Z',
    '9', '10', '11', '12', '13', '14', '15', '17', '18', '19', '20', '21', '22', '23', '24', '26', '27', '28',
    '16', '25', '29', '30', '31', '32', '33', '34', '35', '36', '37'
]

# =============================================================================
# Helper Utilities & Loaders
# =============================================================================

def _merge_intervals(intervals, min_gap=5_000):
    if not intervals:
        return []
    s_int = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [s_int[0]]
    for s, e in s_int[1:]:
        ls, le = merged[-1]
        if s <= le + min_gap:
            merged[-1] = (ls, max(le, e))
        else:
            merged.append((s, e))
    return merged

def load_seq_sizes(tsv_path):
    sizes = {}
    if not tsv_path or not os.path.exists(tsv_path):
        return sizes
    with open(tsv_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'): continue
            p = line.split()
            if len(p) >= 2 and p[0] != 'name' and p[1] != 'size':
                try:
                    sizes[p[0]] = int(p[1])
                except ValueError:
                    continue
    return sizes

def load_chrom_pairs(path):
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
                    ocs, oce = centromeres[key]
                    centromeres[key] = (min(ocs, cs), max(oce, ce))
    return centromeres

def load_p_arm_bed(bed_path):
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

def load_t2t_terminal_telomeres(bed_path):
    telo = defaultdict(set)
    if not bed_path or not os.path.exists(bed_path):
        return telo
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split('\t')
            if len(parts) >= 4:
                chrom, arm = parts[0], parts[3].lower()
                if arm in ('p', 'q'):
                    telo[chrom].add(arm)
    return telo

def build_flip_set_from_centromeres_and_telomeres(seq_sizes, centromeres, p_arm_dict):
    flip_set = set()
    for chrom, size in seq_sizes.items():
        m = re.match(r'^(Pat|Mat)_.*_chromosome_([0-9]+[A-Za-z]*|[A-Za-z]+)$', chrom)
        if not m:
            continue
        side, tok = m.group(1).lower(), m.group(2).upper()
        cen_span = centromeres.get((tok, side))
        if cen_span:
            cen_mid = (cen_span[0] + cen_span[1]) / 2
            if cen_mid > size / 2:
                flip_set.add(chrom)
        elif chrom in p_arm_dict:
            if p_arm_dict[chrom] > size / 2:
                flip_set.add(chrom)
    return flip_set

def load_gaps_bed(gaps_bed):
    gaps = defaultdict(list)
    if not gaps_bed or not os.path.exists(gaps_bed):
        return gaps
    with open(gaps_bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            p = line.split('\t')
            if len(p) >= 3:
                scaf, s, e = p[0], int(p[1]), int(p[2])
                gtype = p[3] if len(p) > 3 else 'GAP'
                gaps[scaf].append((s, e, gtype))
    return gaps

def load_switch_blocks(bed_path):
    sw = defaultdict(list)
    if not bed_path or not os.path.exists(bed_path):
        return sw
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            p = line.split('\t')
            if len(p) >= 3:
                sw[p[0]].append((int(p[1]), int(p[2])))
    return sw

def load_telomere_presence_tsv(tsv_path):
    telo_coll = defaultdict(set)
    telo_nc   = defaultdict(set)
    if not tsv_path or not os.path.exists(tsv_path):
        return telo_coll, telo_nc
    with open(tsv_path) as f:
        hdr = None
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            p = line.split('\t')
            if not hdr:
                hdr = [x.lower() for x in p]
                continue
            row = dict(zip(hdr, p))
            chrom = row.get('chromosome') or row.get('chrom') or row.get('name')
            scaff = row.get('scaffold')
            col_val = (row.get('collinear') or '').lower()
            nc_val = (row.get('non-collinear') or row.get('noncollinear') or '').lower()

            arms_c = set()
            if col_val and col_val != 'none':
                if 'p' in col_val: arms_c.add('p')
                if 'q' in col_val: arms_c.add('q')

            arms_nc = set()
            if nc_val and nc_val != 'none':
                if 'p' in nc_val: arms_nc.add('p')
                if 'q' in nc_val: arms_nc.add('q')

            for key in (chrom, scaff):
                if key:
                    telo_coll[key].update(arms_c)
                    telo_nc[key].update(arms_nc)
    return telo_coll, telo_nc

def parse_chain_detailed(chain_path, is_collinear=True, min_nc_size=MIN_NC_RIBBON_BP):
    ribbons = defaultdict(list)
    query_spans = defaultdict(list)
    insertions = defaultdict(list)

    if not chain_path or not os.path.exists(chain_path):
        return ribbons, query_spans, insertions

    def _add_ribbon(t_n, ts, te, q_n, qs, qe, qstr):
        if te <= ts or qe <= qs:
            return
        if not is_collinear and max(te - ts, qe - qs) < min_nc_size:
            return
        ribbons[t_n].append((t_n, ts, te, q_n, qs, qe, qstr, is_collinear))
        query_spans[q_n].append((qs, qe))

    with open(chain_path) as f:
        t_name = q_name = q_strand = None
        t_pos = q_pos = q_size = 0
        cur_t_start = cur_t_end = cur_q_start = cur_q_end = 0
        in_ribbon = False

        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            if line.startswith('chain'):
                if in_ribbon:
                    _add_ribbon(t_name, cur_t_start, cur_t_end, q_name, cur_q_start, cur_q_end, q_strand)
                in_ribbon = False

                parts = line.split()
                t_name, q_name = parts[2], parts[7]
                q_strand = parts[9]
                t_pos = int(parts[5])
                q_size = int(parts[8])
                q_raw_start = int(parts[10])
                q_raw_end   = int(parts[11])

                if q_strand == '-':
                    q_pos = q_size - q_raw_end
                else:
                    q_pos = q_raw_start
            else:
                parts = line.split()
                size = int(parts[0])
                dt = int(parts[1]) if len(parts) > 1 else 0
                dq = int(parts[2]) if len(parts) > 2 else 0

                blk_t_s = t_pos
                blk_t_e = t_pos + size
                blk_q_s = q_pos
                blk_q_e = q_pos + size

                if not in_ribbon:
                    cur_t_start, cur_t_end = blk_t_s, blk_t_e
                    cur_q_start, cur_q_end = blk_q_s, blk_q_e
                    in_ribbon = True
                else:
                    cur_t_end = blk_t_e
                    cur_q_end = blk_q_e

                t_pos += size + dt
                q_pos += size + dq

                if (dq >= INSERTION_GAP_THRESHOLD or dt >= INSERTION_GAP_THRESHOLD) and in_ribbon:
                    _add_ribbon(t_name, cur_t_start, cur_t_end, q_name, cur_q_start, cur_q_end, q_strand)
                    if dq >= INSERTION_GAP_THRESHOLD:
                        insertions[q_name].append((blk_q_e, blk_q_e + dq, dq, blk_t_e))
                    in_ribbon = False

        if in_ribbon:
            _add_ribbon(t_name, cur_t_start, cur_t_end, q_name, cur_q_start, cur_q_end, q_strand)

    return ribbons, query_spans, insertions

def load_recovered_chains(rec_chain_path):
    rec_ribbons = defaultdict(list)
    rec_spans = defaultdict(list)
    sizes = {}
    if not rec_chain_path or not os.path.exists(rec_chain_path):
        return rec_ribbons, rec_spans, sizes
    with open(rec_chain_path) as f:
        t_name = q_name = q_strand = None
        t_pos = q_pos = 0
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            if line.startswith('chain'):
                p = line.split()
                t_name, q_name, q_strand = p[2], p[7], p[9]
                t_pos, q_pos = int(p[5]), int(p[10])
                sizes[q_name] = int(p[8])
            else:
                p = line.split()
                size = int(p[0])
                dt = int(p[1]) if len(p) > 1 else 0
                dq = int(p[2]) if len(p) > 2 else 0
                ts, te = t_pos, t_pos + size
                qs, qe = q_pos, q_pos + size
                rec_ribbons[t_name].append((t_name, ts, te, q_name, qs, qe, q_strand, False))
                rec_spans[q_name].append((qs, qe))
                t_pos += size + dt
                q_pos += size + dq
    return rec_ribbons, rec_spans, sizes

def classify_chrom_group(token):
    tok = str(token).upper().replace('CHR', '')
    if tok in MACRO_SET: return 'Macrochromosome'
    if tok in MICRO_SET: return 'Microchromosome'
    if tok in DOT_SET:   return 'Dot chromosome'
    return 'Dot chromosome'

# =============================================================================
# Path & Drawing Primitives
# =============================================================================

def _rect_path(x0, x1, yc, h):
    y0, y1 = yc - h / 2, yc + h / 2
    verts = [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    return Path(verts, codes)

def _hourglass_path(x0, x1, yc, h, cs, ce):
    y0, y1 = yc - h / 2, yc + h / 2
    waist_h = h * 0.22
    y_top_in = yc + waist_h / 2
    y_bot_in = yc - waist_h / 2
    verts = [
        (x0, y0), (cs, y0), (cs, y_bot_in), (ce, y_bot_in), (ce, y0),
        (x1, y0), (x1, y1), (ce, y1), (ce, y_top_in), (cs, y_top_in),
        (cs, y1), (x0, y1), (x0, y0)
    ]
    codes = [
        Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
        Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
        Path.LINETO, Path.LINETO, Path.CLOSEPOLY
    ]
    return Path(verts, codes)

def _draw_telo_semi(ax, x_edge, tip_type, y_center, h_bar, rx, fc='black', ec='black', lw=0.4):
    y_top = y_center + h_bar / 2
    y_bot = y_center - h_bar / 2
    x_ext = (x_edge - rx) if tip_type == 'p' else (x_edge + rx)
    verts = [
        (x_edge, y_top),
        (x_ext, y_top), (x_ext, y_bot), (x_edge, y_bot),
        (x_edge, y_top)
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4, Path.CURVE4, Path.LINETO,
        Path.CLOSEPOLY
    ]
    patch = PathPatch(Path(verts, codes), facecolor=fc, edgecolor=ec, linewidth=lw, zorder=9)
    ax.add_patch(patch)

def draw_bezier_ribbon(ax, x0_top, x1_top, y_top, x0_bot, x1_bot, y_bot,
                       color, alpha=0.38, edge_color=None, edge_width=0.4):
    if x1_top < x0_top: x0_top, x1_top = x1_top, x0_top
    if x1_bot < x0_bot: x0_bot, x1_bot = x1_bot, x0_bot
    if x1_top == x0_top and x1_bot == x0_bot: return

    ymid = (y_top + y_bot) / 2.0
    verts = [
        (x0_top, y_top),
        (x0_top, ymid), (x0_bot, ymid), (x0_bot, y_bot),
        (x1_bot, y_bot),
        (x1_bot, ymid), (x1_top, ymid), (x1_top, y_top),
        (x0_top, y_top)
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.CLOSEPOLY
    ]
    patch = PathPatch(Path(verts, codes), facecolor=color, edgecolor=edge_color if edge_color else 'none',
                      linewidth=edge_width if edge_color else 0, alpha=alpha, zorder=3)
    ax.add_patch(patch)

# =============================================================================
# Full Synteny Ideogram Renderer
# =============================================================================

def plot_single_chromosome_synteny(tok, max_global_size, t2t_sizes, centromeres, t2t_telo,
                                   p_arm_dict, flip_set, single_data, dual_data, ont_data,
                                   output_path, dpi=300):
    group_name = classify_chrom_group(tok)

    # Dynamic tick step based on chromosome size
    if max_global_size >= 40_000_000:
        tick_step = 10_000_000
    elif max_global_size >= 15_000_000:
        tick_step = 5_000_000
    elif max_global_size >= 5_000_000:
        tick_step = 2_000_000
    elif max_global_size >= 2_000_000:
        tick_step = 1_000_000
    else:
        tick_step = 500_000

    fig_w = 15.0
    fig_h = 7.5
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')

    outer = GridSpec(3, 1, figure=fig, height_ratios=[0.14, 0.74, 0.12], hspace=0.18,
                     left=0.04, right=0.96, top=0.95, bottom=0.06)

    # ── Header Banner ────────────────────────────────────────────────────────
    ax_head = fig.add_subplot(outer[0, 0])
    ax_head.axis('off')
    
    head_rect = mpatches.FancyBboxPatch((0.0, 0.05), 1.0, 0.90, boxstyle="round,pad=0.02,rounding_size=0.04",
                                        facecolor='#F8FAFC', edgecolor='#E2E8F0', linewidth=1.2,
                                        transform=ax_head.transAxes, zorder=1)
    ax_head.add_patch(head_rect)

    t2t_m_chrom = next((k for k in t2t_sizes.keys() if f'chromosome_{tok}' in k and ('Mat' in k or tok == 'W')), None)
    t2t_p_chrom = next((k for k in t2t_sizes.keys() if f'chromosome_{tok}' in k and ('Pat' in k or tok == 'Z')), None)
    t2t_m_len = (t2t_sizes.get(t2t_m_chrom, 0) / 1e6) if t2t_m_chrom else 0.0
    t2t_p_len = (t2t_sizes.get(t2t_p_chrom, 0) / 1e6) if t2t_p_chrom else 0.0

    ax_head.text(0.02, 0.65, f"Chromosome {tok} Synteny & Assembly QC Ideogram",
                 fontsize=15, fontweight='bold', color='#0F172A', va='center', transform=ax_head.transAxes)
    ax_head.text(0.02, 0.28, f"Classification: {group_name}   |   T2T Reference: Maternal {t2t_m_len:.2f} Mb, Paternal {t2t_p_len:.2f} Mb",
                 fontsize=11, color='#475569', va='center', transform=ax_head.transAxes)

    # ── Synteny Butterfly Area ───────────────────────────────────────────────
    inner_syn = GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[1, 0],
                                        width_ratios=[0.08, 0.46, 0.46], wspace=0.015)

    ax_chr = fig.add_subplot(inner_syn[0, 0])
    ax_mat = fig.add_subplot(inner_syn[0, 1])
    ax_pat = fig.add_subplot(inner_syn[0, 2])

    ax_chr.axis('off')

    H_BAR_T2T  = 0.038
    H_BAR_ASM  = 0.038
    H_RIBBON   = 0.088
    BLOCK_GAP  = 0.055
    BLOCK_H    = H_BAR_T2T + H_RIBBON + H_BAR_ASM + BLOCK_GAP

    total_h = 3 * BLOCK_H

    for ax in (ax_mat, ax_pat):
        ax.set_ylim(0, total_h)
        ax.set_facecolor('white')

    ax_mat.set_title('Maternal', fontsize=12.5, fontweight='bold', pad=8, color='#0F172A')
    ax_pat.set_title('Paternal', fontsize=12.5, fontweight='bold', pad=8, color='#0F172A')

    axis_max = int(np.ceil(max_global_size / tick_step)) * tick_step
    eff_max = max(max_global_size, axis_max)
    x_margin = eff_max * 0.008

    ax_chr.set_xlim(0, 1)
    ax_mat.set_xlim(eff_max + x_margin, -x_margin)
    ax_pat.set_xlim(-x_margin, eff_max + x_margin)

    telo_rx = eff_max * 0.0055

    # Shading background
    ax_mat.axhspan(0, total_h, color='#F8FAFC', zorder=0)
    ax_pat.axhspan(0, total_h, color='#F8FAFC', zorder=0)

    assemblies = [
        ('HiFi Single', single_data, COL_COV_S_DARK, COL_COV_S_LIGHT),
        ('HiFi Dual',   dual_data,   COL_COV_D_DARK, COL_COV_D_LIGHT),
        ('ONT Dual',    ont_data,    COL_COV_O_DARK, COL_COV_O_LIGHT)
    ]

    for b_idx, (asm_title, asm, c_dark, c_light) in enumerate(assemblies):
        y_block_top = total_h - 0.08 - b_idx * BLOCK_H
        y_t2t = y_block_top - H_BAR_T2T / 2
        y_asm = y_t2t - (H_BAR_T2T / 2 + H_RIBBON + H_BAR_ASM / 2)

        # Explicit track labels for T2T Ref and Assembly
        ax_chr.text(0.95, y_t2t, "T2T Ref", ha='right', va='center',
                    fontsize=8.5, fontweight='bold', color='#64748B')
        ax_chr.text(0.95, y_asm, asm_title, ha='right', va='center',
                    fontsize=9.0, fontweight='bold', color=c_dark)

        for side, ax in [('mat', ax_mat), ('pat', ax_pat)]:
            # Resolve T2T chromosome
            if tok == 'W':
                if side != 'mat': continue
                t2t_chrom = 'Mat_NC_133064.1_chromosome_W'
            elif tok == 'Z':
                if side != 'pat': continue
                t2t_chrom = 'Mat_NC_133063.1_chromosome_Z'
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

            # Primary scaffold
            q_prim = asm['pairs'].get(t2t_chrom)
            q_prim_len = asm['sizes'].get(q_prim, t2t_len) if q_prim else t2t_len

            # Unlinked scaffolds
            all_asm_ribs = asm['col_ribbons'].get(t2t_chrom, []) + asm['nc_ribbons'].get(t2t_chrom, [])
            unloc_aligned = defaultdict(int)
            for r in all_asm_ribs:
                qn = r[3]
                if qn != q_prim:
                    unloc_aligned[qn] += abs(r[5] - r[4])

            unloc_candidates = set(unloc_aligned.keys())
            for r in asm.get('rec_ribbons', {}).get(t2t_chrom, []):
                qn = r[3]
                if qn != q_prim:
                    unloc_candidates.add(qn)

            if q_prim:
                parts = q_prim.rsplit('.', 1)
                stem = parts[0]
                suffix = f".{parts[1]}" if len(parts) > 1 else ""
                for s in asm['sizes']:
                    if s.startswith(f"{stem}_unloc_") and (not suffix or s.endswith(suffix)):
                        unloc_candidates.add(s)

            def _unloc_sort_key(u):
                m = re.search(r'unloc_(\d+)', u)
                u_num = int(m.group(1)) if m else 999999
                align_len = unloc_aligned.get(u, 0)
                tot_len = asm['sizes'].get(u, 0)
                return (0 if align_len > 0 else 1, u_num, -tot_len)

            unloc_sorted = sorted(unloc_candidates, key=_unloc_sort_key)

            gap_unloc = max(eff_max * 0.012, 120_000)
            scaffold_offsets = {q_prim: 0}
            unloc_offsets = {}
            curr_x = q_prim_len

            for u in unloc_sorted:
                curr_x += gap_unloc
                u_len = asm['sizes'].get(u, 0)
                scaffold_offsets[u] = curr_x
                unloc_offsets[u] = curr_x
                curr_x += u_len

            # ── A. T2T Reference Track ───────────────────────────────────────
            c_key = (tok.upper(), 'pat' if (side == 'pat' or tok == 'Z') else 'mat')
            cen = centromeres.get(c_key)
            if cen:
                cs, ce = (t2t_len - cen[1], t2t_len - cen[0]) if flip else (cen[0], cen[1])
                cs, ce = max(0, min(cs, ce)), min(t2t_len, max(cs, ce))
                t2t_path = _hourglass_path(0, t2t_len, y_t2t, H_BAR_T2T, cs, ce)
            else:
                t2t_path = _rect_path(0, t2t_len, y_t2t, H_BAR_T2T)

            ax.add_patch(PathPatch(t2t_path, fc=COL_T2T, ec=COL_BORDER, lw=0.55, zorder=5))

            # Collinear coverage on T2T
            t2t_col = _merge_intervals(asm['col_spans_t2t'].get(t2t_chrom, []))
            for s, e in t2t_col:
                x0 = (t2t_len - e) if flip else s
                x1 = (t2t_len - s) if flip else e
                r = mpatches.Rectangle((min(x0, x1), y_t2t - H_BAR_T2T / 2),
                                       abs(x1 - x0), H_BAR_T2T,
                                       fc=c_dark, ec='none', alpha=1.00, zorder=7)
                r.set_clip_path(t2t_path, transform=ax.transData)
                ax.add_patch(r)

            # T2T Telomeres
            t2t_arms = t2t_telo.get(t2t_chrom, set())
            if p_arm in t2t_arms or not t2t_arms:
                _draw_telo_semi(ax, 0, 'p', y_t2t, H_BAR_T2T, telo_rx, fc='#000000', ec='#000000')
            if q_arm in t2t_arms or not t2t_arms:
                _draw_telo_semi(ax, t2t_len, 'q', y_t2t, H_BAR_T2T, telo_rx, fc='#000000', ec='#000000')

            # ── B. Assembly Ideogram Track ───────────────────────────────────
            if q_prim:
                prim_path = _rect_path(0, q_prim_len, y_asm, H_BAR_ASM)
                ax.add_patch(PathPatch(prim_path, fc=COL_UNALIGNED, ec=COL_BORDER, lw=0.55, zorder=5))

                # Non-collinear coverage
                m_nc = _merge_intervals(asm['nc_spans'].get(q_prim, []))
                for s, e in m_nc:
                    x0 = (q_prim_len - e) if flip else s
                    x1 = (q_prim_len - s) if flip else e
                    r = mpatches.Rectangle((min(x0, x1), y_asm - H_BAR_ASM / 2),
                                           abs(x1 - x0), H_BAR_ASM,
                                           fc=c_dark, ec='none', alpha=0.25, zorder=6)
                    r.set_clip_path(prim_path, transform=ax.transData)
                    ax.add_patch(r)

                # Collinear coverage
                m_col = _merge_intervals(asm['col_spans'].get(q_prim, []))
                for s, e in m_col:
                    x0 = (q_prim_len - e) if flip else s
                    x1 = (q_prim_len - s) if flip else e
                    r = mpatches.Rectangle((min(x0, x1), y_asm - H_BAR_ASM / 2),
                                           abs(x1 - x0), H_BAR_ASM,
                                           fc=c_dark, ec='none', alpha=1.00, zorder=7)
                    r.set_clip_path(prim_path, transform=ax.transData)
                    ax.add_patch(r)

                # Switch errors
                for sw_s, sw_e in asm['sw'].get(q_prim, []):
                    x0 = (q_prim_len - sw_e) if flip else sw_s
                    x1 = (q_prim_len - sw_s) if flip else sw_e
                    w_sw = max(abs(x1 - x0), 80_000)
                    r = mpatches.Rectangle((min(x0, x1), y_asm - H_BAR_ASM / 2),
                                           w_sw, H_BAR_ASM,
                                           fc=COL_SWITCH, ec='none', zorder=8)
                    r.set_clip_path(prim_path, transform=ax.transData)
                    ax.add_patch(r)

                # Gaps & Joins
                for gs, ge, gtype in asm['gaps'].get(q_prim, []):
                    gx = (q_prim_len - (gs + ge) / 2) if flip else ((gs + ge) / 2)
                    gap_col = COL_CUR_GAP if 'CURATION' in gtype else COL_ASM_GAP
                    ax.plot([gx, gx], [y_asm - H_BAR_ASM / 2, y_asm + H_BAR_ASM / 2],
                            color=gap_col, lw=0.85, zorder=9)

                # Telomeres
                t_col_p, t_col_q = asm['telo_col'].get(q_prim, set()), asm['telo_col'].get(q_prim, set())
                t_nc_p,  t_nc_q  = asm['telo_nc'].get(q_prim, set()),  asm['telo_nc'].get(q_prim, set())
                if p_arm in t_col_p:
                    _draw_telo_semi(ax, 0, 'p', y_asm, H_BAR_ASM, telo_rx, fc='#000000', ec='#000000')
                elif p_arm in t_nc_p:
                    _draw_telo_semi(ax, 0, 'p', y_asm, H_BAR_ASM, telo_rx, fc='#FFFFFF', ec='#000000', lw=0.6)

                if q_arm in t_col_q:
                    _draw_telo_semi(ax, q_prim_len, 'q', y_asm, H_BAR_ASM, telo_rx, fc='#000000', ec='#000000')
                elif q_arm in t_nc_q:
                    _draw_telo_semi(ax, q_prim_len, 'q', y_asm, H_BAR_ASM, telo_rx, fc='#FFFFFF', ec='#000000', lw=0.6)

            # ── C. Unlocalized Scaffolds ─────────────────────────────────────
            for u in unloc_sorted:
                ux0 = unloc_offsets[u]
                u_len = asm['sizes'].get(u, 0)
                ux1 = ux0 + u_len
                u_path = _rect_path(ux0, ux1, y_asm, H_BAR_ASM)
                ax.add_patch(PathPatch(u_path, fc=COL_UNALIGNED, ec=COL_UNLOC_EDGE, ls='--', lw=0.6, zorder=5))

                # Unloc coverage & features
                for s, e in _merge_intervals(asm['col_spans'].get(u, [])):
                    r = mpatches.Rectangle((ux0 + s, y_asm - H_BAR_ASM / 2), abs(e - s), H_BAR_ASM,
                                           fc=c_dark, ec='none', alpha=1.00, zorder=7)
                    r.set_clip_path(u_path, transform=ax.transData)
                    ax.add_patch(r)

                for sw_s, sw_e in asm['sw'].get(u, []):
                    w_sw = max(abs(sw_e - sw_s), 40_000)
                    r = mpatches.Rectangle((ux0 + sw_s, y_asm - H_BAR_ASM / 2), w_sw, H_BAR_ASM,
                                           fc=COL_SWITCH, ec='none', zorder=8)
                    r.set_clip_path(u_path, transform=ax.transData)
                    ax.add_patch(r)

                # Unloc label
                m_u = re.search(r'unloc_(\d+)', u)
                u_tag = f"u{m_u.group(1)}" if m_u else "u"
                ax.text((ux0 + ux1) / 2, y_asm + H_BAR_ASM / 2 + 0.012, u_tag,
                        fontsize=7.5, color='#64748B', ha='center', va='bottom', rotation=45)

            # ── D. Ribbons ───────────────────────────────────────────────────
            y_rib_top = y_t2t - H_BAR_T2T / 2
            y_rib_bot = y_asm + H_BAR_ASM / 2

            # Collinear Ribbons
            for r in asm['col_ribbons'].get(t2t_chrom, []):
                _, ts, te, qn, qs, qe, qstr, _ = r
                if qn not in scaffold_offsets: continue
                q_off = scaffold_offsets[qn]
                q_len = asm['sizes'].get(qn, q_prim_len)
                x0_top = (t2t_len - te) if flip else ts
                x1_top = (t2t_len - ts) if flip else te
                x0_bot = q_off + ((q_len - qe) if flip else qs)
                x1_bot = q_off + ((q_len - qs) if flip else qe)
                draw_bezier_ribbon(ax, x0_top, x1_top, y_rib_top, x0_bot, x1_bot, y_rib_bot,
                                   color=c_light, alpha=0.60, edge_color=c_dark, edge_width=0.3)

            # Non-collinear Ribbons
            for r in asm['nc_ribbons'].get(t2t_chrom, []):
                _, ts, te, qn, qs, qe, qstr, _ = r
                if qn not in scaffold_offsets: continue
                q_off = scaffold_offsets[qn]
                q_len = asm['sizes'].get(qn, q_prim_len)
                x0_top = (t2t_len - te) if flip else ts
                x1_top = (t2t_len - ts) if flip else te
                x0_bot = q_off + ((q_len - qe) if flip else qs)
                x1_bot = q_off + ((q_len - qs) if flip else qe)
                draw_bezier_ribbon(ax, x0_top, x1_top, y_rib_top, x0_bot, x1_bot, y_rib_bot,
                                   color=COL_RIB_NC, alpha=0.55, edge_color='#D97706', edge_width=0.4)

            # Uncovered Synteny Ribbons
            for r in asm.get('rec_ribbons', {}).get(t2t_chrom, []):
                _, ts, te, qn, qs, qe, qstr, _ = r
                if qn not in scaffold_offsets: continue
                q_off = scaffold_offsets[qn]
                q_len = asm['sizes'].get(qn, q_prim_len)
                x0_top = (t2t_len - te) if flip else ts
                x1_top = (t2t_len - ts) if flip else te
                x0_bot = q_off + ((q_len - qe) if flip else qs)
                x1_bot = q_off + ((q_len - qs) if flip else qe)
                draw_bezier_ribbon(ax, x0_top, x1_top, y_rib_top, x0_bot, x1_bot, y_rib_bot,
                                   color=COL_RIB_REC, alpha=0.55, edge_color='#F43F5E', edge_width=0.4)

    # Format Axes
    for ax in (ax_mat, ax_pat):
        major_ticks = np.arange(0, eff_max + tick_step, tick_step)
        ax.set_xticks(major_ticks)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, p: f"{v / 1e6:.0f}" if v % 1e6 == 0 else f"{v / 1e6:.1f}"))
        ax.tick_params(axis='x', labelsize=10.0, colors='#334155', length=4.5, width=1.0)
        ax.tick_params(axis='y', left=False, labelleft=False)
        for spine in ['top', 'left', 'right']:
            ax.spines[spine].set_visible(False)
        ax.spines['bottom'].set_color('#94A3B8')

    # ── Bottom Legend Card ───────────────────────────────────────────────────
    ax_leg = fig.add_subplot(outer[2, 0])
    ax_leg.axis('off')

    leg_rect = mpatches.FancyBboxPatch((0.0, 0.0), 1.0, 1.0, boxstyle="round,pad=0.02,rounding_size=0.04",
                                       facecolor='#F8FAFC', edgecolor='#E2E8F0', linewidth=1.0,
                                       transform=ax_leg.transAxes, zorder=1)
    ax_leg.add_patch(leg_rect)

    patches = [
        mpatches.Patch(facecolor=COL_COV_S_DARK, edgecolor=COL_COV_S_DARK, label='HiFi Single'),
        mpatches.Patch(facecolor=COL_COV_D_DARK, edgecolor=COL_COV_D_DARK, label='HiFi Dual'),
        mpatches.Patch(facecolor=COL_COV_O_DARK, edgecolor=COL_COV_O_DARK, label='ONT Dual'),
        mpatches.Patch(facecolor='#E2E8F0', label='Collinear Synteny'),
        mpatches.Patch(facecolor=COL_RIB_NC, label='Non-Collinear'),
        mpatches.Patch(facecolor=COL_RIB_REC, label='Uncovered Synteny'),
        mpatches.Patch(facecolor=COL_SWITCH, label='Switch Block'),
        plt.Line2D([0], [0], color=COL_ASM_GAP, lw=2, label='Assembly Gap'),
        plt.Line2D([0], [0], color=COL_CUR_GAP, lw=2, label='Curation Join'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#000000', markersize=6, label='Telomere (Coll)'),
        plt.Line2D([0], [0], marker='o', color='k', markerfacecolor='#FFFFFF', markersize=6, label='Telomere (NC)'),
    ]

    ax_leg.legend(handles=patches, loc='center', ncol=6, frameon=False, fontsize=9.0, columnspacing=1.3)
    ax_leg.text(0.015, 0.5, "Scale: Mb", fontsize=10.0, fontweight='bold', color='#475569',
                va='center', transform=ax_leg.transAxes)

    png_path = output_path.replace('.pdf', '.png')
    pdf_path = output_path.replace('.png', '.pdf')

    plt.savefig(png_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.savefig(pdf_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    return png_path, pdf_path

# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Generate ultra-detailed individual chromosome synteny plots.")
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
    parser.add_argument('--single-rec-chain', required=True)
    parser.add_argument('--dual-rec-chain', required=True)
    parser.add_argument('--ont-rec-chain', required=True)
    parser.add_argument('--single-gaps', required=True)
    parser.add_argument('--dual-gaps', required=True)
    parser.add_argument('--ont-gaps', required=True)
    parser.add_argument('--single-bed', required=True)
    parser.add_argument('--dual-bed', required=True)
    parser.add_argument('--ont-bed', required=True)
    parser.add_argument('--single-telomeres', required=True)
    parser.add_argument('--dual-telomeres', required=True)
    parser.add_argument('--ont-telomeres', required=True)
    parser.add_argument('--telo-p-bed', required=True)
    parser.add_argument('--centromeres', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--chromosomes', default='all')
    parser.add_argument('--dpi', type=int, default=300)

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print("Loading comprehensive genomic datasets...")
    t2t_sizes = load_seq_sizes(args.t2t_tsv)
    p_arm_dict = load_p_arm_bed(args.telo_p_bed)
    centromeres = load_centromeres(args.centromeres)
    t2t_telo = load_t2t_terminal_telomeres(args.telo_p_bed)
    flip_set = build_flip_set_from_centromeres_and_telomeres(t2t_sizes, centromeres, p_arm_dict)

    def _load_asm(sizes_tsv, pairs_file, col_chain, nc_chain, rec_chain, gaps_bed, sw_bed, telo_tsv):
        sizes = load_seq_sizes(sizes_tsv)
        pairs = load_chrom_pairs(pairs_file)
        c_rib, c_spans, c_ins = parse_chain_detailed(col_chain, is_collinear=True)
        nc_rib, nc_spans, _ = parse_chain_detailed(nc_chain, is_collinear=False)
        rec_rib, rec_spans, rec_sz = load_recovered_chains(rec_chain)
        sizes.update(rec_sz)
        
        c_spans_t2t = defaultdict(list)
        for t_n, r_list in c_rib.items():
            for r in r_list:
                c_spans_t2t[t_n].append((r[1], r[2]))

        gaps = load_gaps_bed(gaps_bed)
        sw = load_switch_blocks(sw_bed)
        telo_c, telo_n = load_telomere_presence_tsv(telo_tsv)

        return {
            'sizes': sizes, 'pairs': pairs, 'col_ribbons': c_rib, 'col_spans': c_spans,
            'col_spans_t2t': c_spans_t2t, 'nc_ribbons': nc_rib, 'nc_spans': nc_spans,
            'rec_ribbons': rec_rib, 'rec_spans': rec_spans, 'gaps': gaps, 'sw': sw,
            'telo_col': telo_c, 'telo_nc': telo_n
        }

    single_data = _load_asm(args.single_tsv, args.single_pairs, args.single_chain, args.single_nc_chain,
                            args.single_rec_chain, args.single_gaps, args.single_bed, args.single_telomeres)
    dual_data   = _load_asm(args.dual_tsv, args.dual_pairs, args.dual_chain, args.dual_nc_chain,
                            args.dual_rec_chain, args.dual_gaps, args.dual_bed, args.dual_telomeres)
    ont_data    = _load_asm(args.ont_tsv, args.ont_pairs, args.ont_chain, args.ont_nc_chain,
                            args.ont_rec_chain, args.ont_gaps, args.ont_bed, args.ont_telomeres)

    if args.chromosomes.lower() == 'all':
        chrom_list = ALL_CHROM_ORDER
    else:
        tokens_req = [c.strip().replace('chr', '') for c in args.chromosomes.split(',')]
        chrom_list = [t for t in tokens_req if t]

    print(f"Generating detailed synteny plots for {len(chrom_list)} chromosomes...")

    for tok in chrom_list:
        # Calculate max size across assemblies
        chrom_stems = []
        if tok == 'W':
            chrom_stems = ['Mat_NC_133064.1_chromosome_W']
        elif tok == 'Z':
            chrom_stems = ['Mat_NC_133063.1_chromosome_Z']
        else:
            chrom_stems = [k for k in t2t_sizes.keys() if f'chromosome_{tok}' in k]

        max_s = max([t2t_sizes.get(c, 0) for c in chrom_stems] or [0])

        for asm in [single_data, dual_data, ont_data]:
            for st in chrom_stems:
                prim = asm['pairs'].get(st)
                if prim and prim in asm['sizes']:
                    tot = asm['sizes'][prim]
                    parts = prim.rsplit('.', 1)
                    stem = parts[0]
                    suffix = f".{parts[1]}" if len(parts) > 1 else ""
                    for s in asm['sizes']:
                        if s.startswith(f"{stem}_unloc_") and (not suffix or s.endswith(suffix)):
                            tot += max(max_s * 0.012, 120_000) + asm['sizes'][s]
                    max_s = max(max_s, tot)

        if max_s <= 0:
            max_s = 2_000_000

        out_base = os.path.join(args.output_dir, f"chr{tok}_synteny_qc.png")
        try:
            png_res, pdf_res = plot_single_chromosome_synteny(
                tok, max_s, t2t_sizes, centromeres, t2t_telo, p_arm_dict, flip_set,
                single_data, dual_data, ont_data, out_base, dpi=args.dpi
            )
            print(f"  [✓] Generated: chr{tok} -> {os.path.basename(png_res)}")
        except Exception as e:
            print(f"  [✗] Error processing chr{tok}: {e}")

    print(f"All individual chromosome synteny plots saved in: {args.output_dir}")

if __name__ == '__main__':
    main()
