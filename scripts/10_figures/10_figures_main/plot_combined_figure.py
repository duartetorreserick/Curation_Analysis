#!/usr/bin/env python3
"""
plot_combined_figure.py

Publication-quality multi-panel figure sized at half-height (fig_h = 12.0 in):
  - Row 0 (Compact Top Stats, ~15% height):
      Panel a: Collinear coverage distribution (reduced height by 40%)
      Panel b: Telomere completeness lollipop
      Panel c: Average genome coverage bars across all chromosomes (compact height, unified non-collinear yellow)
  - Rows 1–3 (Extended Ideograms, ~70% height total, ~23% each):
      Panel d: Macrochromosomes (Chr 4, 5, W maternal; Chr 4, 5 paternal)
      Panel e: Microchromosomes (Chr 11, 14, 18 butterfly)
      Panel f: Dot chromosomes (Chr 32, 34, 35 butterfly)
  - Row 4 (Half-height Bottom Row, ~15% height):
      Panel g (Left ~44% width): Gene-model error rates (Frameshifts % & Premature stop codons % with PCG numbers on top)
      Legend  (Right ~56% width): Neater, tidier publication legend card with 'Curation Join' and unified colors
"""

import argparse
import re
import os
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy.stats import gaussian_kde
import matplotlib.ticker as mticker

# =============================================================================
# Constants & Palette (Exact Supplementary Standard)
# =============================================================================

INSERTION_GAP_THRESHOLD = 20_000   # 20 kb gap in query (dq) splits collinear ribbon
MIN_NC_RIBBON_BP        = 20_000   # 20 kb minimum alignment to display non-collinear synteny ribbon

# Ideogram Base Colors
COL_BORDER       = '#2c3e50'
COL_T2T_FILL     = '#FFFFFF'       # T2T Reference is pure white / blank
COL_UNALIGNED    = '#FFFFFF'       # Unaligned assembly background is pure white

# Assembly Coverage Colors (Dark = Collinear, Light = Non-Collinear)
COL_COV_S_DARK   = '#00BCCC'       # Single Collinear (cyan)
COL_COV_S_LIGHT  = '#A0D4DC'       # Single Non-Collinear (light cyan)
COL_COV_D_DARK   = '#6445B0'       # Dual Collinear (purple)
COL_COV_D_LIGHT  = '#C4B8E8'       # Dual Non-Collinear (light purple)
COL_COV_O_DARK   = '#2E8B52'       # ONT Collinear (green)
COL_COV_O_LIGHT  = '#B8E8CC'       # ONT Non-Collinear (light green)
COL_UNCOV        = '#FFFFFF'
COL_ROW_ODD      = '#F2F9FD'

# Native Feature Markers
COL_SWITCH       = '#D946EF'       # Hap-mer switch error block (magenta)
COL_CUR_GAP      = '#E63946'       # Curation join (red)
COL_ASM_GAP      = '#111111'       # Assembly gap (black)

# Synteny Ribbons (Matching Supplementary Standard)
COL_RIB_COLL     = '#38BDF8'       # Sky blue
COL_RIB_COLL_E   = '#0284C7'       # Defined border for collinear ribbons
COL_RIB_NC       = '#F59E0B'       # Amber yellow for non-collinear ribbons (#F59E0B)
COL_RIB_NC_E     = '#B45309'       # Defining border for amber ribbons
COL_RIB_REC      = '#E11D48'       # Red for recovered query-1x (uncovered synteny)
COL_RIB_REC_E    = '#9F1239'
ALPHA_RIB_REC    = 0.25            # 25% opacity for uncovered synteny

# EF panel constants
_EF_ASM_ORDER  = ['Single', 'Dual', 'ONT_Dual']
_EF_CAT_ORDER  = ['macro', 'micro', 'dot']
_EF_CAT_LABELS = {'macro': 'Macro', 'micro': 'Micro', 'dot': 'Dot'}
_EF_ASM_COLORS = {
    'Single':   COL_COV_S_DARK,
    'Dual':     COL_COV_D_DARK,
    'ONT_Dual': COL_COV_O_DARK,
}

# Chromosome subsets for Main Figure
MACRO_TOKENS = ['3', '4', 'W']
MICRO_TOKENS = ['11', '14', '16']
DOT_TOKENS   = ['25', '30', '35', '37']

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
                parts = line.split()
                t_name, q_name = parts[2], parts[7]
                t_size, q_size = int(parts[3]), int(parts[8])
                sizes[t_name] = t_size
                sizes[q_name] = q_size
                q_strand = parts[9]
                t_pos = int(parts[5])
                q_raw_start = int(parts[10])
                q_raw_end   = int(parts[11])
                q_pos = (q_size - q_raw_end) if q_strand == '-' else q_raw_start
            else:
                parts = line.split()
                size = int(parts[0])
                dt = int(parts[1]) if len(parts) > 1 else 0
                dq = int(parts[2]) if len(parts) > 2 else 0

                rec_ribbons[t_name].append((t_name, t_pos, t_pos + size, q_name, q_pos, q_pos + size, q_strand, 'rec'))
                rec_spans[q_name].append((q_pos, q_pos + size))
                t_pos += size + dt
                q_pos += size + dq

    return rec_ribbons, rec_spans, sizes

def load_coverage_summary(path_s=None, path_d=None, path_ont=None):
    def _load(p):
        if not p or not os.path.exists(p):
            return {}
        df = pd.read_csv(p, sep='\t')
        df.columns = [c.strip() for c in df.columns]
        cov = {}
        for _, r in df.iterrows():
            nm = str(r['name'])
            c_pct = float(r.get('collinear_pct', 0.0))
            nc_pct = float(r.get('noncollinear_pct', 0.0))
            u_pct = float(r.get('uncovered_pct', max(0.0, 100.0 - c_pct - nc_pct)))
            cov[nm] = (c_pct, nc_pct, u_pct)
        return cov
    return _load(path_s), _load(path_d), _load(path_ont)

def chrom_token(name):
    m = re.search(r'chromosome_([0-9A-Za-z]+)', name)
    return m.group(1) if m else name

def chrom_sort_key(token):
    if token.isdigit():
        return (0, int(token), '')
    m = re.match(r'(\d+)([A-Za-z]+)', token)
    if m:
        return (0, int(m.group(1)), m.group(2))
    return (1, 0, token)

def is_pat(name):
    return 'Pat' in name

def classify_chrom_group(token):
    m = re.match(r'^(\d+)', token)
    if m:
        n = int(m.group(1))
        if n in (1, 2, 3, 4, 5):
            return 'macro'
        if 6 <= n <= 20:
            return 'micro'
        return 'dot'
    if token in ('Z', 'W'):
        return 'macro'
    return 'dot'

# =============================================================================
# Cytogenetic Drawing Utilities
# =============================================================================

def _rect_path(x0, x1, yc, h):
    y0, y1 = yc - h / 2, yc + h / 2
    verts = [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    return Path(verts, codes)

def _constriction_path(x0, x1, yc, h, c_start, c_end):
    y0 = yc - h / 2
    y1 = yc + h / 2
    w_min = max(0.015 * (x1 - x0), 180_000)
    cen_mid = (c_start + c_end) / 2
    cs = max(x0 + 50_000, cen_mid - w_min / 2)
    ce = min(x1 - 50_000, cen_mid + w_min / 2)

    if cs >= ce or cs <= x0 or ce >= x1:
        return _rect_path(x0, x1, yc, h)

    pinch_d = 0.22 * h
    y_bot_pinch = yc - pinch_d
    y_top_pinch = yc + pinch_d

    verts = [
        (x0, y0),
        (cs, y0),
        (cen_mid, y_bot_pinch),
        (ce, y0),
        (x1, y0),
        (x1, y1),
        (ce, y1),
        (cen_mid, y_top_pinch),
        (cs, y1),
        (x0, y1),
        (x0, y0)
    ]
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY
    ]
    return Path(verts, codes)

def draw_synteny_ribbon(ax, x0_t, x1_t, y_t, x0_q, x1_q, y_q, col, edge_col, alpha=0.45, inverted=False, lw=0.35):
    if abs(x1_t - x0_t) < 1: return
    if abs(x1_q - x0_q) < 1: return

    dy = (y_t - y_q) * 0.45
    c1_left  = (x0_t, y_t - dy)
    c2_left  = (x1_q if inverted else x0_q, y_q + dy)
    p_left_q = (x1_q if inverted else x0_q, y_q)

    c1_right  = (x0_q if inverted else x1_q, y_q + dy)
    c2_right  = (x1_t, y_t - dy)
    p_right_t = (x1_t, y_t)

    verts = [
        (x0_t, y_t),
        c1_left,
        c2_left,
        p_left_q,
        (x0_q if inverted else x1_q, y_q),
        c1_right,
        c2_right,
        p_right_t,
        (x0_t, y_t)
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.LINETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.CLOSEPOLY
    ]
    poly = PathPatch(Path(verts, codes), facecolor=col, edgecolor=edge_col,
                     linewidth=lw, alpha=alpha, zorder=3)
    ax.add_patch(poly)

def _draw_telo_semi(ax, x, arm, yc, bar_h, rx, fc='#000000', ec='#000000', lw=0.6, hollow=False):
    theta = np.linspace(np.pi / 2, 3 * np.pi / 2, 40) if arm == 'p' else np.linspace(-np.pi / 2, np.pi / 2, 40)
    ry = bar_h / 2
    verts = [(x + rx * np.cos(t), yc + ry * np.sin(t)) for t in theta]
    verts.append((x, yc - ry))
    verts.append((x, yc + ry))
    p = mpatches.Polygon(verts, closed=True,
                         fc='none' if hollow else fc,
                         ec=ec, lw=lw, zorder=12)
    ax.add_patch(p)

# =============================================================================
# Synteny Ideogram Panel Drawer for Macro / Micro / Dot
# =============================================================================

def draw_synteny_panel(ax_chr, ax_mat, ax_pat, tokens, max_global_size,
                       t2t_sizes, centromeres, t2t_telo, flip_set,
                       single_data, dual_data, ont_data,
                       tick_step=10_000_000, tick_unit_mb=True):
    """
    Renders butterfly synteny ideograms extending fully across the canvas with increased font sizes.
    """
    H_BAR_T2T = 0.16   # Bar thickness
    H_BAR_ASM = 0.16   # Bar thickness
    H_RIBBON  = 0.38   # Expanded ribbon span for clear visibility
    BLOCK_GAP = 0.08
    BLOCK_H   = H_BAR_T2T + H_RIBBON + H_BAR_ASM + BLOCK_GAP   # 0.78
    ROW_GAP   = 0.28
    ROW_H     = 3 * BLOCK_H + ROW_GAP                           # 2.62
    total_h   = len(tokens) * ROW_H + 0.35

    for ax in (ax_chr, ax_mat, ax_pat):
        ax.set_ylim(-0.35, total_h)
        ax.set_yticks([])
        ax.set_xticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    # Tight margin so ideograms span horizontally
    axis_max = int(np.ceil(max_global_size / tick_step)) * tick_step
    eff_max = max(max_global_size, axis_max)
    x_margin = eff_max * 0.008

    ax_chr.set_xlim(0, 1)
    ax_mat.set_xlim(eff_max + x_margin, -x_margin)
    ax_pat.set_xlim(-x_margin, eff_max + x_margin)

    telo_rx = eff_max * 0.0055

    for row_idx, tok in enumerate(tokens):
        y_row_top = total_h - 0.22 - row_idx * ROW_H

        # Zebra background shading
        if row_idx % 2 == 0:
            y_bg_top = y_row_top + 0.06
            y_bg_bot = y_row_top - (3 * BLOCK_H) - 0.05
            for ax in (ax_chr, ax_mat, ax_pat):
                ax.axhspan(y_bg_bot, y_bg_top, color='#F4F8FA', zorder=0)

        # Row divider line
        if row_idx > 0:
            y_div = y_row_top + ROW_GAP / 2
            ax_mat.axhline(y_div, color='#E2E8F0', lw=0.5, ls=':')
            ax_pat.axhline(y_div, color='#E2E8F0', lw=0.5, ls=':')

        # Chromosome label in Col 0 (ax_chr) - increased font 2 pt
        lbl_text = f"chr{tok}"
        y_mid_row = y_row_top - 1.5 * BLOCK_H + BLOCK_GAP / 2
        ax_chr.text(0.60, y_mid_row, lbl_text, ha='center', va='center',
                    fontsize=10.5, fontweight='bold', color='#0F172A')

        assemblies = [
            ('HiFi Single', single_data, COL_COV_S_DARK, COL_RIB_NC),
            ('HiFi Dual',   dual_data,   COL_COV_D_DARK, COL_RIB_NC),
            ('ONT Dual',    ont_data,    COL_COV_O_DARK, COL_RIB_NC)
        ]

        for b_idx, (asm_title, asm, c_dark, c_light) in enumerate(assemblies):
            y_block_top = y_row_top - b_idx * BLOCK_H
            y_t2t = y_block_top - H_BAR_T2T / 2
            y_asm = y_t2t - (H_BAR_T2T / 2 + H_RIBBON + H_BAR_ASM / 2)

            for side, ax in [('mat', ax_mat), ('pat', ax_pat)]:
                # Resolve T2T chromosome
                if tok == 'W':
                    if side != 'mat':
                        continue
                    t2t_chrom = 'Mat_NC_133064.1_chromosome_W'
                elif tok == 'Z':
                    if side != 'pat':
                        continue
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

                # Unlinked scaffolds (Matching supplementary standard logic)
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
                    has_nc = any(r[3] == u for r in asm['nc_ribbons'].get(t2t_chrom, []))
                    has_other = any(r[3] == u for r in asm['col_ribbons'].get(t2t_chrom, []) + asm.get('rec_ribbons', {}).get(t2t_chrom, []))
                    if has_nc:
                        tier = 0
                    elif has_other:
                        tier = 1
                    else:
                        tier = 2
                    return (tier, u_num, u)

                unloc_qual = sorted(list(unloc_candidates), key=_unloc_sort_key)

                # Compute X-offsets for unlinked scaffolds
                SPACER = max(eff_max * 0.012, 150_000)
                offsets = {q_prim: 0}
                cur_x = q_prim_len
                for u in unloc_qual:
                    cur_x += SPACER
                    offsets[u] = cur_x
                    cur_x += asm['sizes'].get(u, 500_000)

                # =============================================================
                # A. Top Bar: T2T Reference Ideogram (Pure White, Black Tips)
                # =============================================================
                tok_cen = 'W' if (tok == 'W' and side == 'mat') else ('Z' if (tok == 'Z' and side == 'pat') else tok)
                cen_span = centromeres.get((tok_cen.upper(), side.lower()))

                if cen_span:
                    cs = (t2t_len - cen_span[1]) if flip else cen_span[0]
                    ce = (t2t_len - cen_span[0]) if flip else cen_span[1]
                    t2t_path = _constriction_path(0, t2t_len, y_t2t, H_BAR_T2T, min(cs, ce), max(cs, ce))
                else:
                    t2t_path = _rect_path(0, t2t_len, y_t2t, H_BAR_T2T)

                ax.add_patch(PathPatch(t2t_path, fc=COL_T2T_FILL, ec=COL_BORDER, lw=0.55, zorder=5))

                # Liftover: Non-collinear coverage bands on T2T ideogram
                t2t_nc = _merge_intervals([(r[1], r[2]) for r in asm['nc_ribbons'].get(t2t_chrom, [])])
                for s, e in t2t_nc:
                    x0 = (t2t_len - e) if flip else s
                    x1 = (t2t_len - s) if flip else e
                    r = mpatches.Rectangle((min(x0, x1), y_t2t - H_BAR_T2T / 2),
                                           abs(x1 - x0), H_BAR_T2T,
                                           fc=c_dark, ec='none', alpha=0.25, zorder=6)
                    r.set_clip_path(t2t_path, transform=ax.transData)
                    ax.add_patch(r)

                # Liftover: Collinear coverage bands on T2T ideogram
                t2t_col = _merge_intervals([(r[1], r[2]) for r in asm['col_ribbons'].get(t2t_chrom, [])])
                for s, e in t2t_col:
                    x0 = (t2t_len - e) if flip else s
                    x1 = (t2t_len - s) if flip else e
                    r = mpatches.Rectangle((min(x0, x1), y_t2t - H_BAR_T2T / 2),
                                           abs(x1 - x0), H_BAR_T2T,
                                           fc=c_dark, ec='none', alpha=1.00, zorder=7)
                    r.set_clip_path(t2t_path, transform=ax.transData)
                    ax.add_patch(r)

                # T2T Telomeres (Tips filled with black)
                t2t_arms = t2t_telo.get(t2t_chrom, set())
                if p_arm in t2t_arms or not t2t_arms:
                    _draw_telo_semi(ax, 0, 'p', y_t2t, H_BAR_T2T, telo_rx, fc='#000000', ec='#000000')
                if q_arm in t2t_arms or not t2t_arms:
                    _draw_telo_semi(ax, t2t_len, 'q', y_t2t, H_BAR_T2T, telo_rx, fc='#000000', ec='#000000')

                # =============================================================
                # B. Bottom Bar: Assembly Ideogram (Primary Scaffold)
                # =============================================================
                if q_prim:
                    prim_path = _rect_path(0, q_prim_len, y_asm, H_BAR_ASM)
                    ax.add_patch(PathPatch(prim_path, fc=COL_UNALIGNED, ec=COL_BORDER, lw=0.55, zorder=5))

                    # Non-collinear coverage bands (under collinear)
                    m_nc = _merge_intervals(asm['nc_spans'].get(q_prim, []))
                    for s, e in m_nc:
                        x0 = (q_prim_len - e) if flip else s
                        x1 = (q_prim_len - s) if flip else e
                        r = mpatches.Rectangle((min(x0, x1), y_asm - H_BAR_ASM / 2),
                                               abs(x1 - x0), H_BAR_ASM,
                                               fc=c_dark, ec='none', alpha=0.25, zorder=6)
                        r.set_clip_path(prim_path, transform=ax.transData)
                        ax.add_patch(r)

                    # Collinear coverage bands
                    m_col = _merge_intervals(asm['col_spans'].get(q_prim, []))
                    for s, e in m_col:
                        x0 = (q_prim_len - e) if flip else s
                        x1 = (q_prim_len - s) if flip else e
                        r = mpatches.Rectangle((min(x0, x1), y_asm - H_BAR_ASM / 2),
                                               abs(x1 - x0), H_BAR_ASM,
                                               fc=c_dark, ec='none', alpha=1.00, zorder=7)
                        r.set_clip_path(prim_path, transform=ax.transData)
                        ax.add_patch(r)

                    # Switch error blocks (Magenta #D946EF)
                    for sw_s, sw_e in asm['sw'].get(q_prim, []):
                        x0 = (q_prim_len - sw_e) if flip else sw_s
                        x1 = (q_prim_len - sw_s) if flip else sw_e
                        w_sw = max(abs(x1 - x0), 80_000)
                        r = mpatches.Rectangle((min(x0, x1), y_asm - H_BAR_ASM / 2),
                                               w_sw, H_BAR_ASM,
                                               fc=COL_SWITCH, ec='none', zorder=8)
                        r.set_clip_path(prim_path, transform=ax.transData)
                        ax.add_patch(r)

                    # Gaps & Curation Joins (lw=0.75)
                    for gs, ge, gtype in asm['gaps'].get(q_prim, []):
                        gx = (q_prim_len - (gs + ge) / 2) if flip else ((gs + ge) / 2)
                        gap_col = COL_CUR_GAP if 'CURATION' in gtype else COL_ASM_GAP
                        ax.plot([gx, gx], [y_asm - H_BAR_ASM / 2, y_asm + H_BAR_ASM / 2],
                                color=gap_col, lw=0.75, zorder=9, solid_capstyle='butt')

                    # Assembly Telomeres
                    asm_coll_arms = asm['telo_coll'].get(t2t_chrom, set()) | asm['telo_coll'].get(q_prim, set())
                    asm_nc_arms   = asm['telo_nc'].get(t2t_chrom, set()) | asm['telo_nc'].get(q_prim, set())

                    if p_arm in asm_coll_arms:
                        _draw_telo_semi(ax, 0, 'p', y_asm, H_BAR_ASM, telo_rx, fc='#000000', ec='#000000', hollow=False)
                    elif p_arm in asm_nc_arms:
                        _draw_telo_semi(ax, 0, 'p', y_asm, H_BAR_ASM, telo_rx, fc='white', ec='#000000', lw=0.6, hollow=True)

                    if q_arm in asm_coll_arms:
                        _draw_telo_semi(ax, q_prim_len, 'q', y_asm, H_BAR_ASM, telo_rx, fc='#000000', ec='#000000', hollow=False)
                    elif q_arm in asm_nc_arms:
                        _draw_telo_semi(ax, q_prim_len, 'q', y_asm, H_BAR_ASM, telo_rx, fc='white', ec='#000000', lw=0.6, hollow=True)

                # =============================================================
                # C. Unlinked Scaffolds at Chromosome End
                # =============================================================
                for u in unloc_qual:
                    u_off = offsets[u]
                    u_len = asm['sizes'].get(u, 500_000)

                    u_path = _rect_path(u_off, u_off + u_len, y_asm, H_BAR_ASM)
                    ax.add_patch(PathPatch(u_path, fc=COL_UNALIGNED, ec='#64748B', lw=0.5, ls='--', zorder=5))

                    u_nc = _merge_intervals(asm['nc_spans'].get(u, []))
                    for s, e in u_nc:
                        r = mpatches.Rectangle((u_off + s, y_asm - H_BAR_ASM / 2), e - s, H_BAR_ASM,
                                               fc=c_dark, ec='none', alpha=0.25, zorder=6)
                        r.set_clip_path(u_path, transform=ax.transData)
                        ax.add_patch(r)

                    u_col = _merge_intervals(asm['col_spans'].get(u, []))
                    for s, e in u_col:
                        r = mpatches.Rectangle((u_off + s, y_asm - H_BAR_ASM / 2), e - s, H_BAR_ASM,
                                               fc=c_dark, ec='none', alpha=1.00, zorder=7)
                        r.set_clip_path(u_path, transform=ax.transData)
                        ax.add_patch(r)

                    for gs, ge, gtype in asm['gaps'].get(u, []):
                        gx = u_off + (gs + ge) / 2
                        gap_col = COL_CUR_GAP if 'CURATION' in gtype else COL_ASM_GAP
                        ax.plot([gx, gx], [y_asm - H_BAR_ASM / 2, y_asm + H_BAR_ASM / 2],
                                color=gap_col, lw=0.75, zorder=9, solid_capstyle='butt')

                    m_lbl = re.search(r'unloc_([0-9]+)', u)
                    s_tag = f"u{m_lbl.group(1)}" if m_lbl else "unloc"
                    ax.text(u_off + u_len / 2, y_asm + H_BAR_ASM * 0.70, s_tag,
                            ha='left', va='bottom', fontsize=7.0, color='#475569',
                            fontweight='bold', rotation=40)

                # =============================================================
                # D. Synteny Ribbons (Exact Supplementary Standard)
                # =============================================================
                y_rib_t = y_t2t - H_BAR_T2T / 2
                y_rib_q = y_asm + H_BAR_ASM / 2

                # 1. Collinear Ribbons (Assembly color alpha=0.35)
                for r in asm['col_ribbons'].get(t2t_chrom, []):
                    _, ts, te, qn, qs, qe, qstr, _ = r
                    if qn not in offsets: continue
                    q_off = offsets[qn]

                    if flip:
                        t0, t1 = t2t_len - te, t2t_len - ts
                        if qn == q_prim:
                            q0, q1 = q_prim_len - qe, q_prim_len - qs
                        else:
                            q0, q1 = q_off + qs, q_off + qe
                    else:
                        t0, t1 = ts, te
                        q0, q1 = q_off + qs, q_off + qe

                    draw_synteny_ribbon(ax, min(t0, t1), max(t0, t1), y_rib_t,
                                        min(q0, q1), max(q0, q1), y_rib_q,
                                        c_dark, c_dark, alpha=0.35,
                                        inverted=(qstr == '-'), lw=0.35)

                # 2. Non-collinear Ribbons (Orange #F59E0B alpha=0.55)
                for r in asm['nc_ribbons'].get(t2t_chrom, []):
                    _, ts, te, qn, qs, qe, qstr, _ = r
                    if qn not in offsets: continue
                    if max(te - ts, abs(qe - qs)) < MIN_NC_RIBBON_BP:
                        continue
                    q_off = offsets[qn]

                    if flip:
                        t0, t1 = t2t_len - te, t2t_len - ts
                        if qn == q_prim:
                            q0, q1 = q_prim_len - qe, q_prim_len - qs
                        else:
                            q0, q1 = q_off + qs, q_off + qe
                    else:
                        t0, t1 = ts, te
                        q0, q1 = q_off + qs, q_off + qe

                    draw_synteny_ribbon(ax, min(t0, t1), max(t0, t1), y_rib_t,
                                        min(q0, q1), max(q0, q1), y_rib_q,
                                        COL_RIB_NC, COL_RIB_NC_E, alpha=0.55,
                                        inverted=(qstr == '-'), lw=0.35)

                # 3. Recovered Uncovered Synteny Ribbons (Red #E11D48)
                for r in asm.get('rec_ribbons', {}).get(t2t_chrom, []):
                    _, ts, te, qn, qs, qe, qstr, _ = r
                    if qn not in offsets: continue
                    if max(te - ts, abs(qe - qs)) < MIN_NC_RIBBON_BP:
                        continue
                    q_off = offsets[qn]

                    if flip:
                        t0, t1 = t2t_len - te, t2t_len - ts
                        if qn == q_prim:
                            q0, q1 = q_prim_len - qe, q_prim_len - qs
                        else:
                            q0, q1 = q_off + qs, q_off + qe
                    else:
                        t0, t1 = ts, te
                        q0, q1 = q_off + qs, q_off + qe

                    draw_synteny_ribbon(ax, min(t0, t1), max(t0, t1), y_rib_t,
                                        min(q0, q1), max(q0, q1), y_rib_q,
                                        COL_RIB_REC, COL_RIB_REC_E, alpha=ALPHA_RIB_REC,
                                        inverted=(qstr == '-'))

    # Bottom coordinate axis line & scale ticks (Increased font 2 pt)
    y_axis = 0.08
    tick_h = 0.08

    for ax, side in [(ax_mat, 'mat'), (ax_pat, 'pat')]:
        ax.plot([0, axis_max], [y_axis, y_axis], color='#1E293B', lw=0.8, zorder=10)
        for tick_val in range(0, axis_max + 1, tick_step):
            num_val = (tick_val // 1_000_000) if tick_unit_mb else tick_val
            ax.plot([tick_val, tick_val], [y_axis, y_axis + tick_h], color='#1E293B', lw=0.8, zorder=10)
            ax.text(tick_val, y_axis - 0.08, f'{num_val}', ha='center', va='top',
                    fontsize=8.8, color='#1E293B', fontweight='medium')

        lbl_x = axis_max + tick_step * 0.15
        ha = 'right' if side == 'mat' else 'left'
        ax.text(lbl_x, y_axis, 'Mb' if tick_unit_mb else 'bp', ha=ha, va='center',
                fontsize=9.5, fontweight='bold', color='#1E293B')

# =============================================================================
# Panel a / b / c / g Drawers
# =============================================================================

def compute_avg_coverage_by_method(group_order, chrom_groups,
                                   cov_summary_s, cov_summary_d, cov_summary_ont=None):
    def _calc(cov_summary):
        if not cov_summary:
            return {}
        avg = {}
        for tok in group_order:
            members = chrom_groups.get(tok, [])
            if not members:
                continue
            c_sum, nc_sum, u_sum, n = 0.0, 0.0, 0.0, 0
            for nm in members:
                if nm in cov_summary:
                    c, nc, u = cov_summary[nm]
                    c_sum += c; nc_sum += nc; u_sum += u; n += 1
            if n > 0:
                avg[tok] = (c_sum / n, nc_sum / n, u_sum / n)
        return avg
    return _calc(cov_summary_s), _calc(cov_summary_d), _calc(cov_summary_ont)

def _draw_full_coverage_panel(ax, all_macro_order, all_micro_order, all_nano_order,
                               chrom_groups, cov_summary_s, cov_summary_d,
                               cov_summary_ont=None):
    all_order = all_macro_order + all_micro_order + all_nano_order
    avg_s, avg_d, avg_ont = compute_avg_coverage_by_method(
        all_order, chrom_groups, cov_summary_s, cov_summary_d, cov_summary_ont)

    bar_w     = 0.60
    pair_gap  = 0.04
    chrom_gap = 0.50
    group_gap = 1.40

    group_specs = [
        ('Macrochromosomes', all_macro_order),
        ('Microchromosomes', all_micro_order),
        ('Dot chromosomes',  all_nano_order),
    ]

    x_pos_s, x_pos_d, x_pos_o, x_tick = {}, {}, {}, {}
    group_xranges = {}
    x = 0.0
    first = True
    for gname, tokens in group_specs:
        present = [t for t in tokens if t in avg_s or t in avg_d or t in avg_ont]
        if not present:
            continue
        if not first:
            x += group_gap
        first = False
        gx_start = x
        for tok in present:
            x_pos_s[tok] = x
            x_pos_d[tok] = x + bar_w + pair_gap
            x_pos_o[tok] = x + 2 * (bar_w + pair_gap)
            x_tick[tok]  = x + (3 * bar_w + 2 * pair_gap) / 2
            x += 3 * bar_w + 2 * pair_gap + chrom_gap
        group_xranges[gname] = (gx_start, x - chrom_gap)

    total_x = max(x - chrom_gap, 1.0)
    ax.set_xlim(-0.6, total_x + 0.6)
    ax.set_ylim(0, 110)
    ax.set_facecolor('white')

    for i, tok in enumerate(x_tick):
        if i % 2 == 0:
            x0 = x_pos_s[tok] - chrom_gap * 0.35
            x1 = x_pos_o.get(tok, x_pos_d[tok]) + bar_w + chrom_gap * 0.35
            ax.axvspan(x0, x1, color=COL_ROW_ODD, zorder=0)

    # In coverage subplot: Collinear uses assembly color, Non-collinear uses unified amber yellow COL_RIB_NC
    bar_specs = [
        (x_pos_s, avg_s,   COL_COV_S_DARK, COL_RIB_NC),
        (x_pos_d, avg_d,   COL_COV_D_DARK, COL_RIB_NC),
        (x_pos_o, avg_ont, COL_COV_O_DARK, COL_RIB_NC),
    ]
    for tok in x_tick:
        for xc_map, avg, col_dark, col_light in bar_specs:
            xc = xc_map.get(tok)
            if xc is None or not avg or tok not in avg:
                continue
            cp, ncp, up = avg[tok]
            ax.bar(xc, cp,  width=bar_w, bottom=0, align='edge',
                   color=col_dark, ec='none', zorder=2)
            ax.bar(xc, ncp, width=bar_w, bottom=cp, align='edge',
                   color=col_dark, alpha=0.25, ec='none', zorder=2)
            ax.bar(xc, up,  width=bar_w, bottom=cp+ncp, align='edge',
                   color=COL_UNCOV, ec='none', zorder=2)
            ax.add_patch(mpatches.Rectangle((xc, 0), bar_w, 100,
                         fc='none', ec=COL_BORDER, lw=0.3, zorder=4))

    ax.set_xticks(list(x_tick.values()))
    ax.set_xticklabels(list(x_tick.keys()), fontsize=7.2, rotation=0, ha='center')
    ax.tick_params(axis='x', length=0, pad=1.5)

    for gname, (gx0, gx1) in group_xranges.items():
        x_data_mid = (gx0 + gx1) / 2
        x_frac = (x_data_mid - ax.get_xlim()[0]) / (ax.get_xlim()[1] - ax.get_xlim()[0])
        ax.text(x_frac, 1.04, gname, ha='center', va='bottom',
                fontsize=9.0, fontweight='bold', color='#333333',
                transform=ax.transAxes, clip_on=False)
        if gx0 > 0:
            ax.axvline(gx0 - group_gap / 2, color='#cccccc', lw=0.6, ls='--', zorder=1)

    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(['0', '25', '50', '75', '100%'], fontsize=8.0)
    ax.set_ylabel('Coverage (%)', fontsize=9.2, labelpad=2)
    ax.tick_params(axis='y', labelsize=8.0, length=2.5)

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_linewidth(0.4)
    ax.spines['bottom'].set_linewidth(0.4)

def _draw_panel_d_violin(axes_by_cat, all_macro_order, all_micro_order, all_nano_order,
                          chrom_groups, cov_summary_s, cov_summary_d,
                          cov_summary_ont=None):
    group_specs = [
        ('Macro', all_macro_order),
        ('Micro', all_micro_order),
        ('Dot',   all_nano_order),
    ]
    method_info = [
        ('Single',   cov_summary_s,   COL_COV_S_DARK),
        ('Dual',     cov_summary_d,   COL_COV_D_DARK),
        ('ONT_Dual', cov_summary_ont, COL_COV_O_DARK),
    ]

    vdata = {}
    for gname, tokens in group_specs:
        for mname, cov_dict, _ in method_info:
            mat_vals, pat_vals = [], []
            for tok in tokens:
                for nm in chrom_groups.get(tok, []):
                    if cov_dict and nm in cov_dict:
                        coll, _, _ = cov_dict[nm]
                        val = coll
                    else:
                        val = 0.0
                    (pat_vals if is_pat(nm) else mat_vals).append((val, tok))
            vdata[(gname, mname)] = {'mat': mat_vals, 'pat': pat_vals}

    half_w   = 0.14
    meth_gap = 0.06
    meth_centers = {}
    x = 0.0
    for mi, (mname, _, _) in enumerate(method_info):
        if mi > 0:
            x += meth_gap
        meth_centers[mname] = x + half_w
        x += half_w * 2
    x_max = x

    np.random.seed(42)

    for gi, (gname, _) in enumerate(group_specs):
        ax = axes_by_cat.get(gname)
        if ax is None:
            continue

        all_vals = [v for mname, _, _ in method_info
                    for side in ('mat', 'pat')
                    for v, _ in vdata.get((gname, mname), {}).get(side, [])]
        if all_vals:
            dmin, dmax = min(all_vals), max(all_vals)
            pad  = max((dmax - dmin) * 0.10, 1.5)
            y_lo = max(0.0,   dmin - pad)
            y_hi = min(100.0, dmax + pad)
            if y_hi - y_lo < 4:
                mid  = (y_lo + y_hi) / 2
                y_lo = max(0.0,   mid - 2.0)
                y_hi = min(100.0, mid + 2.0)
        else:
            y_lo, y_hi = 0.0, 100.0

        ax.set_xlim(-half_w - 0.08, x_max + half_w + 0.08)
        ax.set_ylim(y_lo, y_hi + (y_hi - y_lo) * 0.06)
        ax.set_facecolor('white')

        def _hv(vals_toks, xc, side, color, _ax=ax, _y_lo=y_lo, _y_hi=y_hi):
            if not vals_toks:
                return []
            vals = np.array([v for v, _ in vals_toks], dtype=float)
            toks = [t for _, t in vals_toks]
            bx0 = xc - half_w if side == 'left' else xc
            bx1 = xc          if side == 'left' else xc + half_w
            wx  = (bx0 + bx1) / 2

            kde_lo = max(_y_lo - 1.0, float(vals.min()) - 1.0)
            kde_hi = min(_y_hi + 1.0, float(vals.max()) + 1.0)
            if kde_hi <= kde_lo:
                kde_hi = kde_lo + 0.5
            y_grid = np.linspace(kde_lo, kde_hi, 200)
            try:
                kde_fn  = gaussian_kde(vals, bw_method='scott')
                density = kde_fn(y_grid)
                scale   = half_w * 0.88 / density.max() if density.max() > 0 else 0
                density = density * scale
            except Exception:
                density = np.zeros(len(y_grid))

            if side == 'left':
                outer_x = xc - density
                poly_x  = np.concatenate([outer_x, [xc, xc]])
                poly_y  = np.concatenate([y_grid,  [y_grid[-1], y_grid[0]]])
            else:
                outer_x = xc + density
                poly_x  = np.concatenate([[xc, xc], outer_x[::-1]])
                poly_y  = np.concatenate([[y_grid[0], y_grid[-1]], y_grid[::-1]])

            _ax.add_patch(mpatches.Polygon(
                list(zip(poly_x, poly_y)), closed=True,
                fc=color, ec=color, lw=0.35, alpha=0.28, zorder=3))
            _ax.plot(outer_x, y_grid, color=color, lw=0.5, alpha=0.75, zorder=4)

            q2 = float(np.median(vals))
            _ax.plot([bx0, bx1], [q2, q2], color=COL_BORDER, lw=0.8, zorder=5)

            jitter = np.random.uniform(-(half_w * 0.38), half_w * 0.38, len(vals))
            xpts = np.clip(wx + jitter, bx0 + 0.005, bx1 - 0.005)
            _ax.scatter(xpts, vals, color=color, s=2.5, zorder=5, ec='none', alpha=0.65)
            return list(zip(xpts, vals, toks, [side] * len(vals)))

        for mname, _, color in method_info:
            xc = meth_centers[mname]
            sides = vdata[(gname, mname)]
            pts_mat = _hv(sides['mat'], xc, 'left',  color)
            pts_pat = _hv(sides['pat'], xc, 'right', color)

            # Label up to 2 lowest outliers below median per assembly
            all_pts = pts_mat + pts_pat
            if all_pts:
                all_vals = [p[1] for p in all_pts]
                med = np.median(all_vals)
                outliers = [p for p in all_pts if p[1] < med - 1e-4]
                outliers.sort(key=lambda p: p[1])
                for xpt, val, tok, s_name in outliers[:2]:
                    ha = 'right' if s_name == 'left' else 'left'
                    dx = -0.012 if s_name == 'left' else 0.012
                    ax.text(xpt + dx, val, f'{tok}', fontsize=6.5, color='#0F172A',
                            fontweight='bold', ha=ha, va='center', zorder=10)

        ax.set_xticks([x_max / 2])
        ax.set_xticklabels([gname], fontsize=8.5, fontweight='medium')
        ax.tick_params(axis='x', length=0, pad=1.5)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=4, integer=True))
        ax.tick_params(axis='y', labelsize=8.0, length=2)
        if gi == 0:
            ax.set_ylabel('Collinear cov. (%)', fontsize=9.0, labelpad=2)
            ax.text(0.02, 1.01, 'Mat◀|▶Pat', transform=ax.transAxes,
                    fontsize=7.0, color='#555555', ha='left', va='bottom')

        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        ax.spines['left'].set_linewidth(0.4)
        ax.spines['bottom'].set_linewidth(0.4)

def _draw_panel_f(ax, df):
    gap = 0.80
    off = 0.18
    slw = 0.9
    ms  = 4.5
    lw  = 0.4

    cat_x   = {cat: i * gap for i, cat in enumerate(_EF_CAT_ORDER)}
    n_asm   = len(_EF_ASM_ORDER)
    offsets = np.linspace(-(n_asm - 1) / 2 * off, (n_asm - 1) / 2 * off, n_asm)

    cat_expected = {}
    for cat in _EF_CAT_ORDER:
        rows = df[df['category'] == cat]
        if not rows.empty:
            cat_expected[cat] = int(rows['n_telo_expected'].iloc[0])

    for ai, asm in enumerate(_EF_ASM_ORDER):
        color = _EF_ASM_COLORS[asm]
        sub   = df[(df['assembly'] == asm) & (df['category'].isin(_EF_CAT_ORDER))]
        sub   = sub.set_index('category')

        for cat in _EF_CAT_ORDER:
            if cat not in sub.index:
                continue
            x         = cat_x[cat] + offsets[ai]
            row       = sub.loc[cat]
            frac      = row['telomere_pct'] / 100.0
            n_present = int(row['n_telo_present'])

            ax.plot([x, x], [0, frac], color=color, lw=slw, solid_capstyle='round', zorder=3)
            ax.scatter(x, frac, color=color, s=ms ** 2 * 0.5, zorder=4, ec='white', linewidths=0.3)
            ax.text(x, frac + 0.05, f'{n_present}', ha='center', va='bottom',
                    fontsize=7.0, color=color, zorder=5)

    ax.axhline(1.0, color='#aaaaaa', lw=0.5, ls='--', zorder=1)
    ax.set_xlim(-gap * 0.6, (len(_EF_CAT_ORDER) - 1) * gap + gap * 0.6)
    ax.set_ylim(0, 1.45)
    ax.set_yticks([0, 0.50, 1.00])
    ax.set_yticklabels(['0', '50', '100%'], fontsize=8.0)
    ax.set_ylabel('Telomere (%)', fontsize=8.5, labelpad=2)

    xtick_pos = [cat_x[c] for c in _EF_CAT_ORDER]
    ax.set_xticks(xtick_pos)
    ax.set_xticklabels([f"{_EF_CAT_LABELS[c]}\n(n={cat_expected.get(c, '?')})" for c in _EF_CAT_ORDER],
                       fontsize=8.0)
    ax.tick_params(axis='x', length=0, pad=2)
    ax.tick_params(axis='y', labelsize=8.0, length=2)

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_linewidth(lw)
    ax.spines['bottom'].set_linewidth(lw)
    ax.set_facecolor('white')

def _load_annotation_tsv(path):
    if not path or not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep='\t')
    for col in ['Frameshifts_errors (%)', 'Premature stop codon(%)']:
        df[col] = df[col].astype(str).str.replace(',', '.').astype(float)
    df['technology'] = df['Assembly'].str.split('+').str[0]
    return df

def _draw_panel_g(ax_left, ax_right, annotation_df):
    """
    Draws Frameshifts (%) and Premature stop codons (%) bar charts spanning the full width
    with PCG counts placed directly at the top of each corresponding bar.
    """
    _G_ROW_SPECS = [('HiFi', 'hap1'), ('HiFi', 'hap2'), ('ONT', 'hap1'), ('ONT', 'hap2')]
    _GX = {('HiFi', 'hap1'): 0.00, ('HiFi', 'hap2'): 0.74, ('ONT', 'hap1'): 1.91, ('ONT', 'hap2'): 2.65}
    xlim_bars = (-0.46, 3.11)
    bar_w = 0.60
    fs_tick = 8.0
    fs_title = 9.0

    col_pcg = next((c for c in ['Nr.protein coding genes', 'PCGs', 'protein_coding_genes'] if c in annotation_df.columns), None)

    panel_specs = [
        (ax_left,  'Frameshifts_errors (%)',  'Frameshifts (%)',           2.0),
        (ax_right, 'Premature stop codon(%)', 'Premature stop codons (%)', 0.2),
    ]

    for ax, col, title, tick_step in panel_specs:
        vals = [
            float(annotation_df.loc[
                (annotation_df['technology'] == t) &
                (annotation_df['Haplotype']  == h), col
            ].iloc[0])
            for t, h in _G_ROW_SPECS
            if not annotation_df.loc[
                (annotation_df['technology'] == t) &
                (annotation_df['Haplotype']  == h)
            ].empty
        ]
        v_max    = max(vals) if vals else 1.0
        ylim_max = float(int(np.ceil(v_max * 1.35 / tick_step)) * tick_step)
        yticks   = list(np.round(np.arange(0, ylim_max + tick_step * 0.5, tick_step), 10))

        ax.set_xlim(*xlim_bars)
        ax.set_ylim(0, ylim_max)
        ax.set_yticks(yticks)
        if tick_step < 1:
            ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: f'{y:.1f}'))
        else:
            ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: f'{y:.0f}'))

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(0.4)
        ax.spines['bottom'].set_linewidth(0.4)
        ax.tick_params(axis='y', length=2.0, width=0.4, labelsize=fs_tick, pad=1.5)
        ax.set_facecolor('none')

        ax.set_title(title, fontsize=fs_title, pad=3, fontweight='bold', loc='center', color='#0F172A')
        ax.plot([0, 1], [1, 1], transform=ax.transAxes, color='#333333', lw=0.25, clip_on=False)

        # Plot bars and place PCG numbers at the top of the bars
        for tech, hap in _G_ROW_SPECS:
            row = annotation_df.loc[
                (annotation_df['technology'] == tech) &
                (annotation_df['Haplotype']  == hap)
            ]
            if row.empty:
                continue
            v = float(row[col].iloc[0])
            x = _GX[(tech, hap)]
            if tech == 'HiFi':
                ax.bar(x, v, width=bar_w, facecolor='black', edgecolor='black', lw=0.4, zorder=3)
            else:
                ax.bar(x, v, width=bar_w, facecolor='white', edgecolor='black', lw=0.5, hatch='////', zorder=3)

            # Add PCG count at the top of the bar
            if col_pcg:
                pcg_val = int(float(str(row[col_pcg].iloc[0]).replace(',', '')))
                ax.text(x, v + ylim_max * 0.03, f"{pcg_val:,}", ha='center', va='bottom',
                        fontsize=7.8, fontweight='bold', color='#1E293B', zorder=5)

        # X-tick labels under bars
        xtick_locs = [_GX[spec] for spec in _G_ROW_SPECS]
        xtick_lbls = [f"{tech} H{hap[-1]}" for tech, hap in _G_ROW_SPECS]
        ax.set_xticks(xtick_locs)
        ax.set_xticklabels(xtick_lbls, fontsize=7.8, rotation=0, ha='center', color='#1E293B')
        ax.tick_params(axis='x', length=0, pad=2)

# =============================================================================
# Neat & Tidy Legend Card (Beside Panel g in Bottom Row)
# =============================================================================

def _draw_legend_card(ax_leg):
    """
    Neat, refined publication legend card with compact symbols and increased font sizes.
    """
    ax_leg.set_xlim(0, 1)
    ax_leg.set_ylim(0, 1)
    ax_leg.axis('off')
    ax_leg.add_patch(mpatches.FancyBboxPatch((0.005, 0.02), 0.99, 0.96,
                                             boxstyle="round,pad=0.015,rounding_size=0.03",
                                             fc='#F8FAFC', ec='#CBD5E1', lw=0.5))

    fsize_hdr = 8.5
    fsize_txt = 7.8

    # ── Section 1: Assemblies (x: 0.02 to 0.23) ────────────────────────────
    ax_leg.text(0.02, 0.85, 'Assemblies (Coll/NC):', fontsize=fsize_hdr, fontweight='bold', color='#1E293B')
    asms = [
        ('HiFi Single', COL_COV_S_DARK),
        ('HiFi Dual',   COL_COV_D_DARK),
        ('ONT Dual',    COL_COV_O_DARK)
    ]
    for i, (name, cd) in enumerate(asms):
        by = 0.63 - i * 0.24
        ax_leg.add_patch(mpatches.Rectangle((0.020, by - 0.03), 0.015, 0.14, fc=cd, ec='none', alpha=1.00))
        ax_leg.add_patch(mpatches.Rectangle((0.035, by - 0.03), 0.015, 0.14, fc=cd, ec='none', alpha=0.25))
        ax_leg.text(0.055, by + 0.04, name, va='center', fontsize=fsize_txt, color='#334155')

    # ── Section 2: Synteny Ribbons & Coverage (x: 0.24 to 0.48) ────────────
    ax_leg.text(0.24, 0.85, 'Synteny Ribbons:', fontsize=fsize_hdr, fontweight='bold', color='#1E293B')
    by0 = 0.63
    ax_leg.add_patch(mpatches.Rectangle((0.240, by0 - 0.03), 0.008, 0.14, fc=COL_COV_S_DARK, ec='none', alpha=0.35))
    ax_leg.add_patch(mpatches.Rectangle((0.248, by0 - 0.03), 0.008, 0.14, fc=COL_COV_D_DARK, ec='none', alpha=0.35))
    ax_leg.add_patch(mpatches.Rectangle((0.256, by0 - 0.03), 0.008, 0.14, fc=COL_COV_O_DARK, ec='none', alpha=0.35))
    ax_leg.add_patch(mpatches.Rectangle((0.240, by0 - 0.03), 0.024, 0.14, fc='none', ec='#64748B', lw=0.4))
    ax_leg.text(0.270, by0 + 0.04, 'Collinear', va='center', fontsize=fsize_txt, color='#334155')

    by1 = 0.63 - 1 * 0.24
    ax_leg.add_patch(mpatches.Rectangle((0.240, by1 - 0.03), 0.024, 0.14, fc=COL_RIB_NC, ec=COL_RIB_NC_E, lw=0.4, alpha=0.55))
    ax_leg.text(0.270, by1 + 0.04, 'Non-collinear', va='center', fontsize=fsize_txt, color='#334155')

    by2 = 0.63 - 2 * 0.24
    ax_leg.add_patch(mpatches.Rectangle((0.240, by2 - 0.03), 0.024, 0.14, fc=COL_RIB_REC, ec=COL_RIB_REC_E, lw=0.4, alpha=ALPHA_RIB_REC))
    ax_leg.text(0.270, by2 + 0.04, 'Uncovered Synteny', va='center', fontsize=fsize_txt, color='#334155')

    # ── Section 3: Feature Markers (x: 0.49 to 0.83) ───────────────────────
    ax_leg.text(0.49, 0.85, 'Feature Markers:', fontsize=fsize_hdr, fontweight='bold', color='#1E293B')
    feats_col1 = [
        ('Switch Block',   'patch',   COL_SWITCH),
        ('Unloc Scaffold', 'unloc',   '#64748B'),
        ('Curation Join',  'line',    COL_CUR_GAP),
    ]
    feats_col2 = [
        ('Assembly Gap',   'line',    COL_ASM_GAP),
        ('Telomere (Coll)','telo_c',  '#000000'),
        ('Telomere (NC)',  'telo_nc', '#000000'),
    ]
    for i, (name, ftype, col) in enumerate(feats_col1):
        by = 0.63 - i * 0.24
        if ftype == 'patch':
            ax_leg.add_patch(mpatches.Rectangle((0.490, by - 0.03), 0.020, 0.14, fc=col, ec='none'))
        elif ftype == 'unloc':
            ax_leg.add_patch(mpatches.Rectangle((0.490, by - 0.03), 0.020, 0.14, fc='none', ec=col, ls='--', lw=0.5))
        elif ftype == 'line':
            ax_leg.plot([0.500, 0.500], [by - 0.03, by + 0.11], color=col, lw=1.2)
        ax_leg.text(0.518, by + 0.04, name, va='center', fontsize=fsize_txt, color='#334155')

    for i, (name, ftype, col) in enumerate(feats_col2):
        by = 0.63 - i * 0.24
        if ftype == 'line':
            ax_leg.plot([0.665, 0.665], [by - 0.03, by + 0.11], color=col, lw=1.2)
        elif ftype == 'telo_c':
            ax_leg.scatter([0.665], [by + 0.04], color=col, s=12, ec='none')
        elif ftype == 'telo_nc':
            ax_leg.scatter([0.665], [by + 0.04], color='white', s=12, ec=col, linewidths=0.7)
        ax_leg.text(0.686, by + 0.04, name, va='center', fontsize=fsize_txt, color='#334155')

    # ── Section 4: Annotation Tech (x: 0.84 to 0.98) ───────────────────────
    ax_leg.text(0.84, 0.85, 'Annotation Tech:', fontsize=fsize_hdr, fontweight='bold', color='#1E293B')
    ax_leg.add_patch(mpatches.Rectangle((0.84, 0.52), 0.026, 0.16, facecolor='black', edgecolor='black', lw=0.4))
    ax_leg.text(0.876, 0.60, 'HiFi', va='center', fontsize=fsize_txt, color='#334155')
    ax_leg.add_patch(mpatches.Rectangle((0.84, 0.20), 0.026, 0.16, facecolor='white', edgecolor='black', lw=0.5, hatch='////'))
    ax_leg.text(0.876, 0.28, 'ONT', va='center', fontsize=fsize_txt, color='#334155')

# =============================================================================
# Main Figure Builder
# =============================================================================

def build_combined_figure(t2t_sizes, centromeres, p_arm_dict, t2t_telo,
                          single_data, dual_data, ont_data,
                          cov_sum_s, cov_sum_d, cov_sum_o,
                          df_stats, annotation_df,
                          out_png, out_pdf=None, dpi=300):
    """
    Constructs the publication combined figure sized with extended ideograms (fig_h = 18.5 in):
      - Panels a, b, c (Row 0): Compact top row
      - Ideograms (Rows 1, 2, 3): Extended ideograms with tall synteny ribbons and T2T liftover
      - Panel g + Legend (Row 4): Compact bottom row with PCG numbers on bars
    """
    flip_set = build_flip_set_from_centromeres_and_telomeres(t2t_sizes, centromeres, p_arm_dict)

    # Classify chromosomes into groups
    chrom_groups = defaultdict(list)
    for chrom in t2t_sizes.keys():
        tok = chrom_token(chrom)
        chrom_groups[tok].append(chrom)

    tokens_sorted = sorted(chrom_groups.keys(), key=chrom_sort_key)
    all_macro = [t for t in tokens_sorted if classify_chrom_group(t) == 'macro']
    all_micro = [t for t in tokens_sorted if classify_chrom_group(t) == 'micro']
    all_dot   = [t for t in tokens_sorted if classify_chrom_group(t) == 'dot']

    # Max global size per panel accommodating scaffolds and unlocs
    def get_max_size(tok_list):
        max_s = 0
        for tok in tok_list:
            for side in ('mat', 'pat'):
                if tok == 'W' and side == 'mat':
                    t2t_c = 'Mat_NC_133064.1_chromosome_W'
                elif tok == 'Z' and side == 'pat':
                    t2t_c = 'Mat_NC_133063.1_chromosome_Z'
                else:
                    t2t_c = next((x for x in t2t_sizes.keys() if f'chromosome_{tok}' in x and (side == 'pat' if 'Pat' in x else side == 'mat')), None)
                if not t2t_c: continue
                max_s = max(max_s, t2t_sizes.get(t2t_c, 0))
                for asm_d in (single_data, dual_data, ont_data):
                    qp = asm_d['pairs'].get(t2t_c)
                    if not qp: continue
                    tot = asm_d['sizes'].get(qp, 0)
                    parts = qp.rsplit('.', 1)
                    stem = parts[0]
                    suf = f".{parts[1]}" if len(parts) > 1 else ""
                    u_set = set()
                    for r in asm_d['col_ribbons'].get(t2t_c, []) + asm_d['nc_ribbons'].get(t2t_c, []) + asm_d.get('rec_ribbons', {}).get(t2t_c, []):
                        if r[3] != qp: u_set.add(r[3])
                    for s in asm_d['sizes']:
                        if s.startswith(f"{stem}_unloc_") and (not suf or s.endswith(suf)):
                            u_set.add(s)
                    for u in u_set:
                        tot += max(max_s * 0.012, 150_000) + asm_d['sizes'].get(u, 500_000)
                    max_s = max(max_s, tot)
        return max_s if max_s > 0 else 10_000_000

    max_mac = get_max_size(MACRO_TOKENS)
    max_mic = get_max_size(MICRO_TOKENS)
    max_dot = get_max_size(DOT_TOKENS)

    fig_w = 14.5
    fig_h = 18.5
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')

    # Vertical Height Ratios:
    # Row 0: Top stats (a, b, c) -> 0.95 (compact)
    # Rows 1–2: Macro & Micro (d, e) -> 2.35 each (3 chromosomes each)
    # Row 3: Dot (f) -> 2.95 (4 chromosomes)
    # Row 4: Panel g + Legend -> 0.55 (compact)
    outer = GridSpec(5, 1, figure=fig,
                     height_ratios=[0.95, 2.35, 2.35, 2.95, 0.55],
                     hspace=0.15,
                     left=0.04, right=0.98, top=0.97, bottom=0.03)

    # ── Row 0: Panels a, b, c (Compact Top Row) ─────────────────────────────
    inner_top = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0, 0],
                                        width_ratios=[0.30, 0.70], wspace=0.10)
    ax_cov = fig.add_subplot(inner_top[0, 1])

    # Panel a reduced by 40% (height ratio 1.3 to 1.0)
    inner_stats = GridSpecFromSubplotSpec(2, 1, subplot_spec=inner_top[0, 0],
                                          height_ratios=[1.3, 1.0], hspace=0.38)
    inner_viol = GridSpecFromSubplotSpec(1, 3, subplot_spec=inner_stats[0], wspace=0.35)
    ax_viol = {
        'Macro': fig.add_subplot(inner_viol[0, 0]),
        'Micro': fig.add_subplot(inner_viol[0, 1]),
        'Dot':   fig.add_subplot(inner_viol[0, 2]),
    }
    ax_f_telo = fig.add_subplot(inner_stats[1])

    _draw_panel_d_violin(ax_viol, all_macro, all_micro, all_dot,
                         chrom_groups, cov_sum_s, cov_sum_d, cov_sum_o)
    if df_stats is not None:
        _draw_panel_f(ax_f_telo, df_stats)
    _draw_full_coverage_panel(ax_cov, all_macro, all_micro, all_dot,
                             chrom_groups, cov_sum_s, cov_sum_d, cov_sum_o)

    # ── Row 1: Macrochromosomes (Panel d, Extended) ──────────────────────────
    inner_mac = GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[1, 0],
                                        width_ratios=[0.038, 0.481, 0.481], wspace=0.015)
    ax_mac_chr = fig.add_subplot(inner_mac[0, 0])
    ax_mac_mat = fig.add_subplot(inner_mac[0, 1])
    ax_mac_pat = fig.add_subplot(inner_mac[0, 2])
    ax_mac_mat.set_title('Maternal', fontsize=12.5, fontweight='bold', pad=4, color='#1E293B')
    ax_mac_pat.set_title('Paternal', fontsize=12.5, fontweight='bold', pad=4, color='#1E293B')

    draw_synteny_panel(ax_mac_chr, ax_mac_mat, ax_mac_pat, MACRO_TOKENS, max_mac,
                       t2t_sizes, centromeres, t2t_telo, flip_set,
                       single_data, dual_data, ont_data,
                       tick_step=10_000_000, tick_unit_mb=True)

    # ── Row 2: Microchromosomes (Panel e, Extended) ──────────────────────────
    inner_mic = GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[2, 0],
                                        width_ratios=[0.038, 0.481, 0.481], wspace=0.015)
    ax_mic_chr = fig.add_subplot(inner_mic[0, 0])
    ax_mic_mat = fig.add_subplot(inner_mic[0, 1])
    ax_mic_pat = fig.add_subplot(inner_mic[0, 2])

    draw_synteny_panel(ax_mic_chr, ax_mic_mat, ax_mic_pat, MICRO_TOKENS, max_mic,
                       t2t_sizes, centromeres, t2t_telo, flip_set,
                       single_data, dual_data, ont_data,
                       tick_step=2_000_000, tick_unit_mb=True)

    # ── Row 3: Dot chromosomes (Panel f, Extended) ───────────────────────────
    inner_dot = GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[3, 0],
                                        width_ratios=[0.038, 0.481, 0.481], wspace=0.015)
    ax_dot_chr = fig.add_subplot(inner_dot[0, 0])
    ax_dot_mat = fig.add_subplot(inner_dot[0, 1])
    ax_dot_pat = fig.add_subplot(inner_dot[0, 2])

    draw_synteny_panel(ax_dot_chr, ax_dot_mat, ax_dot_pat, DOT_TOKENS, max_dot,
                       t2t_sizes, centromeres, t2t_telo, flip_set,
                       single_data, dual_data, ont_data,
                       tick_step=1_000_000, tick_unit_mb=True)

    # ── Row 4: Panel g (Left 44%) + Legend Card (Right 56%) (Compact) ────────
    inner_bot = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[4, 0],
                                        width_ratios=[0.44, 0.56], wspace=0.06)

    # Left: Panel g with 2 full-width barplots
    inner_g = GridSpecFromSubplotSpec(1, 2, subplot_spec=inner_bot[0, 0],
                                      width_ratios=[0.50, 0.50], wspace=0.28)
    ax_g_left  = fig.add_subplot(inner_g[0, 0])
    ax_g_right = fig.add_subplot(inner_g[0, 1])

    if annotation_df is not None:
        _draw_panel_g(ax_g_left, ax_g_right, annotation_df)

    # Right: Neater, Tidier Publication Legend Card
    ax_leg = fig.add_subplot(inner_bot[0, 1])
    _draw_legend_card(ax_leg)

    # ── Panel Lettering a–g ──────────────────────────────────────────────────
    panel_letters = [
        ('a', ax_viol['Macro'], -0.25, 1.08),
        ('b', ax_f_telo,        -0.08, 1.08),
        ('c', ax_cov,           -0.03, 1.04),
        ('d', ax_mac_chr,       -0.45, 1.02),
        ('e', ax_mic_chr,       -0.45, 1.02),
        ('f', ax_dot_chr,       -0.45, 1.02),
    ]
    if annotation_df is not None:
        panel_letters.append(('g', ax_g_left, -0.25, 1.08))

    for tag, ax_ref, lx, ly in panel_letters:
        ax_ref.text(lx, ly, tag, transform=ax_ref.transAxes,
                    fontsize=12.5, fontweight='bold', va='bottom', ha='right', color='#0F172A')

    # ── Save Outputs ─────────────────────────────────────────────────────────
    print(f"Saving PNG: {out_png}")
    fig.savefig(out_png, dpi=dpi, facecolor='white', bbox_inches='tight')
    if out_pdf:
        print(f"Saving PDF: {out_pdf}")
        fig.savefig(out_pdf, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    print("Combined multi-panel figure generated successfully!")

# =============================================================================
# CLI Entry Point
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Multi-panel assembly QC figure generator with synteny ideograms.")
    parser.add_argument('--t2t-tsv', required=True)
    parser.add_argument('--single-tsv', required=True)
    parser.add_argument('--dual-tsv', required=True)
    parser.add_argument('--ont-tsv', required=True)

    parser.add_argument('--single-pairs', required=True)
    parser.add_argument('--dual-pairs', required=True)
    parser.add_argument('--ont-pairs', '--ont-dual-pairs', dest='ont_pairs', required=True)

    parser.add_argument('--single-chain', required=True)
    parser.add_argument('--single-nc-chain', required=True)
    parser.add_argument('--dual-chain', required=True)
    parser.add_argument('--dual-nc-chain', required=True)
    parser.add_argument('--ont-chain', '--ont-dual-chain', dest='ont_chain', required=True)
    parser.add_argument('--ont-nc-chain', '--ont-dual-nc-chain', dest='ont_nc_chain', required=True)

    parser.add_argument('--single-rec-chain', default=None)
    parser.add_argument('--dual-rec-chain', default=None)
    parser.add_argument('--ont-rec-chain', default=None)

    parser.add_argument('--single-gaps', '--single-annotated-gaps', dest='single_gaps', required=True)
    parser.add_argument('--dual-gaps', '--dual-annotated-gaps', dest='dual_gaps', required=True)
    parser.add_argument('--ont-gaps', '--ont-annotated-gaps', dest='ont_gaps', required=True)

    parser.add_argument('--single-bed', required=True)
    parser.add_argument('--dual-bed', required=True)
    parser.add_argument('--ont-bed', required=True)

    parser.add_argument('--single-telomeres', '--single-telomere-presence', dest='single_telomeres', required=True)
    parser.add_argument('--dual-telomeres', '--dual-telomere-presence', dest='dual_telomeres', required=True)
    parser.add_argument('--ont-telomeres', '--ont-telomere-presence', dest='ont_telomeres', required=True)

    parser.add_argument('--centromeres', required=True)
    parser.add_argument('--telo-p-bed', required=True)

    parser.add_argument('--single-cov-sum', '--single-coverage-summary', dest='single_cov_sum', default=None)
    parser.add_argument('--dual-cov-sum', '--dual-coverage-summary', dest='dual_cov_sum', default=None)
    parser.add_argument('--ont-cov-sum', '--ont-coverage-summary', dest='ont_cov_sum', default=None)

    parser.add_argument('--stats-tsv', default=None)
    parser.add_argument('--annotation-tsv', default=None)

    parser.add_argument('--output', required=True)
    parser.add_argument('--format', default='png')
    parser.add_argument('--dpi', type=int, default=300)
    parser.add_argument('--no-timestamp', action='store_true')
    parser.add_argument('--style', default='paper')
    parser.add_argument('--simplify', action='store_true')

    args = parser.parse_args()

    print("Loading data for Multi-Panel Combined QC Figure...")

    # Load sequence sizes
    t2t_sizes = load_seq_sizes(args.t2t_tsv)
    single_sizes = load_seq_sizes(args.single_tsv)
    dual_sizes = load_seq_sizes(args.dual_tsv)
    ont_sizes = load_seq_sizes(args.ont_tsv)

    # Load chains and ribbons
    s_col_rib, s_col_sp, s_col_ins = parse_chain_detailed(args.single_chain, is_collinear=True)
    s_nc_rib, s_nc_sp, s_nc_ins    = parse_chain_detailed(args.single_nc_chain, is_collinear=False)
    s_rec_rib, s_rec_sp, s_rec_sz  = load_recovered_chains(args.single_rec_chain)
    single_sizes.update(s_rec_sz)
    s_telo_coll, s_telo_nc = load_telomere_presence_tsv(args.single_telomeres)

    single_data = {
        'col_ribbons': s_col_rib, 'nc_ribbons': s_nc_rib, 'rec_ribbons': s_rec_rib,
        'col_spans': s_col_sp, 'nc_spans': s_nc_sp, 'rec_spans': s_rec_sp,
        'sizes': single_sizes, 'pairs': load_chrom_pairs(args.single_pairs),
        'sw': load_switch_blocks(args.single_bed),
        'gaps': load_gaps_bed(args.single_gaps),
        'telo_coll': s_telo_coll, 'telo_nc': s_telo_nc
    }

    d_col_rib, d_col_sp, d_col_ins = parse_chain_detailed(args.dual_chain, is_collinear=True)
    d_nc_rib, d_nc_sp, d_nc_ins    = parse_chain_detailed(args.dual_nc_chain, is_collinear=False)
    d_rec_rib, d_rec_sp, d_rec_sz  = load_recovered_chains(args.dual_rec_chain)
    dual_sizes.update(d_rec_sz)
    d_telo_coll, d_telo_nc = load_telomere_presence_tsv(args.dual_telomeres)

    dual_data = {
        'col_ribbons': d_col_rib, 'nc_ribbons': d_nc_rib, 'rec_ribbons': d_rec_rib,
        'col_spans': d_col_sp, 'nc_spans': d_nc_sp, 'rec_spans': d_rec_sp,
        'sizes': dual_sizes, 'pairs': load_chrom_pairs(args.dual_pairs),
        'sw': load_switch_blocks(args.dual_bed),
        'gaps': load_gaps_bed(args.dual_gaps),
        'telo_coll': d_telo_coll, 'telo_nc': d_telo_nc
    }

    o_col_rib, o_col_sp, o_col_ins = parse_chain_detailed(args.ont_chain, is_collinear=True)
    o_nc_rib, o_nc_sp, o_nc_ins    = parse_chain_detailed(args.ont_nc_chain, is_collinear=False)
    o_rec_rib, o_rec_sp, o_rec_sz  = load_recovered_chains(args.ont_rec_chain)
    ont_sizes.update(o_rec_sz)
    o_telo_coll, o_telo_nc = load_telomere_presence_tsv(args.ont_telomeres)

    ont_data = {
        'col_ribbons': o_col_rib, 'nc_ribbons': o_nc_rib, 'rec_ribbons': o_rec_rib,
        'col_spans': o_col_sp, 'nc_spans': o_nc_sp, 'rec_spans': o_rec_sp,
        'sizes': ont_sizes, 'pairs': load_chrom_pairs(args.ont_pairs),
        'sw': load_switch_blocks(args.ont_bed),
        'gaps': load_gaps_bed(args.ont_gaps),
        'telo_coll': o_telo_coll, 'telo_nc': o_telo_nc
    }

    # Merge sizes into t2t_sizes if needed
    for sz_map in (single_sizes, dual_sizes, ont_sizes):
        for k, v in sz_map.items():
            if k not in t2t_sizes and ('Mat_NC' in k or 'Pat_NC' in k):
                t2t_sizes[k] = v

    centromeres = load_centromeres(args.centromeres)
    p_arm_dict  = load_p_arm_bed(args.telo_p_bed)
    t2t_telo    = load_t2t_terminal_telomeres(args.telo_p_bed)

    cov_sum_s, cov_sum_d, cov_sum_o = load_coverage_summary(
        args.single_cov_sum, args.dual_cov_sum, args.ont_cov_sum
    )

    df_stats = None
    if args.stats_tsv and os.path.exists(args.stats_tsv):
        df_stats = pd.read_csv(args.stats_tsv, sep='\t')

    annotation_df = _load_annotation_tsv(args.annotation_tsv) if args.annotation_tsv else None

    # Determine output file paths
    base_out = args.output
    if base_out.endswith('.png') or base_out.endswith('.pdf'):
        base_out = os.path.splitext(base_out)[0]

    out_png = f"{base_out}.png"
    out_pdf = f"{base_out}.pdf"

    build_combined_figure(
        t2t_sizes, centromeres, p_arm_dict, t2t_telo,
        single_data, dual_data, ont_data,
        cov_sum_s, cov_sum_d, cov_sum_o,
        df_stats, annotation_df,
        out_png=out_png, out_pdf=out_pdf, dpi=args.dpi
    )

if __name__ == '__main__':
    main()
