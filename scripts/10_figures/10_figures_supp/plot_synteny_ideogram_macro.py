#!/usr/bin/env python3
"""
plot_synteny_ideogram_macro.py

Generates compact publication-quality supplementary butterfly ideograms with synteny maps
between the T2T reference and assembled haplotypes (HiFi Single, HiFi Dual, and ONT Dual)
for macrochromosomes (Chr 1, 1A, 2, 3, 4, 4A, 5, 6, 7, 8, and W/Z), plus an averaged
coverage barplot panel matching the classical supplementary figure composition.

Features:
  1. Classical cytogenetic ideogram format:
     - Hourglass centromere constriction waist
     - Rounded semicircular telomere caps
     - Collinear coverage bands (dark saturated: Single #00BCCC, Dual #6445B0, ONT #2E8B52)
     - Non-collinear coverage bands (light pastel: Single #A0D4DC, Dual #C4B8E8, ONT #B8E8CC)
     - Clean white uncovered / insertion background
  2. Chain gap / insertion detection (dq >= 20 kb):
     - Collinear synteny ribbons split around large query insertions
     - Insertions in assembly sequence absent from T2T remain unaligned/white
  3. Appending unlinked scaffolds (>= 20 kb alignment) at chromosome ends:
     - Unlinked scaffolds receive their own ideogram blocks at chromosome ends
     - Non-collinear ribbons connect from T2T directly to unlinked scaffolds
  4. Native feature injection:
     - Curation gaps (red tick marks)
     - Assembly contig gaps (charcoal tick marks)
     - Hap-mer switch error blocks (red blocks)
     - Verified telomeres (solid for collinear, hollow for non-collinear)
  5. Compact 4-column composition:
     - Col 0: Chromosome label strip (1, 1A, 2, 3, ..., W/Z)
     - Col 1: Maternal Haplotype (butterfly mirrored: 0 at right)
     - Col 2: Paternal Haplotype (0 at left)
     - Col 3: Coverage panel (0 to 100% horizontal stacked bars per assembly)
     - Alternating subtle zebra row shading (#F4F8FA) across all columns
"""

import argparse
import os
import re
import sys
from collections import defaultdict
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.path as mpath
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import numpy as np
import pandas as pd


# =============================================================================
# Constants & Colors
# =============================================================================

MACRO_TOKENS = ['1', '1A', '2', '3', '4', '4A', '5', '6', '7', '8', 'ZW']

INSERTION_GAP_THRESHOLD = 20_000   # 20 kb gap in query (dq) splits collinear ribbon
MIN_NC_RIBBON_BP        = 20_000   # 20 kb minimum alignment to display non-collinear synteny ribbon

# Palette
COL_BORDER       = '#2c3e50'
COL_T2T_FILL     = '#F8FAFC'
COL_UNALIGNED    = '#FFFFFF'

# Collinear (Dark) & Non-Collinear (Light) Coverage Fills
COL_COV_S_DARK   = '#00BCCC'   # Single assembly
COL_COV_S_LIGHT  = '#A0D4DC'

COL_COV_D_DARK   = '#6445B0'   # Dual assembly
COL_COV_D_LIGHT  = '#C4B8E8'

COL_COV_O_DARK   = '#2E8B52'   # ONT Dual assembly
COL_COV_O_LIGHT  = '#B8E8CC'

# Synteny Ribbons
COL_RIB_COLL     = '#38BDF8'   # Sky blue
COL_RIB_COLL_E   = '#0284C7'
COL_RIB_NC       = '#F59E0B'   # Amber orange
COL_RIB_NC_E     = '#B45309'
COL_RIB_REC      = '#E11D48'   # Red for recovered query-1x chains
COL_RIB_REC_E    = '#9F1239'   # Darker red border
COL_COV_REC      = '#E11D48'   # Red band on assembly ideogram
ALPHA_RIB_REC    = 0.35        # 35% opacity as requested

# Native Feature Markers
COL_CUR_GAP      = '#E63946'   # Curation gap tick mark
COL_ASM_GAP      = '#111111'   # Assembly contig gap tick mark
COL_SWITCH       = '#D946EF'   # Hap-mer switch error block (magenta)
COL_SWITCH_EDGE  = '#701A75'   # Hap-mer switch border (deep purple/magenta)
COL_TELO         = '#C4426A'   # Cytogenetic telomere cap


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
            if len(p) >= 2 and p[0] != 'name' and p[1] != 'size':
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
                    ocs, oce = centromeres[key]
                    centromeres[key] = (min(ocs, cs), max(oce, ce))
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
    for chrom, size in seq_sizes.items():
        m = re.match(r'^(Pat|Mat)_.*_chromosome_([0-9]+[A-Za-z]*|[A-Za-z]+)$', chrom)
        if not m:
            continue
        side, tok = m.group(1).lower(), m.group(2).upper()
        if chrom in p_starts:
            if p_starts[chrom] > size / 2:
                flip_set.add(chrom)
        elif (tok, side) in centromeres:
            cs, ce = centromeres[(tok, side)]
            cen_mid = (cs + ce) / 2
            if cen_mid > size / 2:
                flip_set.add(chrom)
    return flip_set


def load_t2t_terminal_telomeres(bed_path):
    """Loads terminal telomeres present on T2T reference chromosomes."""
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


def load_gaps_bed(gaps_bed):
    """Loads annotated gaps (curation joins and contig gaps)."""
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
    """Loads hap-mer switch error blocks."""
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
    """Loads telomere presence annotations (collinear and non-collinear)."""
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
            if col_val != 'none':
                if 'p' in col_val: arms_c.add('p')
                if 'q' in col_val: arms_c.add('q')

            arms_nc = set()
            if nc_val != 'none':
                if 'p' in nc_val: arms_nc.add('p')
                if 'q' in nc_val: arms_nc.add('q')

            for key in (chrom, scaff):
                if key:
                    telo_coll[key].update(arms_c)
                    telo_nc[key].update(arms_nc)
    return telo_coll, telo_nc


def load_coverage_summary(path):
    """Loads per-chromosome coverage summary mapping chromosome token to (col_pct, nc_pct, unc_pct)."""
    cov = defaultdict(list)
    if not path or not os.path.exists(path):
        return cov
    df = pd.read_csv(path, sep='\t', comment='#')
    for _, row in df.iterrows():
        chrom_tok = str(row['chrom']).strip()
        c = float(row['collinear_pct'])
        nc = float(row['noncollinear_pct'])
        u = float(row['uncovered_pct'])
        cov[chrom_tok].append((c, nc, u))
    return cov


def get_avg_cov_for_chrom(tok, cov_dict):
    """Computes average (col_pct, nc_pct, unc_pct) for chromosome token across haplotypes."""
    if not cov_dict:
        return (0.0, 0.0, 100.0)
    if tok == 'ZW':
        vals = cov_dict.get('W', []) + cov_dict.get('Z', [])
    else:
        vals = cov_dict.get(tok, [])
    if not vals:
        return (0.0, 0.0, 100.0)
    avg_c = sum(v[0] for v in vals) / len(vals)
    avg_nc = sum(v[1] for v in vals) / len(vals)
    avg_u = max(0.0, 100.0 - avg_c - avg_nc)
    return (avg_c, avg_nc, avg_u)


def parse_chain_detailed(chain_path, is_collinear=True, min_nc_size=MIN_NC_RIBBON_BP):
    """
    Parses UCSC chain file into ribbons, splitting across query gaps (dq >= INSERTION_GAP_THRESHOLD).
    Filters out spurious micro-alignments (< MIN_NC_RIBBON_BP) for non-collinear chains.
    Returns:
      ribbons: list of (t_name, t_start, t_end, q_name, q_start, q_end, strand, is_collinear)
      query_spans: dict of q_name -> list of (q_start, q_end)
      insertions: dict of q_name -> list of (q_ins_start, q_ins_end, dq_size, t_pos)
    """
    ribbons = []
    query_spans = defaultdict(list)
    insertions = defaultdict(list)

    if not chain_path or not os.path.exists(chain_path):
        return ribbons, query_spans, insertions

    def _add_ribbon(t_n, ts, te, q_n, qs, qe, qstr):
        if te <= ts or qe <= qs:
            return
        if not is_collinear and max(te - ts, qe - qs) < min_nc_size:
            return
        ribbons.append((t_n, ts, te, q_n, qs, qe, qstr, is_collinear))
        query_spans[q_n].append((qs, qe))

    with open(chain_path) as f:
        t_name = q_name = q_strand = None
        t_start = t_end = q_start = q_end = 0
        q_size = 0
        cur_t_start = cur_t_end = 0
        cur_q_start = cur_q_end = 0
        in_ribbon = False

        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
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
                    cur_t_start = blk_t_s
                    cur_t_end   = blk_t_e
                    cur_q_start = blk_q_s
                    cur_q_end   = blk_q_e
                    in_ribbon = True
                else:
                    cur_t_end = blk_t_e
                    cur_q_end = blk_q_e

                t_pos += size + dt
                q_pos += size + dq

                # If query gap (insertion) or target gap (deletion/unaligned) is >= threshold, finish current ribbon
                if (dq >= INSERTION_GAP_THRESHOLD or dt >= INSERTION_GAP_THRESHOLD) and in_ribbon:
                    _add_ribbon(t_name, cur_t_start, cur_t_end, q_name, cur_q_start, cur_q_end, q_strand)
                    if dq >= INSERTION_GAP_THRESHOLD:
                        insertions[q_name].append((blk_q_e, blk_q_e + dq, dq, blk_t_e))
                    in_ribbon = False

        if in_ribbon:
            _add_ribbon(t_name, cur_t_start, cur_t_end, q_name, cur_q_start, cur_q_end, q_strand)

    return ribbons, query_spans, insertions


def _merge_intervals(intervals, min_gap=5_000):
    """Merges adjacent or overlapping intervals."""
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


# =============================================================================
# Cytogenetic Drawing Utilities
# =============================================================================

def _rect_path(x0, x1, yc, h):
    """Straight rectangular path."""
    y0, y1 = yc - h / 2, yc + h / 2
    verts = [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    return Path(verts, codes)


def _constriction_path(x0, x1, yc, h, c_start, c_end):
    """Cytogenetic ideogram path with an hourglass centromere constriction waist."""
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


def _draw_telo_semi(ax, x_edge, side, yc, h, rx, fc=None, ec=None, lw=0.6, hollow=False):
    """Draws a rounded semicircular telomere cap at chromosome terminus."""
    y0 = yc - h / 2
    y1 = yc + h / 2
    n_pts = 16
    if side == 'p':
        theta = np.linspace(np.pi / 2, 3 * np.pi / 2, n_pts)
        pts = [(x_edge + rx * np.cos(t), yc + (h / 2) * np.sin(t)) for t in theta]
        pts = [(x_edge, y1)] + pts + [(x_edge, y0)]
    else:
        theta = np.linspace(-np.pi / 2, np.pi / 2, n_pts)
        pts = [(x_edge + rx * np.cos(t), yc + (h / 2) * np.sin(t)) for t in theta]
        pts = [(x_edge, y0)] + pts + [(x_edge, y1)]

    verts = pts + [pts[0]]
    codes = [Path.MOVETO] + [Path.LINETO] * (len(pts) - 1) + [Path.CLOSEPOLY]
    path = Path(verts, codes)

    if fc is None:
        fc = 'white' if hollow else '#000000'
    if ec is None:
        ec = '#000000'
    ax.add_patch(PathPatch(path, facecolor=fc, edgecolor=ec, lw=lw, zorder=15))


def draw_synteny_ribbon(ax, x0_t, x1_t, y_t, x0_q, x1_q, y_q, col, edge_col, alpha=0.45, inverted=False):
    """Draws a cubic Bezier ribbon connecting T2T interval to assembly interval."""
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
                     linewidth=0.4, alpha=alpha, zorder=3)
    ax.add_patch(poly)


# =============================================================================
# Main Multi-Panel Plotting Engine
# =============================================================================

def plot_butterfly_macro(
    t2t_sizes, flip_set,
    single_data, dual_data, ont_data,
    centromeres, t2t_telo,
    single_cov_sum, dual_cov_sum, ont_cov_sum,
    out_png, out_pdf, dpi=300
):
    """
    Renders the compact publication-quality macrochromosome synteny ideogram figure.
    Includes Maternal, Paternal, and Coverage barplot panels with classic zebra row striping.
    """
    tokens = MACRO_TOKENS
    n_tokens = len(tokens)

    # 1. Global X scale
    max_global_size = max(t2t_sizes.values()) if t2t_sizes else 160_000_000
    for asm_d in (single_data, dual_data, ont_data):
        if asm_d['sizes']:
            max_global_size = max(max_global_size, max(asm_d['sizes'].values()))

    # Ensure max_global_size accommodates primary scaffold + unloc scaffolds + spacers
    SPACER = 2_200_000
    for tok in tokens:
        for side in ('mat', 'pat'):
            if tok == 'ZW':
                t2t_c = 'Mat_NC_133064.1_chromosome_W' if side == 'mat' else 'Mat_NC_133063.1_chromosome_Z'
            else:
                t2t_c = next((k for k in t2t_sizes.keys() if f'chromosome_{tok}' in k and (side == 'pat' if ('Pat' in k or tok == 'Z') else side == 'mat')), None)
            if not t2t_c: continue
            for asm_d in (single_data, dual_data, ont_data):
                qp = asm_d['pairs'].get(t2t_c)
                if not qp: continue
                tot = asm_d['sizes'].get(qp, 0)
                parts = qp.rsplit('.', 1)
                stem = parts[0]
                suf = f".{parts[1]}" if len(parts) > 1 else ""
                u_set = set()
                for r in asm_d['ribbons'].get(t2t_c, []) + asm_d.get('rec_ribbons', {}).get(t2t_c, []):
                    if r[3] != qp: u_set.add(r[3])
                for s in asm_d['sizes']:
                    if s.startswith(f"{stem}_unloc_") and (not suf or s.endswith(suf)):
                        u_set.add(s)
                for u in u_set:
                    tot += SPACER + asm_d['sizes'].get(u, 500_000)
                max_global_size = max(max_global_size, tot)

    max_global_size = max_global_size + 6_000_000
    x_margin = max_global_size * 0.02

    # 2. Compact Geometry Layout
    H_BAR_T2T = 0.10
    H_BAR_ASM = 0.10
    H_RIBBON  = 0.14
    BLOCK_GAP = 0.05
    BLOCK_H   = H_BAR_T2T + H_RIBBON + H_BAR_ASM + BLOCK_GAP  # 0.39
    ROW_GAP   = 0.20
    ROW_H     = 3 * BLOCK_H + ROW_GAP                         # 1.37

    total_plot_height = n_tokens * ROW_H + 0.9

    # 3. Figure & GridSpec
    fig_w, fig_h = 22.0, 16.0
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')

    # GridSpec: [ Chr (Left) | Maternal | Paternal | Right Column ]
    gs = GridSpec(1, 4, figure=fig, width_ratios=[0.26, 3.65, 3.65, 1.44],
                  wspace=0.035, left=0.035, right=0.98, top=0.94, bottom=0.06)

    ax_chr = fig.add_subplot(gs[0, 0])
    ax_mat = fig.add_subplot(gs[0, 1])
    ax_pat = fig.add_subplot(gs[0, 2])

    # Right column split: [ Compact Coverage (rows 1-4) | Legend (rows 5-11) ]
    gs_right = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0, 3], height_ratios=[4.2, 6.8], hspace=0.20)
    ax_cov = fig.add_subplot(gs_right[0, 0])
    ax_leg = fig.add_subplot(gs_right[1, 0])
    ax_leg.axis('off')

    y_top = total_plot_height
    y_bot = 0.0

    for ax in (ax_chr, ax_mat, ax_pat):
        ax.set_ylim(y_bot, y_top)
        ax.set_yticks([])
        ax.set_xticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    ax_chr.set_xlim(0, 1)
    ax_mat.set_xlim(max_global_size + x_margin, -x_margin)
    ax_pat.set_xlim(-x_margin, max_global_size + x_margin)

    telo_rx = max_global_size * 0.0050

    # 4. Loop over Chromosome Rows
    for row_idx, tok in enumerate(tokens):
        y_row_top = y_top - 0.45 - row_idx * ROW_H

        # Zebra background shading for even rows
        if row_idx % 2 == 0:
            y_bg_top = y_row_top + 0.07
            y_bg_bot = y_row_top - (3 * BLOCK_H) - 0.05
            for ax in (ax_chr, ax_mat, ax_pat):
                ax.axhspan(y_bg_bot, y_bg_top, color='#F4F8FA', zorder=0)

        # Row divider line for subsequent rows
        if row_idx > 0:
            y_div = y_row_top + ROW_GAP / 2
            ax_mat.axhline(y_div, color='#E2E8F0', lw=0.6, ls=':')
            ax_pat.axhline(y_div, color='#E2E8F0', lw=0.6, ls=':')

        # Chromosome label in Col 0 (ax_chr)
        lbl_text = 'W/Z' if tok == 'ZW' else tok
        y_mid_row = y_row_top - 1.5 * BLOCK_H + BLOCK_GAP / 2
        ax_chr.text(0.70, y_mid_row, lbl_text, ha='center', va='center',
                    fontsize=11.5, fontweight='bold', color='#0F172A')

        # 3 Assembly Blocks per Chromosome
        assemblies = [
            ('HiFi Single', single_data, single_cov_sum, COL_COV_S_DARK, COL_COV_S_LIGHT, '#1D4ED8'),
            ('HiFi Dual',   dual_data,   dual_cov_sum,   COL_COV_D_DARK, COL_COV_D_LIGHT, '#6445B0'),
            ('ONT Dual',    ont_data,    ont_cov_sum,    COL_COV_O_DARK, COL_COV_O_LIGHT, '#2E8B52')
        ]

        for b_idx, (asm_title, asm, cov_sum_dict, c_dark, c_light, c_lbl) in enumerate(assemblies):
            y_block_top = y_row_top - b_idx * BLOCK_H
            y_t2t = y_block_top - H_BAR_T2T / 2
            y_asm = y_t2t - (H_BAR_T2T / 2 + H_RIBBON + H_BAR_ASM / 2)

            # Draw Coverage bar in compact ax_cov
            cp, ncp, up = get_avg_cov_for_chrom(tok, cov_sum_dict)
            y_cov_c = (len(tokens) - 1) - row_idx
            bar_h = 0.22
            y_cov_bar = y_cov_c + (1 - b_idx) * bar_h
            ax_cov.barh(y_cov_bar, cp, height=bar_h * 0.90, left=0, color=c_dark, ec='none', zorder=4)
            ax_cov.barh(y_cov_bar, ncp, height=bar_h * 0.90, left=cp, color=c_light, ec='none', zorder=4)
            ax_cov.barh(y_cov_bar, up, height=bar_h * 0.90, left=cp + ncp, color='#FFFFFF', ec='none', zorder=4)
            ax_cov.add_patch(mpatches.Rectangle((0, y_cov_bar - (bar_h * 0.90) / 2), 100, bar_h * 0.90,
                                                fc='none', ec='#94A3B8', lw=0.4, zorder=5))

            # Render Maternal (left) and Paternal (right)
            for side, ax in [('mat', ax_mat), ('pat', ax_pat)]:
                # Resolve T2T chromosome
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

                # Primary scaffold
                q_prim = asm['pairs'].get(t2t_chrom)
                q_prim_len = asm['sizes'].get(q_prim, t2t_len) if q_prim else t2t_len

                # Find unlinked scaffolds aligning to this chromosome or associated with q_prim
                all_asm_ribs = asm['ribbons'].get(t2t_chrom, [])
                unloc_aligned = defaultdict(int)
                for r in all_asm_ribs:
                    qn = r[3]
                    if qn != q_prim:
                        unloc_aligned[qn] += abs(r[5] - r[4])

                # Collect unlinked scaffolds:
                # 1. Any scaffold aligning to this chromosome (no minimum threshold)
                # 2. Any unlocalized scaffold belonging to q_prim by naming pattern
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

                def _unloc_sort_key(u_name):
                    align_bp = unloc_aligned.get(u_name, 0)
                    return (-align_bp, u_name)

                unloc_qual = sorted(unloc_candidates, key=_unloc_sort_key)

                # Compute X-offsets for primary and unlinked scaffolds
                SPACER = 2_200_000  # 2.2 Mb spacer between primary and unlinked
                offsets = {q_prim: 0}
                cur_x = q_prim_len
                for u in unloc_qual:
                    cur_x += SPACER
                    offsets[u] = cur_x
                    cur_x += asm['sizes'].get(u, 500_000)

                # =============================================================
                # A. Top Bar: T2T Reference Ideogram
                # =============================================================
                tok_cen = 'W' if (tok == 'ZW' and side == 'mat') else ('Z' if (tok == 'ZW' and side == 'pat') else tok)
                cen_span = centromeres.get((tok_cen.upper(), side.lower()))

                if cen_span:
                    cs = (t2t_len - cen_span[1]) if flip else cen_span[0]
                    ce = (t2t_len - cen_span[0]) if flip else cen_span[1]
                    t2t_path = _constriction_path(0, t2t_len, y_t2t, H_BAR_T2T, min(cs, ce), max(cs, ce))
                else:
                    t2t_path = _rect_path(0, t2t_len, y_t2t, H_BAR_T2T)

                ax.add_patch(PathPatch(t2t_path, fc=COL_T2T_FILL, ec=COL_BORDER, lw=0.6, zorder=5))

                # T2T Telomeres (colored black)
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
                    ax.add_patch(PathPatch(prim_path, fc=COL_UNALIGNED, ec=COL_BORDER, lw=0.6, zorder=5))

                    # Non-collinear coverage bands (under collinear)
                    m_nc = _merge_intervals(asm['nc_spans'].get(q_prim, []))
                    for s, e in m_nc:
                        x0 = (q_prim_len - e) if flip else s
                        x1 = (q_prim_len - s) if flip else e
                        r = mpatches.Rectangle((min(x0, x1), y_asm - H_BAR_ASM / 2),
                                               abs(x1 - x0), H_BAR_ASM,
                                               fc=c_light, ec='none', alpha=0.85, zorder=6)
                        r.set_clip_path(prim_path, transform=ax.transData)
                        ax.add_patch(r)

                    # Recovered query-1x coverage bands (red)
                    m_rec = _merge_intervals(asm.get('rec_spans', {}).get(q_prim, []))
                    for s, e in m_rec:
                        x0 = (q_prim_len - e) if flip else s
                        x1 = (q_prim_len - s) if flip else e
                        r = mpatches.Rectangle((min(x0, x1), y_asm - H_BAR_ASM / 2),
                                               abs(x1 - x0), H_BAR_ASM,
                                               fc=COL_COV_REC, ec='none', alpha=0.85, zorder=6.5)
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

                    # Native hap-mer switch blocks
                    for sw_s, sw_e in asm['sw'].get(q_prim, []):
                        x0 = (q_prim_len - sw_e) if flip else sw_s
                        x1 = (q_prim_len - sw_s) if flip else sw_e
                        w_sw = max(abs(x1 - x0), 80_000)
                        r = mpatches.Rectangle((min(x0, x1), y_asm - H_BAR_ASM / 2),
                                               w_sw, H_BAR_ASM,
                                               fc=COL_SWITCH, ec='none', zorder=8)
                        r.set_clip_path(prim_path, transform=ax.transData)
                        ax.add_patch(r)

                    # Native gaps (curation joins and contig gaps)
                    for gs, ge, gtype in asm['gaps'].get(q_prim, []):
                        gx = (q_prim_len - (gs + ge) / 2) if flip else ((gs + ge) / 2)
                        if 'CURATION' in gtype:
                            ax.plot([gx, gx], [y_asm - H_BAR_ASM / 2, y_asm + H_BAR_ASM / 2],
                                    color=COL_CUR_GAP, lw=1.0, zorder=9, solid_capstyle='butt')
                        else:
                            ax.plot([gx, gx], [y_asm - H_BAR_ASM / 2, y_asm + H_BAR_ASM / 2],
                                    color=COL_ASM_GAP, lw=0.7, zorder=8, solid_capstyle='butt')

                    # Assembly Telomeres
                    asm_coll_arms = asm['telo_coll'].get(t2t_chrom, set()) | asm['telo_coll'].get(q_prim, set())
                    asm_nc_arms   = asm['telo_nc'].get(t2t_chrom, set()) | asm['telo_nc'].get(q_prim, set())

                    if p_arm in asm_coll_arms:
                        _draw_telo_semi(ax, 0, 'p', y_asm, H_BAR_ASM, telo_rx, fc=c_dark, ec=c_dark, hollow=False)
                    elif p_arm in asm_nc_arms:
                        _draw_telo_semi(ax, 0, 'p', y_asm, H_BAR_ASM, telo_rx, fc='white', ec=c_dark, lw=0.9, hollow=True)

                    if q_arm in asm_coll_arms:
                        _draw_telo_semi(ax, q_prim_len, 'q', y_asm, H_BAR_ASM, telo_rx, fc=c_dark, ec=c_dark, hollow=False)
                    elif q_arm in asm_nc_arms:
                        _draw_telo_semi(ax, q_prim_len, 'q', y_asm, H_BAR_ASM, telo_rx, fc='white', ec=c_dark, lw=0.9, hollow=True)

                # =============================================================
                # C. Unlinked Scaffolds at Chromosome End
                # =============================================================
                for u in unloc_qual:
                    u_off = offsets[u]
                    u_len = asm['sizes'].get(u, 500_000)

                    u_path = _rect_path(u_off, u_off + u_len, y_asm, H_BAR_ASM)
                    ax.add_patch(PathPatch(u_path, fc=COL_UNALIGNED, ec='#64748B', lw=0.6, ls='--', zorder=5))

                    # Non-collinear bands
                    u_nc = _merge_intervals(asm['nc_spans'].get(u, []))
                    for s, e in u_nc:
                        r = mpatches.Rectangle((u_off + s, y_asm - H_BAR_ASM / 2), e - s, H_BAR_ASM,
                                               fc=c_light, ec='none', alpha=0.85, zorder=6)
                        r.set_clip_path(u_path, transform=ax.transData)
                        ax.add_patch(r)

                    # Recovered bands on unlinked
                    u_rec = _merge_intervals(asm.get('rec_spans', {}).get(u, []))
                    for s, e in u_rec:
                        r = mpatches.Rectangle((u_off + s, y_asm - H_BAR_ASM / 2), e - s, H_BAR_ASM,
                                               fc=COL_COV_REC, ec='none', alpha=0.85, zorder=6.5)
                        r.set_clip_path(u_path, transform=ax.transData)
                        ax.add_patch(r)

                    # Collinear bands on unlinked
                    u_col = _merge_intervals(asm['col_spans'].get(u, []))
                    for s, e in u_col:
                        r = mpatches.Rectangle((u_off + s, y_asm - H_BAR_ASM / 2), e - s, H_BAR_ASM,
                                               fc=c_dark, ec='none', alpha=1.00, zorder=7)
                        r.set_clip_path(u_path, transform=ax.transData)
                        ax.add_patch(r)

                    # Gaps on unlinked
                    for gs, ge, gtype in asm['gaps'].get(u, []):
                        gx = u_off + (gs + ge) / 2
                        if 'CURATION' in gtype:
                            ax.plot([gx, gx], [y_asm - H_BAR_ASM / 2, y_asm + H_BAR_ASM / 2],
                                    color=COL_CUR_GAP, lw=1.0, zorder=9)
                        else:
                            ax.plot([gx, gx], [y_asm - H_BAR_ASM / 2, y_asm + H_BAR_ASM / 2],
                                    color=COL_ASM_GAP, lw=0.6, zorder=8)

                    # Label above unlinked scaffold
                    m_lbl = re.search(r'unloc_([0-9]+)', u)
                    s_tag = f"u{m_lbl.group(1)}" if m_lbl else "unloc"
                    ax.text(u_off + u_len / 2, y_asm + H_BAR_ASM * 0.70, s_tag,
                            ha='left', va='bottom', fontsize=5.0, color='#475569',
                            fontweight='bold', rotation=40)

                # =============================================================
                # D. Synteny Ribbons (Collinear & Non-Collinear)
                # =============================================================
                y_rib_t = y_t2t - H_BAR_T2T / 2
                y_rib_q = y_asm + H_BAR_ASM / 2

                # 1. Collinear Ribbons
                for r in asm['col_ribbons'].get(t2t_chrom, []):
                    _, ts, te, qn, qs, qe, qstr, _ = r
                    if qn not in offsets: continue
                    q_off = offsets[qn]
                    q_len_cur = asm['sizes'].get(qn, q_prim_len)

                    if qn == q_prim and flip:
                        t0, t1 = t2t_len - te, t2t_len - ts
                        q0, q1 = q_prim_len - qe, q_prim_len - qs
                    else:
                        t0, t1 = ts, te
                        q0, q1 = q_off + qs, q_off + qe

                    draw_synteny_ribbon(ax, min(t0, t1), max(t0, t1), y_rib_t,
                                        min(q0, q1), max(q0, q1), y_rib_q,
                                        c_dark, c_dark, alpha=0.25,
                                        inverted=(qstr == '-'))

                # 2. Non-collinear Ribbons (on top of collinear)
                for r in asm['nc_ribbons'].get(t2t_chrom, []):
                    _, ts, te, qn, qs, qe, qstr, _ = r
                    if qn not in offsets: continue
                    if max(te - ts, abs(qe - qs)) < MIN_NC_RIBBON_BP:
                        continue
                    q_off = offsets[qn]
                    q_len_cur = asm['sizes'].get(qn, q_prim_len)

                    if qn == q_prim and flip:
                        t0, t1 = t2t_len - te, t2t_len - ts
                        q0, q1 = q_prim_len - qe, q_prim_len - qs
                    else:
                        t0, t1 = ts, te
                        q0, q1 = q_off + qs, q_off + qe

                    draw_synteny_ribbon(ax, min(t0, t1), max(t0, t1), y_rib_t,
                                        min(q0, q1), max(q0, q1), y_rib_q,
                                        COL_RIB_NC, COL_RIB_NC_E, alpha=0.55,
                                        inverted=(qstr == '-'))

                # 3. Recovered 1x Ribbons (Red, alpha=0.35)
                for r in asm.get('rec_ribbons', {}).get(t2t_chrom, []):
                    _, ts, te, qn, qs, qe, qstr, _ = r
                    if qn not in offsets: continue
                    if max(te - ts, abs(qe - qs)) < MIN_NC_RIBBON_BP:
                        continue
                    q_off = offsets[qn]
                    q_len_cur = asm['sizes'].get(qn, q_prim_len)

                    if qn == q_prim and flip:
                        t0, t1 = t2t_len - te, t2t_len - ts
                        q0, q1 = q_prim_len - qe, q_prim_len - qs
                    else:
                        t0, t1 = ts, te
                        q0, q1 = q_off + qs, q_off + qe

                    draw_synteny_ribbon(ax, min(t0, t1), max(t0, t1), y_rib_t,
                                        min(q0, q1), max(q0, q1), y_rib_q,
                                        COL_RIB_REC, COL_RIB_REC_E, alpha=ALPHA_RIB_REC,
                                        inverted=(qstr == '-'))

    # 5. Column Headers & Scale Bars
    ax_mat.set_title('Maternal', fontsize=14, fontweight='bold', pad=12, color='#1E293B')
    ax_pat.set_title('Paternal', fontsize=14, fontweight='bold', pad=12, color='#1E293B')
    # Configure compact coverage panel in upper right column
    n_tok = len(tokens)
    ax_cov.set_ylim(-0.6, n_tok - 0.4)
    ax_cov.set_yticks(range(n_tok))
    y_labels = [('W/Z' if t == 'ZW' else t) for t in reversed(tokens)]
    ax_cov.set_yticklabels(y_labels, fontsize=8.0, fontweight='bold', color='#1E293B')
    ax_cov.tick_params(axis='y', length=2.5, width=0.6, color='#94A3B8')

    ax_cov.set_xlim(0, 100)
    ax_cov.set_xticks([0, 25, 50, 75, 100])
    ax_cov.set_xticklabels(['0', '25', '50', '75', '100%'], fontsize=7.5, fontweight='medium', color='#1E293B')
    ax_cov.set_xlabel('Coverage (%)', fontsize=8.5, fontweight='bold', labelpad=4, color='#1E293B')
    ax_cov.grid(axis='x', color='#E2E8F0', linestyle='--', linewidth=0.6, alpha=0.9, zorder=0)
    ax_cov.set_axisbelow(True)
    ax_cov.spines['top'].set_visible(False)
    ax_cov.spines['right'].set_visible(False)
    ax_cov.spines['left'].set_visible(True)
    ax_cov.spines['left'].set_color('#CBD5E1')
    ax_cov.spines['left'].set_linewidth(0.6)
    ax_cov.spines['bottom'].set_visible(True)
    ax_cov.spines['bottom'].set_color('#CBD5E1')
    ax_cov.spines['bottom'].set_linewidth(0.6)
    ax_cov.tick_params(axis='x', length=2.5, width=0.6, which='major', labelsize=7.5, color='#1E293B')
    ax_cov.set_title('Coverage', fontsize=11.5, fontweight='bold', pad=8, color='#1E293B')

    # Alternating row background shading in ax_cov
    for r_idx in range(n_tok):
        if r_idx % 2 == 0:
            y_c = (n_tok - 1) - r_idx
            ax_cov.axhspan(y_c - 0.45, y_c + 0.45, color='#F8FAFC', zorder=1)

    # Extended Coordinate Axis with Ticks every 20 Mb at Bottom
    axis_max = int(np.ceil(max_global_size / 20_000_000)) * 20_000_000
    y_axis = 0.28
    tick_h = 0.08

    for ax, side in [(ax_mat, 'mat'), (ax_pat, 'pat')]:
        # Continuous horizontal baseline
        ax.plot([0, axis_max], [y_axis, y_axis], color='#1E293B', lw=1.2, zorder=10)

        # Periodic ticks and labels every 20 Mb
        for tick_val in range(0, axis_max + 1, 20_000_000):
            mb_val = tick_val // 1_000_000
            ax.plot([tick_val, tick_val], [y_axis, y_axis + tick_h], color='#1E293B', lw=1.0, zorder=10)
            ax.text(tick_val, y_axis - 0.08, f'{mb_val}', ha='center', va='top',
                    fontsize=7.8, color='#1E293B', fontweight='medium')

        # Unit 'Mb' at outer margin
        lbl_x = axis_max + 4_000_000
        ha = 'right' if side == 'mat' else 'left'
        ax.text(lbl_x, y_axis, 'Mb', ha=ha, va='center', fontsize=8.2, fontweight='bold', color='#1E293B')

    # 6. Comprehensive Publication Legend (2 columns aligned with ax_cov on both sides)
    legend_elements = [
        # Row 1: Single & Dual Collinear
        mpatches.Patch(facecolor=COL_COV_S_DARK, edgecolor='none', label='Single Collinear'),
        mpatches.Patch(facecolor=COL_COV_D_DARK, edgecolor='none', label='Dual Collinear'),
        # Row 2: ONT Collinear & Non-collinear Coverage
        mpatches.Patch(facecolor=COL_COV_O_DARK, edgecolor='none', label='ONT Collinear'),
        mpatches.Patch(facecolor='#CBD5E1', edgecolor='none', label='Non-collinear Cov.'),
        # Row 3: Collinear & Non-collinear Ribbons
        mpatches.Patch(facecolor='#64748B', edgecolor='none', alpha=0.25, label='Collinear Ribbon'),
        mpatches.Patch(facecolor=COL_RIB_NC, edgecolor=COL_RIB_NC_E, lw=0.4, alpha=0.65, label='Non-collinear Ribbon'),
        # Row 4: Recovered 1x Ribbons & Coverage
        mpatches.Patch(facecolor=COL_RIB_REC, edgecolor=COL_RIB_REC_E, lw=0.4, alpha=ALPHA_RIB_REC, label='Recovered 1x Ribbon'),
        mpatches.Patch(facecolor=COL_COV_REC, edgecolor='none', label='Recovered 1x Cov.'),
        # Row 5: Unaligned insertion & Unlinked Scaffold
        mpatches.Patch(facecolor=COL_UNALIGNED, edgecolor=COL_BORDER, lw=0.7, label='Unaligned (≥20kb)'),
        mpatches.Patch(facecolor='none', edgecolor='#64748B', lw=0.8, ls='--', label='Unlinked Scaffold'),
        # Row 5: Curation Join & Assembly Gap
        plt.Line2D([0], [0], color=COL_CUR_GAP, lw=1.8, label='Curation Gap (Join)'),
        plt.Line2D([0], [0], color=COL_ASM_GAP, lw=1.2, label='Assembly Gap'),
        # Row 6: Switch Block & T2T Reference
        mpatches.Patch(facecolor=COL_SWITCH, edgecolor='none', label='Switch Block'),
        mpatches.Patch(facecolor=COL_T2T_FILL, edgecolor=COL_BORDER, lw=0.7, label='T2T Reference'),
        # Row 7: Assembly Telomeres (Collinear vs Non-collinear)
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#475569', markeredgecolor='#475569',
                   markersize=6, label='Collinear Telomere'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='white', markeredgecolor='#475569',
                   markeredgewidth=1.1, markersize=6, label='Non-collinear Telo.')
    ]

    leg = ax_leg.legend(handles=legend_elements, loc='upper left',
                        bbox_to_anchor=(0.0, 0.68, 1.0, 0.32), mode='expand',
                        ncol=2, frameon=True,
                        facecolor='#F8FAFC', edgecolor='#CBD5E1', fontsize=7.2,
                        handlelength=1.0, handleheight=0.7, handletextpad=0.5,
                        borderpad=0.7, labelspacing=0.6)
    leg.get_frame().set_linewidth(0.8)

    # 7. Save Outputs
    print(f"Saving PNG: {out_png}")
    fig.savefig(out_png, dpi=dpi, facecolor='white', bbox_inches='tight')
    if out_pdf:
        print(f"Saving PDF: {out_pdf}")
        fig.savefig(out_pdf, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    print("Butterfly synteny ideogram successfully generated!")


# =============================================================================
# CLI Main
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description='Publication-quality Macrochromosome Butterfly Synteny Ideogram')
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

    parser.add_argument('--single-rec-chain', default=None, help='Recovered uncovered chain file for Single')
    parser.add_argument('--dual-rec-chain', default=None, help='Recovered uncovered chain file for Dual')
    parser.add_argument('--ont-rec-chain', default=None, help='Recovered uncovered chain file for ONT')

    parser.add_argument('--single-gaps', required=True)
    parser.add_argument('--dual-gaps', required=True)
    parser.add_argument('--ont-gaps', required=True)

    parser.add_argument('--single-bed', required=True)
    parser.add_argument('--dual-bed', required=True)
    parser.add_argument('--ont-bed', required=True)

    parser.add_argument('--single-telomeres', required=True)
    parser.add_argument('--dual-telomeres', required=True)
    parser.add_argument('--ont-telomeres', required=True)

    parser.add_argument('--single-cov-sum', default=None)
    parser.add_argument('--dual-cov-sum', default=None)
    parser.add_argument('--ont-cov-sum', default=None)

    parser.add_argument('--centromeres', required=True)
    parser.add_argument('--telo-p-bed', required=True)

    parser.add_argument('--output-png', required=True)
    parser.add_argument('--output-pdf', default=None)
    parser.add_argument('--dpi', type=int, default=300)

    return parser.parse_args()


def load_assembly_package(name, tsv, pairs_file, col_chain, nc_chain, gaps_bed, sw_bed, telo_tsv, rec_chain=None):
    """Loads all data layers for a single assembly."""
    print(f"Loading data package for {name}...")
    sizes = load_seq_sizes(tsv)
    pairs = load_chrom_pairs(pairs_file)
    gaps  = load_gaps_bed(gaps_bed)
    sw    = load_switch_blocks(sw_bed)
    telo_coll, telo_nc = load_telomere_presence_tsv(telo_tsv)

    col_ribbons_list, col_spans, col_ins = parse_chain_detailed(col_chain, is_collinear=True)
    nc_ribbons_list, nc_spans, nc_ins   = parse_chain_detailed(nc_chain, is_collinear=False)

    col_ribbons_by_t = defaultdict(list)
    for r in col_ribbons_list:
        col_ribbons_by_t[r[0]].append(r)

    nc_ribbons_by_t = defaultdict(list)
    for r in nc_ribbons_list:
        nc_ribbons_by_t[r[0]].append(r)

    rec_ribbons_by_t = defaultdict(list)
    rec_spans = defaultdict(list)
    if rec_chain and os.path.exists(rec_chain):
        rec_ribbons_list, rec_spans, _ = parse_chain_detailed(rec_chain, is_collinear=True)
        for r in rec_ribbons_list:
            rec_ribbons_by_t[r[0]].append(r)

    all_ribbons_by_t = defaultdict(list)
    for r in (col_ribbons_list + nc_ribbons_list):
        all_ribbons_by_t[r[0]].append(r)

    return {
        'name': name,
        'sizes': sizes,
        'pairs': pairs,
        'gaps': gaps,
        'sw': sw,
        'telo_coll': telo_coll,
        'telo_nc': telo_nc,
        'col_spans': col_spans,
        'nc_spans': nc_spans,
        'rec_spans': rec_spans,
        'insertions': col_ins,
        'col_ribbons': col_ribbons_by_t,
        'nc_ribbons': nc_ribbons_by_t,
        'rec_ribbons': rec_ribbons_by_t,
        'ribbons': all_ribbons_by_t
    }


def main():
    args = parse_args()

    print("Loading T2T sequence sizes, centromeres, and telomeres...")
    t2t_sizes    = load_seq_sizes(args.t2t_tsv)
    centromeres  = load_centromeres(args.centromeres)
    p_starts     = load_p_arm_bed(args.telo_p_bed)
    flip_set     = build_flip_set_from_centromeres_and_telomeres(t2t_sizes, centromeres, p_starts)
    t2t_telo     = load_t2t_terminal_telomeres(args.telo_p_bed)

    print("Loading coverage summaries...")
    single_cov_sum = load_coverage_summary(args.single_cov_sum)
    dual_cov_sum   = load_coverage_summary(args.dual_cov_sum)
    ont_cov_sum    = load_coverage_summary(args.ont_cov_sum)

    # Load 3 assemblies
    single_data = load_assembly_package('HiFi Single', args.single_tsv, args.single_pairs,
                                        args.single_chain, args.single_nc_chain,
                                        args.single_gaps, args.single_bed, args.single_telomeres,
                                        rec_chain=args.single_rec_chain)

    dual_data   = load_assembly_package('HiFi Dual', args.dual_tsv, args.dual_pairs,
                                        args.dual_chain, args.dual_nc_chain,
                                        args.dual_gaps, args.dual_bed, args.dual_telomeres,
                                        rec_chain=args.dual_rec_chain)

    ont_data    = load_assembly_package('ONT Dual', args.ont_tsv, args.ont_pairs,
                                        args.ont_chain, args.ont_nc_chain,
                                        args.ont_gaps, args.ont_bed, args.ont_telomeres,
                                        rec_chain=args.ont_rec_chain)

    print("Rendering publication-quality butterfly synteny ideogram...")
    plot_butterfly_macro(
        t2t_sizes=t2t_sizes,
        flip_set=flip_set,
        single_data=single_data,
        dual_data=dual_data,
        ont_data=ont_data,
        centromeres=centromeres,
        t2t_telo=t2t_telo,
        single_cov_sum=single_cov_sum,
        dual_cov_sum=dual_cov_sum,
        ont_cov_sum=ont_cov_sum,
        out_png=args.output_png,
        out_pdf=args.output_pdf,
        dpi=args.dpi
    )


if __name__ == '__main__':
    main()
