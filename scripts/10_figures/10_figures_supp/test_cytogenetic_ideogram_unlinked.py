#!/usr/bin/env python3
"""
test_cytogenetic_ideogram_unlinked.py

Benchmark prototype for:
  1. Classical cytogenetic ideogram styling (hourglass centromeres, rounded telomere caps,
     dark/light coverage bands on white background).
  2. Chain gap / insertion detection (dq >= 20 kb) splitting collinear ribbons and
     leaving insertions unaligned/white on the assembly ideogram.
  3. Appending unlinked scaffolds (>= 20 kb) at chromosome ends with non-collinear synteny ribbons.
"""

import os
import re
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.path as mpath
from matplotlib.patches import PathPatch
from matplotlib.path import Path

# Paths
BASE = "/lustre/fs5/vgl/scratch/eduarte/curations/birds/bTaeGut7/redo_curation_analysis"
T2T_TSV = f"{BASE}/results/dual/03_coverage_dual/chain_cov.target_cov.tsv"
DUAL_TSV = f"{BASE}/results/dual/03_coverage_dual/chain_cov.query_cov.tsv"
DUAL_PAIRS = f"{BASE}/results/dual/02_chainpipeline_dual/t2t.vs.dual.best_chrom_pairs.tsv"
DUAL_COLL = f"{BASE}/results/dual/02_chainpipeline_dual/t2t.vs.dual.T2T.vs.ASM.target.collinear.chain"
DUAL_NC = f"{BASE}/results/dual/02_chainpipeline_dual/t2t.vs.dual.T2T.vs.ASM.target.non-collinear.chain"
DUAL_GAPS = f"{BASE}/results/dual/06_gaps_dual/dual_combined.renamed.sorted.reoriented.annotated.gaps.bed"
DUAL_SW = f"{BASE}/results/dual/05_hapmers_dual/dual.switch_blocks.final.bed"
DUAL_TELO = f"{BASE}/results/dual/04_telomeres_dual/dual_telomere_presence.tsv"
CENTROMERES_GFF = f"{BASE}/data/t2t/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff"
TELO_P_BED = f"{BASE}/data/t2t/bTaeGut7.T2T.fasta_terminal_telomeres.bed"

OUT_PNG = f"{BASE}/results/figures/10_figures_supp/test_chrZ_cytogenetic.png"

# Palette
COL_BORDER     = '#2c3e50'
COL_T2T_FILL   = '#F8FAFC'
COL_COV_D_DARK = '#6445B0'   # Dual collinear coverage
COL_COV_D_LIGHT= '#C4B8E8'   # Dual non-collinear coverage
COL_UNALIGNED  = '#FFFFFF'
COL_RIB_COLL   = '#38BDF8'
COL_RIB_COLL_E = '#0284C7'
COL_RIB_NC     = '#F59E0B'
COL_RIB_NC_E   = '#B45309'
COL_CUR_GAP    = '#E63946'
COL_ASM_GAP    = '#111111'
COL_SWITCH     = '#EF4444'
COL_TELO       = '#C4426A'


# -----------------------------------------------------------------------------
# Path Helpers
# -----------------------------------------------------------------------------

def _rect_path(x0, x1, yc, h):
    hb = h / 2
    verts = [(x0, yc - hb), (x1, yc - hb), (x1, yc + hb), (x0, yc + hb), (x0, yc - hb)]
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    return Path(verts, codes)


def _constriction_path(x0, x1, yc, h, cs, ce):
    """Draws an hourglass centromere constriction waist."""
    hb = h / 2
    wh = hb * 0.40  # inward pinch
    cx = (cs + ce) / 2
    w_cen = max(ce - cs, (x1 - x0) * 0.02)
    s_cen = max(x0, cx - w_cen / 2)
    e_cen = min(x1, cx + w_cen / 2)

    verts = [
        (x0, yc - hb),
        (s_cen, yc - hb),
        (cx, yc - hb + wh),
        (e_cen, yc - hb),
        (x1, yc - hb),
        (x1, yc + hb),
        (e_cen, yc + hb),
        (cx, yc + hb - wh),
        (s_cen, yc + hb),
        (x0, yc + hb),
        (x0, yc - hb)
    ]
    codes = [
        Path.MOVETO,
        Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
        Path.LINETO,
        Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
        Path.CLOSEPOLY
    ]
    return Path(verts, codes)


def _telo_semi(ax, x_base, side, yc, h, rx, hollow=False):
    """Draws a semicircular telomere cap."""
    hb = h / 2
    theta = np.linspace(np.pi / 2, 3 * np.pi / 2, 40) if side == 'p' else np.linspace(-np.pi / 2, np.pi / 2, 40)
    xs = x_base + rx * np.cos(theta)
    ys = yc + hb * np.sin(theta)
    fc = 'none' if hollow else COL_TELO
    ax.add_patch(mpatches.Polygon(list(zip(xs, ys)), closed=True,
                                  fc=fc, ec=COL_TELO, lw=0.8, zorder=10, clip_on=False))


def _merge_intervals(intervals):
    if not intervals: return []
    ivs = sorted(intervals)
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


# -----------------------------------------------------------------------------
# Chain parser with dq >= 20 kb gap splitting
# -----------------------------------------------------------------------------

def parse_chain_detailed(chain_path, is_collinear=True, insertion_threshold=20_000, min_block_size=5_000):
    ribbons = []
    aligned_query = defaultdict(list)
    insertions = defaultdict(list)

    if not chain_path or not os.path.exists(chain_path):
        return ribbons, aligned_query, insertions

    def _proc(hdr, blks):
        tName, qName = hdr['tName'], hdr['qName']
        qSize = hdr['qSize']
        qStrand = hdr['qStrand']
        t_pos = hdr['tStart']
        q_pos = hdr['qStart']

        cur_ts, cur_te = None, None
        cur_qs, cur_qe = None, None

        for blk in blks:
            sz = blk[0]
            dt = blk[1] if len(blk) > 1 else 0
            dq = blk[2] if len(blk) > 2 else 0

            blk_ts = t_pos
            blk_te = t_pos + sz
            if qStrand == '-':
                blk_qs = qSize - (q_pos + sz)
                blk_qe = qSize - q_pos
            else:
                blk_qs = q_pos
                blk_qe = q_pos + sz

            # Record aligned query span
            aligned_query[qName].append((min(blk_qs, blk_qe), max(blk_qs, blk_qe)))

            # Check if query gap is huge (>= 20 kb)
            if dq >= insertion_threshold:
                # Query insertion
                if qStrand == '-':
                    ins_s = qSize - (q_pos + sz + dq)
                    ins_e = qSize - (q_pos + sz)
                else:
                    ins_s = q_pos + sz
                    ins_e = q_pos + sz + dq
                insertions[qName].append((min(ins_s, ins_e), max(ins_s, ins_e)))

                # Close current ribbon
                if cur_ts is not None:
                    if (cur_te - cur_ts >= min_block_size) and (abs(cur_qe - cur_qs) >= min_block_size):
                        ribbons.append((tName, cur_ts, cur_te, qName,
                                        min(cur_qs, cur_qe), max(cur_qs, cur_qe),
                                        qStrand, is_collinear))
                    cur_ts, cur_te = None, None
                    cur_qs, cur_qe = None, None

                t_pos += sz + dt
                q_pos += sz + dq
                continue

            # Standard ribbon extension
            if cur_ts is None:
                cur_ts, cur_te = blk_ts, blk_te
                cur_qs, cur_qe = blk_qs, blk_qe
            else:
                if dt <= 25000 and dq <= 25000:
                    cur_te = blk_te
                    if qStrand == '-':
                        cur_qs = min(cur_qs, blk_qs)
                    else:
                        cur_qe = max(cur_qe, blk_qe)
                else:
                    if (cur_te - cur_ts >= min_block_size) and (abs(cur_qe - cur_qs) >= min_block_size):
                        ribbons.append((tName, cur_ts, cur_te, qName,
                                        min(cur_qs, cur_qe), max(cur_qs, cur_qe),
                                        qStrand, is_collinear))
                    cur_ts, cur_te = blk_ts, blk_te
                    cur_qs, cur_qe = blk_qs, blk_qe

            t_pos += sz + dt
            q_pos += sz + dq

        if cur_ts is not None and (cur_te - cur_ts >= min_block_size) and (abs(cur_qe - cur_qs) >= min_block_size):
            ribbons.append((tName, cur_ts, cur_te, qName,
                            min(cur_qs, cur_qe), max(cur_qs, cur_qe),
                            qStrand, is_collinear))

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

    return ribbons, aligned_query, insertions


def draw_bezier_ribbon(ax, t0, t1, yt, q0, q1, yq, col_fill, col_edge, alpha=0.35, inverted=False):
    if t1 <= t0 or q1 <= q0: return
    dy = yt - yq
    cy_top = yt - dy * 0.45
    cy_bot = yq + dy * 0.45
    if not inverted:
        pts = [(t0, yt), (t1, yt), (t1, cy_top), (q1, cy_bot), (q1, yq),
               (q0, yq), (q0, cy_bot), (t0, cy_top), (t0, yt)]
    else:
        pts = [(t0, yt), (t1, yt), (t1, cy_top), (q0, cy_bot), (q0, yq),
               (q1, yq), (q1, cy_bot), (t0, cy_top), (t0, yt)]
    codes = [Path.MOVETO, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    ax.add_patch(mpatches.PathPatch(Path(pts, codes), facecolor=col_fill,
                                    edgecolor=col_edge, alpha=alpha, lw=0.4, zorder=3))


# -----------------------------------------------------------------------------
# Main Test
# -----------------------------------------------------------------------------

def main():
    print("Testing cytogenetic ideogram with unlinked scaffolds and chain gap splitting on Chr Z...")
    # Load sizes
    t2t_sizes = {}
    with open(T2T_TSV) as fh:
        for line in fh:
            if line.startswith('#'): continue
            p = line.strip().split()
            if len(p) >= 2 and p[0] != 'name' and p[1] != 'size': t2t_sizes[p[0]] = int(p[1])

    dual_sizes = {}
    with open(DUAL_TSV) as fh:
        for line in fh:
            if line.startswith('#'): continue
            p = line.strip().split()
            if len(p) >= 2 and p[0] != 'name' and p[1] != 'size': dual_sizes[p[0]] = int(p[1])

    # Parse collinear and non-collinear chains
    col_ribs, col_q_spans, col_insertions = parse_chain_detailed(DUAL_COLL, is_collinear=True)
    nc_ribs, nc_q_spans, nc_insertions = parse_chain_detailed(DUAL_NC, is_collinear=False)

    t2t_z = "Mat_NC_133063.1_chromosome_Z"
    t2t_len = t2t_sizes[t2t_z]
    q_primary = "Pat_SUPER_Z.H1"
    q_prim_len = dual_sizes[q_primary]

    # Find unlinked scaffolds aligning to Chr Z with >= 20 kb
    all_z_ribs = [r for r in (col_ribs + nc_ribs) if r[0] == t2t_z]
    unloc_scaffs = {}
    for r in all_z_ribs:
        qn = r[3]
        if qn != q_primary:
            al_len = abs(r[5] - r[4])
            unloc_scaffs[qn] = unloc_scaffs.get(qn, 0) + al_len

    unloc_qual = [k for k, v in unloc_scaffs.items() if v >= 20_000]
    print(f"Qualified unlinked scaffolds for Chr Z: {unloc_qual} ({unloc_scaffs})")

    # Layout figure
    fig, ax = plt.subplots(figsize=(16, 6), facecolor='white')
    ax.set_facecolor('white')

    # Coordinates
    SPACER = 1_500_000 # 1.5 Mb spacer between primary and unlinked
    offsets = {q_primary: 0}
    cur_x = q_prim_len
    for u in unloc_qual:
        cur_x += SPACER
        offsets[u] = cur_x
        cur_x += dual_sizes.get(u, 500_000)

    max_x = max(t2t_len, cur_x) * 1.03
    ax.set_xlim(-max_x * 0.02, max_x)
    ax.set_ylim(-0.5, 3.0)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_yticks([])

    # Vertical levels
    BAR_H = 0.35
    y_t2t = 2.0
    y_asm = 0.6

    # 1. T2T Reference Bar
    # Centromere
    cs_z, ce_z = 41_000_000, 43_000_000 # rough for test
    t2t_path = _constriction_path(0, t2t_len, y_t2t, BAR_H, cs_z, ce_z)
    ax.add_patch(PathPatch(t2t_path, fc=COL_T2T_FILL, ec=COL_BORDER, lw=0.9, zorder=4))
    _telo_semi(ax, 0, 'p', y_t2t, BAR_H, max_x * 0.005)
    _telo_semi(ax, t2t_len, 'q', y_t2t, BAR_H, max_x * 0.005)
    ax.text(-max_x * 0.015, y_t2t, 'T2T Ref', ha='right', va='center', fontweight='bold', fontsize=10)

    # 2. Assembly Primary Bar
    prim_path = _rect_path(0, q_prim_len, y_asm, BAR_H)
    ax.add_patch(PathPatch(prim_path, fc=COL_UNALIGNED, ec=COL_BORDER, lw=0.9, zorder=4))

    # Bands inside Primary bar
    # Collinear bands
    m_col = _merge_intervals(col_q_spans[q_primary])
    for s, e in m_col:
        r = mpatches.Rectangle((s, y_asm - BAR_H/2), e - s, BAR_H, fc=COL_COV_D_DARK, ec='none', zorder=5)
        r.set_clip_path(prim_path, transform=ax.transData)
        ax.add_patch(r)

    # Non-collinear bands
    m_nc = _merge_intervals(nc_q_spans[q_primary])
    for s, e in m_nc:
        r = mpatches.Rectangle((s, y_asm - BAR_H/2), e - s, BAR_H, fc=COL_COV_D_LIGHT, ec='none', zorder=5)
        r.set_clip_path(prim_path, transform=ax.transData)
        ax.add_patch(r)

    _telo_semi(ax, 0, 'p', y_asm, BAR_H, max_x * 0.005)
    _telo_semi(ax, q_prim_len, 'q', y_asm, BAR_H, max_x * 0.005)
    ax.text(-max_x * 0.015, y_asm, 'HiFi Dual (Chr Z)', ha='right', va='center', fontweight='bold', fontsize=10)

    # 3. Unlinked Scaffolds Bars
    for u in unloc_qual:
        u_off = offsets[u]
        u_len = dual_sizes[u]
        u_path = _rect_path(u_off, u_off + u_len, y_asm, BAR_H)
        ax.add_patch(PathPatch(u_path, fc=COL_UNALIGNED, ec='#94A3B8', lw=0.8, ls='--', zorder=4))

        # Non-collinear bands on unlinked
        u_nc = _merge_intervals(nc_q_spans[u])
        for s, e in u_nc:
            r = mpatches.Rectangle((u_off + s, y_asm - BAR_H/2), e - s, BAR_H, fc=COL_COV_D_LIGHT, ec='none', zorder=5)
            r.set_clip_path(u_path, transform=ax.transData)
            ax.add_patch(r)

        # Label above unlinked
        lbl = u.replace('Pat_SUPER_', '').replace('.H1', '')
        ax.text(u_off + u_len/2, y_asm + BAR_H*0.7, lbl, ha='center', va='bottom', fontsize=7.5, color='#475569')

    # 4. Ribbons
    for r in all_z_ribs:
        tn, ts, te, qn, qs, qe, qstr, is_col = r
        if qn not in offsets: continue
        q_off = offsets[qn]
        q0 = q_off + qs
        q1 = q_off + qe
        c_fill = COL_RIB_COLL if is_col else COL_RIB_NC
        c_edge = COL_RIB_COLL_E if is_col else COL_RIB_NC_E
        draw_bezier_ribbon(ax, ts, te, y_t2t - BAR_H/2, q0, q1, y_asm + BAR_H/2,
                           c_fill, c_edge, alpha=0.35 if is_col else 0.55, inverted=(qstr == '-'))

    # Scale Bar
    ax.plot([0, 20_000_000], [0.1, 0.1], color='#0F172A', lw=2)
    ax.text(10_000_000, -0.05, '20 Mb', ha='center', va='top', fontsize=9, fontweight='bold')

    plt.tight_layout()
    fig.savefig(OUT_PNG, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved test figure to {OUT_PNG}")

if __name__ == '__main__':
    main()
