#!/usr/bin/env python3
"""
plot_ideogram_haplotype_avgcov_v6.py

Five-column butterfly ideogram:
  Col 0: Maternal ideogram — HiFi (mirrored butterfly; p-arm toward center)
  Col 1: Paternal ideogram — HiFi
  Col 2: Averaged coverage bar (Single + Dual + ONT Dual, averaged over mat+pat)
  Col 3: Maternal ideogram — ONT Dual (mirrored butterfly)
  Col 4: Paternal ideogram — ONT Dual

Each chromosome slot displays three bar positions:
  top bar    — Single assembly (teal)
  middle bar — Dual assembly   (purple)
  bottom bar — ONT Dual        (green)

Features:
- Collinear and non-collinear chain coverage
- Switch-error blocks (hapmers)
- Assembly vs Curation gaps
- Telomere presence markers (p-arm and q-arm, solid for collinear, hollow for non-collinear)
- Centromere constriction waist
- P-arm-first reorientation relative to T2T reference
"""

import os
import re
import argparse
from datetime import datetime
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.path as mpath
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import blended_transform_factory
from matplotlib.legend_handler import HandlerPatch
from matplotlib.lines import Line2D


# ---- Palette ----------------------------------------------------------------
COL_COV_S_DARK = '#00BCCC'
COL_COV_D_DARK = '#6445B0'
COL_COV_O_DARK = '#2E8B52'
COL_COV_S      = '#A0D4DC'
COL_COV_D      = '#C4B8E8'
COL_COV_O      = '#B8E8CC'
COL_UNCOV      = '#ffffff'
COL_BORDER     = '#2c3e50'
COL_ROW_ODD    = '#f2f9fd'
COL_SW         = '#C49060'
COL_TELO       = '#C4426A'
COL_PAT        = '#B85C3E'
COL_MAT        = '#3E6EB8'
COL_CEN        = '#D62728'
COL_GAP        = '#E63946'   # curation gap markers
COL_ASM_GAP    = '#111111'   # assembly gap tick marks


# ---- Style presets ----------------------------------------------------------
STYLES = {
    'paper': dict(
        fig_w=22.0, fig_h=14.0, BAR_H=0.48,
        border_lw=0.4,  telo_lw_h=0.5,  telo_lw_f=0.3,
        font_base=7,   font_chrom=9,  font_pm=6,
        font_tick=6,   font_xlabel=7,  font_title=9,
        font_n=6,      font_legend=7,  ncol_legend=6,
        dpi_default=300,
    ),
    'poster': dict(
        fig_w=30.0, fig_h=32.0, BAR_H=0.85,
        border_lw=0.9,  telo_lw_h=1.4,  telo_lw_f=0.8,
        font_base=16,  font_chrom=22, font_pm=14,
        font_tick=16,  font_xlabel=18, font_title=22,
        font_n=14,     font_legend=17, ncol_legend=6,
        dpi_default=150,
    ),
}


class _SemiHandler(HandlerPatch):
    def __init__(self, hollow=False, **kw):
        self._hollow = hollow
        super().__init__(**kw)

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        r  = min(width, height) * 0.46
        cx = width * 0.5 - xdescent
        cy = height * 0.5 - ydescent
        theta = np.linspace(np.pi / 2, 3 * np.pi / 2, 40)
        verts = list(zip(cx + r * np.cos(theta), cy + r * np.sin(theta)))
        p = mpatches.Polygon(verts, closed=True,
                             fc='none' if self._hollow else COL_TELO,
                             ec=COL_TELO, lw=0.9, transform=trans)
        return [p]


class _ConstrictionHandler(HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        w, h = width, height
        cx   = w / 2 - xdescent
        cy   = h / 2 - ydescent
        wh   = h * 0.13
        cs   = cx - w * 0.13
        ce   = cx + w * 0.13
        verts = np.array([
            [-xdescent,    -ydescent],
            [cs,           -ydescent],
            [cx,           cy - wh],
            [ce,           -ydescent],
            [w - xdescent, -ydescent],
            [w - xdescent, h - ydescent],
            [ce,           h - ydescent],
            [cx,           cy + wh],
            [cs,           h - ydescent],
            [-xdescent,    h - ydescent],
        ])
        codes = [mpath.Path.MOVETO] + [mpath.Path.LINETO] * 9
        p = PathPatch(mpath.Path(verts, np.array(codes), closed=True),
                      fc=orig_handle.get_facecolor(),
                      ec=orig_handle.get_edgecolor(),
                      lw=orig_handle.get_linewidth(),
                      transform=trans)
        return [p]


# =============================================================================
# Chain / interval utilities
# =============================================================================

def _merge(intervals):
    if not intervals:
        return []
    ivs = sorted(intervals)
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def _parse_chain_target(path):
    if not path or not os.path.exists(path):
        return
    with open(path) as fh:
        header = None
        blocks = []
        for raw in fh:
            line = raw.rstrip('\n')
            if not line or line.startswith('#'):
                if header is not None and not line:
                    yield _emit(header, blocks)
                    header = None; blocks = []
                continue
            if line.startswith('chain'):
                if header is not None:
                    yield _emit(header, blocks)
                p = line.split()
                header = dict(tName=p[2], tSize=int(p[3]), tStart=int(p[5]))
                blocks = []
            else:
                parts = line.split()
                if parts:
                    blocks.append(tuple(int(x) for x in parts))
        if header is not None:
            yield _emit(header, blocks)


def _emit(header, blocks):
    pos = header['tStart']
    ivs = []
    for block in blocks:
        size = block[0]
        ivs.append((pos, pos + size))
        pos += size + (block[1] if len(block) == 3 else 0)
    return header['tName'], header['tSize'], ivs


def compute_covered(chain_path, seq_sizes):
    covered = defaultdict(list)
    for tname, _, ivs in _parse_chain_target(chain_path):
        covered[tname].extend(ivs)
    return {name: _merge(covered.get(name, [])) for name in seq_sizes}


def merge_coverage(cov_a, cov_b):
    all_names = set(cov_a) | set(cov_b)
    return {name: _merge(cov_a.get(name, []) + cov_b.get(name, []))
            for name in all_names}


# =============================================================================
# Loaders
# =============================================================================

def load_tsv(path):
    df = pd.read_csv(path, sep='\t', comment='#')
    df.columns = [c.strip() for c in df.columns]
    return df.set_index('name')


def load_fai(path):
    sizes = {}
    if not path or not os.path.exists(path):
        return sizes
    with open(path) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= 2:
                sizes[p[0]] = int(p[1])
    return sizes


def load_chrom_pairs(path):
    mapping = {}
    if not path or not os.path.exists(path):
        return mapping
    with open(path) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= 2:
                mapping[p[0]] = p[1]
    return mapping


def load_switch_blocks(bed_path):
    blocks = defaultdict(list)
    if not bed_path or not os.path.exists(bed_path):
        return dict(blocks)
    with open(bed_path) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= 4:
                blocks[p[0]].append((int(p[1]), int(p[2]), p[3]))
    return dict(blocks)


def load_telomere_presence_tsv(path):
    """
    Loads telomere presence from <asm>_telomere_presence.tsv:
      chromosome  scaffold  collinear  non-collinear  teloscope_rescued  type_of_tele_rescued
    Returns:
      (telo_coll, telo_nc) dicts of {t2t_chrom: set('p', 'q')}
    """
    telo_coll = defaultdict(set)
    telo_nc   = defaultdict(set)
    if not path or not os.path.exists(path):
        return dict(telo_coll), dict(telo_nc)
    with open(path) as f:
        header = None
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            p = line.split('\t')
            if header is None:
                header = [c.lower() for c in p]
                continue
            chrom = p[0]
            coll_val = p[2].lower() if len(p) > 2 else 'none'
            nc_val   = p[3].lower() if len(p) > 3 else 'none'
            for arm in ['p', 'q']:
                if arm in coll_val:
                    telo_coll[chrom].add(arm)
                elif arm in nc_val:
                    telo_nc[chrom].add(arm)
    return dict(telo_coll), dict(telo_nc)


def load_p_arm_bed(bed_path):
    """
    Reads p-arm telomere positions.
    Can accept a 2-column BED (chrom, start) or a 5-column BED where col 4/5 is 'p'.
    """
    p_starts = {}
    if not bed_path or not os.path.exists(bed_path):
        return p_starts
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            p = line.split('\t')
            if len(p) >= 5:
                if p[4].lower() == 'p' or p[3].lower() == 'p':
                    p_starts[p[0]] = int(p[1])
            elif len(p) >= 2:
                p_starts[p[0]] = int(p[1])
    return p_starts


def build_flip_set(p_arm_positions, seq_sizes):
    return {nm for nm, p_start in p_arm_positions.items()
            if nm in seq_sizes and p_start > seq_sizes[nm] / 2}


def load_centromeres(gff_path):
    centromeres = {}
    if not gff_path or not os.path.exists(gff_path):
        return centromeres
    with open(gff_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 5:
                continue
            m = re.match(r'^chr([0-9]+[A-Za-z]*|[A-Za-z]+)_(pat|mat)$', parts[0])
            if m:
                centromeres[(m.group(1), m.group(2))] = (int(parts[3]), int(parts[4]))
    return centromeres


def load_annotated_gaps_bed(path):
    """
    Reads annotated gaps BED (chrom, start, end, label) where label is ASSEMBLY or CURATION.
    Returns (assembly_raw, curation_raw) as {chrom: [(start, end), ...]}.
    """
    assembly = defaultdict(list)
    curation = defaultdict(list)
    if not path or not os.path.exists(path):
        return dict(assembly), dict(curation)
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            p = line.split('\t')
            if len(p) >= 4:
                s, e, label = int(p[1]), int(p[2]), p[3].upper()
                if 'CURATION' in label:
                    curation[p[0]].append((s, e))
                else:
                    assembly[p[0]].append((s, e))
            elif len(p) >= 3:
                assembly[p[0]].append((int(p[1]), int(p[2])))
    return dict(assembly), dict(curation)


def load_single_coverage_summary(path):
    """Loads per-chromosome coverage percentages from an assembly coverage summary TSV."""
    cov_dict = {}
    if not path or not os.path.exists(path):
        return cov_dict
    df = pd.read_csv(path, sep='\t', comment='#')
    for _, row in df.iterrows():
        cov_dict[row['name']] = (
            float(row['collinear_pct']),
            float(row['noncollinear_pct']),
            float(row['uncovered_pct'])
        )
    return cov_dict


def load_combined_coverage_summary(path):
    """Reads legacy combined coverage summary TSV with an 'assembly' column."""
    if not path or not os.path.exists(path):
        return {}, {}, {}
    df = pd.read_csv(path, sep='\t', comment='#')
    cov_s, cov_d, cov_ont = {}, {}, {}
    for _, row in df.iterrows():
        t = (float(row['collinear_pct']),
             float(row['noncollinear_pct']),
             float(row['uncovered_pct']))
        asm_name = str(row['assembly']).lower()
        if 'single' in asm_name:
            cov_s[row['name']] = t
        elif 'dual' in asm_name and 'ont' not in asm_name:
            cov_d[row['name']] = t
        elif 'ont' in asm_name:
            cov_ont[row['name']] = t
    return cov_s, cov_d, cov_ont


# =============================================================================
# Chain liftover for Gaps and Switch blocks
# =============================================================================

def build_chain_query_liftover(chain_path):
    liftover = defaultdict(lambda: defaultdict(list))
    if not chain_path or not os.path.exists(chain_path):
        return {}

    def _process(hdr, blks):
        t_pos = hdr['tStart']; q_pos = hdr['qStart']
        q_strand = hdr['qStrand']; q_size = hdr['qSize']
        for blk in blks:
            size = blk[0]; dt = blk[1] if len(blk) > 1 else 0; dq = blk[2] if len(blk) > 2 else 0
            qs = (q_size - (q_pos + size)) if q_strand == '-' else q_pos
            qe = (q_size - q_pos)          if q_strand == '-' else q_pos + size
            liftover[hdr['qName']][hdr['tName']].append((qs, qe, t_pos, q_strand))
            t_pos += size + dt; q_pos += size + dq

    with open(chain_path) as fh:
        header = None; blocks = []
        for raw in fh:
            line = raw.rstrip('\n')
            if not line or line.startswith('#'):
                if header is not None and not line:
                    _process(header, blocks); header = None; blocks = []
                continue
            if line.startswith('chain'):
                if header is not None:
                    _process(header, blocks)
                p = line.split()
                header = {'tName': p[2], 'tSize': int(p[3]), 'tStrand': p[4],
                          'tStart': int(p[5]), 'qName': p[7], 'qSize': int(p[8]),
                          'qStrand': p[9], 'qStart': int(p[10])}
                blocks = []
            else:
                parts = line.split()
                if parts:
                    blocks.append(tuple(int(x) for x in parts))
        if header is not None:
            _process(header, blocks)

    return {qn: {tn: sorted(blks, key=lambda x: x[0]) for tn, blks in td.items()}
            for qn, td in liftover.items()}


def _map_gap_pos(gs, col_blks, nc_blks):
    """
    Lifts a gap coordinate gs to T2T coordinates using collinear and non-collinear blocks.
    Prioritizes direct overlap with collinear blocks, then non-collinear blocks,
    and finally clamps to the nearest flanking block boundary.
    """
    # 1. Direct overlap with collinear block
    for blk in col_blks:
        qs, qe, ts = blk[0], blk[1], blk[2]
        strand = blk[3] if len(blk) > 3 else '+'
        if qs <= gs <= qe:
            return ts + (qe - gs if strand == '-' else gs - qs)
    # 2. Direct overlap with non-collinear block
    for blk in nc_blks:
        qs, qe, ts = blk[0], blk[1], blk[2]
        strand = blk[3] if len(blk) > 3 else '+'
        if qs <= gs <= qe:
            return ts + (qe - gs if strand == '-' else gs - qs)
    # 3. Flanking/nearest block across all blocks
    all_blks = col_blks + nc_blks
    if not all_blks:
        return None
    preceding = [b for b in all_blks if b[1] < gs]
    following = [b for b in all_blks if b[0] > gs]
    b_prec = max(preceding, key=lambda x: x[1]) if preceding else None
    b_foll = min(following, key=lambda x: x[0]) if following else None
    if b_prec and b_foll:
        d_p = gs - b_prec[1]
        d_f = b_foll[0] - gs
        chosen = b_prec if d_p <= d_f else b_foll
        chosen_edge = 'end' if d_p <= d_f else 'start'
    elif b_prec:
        chosen = b_prec; chosen_edge = 'end'
    elif b_foll:
        chosen = b_foll; chosen_edge = 'start'
    else:
        return None
    qs, qe, ts = chosen[0], chosen[1], chosen[2]
    strand = chosen[3] if len(chosen) > 3 else '+'
    if chosen_edge == 'end':
        return ts if strand == '-' else ts + (qe - qs)
    else:
        return ts + (qe - qs) if strand == '-' else ts


def separate_curation_gaps(positions, min_sep=350_000, chrom_size=None):
    """
    Separates curation gaps that map to identical or near-identical coordinates
    (e.g., due to falling within unaligned/introduced contigs that clamped to the
    same flanking chain block). Centers clusters around their anchor point and
    ensures adjacent gaps are separated by at least min_sep so introduced sequences
    are clearly distinguishable.
    """
    if not positions:
        return []
    sorted_items = sorted(positions, key=lambda x: (x[1], x[0]))

    clusters = []
    curr_cluster = [sorted_items[0]]
    for item in sorted_items[1:]:
        if item[1] - curr_cluster[-1][1] < min_sep:
            curr_cluster.append(item)
        else:
            clusters.append(curr_cluster)
            curr_cluster = [item]
    clusters.append(curr_cluster)

    adjusted = []
    for cl in clusters:
        k = len(cl)
        if k == 1:
            adjusted.append(cl[0][1])
        else:
            center = sum(x[1] for x in cl) / k
            span = (k - 1) * min_sep
            start = center - span / 2
            for i in range(k):
                adjusted.append(int(start + i * min_sep))

    for i in range(1, len(adjusted)):
        if adjusted[i] < adjusted[i-1] + min_sep:
            adjusted[i] = adjusted[i-1] + min_sep
    if chrom_size and adjusted[-1] > chrom_size:
        shift = adjusted[-1] - chrom_size
        adjusted = [max(0, p - shift) for p in adjusted]
    for i in range(len(adjusted)-2, -1, -1):
        if adjusted[i] > adjusted[i+1] - min_sep:
            adjusted[i] = max(0, adjusted[i+1] - min_sep)

    return adjusted


def lift_gaps_to_t2t(gaps_raw, pairs_path, chain_path, nc_chain_path=None, seq_sizes=None, is_curation=False):
    """
    Unified liftover for gaps from assembly coordinates to T2T coordinates,
    supporting both collinear and non-collinear chain alignments.
    When is_curation=True, separates clustered curation gaps so introduced sequences
    are clearly distinguishable.
    """
    if not gaps_raw or not pairs_path or not chain_path:
        return {}
    pairs = load_chrom_pairs(pairs_path)
    super_to_t2t = {v: k for k, v in pairs.items()}
    liftover_col = build_chain_query_liftover(chain_path)
    liftover_nc  = build_chain_query_liftover(nc_chain_path) if nc_chain_path else {}

    result = defaultdict(list)
    for super_name, intervals in gaps_raw.items():
        t2t_name = super_to_t2t.get(super_name)
        if not t2t_name:
            continue
        col_blocks = liftover_col.get(super_name, {}).get(t2t_name, [])
        nc_blocks  = liftover_nc.get(super_name, {}).get(t2t_name, [])
        if not col_blocks and not nc_blocks:
            continue
        mapped = []
        for gs, ge in intervals:
            pos = _map_gap_pos(gs, col_blocks, nc_blocks)
            if pos is not None:
                mapped.append((gs, pos))
        if not mapped:
            continue
        if is_curation:
            sz = seq_sizes.get(t2t_name) if seq_sizes else None
            min_sep = max(int(sz * 0.005), 80_000) if sz else 350_000
            sep_positions = separate_curation_gaps(mapped, min_sep=min_sep, chrom_size=sz)
            for pos in sep_positions:
                result[t2t_name].append((pos - 50, pos + 50))
        else:
            for _, pos in mapped:
                result[t2t_name].append((pos - 50, pos + 50))
    return dict(result)


def extract_insertions_and_gaps(asm_gaps_raw, cur_gaps_raw, pairs_path, chain_path, nc_chain_path=None, seq_sizes=None, min_ins_size=30_000):
    """
    Extracts query insertions (dq >= min_ins_size) from collinear and non-collinear chains,
    and classifies curation and assembly gaps into:
      1. Aligned gaps (falling in aligned chain blocks) -> lifted to T2T coordinates.
      2. Insertion gaps (falling inside query insertions) -> stored within their insertion event.
    """
    if not pairs_path or not chain_path:
        return {}, {}, {}
    pairs = load_chrom_pairs(pairs_path)
    super_to_t2t = {v: k for k, v in pairs.items()}

    def parse_chain_detailed(cp):
        chains = {}
        if not cp: return chains
        with open(cp) as f:
            hdr = None; blks = []
            for line in f:
                line = line.strip()
                if line.startswith('chain'):
                    if hdr: chains.setdefault((hdr['qName'], hdr['tName']), []).append((hdr, blks))
                    p = line.split()
                    hdr = {'tName': p[2], 'tStart': int(p[5]), 'tEnd': int(p[6]),
                           'qName': p[7], 'qStart': int(p[10]), 'qEnd': int(p[11]),
                           'strand': p[9]}
                    blks = []
                elif hdr and line and not line.startswith('#'):
                    p = line.split()
                    blks.append((int(p[0]), int(p[1]) if len(p)>1 else 0, int(p[2]) if len(p)>2 else 0))
            if hdr: chains.setdefault((hdr['qName'], hdr['tName']), []).append((hdr, blks))
        return chains

    chains_col = parse_chain_detailed(chain_path)
    chains_nc  = parse_chain_detailed(nc_chain_path) if nc_chain_path else {}
    liftover_col = build_chain_query_liftover(chain_path)
    liftover_nc  = build_chain_query_liftover(nc_chain_path) if nc_chain_path else {}

    aligned_cur = {}
    aligned_asm = {}
    insertions_map = {}

    all_q_names = set(list((cur_gaps_raw or {}).keys()) + list((asm_gaps_raw or {}).keys()))
    for super_name in all_q_names:
        t2t_name = super_to_t2t.get(super_name)
        if not t2t_name: continue
        all_chains = chains_col.get((super_name, t2t_name), []) + chains_nc.get((super_name, t2t_name), [])
        
        # 1. Discover all query insertions >= min_ins_size
        ins_list = []
        for h, blks in all_chains:
            t_curr = h['tStart']; q_curr = h['qStart']
            for size, dt, dq in blks:
                if dq >= min_ins_size:
                    t_anchor = t_curr + size if h['strand'] == '+' else t_curr + dt
                    ins_list.append({
                        't_anchor': t_anchor,
                        'dt': dt,
                        'dq': dq,
                        'q_start': q_curr + size,
                        'q_end': q_curr + size + dq,
                        'strand': h['strand'],
                        'cur_gaps': [],
                        'asm_gaps': []
                    })
                t_curr += size + dt
                q_curr += size + dq
                
        col_blks = liftover_col.get(super_name, {}).get(t2t_name, [])
        nc_blks  = liftover_nc.get(super_name, {}).get(t2t_name, [])
        
        # 2. Process curation gaps
        mapped_cur = []
        for gs, ge in (cur_gaps_raw or {}).get(super_name, []):
            mid = (gs + ge) / 2
            inside_ins = None
            for ins in ins_list:
                if ins['q_start'] <= mid <= ins['q_end']:
                    inside_ins = ins
                    break
            if inside_ins:
                inside_ins['cur_gaps'].append((mid - inside_ins['q_start'], gs, ge))
            else:
                pos = _map_gap_pos(gs, col_blks, nc_blks)
                if pos is not None:
                    mapped_cur.append((gs, pos))
        if mapped_cur:
            sz = seq_sizes.get(t2t_name) if seq_sizes else None
            min_sep = max(int(sz * 0.005), 80_000) if sz else 350_000
            sep_positions = separate_curation_gaps(mapped_cur, min_sep=min_sep, chrom_size=sz)
            aligned_cur[t2t_name] = [(p - 50, p + 50) for p in sep_positions]
            
        # 3. Process assembly gaps
        for gs, ge in (asm_gaps_raw or {}).get(super_name, []):
            mid = (gs + ge) / 2
            inside_ins = None
            for ins in ins_list:
                if ins['q_start'] <= mid <= ins['q_end']:
                    inside_ins = ins
                    break
            if inside_ins:
                inside_ins['asm_gaps'].append((mid - inside_ins['q_start'], gs, ge))
            else:
                pos = _map_gap_pos(gs, col_blks, nc_blks)
                if pos is not None:
                    aligned_asm.setdefault(t2t_name, []).append((pos - 50, pos + 50))

        if ins_list:
            insertions_map[t2t_name] = ins_list

    return aligned_asm, aligned_cur, insertions_map


_LIFTOVER_GAP         = 50_000
_CHAIN_EDGE_TOLERANCE = 10

def liftover_blocks(raw_blocks, chain_liftover, scaffold_name, t2t_name):
    chain_blocks = chain_liftover.get(scaffold_name, {}).get(t2t_name, [])
    if not chain_blocks:
        return []
    result = []
    for h_start, h_end, _ in raw_blocks:
        overlapping = [blk for blk in chain_blocks
                       if max(h_start, blk[0]) < min(h_end, blk[1])]
        if not overlapping:
            continue
        first_qs = overlapping[0][0];  last_qe = overlapping[-1][1]
        if h_start < first_qs - _CHAIN_EDGE_TOLERANCE or h_end > last_qe + _CHAIN_EDGE_TOLERANCE:
            continue
        ivs = []
        for blk in overlapping:
            qs, qe, ts = blk[0], blk[1], blk[2]
            strand = blk[3] if len(blk) > 3 else '+'
            ov_s = max(h_start, qs);  ov_e = min(h_end, qe)
            if strand == '-':
                p_s  = ts + (qe - ov_e);  p_e  = ts + (qe - ov_s)
            else:
                p_s  = ts + (ov_s - qs);  p_e  = ts + (ov_e - qs)
            if p_s < p_e:
                ivs.append((p_s, p_e))
        if not ivs:
            continue
        merged_s, merged_e = ivs[0]
        for p_s, p_e in ivs[1:]:
            if p_s - merged_e <= _LIFTOVER_GAP:
                merged_e = max(merged_e, p_e)
            else:
                result.append((merged_s, merged_e, ''))
                merged_s, merged_e = p_s, p_e
        result.append((merged_s, merged_e, ''))
    return sorted(result, key=lambda x: x[0])


def build_switch_lookup(pairs_path, bed_path, chain_path):
    if not pairs_path or not bed_path or not chain_path:
        return {}
    pairs = load_chrom_pairs(pairs_path)
    switches = load_switch_blocks(bed_path)
    liftover = build_chain_query_liftover(chain_path)
    result = {}
    for t2t_name, super_full in pairs.items():
        raw = switches.get(super_full, [])
        result[t2t_name] = liftover_blocks(raw, liftover, super_full, t2t_name)
    return result


# =============================================================================
# Sort / classification helpers
# =============================================================================

def chrom_token(name):
    m = re.search(r'(?:chromosome_|SUPER_)([0-9]+[A-Za-z]*|[A-Za-z]+)', name)
    return m.group(1) if m else name


def chrom_sort_key(token):
    SEX = {'Z': 900, 'W': 901}
    m = re.match(r'^(\d+)([A-Za-z]*)$', token)
    if m:
        return (0, int(m.group(1)), m.group(2).upper())
    return (1, SEX.get(token.upper(), 999), token.upper())


def is_pat(name):
    tok = chrom_token(name)
    if tok.upper() == 'Z': return True
    if tok.upper() == 'W': return False
    return name.startswith('Pat') or '.pat' in name.lower()


DOT_CHROMS = {'16', '25', '29', '30', '31', '32', '33', '34', '35', '36', '37'}

def classify_chrom_group(token):
    if token.upper() in ('Z', 'W'): return 'macro'
    if token in DOT_CHROMS:         return 'nano'
    m = re.match(r'^(\d+)', token)
    if m:
        n = int(m.group(1))
        if n <= 8:    return 'macro'
        elif n <= 28: return 'micro'
        else:         return 'nano'
    return 'macro'


# =============================================================================
# Geometry
# =============================================================================

def _rect_path(x0, x1, yc, bar_h):
    h = bar_h / 2
    xs = np.array([x0, x1, x1, x0])
    ys = np.array([yc - h, yc - h, yc + h, yc + h])
    verts = np.column_stack([xs, ys])
    codes = np.array([mpath.Path.MOVETO, mpath.Path.LINETO,
                      mpath.Path.LINETO, mpath.Path.LINETO])
    return mpath.Path(verts, codes, closed=True)


def _constriction_path(x0, x1, yc, bar_h, cs, ce, waist_frac=0.55):
    h  = bar_h / 2
    wh = bar_h * waist_frac / 2
    cx = (cs + ce) / 2
    verts = np.array([
        [x0,  yc - h],
        [cs,  yc - h],
        [cx,  yc - wh],
        [ce,  yc - h],
        [x1,  yc - h],
        [x1,  yc + h],
        [ce,  yc + h],
        [cx,  yc + wh],
        [cs,  yc + h],
        [x0,  yc + h],
    ])
    codes = np.array([mpath.Path.MOVETO] + [mpath.Path.LINETO] * 9)
    return mpath.Path(verts, codes, closed=True)


# =============================================================================
# Y-layout
# =============================================================================

def make_y_layout_triple(group_order, chrom_groups, hap_filter,
                         BAR_H, HAP_GAP, GROUP_GAP,
                         SUB_H=0.0, SUB_GAP=0.0):
    """Three bar positions per chromosome: Single (top), Dual (middle), ONT (bottom)."""
    h = BAR_H / 2
    yp_s    = {}
    yp_d    = {}
    yp_o    = {}
    y_sub_s = {}
    y_sub_d = {}
    y_sub_o = {}
    y_lbl   = {}
    y_spn   = {}
    y = 0.0
    for tok in group_order:
        nm = next((n for n in chrom_groups.get(tok, []) if hap_filter(n)), None)
        if nm is None:
            continue
        span_top   = y
        yp_s[nm]   = y + h;  y += BAR_H
        if SUB_H > 0:
            y += SUB_GAP
            y_sub_s[nm] = y + SUB_H / 2; y += SUB_H
        y          += HAP_GAP
        yp_d[nm]   = y + h;  y += BAR_H
        if SUB_H > 0:
            y += SUB_GAP
            y_sub_d[nm] = y + SUB_H / 2; y += SUB_H
        y          += HAP_GAP
        yp_o[nm]   = y + h;  y += BAR_H
        if SUB_H > 0:
            y += SUB_GAP
            y_sub_o[nm] = y + SUB_H / 2; y += SUB_H
        span_bot   = y
        y_lbl[tok] = (span_top + span_bot) / 2
        y_spn[tok] = (span_top, span_bot)
        y          += GROUP_GAP
    total_y = max(y - GROUP_GAP, BAR_H)
    return yp_s, yp_d, yp_o, y_sub_s, y_sub_d, y_sub_o, y_lbl, y_spn, total_y


# =============================================================================
# Averaged coverage computation
# =============================================================================

def compute_avg_coverage_by_method(group_order, chrom_groups,
                                   cov_summary_s, cov_summary_d,
                                   cov_summary_ont=None):
    avg_s, avg_d, avg_ont = {}, {}, {}
    for tok in group_order:
        for cov_in, cov_out in [(cov_summary_s, avg_s),
                                (cov_summary_d, avg_d),
                                (cov_summary_ont, avg_ont)]:
            if not cov_in:
                continue
            vals = [cov_in[nm] for nm in chrom_groups.get(tok, [])
                    if nm in cov_in]
            if vals:
                c  = sum(v[0] for v in vals) / len(vals)
                nc = sum(v[1] for v in vals) / len(vals)
                cov_out[tok] = (c, nc, max(0.0, 100.0 - c - nc))
    return avg_s, avg_d, avg_ont


# =============================================================================
# Drawing functions
# =============================================================================

def _draw_avg_cov_panel(ax_b, hap_names, yp_single, yp_dual, yp_ont, y_spans,
                        group_order, total_y, margin_y, BAR_H, GROUP_GAP,
                        avg_s, avg_d, avg_ont=None,
                        style=None, show_xlabel=False, rasterized=False):
    st = style or STYLES['paper']
    hb = BAR_H / 2

    ax_b.set_xlim(-2, 102)
    ax_b.set_ylim(total_y + margin_y, -margin_y)
    ax_b.set_facecolor('white')

    for i, tok in enumerate(group_order):
        if i % 2 == 0 and tok in y_spans:
            y0, y1 = y_spans[tok]
            ax_b.axhspan(y0 - GROUP_GAP * 0.45, y1 + GROUP_GAP * 0.45,
                         color=COL_ROW_ODD, zorder=0)

    for nm in hap_names:
        tok = chrom_token(nm)
        bar_specs = [
            (yp_single.get(nm), avg_s,   COL_COV_S_DARK, COL_COV_S),
            (yp_dual.get(nm),   avg_d,   COL_COV_D_DARK, COL_COV_D),
            (yp_ont.get(nm),    avg_ont, COL_COV_O_DARK, COL_COV_O),
        ]
        for yc, cov_avg, col_dark, col_light in bar_specs:
            if yc is None or not cov_avg or tok not in cov_avg:
                continue
            cp, ncp, up = cov_avg[tok]
            ax_b.barh(yc - hb, cp,  height=BAR_H, left=0,        align='edge',
                      color=col_dark,  ec='none', zorder=2)
            ax_b.barh(yc - hb, ncp, height=BAR_H, left=cp,       align='edge',
                      color=col_light, alpha=0.8,  ec='none', zorder=2)
            ax_b.barh(yc - hb, up,  height=BAR_H, left=cp + ncp, align='edge',
                      color=COL_UNCOV, ec='none', zorder=2)
            ax_b.add_patch(mpatches.Rectangle((0, yc - hb), 100, BAR_H,
                           fc='none', ec=COL_BORDER, lw=st['border_lw'], zorder=4))

    ax_b.set_yticks([])
    ax_b.set_xticks([0, 50, 100])
    ax_b.set_xticklabels(['0', '50', '100'])
    if show_xlabel:
        ax_b.set_xlabel('Coverage (%)', fontsize=st['font_xlabel'], labelpad=3)
    ax_b.tick_params(axis='x', labelsize=st['font_tick'], which='major', length=2)
    ax_b.spines['top'].set_visible(False)
    ax_b.spines['right'].set_visible(False)
    ax_b.spines['left'].set_visible(False)
    ax_b.spines['bottom'].set_linewidth(0.4)
    ax_b.tick_params(left=False)
    if rasterized:
        ax_b.set_rasterization_zorder(5)


def _draw_haplotype_panel(ax, fig_w, fig_h,
                          hap_names,
                          yp_single, yp_dual, y_spans, y_label,
                          group_order, total_y, margin_x, margin_y,
                          max_size, size_of,
                          cov_s, nc_cov_s, switch_s,
                          cov_d, nc_cov_d, switch_d,
                          telo_s, telo_nc_s,
                          telo_d, telo_nc_d,
                          BAR_H, GROUP_GAP,
                          show_chrom_labels=True,
                          label_side='left',
                          show_xlabel=False,
                          mirror_x=False,
                          centromeres=None,
                          flip_set=None,
                          style=None,
                          rasterized=False,
                          simplify=False,
                          yp_ont=None,
                          cov_ont=None,
                          nc_cov_ont=None,
                          switch_ont=None,
                          telo_ont=None,
                          telo_nc_ont=None,
                          gaps_dual=None,
                          gaps_single=None,
                          gaps_ont=None,
                          gaps_dual_asm=None,
                          gaps_single_asm=None,
                          gaps_ont_asm=None,
                          yp_sub_single=None,
                          yp_sub_dual=None,
                          yp_sub_ont=None,
                          SUB_H=0.0,
                          insertions_single=None,
                          insertions_dual=None,
                          insertions_ont=None,
                          insertion_track='none'):
    st = style or STYLES['paper']
    h  = BAR_H / 2

    if mirror_x:
        ax.set_xlim(max_size + margin_x, -margin_x)
    else:
        ax.set_xlim(-margin_x, max_size + margin_x)
    ax.set_ylim(total_y + margin_y, -margin_y)
    ax.set_facecolor('white')

    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{v/1e6:.0f}'))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(10e6))
    if show_xlabel:
        ax.set_xlabel('Position (Mb)', fontsize=st['font_xlabel'], labelpad=3)
    ax.tick_params(axis='x', labelsize=st['font_tick'], which='major', length=2.5)
    ax.tick_params(axis='x', which='minor', length=1.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_linewidth(st['border_lw'] * 0.7)
    ax.tick_params(left=False)

    for i, tok in enumerate(group_order):
        if i % 2 == 0 and tok in y_spans:
            y0, y1 = y_spans[tok]
            ax.axhspan(y0 - GROUP_GAP * 0.45, y1 + GROUP_GAP * 0.45,
                       color=COL_ROW_ODD, zorder=0)

    ax_data_trans = blended_transform_factory(ax.transAxes, ax.transData)

    _pos    = ax.get_position()
    _ax_w   = _pos.width  * fig_w
    _ax_h   = _pos.height * fig_h
    _x_span = abs(ax.get_xlim()[1] - ax.get_xlim()[0])
    _y_span = abs(ax.get_ylim()[1] - ax.get_ylim()[0])
    telo_ry = h
    _natural_rx = telo_ry * (_ax_h / _ax_w) * (_x_span / _y_span)
    telo_rx = max(_natural_rx, max_size * 0.006)

    def _telo_semi(x_base, side, y_low, y_high, hollow=False):
        ymid  = (y_low + y_high) / 2
        theta = (np.linspace(np.pi / 2, 3 * np.pi / 2, 60) if side == 'p'
                 else np.linspace(-np.pi / 2, np.pi / 2, 60))
        xs = x_base + telo_rx * np.cos(theta)
        ys = ymid   + telo_ry * np.sin(theta)
        fc = 'none' if hollow else COL_TELO
        lw = st['telo_lw_h'] if hollow else st['telo_lw_f']
        ax.add_patch(mpatches.Polygon(list(zip(xs, ys)), closed=True,
                                      fc=fc, ec=COL_TELO, lw=lw,
                                      zorder=10, clip_on=False))

    def _draw_bar(nm, yc, col_dark, col_light, cov_, nc_cov_, sw_, telo_, telo_nc_,
                  gaps_=None, gaps_asm_=None,
                  y_sub=None, SUB_H=0.0, insertions_=None, insertion_track='none'):
        size = size_of(nm)
        flip = flip_set is not None and nm in flip_set

        def _fi(ivs):
            return [(size - e, size - s) for s, e in ivs] if flip else list(ivs)

        def _fs(ivs):
            raw = [(size - e, size - s, t) for s, e, t in ivs] if flip else list(ivs)
            min_w = max_size * 0.004
            result = []
            for s, e, t in raw:
                if e - s < min_w:
                    mid = (s + e) / 2
                    s = max(0, mid - min_w / 2)
                    e = min(size, mid + min_w / 2)
                result.append((s, e, t))
            return result

        def _fa(arm_set):
            return {('q' if a == 'p' else 'p') for a in arm_set} if flip else arm_set

        cen_cs = cen_ce = None
        if centromeres:
            tok = chrom_token(nm)
            hap = 'pat' if is_pat(nm) else 'mat'
            cen = centromeres.get((tok, hap))
            if cen:
                cs, ce = cen
                if flip:
                    cs, ce = size - ce, size - cs
                cx = (cs + ce) / 2
                cw = max_size * 0.03
                cen_cs = max(0, cx - cw)
                cen_ce = min(size, cx + cw)

        if cen_cs is not None:
            chrom_p = _constriction_path(0, size, yc, BAR_H, cen_cs, cen_ce)
        else:
            chrom_p = _rect_path(0, size, yc, BAR_H)

        ax.add_patch(PathPatch(chrom_p, fc='white', ec=COL_BORDER,
                               lw=st['border_lw'], zorder=2))

        def _band(intervals, color, alpha, z, _cp=chrom_p):
            if not intervals:
                return
            if simplify:
                verts, codes = [], []
                for s, e in intervals:
                    verts.extend([(s, yc-h), (e, yc-h), (e, yc+h), (s, yc+h)])
                    codes.extend([mpath.Path.MOVETO, mpath.Path.LINETO,
                                  mpath.Path.LINETO, mpath.Path.LINETO])
                p = PathPatch(mpath.Path(np.array(verts), np.array(codes), closed=True),
                              fc=color, ec='none', alpha=alpha, zorder=z)
                p.set_clip_path(_cp, transform=ax.transData)
                ax.add_patch(p)
            else:
                for s, e in intervals:
                    r = mpatches.Rectangle((s, yc - h), e - s, BAR_H,
                                           fc=color, ec='none', alpha=alpha, zorder=z)
                    r.set_clip_path(_cp, transform=ax.transData)
                    ax.add_patch(r)

        _band(_fi(nc_cov_.get(nm, [])), col_light, 0.80, 3)
        _band(_fi(cov_.get(nm,    [])), col_dark,  1.00, 4)

        if cov_.get(nm) or nc_cov_.get(nm):
            sw_ivs = _fs(sw_.get(nm, []))
            if simplify and sw_ivs:
                verts, codes = [], []
                for s, e, _ in sw_ivs:
                    verts.extend([(s, yc-h), (e, yc-h), (e, yc+h), (s, yc+h)])
                    codes.extend([mpath.Path.MOVETO, mpath.Path.LINETO,
                                  mpath.Path.LINETO, mpath.Path.LINETO])
                p = PathPatch(mpath.Path(np.array(verts), np.array(codes), closed=True),
                              fc=COL_SW, ec='none', alpha=0.85, zorder=7)
                p.set_clip_path(chrom_p, transform=ax.transData)
                ax.add_patch(p)
            else:
                for s, e, _ in sw_ivs:
                    r = mpatches.Rectangle((s, yc - h), e - s, BAR_H,
                                           fc=COL_SW, ec='none', alpha=0.85, zorder=7)
                    r.set_clip_path(chrom_p, transform=ax.transData)
                    ax.add_patch(r)

        cur_gap_lw = max(0.8, st['border_lw'] * 1.8)
        asm_gap_lw = max(0.5, st['border_lw'] * 1.2)

        if gaps_asm_:
            for gps, gpe in gaps_asm_.get(nm, []):
                if flip:
                    gps, gpe = size - gpe, size - gps
                mid_g = (gps + gpe) / 2
                line, = ax.plot([mid_g, mid_g], [yc - h, yc + h],
                                color=COL_ASM_GAP, lw=asm_gap_lw, zorder=8,
                                solid_capstyle='butt', clip_on=True)
                line.set_clip_path(chrom_p, transform=ax.transData)

        if gaps_:
            for gps, gpe in gaps_.get(nm, []):
                if flip:
                    gps, gpe = size - gpe, size - gps
                mid_g = (gps + gpe) / 2
                line, = ax.plot([mid_g, mid_g], [yc - h, yc + h],
                                color=COL_GAP, lw=cur_gap_lw, zorder=9,
                                solid_capstyle='butt', clip_on=True)
                line.set_clip_path(chrom_p, transform=ax.transData)

        ax.add_patch(PathPatch(chrom_p, fc='none', ec=COL_BORDER,
                               lw=st['border_lw'], zorder=13))
        ax.plot([0, 0],       [yc - h, yc + h], color=COL_BORDER,
                lw=st['border_lw'], solid_capstyle='butt', zorder=14, clip_on=False)
        ax.plot([size, size], [yc - h, yc + h], color=COL_BORDER,
                lw=st['border_lw'], solid_capstyle='butt', zorder=14, clip_on=False)

        if telo_:
            arms = _fa(telo_.get(nm, set()))
            if 'p' in arms: _telo_semi(0,    'p', yc - h, yc + h)
            if 'q' in arms: _telo_semi(size, 'q', yc - h, yc + h)
        if telo_nc_:
            arms = _fa(telo_nc_.get(nm, set()))
            if 'p' in arms: _telo_semi(0,    'p', yc - h, yc + h, hollow=True)
            if 'q' in arms: _telo_semi(size, 'q', yc - h, yc + h, hollow=True)

        # ---- Draw Query Insertion Sub-track (Option A or Option B) ----------
        if y_sub is not None and SUB_H > 0 and insertion_track in ('optionA', 'optionB'):
            h_sub = SUB_H / 2
            ax.add_patch(mpatches.Rectangle((0, y_sub - h_sub), size, SUB_H,
                                           fc='#F8F9FA', ec='#D0D7DE', lw=st['border_lw'] * 0.5, zorder=2))
            ins_list = (insertions_ or {}).get(nm, [])
            if insertion_track == 'optionA':
                for ins in ins_list:
                    t_anchor = ins['t_anchor']
                    if flip:
                        t_anchor = size - t_anchor
                    w = max(max_size * 0.007, 200_000)
                    if ins['dt'] > w:
                        w = ins['dt']
                    x0 = max(0, t_anchor - w / 2)
                    x1 = min(size, t_anchor + w / 2)
                    ins_patch = mpatches.Rectangle((x0, y_sub - h_sub), x1 - x0, SUB_H,
                                                   fc='#E76F51', ec='#2C3E50', lw=st['border_lw'] * 0.8,
                                                   alpha=0.9, zorder=5)
                    ax.add_patch(ins_patch)
                    
                    cgaps = ins.get('cur_gaps', [])
                    if len(cgaps) == 1:
                        mid_x = (x0 + x1) / 2
                        ax.plot([mid_x, mid_x], [y_sub - h_sub, y_sub + h_sub],
                                color='#B22222', lw=max(0.7, st['border_lw'] * 1.6), zorder=8, solid_capstyle='butt')
                    elif len(cgaps) > 1:
                        for rel_pos, gs, ge in cgaps:
                            frac = min(1.0, max(0.0, rel_pos / max(ins['dq'], 1)))
                            tick_x = x0 + (x1 - x0) * (0.12 + 0.76 * frac)
                            ax.plot([tick_x, tick_x], [y_sub - h_sub, y_sub + h_sub],
                                    color='#B22222', lw=max(0.7, st['border_lw'] * 1.6), zorder=8, solid_capstyle='butt')

                    agaps = ins.get('asm_gaps', [])
                    for rel_pos, gs, ge in agaps:
                        frac = min(1.0, max(0.0, rel_pos / max(ins['dq'], 1)))
                        tick_x = x0 + (x1 - x0) * (0.12 + 0.76 * frac)
                        ax.plot([tick_x, tick_x], [y_sub - h_sub, y_sub + h_sub],
                                color=COL_ASM_GAP, lw=max(0.5, st['border_lw'] * 1.2), zorder=7, solid_capstyle='butt')

                    if ins['dq'] >= 500_000 or len(cgaps) >= 3:
                        n_g = len(cgaps)
                        txt = f"+{ins['dq']/1e6:.1f}M ({n_g})" if n_g > 0 else f"+{ins['dq']/1e6:.1f}M"
                        ax.text(t_anchor, y_sub + h_sub + 0.04 * BAR_H, txt,
                                fontsize=max(3.8, st['font_n'] * 0.75), color='#8B2500', fontweight='bold',
                                ha='center', va='top', zorder=12)

            elif insertion_track == 'optionB':
                for ins in ins_list:
                    t_anchor = ins['t_anchor']
                    if flip:
                        t_anchor = size - t_anchor
                    seg_w = max(max_size * 0.015, min(max_size * 0.09, 1_200_000 * np.log10(ins['dq'] / 10_000 + 1)))
                    x_left  = max(0, t_anchor - seg_w / 2)
                    x_right = min(size, t_anchor + seg_w / 2)
                    
                    ax.plot([t_anchor, x_left],  [yc + h, y_sub - h_sub], color='#E76F51', lw=st['border_lw'] * 0.8, alpha=0.6, zorder=3)
                    ax.plot([t_anchor, x_right], [yc + h, y_sub - h_sub], color='#E76F51', lw=st['border_lw'] * 0.8, alpha=0.6, zorder=3)
                    ax.fill([t_anchor, x_right, x_left], [yc + h, y_sub - h_sub, y_sub - h_sub],
                            color='#FCEADE', alpha=0.25, zorder=2)
                    
                    seg_patch = mpatches.Rectangle((x_left, y_sub - h_sub), x_right - x_left, SUB_H,
                                                   fc='#FCEADE', ec='#E76F51', lw=st['border_lw'] * 0.8, zorder=5)
                    ax.add_patch(seg_patch)
                    
                    cgaps = ins.get('cur_gaps', [])
                    for rel_pos, gs, ge in cgaps:
                        frac = min(1.0, max(0.0, rel_pos / max(ins['dq'], 1)))
                        tick_x = x_left + frac * (x_right - x_left)
                        ax.plot([tick_x, tick_x], [y_sub - h_sub, y_sub + h_sub],
                                color='#B22222', lw=max(0.7, st['border_lw'] * 1.6), zorder=8, solid_capstyle='butt')
                                
                    agaps = ins.get('asm_gaps', [])
                    for rel_pos, gs, ge in agaps:
                        frac = min(1.0, max(0.0, rel_pos / max(ins['dq'], 1)))
                        tick_x = x_left + frac * (x_right - x_left)
                        ax.plot([tick_x, tick_x], [y_sub - h_sub, y_sub + h_sub],
                                color=COL_ASM_GAP, lw=max(0.5, st['border_lw'] * 1.2), zorder=7, solid_capstyle='butt')

                    if ins['dq'] >= 500_000 or len(cgaps) >= 3:
                        n_g = len(cgaps)
                        txt = f"+{ins['dq']/1e6:.1f}Mb ({n_g} joins)" if n_g > 0 else f"+{ins['dq']/1e6:.1f}Mb"
                        ax.text((x_left + x_right) / 2, y_sub + h_sub + 0.04 * BAR_H, txt,
                                fontsize=max(3.8, st['font_n'] * 0.75), color='#8B2500', fontweight='bold',
                                ha='center', va='top', zorder=12)

    for nm in hap_names:
        if nm in yp_single:
            _draw_bar(nm, yp_single[nm],
                      COL_COV_S_DARK, COL_COV_S,
                      cov_s, nc_cov_s, switch_s, telo_s, telo_nc_s,
                      gaps_=gaps_single, gaps_asm_=gaps_single_asm,
                      y_sub=yp_sub_single.get(nm) if yp_sub_single else None,
                      SUB_H=SUB_H, insertions_=insertions_single, insertion_track=insertion_track)
        if nm in yp_dual:
            _draw_bar(nm, yp_dual[nm],
                      COL_COV_D_DARK, COL_COV_D,
                      cov_d, nc_cov_d, switch_d, telo_d, telo_nc_d,
                      gaps_=gaps_dual, gaps_asm_=gaps_dual_asm,
                      y_sub=yp_sub_dual.get(nm) if yp_sub_dual else None,
                      SUB_H=SUB_H, insertions_=insertions_dual, insertion_track=insertion_track)
        if yp_ont and nm in yp_ont:
            _draw_bar(nm, yp_ont[nm],
                      COL_COV_O_DARK, COL_COV_O,
                      cov_ont or {}, nc_cov_ont or {}, switch_ont or {},
                      telo_ont, telo_nc_ont,
                      gaps_=gaps_ont, gaps_asm_=gaps_ont_asm,
                      y_sub=yp_sub_ont.get(nm) if yp_sub_ont else None,
                      SUB_H=SUB_H, insertions_=insertions_ont, insertion_track=insertion_track)

    ax.set_yticks([])

    _lbl_remap = {'W': 'W/Z'}
    if show_chrom_labels:
        for tok in group_order:
            if tok not in y_label:
                continue
            disp = _lbl_remap.get(tok, tok)
            if label_side == 'right':
                ax.text(1.045, y_label[tok], disp,
                        transform=ax_data_trans, ha='left', va='center',
                        fontsize=st['font_chrom'], fontweight='bold',
                        clip_on=False, zorder=20)
            else:
                ax.text(-0.045, y_label[tok], disp,
                        transform=ax_data_trans, ha='right', va='center',
                        fontsize=st['font_chrom'], fontweight='bold',
                        clip_on=False, zorder=20)

    if rasterized:
        ax.set_rasterization_zorder(12)


# =============================================================================
# Main figure builder
# =============================================================================

def build_haplotype_figure(df_s, df_d,
                           cov_s, cov_d, nc_cov_s, nc_cov_d,
                           switch_s, switch_d,
                           telo_s, telo_d, telo_nc_s, telo_nc_d,
                           cov_summary_s, cov_summary_d, cov_summary_ont,
                           cov_ont, nc_cov_ont, switch_ont,
                           telo_ont, telo_nc_ont,
                           out_path, dpi, fmt,
                           style_name='paper', layout='butterfly',
                           rasterized=False,
                           simplify=False,
                           centromeres=None,
                           flip_set=None,
                           gaps_dual=None,
                           gaps_single=None,
                           gaps_ont=None,
                            gaps_dual_asm=None,
                            gaps_single_asm=None,
                            gaps_ont_asm=None,
                            insertions_dual=None,
                            insertions_single=None,
                            insertions_ont=None,
                            insertion_track='none'):

    all_names = sorted(
        set(df_s.index) | set(df_d.index),
        key=lambda nm: (chrom_sort_key(chrom_token(nm)), 0 if is_pat(nm) else 1)
    )
    chrom_groups = {}
    chrom_order  = []
    for nm in all_names:
        tok = chrom_token(nm)
        if tok not in chrom_groups:
            chrom_groups[tok] = []
            chrom_order.append(tok)
        chrom_groups[tok].append(nm)

    macro_order = [t for t in chrom_order if classify_chrom_group(t) == 'macro']
    micro_order = [t for t in chrom_order if classify_chrom_group(t) == 'micro']
    nano_order  = [t for t in chrom_order if classify_chrom_group(t) == 'nano']

    def size_of(nm):
        return int(df_s.loc[nm, 'size'] if nm in df_s.index else df_d.loc[nm, 'size'])

    st = STYLES.get(style_name, STYLES['paper'])
    plt.rcParams.update({
        'font.family':       'sans-serif',
        'font.sans-serif':   ['Arial', 'Helvetica Neue', 'Helvetica', 'DejaVu Sans'],
        'font.size':         st['font_base'],
        'axes.linewidth':    st['border_lw'] * 0.7,
        'xtick.major.width': st['border_lw'] * 0.7,
        'xtick.minor.width': st['border_lw'] * 0.5,
        'xtick.major.size':  3.0,
        'xtick.minor.size':  1.8,
    })

    BAR_H     = st['BAR_H']
    SUB_H     = BAR_H * 0.32 if insertion_track != 'none' else 0.0
    SUB_GAP   = BAR_H * 0.12 if insertion_track != 'none' else 0.0
    HAP_GAP   = BAR_H * 0.14
    GROUP_GAP = BAR_H * (0.58 if insertion_track != 'none' else 0.52)
    is_mat    = lambda nm: not is_pat(nm)

    yps_pat_mac, ypd_pat_mac, ypo_pat_mac, y_sub_s_pat_mac, y_sub_d_pat_mac, y_sub_o_pat_mac, yl_pat_mac, ys_pat_mac, _ty_p_mac = \
        make_y_layout_triple(macro_order, chrom_groups, is_pat, BAR_H, HAP_GAP, GROUP_GAP, SUB_H, SUB_GAP)
    yps_mat_mac, ypd_mat_mac, ypo_mat_mac, y_sub_s_mat_mac, y_sub_d_mat_mac, y_sub_o_mat_mac, yl_mat_mac, ys_mat_mac, _ty_m_mac = \
        make_y_layout_triple(macro_order, chrom_groups, is_mat, BAR_H, HAP_GAP, GROUP_GAP, SUB_H, SUB_GAP)
    ty_mac = max(_ty_p_mac, _ty_m_mac)

    yps_pat_mic, ypd_pat_mic, ypo_pat_mic, y_sub_s_pat_mic, y_sub_d_pat_mic, y_sub_o_pat_mic, yl_pat_mic, ys_pat_mic, _ty_p_mic = \
        make_y_layout_triple(micro_order, chrom_groups, is_pat, BAR_H, HAP_GAP, GROUP_GAP, SUB_H, SUB_GAP)
    yps_mat_mic, ypd_mat_mic, ypo_mat_mic, y_sub_s_mat_mic, y_sub_d_mat_mic, y_sub_o_mat_mic, yl_mat_mic, ys_mat_mic, _ty_m_mic = \
        make_y_layout_triple(micro_order, chrom_groups, is_mat, BAR_H, HAP_GAP, GROUP_GAP, SUB_H, SUB_GAP)
    ty_mic = max(_ty_p_mic, _ty_m_mic)

    yps_pat_nan, ypd_pat_nan, ypo_pat_nan, y_sub_s_pat_nan, y_sub_d_pat_nan, y_sub_o_pat_nan, yl_pat_nan, ys_pat_nan, _ty_p_nan = \
        make_y_layout_triple(nano_order,  chrom_groups, is_pat, BAR_H, HAP_GAP, GROUP_GAP, SUB_H, SUB_GAP)
    yps_mat_nan, ypd_mat_nan, ypo_mat_nan, y_sub_s_mat_nan, y_sub_d_mat_nan, y_sub_o_mat_nan, yl_mat_nan, ys_mat_nan, _ty_m_nan = \
        make_y_layout_triple(nano_order,  chrom_groups, is_mat, BAR_H, HAP_GAP, GROUP_GAP, SUB_H, SUB_GAP)
    ty_nan = max(_ty_p_nan, _ty_m_nan)

    def _hap_names(grp_order, filt):
        return [nm for tok in grp_order
                for nm in chrom_groups.get(tok, []) if filt(nm)]

    names_pat_mac = _hap_names(macro_order, is_pat)
    names_pat_mic = _hap_names(micro_order, is_pat)
    names_pat_nan = _hap_names(nano_order,  is_pat)
    names_mat_mac = _hap_names(macro_order, is_mat)
    names_mat_mic = _hap_names(micro_order, is_mat)
    names_mat_nan = _hap_names(nano_order,  is_mat)

    def _max_size(names):
        return max((size_of(nm) for nm in names), default=1)

    max_mac  = _max_size(names_pat_mac + names_mat_mac)
    max_mic  = _max_size(names_pat_mic + names_mat_mic)
    max_nan  = _max_size(names_pat_nan + names_mat_nan)
    margin_y = GROUP_GAP * 0.6

    has_avg = (cov_summary_s is not None) or (cov_summary_d is not None) or \
              (cov_summary_ont is not None)
    avg_s_mac, avg_d_mac, avg_ont_mac = compute_avg_coverage_by_method(
        macro_order, chrom_groups, cov_summary_s, cov_summary_d, cov_summary_ont)
    avg_s_mic, avg_d_mic, avg_ont_mic = compute_avg_coverage_by_method(
        micro_order, chrom_groups, cov_summary_s, cov_summary_d, cov_summary_ont)
    avg_s_nan, avg_d_nan, avg_ont_nan = compute_avg_coverage_by_method(
        nano_order,  chrom_groups, cov_summary_s, cov_summary_d, cov_summary_ont)

    fig_w = st['fig_w']
    fig_h = st['fig_h'] * (1.18 if insertion_track != 'none' else 1.0)
    fig   = plt.figure(figsize=(fig_w, fig_h), facecolor='white')

    hr = [ty_mac, ty_mic, ty_nan]

    if has_avg:
        ncols = 3
        wr    = [3, 3, 1]
        mat_col, pat_col, avg_col = 0, 1, 2
    else:
        ncols = 2
        wr    = [1, 1]
        mat_col, pat_col = 0, 1
        avg_col = None

    gs = GridSpec(3, ncols, figure=fig,
                  width_ratios=wr, height_ratios=hr,
                  hspace=0.14, wspace=0.06,
                  left=0.07, right=0.98, top=0.93, bottom=0.06)

    ax_mm  = fig.add_subplot(gs[0, mat_col])
    ax_mu  = fig.add_subplot(gs[1, mat_col])
    ax_mn  = fig.add_subplot(gs[2, mat_col])
    ax_pm  = fig.add_subplot(gs[0, pat_col])
    ax_pu  = fig.add_subplot(gs[1, pat_col])
    ax_pn  = fig.add_subplot(gs[2, pat_col])
    ax_avg_mac = fig.add_subplot(gs[0, avg_col]) if has_avg else None
    ax_avg_mic = fig.add_subplot(gs[1, avg_col]) if has_avg else None
    ax_avg_nan = fig.add_subplot(gs[2, avg_col]) if has_avg else None

    mirror_mat_bf = (layout == 'butterfly')

    row_specs = [
        (ax_mm, ax_pm, ax_avg_mac,
         names_mat_mac, yps_mat_mac, ypd_mat_mac, ypo_mat_mac, y_sub_s_mat_mac, y_sub_d_mat_mac, y_sub_o_mat_mac, yl_mat_mac, ys_mat_mac,
         names_pat_mac, yps_pat_mac, ypd_pat_mac, ypo_pat_mac, y_sub_s_pat_mac, y_sub_d_pat_mac, y_sub_o_pat_mac, yl_pat_mac, ys_pat_mac,
         ty_mac, max_mac, macro_order, avg_s_mac, avg_d_mac, avg_ont_mac),
        (ax_mu, ax_pu, ax_avg_mic,
         names_mat_mic, yps_mat_mic, ypd_mat_mic, ypo_mat_mic, y_sub_s_mat_mic, y_sub_d_mat_mic, y_sub_o_mat_mic, yl_mat_mic, ys_mat_mic,
         names_pat_mic, yps_pat_mic, ypd_pat_mic, ypo_pat_mic, y_sub_s_pat_mic, y_sub_d_pat_mic, y_sub_o_pat_mic, yl_pat_mic, ys_pat_mic,
         ty_mic, max_mic, micro_order, avg_s_mic, avg_d_mic, avg_ont_mic),
        (ax_mn, ax_pn, ax_avg_nan,
         names_mat_nan, yps_mat_nan, ypd_mat_nan, ypo_mat_nan, y_sub_s_mat_nan, y_sub_d_mat_nan, y_sub_o_mat_nan, yl_mat_nan, ys_mat_nan,
         names_pat_nan, yps_pat_nan, ypd_pat_nan, ypo_pat_nan, y_sub_s_pat_nan, y_sub_d_pat_nan, y_sub_o_pat_nan, yl_pat_nan, ys_pat_nan,
         ty_nan, max_nan, nano_order, avg_s_nan, avg_d_nan, avg_ont_nan),
    ]

    n_rows = len(row_specs)
    for row_i, (ax_m, ax_p, ax_avg,
                m_nms, yps_m, ypd_m, ypo_m, y_sub_s_m, y_sub_d_m, y_sub_o_m, yl_m, ys_m,
                p_nms, yps_p, ypd_p, ypo_p, y_sub_s_p, y_sub_d_p, y_sub_o_p, yl_p, ys_p,
                ty, max_sz, grp_order,
                avg_s_grp, avg_d_grp, avg_ont_grp) in enumerate(row_specs):

        if not p_nms and not m_nms:
            for ax in [ax_m, ax_p]:
                ax.set_visible(False)
            if ax_avg:
                ax_avg.set_visible(False)
            continue

        mx        = max_sz * 0.008
        is_bottom = (row_i == n_rows - 1)

        # Maternal panel
        _draw_haplotype_panel(ax_m, fig_w, fig_h,
                              m_nms, yps_m, ypd_m, ys_m, yl_m,
                              grp_order, ty, mx, margin_y,
                              max_sz, size_of,
                              cov_s, nc_cov_s, switch_s,
                              cov_d, nc_cov_d, switch_d,
                              telo_s, telo_nc_s, telo_d, telo_nc_d,
                              BAR_H, GROUP_GAP,
                              show_chrom_labels=True,
                              label_side='left',
                              show_xlabel=is_bottom,
                              mirror_x=mirror_mat_bf,
                              centromeres=centromeres,
                              flip_set=flip_set,
                              style=st, rasterized=rasterized,
                              simplify=simplify,
                              yp_ont=ypo_m,
                              cov_ont=cov_ont, nc_cov_ont=nc_cov_ont,
                              switch_ont=switch_ont,
                              telo_ont=telo_ont, telo_nc_ont=telo_nc_ont,
                              gaps_dual=gaps_dual,
                              gaps_single=gaps_single,
                              gaps_ont=gaps_ont,
                              gaps_dual_asm=gaps_dual_asm,
                              gaps_single_asm=gaps_single_asm,
                              gaps_ont_asm=gaps_ont_asm,
                              yp_sub_single=y_sub_s_m,
                              yp_sub_dual=y_sub_d_m,
                              yp_sub_ont=y_sub_o_m,
                              SUB_H=SUB_H,
                              insertions_single=insertions_single,
                              insertions_dual=insertions_dual,
                              insertions_ont=insertions_ont,
                              insertion_track=insertion_track)

        # Paternal panel
        _draw_haplotype_panel(ax_p, fig_w, fig_h,
                              p_nms, yps_p, ypd_p, ys_p, yl_p,
                              grp_order, ty, mx, margin_y,
                              max_sz, size_of,
                              cov_s, nc_cov_s, switch_s,
                              cov_d, nc_cov_d, switch_d,
                              telo_s, telo_nc_s, telo_d, telo_nc_d,
                              BAR_H, GROUP_GAP,
                              show_chrom_labels=False,
                              label_side='right',
                              show_xlabel=is_bottom,
                              mirror_x=False,
                              centromeres=centromeres,
                              flip_set=flip_set,
                              style=st, rasterized=rasterized,
                              simplify=simplify,
                              yp_ont=ypo_p,
                              cov_ont=cov_ont, nc_cov_ont=nc_cov_ont,
                              switch_ont=switch_ont,
                              telo_ont=telo_ont, telo_nc_ont=telo_nc_ont,
                              gaps_dual=gaps_dual,
                              gaps_single=gaps_single,
                              gaps_ont=gaps_ont,
                              gaps_dual_asm=gaps_dual_asm,
                              gaps_single_asm=gaps_single_asm,
                              gaps_ont_asm=gaps_ont_asm,
                              yp_sub_single=y_sub_s_p,
                              yp_sub_dual=y_sub_d_p,
                              yp_sub_ont=y_sub_o_p,
                              SUB_H=SUB_H,
                              insertions_single=insertions_single,
                              insertions_dual=insertions_dual,
                              insertions_ont=insertions_ont,
                              insertion_track=insertion_track)

        # Coverage panel
        if ax_avg is not None:
            _draw_avg_cov_panel(ax_avg, m_nms,
                                yps_m, ypd_m, ypo_m, ys_m,
                                grp_order, ty, margin_y, BAR_H, GROUP_GAP,
                                avg_s_grp, avg_d_grp, avg_ont_grp,
                                style=st, show_xlabel=is_bottom,
                                rasterized=rasterized)

    # ---- Legend -------------------------------------------------------------
    _telo_coll = mpatches.Patch(fc=COL_TELO, ec=COL_TELO, label='Telomere — collinear')
    _telo_nc   = mpatches.Patch(fc='none',   ec=COL_TELO, label='Telomere — non-collinear')
    legend_handles = [
        mpatches.Patch(fc=COL_COV_S_DARK, ec='none',             label='Single — collinear'),
        mpatches.Patch(fc=COL_COV_D_DARK, ec='none',             label='Dual — collinear'),
        mpatches.Patch(fc=COL_COV_O_DARK, ec='none',             label='ONT Dual — collinear'),
        mpatches.Patch(fc=COL_COV_S,      ec='none', alpha=0.80, label='Single — non-collinear'),
        mpatches.Patch(fc=COL_COV_D,      ec='none', alpha=0.80, label='Dual — non-collinear'),
        mpatches.Patch(fc=COL_COV_O,      ec='none', alpha=0.80, label='ONT Dual — non-collinear'),
        mpatches.Patch(fc=COL_SW,         ec='none',             label='Hap. switch'),
        Line2D([0], [0], color=COL_GAP, lw=1.5,                  label='Curation gap'),
        Line2D([0], [0], color=COL_ASM_GAP, lw=1.2,              label='Assembly gap'),
        mpatches.Patch(fc=COL_BORDER, ec=COL_BORDER, lw=0.5,    label='Centromere'),
        _telo_coll, _telo_nc,
    ]
    if insertion_track in ('optionA', 'optionB'):
        legend_handles.append(mpatches.Patch(fc='#E76F51', ec='#2C3E50', lw=0.5, label='Assembly insertion (absent in T2T)'))
    _cen_handle = legend_handles[9]
    fig.legend(handles=legend_handles, loc='lower center',
               fontsize=st['font_legend'], ncol=st['ncol_legend'],
               framealpha=0.9, bbox_to_anchor=(0.50, 0.002),
               handlelength=1.2, handleheight=1.0,
               columnspacing=0.8, borderpad=0.6,
               frameon=True, edgecolor='#cccccc',
               handler_map={_telo_coll:  _SemiHandler(hollow=False),
                            _telo_nc:    _SemiHandler(hollow=True),
                            _cen_handle: _ConstrictionHandler()})

    # ---- Column headers -----------------------------------------------------
    lw_under = st['border_lw'] * 1.8
    bb_mat   = ax_mm.get_position()
    bb_pat   = ax_pm.get_position()
    header_y = bb_mat.y1 + 0.018

    mat_cx = (bb_mat.x0 + bb_mat.x1) / 2
    pat_cx = (bb_pat.x0 + bb_pat.x1) / 2
    fig.text(mat_cx, header_y, 'Maternal', ha='center', va='bottom',
             fontsize=st['font_title'], fontweight='bold', color='black',
             transform=fig.transFigure)
    fig.text(pat_cx, header_y, 'Paternal', ha='center', va='bottom',
             fontsize=st['font_title'], fontweight='bold', color='black',
             transform=fig.transFigure)
    fig.add_artist(Line2D([bb_mat.x0, bb_mat.x1],
                           [header_y - 0.004, header_y - 0.004],
                           color='black', lw=lw_under,
                           transform=fig.transFigure, solid_capstyle='butt'))
    fig.add_artist(Line2D([bb_pat.x0, bb_pat.x1],
                           [header_y - 0.004, header_y - 0.004],
                           color='black', lw=lw_under,
                           transform=fig.transFigure, solid_capstyle='butt'))

    if has_avg and ax_avg_mac is not None:
        bb_avg = ax_avg_mac.get_position()
        avg_cx = (bb_avg.x0 + bb_avg.x1) / 2
        fig.text(avg_cx, header_y, 'Coverage', ha='center', va='bottom',
                 fontsize=st['font_title'] * 0.75, fontweight='bold', color='black',
                 transform=fig.transFigure)
        fig.add_artist(Line2D([bb_avg.x0, bb_avg.x1],
                               [header_y - 0.004, header_y - 0.004],
                               color='black', lw=lw_under,
                               transform=fig.transFigure, solid_capstyle='butt'))

    out_dir = os.path.dirname(os.path.abspath(out_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches='tight', pad_inches=0.05)
    print(f'Saved: {out_path}')
    plt.close(fig)


# =============================================================================
# Entry point
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Butterfly haplotype ideogram (HiFi Single + Dual and ONT Dual) with coverage panel.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # HiFi inputs
    parser.add_argument('--single-tsv',           required=True,  metavar='FILE')
    parser.add_argument('--dual-tsv',             required=True,  metavar='FILE')
    parser.add_argument('--single-chain',         required=True,  metavar='FILE')
    parser.add_argument('--dual-chain',           required=True,  metavar='FILE')
    parser.add_argument('--single-nc-chain',      required=True,  metavar='FILE')
    parser.add_argument('--dual-nc-chain',        required=True,  metavar='FILE')
    parser.add_argument('--single-unlocs-chain',  required=False, metavar='FILE', default=None)
    parser.add_argument('--dual-unlocs-chain',    required=False, metavar='FILE', default=None)
    parser.add_argument('--single-pairs',         required=True,  metavar='FILE')
    parser.add_argument('--dual-pairs',           required=True,  metavar='FILE')
    parser.add_argument('--single-bed',           required=True,  metavar='FILE')
    parser.add_argument('--dual-bed',             required=True,  metavar='FILE')
    # ONT Dual inputs
    parser.add_argument('--ont-dual-chain',       required=True,  metavar='FILE')
    parser.add_argument('--ont-dual-nc-chain',    required=True,  metavar='FILE')
    parser.add_argument('--ont-dual-unlocs-chain', required=False, metavar='FILE', default=None)
    parser.add_argument('--ont-dual-pairs',       required=False, metavar='FILE', default=None)
    parser.add_argument('--ont-dual-bed',         required=False, metavar='FILE', default=None)
    # Telomeres (clean *_telomere_presence.tsv)
    parser.add_argument('--single-telomere-presence', required=False, metavar='FILE', default=None)
    parser.add_argument('--dual-telomere-presence',   required=False, metavar='FILE', default=None)
    parser.add_argument('--ont-telomere-presence',    required=False, metavar='FILE', default=None)
    # Coverage summaries
    parser.add_argument('--single-coverage-summary', required=False, metavar='FILE', default=None)
    parser.add_argument('--dual-coverage-summary',   required=False, metavar='FILE', default=None)
    parser.add_argument('--ont-coverage-summary',    required=False, metavar='FILE', default=None)
    parser.add_argument('--coverage-summary',        required=False, metavar='FILE', default=None)
    # Annotated Gaps BEDs
    parser.add_argument('--single-annotated-gaps', required=False, metavar='FILE', default=None)
    parser.add_argument('--dual-annotated-gaps',   required=False, metavar='FILE', default=None)
    parser.add_argument('--ont-annotated-gaps',    required=False, metavar='FILE', default=None)
    # Centromeres & Telomeres P-arm
    parser.add_argument('--centromeres',          required=False, metavar='FILE', default=None)
    parser.add_argument('--telo-p-bed',           required=False, metavar='FILE', default=None)
    # Output & Styling
    parser.add_argument('--insertion-track',      default='none', choices=['none', 'optionA', 'optionB'],
                        help='Display unaligned assembly sequence insertions with curation gaps as sub-tracks: '
                             'none (default), optionA (compact block footprint), optionB (expanded segment loops)')
    parser.add_argument('--output',               required=True,  metavar='FILE')
    parser.add_argument('--layout',    default='butterfly', choices=['butterfly', 'columns'])
    parser.add_argument('--style',     default='paper',     choices=['paper', 'poster'])
    parser.add_argument('--format',    default='png',       choices=['png', 'pdf', 'svg'])
    parser.add_argument('--dpi',       type=int, default=None)
    parser.add_argument('--rasterize', action='store_true')
    parser.add_argument('--simplify',  action='store_true')
    parser.add_argument('--datestamp', action='store_true',
                        help='Append YYYYMMDD datestamp to output filename')

    args = parser.parse_args()

    out = args.output
    if not out.endswith(f'.{args.format}'):
        out = f'{out}.{args.format}'
    if args.datestamp:
        base, ext = os.path.splitext(out)
        today = datetime.now().strftime("%Y%m%d")
        if not base.endswith(today):
            out = f'{base}_{today}{ext}'

    print('Loading sequence sizes from TSVs...')
    df_s = load_tsv(args.single_tsv)
    df_d = load_tsv(args.dual_tsv)
    seq_sizes = {nm: int(df_s.loc[nm, 'size'] if nm in df_s.index else df_d.loc[nm, 'size'])
                 for nm in set(df_s.index) | set(df_d.index)}

    print('Computing HiFi collinear coverage...')
    cov_s = compute_covered(args.single_chain, seq_sizes)
    cov_d = compute_covered(args.dual_chain,   seq_sizes)

    print('Computing HiFi non-collinear coverage...')
    nc_cov_s = compute_covered(args.single_nc_chain, seq_sizes)
    nc_cov_d = compute_covered(args.dual_nc_chain,   seq_sizes)

    if args.single_unlocs_chain:
        nc_cov_s = merge_coverage(nc_cov_s, compute_covered(args.single_unlocs_chain, seq_sizes))
    if args.dual_unlocs_chain:
        nc_cov_d = merge_coverage(nc_cov_d, compute_covered(args.dual_unlocs_chain,   seq_sizes))

    print('Computing ONT Dual collinear & non-collinear coverage...')
    cov_ont = compute_covered(args.ont_dual_chain, seq_sizes)
    nc_cov_ont = compute_covered(args.ont_dual_nc_chain, seq_sizes)
    if args.ont_dual_unlocs_chain:
        nc_cov_ont = merge_coverage(nc_cov_ont, compute_covered(args.ont_dual_unlocs_chain, seq_sizes))

    print('Building switch-error lookups...')
    switch_s = build_switch_lookup(args.single_pairs, args.single_bed, args.single_chain)
    switch_d = build_switch_lookup(args.dual_pairs,   args.dual_bed,   args.dual_chain)
    switch_ont = {}
    if args.ont_dual_pairs and args.ont_dual_bed:
        switch_ont = build_switch_lookup(args.ont_dual_pairs, args.ont_dual_bed, args.ont_dual_chain)

    print('Loading telomere presence records...')
    telo_s, telo_nc_s = load_telomere_presence_tsv(args.single_telomere_presence)
    telo_d, telo_nc_d = load_telomere_presence_tsv(args.dual_telomere_presence)
    telo_ont, telo_nc_ont = load_telomere_presence_tsv(args.ont_telomere_presence)

    print('Loading coverage summaries...')
    if args.coverage_summary:
        cov_summary_s, cov_summary_d, cov_summary_ont = load_combined_coverage_summary(args.coverage_summary)
    else:
        cov_summary_s = load_single_coverage_summary(args.single_coverage_summary)
        cov_summary_d = load_single_coverage_summary(args.dual_coverage_summary)
        cov_summary_ont = load_single_coverage_summary(args.ont_coverage_summary)

    centromeres = None
    if args.centromeres:
        print('Loading centromere annotations...')
        centromeres = load_centromeres(args.centromeres)

    flip_set = None
    if args.telo_p_bed:
        print('Building orientation flip set from p-arm telomere positions...')
        p_arm_pos = load_p_arm_bed(args.telo_p_bed)
        flip_set  = build_flip_set(p_arm_pos, seq_sizes)
        if flip_set:
            print(f'  {len(flip_set)} chromosome(s) flipped p-arm-first: {sorted(flip_set)}')
        else:
            print('  All chromosomes already p-arm-first.')

    if centromeres:
        cen_flip = set()
        for nm, size in seq_sizes.items():
            tok = chrom_token(nm)
            hap = 'pat' if is_pat(nm) else 'mat'
            cen = centromeres.get((tok, hap))
            if cen:
                midpoint = (cen[0] + cen[1]) / 2
                if midpoint > size / 2:
                    cen_flip.add(nm)
        if cen_flip:
            print(f'  Centromere flip: {len(cen_flip)} chromosome(s) flipped toward center: {sorted(cen_flip)}')
            flip_set = (flip_set or set()) | cen_flip

    gaps_single_asm = gaps_single = None
    gaps_dual_asm   = gaps_dual   = None
    gaps_ont_asm    = gaps_ont    = None
    ins_single = ins_dual = ins_ont = None

    if args.insertion_track != 'none':
        print(f'Extracting query insertions and classifying gaps (track mode: {args.insertion_track})...')
        if args.single_annotated_gaps:
            asm_s, cur_s = load_annotated_gaps_bed(args.single_annotated_gaps)
            gaps_single_asm, gaps_single, ins_single = extract_insertions_and_gaps(
                asm_s, cur_s, args.single_pairs, args.single_chain, args.single_nc_chain, seq_sizes=seq_sizes)
            print(f'  Single: {sum(len(v) for v in gaps_single.values())} aligned curation gaps, '
                  f'{sum(len(ins["cur_gaps"]) for v in ins_single.values() for ins in v)} inside insertions.')
        if args.dual_annotated_gaps:
            asm_d, cur_d = load_annotated_gaps_bed(args.dual_annotated_gaps)
            gaps_dual_asm, gaps_dual, ins_dual = extract_insertions_and_gaps(
                asm_d, cur_d, args.dual_pairs, args.dual_chain, args.dual_nc_chain, seq_sizes=seq_sizes)
            print(f'  Dual: {sum(len(v) for v in gaps_dual.values())} aligned curation gaps, '
                  f'{sum(len(ins["cur_gaps"]) for v in ins_dual.values() for ins in v)} inside insertions.')
        if args.ont_annotated_gaps:
            asm_o, cur_o = load_annotated_gaps_bed(args.ont_annotated_gaps)
            gaps_ont_asm, gaps_ont, ins_ont = extract_insertions_and_gaps(
                asm_o, cur_o, args.ont_dual_pairs, args.ont_dual_chain, args.ont_dual_nc_chain, seq_sizes=seq_sizes)
            print(f'  ONT: {sum(len(v) for v in gaps_ont.values())} aligned curation gaps, '
                  f'{sum(len(ins["cur_gaps"]) for v in ins_ont.values() for ins in v)} inside insertions.')
    else:
        print('Loading and lifting annotated gaps...')
        if args.single_annotated_gaps:
            asm_s, cur_s = load_annotated_gaps_bed(args.single_annotated_gaps)
            gaps_single_asm = lift_gaps_to_t2t(asm_s, args.single_pairs, args.single_chain, args.single_nc_chain, seq_sizes=seq_sizes, is_curation=False)
            gaps_single     = lift_gaps_to_t2t(cur_s, args.single_pairs, args.single_chain, args.single_nc_chain, seq_sizes=seq_sizes, is_curation=True)
            print(f'  Single: {sum(len(v) for v in gaps_single.values())} curation, {sum(len(v) for v in gaps_single_asm.values())} assembly gaps lifted.')

        if args.dual_annotated_gaps:
            asm_d, cur_d = load_annotated_gaps_bed(args.dual_annotated_gaps)
            gaps_dual_asm = lift_gaps_to_t2t(asm_d, args.dual_pairs, args.dual_chain, args.dual_nc_chain, seq_sizes=seq_sizes, is_curation=False)
            gaps_dual     = lift_gaps_to_t2t(cur_d, args.dual_pairs, args.dual_chain, args.dual_nc_chain, seq_sizes=seq_sizes, is_curation=True)
            print(f'  Dual: {sum(len(v) for v in gaps_dual.values())} curation, {sum(len(v) for v in gaps_dual_asm.values())} assembly gaps lifted.')

        if args.ont_annotated_gaps:
            asm_o, cur_o = load_annotated_gaps_bed(args.ont_annotated_gaps)
            gaps_ont_asm = lift_gaps_to_t2t(asm_o, args.ont_dual_pairs, args.ont_dual_chain, args.ont_dual_nc_chain, seq_sizes=seq_sizes, is_curation=False)
            gaps_ont     = lift_gaps_to_t2t(cur_o, args.ont_dual_pairs, args.ont_dual_chain, args.ont_dual_nc_chain, seq_sizes=seq_sizes, is_curation=True)
            print(f'  ONT: {sum(len(v) for v in gaps_ont.values())} curation, {sum(len(v) for v in gaps_ont_asm.values())} assembly gaps lifted.')

    dpi = args.dpi or STYLES.get(args.style, STYLES['paper'])['dpi_default']
    print(f'Building figure (layout={args.layout}, style={args.style}, dpi={dpi})...')
    build_haplotype_figure(
        df_s, df_d,
        cov_s, cov_d, nc_cov_s, nc_cov_d,
        switch_s, switch_d,
        telo_s, telo_d, telo_nc_s, telo_nc_d,
        cov_summary_s, cov_summary_d, cov_summary_ont,
        cov_ont, nc_cov_ont, switch_ont,
        telo_ont, telo_nc_ont,
        out_path=out, dpi=dpi, fmt=args.format,
        style_name=args.style, layout=args.layout,
        rasterized=args.rasterize,
        simplify=args.simplify,
        centromeres=centromeres,
        flip_set=flip_set,
        gaps_dual=gaps_dual,
        gaps_single=gaps_single,
        gaps_ont=gaps_ont,
        gaps_dual_asm=gaps_dual_asm,
        gaps_single_asm=gaps_single_asm,
        gaps_ont_asm=gaps_ont_asm,
        insertions_dual=ins_dual,
        insertions_single=ins_single,
        insertions_ont=ins_ont,
        insertion_track=args.insertion_track
    )
    print(f'SUCCESS! Figure saved to {out}')


if __name__ == '__main__':
    main()
