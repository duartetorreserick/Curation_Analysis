#!/usr/bin/env python3
"""
plot_combined_figure.py

Combined 6-panel figure:

  Panel A — Macrochromosomes Z and W, maternal only
  Panel B — Microchromosomes 7 and 10, butterfly (mat | pat)
  Panel C — Dot chromosomes 16 and 32, butterfly (mat | pat)
  Panel D — Coverage average per chromosome (horizontal stacked bars)
             occupies 70 % of row D width
  Panel e — Grouped stacked bar: collinear / non-collinear / switch-error
             fractions per category (Macro / Micro / Dot)  [right 30 %]
  Panel f — Lollipop: telomere completeness per category   [right 30 %]

Usage:
  python3 plot_combined_figure.py \\
      --single-tsv FILE --dual-tsv FILE \\
      --single-chain FILE --dual-chain FILE \\
      --single-nc-chain FILE --dual-nc-chain FILE \\
      --ont-dual-chain FILE --ont-dual-nc-chain FILE \\
      [--single-unlocs-chain FILE] [--dual-unlocs-chain FILE] \\
      [--telo-combined-report FILE] [--ont-telo-combined-report FILE] \\
      --dual-pairs FILE --single-pairs FILE \\
      --single-bed FILE --dual-bed FILE \\
      --dual-fai FILE --single-fai FILE \\
      --stats-tsv FILE \\
      [--coverage-summary FILE] [--centromeres FILE] [--telo-p-bed FILE] \\
      --output FILE \\
      [--style paper|poster] [--format png|pdf|svg] [--dpi INT]
      [--rasterize] [--simplify]
"""

import argparse
import re
import os
from collections import defaultdict
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.path as mpath
from matplotlib.patches import PathPatch
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.transforms import blended_transform_factory
from matplotlib.legend_handler import HandlerPatch
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde
import matplotlib.ticker as mticker

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
        fig_w=3600/300, fig_h=2250/300, BAR_H=0.22,
        border_lw=0.25, telo_lw_h=0.4,  telo_lw_f=0.25,
        font_base=9,   font_chrom=11, font_pm=7,
        font_tick=8,   font_xlabel=9,  font_title=11,
        font_n=8,      font_legend=9,  ncol_legend=5,
        dpi_default=300,
        ef_bar_w=0.15, ef_group_gap=0.80, ef_lollipop_offset=0.18,
        ef_stem_lw=0.9, ef_dot_ms=4,
    ),
    'poster': dict(
        fig_w=20.0, fig_h=36.0, BAR_H=1.1,
        border_lw=0.9,  telo_lw_h=1.4,  telo_lw_f=0.8,
        font_base=16,  font_chrom=22, font_pm=13,
        font_tick=16,  font_xlabel=18, font_title=22,
        font_n=14,     font_legend=17, ncol_legend=5,
        dpi_default=150,
        ef_bar_w=0.30, ef_group_gap=1.50, ef_lollipop_offset=0.35,
        ef_stem_lw=2.0, ef_dot_ms=10,
    ),
}

# ---- EF panel constants -----------------------------------------------------
_EF_ASM_ORDER  = ['Single', 'Dual', 'ONT_Dual']
_EF_CAT_ORDER  = ['macro', 'micro', 'dot']
_EF_CAT_LABELS = {'macro': 'Macro', 'micro': 'Micro', 'dot': 'Dot'}
_EF_ASM_COLORS = {
    'Single':   COL_COV_S_DARK,
    'Dual':     COL_COV_D_DARK,
    'ONT_Dual': COL_COV_O_DARK,
}
_EF_ASM_NC_COLORS = {
    'Single':   COL_COV_S,
    'Dual':     COL_COV_D,
    'ONT_Dual': COL_COV_O,
}

# ---- Selected chromosomes per panel -----------------------------------------
MACRO_SELECTED = {'Z', 'W', '4'}
MICRO_SELECTED = {'10', '15'}
NANO_SELECTED  = {'31', '32'}


# =============================================================================
# Legend handlers
# =============================================================================

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
    with open(path) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= 2:
                sizes[p[0]] = int(p[1])
    return sizes


def load_chrom_pairs(path):
    mapping = {}
    with open(path) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= 2:
                mapping[p[0]] = p[1]
    return mapping


def load_switch_blocks(bed_path):
    blocks = defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= 4:
                blocks[p[0]].append((int(p[1]), int(p[2]), p[3]))
    return dict(blocks)


def load_telomere_ends_from_combined_report(path):
    # Format: type(collinear/non-collinear)  chrom  arm  SINGLE(ok/MISSING)  DUAL(ok/MISSING)
    telo_s = defaultdict(set); telo_d = defaultdict(set)
    telo_nc_s = defaultdict(set); telo_nc_d = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 5:
                continue
            kind, chrom, arm, single_val, dual_val = cols[0], cols[1], cols[2], cols[3], cols[4]
            is_coll = kind.strip().lower() == 'collinear'
            if is_coll:
                if single_val.strip().lower() == 'ok':
                    telo_s[chrom].add(arm)
                if dual_val.strip().lower() == 'ok':
                    telo_d[chrom].add(arm)
            else:
                if single_val.strip().lower() == 'ok' and arm not in telo_s.get(chrom, set()):
                    telo_nc_s[chrom].add(arm)
                if dual_val.strip().lower() == 'ok' and arm not in telo_d.get(chrom, set()):
                    telo_nc_d[chrom].add(arm)
    return dict(telo_s), dict(telo_d), dict(telo_nc_s), dict(telo_nc_d)


def load_ont_telomere_report(path):
    telo_ont    = defaultdict(set)
    telo_nc_ont = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 4:
                continue
            chrom, arm, nc_status, coll_status = cols[:4]
            if coll_status.strip().lower() == 'ok':
                telo_ont[chrom].add(arm)
            if nc_status.strip().lower() == 'ok':
                telo_nc_ont[chrom].add(arm)
    return dict(telo_ont), dict(telo_nc_ont)


def load_telomere_presence_tsv(path):
    """
    Parses modern <asm>_telomere_presence.tsv:
    Columns: chromosome, scaffold, collinear, non-collinear, teloscope_rescued, type_of_tele_rescued
    Returns (telo_coll, telo_nc): {chrom: set(arm)}.
    """
    telo_coll = defaultdict(set)
    telo_nc   = defaultdict(set)
    if not path or not os.path.exists(path):
        return dict(telo_coll), dict(telo_nc)
    df = pd.read_csv(path, sep='\t')
    for _, row in df.iterrows():
        chrom = str(row['chromosome']).strip()
        coll = str(row.get('collinear', '')).strip().lower()
        nc = str(row.get('non-collinear', '')).strip().lower()
        rescued = str(row.get('teloscope_rescued', '')).strip().lower() == 'yes'
        rescued_type = str(row.get('type_of_tele_rescued', '')).strip().lower()

        for arm in ['p', 'q']:
            if arm in coll:
                telo_coll[chrom].add(arm)
            elif arm in nc:
                telo_nc[chrom].add(arm)
            elif rescued and arm in rescued_type and arm not in telo_coll[chrom]:
                telo_nc[chrom].add(arm)
    return dict(telo_coll), dict(telo_nc)


def load_p_arm_bed(bed_path):
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


def load_coverage_summary(path_s=None, path_d=None, path_ont=None, combined=None):
    cov_s, cov_d, cov_ont = {}, {}, {}
    if combined and os.path.exists(combined):
        df = pd.read_csv(combined, sep='\t', comment='#')
        for _, row in df.iterrows():
            t = (float(row['collinear_pct']),
                 float(row['noncollinear_pct']),
                 float(row['uncovered_pct']))
            asm = str(row['assembly']).strip()
            if asm == 'Single':     cov_s[row['name']] = t
            elif asm == 'Dual':     cov_d[row['name']] = t
            elif asm in ('ONT_Dual', 'ONT'): cov_ont[row['name']] = t
        return cov_s, cov_d, cov_ont

    def _load_single(p):
        res = {}
        if p and os.path.exists(p):
            df = pd.read_csv(p, sep='\t', comment='#')
            for _, row in df.iterrows():
                res[row['name']] = (float(row['collinear_pct']),
                                    float(row['noncollinear_pct']),
                                    float(row['uncovered_pct']))
        return res

    return _load_single(path_s), _load_single(path_d), _load_single(path_ont)


# =============================================================================
# Switch-error liftover
# =============================================================================

def build_chain_query_liftover(chain_path):
    liftover = defaultdict(lambda: defaultdict(list))

    def _process(hdr, blks):
        t_pos = hdr['tStart']; q_pos = hdr['qStart']
        q_strand = hdr['qStrand']; q_size = hdr['qSize']
        for blk in blks:
            size = blk[0]; dt = blk[1] if len(blk) > 1 else 0; dq = blk[2] if len(blk) > 2 else 0
            qs = (q_size - (q_pos + size)) if q_strand == '-' else q_pos
            qe = (q_size - q_pos)          if q_strand == '-' else q_pos + size
            liftover[hdr['qName']][hdr['tName']].append((qs, qe, t_pos))
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
                if header is not None: _process(header, blocks)
                p = line.split()
                header = {'tName': p[2], 'tSize': int(p[3]), 'tStrand': p[4],
                          'tStart': int(p[5]), 'qName': p[7], 'qSize': int(p[8]),
                          'qStrand': p[9], 'qStart': int(p[10])}
                blocks = []
            else:
                parts = line.split()
                if parts: blocks.append(tuple(int(x) for x in parts))
        if header is not None: _process(header, blocks)

    return {qn: {tn: sorted(blks, key=lambda x: x[0]) for tn, blks in td.items()}
            for qn, td in liftover.items()}


_LIFTOVER_GAP         = 50_000  # split switch block display at T2T gaps larger than this
_CHAIN_EDGE_TOLERANCE = 10      # allow blocks that start/end this many bp outside the chain boundary

def liftover_blocks(raw_blocks, chain_liftover, scaffold_name, t2t_name):
    chain_blocks = chain_liftover.get(scaffold_name, {}).get(t2t_name, [])
    if not chain_blocks:
        return []
    result = []
    for h_start, h_end, _ in raw_blocks:
        overlapping = [(qs, qe, ts) for qs, qe, ts in chain_blocks
                       if max(h_start, qs) < min(h_end, qe)]
        if not overlapping:
            continue
        first_qs = overlapping[0][0];  last_qe = overlapping[-1][1]
        # Drop blocks that extend outside the chain's covered span; allow a
        # small tolerance for chain-boundary edge artifacts (e.g. off-by-one).
        if h_start < first_qs - _CHAIN_EDGE_TOLERANCE or h_end > last_qe + _CHAIN_EDGE_TOLERANCE:
            continue
        # Project each overlapping chain block; merge runs within 50 kb,
        # split at large T2T gaps (e.g. T2T insertions relative to scaffold)
        ivs = []
        for qs, qe, ts in overlapping:
            ov_s = max(h_start, qs);  ov_e = min(h_end, qe)
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


def build_switch_lookup(pairs_path, bed_path, fai_path, chain_path):
    pairs    = load_chrom_pairs(pairs_path)
    switches = load_switch_blocks(bed_path)
    liftover = build_chain_query_liftover(chain_path)
    result = {}
    for t2t_name, super_full in pairs.items():
        raw = switches.get(super_full, [])
        result[t2t_name] = liftover_blocks(raw, liftover, super_full, t2t_name)
    return result


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


def _lift_gap_start(gs, ge, chain_blocks):
    """Lifts a gap position to T2T coordinates via chain blocks."""
    for qs, qe, ts in chain_blocks:
        if qs <= gs <= qe:
            return ts + (gs - qs)
    preceding = [(qs, qe, ts) for qs, qe, ts in chain_blocks if qe < gs]
    if preceding:
        qs, qe, ts = max(preceding, key=lambda x: x[1])
        return ts + (qe - qs)
    following = sorted([(qs, qe, ts) for qs, qe, ts in chain_blocks if qs > gs],
                       key=lambda x: x[0])
    if following:
        return following[0][2]
    return None


def lift_gaps_to_t2t(gaps_raw, pairs_path, chain_path):
    """
    Unified liftover for gaps from assembly coordinates to T2T coordinates.
    gaps_raw keys are assembly sequence names (e.g. Pat_SUPER_2.H2) matching best_chrom_pairs.
    Returns:
      {t2t_name: [(t2t_start, t2t_end), ...]}
    """
    if not gaps_raw or not pairs_path or not chain_path:
        return {}
    pairs = load_chrom_pairs(pairs_path)
    super_to_t2t = {v: k for k, v in pairs.items()}
    liftover = build_chain_query_liftover(chain_path)

    result = defaultdict(list)
    for super_name, intervals in gaps_raw.items():
        t2t_name = super_to_t2t.get(super_name)
        if not t2t_name:
            continue
        chain_blocks = liftover.get(super_name, {}).get(t2t_name, [])
        if not chain_blocks:
            continue
        for gs, ge in intervals:
            pos = _lift_gap_start(gs, ge, chain_blocks)
            if pos is not None:
                result[t2t_name].append((pos - 50, pos + 50))
    return dict(result)


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
    if tok.upper() == 'Z': return True   # Z always on paternal (male) side
    if tok.upper() == 'W': return False  # W always on maternal (female) side
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
# Geometry helpers
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
        [x0,  yc - h], [cs,  yc - h], [cx,  yc - wh], [ce,  yc - h], [x1,  yc - h],
        [x1,  yc + h], [ce,  yc + h], [cx,  yc + wh], [cs,  yc + h], [x0,  yc + h],
    ])
    codes = np.array([mpath.Path.MOVETO] + [mpath.Path.LINETO] * 9)
    return mpath.Path(verts, codes, closed=True)


# =============================================================================
# Y-layout — triple bars per chromosome (Single / Dual / ONT)
# =============================================================================

def make_y_layout_triple(group_order, chrom_groups, hap_filter, BAR_H, HAP_GAP, GROUP_GAP):
    h = BAR_H / 2
    yp_s  = {}
    yp_d  = {}
    yp_o  = {}
    y_lbl = {}
    y_spn = {}
    y = 0.0
    for tok in group_order:
        nm = next((n for n in chrom_groups.get(tok, []) if hap_filter(n)), None)
        if nm is None:
            continue
        span_top   = y
        yp_s[nm]   = y + h;  y += BAR_H
        y          += HAP_GAP
        yp_d[nm]   = y + h;  y += BAR_H
        y          += HAP_GAP
        yp_o[nm]   = y + h;  y += BAR_H
        span_bot   = y
        y_lbl[tok] = (span_top + span_bot) / 2
        y_spn[tok] = (span_top, span_bot)
        y          += GROUP_GAP
    total_y = max(y - GROUP_GAP, BAR_H)
    return yp_s, yp_d, yp_o, y_lbl, y_spn, total_y


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
            vals = [cov_in[nm] for nm in chrom_groups.get(tok, []) if nm in cov_in]
            if vals:
                c  = sum(v[0] for v in vals) / len(vals)
                nc = sum(v[1] for v in vals) / len(vals)
                cov_out[tok] = (c, nc, max(0.0, 100.0 - c - nc))
    return avg_s, avg_d, avg_ont


# =============================================================================
# Draw one haplotype ideogram panel
# =============================================================================

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
                           label_pad=2,
                           show_xlabel=False,
                           mirror_x=False,
                           show_sd=True,
                           centromeres=None,
                           flip_set=None,
                           style=None,
                           rasterized=False,
                           simplify=False,
                           yp_ont=None,
                           cov_ont=None,
                           nc_cov_ont=None,
                            telo_ont=None,
                            telo_nc_ont=None,
                            gaps_s=None,
                            gaps_s_asm=None,
                            gaps_d=None,
                            gaps_d_asm=None,
                            gaps_ont=None,
                            gaps_ont_asm=None):
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
                  gaps_=None, gaps_asm_=None):
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
            cen = centromeres.get((tok, hap)) or centromeres.get((tok, 'pat' if hap == 'mat' else 'mat'))
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

        if gaps_asm_:
            for gps, gpe in gaps_asm_.get(nm, []):
                if flip:
                    gps, gpe = size - gpe, size - gps
                ax.plot([gps, gps], [yc - h, yc + h],
                        color=COL_ASM_GAP, lw=0.6, zorder=8,
                        solid_capstyle='butt')

        if gaps_:
            gap_ivs = gaps_.get(nm, [])
            min_gap_w = max_size * 0.004
            for gps, gpe in gap_ivs:
                if flip:
                    gps, gpe = size - gpe, size - gps
                gw = max(gpe - gps, min_gap_w)
                rect = mpatches.Rectangle((gps, yc - h), gw, BAR_H,
                                          fc=COL_GAP, ec='none', alpha=0.9, zorder=9)
                rect.set_clip_path(chrom_p, transform=ax.transData)
                ax.add_patch(rect)

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

    for nm in hap_names:
        if nm in yp_single:
            _draw_bar(nm, yp_single[nm],
                      COL_COV_S_DARK, COL_COV_S,
                      cov_s, nc_cov_s, switch_s, telo_s, telo_nc_s,
                      gaps_=gaps_s, gaps_asm_=gaps_s_asm)
        if nm in yp_dual:
            _draw_bar(nm, yp_dual[nm],
                      COL_COV_D_DARK, COL_COV_D,
                      cov_d, nc_cov_d, switch_d, telo_d, telo_nc_d,
                      gaps_=gaps_d, gaps_asm_=gaps_d_asm)
        if yp_ont and nm in yp_ont:
            _draw_bar(nm, yp_ont[nm],
                      COL_COV_O_DARK, COL_COV_O,
                      cov_ont or {}, nc_cov_ont or {}, {},
                      telo_ont, telo_nc_ont,
                      gaps_=gaps_ont, gaps_asm_=gaps_ont_asm)

    if show_chrom_labels:
        chrom_ticks = [(y_label[tok], tok) for tok in group_order if tok in y_label]
        ax.set_yticks([yc for yc, _ in chrom_ticks])
        _lbl_remap = {'W': 'W/Z'}
        ax.set_yticklabels([_lbl_remap.get(lbl, lbl) for _, lbl in chrom_ticks],
                           fontsize=st['font_pm'] + 1, color='#333333')
        ax.yaxis.set_tick_params(length=0, pad=label_pad)
        if label_side == 'right':
            ax.yaxis.tick_right()
    else:
        ax.set_yticks([])

    if rasterized:
        ax.set_rasterization_zorder(12)


# =============================================================================
# Panel f: split half-violin — coverage per chromosome by category (mat | pat)
# =============================================================================

def _draw_full_coverage_panel(ax, all_macro_order, all_micro_order, all_nano_order,
                               chrom_groups, cov_summary_s, cov_summary_d,
                               cov_summary_ont=None,
                               BAR_H=0.5, style=None, rasterized=False):
    st = style or STYLES['paper']

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

    x_pos_s = {}
    x_pos_d = {}
    x_pos_o = {}
    x_tick  = {}
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

    bar_specs = [
        (x_pos_s, avg_s,   COL_COV_S_DARK, COL_COV_S),
        (x_pos_d, avg_d,   COL_COV_D_DARK, COL_COV_D),
        (x_pos_o, avg_ont, COL_COV_O_DARK, COL_COV_O),
    ]
    for tok in x_tick:
        for xc_map, avg, col_dark, col_light in bar_specs:
            xc = xc_map.get(tok)
            if xc is None or not avg or tok not in avg:
                continue
            cp, ncp, up = avg[tok]
            ax.bar(xc, cp,  width=bar_w, bottom=0,      align='edge',
                   color=col_dark,  ec='none', zorder=2)
            ax.bar(xc, ncp, width=bar_w, bottom=cp,     align='edge',
                   color=col_light, alpha=0.8, ec='none', zorder=2)
            ax.bar(xc, up,  width=bar_w, bottom=cp+ncp, align='edge',
                   color=COL_UNCOV, ec='none', zorder=2)
            ax.add_patch(mpatches.Rectangle((xc, 0), bar_w, 100,
                         fc='none', ec=COL_BORDER, lw=st['border_lw'], zorder=4))

    ax.set_xticks(list(x_tick.values()))
    ax.set_xticklabels(list(x_tick.keys()),
                       fontsize=6,
                       rotation=0, ha='center')
    ax.tick_params(axis='x', length=0, pad=2)

    for gname, (gx0, gx1) in group_xranges.items():
        x_data_mid = (gx0 + gx1) / 2
        x_frac = (x_data_mid - ax.get_xlim()[0]) / (ax.get_xlim()[1] - ax.get_xlim()[0])
        ax.text(x_frac, 1.04, gname,
                ha='center', va='bottom',
                fontsize=st['font_pm'] + 1, fontweight='bold', color='#333333',
                transform=ax.transAxes, clip_on=False)
        if gx0 > 0:
            ax.axvline(gx0 - group_gap / 2,
                       color='#cccccc', lw=0.8, ls='--', zorder=1)

    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(['0', '25', '50', '75', '100'])
    ax.set_ylabel('Coverage (%)', fontsize=st['font_xlabel'], labelpad=3)
    ax.tick_params(axis='y', labelsize=st['font_tick'], length=2)

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_linewidth(st['border_lw'] * 0.7)
    ax.spines['bottom'].set_linewidth(st['border_lw'] * 0.7)

    if rasterized:
        ax.set_rasterization_zorder(5)


# =============================================================================
# Panel d — split half-violin: coverage per chromosome by category (mat | pat)
# =============================================================================

def _draw_panel_d_violin(axes_by_cat, all_macro_order, all_micro_order, all_nano_order,
                          chrom_groups, cov_summary_s, cov_summary_d,
                          cov_summary_ont=None, style=None):
    """Split half-violin per category, each on its own axes with zoomed y-range.

    axes_by_cat: dict {'Macro': ax, 'Micro': ax, 'Dot': ax}
    Left half = maternal, right half = paternal chromosomes.
    Y-value = collinear coverage only.  Zoomed per group.
    All chromosomes shown as dots; IQR outliers labelled with chromosome token.
    """
    st = style or STYLES['paper']

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

    # Collect per (group, method) → {mat: [(val, tok)], pat: [(val, tok)]}
    # val = collinear coverage only
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

    # Shared x layout: 3 method split-violins per axes
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

    # Pre-compute which Dot tokens to label per (method, side):
    #   ONT_Dual: label '30' and '25' on maternal side only
    #   Single/Dual: label '30' and '25' on whichever haplotype has higher coverage
    _DOT_LABEL_TOKS = ('30', '25')
    _dot_force = {}
    for mname, cov_dict, _ in method_info:
        if mname == 'ONT_Dual':
            _dot_force[(mname, 'mat')] = set(_DOT_LABEL_TOKS)
            _dot_force[(mname, 'pat')] = set()
        else:
            mat_force, pat_force = set(), set()
            for tok_check in _DOT_LABEL_TOKS:
                mat_v = next((v for v, t in vdata[('Dot', mname)]['mat'] if t == tok_check), 0.0)
                pat_v = next((v for v, t in vdata[('Dot', mname)]['pat'] if t == tok_check), 0.0)
                (mat_force if mat_v >= pat_v else pat_force).add(tok_check)
            _dot_force[(mname, 'mat')] = mat_force
            _dot_force[(mname, 'pat')] = pat_force

    for gi, (gname, _) in enumerate(group_specs):
        ax = axes_by_cat.get(gname)
        if ax is None:
            continue

        # Determine y-range from this category's data
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

        def _hv(vals_toks, xc, side, color, force_toks=None,
                _ax=ax, _y_lo=y_lo, _y_hi=y_hi):
            if not vals_toks:
                return
            vals = np.array([v for v, _ in vals_toks], dtype=float)
            toks = [t for _, t in vals_toks]

            bx0 = xc - half_w if side == 'left' else xc
            bx1 = xc          if side == 'left' else xc + half_w
            wx  = (bx0 + bx1) / 2

            # KDE over the zoomed y range
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

            # Half-violin polygon
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
                fc=color, ec=color, lw=0.3, alpha=0.28, zorder=3))
            _ax.plot(outer_x, y_grid, color=color, lw=0.5, alpha=0.75, zorder=4)

            # Median tick
            q2 = float(np.median(vals))
            _ax.plot([bx0, bx1], [q2, q2], color=COL_BORDER, lw=0.8, zorder=5)

            # Scatter with jitter; label specific forced tokens
            jitter = np.random.uniform(-(half_w * 0.38), half_w * 0.38, len(vals))
            xpts      = np.clip(wx + jitter, bx0 + 0.005, bx1 - 0.005)
            lbl_set   = force_toks or set()
            lbl_mask  = [tok in lbl_set for tok in toks]

            for xi, yi, is_lbl in zip(xpts, vals, lbl_mask):
                if not is_lbl:
                    _ax.scatter(xi, yi, color=color, s=2,
                                zorder=5, ec='none', alpha=0.55)
            for xi, yi, tok, is_lbl in zip(xpts, vals, toks, lbl_mask):
                if is_lbl:
                    _ax.scatter(xi, yi, color=color, s=7,
                                zorder=7, ec='white', linewidths=0.25)
                    y_off = (_y_hi - _y_lo) * 0.04
                    _ax.text(xi, yi - y_off, tok,
                             ha='center', va='top',
                             fontsize=st['font_tick'] - 3,
                             color=color, zorder=8)

        for mname, _, color in method_info:
            xc    = meth_centers[mname]
            sides = vdata[(gname, mname)]
            if gname == 'Macro':
                mat_force = {'Z', 'W', '4'}
                pat_force = {'Z', 'W', '4'}
            elif gname == 'Micro':
                mat_force = {'22'}
                pat_force = {'22'}
            elif gname == 'Dot':
                mat_force = _dot_force[(mname, 'mat')]
                pat_force = _dot_force[(mname, 'pat')]
            else:
                mat_force = set()
                pat_force = set()
            _hv(sides['mat'], xc, 'left',  color, force_toks=mat_force)
            _hv(sides['pat'], xc, 'right', color, force_toks=pat_force)

        # X-tick: category name centred
        ax.set_xticks([x_max / 2])
        ax.set_xticklabels([gname], fontsize=st['font_tick'] - 1)
        ax.tick_params(axis='x', length=0, pad=3)

        # Y-axis: auto-ticks based on zoomed range; label only on first panel
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=4, integer=True))
        ax.tick_params(axis='y', labelsize=st['font_tick'] - 1, length=2)
        if gi == 0:
            ax.set_ylabel('Collinear coverage (%)',
                          fontsize=st['font_tick'] - 1, labelpad=2)
            ax.text(0.02, 1.01, 'Mat◀|▶Pat',
                    transform=ax.transAxes,
                    fontsize=st['font_tick'] - 3, color='#555555',
                    ha='left', va='bottom')
        else:
            ax.yaxis.set_tick_params(labelleft=True)   # keep ticks, show labels

        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        ax.spines['left'].set_linewidth(st['border_lw'] * 0.5)
        ax.spines['bottom'].set_linewidth(st['border_lw'] * 0.5)


# =============================================================================
# Panel e — stacked bar (coverage fractions per category)
# =============================================================================

def _draw_panel_e(ax, df, st):
    bw  = st['ef_bar_w']
    gap = st['ef_group_gap']
    lw  = st['border_lw']

    cat_x   = {cat: i * gap for i, cat in enumerate(_EF_CAT_ORDER)}
    n_asm   = len(_EF_ASM_ORDER)
    offsets = np.linspace(-(n_asm - 1) / 2 * bw, (n_asm - 1) / 2 * bw, n_asm)

    for ai, asm in enumerate(_EF_ASM_ORDER):
        color    = _EF_ASM_COLORS[asm]
        nc_color = _EF_ASM_NC_COLORS[asm]
        sub = df[(df['assembly'] == asm) & (df['category'].isin(_EF_CAT_ORDER))]
        sub = sub.set_index('category')

        for cat in _EF_CAT_ORDER:
            if cat not in sub.index:
                continue
            x    = cat_x[cat] + offsets[ai]
            row  = sub.loc[cat]
            coll = row['collinear_pct']    / 100.0
            nc   = row['noncollinear_pct'] / 100.0
            sw   = row['switch_pct']       / 100.0

            ax.bar(x, coll, width=bw, bottom=0,
                   color=color, ec=COL_BORDER, lw=lw * 0.4, zorder=3)
            if nc > 0:
                ax.bar(x, nc, width=bw, bottom=coll,
                       color=nc_color, ec=COL_BORDER, lw=lw * 0.4, zorder=3)
            if sw > 1e-6:
                ax.bar(x, sw, width=bw, bottom=0,
                       color=COL_SW, ec=COL_BORDER, lw=lw * 0.4,
                       alpha=0.75, zorder=4)

    ax.set_xlim(-gap * 0.6, (len(_EF_CAT_ORDER) - 1) * gap + gap * 0.6)
    ax.set_ylim(0, 1.05)
    ax.set_yticks([0, 0.50, 1.00])
    ax.set_yticklabels(['0', '50', '100'], fontsize=st['font_tick'] - 1)
    ax.set_ylabel('Coverage (%)', fontsize=st['font_tick'] - 1, labelpad=2)

    xtick_pos = [cat_x[c] for c in _EF_CAT_ORDER]
    ax.set_xticks(xtick_pos)
    ax.set_xticklabels([_EF_CAT_LABELS[c] for c in _EF_CAT_ORDER],
                        fontsize=st['font_tick'] - 1)
    ax.tick_params(axis='x', length=0, pad=3)
    ax.tick_params(axis='y', labelsize=st['font_tick'] - 1, length=2)

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_linewidth(lw * 0.5)
    ax.spines['bottom'].set_linewidth(lw * 0.5)
    ax.set_facecolor('white')


# =============================================================================
# Panel f — lollipop (telomere completeness per category)
# =============================================================================

def _draw_panel_f(ax, df, st):
    gap = st['ef_group_gap']
    off = st['ef_lollipop_offset']
    slw = st['ef_stem_lw']
    ms  = st['ef_dot_ms']
    lw  = st['border_lw']

    cat_x   = {cat: i * gap for i, cat in enumerate(_EF_CAT_ORDER)}
    n_asm   = len(_EF_ASM_ORDER)
    offsets = np.linspace(-(n_asm - 1) / 2 * off, (n_asm - 1) / 2 * off, n_asm)

    # Expected telomere counts per category (same across assemblies)
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

            ax.plot([x, x], [0, frac], color=color,
                    lw=slw, solid_capstyle='round', zorder=3)
            ax.scatter(x, frac, color=color, s=ms ** 2 * 0.5,
                       zorder=4, ec='white', linewidths=0.4)
            ax.text(x, frac + 0.05, f'{n_present}',
                    ha='center', va='bottom', fontsize=st['font_tick'] - 2,
                    color=color, zorder=5)

    ax.axhline(1.0, color='#aaaaaa', lw=0.6, ls='--', zorder=1)

    ax.set_xlim(-gap * 0.6, (len(_EF_CAT_ORDER) - 1) * gap + gap * 0.6)
    ax.set_ylim(0, 1.40)
    ax.set_yticks([0, 0.50, 1.00])
    ax.set_yticklabels(['0', '50', '100'], fontsize=st['font_tick'] - 1)
    ax.set_ylabel('Telomere (%)', fontsize=st['font_tick'] - 1, labelpad=2)

    xtick_pos = [cat_x[c] for c in _EF_CAT_ORDER]
    ax.set_xticks(xtick_pos)
    ax.set_xticklabels(
        [f"{_EF_CAT_LABELS[c]}\n(n={cat_expected.get(c, '?')})" for c in _EF_CAT_ORDER],
        fontsize=st['font_tick'] - 1)
    ax.tick_params(axis='x', length=0, pad=3)
    ax.tick_params(axis='y', labelsize=st['font_tick'] - 1, length=2)

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_linewidth(lw * 0.5)
    ax.spines['bottom'].set_linewidth(lw * 0.5)
    ax.set_facecolor('white')


# =============================================================================
# Panel g — gene-model error rates dot plot
# =============================================================================

_COL_HIFI_G = '#26518F'
_COL_ONT_G  = '#6B8E23'

_G_ROW_SPECS = [
    ('HiFi', 'hap1'),
    ('HiFi', 'hap2'),
    ('ONT',  'hap1'),
    ('ONT',  'hap2'),
]
# Bar geometry: slot=1, bar_h=0.62*slot, gap_within(hap1-hap2)=0.12*slot edge-to-edge,
# gap_between(HiFi-ONT)=0.55*slot edge-to-edge → center-to-center = bar_h + gap.
_G_BAR_H = 0.62
_G_Y = {
    ('HiFi', 'hap1'):  2.65,
    ('HiFi', 'hap2'):  1.91,
    ('ONT',  'hap1'):  0.74,
    ('ONT',  'hap2'):  0.00,
}
_G_YLIM = (-0.46, 3.11)


def _load_annotation_tsv(path):
    df = pd.read_csv(path, sep='\t')
    for col in ['Frameshifts_errors (%)', 'Premature stop codon(%)']:
        df[col] = df[col].astype(str).str.replace(',', '.').astype(float)
    df['technology'] = df['Assembly'].str.split('+').str[0]
    return df


def _clean_xlim(vmin, vmax, n_breaks=4):
    """Return (lo, hi, ticks) with 15 % padding, rounded to a clean interval."""
    pad  = max((vmax - vmin) * 0.15, 1e-6)
    lo   = vmin - pad
    hi   = vmax + pad
    rng  = hi - lo
    step_raw = rng / (n_breaks - 1)
    mag  = 10 ** np.floor(np.log10(step_raw))
    step = np.ceil(step_raw / mag) * mag
    lo   = np.floor(lo / step) * step
    hi   = np.ceil(hi  / step) * step
    ticks = np.round(np.arange(lo, hi + step * 0.01, step), 10)
    return float(lo), float(hi), ticks.tolist()


def _draw_panel_g(ax_left, ax_right, ax_pcg, annotation_df, style=None):
    """Two vertical bar subpanels (Frameshifts | Premature stop) + PCG text column.
    B&W encoding: HiFi = solid black; ONT = white fill + diagonal hatch + black edge.
    """
    st       = style or STYLES['paper']
    fs_tick  = st['font_tick'] - 1   # 7 pt
    fs_title = st['font_tick']       # 8 pt
    bar_w    = _G_BAR_H              # 0.62  (bar width in x-data units)

    # x-centres of the four vertical bars (same spacing as old y-positions)
    _GX = {
        ('HiFi', 'hap1'): 0.00,
        ('HiFi', 'hap2'): 0.74,   # 0.62 bar + 0.12 within-group gap
        ('ONT',  'hap1'): 1.91,   # + 0.62 bar + 0.55 between-group gap
        ('ONT',  'hap2'): 2.65,   # + 0.62 bar + 0.12 within-group gap
    }
    xlim_bars = _G_YLIM            # (-0.46, 3.11) — same span as before

    # B&W fill styles: HiFi solid black, ONT white + hatch
    _FC = {
        'HiFi': dict(facecolor='black', edgecolor='black', lw=0.5, hatch=None),
        'ONT':  dict(facecolor='white', edgecolor='black', lw=0.7, hatch='////'),
    }

    # ── bar subpanels ──────────────────────────────────────────────────────────
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
        v_max    = max(vals)
        ylim_max = float(int(np.ceil(v_max * 1.15 / tick_step)) * tick_step)
        yticks   = list(np.round(np.arange(0, ylim_max + tick_step * 0.5, tick_step), 10))

        ax.set_xlim(*xlim_bars)
        ax.set_ylim(0, ylim_max)
        ax.set_xticks([])
        ax.set_yticks(yticks)
        if tick_step < 1:
            ax.yaxis.set_major_formatter(
                mticker.FuncFormatter(lambda y, _: f'{y:.1f}'))
        else:
            ax.yaxis.set_major_formatter(
                mticker.FuncFormatter(lambda y, _: f'{y:.0f}'))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(0.4)
        ax.spines['bottom'].set_linewidth(0.4)
        ax.tick_params(axis='y', length=2.5, width=0.4, labelsize=fs_tick, pad=2)
        ax.set_facecolor('none')

        # Subpanel title with hairline rule
        ax.set_title(title, fontsize=fs_title, pad=3, fontweight='normal', loc='center')
        ax.plot([0, 1], [1, 1], transform=ax.transAxes,
                color='#333333', lw=0.25, clip_on=False, solid_capstyle='butt')

        for tech, hap in _G_ROW_SPECS:
            row = annotation_df.loc[
                (annotation_df['technology'] == tech) &
                (annotation_df['Haplotype']  == hap)
            ]
            if row.empty:
                continue
            v  = float(row[col].iloc[0])
            x  = _GX[(tech, hap)]
            fc = _FC[tech]
            ax.bar(x, v, width=bar_w, bottom=0,
                   facecolor=fc['facecolor'], edgecolor=fc['edgecolor'],
                   linewidth=fc['lw'], hatch=fc['hatch'],
                   zorder=4, align='center')
            ax.text(x, v + ylim_max * 0.025, f'{v:.2f}',
                    ha='center', va='bottom', fontsize=fs_tick,
                    color='#333333', clip_on=False)

        # x-tick labels: Hap 1 / Hap 2 below each bar
        ax.set_xticks([_GX[(t, h)] for t, h in _G_ROW_SPECS])
        ax.set_xticklabels(['Hap 1', 'Hap 2', 'Hap 1', 'Hap 2'], fontsize=fs_tick)
        ax.tick_params(axis='x', length=0, pad=2)

        # HiFi / ONT group labels with bracket line below the Hap labels
        lo, hi = xlim_bars
        span   = hi - lo
        for tech, hap_pair in [('HiFi', ['hap1', 'hap2']), ('ONT', ['hap1', 'hap2'])]:
            xs     = [_GX[(tech, h)] for h in hap_pair]
            a_l    = (min(xs) - bar_w / 2 - lo) / span   # axes fraction
            a_r    = (max(xs) + bar_w / 2 - lo) / span
            a_mid  = (a_l + a_r) / 2
            y_spi  = -0.22   # bracket spine y in axes coords
            y_arm  = -0.19   # arm tips (shorter than spine — arms point up)
            # horizontal spine
            ax.plot([a_l, a_r], [y_spi, y_spi], transform=ax.transAxes,
                    color='#333333', lw=0.5, clip_on=False, solid_capstyle='butt')
            # left arm up
            ax.plot([a_l, a_l], [y_spi, y_arm], transform=ax.transAxes,
                    color='#333333', lw=0.5, clip_on=False, solid_capstyle='butt')
            # right arm up
            ax.plot([a_r, a_r], [y_spi, y_arm], transform=ax.transAxes,
                    color='#333333', lw=0.5, clip_on=False, solid_capstyle='butt')
            ax.text(a_mid, y_spi - 0.06, tech, transform=ax.transAxes,
                    ha='center', va='top', fontsize=fs_title,
                    fontweight='bold', color='#333333', clip_on=False)

    # ── PCG text column ────────────────────────────────────────────────────────
    # Four rows top-to-bottom: HiFi Hap1, HiFi Hap2, ONT Hap1, ONT Hap2
    _PCG_Y = {
        ('HiFi', 'hap1'): 3.0,
        ('HiFi', 'hap2'): 2.0,
        ('ONT',  'hap1'): 1.0,
        ('ONT',  'hap2'): 0.0,
    }
    pcg_ylim = (-0.5, 3.5)

    ax_pcg.set_xlim(0, 1)
    ax_pcg.set_ylim(*pcg_ylim)
    ax_pcg.set_xticks([])
    ax_pcg.set_yticks([])
    for sp in ax_pcg.spines.values():
        sp.set_visible(False)
    ax_pcg.set_facecolor('none')

    ax_pcg.set_title('Protein-coding genes', fontsize=fs_title, pad=3,
                     fontweight='normal', loc='right')
    ax_pcg.plot([0, 1], [1, 1], transform=ax_pcg.transAxes,
                color='#333333', lw=0.25, clip_on=False, solid_capstyle='butt')

    trans_pcg = blended_transform_factory(ax_pcg.transAxes, ax_pcg.transData)
    for tech, hap in _G_ROW_SPECS:
        row = annotation_df.loc[
            (annotation_df['technology'] == tech) &
            (annotation_df['Haplotype']  == hap)
        ]
        if row.empty:
            continue
        n = int(row['Nr.protein coding genes'].iloc[0])
        y = _PCG_Y[(tech, hap)]
        ax_pcg.text(0.95, y, f'{n:,}',
                    transform=trans_pcg, ha='right', va='center',
                    fontsize=fs_tick, color='#333333', clip_on=False)


# =============================================================================
# Custom legend band
# =============================================================================

def _draw_legend_band(fig, annotation_df=None, style=None):
    """
    Horizontal legend band matching the standalone SVG spec (400 × 60 pt canvas).

    Coordinate mapping: SVG x / 400 → matplotlib x [0,1];
                        1 − SVG y / 60 → matplotlib y [0,1] (SVG y-down → mpl y-up).

    Named SVG groups (via gid): swatch-matrix, platform-brackets,
    platform-markers, feature-glyphs.
    """
    st    = style or STYLES['paper']
    fsize = st['font_legend']

    # ── Legend axes ────────────────────────────────────────────────────────
    lax_l, lax_b, lax_w, lax_h = 0.04, 0.030, 0.92, 0.082
    lax = fig.add_axes([lax_l, lax_b, lax_w, lax_h],
                       label='legend-band', facecolor='white', zorder=8)
    lax.set_xlim(0, 1)
    lax.set_ylim(0, 1)
    lax.axis('off')

    fig_w_in, fig_h_in = fig.get_size_inches()
    ax_w_in = lax_w * fig_w_in
    ax_h_in = lax_h * fig_h_in

    # Coordinate converters (SVG canvas: 420 × 60 pt)
    W_CANVAS = 420.0
    def sx(x): return x / W_CANVAS         # SVG x → mpl x
    def sy(y): return 1.0 - y / 60.0       # SVG y (down) → mpl y (up)
    def sw(w): return w / W_CANVAS         # SVG width → mpl x-width
    def sh(h): return h / 60.0             # SVG height → mpl y-height (magnitude)

    # For a physical pt size → mpl data unit (aspect-correct circles/discs)
    def ptx(pt): return pt / 72.0 / ax_w_in
    def pty(pt): return pt / 72.0 / ax_h_in

    # Row y-centres from SVG spec (y = 12, 30, 48 of 60)
    ry = [sy(12), sy(30), sy(48)]   # [0.800, 0.500, 0.200]

    # Swatch: 9 × 9 pt rectangles, no stroke
    swatch_xw = sw(9)    # x-width  = 9 / W_CANVAS
    swatch_yh = sh(9)    # y-height = 9 / 60

    # ══════════════════════════════════════════════════════════════════════
    # GROUP — swatch-matrix   (col A: collinear; col B: non-collinear)
    # ══════════════════════════════════════════════════════════════════════
    col_a_sw  = sx(6)
    col_a_lbl = sx(18)
    col_b_sw  = sx(90)
    col_b_lbl = sx(102)

    g1_rows = [
        (COL_COV_S_DARK, 'Single — collinear',    COL_COV_S, 'Single — non-collinear'),
        (COL_COV_D_DARK, 'Dual — collinear',      COL_COV_D, 'Dual — non-collinear'),
        (COL_COV_O_DARK, 'ONT Dual — collinear',  COL_COV_O, 'ONT Dual — non-collinear'),
    ]

    for i, (c_c, lbl_c, c_n, lbl_n) in enumerate(g1_rows):
        y_ctr = ry[i]
        y_bot = y_ctr - swatch_yh / 2   # bottom of 9pt rect in mpl coords

        lax.add_patch(mpatches.Rectangle(
            (col_a_sw, y_bot), swatch_xw, swatch_yh,
            fc=c_c, ec='none', transform=lax.transData,
            gid=f'swatch-matrix-col-sw-{i+1}'))
        lax.text(col_a_lbl, y_ctr, lbl_c,
                 va='center', ha='left', fontsize=fsize,
                 transform=lax.transData, gid=f'swatch-matrix-col-lbl-{i+1}')

        lax.add_patch(mpatches.Rectangle(
            (col_b_sw, y_bot), swatch_xw, swatch_yh,
            fc=c_n, ec='none', transform=lax.transData,
            gid=f'swatch-matrix-nc-sw-{i+1}'))
        lax.text(col_b_lbl, y_ctr, lbl_n,
                 va='center', ha='left', fontsize=fsize,
                 transform=lax.transData, gid=f'swatch-matrix-nc-lbl-{i+1}')

    # ══════════════════════════════════════════════════════════════════════
    # GROUP — platform-brackets
    # ══════════════════════════════════════════════════════════════════════
    brk_arm = sx(184)
    brk_spi = sx(188)

    # Upper bracket  (mpl y: high → low because SVG y=5 is near top)
    lax.plot([brk_arm, brk_spi, brk_spi, brk_arm],
             [sy(5),   sy(5),   sy(36.5), sy(36.5)],
             color='#000000', lw=0.5, solid_capstyle='butt',
             transform=lax.transData, gid='platform-brackets-hifi', clip_on=False)

    # Lower bracket
    lax.plot([brk_arm, brk_spi, brk_spi, brk_arm],
             [sy(41.5), sy(41.5), sy(55), sy(55)],
             color='#000000', lw=0.5, solid_capstyle='butt',
             transform=lax.transData, gid='platform-brackets-ont', clip_on=False)

    # ══════════════════════════════════════════════════════════════════════
    # GROUP — platform-markers
    # ══════════════════════════════════════════════════════════════════════
    g2_cx   = sx(197)
    lbl2_x  = sx(205)
    dot_s   = np.pi * 4.0 ** 2

    # HiFi: midpoint of upper bracket (SVG cy=20.75); ONT: lower (cy=48.25)
    for fc, ec, lbl, cy_svg in [
        (_COL_HIFI_G, 'none',    'HiFi', 20.75),
        (_COL_ONT_G,  '#4A4A10', 'ONT',  48.25),
    ]:
        y = sy(cy_svg)
        lax.scatter([g2_cx], [y], s=dot_s, color=fc, ec=ec, linewidths=0.75,
                    transform=lax.transData, clip_on=False, zorder=5,
                    gid=f'platform-markers-{lbl.lower()}')
        lax.text(lbl2_x, y, lbl,
                 va='center', ha='left', fontsize=fsize,
                 transform=lax.transData, gid=f'platform-markers-{lbl.lower()}-lbl')

    # ── Divider — full-height rule, SVG x=230 ─────────────────────────────
    lax.plot([sx(230), sx(230)], [sy(4), sy(56)], '-',
             color='#000000', lw=0.5, solid_capstyle='butt',
             transform=lax.transData, gid='divider', clip_on=False)

    # ══════════════════════════════════════════════════════════════════════
    # GROUP — feature-glyphs (2 columns: switch/gaps & telomeres)
    # ══════════════════════════════════════════════════════════════════════
    col1_sw_x  = sx(240)
    col1_lbl_x = sx(249)

    # Row 1 (y=12): Hap. switch — 5 pt wide × 9 pt tall orange rect
    lax.add_patch(mpatches.Rectangle(
        (col1_sw_x - sw(2.5), sy(16.5)), sw(5), sh(9),
        fc=COL_SW, ec='none', transform=lax.transData,
        gid='feature-glyphs-hapswitch'))
    lax.text(col1_lbl_x, ry[0], 'Hap. switch',
             va='center', ha='left', fontsize=fsize,
             transform=lax.transData, gid='feature-glyphs-hapswitch-lbl')

    # Row 2 (y=30): Assembly gap — vertical tick mark
    lax.plot([col1_sw_x, col1_sw_x], [sy(30) - sh(4.5), sy(30) + sh(4.5)],
             color=COL_ASM_GAP, lw=1.2, solid_capstyle='butt',
             transform=lax.transData, gid='feature-glyphs-asmgap')
    lax.text(col1_lbl_x, ry[1], 'Assembly gap',
             va='center', ha='left', fontsize=fsize,
             transform=lax.transData, gid='feature-glyphs-asmgap-lbl')

    # Row 3 (y=48): Curation gap — red rectangle
    lax.add_patch(mpatches.Rectangle(
        (col1_sw_x - sw(2.5), sy(52.5)), sw(5), sh(9),
        fc=COL_GAP, ec='none', alpha=0.9, transform=lax.transData,
        gid='feature-glyphs-curgap'))
    lax.text(col1_lbl_x, ry[2], 'Curation gap',
             va='center', ha='left', fontsize=fsize,
             transform=lax.transData, gid='feature-glyphs-curgap-lbl')

    # Column 2 (x=314): Telomeres
    col2_flat_x = sx(314)
    col2_lbl_x  = sx(322)
    R_in        = 4.5 / 72.0          # 4.5 pt radius in inches
    telo_rx     = R_in / ax_w_in      # x-radius in data units
    telo_ry     = R_in / ax_h_in      # y-radius in data units
    theta       = np.linspace(np.pi / 2, 3 * np.pi / 2, 60)   # left-half arc

    for j, (filled, lbl, gid_pre) in enumerate([
        (True,  'Telomere — collinear',     'feature-glyphs-telo-col'),
        (False, 'Telomere — non-collinear', 'feature-glyphs-telo-noncol'),
    ]):
        y     = ry[j]
        arc_x = col2_flat_x + telo_rx * np.cos(theta)   # bulges left of flat_x
        arc_y = y           + telo_ry * np.sin(theta)
        disc  = mpatches.Polygon(list(zip(arc_x, arc_y)), closed=True,
                                  fc=COL_TELO if filled else 'none',
                                  ec=COL_TELO, lw=0.75,
                                  transform=lax.transData, gid=f'{gid_pre}-glyph')
        lax.add_patch(disc)
        lax.text(col2_lbl_x, y, lbl,
                 va='center', ha='left', fontsize=fsize,
                 transform=lax.transData, gid=f'{gid_pre}-lbl')

    # ── Invisible bounding rects → SVG group anchors for Illustrator ──────
    for gid, x0, x1 in [
        ('legend-swatch-matrix',     0.0,       brk_arm),
        ('legend-platform-brackets', brk_arm,   sx(230)),
        ('legend-feature-glyphs',    sx(230),   1.0),
    ]:
        lax.add_patch(mpatches.Rectangle(
            (x0, 0.0), x1 - x0, 1.0,
            fc='none', ec='none', lw=0, alpha=0, zorder=0,
            transform=lax.transData, gid=gid))

    # ── Border box around legend ───────────────────────────────────────────
    lax.add_patch(mpatches.Rectangle(
        (0.0, 0.0), 1.0, 1.0,
        fc='none', ec='#888888', lw=0.6,
        transform=lax.transAxes, clip_on=False, zorder=11,
        gid='legend-border'))

    return lax


# =============================================================================
# Main figure builder
# =============================================================================

def build_combined_figure(df_s, df_d,
                           cov_s, cov_d, nc_cov_s, nc_cov_d,
                           switch_s, switch_d,
                           telo_s, telo_d, telo_nc_s, telo_nc_d,
                           cov_summary_s, cov_summary_d, cov_summary_ont,
                           cov_ont, nc_cov_ont,
                           telo_ont, telo_nc_ont,
                           df_stats,
                           out_path, dpi, fmt,
                           style_name='paper',
                           rasterized=False,
                           simplify=False,
                           centromeres=None,
                           flip_set=None,
                           annotation_df=None):

    # ---- Chromosome ordering ------------------------------------------------
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

    all_macro_order = [t for t in chrom_order if classify_chrom_group(t) == 'macro']
    all_micro_order = [t for t in chrom_order if classify_chrom_group(t) == 'micro']
    all_nano_order  = [t for t in chrom_order if classify_chrom_group(t) == 'nano']

    macro_order = [t for t in all_macro_order if t in MACRO_SELECTED]
    micro_order = [t for t in all_micro_order if t in MICRO_SELECTED]
    nano_order  = [t for t in all_nano_order  if t in NANO_SELECTED]

    def size_of(nm):
        return int(df_s.loc[nm, 'size'] if nm in df_s.index else df_d.loc[nm, 'size'])

    # ---- Style & layout params ----------------------------------------------
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
    HAP_GAP   = BAR_H * 0.14
    GROUP_GAP = BAR_H * 0.52
    is_mat    = lambda nm: not is_pat(nm)

    def _hap_names(grp_order, filt):
        return [nm for tok in grp_order
                for nm in chrom_groups.get(tok, []) if filt(nm)]

    # ---- Triple y-layouts per panel group -----------------------------------
    yps_mac_mat, ypd_mac_mat, ypo_mac_mat, yl_mac_mat, ys_mac_mat, _ty_mac_mat = \
        make_y_layout_triple(macro_order, chrom_groups, is_mat, BAR_H, HAP_GAP, GROUP_GAP)
    yps_mac_pat, ypd_mac_pat, ypo_mac_pat, yl_mac_pat, ys_mac_pat, _ty_mac_pat = \
        make_y_layout_triple(macro_order, chrom_groups, is_pat, BAR_H, HAP_GAP, GROUP_GAP)
    ty_mac = max(_ty_mac_mat, _ty_mac_pat)

    yps_mic_mat, ypd_mic_mat, ypo_mic_mat, yl_mic_mat, ys_mic_mat, _ty_mic_mat = \
        make_y_layout_triple(micro_order, chrom_groups, is_mat, BAR_H, HAP_GAP, GROUP_GAP)
    yps_mic_pat, ypd_mic_pat, ypo_mic_pat, yl_mic_pat, ys_mic_pat, _ty_mic_pat = \
        make_y_layout_triple(micro_order, chrom_groups, is_pat, BAR_H, HAP_GAP, GROUP_GAP)
    ty_mic = max(_ty_mic_mat, _ty_mic_pat)

    yps_nan_mat, ypd_nan_mat, ypo_nan_mat, yl_nan_mat, ys_nan_mat, _ty_nan_mat = \
        make_y_layout_triple(nano_order,  chrom_groups, is_mat, BAR_H, HAP_GAP, GROUP_GAP)
    yps_nan_pat, ypd_nan_pat, ypo_nan_pat, yl_nan_pat, ys_nan_pat, _ty_nan_pat = \
        make_y_layout_triple(nano_order,  chrom_groups, is_pat, BAR_H, HAP_GAP, GROUP_GAP)
    ty_nan = max(_ty_nan_mat, _ty_nan_pat)

    margin_y  = GROUP_GAP * 0.6

    names_mac_mat = _hap_names(macro_order, is_mat)
    names_mac_pat = _hap_names(macro_order, is_pat)
    names_mic_mat = _hap_names(micro_order, is_mat)
    names_mic_pat = _hap_names(micro_order, is_pat)
    names_nan_mat = _hap_names(nano_order,  is_mat)
    names_nan_pat = _hap_names(nano_order,  is_pat)

    def _max_size(names):
        return max((size_of(nm) for nm in names), default=1)

    max_mac = _max_size(names_mac_mat + names_mac_pat)
    max_mic = _max_size(names_mic_mat + names_mic_pat)
    max_nan = _max_size(names_nan_mat + names_nan_pat)

    has_avg    = (cov_summary_s is not None) or (cov_summary_d is not None) or \
                 (cov_summary_ont is not None)
    has_panelg = annotation_df is not None
    cov_panel_h = max(ty_mac, ty_mic, ty_nan) * 1.6

    # ---- Figure: 4-5 rows (A+B+C / D / E / F [/ G]) — stats row first ---------
    fig_w = st['fig_w']
    if has_avg:
        base_ratios = [cov_panel_h, ty_mac, ty_mic, ty_nan]
    else:
        base_ratios = [ty_mac, ty_mic, ty_nan]

    height_ratios = base_ratios + ([ty_mic] if has_panelg else [])
    n_outer_rows  = len(height_ratios)

    # Scale figure height proportionally to accommodate panel g
    total_base = sum(base_ratios)
    total_new  = sum(height_ratios)
    fig_h = st['fig_h'] * (total_new / total_base)

    fig   = plt.figure(figsize=(fig_w, fig_h), facecolor='white')

    outer = GridSpec(n_outer_rows, 1, figure=fig,
                     height_ratios=height_ratios,
                     hspace=1.10,
                     left=0.08, right=0.97, top=0.93, bottom=0.14)

    # Row 0 — panels a/b/c (violin + lollipop + coverage) when has_avg
    ax_cov = ax_E = ax_F = None
    ax_viol = {}   # {'Macro': ax, 'Micro': ax, 'Dot': ax}
    if has_avg:
        inner_d = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0, 0],
                                          width_ratios=[0.36, 0.64], wspace=0.12)
        ax_cov  = fig.add_subplot(inner_d[0, 1])   # panel c — coverage on RIGHT
        if df_stats is not None:
            inner_de = GridSpecFromSubplotSpec(2, 1, subplot_spec=inner_d[0, 0],
                                               height_ratios=[3.0, 0.9], hspace=0.42)
            # Panel a: 3 sub-axes side-by-side, one per category
            inner_viol = GridSpecFromSubplotSpec(1, 3, subplot_spec=inner_de[0],
                                                 wspace=0.45)
            ax_viol = {
                'Macro': fig.add_subplot(inner_viol[0, 0]),
                'Micro': fig.add_subplot(inner_viol[0, 1]),
                'Dot':   fig.add_subplot(inner_viol[0, 2]),
            }
            ax_E = ax_viol['Macro']   # used for label placement
            ax_F = fig.add_subplot(inner_de[1])    # panel b — lollipop

    # Row 1 (or 0 when no stats) — Panel D: Macrochromosomes, nested 2-col (mat | pat)
    _mac_row = 1 if has_avg else 0
    inner_a    = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[_mac_row, 0], wspace=0.05)
    ax_mac_mat = fig.add_subplot(inner_a[0, 0])
    ax_mac_pat = fig.add_subplot(inner_a[0, 1])

    # Row 2 (or 1) — Panel E: Microchromosomes, nested 2-col
    _mic_row = 2 if has_avg else 1
    inner_b    = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[_mic_row, 0], wspace=0.05)
    ax_mic_mat = fig.add_subplot(inner_b[0, 0])
    ax_mic_pat = fig.add_subplot(inner_b[0, 1])

    # Row 3 (or 2) — Panel F: Dot chromosomes, nested 2-col
    _nan_row = 3 if has_avg else 2
    inner_c    = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[_nan_row, 0], wspace=0.05)
    ax_nan_mat = fig.add_subplot(inner_c[0, 0])
    ax_nan_pat = fig.add_subplot(inner_c[0, 1])

    # Row 4 (or 3) — Panel G: Gene-model error dot plot
    # 4-column layout: [label-space | left-plot | right-plot | pcg-space]
    # Columns 0 and 3 are reserved blank space (no axes created) for the
    # left-side row labels and the right-side protein-coding-genes column.
    ax_g_left = ax_g_right = ax_g_pcg = None
    if has_panelg:
        _g_row  = _nan_row + 1
        # 6-column layout aligned with panels d/e/f (outer left=0.08, right=0.97).
        # Width budget: gutter 10% | Frameshifts 37% | gap 2% | stop 37% | gap 2% | PCG 12%
        inner_g = GridSpecFromSubplotSpec(
            1, 6, subplot_spec=outer[_g_row, 0],
            width_ratios=[10, 37, 2, 37, 2, 12],
            wspace=0.0)
        ax_g_left  = fig.add_subplot(inner_g[0, 1])   # Frameshifts
        ax_g_right = fig.add_subplot(inner_g[0, 3])   # Premature stop codons
        ax_g_pcg   = fig.add_subplot(inner_g[0, 5])   # Protein-coding genes

    # ---- Common ONT kwargs --------------------------------------------------
    _common_ont = dict(
        yp_ont=None,
        cov_ont=cov_ont,
        nc_cov_ont=nc_cov_ont,
        telo_ont=telo_ont,
        telo_nc_ont=telo_nc_ont,
    )

    # ---- Draw Panel D — maternal (W, mirrored) ------------------------------
    _draw_haplotype_panel(ax_mac_mat, fig_w, fig_h,
                          names_mac_mat,
                          yps_mac_mat, ypd_mac_mat, ys_mac_mat, yl_mac_mat,
                          macro_order, ty_mac, max_mac * 0.008, margin_y,
                          max_mac, size_of,
                          cov_s, nc_cov_s, switch_s,
                          cov_d, nc_cov_d, switch_d,
                          telo_s, telo_nc_s, telo_d, telo_nc_d,
                          BAR_H, GROUP_GAP,
                          show_chrom_labels=True, label_side='left', label_pad=-5,
                          show_xlabel=True, mirror_x=True, show_sd=True,
                          centromeres=centromeres, flip_set=flip_set,
                          style=st, rasterized=rasterized, simplify=simplify,
                          **{**_common_ont, 'yp_ont': ypo_mac_mat})

    # ---- Draw Panel D — paternal (Z) ----------------------------------------
    _draw_haplotype_panel(ax_mac_pat, fig_w, fig_h,
                          names_mac_pat,
                          yps_mac_pat, ypd_mac_pat, ys_mac_pat, yl_mac_pat,
                          macro_order, ty_mac, max_mac * 0.008, margin_y,
                          max_mac, size_of,
                          cov_s, nc_cov_s, switch_s,
                          cov_d, nc_cov_d, switch_d,
                          telo_s, telo_nc_s, telo_d, telo_nc_d,
                          BAR_H, GROUP_GAP,
                          show_chrom_labels=False, label_side='right',
                          show_xlabel=True, mirror_x=False, show_sd=False,
                          centromeres=centromeres, flip_set=flip_set,
                          style=st, rasterized=rasterized, simplify=simplify,
                          **{**_common_ont, 'yp_ont': ypo_mac_pat})

    # ---- Draw Panel B — maternal (mirrored) ---------------------------------
    _draw_haplotype_panel(ax_mic_mat, fig_w, fig_h,
                          names_mic_mat,
                          yps_mic_mat, ypd_mic_mat, ys_mic_mat, yl_mic_mat,
                          micro_order, ty_mic, max_mic * 0.008, margin_y,
                          max_mic, size_of,
                          cov_s, nc_cov_s, switch_s,
                          cov_d, nc_cov_d, switch_d,
                          telo_s, telo_nc_s, telo_d, telo_nc_d,
                          BAR_H, GROUP_GAP,
                          show_chrom_labels=True, label_side='left',
                          show_xlabel=True, mirror_x=True, show_sd=True,
                          centromeres=centromeres, flip_set=flip_set,
                          style=st, rasterized=rasterized, simplify=simplify,
                          **{**_common_ont, 'yp_ont': ypo_mic_mat})

    # ---- Draw Panel B — paternal --------------------------------------------
    _draw_haplotype_panel(ax_mic_pat, fig_w, fig_h,
                          names_mic_pat,
                          yps_mic_pat, ypd_mic_pat, ys_mic_pat, yl_mic_pat,
                          micro_order, ty_mic, max_mic * 0.008, margin_y,
                          max_mic, size_of,
                          cov_s, nc_cov_s, switch_s,
                          cov_d, nc_cov_d, switch_d,
                          telo_s, telo_nc_s, telo_d, telo_nc_d,
                          BAR_H, GROUP_GAP,
                          show_chrom_labels=False, label_side='right',
                          show_xlabel=True, mirror_x=False, show_sd=False,
                          centromeres=centromeres, flip_set=flip_set,
                          style=st, rasterized=rasterized, simplify=simplify,
                          **{**_common_ont, 'yp_ont': ypo_mic_pat})

    # ---- Draw Panel C — maternal (mirrored) ---------------------------------
    _draw_haplotype_panel(ax_nan_mat, fig_w, fig_h,
                          names_nan_mat,
                          yps_nan_mat, ypd_nan_mat, ys_nan_mat, yl_nan_mat,
                          nano_order, ty_nan, max_nan * 0.008, margin_y,
                          max_nan, size_of,
                          cov_s, nc_cov_s, switch_s,
                          cov_d, nc_cov_d, switch_d,
                          telo_s, telo_nc_s, telo_d, telo_nc_d,
                          BAR_H, GROUP_GAP,
                          show_chrom_labels=True, label_side='left',
                          show_xlabel=True, mirror_x=True, show_sd=True,
                          centromeres=centromeres, flip_set=flip_set,
                          style=st, rasterized=rasterized, simplify=simplify,
                          **{**_common_ont, 'yp_ont': ypo_nan_mat})

    # ---- Draw Panel C — paternal --------------------------------------------
    _draw_haplotype_panel(ax_nan_pat, fig_w, fig_h,
                          names_nan_pat,
                          yps_nan_pat, ypd_nan_pat, ys_nan_pat, yl_nan_pat,
                          nano_order, ty_nan, max_nan * 0.008, margin_y,
                          max_nan, size_of,
                          cov_s, nc_cov_s, switch_s,
                          cov_d, nc_cov_d, switch_d,
                          telo_s, telo_nc_s, telo_d, telo_nc_d,
                          BAR_H, GROUP_GAP,
                          show_chrom_labels=False, label_side='right',
                          show_xlabel=True, mirror_x=False, show_sd=False,
                          centromeres=centromeres, flip_set=flip_set,
                          style=st, rasterized=rasterized, simplify=simplify,
                          **{**_common_ont, 'yp_ont': ypo_nan_pat})

    # ---- Draw Panel f (coverage per chromosome) -----------------------------
    if has_avg and ax_cov is not None:
        _draw_full_coverage_panel(ax_cov,
                                   all_macro_order, all_micro_order, all_nano_order,
                                   chrom_groups,
                                   cov_summary_s, cov_summary_d, cov_summary_ont,
                                   BAR_H=BAR_H, style=st, rasterized=rasterized)

    # ---- Draw Panel d (half-violin per category) and e (lollipop) ----------
    if ax_viol:
        _draw_panel_d_violin(ax_viol,
                              all_macro_order, all_micro_order, all_nano_order,
                              chrom_groups,
                              cov_summary_s, cov_summary_d, cov_summary_ont,
                              style=st)
    if ax_F is not None and df_stats is not None:
        _draw_panel_f(ax_F, df_stats, st)

    if has_panelg and ax_g_left is not None:
        _draw_panel_g(ax_g_left, ax_g_right, ax_g_pcg, annotation_df, style=st)

    # ---- Panel labels -------------------------------------------------------
    lw_under = st['border_lw'] * 1.8
    label_kw = dict(transform=fig.transFigure,
                    fontsize=st['font_title'] + 2, fontweight='bold',
                    va='bottom', ha='left', color='black')

    title_fs = st['font_title'] * 0.85
    ax_mac_mat.set_title('Maternal', fontsize=title_fs, color='black',
                          fontweight='bold', pad=4)
    ax_mac_pat.set_title('Paternal', fontsize=title_fs, color='black',
                          fontweight='bold', pad=4)
    ax_mic_mat.set_title('Maternal', fontsize=title_fs, color='black',
                          fontweight='bold', pad=4)
    ax_mic_pat.set_title('Paternal', fontsize=title_fs, color='black',
                          fontweight='bold', pad=4)
    ax_nan_mat.set_title('Maternal', fontsize=title_fs, color='black',
                          fontweight='bold', pad=4)
    ax_nan_pat.set_title('Paternal', fontsize=title_fs, color='black',
                          fontweight='bold', pad=4)

    def _row_header(ax_left, ax_right, label_letter, title_text, header_gap=0.038):
        bb_l = ax_left.get_position()
        bb_r = ax_right.get_position() if ax_right is not None else bb_l
        hy   = bb_l.y1 + header_gap
        fig.text(max(bb_l.x0 - 0.06, 0.01), hy, label_letter, **label_kw)
        cx = (bb_l.x0 + bb_r.x1) / 2
        fig.text(cx, hy, title_text, ha='center', va='bottom',
                 fontsize=st['font_title'], fontweight='bold', color='black',
                 transform=fig.transFigure)
        fig.add_artist(Line2D([bb_l.x0, bb_r.x1],
                               [hy - 0.004, hy - 0.004],
                               color='black', lw=lw_under,
                               transform=fig.transFigure, solid_capstyle='butt'))

    _mac_lbl = 'd' if has_avg else 'a'
    _mic_lbl = 'e' if has_avg else 'b'
    _nan_lbl = 'f' if has_avg else 'c'
    _row_header(ax_mac_mat, ax_mac_pat,  _mac_lbl, 'Macrochromosomes')
    _row_header(ax_mic_mat, ax_mic_pat,  _mic_lbl, 'Microchromosomes')
    _row_header(ax_nan_mat, ax_nan_pat,  _nan_lbl, 'Dot chromosomes')
    if has_panelg and ax_g_left is not None:
        # Panel g header spans the full outer width (0.08–0.97), matching d/e/f.
        # "g" is placed at the same left-margin x as d/e/f labels.
        hy_g  = ax_g_left.get_position().y1 + 0.038
        g_lx  = max(ax_mac_mat.get_position().x0 - 0.06, 0.01)
        fig.text(g_lx, hy_g, 'g', **label_kw)
        fig.text(0.525, hy_g, 'Gene model errors', ha='center', va='bottom',
                 fontsize=st['font_title'], fontweight='bold', color='black',
                 transform=fig.transFigure)
        fig.add_artist(Line2D([0.08, 0.97], [hy_g - 0.004, hy_g - 0.004],
                              color='black', lw=lw_under,
                              transform=fig.transFigure, solid_capstyle='butt'))

    if has_avg and ax_cov:
        bb_f  = ax_cov.get_position()
        hy_f  = bb_f.y1 + 0.045
        fig.text(max(bb_f.x0 - 0.04, 0.01), hy_f, 'c', **label_kw)
        fig.text((bb_f.x0 + bb_f.x1) / 2, hy_f,
                 'Coverage average (2n=80)',
                 ha='center', va='bottom',
                 fontsize=st['font_title'], fontweight='bold', color='black',
                 transform=fig.transFigure)
        fig.add_artist(Line2D([bb_f.x0, bb_f.x1],
                               [hy_f - 0.004, hy_f - 0.004],
                               color='black', lw=lw_under,
                               transform=fig.transFigure, solid_capstyle='butt'))

        if ax_E is not None:
            bb_d = ax_E.get_position()
            bb_e = ax_F.get_position()
            fig.text(max(bb_d.x0 - 0.06, 0.01), bb_d.y1 + 0.038, 'a', **label_kw)
            fig.text(max(bb_e.x0 - 0.06, 0.01), bb_e.y1 + 0.038, 'b', **label_kw)

    # ---- Legend -------------------------------------------------------------
    _draw_legend_band(fig, annotation_df=annotation_df, style=st)

    out_dir = os.path.dirname(os.path.abspath(out_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    base = os.path.splitext(out_path)[0]
    fig.savefig(base + '.png', dpi=dpi)
    print(f'Saved: {base}.png')
    fig.savefig(base + '.pdf', format='pdf')
    print(f'Saved: {base}.pdf')
    fig.savefig(base + '.svg', format='svg')
    print(f'Saved: {base}.svg')
    plt.close(fig)


# =============================================================================
# Entry point
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Combined 6-panel ideogram figure (A–D ideograms + e/f summary stats).',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--single-tsv',               required=True,  metavar='FILE')
    parser.add_argument('--dual-tsv',                 required=True,  metavar='FILE')
    parser.add_argument('--single-chain',             required=True,  metavar='FILE')
    parser.add_argument('--dual-chain',               required=True,  metavar='FILE')
    parser.add_argument('--single-nc-chain',          required=True,  metavar='FILE')
    parser.add_argument('--dual-nc-chain',            required=True,  metavar='FILE')
    parser.add_argument('--single-unlocs-chain',      required=False, metavar='FILE', default=None)
    parser.add_argument('--dual-unlocs-chain',        required=False, metavar='FILE', default=None)
    parser.add_argument('--telo-combined-report',     required=False, metavar='FILE', default=None)
    parser.add_argument('--single-telomere-presence', required=False, metavar='FILE', default=None)
    parser.add_argument('--dual-telomere-presence',   required=False, metavar='FILE', default=None)
    parser.add_argument('--ont-telomere-presence',    required=False, metavar='FILE', default=None)
    parser.add_argument('--dual-pairs',               required=True,  metavar='FILE')
    parser.add_argument('--single-pairs',             required=True,  metavar='FILE')
    parser.add_argument('--single-bed',               required=True,  metavar='FILE')
    parser.add_argument('--dual-bed',                 required=True,  metavar='FILE')
    parser.add_argument('--dual-fai',                 required=False, metavar='FILE', default=None)
    parser.add_argument('--single-fai',               required=False, metavar='FILE', default=None)
    parser.add_argument('--ont-dual-chain',           required=True,  metavar='FILE')
    parser.add_argument('--ont-dual-nc-chain',        required=True,  metavar='FILE')
    parser.add_argument('--ont-telo-combined-report', required=False, metavar='FILE', default=None)
    parser.add_argument('--coverage-summary',         required=False, metavar='FILE', default=None)
    parser.add_argument('--single-coverage-summary',  required=False, metavar='FILE', default=None)
    parser.add_argument('--dual-coverage-summary',    required=False, metavar='FILE', default=None)
    parser.add_argument('--ont-coverage-summary',     required=False, metavar='FILE', default=None)
    parser.add_argument('--centromeres',              required=False, metavar='FILE', default=None)
    parser.add_argument('--telo-p-bed',               required=False, metavar='FILE', default=None)
    parser.add_argument('--stats-tsv',                required=False, metavar='FILE', default=None,
                        help='Wide-form TSV from summarize_stats_by_category.py (for panels e/f)')
    parser.add_argument('--annotation-tsv',           required=False, metavar='FILE', default=None,
                        help='Gene-model error rates TSV (for panel g)')
    parser.add_argument('--output',                   required=True,  metavar='FILE')
    parser.add_argument('--no-timestamp',             action='store_true',
                        help='Do not append timestamp to output filename')
    parser.add_argument('--style',     default='paper', choices=['paper', 'poster'])
    parser.add_argument('--format',    default='png',   choices=['png', 'pdf', 'svg'])
    parser.add_argument('--dpi',       type=int, default=None)
    parser.add_argument('--rasterize', action='store_true')
    parser.add_argument('--simplify',  action='store_true')
    args = parser.parse_args()

    out = args.output
    if not out.endswith(f'.{args.format}'):
        out = f'{out}.{args.format}'
    if not args.no_timestamp:
        base, ext = os.path.splitext(out)
        out = f'{base}_{datetime.now().strftime("%Y%m%d_%H%M%S")}{ext}'

    print('Loading TSVs...')
    df_s = load_tsv(args.single_tsv)
    df_d = load_tsv(args.dual_tsv)
    seq_sizes = {nm: int(df_s.loc[nm, 'size'] if nm in df_s.index else df_d.loc[nm, 'size'])
                 for nm in set(df_s.index) | set(df_d.index)}

    df_stats = None
    if args.stats_tsv:
        print('Loading stats TSV (panels e/f)...')
        df_stats = pd.read_csv(args.stats_tsv, sep='\t')

    annotation_df = None
    if args.annotation_tsv:
        print('Loading annotation TSV (panel g)...')
        annotation_df = _load_annotation_tsv(args.annotation_tsv)

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

    print('Computing ONT Dual collinear coverage...')
    cov_ont = compute_covered(args.ont_dual_chain, seq_sizes)

    print('Computing ONT Dual non-collinear coverage...')
    nc_cov_ont = compute_covered(args.ont_dual_nc_chain, seq_sizes)

    print('Building HiFi switch-error lookups...')
    switch_s = build_switch_lookup(args.single_pairs, args.single_bed,
                                   args.single_fai,   args.single_chain)
    switch_d = build_switch_lookup(args.dual_pairs,   args.dual_bed,
                                   args.dual_fai,     args.dual_chain)

    telo_s = telo_d = telo_nc_s = telo_nc_d = None
    if args.single_telomere_presence and os.path.exists(args.single_telomere_presence):
        print('Loading Single telomere presence...')
        telo_s, telo_nc_s = load_telomere_presence_tsv(args.single_telomere_presence)
    if args.dual_telomere_presence and os.path.exists(args.dual_telomere_presence):
        print('Loading Dual telomere presence...')
        telo_d, telo_nc_d = load_telomere_presence_tsv(args.dual_telomere_presence)

    if (telo_s is None or telo_d is None) and args.telo_combined_report:
        print('Loading HiFi telomere report...')
        t_s, t_d, t_nc_s, t_nc_d = load_telomere_ends_from_combined_report(args.telo_combined_report)
        if telo_s is None: telo_s, telo_nc_s = t_s, t_nc_s
        if telo_d is None: telo_d, telo_nc_d = t_d, t_nc_d

    telo_ont = telo_nc_ont = None
    if args.ont_telomere_presence and os.path.exists(args.ont_telomere_presence):
        print('Loading ONT telomere presence...')
        telo_ont, telo_nc_ont = load_telomere_presence_tsv(args.ont_telomere_presence)
    elif args.ont_telo_combined_report:
        print('Loading ONT telomere report...')
        telo_ont, telo_nc_ont = load_ont_telomere_report(args.ont_telo_combined_report)

    cov_summary_s = cov_summary_d = cov_summary_ont = None
    if args.coverage_summary or args.single_coverage_summary or args.dual_coverage_summary or args.ont_coverage_summary:
        print('Loading coverage summary...')
        cov_summary_s, cov_summary_d, cov_summary_ont = load_coverage_summary(
            path_s=args.single_coverage_summary,
            path_d=args.dual_coverage_summary,
            path_ont=args.ont_coverage_summary,
            combined=args.coverage_summary
        )

    centromeres = None
    if args.centromeres:
        print('Loading centromere annotations...')
        centromeres = load_centromeres(args.centromeres)

    flip_set = None
    if args.telo_p_bed:
        print('Building orientation flip set...')
        p_arm_pos = load_p_arm_bed(args.telo_p_bed)
        flip_set  = build_flip_set(p_arm_pos, seq_sizes)
        if flip_set:
            print(f'  {len(flip_set)} chromosome(s) flipped: {sorted(flip_set)}')
        else:
            print('  No flips needed.')

    if centromeres:
        print('Building centromere-based orientation...')
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
            print(f'  Centromere flip: {len(cen_flip)} chromosome(s)')
            flip_set = (flip_set or set()) | cen_flip

    dpi = args.dpi or STYLES.get(args.style, STYLES['paper'])['dpi_default']
    print(f'Building combined figure (style={args.style}, dpi={dpi})...')
    build_combined_figure(df_s, df_d,
                          cov_s, cov_d, nc_cov_s, nc_cov_d,
                          switch_s, switch_d,
                          telo_s, telo_d, telo_nc_s, telo_nc_d,
                          cov_summary_s, cov_summary_d, cov_summary_ont,
                          cov_ont, nc_cov_ont,
                          telo_ont, telo_nc_ont,
                          df_stats=df_stats,
                          out_path=out, dpi=dpi, fmt=args.format,
                          style_name=args.style,
                          rasterized=args.rasterize,
                          simplify=args.simplify,
                          centromeres=centromeres,
                          flip_set=flip_set,
                          annotation_df=annotation_df)
    print('Done.')


if __name__ == '__main__':
    main()
