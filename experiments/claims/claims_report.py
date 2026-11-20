"""Table and figures of the claims test (reads the CSVs written by claims_test.ipynb / claims_lib)."""
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle, Polygon, Ellipse
from matplotlib.lines import Line2D

from claims_lib import *

# ------------------------------------------------------------------ style (dataviz skill: thin marks, recessive chrome)
INK, INK2, MUTED, GRID, SURFACE = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#fcfcfb'
SEQ_BLUE = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b']   # palette.md steps 100-700
CM_SEQ = LinearSegmentedColormap.from_list('seq_blue', SEQ_BLUE)
SEQ_RED = ['#fbdcda', '#f5b4b1', '#ee8a86', '#e34948', '#c93a39', '#a62d2c', '#7f1f1e']   # one-hue red ramp for b-positive cells
CM_SEQ_RED = LinearSegmentedColormap.from_list('seq_red', SEQ_RED)
CM_DIV = LinearSegmentedColormap.from_list('div_rb', ['#c0392b', '#e34948', '#f0efec', '#2a78d6', '#184f95'])  # red-gray-blue
EST_COL = {'classical': C_GR, 'positive': C_BPOS}   # the paper's estimator colours (Figs 4-10)
STYLE_MK = {'thrust': ('^', '#e34948'), 'normal': ('v', '#2a78d6'), 'strike-slip': ('s', '#898781'), 'oblique': ('o', '#c3c2b7')}
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 8, 'axes.edgecolor': '#c3c2b7', 'axes.linewidth': 0.6,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.labelcolor': INK2, 'text.color': INK,
                     'axes.titlesize': 9, 'axes.titleweight': 'semibold', 'legend.frameon': False,
                     'pdf.fonttype': 42, 'savefig.facecolor': 'white'})

ROWS = [  # spatial rows (key, label, claim group, base kind)
    # 2026-10-29: the two strain-rate rows were merged into one (low-b & high rate + high-b & low rate over all nodes with
    # strain data, chance = 50% because the strain split is at the node median); the two GCMT mechanism rows were removed
    # (only LMS thrust, WYN normal and NW reached 10 mechanisms). The old keys stay in build_table for events_*.csv readers.
    ('C12_rate', 'b and strain rate on opposite sides\n(low-b & high rate, high-b & low rate)', 'C1-C2', 'half'),
    ('C1_regime', 'Contractional nodes that are low-b', 'C1', 'low'),
    ('C2_regime', 'Extensional nodes that are high-b', 'C2', 'high'),
    ('C3_hf', 'High-heat-flow nodes that are high-b', 'C3', 'high'),
]
C4_COLS = [  # temporal columns (key, label)
    ('drop', 'Marked events preceded by a b drop'),
    ('rise', 'Marked events followed by a b rise'),
    # ('hit', 'Drops followed by an M>=4 within 7 d'),  removed 2026-10-29: a false-alarm count needs an a priori target
    # and validity window that regional windows cannot supply; alarm episodes are still written to alarms_C4.csv
]


def _pct(k, n):
    return 100.0 * k / n if n else np.nan


def utsu_p(b1, n1, b2, n2):
    """Utsu (1992) probability that two samples (b1, n1) and (b2, n2) share one b value, as written in
    Wiemer & Wyss (1997, JGR 102, eq. after Utsu 1992): P = exp(-dA/2 - 2),
    dA = -2 n ln n + 2 n1 ln(n1 + n2 b1/b2) + 2 n2 ln(n1 b2/b1 + n2) - 2, n = n1 + n2."""
    b1, n1, b2, n2 = (np.asarray(x, float) for x in (b1, n1, b2, n2))
    n = n1 + n2
    dA = -2 * n * np.log(n) + 2 * n1 * np.log(n1 + n2 * b1 / b2) + 2 * n2 * np.log(n1 * b2 / b1 + n2) - 2
    return np.exp(-dA / 2 - 2)


P_SIG = 0.05   # Utsu-test level of the sensitivity rows (*_utsu in claims_table.csv): low-b / high-b only if the test rejects
               # b_node = b0 (Wiemer & Wyss 1997; Schorlemmer et al. 2004, 2005). Main rows: plain below/above-b0 split
               # (Schorlemmer et al. 2005, classes relative to the regional average b); Yuxin's choice 2026-10-30 after the
               # Utsu version showed the same direction in every cell.


def node_classes(v, bcol, ncol, b0, n0, sig=True):
    """Boolean (low, high) node masks: below/above the regional b0; with sig=True also Utsu p < P_SIG."""
    low, high = v[bcol] < b0, v[bcol] > b0
    if sig:
        p = utsu_p(v[bcol].to_numpy(), v[ncol].to_numpy(), b0, n0)
        low, high = low & (p < P_SIG), high & (p < P_SIG)
    return low, high


def build_table(thr, gcmt, min_cov_nodes=20, min_mech=10):
    """Long table of the spatial rows: one row per (claim row, region, estimator) with k, n, pct, base.

    Rows of ROWS use the plain below/above-b0 split; the same rows restricted to nodes whose b differs from b0 at
    Utsu p < P_SIG are appended with the suffix _utsu (sensitivity, CSV only)."""
    out = []
    extra = [('C1_rate', 'Low-b nodes with high strain rate', 'C1'), ('C2_rate', 'High-b nodes with low strain rate', 'C2'),
             ('C1_mech', 'Thrust mechanisms in low-b nodes', 'C1'), ('C2_mech', 'Normal mechanisms in high-b nodes', 'C2')]
    for region in REGIONS:
        nd = pd.read_csv(OUT / f'nodes_{region}.csv')
        g = gcmt[in_box(gcmt.lon, gcmt.lat, BOXES[region]) & (gcmt.mw >= 5.0)]
        for est, bcol in (('classical', 'b_c'), ('positive', 'b_p')):
            b0 = thr.loc[region, bcol]
            ncol = 'n_c' if est == 'classical' else 'n_p'
            n0 = thr.loc[region, ncol] if ncol in thr.columns else thr.loc[region, 'n'] * (1.0 if est == 'classical' else 0.5)
            v = nd[np.isfinite(nd[bcol])].copy()
            for sig, suffix in ((False, ''), (True, '_utsu')):
                low, high = node_classes(v, bcol, ncol, b0, n0, sig=sig)
                v['low'], v['high'] = low.to_numpy(), high.to_numpy()
                base_low, base_high = 100 * low.mean(), 100 * high.mean()
                share = (base_low, base_high, 100 * (~low & ~high).mean())
                rec = {}
                s_ = v[np.isfinite(v.ws_e2inv)]
                sl, sh = s_.low, s_.high
                hs = s_.ws_e2inv > s_.ws_e2inv.median()
                rec['C1_rate'] = ((sl & hs).sum(), sl.sum(), 50.0)
                rec['C2_rate'] = ((sh & ~hs).sum(), sh.sum(), 50.0)
                rec['C12_rate'] = ((sl & hs).sum() + (sh & ~hs).sum(), sl.sum() + sh.sum(), 50.0)
                q1, q2 = s_.ws_dilat.quantile([1 / 3, 2 / 3])
                contr, ext = s_.ws_dilat <= q1, s_.ws_dilat >= q2
                rec['C1_regime'] = ((contr & sl).sum(), contr.sum(), base_low)
                rec['C2_regime'] = ((ext & sh).sum(), ext.sum(), base_high)
                nl = nearest_node_value(v, g.lon.to_numpy(), g.lat.to_numpy(), col='low') if len(g) else np.array([])
                nh = nearest_node_value(v, g.lon.to_numpy(), g.lat.to_numpy(), col='high') if len(g) else np.array([])
                th = (g['style'] == 'thrust').to_numpy() & np.isfinite(nl)
                no = (g['style'] == 'normal').to_numpy() & np.isfinite(nh)
                rec['C1_mech'] = ((nl[th] > 0.5).sum(), th.sum(), base_low) if th.sum() >= min_mech else (np.nan, th.sum(), base_low)
                rec['C2_mech'] = ((nh[no] > 0.5).sum(), no.sum(), base_high) if no.sum() >= min_mech else (np.nan, no.sum(), base_high)
                vh = v[np.isfinite(v.hf)]
                hh = vh.hf > vh.hf.median()
                rec['C3_hf'] = ((vh.high & hh).sum(), hh.sum(), base_high) if (len(vh) >= min_cov_nodes and region != 'XSH') else (np.nan, 0, base_high)
                for key, label, grp in [(k_, l_, g_) for k_, l_, g_, _ in ROWS] + extra:   # extra rows go to the CSV only
                    k, n, base = rec[key]
                    out.append(dict(row=key + suffix, label=label, claim=grp, region=region, estimator=est,
                                    k=k, n=n, pct=_pct(k, n) if np.isfinite(k) else np.nan, base=base,
                                    share_low=share[0], share_high=share[1], share_neutral=share[2]))
    return pd.DataFrame(out)


def build_c4_table(summaries):
    """Long table of the temporal columns: one row per (column, window, estimator) with k, n, pct; plus totals."""
    sm = pd.DataFrame(summaries); out = []
    for est in ('classical', 'positive'):
        u = sm[sm.estimator == est]
        for key, label in C4_COLS:
            kk, nn = {'drop': ('drop_k', 'drop_n'), 'rise': ('rise_k', 'rise_n'), 'hit': ('n_hits', 'n_episodes')}[key]
            for _, r in u.iterrows():
                out.append(dict(col=key, label=label, window=r.window, fig=r.fig, estimator=est, k=r[kk], n=r[nn],
                                pct=_pct(r[kk], r[nn]) if r[nn] else np.nan, alarm_fraction=r.alarm_fraction,
                                hit_ref_share=r.hit_ref_share if key == 'hit' else np.nan,
                                n_episodes=r.n_episodes, n_marked=r.n_marked))
            out.append(dict(col=key, label=label, window='ALL', fig='Total', estimator=est, k=u[kk].sum(), n=u[nn].sum(),
                            pct=_pct(u[kk].sum(), u[nn].sum()) if u[nn].sum() else np.nan, alarm_fraction=np.nan, hit_ref_share=np.nan,
                            n_episodes=u.n_episodes.sum(), n_marked=u.n_marked.sum()))
    return pd.DataFrame(out)


def table_markdown(tab, c4):
    """Two wide markdown tables: spatial rows x regions, and temporal columns x windows."""
    lines = ['**Spatial claims** (cell = k/n (%), classical | positive; ref. = share expected with no relation)', '',
             '| Claim | Row | ' + ' | '.join(REGIONS) + ' |', '|' + '---|' * (len(REGIONS) + 2)]
    for key, label, grp, basek in ROWS:
        label = label.replace('\n', ' '); cells = []
        for region in REGIONS:
            parts = []
            for est in ('classical', 'positive'):
                r = tab[(tab.row == key) & (tab.region == region) & (tab.estimator == est)].iloc[0]
                parts.append(f'{int(r.k)}/{int(r.n)} ({r.pct:.0f}%)' if np.isfinite(r.pct)
                             else (f'n/a ({int(r.n)})' if key.endswith('_mech') else 'n/a'))   # mech: n/a shows the mechanism count
            b2 = tab[(tab.row == key) & (tab.region == region)].base.tolist()
            ref = 'ref. 50%' if basek == 'half' else f'ref. {b2[0]:.0f}/{b2[1]:.0f}%'
            cells.append(' \\| '.join(parts) + ('' if all(p.startswith('n/a') for p in parts) else f'<br>{ref}'))
        lines.append(f'| {grp} | {label} | ' + ' | '.join(cells) + ' |')
    lines += ['', '**Temporal claims (C4), main-text windows** (cell = k/n (%), classical | positive)', '',
              '| Window | Marked events | ' + ' | '.join(l for _, l in C4_COLS) + ' |',
              '|' + '---|' * (len(C4_COLS) + 2)]
    for fig in [FIG_OF[l] for l in MAIN_LABELS]:
        cells = []
        for key, _ in C4_COLS:
            parts = []
            for est in ('classical', 'positive'):
                r = c4[(c4.col == key) & (c4.fig == fig) & (c4.estimator == est)].iloc[0]
                parts.append(f'{int(r.k)}/{int(r.n)} ({r.pct:.0f}%)' if np.isfinite(r.pct) else 'n/a')
            cells.append(' \\| '.join(parts))
        r0 = c4[(c4.col == 'drop') & (c4.fig == fig)]
        lines.append(f'| {fig} {r0.window.iloc[0]} | {int(r0.n_marked.iloc[0])} | ' + ' | '.join(cells) + ' |')
    sh = []
    for region in REGIONS:
        r2 = [tab[(tab.row == 'C12_rate') & (tab.region == region) & (tab.estimator == e)].iloc[0] for e in ('classical', 'positive')]
        sh.append(f'{region} ' + ' \\| '.join(f'{r.share_low:.0f}/{r.share_high:.0f}' for r in r2))
    lines.append('')
    lines.append('Node classes low / high (% of valid nodes, classical | positive): ' + '; '.join(sh) + '.')
    lines.append('')
    lines.append('Drop = last window before the event >= 10% below the pre-sequence reference (median of the earlier windows of the panel); '
                 'rise = first post-event window >= 10% above the same reference; n = marked events with a valid window. '
                 'Low-b / high-b node = below / above the regional whole-period b (b0; Schorlemmer et al. 2005, classes relative '
                 'to the regional average). Restricting the classes to nodes that differ from b0 at Utsu (1992) p < 0.05 '
                 '(Wiemer & Wyss 1997; Schorlemmer et al. 2004) leaves the direction of every cell unchanged (rows *_utsu in '
                 'claims_table.csv). Strain-rate row: n = all nodes with strain data; chance = 50% because the strain split is '
                 'at the node median.')
    return '\n'.join(lines)


# ------------------------------------------------------------------ C4 per-event table (Form B, chosen 2026-10-29)
SYM_YES, SYM_NO, SYM_NA = '●', '○', '·'   # filled circle, open circle, middle dot


def _sym(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return SYM_NA
    return SYM_YES if bool(v) else SYM_NO


def c4_event_table(events):
    """Form B of the C4 result: one row per main-text window, one cell per marked event in time order.

    Cell = magnitude, '*' if the event is an aftershock of an earlier marked event (not a sequence initiator),
    then two symbols: left = classical, right = b-positive; filled = criterion met, open = not met, dot = no valid window.
    Returns (long DataFrame, markdown string). The markdown table shows the drop symbols and the drop k/n of both estimators;
    the rise results are summarised in the footnote because almost no event meets the rise criterion.
    """
    ev = events.copy()
    ev['time'] = pd.to_datetime(ev['time'])
    rows, n_cols = [], 0
    for lab in MAIN_LABELS:
        g = ev[ev.window == lab]
        if g.empty:
            continue
        keys = sorted(set(zip(g.time, g.mag)))
        n_cols = max(n_cols, len(keys))
        for i, (t, m) in enumerate(keys, 1):
            gg = g[(g.time == t) & (g.mag == m)]
            r = {x.estimator: x for x in gg.itertuples()}
            c, p = r.get('classical'), r.get('positive')
            a = c if c is not None else p
            rows.append(dict(window=lab, fig=FIG_OF[lab], idx=i, time=t, mag=m, initiator=bool(a.initiator),
                             drop_c=getattr(c, 'drop', np.nan), drop_p=getattr(p, 'drop', np.nan),
                             rise_c=getattr(c, 'rise', np.nan), rise_p=getattr(p, 'rise', np.nan),
                             drop_rel_c=getattr(c, 'drop_rel', np.nan), drop_rel_p=getattr(p, 'drop_rel', np.nan),
                             rise_rel_c=getattr(c, 'rise_rel', np.nan), rise_rel_p=getattr(p, 'rise_rel', np.nan)))
    long = pd.DataFrame(rows)

    def kn(col, g):
        v = g[col]
        ok = v.notna()
        return f'{int(v[ok].astype(bool).sum())}/{int(ok.sum())}'

    lines = ['**C4 per marked event, main-text windows** (cell = magnitude, * = aftershock of an earlier marked event, '
             'then classical | positive: ' + SYM_YES + ' drop before, ' + SYM_NO + ' no drop, ' + SYM_NA + ' no valid window)', '',
             '| Window | ' + ' | '.join(f'{i}' for i in range(1, n_cols + 1)) + ' | Drop, c \\| p |',
             '|' + '---|' * (n_cols + 2)]
    rise_note = []
    for lab in MAIN_LABELS:
        g = long[long.window == lab]
        if g.empty:
            continue
        cells = [f'{r.mag:.1f}{"" if r.initiator else "*"} {_sym(r.drop_c)}{_sym(r.drop_p)}' for r in g.itertuples()]
        cells += [''] * (n_cols - len(cells))
        lines.append(f'| {FIG_OF[lab]} {lab} | ' + ' | '.join(cells) + f' | {kn("drop_c", g)} \\| {kn("drop_p", g)} |')
        rise_note.append(f'{FIG_OF[lab]} {kn("rise_c", g)} | {kn("rise_p", g)}')
    lines += ['', 'Rise after (first post-event window >= 10% above the same reference), classical | positive: ' + '; '.join(rise_note) + '.',
              'Drop = last window ending before the event >= 10% below the pre-sequence reference (median of the earlier '
              'windows of the panel outside earlier sequences); events in time order within each window.']
    return long, '\n'.join(lines)


# ------------------------------------------------------------------------------- figures
def _cell(ax, x0, y0, w, h, pct, est, k, n, invert=False):
    norm = Normalize(0, 100)
    if np.isfinite(pct):
        shade = 100 - pct if invert else pct
        col = (CM_SEQ if est == 'classical' else CM_SEQ_RED)(norm(shade)); txtc = 'white' if shade > 55 else INK
        ax.add_patch(Rectangle((x0 + 0.02, y0 + 0.04), w - 0.04, h - 0.08, fc=col, ec='none'))
        ax.text(x0 + w / 2, y0 + 0.40, f'{pct:.0f}%', ha='center', va='center', fontsize=7.5, color=txtc, fontweight='semibold')
        kn = f'{int(k)}/{int(n)}'
        ax.text(x0 + w / 2, y0 + 0.75, kn, ha='center', va='center', fontsize=5.6 if len(kn) <= 7 else 4.6, color=txtc)
    else:
        ax.add_patch(Rectangle((x0 + 0.02, y0 + 0.04), w - 0.04, h - 0.08, fc='#f3f2ee', ec='none'))
        ax.text(x0 + w / 2, y0 + 0.55, 'n/a', ha='center', va='center', fontsize=6.5, color=MUTED)


def fig_summary(tab, c4, path, events=None):
    """(a) spatial claims x regions, cells split classical | positive, colour = % supporting, text = k/n;
    (b) temporal claims: one glyph per marked event of each main-text window (left half classical, right half
    positive; filled = criterion met), size = magnitude, dashed outline = aftershock of an earlier marked event,
    and the k/n of both estimators at the end of each strip."""
    nrow, ncol = len(ROWS), len(REGIONS)
    wins = [FIG_OF[l] for l in MAIN_LABELS]
    nrow2, ncol2 = len(wins), len(C4_COLS)
    fig = plt.figure(figsize=(7.2, 0.42 * (nrow + nrow2) + 3.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[nrow + 1.2, nrow2 + 1.9], hspace=0.42, left=0.40, right=0.955, top=0.95, bottom=0.12)
    REF = {'half': 'ref. 50%', 'low': 'ref. share of low-b nodes', 'high': 'ref. share of high-b nodes'}
    # (a)
    ax = fig.add_subplot(gs[0]); ax.set_xlim(0, ncol); ax.set_ylim(nrow, 0); ax.axis('off')
    for i, (key, label, grp, basek) in enumerate(ROWS):
        for j, region in enumerate(REGIONS):
            for h, est in enumerate(('classical', 'positive')):
                r = tab[(tab.row == key) & (tab.region == region) & (tab.estimator == est)].iloc[0]
                _cell(ax, j + 0.5 * h, i, 0.5, 1, r.pct, est, r.k, r.n)
        ax.text(-0.08, i + 0.36, label, ha='right', va='center', fontsize=7, color=INK, linespacing=0.95)
        ax.text(-0.08, i + 0.76, REF[basek], ha='right', va='center', fontsize=5.5, color=MUTED)
        if i > 0 and grp.startswith('C3') and not ROWS[i - 1][2].startswith('C3'):
            ax.plot([0, ncol], [i, i], color='white', lw=3)   # gap between the stress rows and the fluid rows
    for j, region in enumerate(REGIONS):
        name = {'LMS': 'Longmen-\nshan', 'XSH': 'Xianshui-\nhe', 'WYN': 'Western\nYunnan', 'XJ': 'Xiaojiang', 'NCP': 'North\nChina Plain', 'NW': 'Northwest\nChina'}[region]
        ax.text(j + 0.5, -0.45, name, ha='center', va='center', fontsize=6.8, fontweight='semibold', color=INK, linespacing=0.95)
        ax.text(j + 0.25, -0.07, 'classical', ha='center', va='center', fontsize=5.3, color=EST_COL['classical'])
        ax.text(j + 0.75, -0.07, 'positive', ha='center', va='center', fontsize=5.3, color=EST_COL['positive'])
    ax.text(-0.08, -0.6, 'a  Spatial claims, 2008\u20132023 maps', ha='right', va='center', fontsize=8, fontweight='bold', color=INK)
    n_stress = sum(not g.startswith('C3') for _, _, g, _ in ROWS)
    for grp, y in (('Stress (C1, C2)', n_stress / 2), ('Fluids (C3)', n_stress + (nrow - n_stress) / 2)):
        ax.text(ncol + 0.08, y, grp, rotation=270, ha='left', va='center', fontsize=7, color=INK2)
    # (b) one glyph per marked event
    ax = fig.add_subplot(gs[1]); ax.set_xlim(0, 2); ax.set_ylim(nrow2, 0); ax.axis('off')
    bb = ax.get_position()
    xpp = 2 / (bb.width * fig.get_figwidth() * 72); ypp = nrow2 / (bb.height * fig.get_figheight() * 72)   # data units per point
    ev = events.copy(); ev['time'] = pd.to_datetime(ev['time'])
    n_max = int(ev.groupby(['window', 'estimator']).size().max()); dx = 0.66 / n_max

    def tri(v):
        if isinstance(v, str):
            return {'True': True, 'False': False}.get(v)
        return None if v is None or (isinstance(v, float) and np.isnan(v)) else bool(v)

    def glyph(cx, cy, m, yes_c, yes_p, dashed):
        r = 2.6 + 2.9 * (min(max(m, 4.4), 7.9) - 4.4) / 3.5          # radius in points
        for half, yes in (('L', yes_c), ('R', yes_p)):
            th = np.linspace(np.pi / 2, 3 * np.pi / 2, 30) if half == 'L' else np.linspace(-np.pi / 2, np.pi / 2, 30)
            col = EST_COL['classical' if half == 'L' else 'positive']
            fc = col if yes is True else ('#e1e0d9' if yes is None else 'white')
            ax.add_patch(Polygon(np.column_stack([cx + r * xpp * np.cos(th), cy + r * ypp * np.sin(th)]), closed=True, fc=fc, ec='none', zorder=3, clip_on=False))
        ax.add_patch(Ellipse((cx, cy), 2 * r * xpp, 2 * r * ypp, fc='none', ec=INK2, lw=0.55, ls=(0, (1.5, 1.2)) if dashed else '-', zorder=4, clip_on=False))
        return r
    for i, (label_w, figname) in enumerate(zip(MAIN_LABELS, wins)):
        ec = ev[(ev.window == label_w) & (ev.estimator == 'classical')].sort_values('time').reset_index(drop=True)
        ep = ev[(ev.window == label_w) & (ev.estimator == 'positive')].sort_values('time').reset_index(drop=True)
        cy = i + 0.58
        for g, key in enumerate(('drop', 'rise')):
            for k in range(len(ec)):
                cx = g + 0.03 + dx * (k + 0.5)
                r = glyph(cx, cy, ec.mag[k], tri(ec[key][k]), tri(ep[key][k]), not bool(ec.initiator[k]))
                if ec.mag[k] >= 6:
                    ax.text(cx, cy - (r + 1.2) * ypp, f'{ec.mag[k]:.1f}', ha='center', va='bottom', fontsize=4.8, color=INK2)
            for h, est in enumerate(('classical', 'positive')):
                rr = c4[(c4.col == key) & (c4.fig == figname) & (c4.estimator == est)].iloc[0]
                ax.text(g + 0.785 + 0.175 * h, cy, f'{int(rr.k)}/{int(rr.n)}' if np.isfinite(rr.pct) else 'n/a',
                        ha='center', va='center', fontsize=6.5, color=EST_COL[est], fontweight='semibold')
        lab = f"{figname}  {REGION_NAME[label_w.split('_')[0]]} {label_w.split('_')[1]}"
        ax.text(-0.03, i + 0.42, lab, ha='right', va='center', fontsize=7, color=INK)
        ax.text(-0.03, i + 0.74, f"{len(ec)} events · {int((ec.mag >= 6).sum())} with M≥6", ha='right', va='center', fontsize=5.2, color=MUTED)
        if i:
            ax.plot([0, 2], [i, i], color=GRID, lw=0.4, zorder=1)
    for g, short in enumerate(('Marked events preceded by a b drop', 'Marked events followed by a b rise')):
        ax.text(g + 0.5, -0.55, short, ha='center', va='center', fontsize=6.6, fontweight='semibold', color=INK)
        ax.text(g + 0.05, -0.12, 'time →', ha='left', va='center', fontsize=5, color=MUTED)
        ax.text(g + 0.785, -0.12, 'classical', ha='center', va='center', fontsize=5.1, color=EST_COL['classical'])
        ax.text(g + 0.96, -0.12, 'positive', ha='center', va='center', fontsize=5.1, color=EST_COL['positive'])
    from matplotlib.transforms import blended_transform_factory
    ax.text(0.03, -0.9, 'b  Temporal claims (C4), main-text windows', transform=blended_transform_factory(fig.transFigure, ax.transData),
            ha='left', va='center', fontsize=8, fontweight='bold', color=INK)
    # legend glyphs
    items = [((True, True, False), 'met by both estimators'), ((True, False, False), 'b-classical only'), ((False, True, False), 'b-positive only'),
             ((False, False, False), 'not met'), ((None, None, False), 'no valid window'), ((False, False, True), 'aftershock of an earlier marked event')]
    xs = [-1.38, -0.86, -0.46, -0.08, 0.20, 0.60]          # the legend runs under the row labels too (axes x < 0 is free)
    for (yc, yp, dsh), lab, x in zip([i[0] for i in items], [i[1] for i in items], xs):
        glyph(x, nrow2 + 0.42, 5.5, yc, yp, dsh)
        ax.text(x + 0.035, nrow2 + 0.42, lab, ha='left', va='center', fontsize=5.4, color=INK2)
    ax.text(-1.40, nrow2 + 0.78, 'Circle size = magnitude (M≥6 labelled); left half b-classical, right half b-positive; filled = ≥ 10% below (drop) / '
            'above (rise) the pre-sequence reference.', ha='left', va='center', fontsize=5.4, color=MUTED)
    norm = Normalize(0, 100)
    for x, cm, lab in ((0.42, CM_SEQ, 'b-classical'), (0.70, CM_SEQ_RED, 'b-positive')):
        cax = fig.add_axes([x, 0.045, 0.22, 0.012])
        cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cm, norm=norm, orientation='horizontal')
        cb.set_label(f'{lab}: % supporting (small text: k/n)', fontsize=6.3, color=INK2); cb.ax.tick_params(labelsize=6, colors=MUTED); cb.outline.set_visible(False)
    fig.savefig(path); plt.close(fig)


def _faults(ax, box, lw=0.4):
    """3D-CFM fault points (Sichuan-Yunnan only) as a faint layer."""
    if not FAULTS_CFM.exists():
        return
    f = pd.read_csv(FAULTS_CFM, sep=r'\s+', header=None, usecols=[2, 3], names=['lon', 'lat'])
    f = f[in_box(f.lon, f.lat, box)]
    if len(f):
        ax.scatter(f.lon, f.lat, s=0.15, c='#52514e', alpha=0.35, lw=0, rasterized=True)


def _map_axes(ax, box):
    x0, x1, y0, y1 = box
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect(1 / np.cos(np.radians((y0 + y1) / 2)))
    ax.tick_params(labelsize=6, length=2)
    for s in ax.spines.values():
        s.set_linewidth(0.5)


def _grid(nd, col):
    lons, lats = np.sort(nd.lon.unique()), np.sort(nd.lat.unique())
    Z = nd.pivot(index='lat', columns='lon', values=col).reindex(index=lats, columns=lons).to_numpy()
    return lons, lats, Z


def figS_C12(thr, gcmt, path, bcol='b_c'):
    """Per region: b anomaly (b - b0), strain-rate magnitude, dilatation; GCMT mechanisms on the b map."""
    fig, axes = plt.subplots(len(REGIONS), 3, figsize=(7.2, 2.35 * len(REGIONS)))
    est = 'b_c' if bcol == 'b_c' else 'b_p'
    for i, region in enumerate(REGIONS):
        nd = pd.read_csv(OUT / f'nodes_{region}.csv'); box = BOXES[region]
        b0 = thr.loc[region, est]
        lons, lats, Zb = _grid(nd, bcol)
        _, _, Ze = _grid(nd, 'ws_e2inv'); _, _, Zd = _grid(nd, 'ws_dilat')
        # b anomaly
        ax = axes[i, 0]; _map_axes(ax, box)
        dv = np.nanpercentile(np.abs(Zb - b0), 95) if np.isfinite(Zb).any() else 0.3
        im = ax.pcolormesh(lons, lats, Zb - b0, cmap=CM_DIV, norm=TwoSlopeNorm(0, -dv, dv), shading='nearest', rasterized=True)
        _faults(ax, box)
        g = gcmt[in_box(gcmt.lon, gcmt.lat, box) & (gcmt.mw >= 5.0)]
        for sty, (mk, c) in STYLE_MK.items():
            gg = g[g['style'] == sty]
            if len(gg):
                ax.scatter(gg.lon, gg.lat, s=9, marker=mk, fc=c, ec='white', lw=0.3, zorder=5)
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.ax.tick_params(labelsize=5); cb.outline.set_visible(False)
        ax.set_title(f'{REGION_NAME[region]}\nb − b₀ (b₀ = {b0:.2f})', loc='left', fontsize=7.5)
        # strain rate magnitude
        ax = axes[i, 1]; _map_axes(ax, box)
        vmax = np.nanpercentile(Ze, 98) if np.isfinite(Ze).any() else 1
        im = ax.pcolormesh(lons, lats, Ze, cmap=CM_SEQ, vmin=0, vmax=vmax, shading='nearest', rasterized=True)
        _faults(ax, box)
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.ax.tick_params(labelsize=5); cb.outline.set_visible(False)
        ax.set_title('\nStrain rate, 2nd invariant\n(nstrain/yr)', loc='left', fontsize=7.5)
        # dilatation
        ax = axes[i, 2]; _map_axes(ax, box)
        dmax = np.nanpercentile(np.abs(Zd), 95) if np.isfinite(Zd).any() else 1
        im = ax.pcolormesh(lons, lats, Zd, cmap=CM_DIV, norm=TwoSlopeNorm(0, -dmax, dmax), shading='nearest', rasterized=True)
        _faults(ax, box)
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.ax.tick_params(labelsize=5); cb.outline.set_visible(False)
        ax.set_title('\nDilatation rate (nstrain/yr)\nred = contraction', loc='left', fontsize=7.5)
    handles = [Line2D([], [], marker=mk, ls='', mfc=c, mec='white', ms=5, label=sty) for sty, (mk, c) in STYLE_MK.items()]
    fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=7, title='GCMT M≥5 mechanisms (1976–2023), on the b maps', title_fontsize=7, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.02, 1, 1), w_pad=1.5)
    fig.savefig(path); plt.close(fig)


def figS_C3(thr, hf, path, bcol='b_c'):
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.2))
    est = 'b_c' if bcol == 'b_c' else 'b_p'
    hnorm = Normalize(30, 100)
    for ax, region in zip(axes.ravel(), REGIONS):
        nd = pd.read_csv(OUT / f'nodes_{region}.csv'); box = BOXES[region]; b0 = thr.loc[region, est]
        lons, lats, Zb = _grid(nd, bcol); _map_axes(ax, box)
        dv = np.nanpercentile(np.abs(Zb - b0), 95) if np.isfinite(Zb).any() else 0.3
        ax.pcolormesh(lons, lats, Zb - b0, cmap=CM_DIV, norm=TwoSlopeNorm(0, -dv, dv), shading='nearest', rasterized=True, alpha=0.9)
        _faults(ax, box)
        h = hf[in_box(hf.lon, hf.lat, box)]
        if len(h):
            ax.scatter(h.lon, h.lat, s=14, c=h.hf, cmap='Oranges', norm=hnorm, ec=INK, lw=0.3, zorder=5)
        ax.set_title(f'{REGION_NAME[region]}\n{len(h)} heat-flow sites', loc='left', fontsize=7.5)
    handles = [Line2D([], [], marker='o', ls='', mfc='#eda100', mec=INK, ms=5, label='heat-flow site (colour = mW/m²)')]
    fig.legend(handles=handles, loc='lower center', ncol=1, fontsize=7, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle('b − b₀ (red low, blue high) with heat-flow sites', fontsize=9, x=0.02, ha='left')
    fig.tight_layout(rect=(0, 0.03, 1, 0.97), w_pad=1.5)
    fig.savefig(path); plt.close(fig)


def figS_C4(events, figrun, path, days=365):
    """Small multiples: b(t) of both estimators around each marked event of the main-text windows, with the
    pre-event reference (dashed), the pre-event and post-event windows (thick) and the event (vertical line)."""
    ev = events[events.estimator == 'classical'].copy()
    ev['order'] = ev['window'].map({l: i for i, l in enumerate(MAIN_LABELS)})
    ev = ev.sort_values(['order', 'time']).reset_index(drop=True)
    n = len(ev); ncol = 4; nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.2, 1.55 * nrow))
    for ax, (_, r) in zip(axes.ravel(), ev.iterrows()):
        t_ev = pd.Timestamp(r.time); lo, hi = t_ev - pd.Timedelta(days=days), t_ev + pd.Timedelta(days=days)
        for est, cls in (('classical', ClassicBValueEstimator), ('positive', BPositiveBValueEstimator)):
            s_, _ = series_for_window(r.window, cls, figrun)
            w = s_[(s_.t_last >= lo) & (s_.t_last <= hi)]
            x = (w.t_last - t_ev).dt.total_seconds() / 86400
            ax.plot(x, w.b, color=EST_COL[est], lw=0.8)
            e2 = events[(events.estimator == est) & (events.window == r.window) & (pd.to_datetime(events.time) == t_ev)].iloc[0]
            if np.isfinite(e2.ref):
                ax.plot([-e2.pre_days - 30, 0], [e2.ref, e2.ref], color=EST_COL[est], lw=0.6, ls=(0, (3, 2)), alpha=0.8)
            if np.isfinite(e2.b_pre):
                ax.plot([-e2.pre_days, 0], [e2.b_pre, e2.b_pre], color=EST_COL[est], lw=2.2, alpha=0.9)
            if np.isfinite(e2.b_post):
                ax.plot([e2.post_start_days, e2.post_start_days + e2.post_days], [e2.b_post, e2.b_post], color=EST_COL[est], lw=2.2, alpha=0.9)
        ax.axvline(0, color=INK, lw=0.6)
        ax.set_xlim(-days, days)
        ax.set_title(f'{r.fig} #{int(r.number)} {t_ev:%Y-%m-%d} M{r.mag:.1f}', fontsize=6.3, loc='left')
        ax.tick_params(labelsize=5.5); ax.grid(axis='y', color=GRID, lw=0.4)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    for ax in axes.ravel()[n:]:
        ax.axis('off')
    handles = [Line2D([], [], color=C_GR, lw=1.2, label='b-classical'), Line2D([], [], color=C_BPOS, lw=1.2, label='b-positive'),
               Line2D([], [], color=INK2, lw=0.8, ls=(0, (3, 2)), label='reference (median of the earlier pre-sequence windows of the panel)'),
               Line2D([], [], color=INK2, lw=2.2, label='pre-event / post-event window')]
    fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=6.5, bbox_to_anchor=(0.5, 0.0))
    fig.text(0.5, 0.025, 'days relative to the event', ha='center', fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path); plt.close(fig)


# ------------------------------------------------------------------ main-text table (spatial rows + C4 Form B) rendered as a figure

TABLE_ROWS = [  # (key, title, sub-line)
    ('C12_rate', 'C1–C2  Low b with high strain rate,\n           high b with low strain rate', 'share of all nodes; chance 50%'),
    ('C1_regime', 'C1  Contractional nodes that are low-b', 'reference: share of low-b nodes'),
    ('C2_regime', 'C2  Extensional nodes that are high-b', 'reference: share of high-b nodes'),
    ('C3_hf', 'C3  High-heat-flow nodes that are high-b', 'reference: share of high-b nodes'),
]
REGION_NAMES = {'LMS': 'Longmenshan', 'XSH': 'Xianshuihe', 'WYN': 'Western\nYunnan', 'XJ': 'Xiaojiang',
                'NCP': 'North China\nPlain', 'NW': 'Northwest\nChina'}
WINDOW_NAMES = {'LMS_2006-2010': 'Longmenshan 2006–2010', 'LMS_2011-2016': 'Longmenshan 2011–2016',
                'XSH_2022-2023': 'Xianshuihe 2022–2023', 'WYN_2020-2023': 'Western Yunnan 2020–2023',
                'XJ_2016-2023': 'Xiaojiang 2016–2023', 'NCP_2008-2023': 'North China Plain 2008–2023',
                'NW_2006-2023': 'Northwest China 2006–2023'}


def fig_table(tab, events, path):
    """One figure with the two main-text tables: (a) spatial claims per region, (b) C4 per marked event (Form B)."""
    long, _ = c4_event_table(events)
    n_ev = int(long.groupby('window').size().max())
    W = 7.2
    ha, hb = 0.60 * len(TABLE_ROWS) + 0.62, 0.40 * long.window.nunique() + 1.05
    fig = plt.figure(figsize=(W, ha + hb + 0.55))
    gs = fig.add_gridspec(2, 1, height_ratios=[ha, hb], hspace=0.08, left=0.01, right=0.99, top=0.985, bottom=0.015)
    # ---------------------------------------------------------------- (a) spatial claims
    ax = fig.add_subplot(gs[0]); ax.set_xlim(0, W); ax.set_ylim(ha, 0); ax.axis('off')
    x_lab, cw = 2.55, (W - 2.55) / len(REGIONS)
    ax.text(0.02, 0.16, 'a  Spatial claims, 2008–2023 b-value maps (cell: % of nodes supporting the claim, b-classical | b-positive)',
            fontsize=7.2, fontweight='bold', color=INK, va='center')
    y0 = 0.62
    for j, region in enumerate(REGIONS):
        cx = x_lab + (j + 0.5) * cw
        ax.text(cx, y0 - 0.25, REGION_NAMES[region], ha='center', va='center', fontsize=6.6, fontweight='semibold', color=INK, linespacing=0.95)
    ax.plot([0.02, W - 0.02], [y0, y0], color=INK2, lw=0.8)
    ax.plot([0.02, W - 0.02], [0.34, 0.34], color=INK2, lw=0.8)
    rh = 0.60
    for i, (key, title, sub) in enumerate(TABLE_ROWS):
        y = y0 + i * rh
        if i % 2 == 1:
            ax.add_patch(Rectangle((0.02, y), W - 0.04, rh, color='#f3f2ee', lw=0, zorder=0))
        ax.text(0.08, y + 0.22, title, fontsize=6.6, color=INK, va='center', linespacing=1.0)
        ax.text(0.08, y + 0.45, sub, fontsize=5.4, color=MUTED, va='center')
        for j, region in enumerate(REGIONS):
            cx = x_lab + (j + 0.5) * cw
            rr = {est: tab[(tab.row == key) & (tab.region == region) & (tab.estimator == est)].iloc[0] for est in ('classical', 'positive')}
            if not np.isfinite(rr['classical'].pct):
                ax.text(cx, y + 0.30, 'n/a', ha='center', va='center', fontsize=6.2, color=MUTED); continue
            pc, pp = rr['classical'], rr['positive']
            ax.text(cx - 0.05, y + 0.19, f'{pc.pct:.0f}%', ha='right', va='center', fontsize=7.2, color=EST_COL['classical'], fontweight='semibold')
            ax.text(cx, y + 0.19, '|', ha='center', va='center', fontsize=6, color=MUTED)
            ax.text(cx + 0.05, y + 0.19, f'{pp.pct:.0f}%', ha='left', va='center', fontsize=7.2, color=EST_COL['positive'], fontweight='semibold')
            ref = 'ref. 50%' if key == 'C12_rate' else f'ref. {pc.base:.0f}% | {pp.base:.0f}%'
            ax.text(cx, y + 0.35, ref, ha='center', va='center', fontsize=4.9, color=MUTED)
            ax.text(cx, y + 0.47, f'n {int(pc.n):,}' + ('' if pc.n == pp.n else f' | {int(pp.n):,}'), ha='center', va='center', fontsize=4.9, color=MUTED)
    ax.plot([0.02, W - 0.02], [y0 + len(TABLE_ROWS) * rh, y0 + len(TABLE_ROWS) * rh], color=INK2, lw=0.8)
    # ---------------------------------------------------------------- (b) C4 per marked event
    ax = fig.add_subplot(gs[1]); ax.set_xlim(0, W); ax.set_ylim(hb, 0); ax.axis('off')
    x_lab, x_res = 1.95, W - 0.95
    ew = (x_res - x_lab) / n_ev
    ax.text(0.02, 0.16, 'b  Temporal claim C4: b drop before each marked event of Figs 4–10 (left half b-classical, right half b-positive)',
            fontsize=7.2, fontweight='bold', color=INK, va='center')
    y0 = 0.62
    ax.plot([0.02, W - 0.02], [y0, y0], color=INK2, lw=0.8)
    ax.text(x_lab + 0.02, y0 - 0.16, 'events in time order  →', fontsize=5.6, color=MUTED, va='center')
    ax.text(x_res + 0.22, y0 - 0.22, 'drop before', ha='center', va='center', fontsize=6.0, color=INK, fontweight='semibold')
    ax.text(x_res + 0.22, y0 - 0.08, 'classical | positive', ha='center', va='center', fontsize=5.0, color=MUTED)
    rh = 0.40
    bb = ax.get_position(); ypp = hb / (bb.height * fig.get_figheight() * 72); xpp = W / (bb.width * fig.get_figwidth() * 72)
    rpt = 4.2   # glyph radius in points
    for i, lab in enumerate(MAIN_LABELS):
        g = long[long.window == lab]
        y = y0 + i * rh
        if i % 2 == 1:
            ax.add_patch(Rectangle((0.02, y), W - 0.04, rh, color='#f3f2ee', lw=0, zorder=0))
        ax.text(0.08, y + 0.14, f'{FIG_OF[lab]}  {WINDOW_NAMES[lab]}', fontsize=6.4, color=INK, va='center')
        nm6 = int((g.mag >= 6).sum())
        ax.text(0.08, y + 0.29, f'{len(g)} marked events · {nm6} with M ≥ 6', fontsize=5.0, color=MUTED, va='center')
        for k, r in enumerate(g.itertuples()):
            cx, cy = x_lab + (k + 0.5) * ew, y + 0.25
            for half, val, col in ((0, r.drop_c, EST_COL['classical']), (1, r.drop_p, EST_COL['positive'])):
                isna = isinstance(val, float) and np.isnan(val)
                fill = '#e1e0d9' if isna else (col if bool(val) else 'white')
                ang = np.linspace(np.pi / 2, 3 * np.pi / 2, 24) if half == 0 else np.linspace(-np.pi / 2, np.pi / 2, 24)
                pts = np.column_stack([cx + rpt * xpp * np.cos(ang), cy + rpt * ypp * np.sin(ang)])
                ax.add_patch(Polygon(pts, closed=True, facecolor=fill, edgecolor='none', lw=0, clip_on=False, zorder=2))
            ax.add_patch(Ellipse((cx, cy), 2 * rpt * xpp, 2 * rpt * ypp, facecolor='none', edgecolor=INK2, lw=0.6,
                                 ls='-' if r.initiator else (0, (1.5, 1.2)), clip_on=False, zorder=3))
            ax.text(cx, cy - 0.145, f'{r.mag:.1f}', ha='center', va='center', fontsize=5.0,
                    color=INK if r.mag >= 6 else MUTED, fontweight='semibold' if r.mag >= 6 else 'normal')
        kc = f'{int(g.drop_c.dropna().astype(bool).sum())}/{int(g.drop_c.notna().sum())}'
        kp = f'{int(g.drop_p.dropna().astype(bool).sum())}/{int(g.drop_p.notna().sum())}'
        ax.text(x_res + 0.17, y + 0.22, kc, ha='right', va='center', fontsize=7.0, color=EST_COL['classical'], fontweight='semibold')
        ax.text(x_res + 0.22, y + 0.22, '|', ha='center', va='center', fontsize=6, color=MUTED)
        ax.text(x_res + 0.27, y + 0.22, kp, ha='left', va='center', fontsize=7.0, color=EST_COL['positive'], fontweight='semibold')
    yb = y0 + len(MAIN_LABELS) * rh
    ax.plot([0.02, W - 0.02], [yb, yb], color=INK2, lw=0.8)
    rise = []
    for lab in MAIN_LABELS:
        g = long[long.window == lab]
        rise.append(f"{FIG_OF[lab]} {int(g.rise_c.dropna().astype(bool).sum())}/{int(g.rise_c.notna().sum())} | {int(g.rise_p.dropna().astype(bool).sum())}/{int(g.rise_p.notna().sum())}")
    ax.text(0.08, yb + 0.13, 'Filled = b of the last window before the event ≥ 10% below the pre-sequence reference (median of the earlier '
            'windows of the panel); dashed outline = aftershock of an earlier marked event; grey = no valid window.', fontsize=4.6, color=MUTED, va='center')
    ax.text(0.08, yb + 0.26, 'Rise after (first post-event window ≥ 10% above the same reference), classical | positive: ' + '; '.join(rise) + '.',
            fontsize=4.8, color=MUTED, va='center')
    fig.savefig(path, bbox_inches='tight', pad_inches=0.04)
    plt.close(fig)
