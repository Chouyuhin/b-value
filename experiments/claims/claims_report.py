"""Table and figures of the claims test (reads the CSVs written by claims_test.ipynb / claims_lib)."""
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
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
    ('C1_rate', 'Low-b nodes with high strain rate', 'C1', 'half'),
    ('C2_rate', 'High-b nodes with low strain rate', 'C2', 'half'),
    ('C1_regime', 'Contractional nodes that are low-b', 'C1', 'low'),
    ('C2_regime', 'Extensional nodes that are high-b', 'C2', 'high'),
    ('C1_mech', 'Thrust mechanisms in low-b nodes', 'C1', 'low'),
    ('C2_mech', 'Normal mechanisms in high-b nodes', 'C2', 'high'),
    ('C3_hf', 'High-heat-flow nodes that are high-b', 'C3', 'high'),
]
C4_COLS = [  # temporal columns (key, label)
    ('drop', 'Marked events preceded by a b drop'),
    ('rise', 'Marked events followed by a b rise'),
    ('hit', 'Drops followed by an M>=4 within 7 d'),
]


def _pct(k, n):
    return 100.0 * k / n if n else np.nan


def build_table(thr, gcmt, min_cov_nodes=20, min_mech=10):
    """Long table of the spatial rows: one row per (claim row, region, estimator) with k, n, pct, base."""
    out = []
    for region in REGIONS:
        nd = pd.read_csv(OUT / f'nodes_{region}.csv')
        g = gcmt[in_box(gcmt.lon, gcmt.lat, BOXES[region]) & (gcmt.mw >= 5.0)]
        for est, bcol in (('classical', 'b_c'), ('positive', 'b_p')):
            b0 = thr.loc[region, 'b_c' if est == 'classical' else 'b_p']
            v = nd[np.isfinite(nd[bcol])]
            low, high = v[bcol] < b0, v[bcol] > b0
            base_low, base_high = 100 * low.mean(), 100 * high.mean()
            rec = {}
            s_ = v[np.isfinite(v.ws_e2inv)]
            sl, sh = s_[bcol] < b0, s_[bcol] > b0
            hs = s_.ws_e2inv > s_.ws_e2inv.median()
            rec['C1_rate'] = ((sl & hs).sum(), sl.sum(), 50.0)
            rec['C2_rate'] = ((sh & ~hs).sum(), sh.sum(), 50.0)
            q1, q2 = s_.ws_dilat.quantile([1 / 3, 2 / 3])
            contr, ext = s_.ws_dilat <= q1, s_.ws_dilat >= q2
            rec['C1_regime'] = ((contr & sl).sum(), contr.sum(), base_low)
            rec['C2_regime'] = ((ext & sh).sum(), ext.sum(), base_high)
            bn = nearest_node_value(nd, g.lon.to_numpy(), g.lat.to_numpy(), col=bcol) if len(g) else np.array([])
            th = (g['style'] == 'thrust').to_numpy() & np.isfinite(bn)
            no = (g['style'] == 'normal').to_numpy() & np.isfinite(bn)
            rec['C1_mech'] = ((bn[th] < b0).sum(), th.sum(), base_low) if th.sum() >= min_mech else (np.nan, th.sum(), base_low)
            rec['C2_mech'] = ((bn[no] > b0).sum(), no.sum(), base_high) if no.sum() >= min_mech else (np.nan, no.sum(), base_high)
            vh = v[np.isfinite(v.hf)]
            hh = vh.hf > vh.hf.median()
            rec['C3_hf'] = (((vh[bcol] > b0) & hh).sum(), hh.sum(), base_high) if (len(vh) >= min_cov_nodes and region != 'XSH') else (np.nan, 0, base_high)
            for key, label, grp, _ in ROWS:
                k, n, base = rec[key]
                out.append(dict(row=key, label=label, claim=grp, region=region, estimator=est,
                                k=k, n=n, pct=_pct(k, n) if np.isfinite(k) else np.nan, base=base))
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
        cells = []
        for region in REGIONS:
            parts = []
            for est in ('classical', 'positive'):
                r = tab[(tab.row == key) & (tab.region == region) & (tab.estimator == est)].iloc[0]
                parts.append(f'{int(r.k)}/{int(r.n)} ({r.pct:.0f}%)' if np.isfinite(r.pct) else 'n/a')
            b2 = tab[(tab.row == key) & (tab.region == region)].base.tolist()
            ref = 'ref. 50%' if basek == 'half' else f'ref. {b2[0]:.0f}/{b2[1]:.0f}%'
            cells.append(' \\| '.join(parts) + ('' if parts == ['n/a', 'n/a'] else f'<br>{ref}'))
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
                parts.append(f'{int(r.k)}/{int(r.n)} ({r.pct:.0f}%)' if np.isfinite(r.pct) else f'n/a (0/{int(r.n)})' if key == 'hit' else 'n/a')
            cell = ' \\| '.join(parts)
            if key == 'hit' and fig != 'Total':
                rs = c4[(c4.col == 'hit') & (c4.fig == fig)].hit_ref_share.iloc[0]
                cell += f'<br>random date: {100 * rs:.0f}%'
            cells.append(cell)
        r0 = c4[(c4.col == 'hit') & (c4.fig == fig)]
        lines.append(f'| {fig} {r0.window.iloc[0]} | {int(r0.n_marked.iloc[0])} | ' + ' | '.join(cells) + ' |')
    lines.append('')
    lines.append('Drops = alarm episodes (runs of windows >= 10% below their reference). A drop counts as a hit when an M>=4 occurs in the '
                 'window region within the following 7 days; false-alarm rate = 100% minus the hit rate; "random date" = share of days in the window '
                 'followed by an M>=4 within 7 days (chance reference).')
    return '\n'.join(lines)


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


def fig_summary(tab, c4, path):
    """(a) spatial claims x regions; (b) temporal claims x main-text windows. Cells split classical | positive;
    colour = % supporting the claim (false-alarm column shaded by its complement), text = k/n."""
    nrow, ncol = len(ROWS), len(REGIONS)
    wins = [FIG_OF[l] for l in MAIN_LABELS]
    nrow2, ncol2 = len(wins), len(C4_COLS)
    fig = plt.figure(figsize=(7.2, 0.42 * (nrow + nrow2) + 3.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[nrow + 1.2, nrow2 + 1.9], hspace=0.42, left=0.40, right=0.955, top=0.95, bottom=0.10)
    REF = {'half': 'ref. 50%', 'low': 'ref. share of low-b nodes', 'high': 'ref. share of high-b nodes'}
    # (a)
    ax = fig.add_subplot(gs[0]); ax.set_xlim(0, ncol); ax.set_ylim(nrow, 0); ax.axis('off')
    for i, (key, label, grp, basek) in enumerate(ROWS):
        for j, region in enumerate(REGIONS):
            for h, est in enumerate(('classical', 'positive')):
                r = tab[(tab.row == key) & (tab.region == region) & (tab.estimator == est)].iloc[0]
                _cell(ax, j + 0.5 * h, i, 0.5, 1, r.pct, est, r.k, r.n)
        ax.text(-0.08, i + 0.40, label, ha='right', va='center', fontsize=7, color=INK)
        ax.text(-0.08, i + 0.74, REF[basek], ha='right', va='center', fontsize=5.5, color=MUTED)
        if i in (2, 6):
            ax.plot([0, ncol], [i, i], color='white', lw=3)
    for j, region in enumerate(REGIONS):
        name = {'LMS': 'Longmen-\nshan', 'XSH': 'Xianshui-\nhe', 'WYN': 'Western\nYunnan', 'XJ': 'Xiaojiang', 'NCP': 'North\nChina Plain', 'NW': 'Northwest\nChina'}[region]
        ax.text(j + 0.5, -0.45, name, ha='center', va='center', fontsize=6.8, fontweight='semibold', color=INK, linespacing=0.95)
        ax.text(j + 0.25, -0.07, 'classical', ha='center', va='center', fontsize=5.3, color=EST_COL['classical'])
        ax.text(j + 0.75, -0.07, 'positive', ha='center', va='center', fontsize=5.3, color=EST_COL['positive'])
    ax.text(-0.08, -0.6, 'a  Spatial claims, 2008\u20132023 maps', ha='right', va='center', fontsize=8, fontweight='bold', color=INK)
    for grp, y in (('Stress (C1, C2)', 3.0), ('Fluids (C3)', 6.5)):
        ax.text(ncol + 0.08, y, grp, rotation=270, ha='left', va='center', fontsize=7, color=INK2)
    # (b)
    ax = fig.add_subplot(gs[1]); ax.set_xlim(0, ncol2); ax.set_ylim(nrow2, 0); ax.axis('off')
    for i, figname in enumerate(wins):
        for j, (key, label) in enumerate(C4_COLS):
            for h, est in enumerate(('classical', 'positive')):
                r = c4[(c4.col == key) & (c4.fig == figname) & (c4.estimator == est)].iloc[0]
                _cell(ax, j + 0.5 * h, i, 0.5, 1, r.pct, est, r.k, r.n)
        r0 = c4[(c4.col == 'hit') & (c4.fig == figname)]
        lab = f"{figname}  {REGION_NAME[r0.window.iloc[0].split('_')[0]]} {r0.window.iloc[0].split('_')[1]}"
        ax.text(-0.08, i + 0.40, lab, ha='right', va='center', fontsize=7, color=INK)
        sub = f"{int(r0.n_marked.iloc[0])} marked events" + (f" \u00b7 random date followed by M\u22654 within 7 d: {100 * r0.hit_ref_share.iloc[0]:.0f}%" if np.isfinite(r0.hit_ref_share.iloc[0]) else '')
        ax.text(-0.08, i + 0.74, sub, ha='right', va='center', fontsize=5.2, color=MUTED)
    for j, (key, label) in enumerate(C4_COLS):
        short = {'drop': 'Marked events\npreceded by\na b drop', 'rise': 'Marked events\nfollowed by\na b rise',
                 'hit': 'Drops followed\nby an M\u22654\nwithin 7 d'}[key]
        ax.text(j + 0.5, -0.55, short, ha='center', va='center', fontsize=6.2, fontweight='semibold', color=INK, linespacing=0.95)
        ax.text(j + 0.25, -0.07, 'classical', ha='center', va='center', fontsize=5.3, color=EST_COL['classical'])
        ax.text(j + 0.75, -0.07, 'positive', ha='center', va='center', fontsize=5.3, color=EST_COL['positive'])
    ax.text(-0.08, -0.9, 'b  Temporal claims (C4), main-text windows', ha='right', va='center', fontsize=8, fontweight='bold', color=INK)
    ax.text(-0.08, nrow2 + 0.35, 'Drops: text = hits / drops; false-alarm rate = 100% \u2212 hit rate.', fontsize=5.5, color=MUTED, va='center')
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
               Line2D([], [], color=INK2, lw=0.8, ls=(0, (3, 2)), label='reference (mean of the earlier windows of the panel)'),
               Line2D([], [], color=INK2, lw=2.2, label='pre-event / post-event window')]
    fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=6.5, bbox_to_anchor=(0.5, 0.0))
    fig.text(0.5, 0.025, 'days relative to the event', ha='center', fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path); plt.close(fig)
