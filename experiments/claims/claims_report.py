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

ROWS = [  # (key, label, claim group, base kind)
    ('C1_rate', 'Low-b nodes with high strain rate', 'C1', 'half'),
    ('C2_rate', 'High-b nodes with low strain rate', 'C2', 'half'),
    ('C1_regime', 'Contractional nodes that are low-b', 'C1', 'low'),
    ('C2_regime', 'Extensional nodes that are high-b', 'C2', 'high'),
    ('C1_mech', 'Thrust mechanisms in low-b nodes', 'C1', 'low'),
    ('C2_mech', 'Normal mechanisms in high-b nodes', 'C2', 'high'),
    ('C3_hf', 'High-heat-flow nodes that are high-b', 'C3', 'high'),
    ('C4a', 'M>=6 preceded by a b drop', 'C4', None),
    ('C4b', 'M>=6 followed by a b rise', 'C4', None),
    ('C4_alarm', 'False-alarm rate: b-drop alarms with no M>=5 within 1 yr', 'C4', None),
]


def _pct(k, n):
    return 100.0 * k / n if n else np.nan


def build_table(thr, gcmt, events, alarm_summary, min_cov_nodes=20, min_mech=10):
    """Long table: one row per (claim row, region, estimator) with k, n, pct, base."""
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
            # strain rate
            s = v[np.isfinite(v.ws_e2inv)]
            sl, sh = s[bcol] < b0, s[bcol] > b0
            hs = s.ws_e2inv > s.ws_e2inv.median()
            rec['C1_rate'] = ((sl & hs).sum(), sl.sum(), 50.0)
            rec['C2_rate'] = ((sh & ~hs).sum(), sh.sum(), 50.0)
            q1, q2 = s.ws_dilat.quantile([1 / 3, 2 / 3])
            contr, ext = s.ws_dilat <= q1, s.ws_dilat >= q2
            rec['C1_regime'] = ((contr & sl).sum(), contr.sum(), base_low)
            rec['C2_regime'] = ((ext & sh).sum(), ext.sum(), base_high)
            # mechanisms
            bn = nearest_node_value(nd, g.lon.to_numpy(), g.lat.to_numpy(), col=bcol) if len(g) else np.array([])
            th = (g['style'] == 'thrust').to_numpy() & np.isfinite(bn)
            no = (g['style'] == 'normal').to_numpy() & np.isfinite(bn)
            rec['C1_mech'] = ((bn[th] < b0).sum(), th.sum(), base_low) if th.sum() >= min_mech else (np.nan, th.sum(), base_low)
            rec['C2_mech'] = ((bn[no] > b0).sum(), no.sum(), base_high) if no.sum() >= min_mech else (np.nan, no.sum(), base_high)
            # fluids
            vh = v[np.isfinite(v.hf)]
            hh = vh.hf > vh.hf.median()
            rec['C3_hf'] = (((vh[bcol] > b0) & hh).sum(), hh.sum(), base_high) if (len(vh) >= min_cov_nodes and region != 'XSH') else (np.nan, 0, base_high)
            # temporal (clean M>=6 since 2008)
            e = events[(events.region == region) & (events.estimator == est) & (events['set'] == 'M6_2008')]
            e = e[~e.in_sequence.astype(bool)]
            d = e['drop'].dropna().astype(bool); r = e['rise'].dropna().astype(bool)
            rec['C4a'] = (d.sum(), len(d), np.nan) if len(d) else (np.nan, 0, np.nan)
            rec['C4b'] = (r.sum(), len(r), np.nan) if len(r) else (np.nan, 0, np.nan)
            a = alarm_summary[(alarm_summary.region == region) & (alarm_summary.estimator == est)].iloc[0]
            n_ep, n_false = a.n_episodes_clean, a.n_false_clean
            rec['C4_alarm'] = (n_false, n_ep, np.nan) if n_ep else (np.nan, 0, np.nan)      # false alarms / clean episodes
            for key, label, grp, _ in ROWS:
                k, n, base = rec[key]
                out.append(dict(row=key, label=label, claim=grp, region=region, estimator=est,
                                k=k, n=n, pct=_pct(k, n) if np.isfinite(k) else np.nan, base=base,
                                alarm_fraction=a.alarm_fraction if key == 'C4_alarm' else np.nan))
    return pd.DataFrame(out)


def table_markdown(tab):
    """Wide markdown: rows = claim rows, columns = regions, cell = 'k/n (p%) | k/n (p%)' with base rate."""
    lines = ['| Claim | Row | ' + ' | '.join(REGIONS) + ' | Mean % (regions tested) |', '|' + '---|' * (len(REGIONS) + 3)]
    for key, label, grp, basek in ROWS:
        cells = []; means = []
        for region in REGIONS:
            parts = []
            for est in ('classical', 'positive'):
                r = tab[(tab.row == key) & (tab.region == region) & (tab.estimator == est)].iloc[0]
                if np.isfinite(r.pct):
                    parts.append(f'{int(r.k)}/{int(r.n)} ({r.pct:.0f}%)'); means.append(r.pct)
                else:
                    parts.append('n/a')
            base = tab[(tab.row == key) & (tab.region == region)].base.iloc[0]
            cell = ' \\| '.join(parts)
            if basek and np.isfinite(base) and 'n/a' != parts[0] or (basek and 'n/a' != parts[1]):
                b2 = tab[(tab.row == key) & (tab.region == region)].base.tolist()
                cell += f'<br>base {b2[0]:.0f}% \\| {b2[1]:.0f}%' if basek != 'half' else '<br>base 50%'
            cells.append(cell)
        n_tested = sum(1 for region in REGIONS if np.isfinite(tab[(tab.row == key) & (tab.region == region)].pct).any())
        lines.append(f'| {grp} | {label} | ' + ' | '.join(cells) + f' | {np.mean(means):.0f}% ({n_tested}) |' if means else f'| {grp} | {label} | ' + ' | '.join(cells) + ' | n/a |')
    af = tab[tab.row == 'C4_alarm'].pivot(index='estimator', columns='region', values='alarm_fraction')
    lines.append('')
    lines.append('Alarm fraction of 2008-2023 (classical | positive): ' + '; '.join(
        f"{r} {100*af.loc['classical', r]:.0f}% | {100*af.loc['positive', r]:.0f}%" for r in REGIONS))
    return '\n'.join(lines)


# ------------------------------------------------------------------------------- figures
def fig_summary(tab, path):
    """Claims x regions matrix; each cell split classical | positive; colour = % supporting, text = k/n."""
    nrow, ncol = len(ROWS), len(REGIONS)
    fig, ax = plt.subplots(figsize=(7.2, 0.42 * nrow + 1.7))
    ax.set_xlim(0, ncol); ax.set_ylim(nrow, 0); ax.axis('off')
    norm = Normalize(0, 100)
    REF = {'half': 'ref. 50%', 'low': 'ref. share of low-b nodes', 'high': 'ref. share of high-b nodes', None: ''}
    for i, (key, label, grp, basek) in enumerate(ROWS):
        for j, region in enumerate(REGIONS):
            for h, est in enumerate(('classical', 'positive')):
                r = tab[(tab.row == key) & (tab.region == region) & (tab.estimator == est)].iloc[0]
                x0 = j + 0.5 * h
                if np.isfinite(r.pct):
                    shade = 100 - r.pct if key == 'C4_alarm' else r.pct      # false-alarm row: colour = share of alarms that were followed by an M>=6
                    col = (CM_SEQ if est == 'classical' else CM_SEQ_RED)(norm(shade)); txtc = 'white' if shade > 55 else INK
                    ax.add_patch(Rectangle((x0 + 0.02, i + 0.04), 0.46, 0.92, fc=col, ec='none'))
                    ax.text(x0 + 0.25, i + 0.40, f'{r.pct:.0f}%', ha='center', va='center', fontsize=7.5, color=txtc, fontweight='semibold')
                    kn = f'{int(r.k)}/{int(r.n)}'
                    ax.text(x0 + 0.25, i + 0.75, kn, ha='center', va='center', fontsize=5.6 if len(kn) <= 7 else 4.6, color=txtc)
                else:
                    ax.add_patch(Rectangle((x0 + 0.02, i + 0.04), 0.46, 0.92, fc='#f3f2ee', ec='none'))
                    ax.text(x0 + 0.25, i + 0.55, 'n/a', ha='center', va='center', fontsize=6.5, color=MUTED)
        ax.text(-0.08, i + 0.40 if (basek or key == 'C4_alarm') else i + 0.52, label, ha='right', va='center', fontsize=7, color=INK)
        if basek:
            ax.text(-0.08, i + 0.74, REF[basek], ha='right', va='center', fontsize=5.5, color=MUTED)
        if key == 'C4_alarm':
            ax.text(-0.08, i + 0.74, 'text: false alarms / clean episodes; colour: share followed by an M\u22655', ha='right', va='center', fontsize=5.5, color=MUTED)
        if i in (2, 6, 7):
            ax.plot([-0.0, ncol], [i, i], color='white', lw=3)
    for j, region in enumerate(REGIONS):
        name = {'LMS': 'Longmen-\nshan', 'XSH': 'Xianshui-\nhe', 'WYN': 'Western\nYunnan', 'XJ': 'Xiaojiang', 'NCP': 'North\nChina Plain', 'NW': 'Northwest\nChina'}[region]
        ax.text(j + 0.5, -0.42, name, ha='center', va='center', fontsize=6.8, fontweight='semibold', color=INK, linespacing=0.95)
        ax.text(j + 0.25, -0.06, 'classical', ha='center', va='center', fontsize=5.3, color=EST_COL['classical'])
        ax.text(j + 0.75, -0.06, 'positive', ha='center', va='center', fontsize=5.3, color=EST_COL['positive'])
    for grp, y in (('Stress (C1, C2)', 3.0), ('Fluids (C3)', 6.5), ('Fore-/aftershocks (C4)', 8.5)):
        ax.text(ncol + 0.08, y, grp, rotation=270, ha='left', va='center', fontsize=7, color=INK2)
    for x, cm, lab in ((0.42, CM_SEQ, 'b-classical'), (0.70, CM_SEQ_RED, 'b-positive')):
        cax = fig.add_axes([x, 0.055, 0.22, 0.016])
        cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cm, norm=norm, orientation='horizontal')
        cb.set_label(f'{lab}: % supporting (small text: k/n)', fontsize=6.3, color=INK2); cb.ax.tick_params(labelsize=6, colors=MUTED); cb.outline.set_visible(False)
    fig.subplots_adjust(left=0.42, right=0.955, top=0.90, bottom=0.13)
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


def figS_C4(events, figrun, path, years=1.0):
    """Small multiples: b(t) of both estimators around each M>=6 (clean and in-sequence), with reference,
    pre-event and post-event windows marked."""
    ev = events[(events['set'] == 'M6_2008') & (events.estimator == 'classical')].sort_values(['region', 'time']).reset_index(drop=True)
    n = len(ev); ncol = 4; nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.2, 1.55 * nrow), sharex=False)
    for ax, (_, r) in zip(axes.ravel(), ev.iterrows()):
        t_ev = pd.Timestamp(r.time); lo, hi = t_ev - pd.Timedelta(days=365.25 * years), t_ev + pd.Timedelta(days=365.25 * years)
        for est, cls in (('classical', ClassicBValueEstimator), ('positive', BPositiveBValueEstimator)):
            s, _ = series_for_window(r.label, cls, figrun); s = add_reference(s)
            w = s[(s.t_last >= lo) & (s.t_last <= hi)]
            x = (w.t_last - t_ev).dt.total_seconds() / 86400
            ax.plot(x, w.b, color=EST_COL[est], lw=0.8)
            e2 = events[(events['set'] == 'M6_2008') & (events.estimator == est) & (events.region == r.region) & (pd.to_datetime(events.time) == t_ev)].iloc[0]
            if np.isfinite(e2.ref):
                ax.plot([-365.25 * years, 0], [e2.ref, e2.ref], color=EST_COL[est], lw=0.6, ls=(0, (3, 2)), alpha=0.8)
            if np.isfinite(e2.b_pre):
                ax.plot([-e2.pre_days, 0], [e2.b_pre, e2.b_pre], color=EST_COL[est], lw=2.2, alpha=0.9)
            if np.isfinite(e2.b_post):
                ax.plot([e2.post_start_days, e2.post_start_days + e2.post_days], [e2.b_post, e2.b_post], color=EST_COL[est], lw=2.2, alpha=0.9)
        ax.axvline(0, color=INK, lw=0.6)
        ax.set_xlim(-365.25 * years, 365.25 * years)
        tag = ' – in seq.' if r.in_sequence else ''
        ax.set_title(f'{r.region} {t_ev:%Y-%m-%d} M{r.mag:.1f}{tag}', fontsize=6.3, loc='left')
        ax.tick_params(labelsize=5.5); ax.grid(axis='y', color=GRID, lw=0.4)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    for ax in axes.ravel()[n:]:
        ax.axis('off')
    handles = [Line2D([], [], color=C_GR, lw=1.2, label='b-classical'), Line2D([], [], color=C_BPOS, lw=1.2, label='b-positive'),
               Line2D([], [], color=INK2, lw=0.8, ls=(0, (3, 2)), label='reference (mean of the preceding year)'),
               Line2D([], [], color=INK2, lw=2.2, label='pre-event / post-event window')]
    fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=6.5, bbox_to_anchor=(0.5, 0.0))
    fig.text(0.5, 0.025, 'days relative to the event', ha='center', fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path); plt.close(fig)
