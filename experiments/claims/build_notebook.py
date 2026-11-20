"""Writes claims_test.ipynb (run once; the notebook itself is then executed with jupyter nbconvert)."""
import json
from pathlib import Path

cells = []


def md(s):
    cells.append({'cell_type': 'markdown', 'metadata': {}, 'source': s.strip('\n')})


def code(s):
    cells.append({'cell_type': 'code', 'metadata': {}, 'execution_count': None, 'outputs': [], 'source': s.strip('\n')})


md("""
# Four literature b-value claims tested on the six sub-regions

Percentages of units (grid nodes, GCMT mechanisms, M≥6 events, alarm episodes) that support each claim, with the paper's own
settings (N, Mc, dmc = 0.1, grid radius 0.5°, N per node). No hypothesis tests. Analysis functions live in `claims_lib.py`
(paper helpers copied verbatim from `bvalue4longmenshan.ipynb`), table and figures in `claims_report.py`.

Claims: **C1** low b ↔ high stress, **C2** high b ↔ low stress, **C3** high b ↔ fluids, **C4** b drops before (C4a) and rises
after (C4b) the marked events.
""")
code("""
import os, time, warnings
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import matplotlib; import matplotlib.pyplot as plt
from IPython.display import Markdown, display
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
from claims_lib import *
from claims_report import *
OUT.mkdir(exist_ok=True)
figrun = load_figrun()
figrun[['region', 't_start', 't_end', 'n_events', 'mc', 'n_m', 'b_classical', 'b_classical_std', 'b_positive', 'b_positive_std']]
""")

md("""
## 1. Validation: reproduce the paper's whole-window values and the SeismoStats sliding window
""")
code("""
from seismostats.plots import plot_b_series_constant_nm
for label in ['LMS_2006-2010', 'XSH_2022-2023', 'NCP_2008-2023']:
    t0, t1 = window_span(label); df = load_region(figrun.loc[label, 'region'], t0, t1)
    mc = float(estimate_mc_maxc(df['mag'].to_numpy(), fmd_bin=0.1)[0])
    ec = b_window(df['mag'], df['time'], figrun.loc[label, 'mc'], ClassicBValueEstimator)
    ep = b_window(df['mag'], df['time'], figrun.loc[label, 'mc'], BPositiveBValueEstimator)
    assert len(df) == figrun.loc[label, 'n_events'] and mc == figrun.loc[label, 'mc']
    assert round(ec.b_value, 3) == figrun.loc[label, 'b_classical'] and round(ep.b_value, 3) == figrun.loc[label, 'b_positive']
    assert round(ec.std, 3) == figrun.loc[label, 'b_classical_std'] and round(ep.std, 3) == figrun.loc[label, 'b_positive_std']
    print(f'{label}: n={len(df)}, Mc={mc}, b_c={ec.b_value:.3f}±{ec.std:.3f}, b_p={ep.b_value:.3f}±{ep.std:.3f}  == figrun')
# sliding window vs seismostats.plots.plot_b_series_constant_nm (XSH 2022-2023, N=300)
label = 'XSH_2022-2023'; t0, t1 = window_span(label); df = load_region('XSH', t0, t1); mc = figrun.loc[label, 'mc']; n_m = int(figrun.loc[label, 'n_m'])
for cls in (ClassicBValueEstimator, BPositiveBValueEstimator):
    s = b_series(df['mag'], df['time'], mc, n_m, cls)
    fig, ax = plt.subplots(); kw = dict(dmc=0.1) if cls is BPositiveBValueEstimator else {}
    plot_b_series_constant_nm(df['mag'].to_numpy(), 0.1, mc, df['time'].to_numpy(), n_m=n_m, x_variable=df['time'].to_numpy(),
                              plot_technique='right', ax=ax, b_method=cls, **kw)
    y = ax.lines[0].get_ydata(); y = y[~np.isnan(y)]; plt.close(fig)
    assert np.allclose(y, s['b'].to_numpy()); print(cls.__name__, 'series reproduced,', len(s), 'windows')
""")
md("""
Sanity figure for the C4 rule on a constant-b synthetic catalogue (b = 1, Mc 1.9, magnitudes binned to 0.1, Poisson times,
N = 500): the share of windows ≥10% below the whole-period b (the alarm rule).
""")
code("""
rng = np.random.default_rng(0)
n = 60000; mags = bin_to_precision(1.9 + rng.exponential(1 / np.log(10), n), 0.1)
times = (pd.Timestamp('2008-01-01') + pd.to_timedelta(np.sort(rng.uniform(0, 15 * 365.25, n)), unit='D')).to_numpy()
for cls in (ClassicBValueEstimator, BPositiveBValueEstimator):
    s = b_series(mags, times, 1.9, 500, cls); b_ref = b_window(mags, times, 1.9, cls).b_value
    print(f'{cls.__name__}: windows ≥10% below the whole-period b: {100 * ((s.b - b_ref) / b_ref <= -0.10).mean():.1f}%')
""")

md("""
## 2. Spatial b-maps 2008–2023, strain rate, heat flow (nodes_<REG>.csv)

Threshold b₀ = whole-window b of the same region, period and estimator (Mc = MAXC + 0.2 of the 2008–2023 catalogue).
Strain-rate tensor from Wang & Shen (2020) velocities (`ws_`), and from the SCB-fixed Sichuan–Yunnan field (`ds_`) for the
four SW regions. Node spacing 0.05°; N per node and the 0.5° radius as in the paper.
""")
code("""
ws = load_gnss_ws(); ds = load_gnss_ds01(); hf = load_heatflow()
thr = {}
for region in REGIONS:
    t = time.time(); df = load_region(region, MAP_T0, MAP_T1)
    mc = float(estimate_mc_maxc(df['mag'].to_numpy(), fmd_bin=0.1)[0])
    bc = b_window(df['mag'], df['time'], mc, ClassicBValueEstimator); bp = b_window(df['mag'], df['time'], mc, BPositiveBValueEstimator)
    thr[region] = dict(mc=mc, n=len(df), b_c=bc.b_value, b_p=bp.b_value, n_c=bc.n, n_p=bp.n)   # n_c / n_p = sample sizes of the two regional fits (Utsu test)
    nc = bmap_nodes(df, mc, N_PER_GRID[region], ClassicBValueEstimator, BOXES[region])
    npos = bmap_nodes(df, mc, N_PER_GRID[region], BPositiveBValueEstimator, BOXES[region])
    nodes = nc.rename(columns={'b': 'b_c', 'std': 'std_c', 'n': 'n_c'}).drop(columns=['radius'])
    nodes['b_p'], nodes['std_p'], nodes['n_p'], nodes['radius'] = npos['b'], npos['std'], npos['n'], nc['radius']
    nodes = pd.concat([nodes, strain_at(nodes, ws).add_prefix('ws_')], axis=1)
    if region in ('LMS', 'XSH', 'WYN', 'XJ'):
        nodes = pd.concat([nodes, strain_at(nodes, ds).add_prefix('ds_')], axis=1)
    h = heatflow_at(nodes, hf); nodes['hf'], nodes['n_hf'] = h['hf'], h['n_hf']
    nodes['region'] = region
    nodes.to_csv(OUT / f'nodes_{region}.csv', index=False)
    v = nodes[np.isfinite(nodes.b_c)]
    print(f"{region}: Mc {mc}, n {len(df)}, b0 {bc.b_value:.3f} | {bp.b_value:.3f}; valid nodes {len(v)}, with strain {np.isfinite(v.ws_e2inv).sum()}, "
          f"median strain {np.nanmedian(v.ws_e2inv):.0f} nstrain/yr, heat-flow nodes {np.isfinite(v.hf).sum()}  ({time.time() - t:.0f}s)")
thr = pd.DataFrame(thr).T; thr.to_csv(OUT / 'thresholds_2008-2023.csv'); thr
""")

md("""
## 3. C4 on the main-text windows (Figs 4–10): marked events, alarm episodes (events_C4.csv, alarms_C4.csv)

Unit = each main-text window with its own N and Mc (the paper's settings); targets = the events marked in panel (d)
(fifth-largest rule). Rules follow the published definitions (revised 2026-10-29):

- classical series: every window re-cut at max(panel Mc, MAXC + 0.2 of the window), ≥ 50 events above it
  (Gulia & Wiemer 2019, Methods steps 5–6a; Gulia et al. 2018); b-positive: panel Mc, dmc = 0.1, no re-cut and no
  exclusion (van der Elst 2021; Lippiello & Petrillo 2024).
- reference of a window = median of the earlier windows of the panel that share no event with it, leaving out windows
  that end inside the Gardner–Knopoff time window of an earlier marked event (Gulia & Wiemer 2019, 6a; Gulia et al. 2018);
  fewer than 5 eligible windows → the region's 2008–2023 whole-period b (Gulia & Wiemer 2019, 6b).
- drop = last window before the event ≥ 10% below its reference; rise = first post-event window ≥ 10% above the same
  reference (Gulia & Wiemer 2019, step 7). The classical post window starts when the per-window Mc is back at the panel
  floor (Gulia et al. 2018: events removed until Mc(t) matches the pre-cut level); the b-positive post window is the
  first window whose events all follow the event (van der Elst 2021).
- alarm episodes (runs of windows ≥ 10% below the reference) and their 7-day M≥4 follow-up are still written to
  alarms_C4.csv, but the false-alarm column was removed from the table on 2026-10-29: a false-alarm count needs an a priori
  target and validity window (Gulia et al. 2024 use a larger event within 3 months and the G-K radius, a question unrelated
  to this section), which regional windows cannot supply without arbitrary choices.
""")
code("""
# the fast MAXC used for the per-window re-cut must equal seismostats' estimate_mc_maxc
rng = np.random.default_rng(0)
for label in MAIN_LABELS:
    t0, t1 = window_span(label); df_ = load_region(figrun.loc[label, 'region'], t0, t1); n_ = int(figrun.loc[label, 'n_m'])
    for i0 in rng.integers(0, max(1, len(df_) - n_), 20):
        w_ = df_['mag'].to_numpy()[i0:i0 + n_]
        assert np.isclose(maxc(w_), estimate_mc_maxc(w_, fmd_bin=0.1)[0]), (label, i0)
print('maxc == estimate_mc_maxc on 140 random windows')
EV, EP, SM = [], [], []
for label in MAIN_LABELS:
    for cls in (ClassicBValueEstimator, BPositiveBValueEstimator):
        ev, ep, sm = window_c4(label, cls, figrun); EV.append(ev); EP.append(ep); SM.append(sm)
events = pd.concat(EV, ignore_index=True); alarms = pd.concat(EP, ignore_index=True); c4_summary = pd.DataFrame(SM)
events.to_csv(OUT / 'events_C4.csv', index=False); alarms.to_csv(OUT / 'alarms_C4.csv', index=False); c4_summary.to_csv(OUT / 'alarms_C4_summary.csv', index=False)
c4_summary[['fig', 'window', 'estimator', 'n_m', 'mc', 'b_ww', 'b_regional', 'n_windows', 'n_windows_valid', 'n_windows_in_seq', 'n_marked', 'n_ref_regional',
            'drop_k', 'drop_n', 'rise_k', 'rise_n', 'n_episodes', 'n_hits', 'hit_ref_share']].round(3)
""")
code("""
events[['fig', 'estimator', 'number', 'time', 'mag', 'initiator', 'ref', 'ref_source', 'n_ref', 'b_pre', 'mc_pre', 'pre_days', 'drop_rel', 'drop',
        'b_post', 'mc_post', 'post_start_days', 'post_days', 'rise_rel', 'rise']].round(3)
""")
md("""
SI sensitivity: b-positive with dmc = 0.2 (van der Elst 2021's default; the paper's panel (c) checks 0.1/0.2/0.4).
""")
code("""
SM2 = [window_c4(label, BPositiveBValueEstimator, figrun, dmc=0.2)[2] for label in MAIN_LABELS]
c4_dmc02 = pd.DataFrame(SM2); c4_dmc02.to_csv(OUT / 'alarms_C4_summary_dmc02.csv', index=False)
cmp = c4_summary[c4_summary.estimator == 'positive'].set_index('fig')[['drop_k', 'drop_n', 'rise_k', 'rise_n', 'n_episodes', 'n_hits']].join(
    c4_dmc02.set_index('fig')[['drop_k', 'drop_n', 'rise_k', 'rise_n', 'n_episodes', 'n_hits']], rsuffix='_dmc02')
cmp
""")

md("""
## 4. GCMT mechanisms and SI rows (epicentre b, depth)
""")
code("""
gcmt = parse_gcmt()
lines = []
for region in REGIONS:
    g = gcmt[in_box(gcmt.lon, gcmt.lat, BOXES[region]) & (gcmt.mw >= 5.0)]
    vc = g['style'].value_counts()
    lines.append(dict(region=region, n=len(g), **{k: int(vc.get(k, 0)) for k in ('thrust', 'normal', 'strike-slip', 'oblique')},
                      dominant=vc.index[0] if len(vc) else '', dominant_pct=100 * vc.iloc[0] / len(g) if len(g) else np.nan))
gcmt_summary = pd.DataFrame(lines); gcmt_summary.to_csv(OUT / 'gcmt_summary.csv', index=False); gcmt_summary
""")
code("""
EPI, DEP = [], []
for region in REGIONS:
    ev = targets(region, 6.0)
    if len(ev): EPI.append(epicentre_rows(region, ev, figrun))
    bins, summ = depth_rows(region, float(thr.loc[region, 'mc'])); DEP.append(bins)
    print(region, 'depth bins with b decreasing (classical | positive):', summ['all_dec_c'], '|', summ['all_dec_p'])
epi = pd.concat(EPI, ignore_index=True); epi.to_csv(OUT / 'epicentres_SI.csv', index=False)
pd.concat(DEP, ignore_index=True).to_csv(OUT / 'depth_SI.csv', index=False)
epi[[c for c in ['region', 'time', 'mag', 'mc', 'n', 'b0_c', 'b_epi_c', 'low_c', 'pct_c', 'b0_p', 'b_epi_p', 'low_p', 'pct_p', 'note'] if c in epi]].round(3)
""")

md("""
## 5. Table and figures
""")
code("""
tab = build_table(thr, gcmt); tab.to_csv(OUT / 'claims_table.csv', index=False)
c4 = build_c4_table(SM); c4.to_csv(OUT / 'claims_table_C4.csv', index=False)
mdtxt = table_markdown(tab, c4); open(OUT / 'claims_table.md', 'w').write(mdtxt); display(Markdown(mdtxt))
c4_long, c4_md = c4_event_table(events)   # Form B (chosen 2026-10-29): one row per window, one cell per marked event
c4_long.to_csv(OUT / 'table_C4_events.csv', index=False); open(OUT / 'table_C4_events.md', 'w').write(c4_md); display(Markdown(c4_md))
""")
code("""
fig_summary(tab, c4, OUT / 'fig_summary.pdf', events=events)
fig_table(tab, events, OUT / 'table_main.pdf')   # main-text table (spatial rows + C4 Form B) rendered as one figure, 2026-10-30
figS_C12(thr, gcmt, OUT / 'figS_C12.pdf')
figS_C3(thr, hf, OUT / 'figS_C3.pdf')
figS_C4(events, figrun, OUT / 'figS_C4.pdf')
sorted(p.name for p in OUT.iterdir())
""")

for i, c in enumerate(cells):
    c['id'] = f'cell-{i:02d}'
nb = {'cells': cells, 'metadata': {'kernelspec': {'display_name': 'gmtenv', 'language': 'python', 'name': 'gmtenv'},
                                   'language_info': {'name': 'python', 'version': '3.12.5'}},
      'nbformat': 4, 'nbformat_minor': 5}
Path(__file__).with_name('claims_test.ipynb').write_text(json.dumps(nb, indent=1, ensure_ascii=False))
print('written claims_test.ipynb with', len(cells), 'cells')
