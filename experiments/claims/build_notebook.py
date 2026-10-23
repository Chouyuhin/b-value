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
after (C4b) mainshocks, plus the alarm/false-alarm view of C4.
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
N = 500): the share of windows flagged as alarms by the rule (≥10% below the preceding-year mean and ≥2 std).
""")
code("""
rng = np.random.default_rng(0)
n = 60000; mags = bin_to_precision(1.9 + rng.exponential(1 / np.log(10), n), 0.1)
times = (pd.Timestamp('2008-01-01') + pd.to_timedelta(np.sort(rng.uniform(0, 15 * 365.25, n)), unit='D')).to_numpy()
for cls in (ClassicBValueEstimator, BPositiveBValueEstimator):
    s = add_reference(b_series(mags, times, 1.9, 500, cls))
    print(f'{cls.__name__}: alarm windows {100 * s.alarm.mean():.1f}%  (≥10% drop alone: {100 * (s.rel <= -0.10).mean():.1f}%)')
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
    thr[region] = dict(mc=mc, n=len(df), b_c=bc.b_value, b_p=bp.b_value)
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
## 3. C4: events, pre/post windows, alarm episodes (events_C4.csv, alarms_C4.csv)

Rules: pre-event window = last window ending before the event; reference = mean of the windows ending in the year before the
pre-event window starts; drop = ≥10% below the reference; post window = first N events after a 7-day
blind time; rise = ≥10% above the pre-event b. Alarm = window ≥10% below the running reference (same rule); contiguous windows form one
episode; alarm duration 1 yr; false-alarm target M≥5 in the box (event rows: M≥6). Events and episodes inside the Gardner–Knopoff time window of a larger
earlier event are flagged `in_sequence` and excluded from the clean counts.
""")
code("""
EV, EP, SM = [], [], []
for region in REGIONS:
    for cls in (ClassicBValueEstimator, BPositiveBValueEstimator):
        et = event_table(region, cls, figrun, min_mag=6.0); et['set'] = 'M6_2008'
        ev5 = targets_main_windows(region, 5.0)
        et5 = event_table(region, cls, figrun, ev=ev5) if len(ev5) else pd.DataFrame()
        if len(et5): et5['set'] = 'M5_main'
        ep, summ = alarm_table(region, cls, figrun, target_mag=ALARM_TARGET); summ.pop('per_tile')
        EV += [et, et5]; EP.append(ep); SM.append(summ)
events = pd.concat(EV, ignore_index=True); alarms = pd.concat(EP, ignore_index=True); alarm_summary = pd.DataFrame(SM)
events.to_csv(OUT / 'events_C4.csv', index=False); alarms.to_csv(OUT / 'alarms_C4.csv', index=False); alarm_summary.to_csv(OUT / 'alarms_C4_summary.csv', index=False)
cols = ['region', 'estimator', 'time', 'mag', 'label', 'n_m', 'mc', 'in_sequence', 'b_pre', 'std_pre', 'pre_days', 'ref', 'drop_rel', 'drop', 'b_post', 'std_post', 'post_start_days', 'post_days', 'rise_rel', 'rise']
events[events['set'] == 'M6_2008'][cols].round(3)
""")
code("""
alarm_summary
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
tab = build_table(thr, gcmt, events, alarm_summary); tab.to_csv(OUT / 'claims_table.csv', index=False)
mdtxt = table_markdown(tab); open(OUT / 'claims_table.md', 'w').write(mdtxt); display(Markdown(mdtxt))
""")
code("""
fig_summary(tab, OUT / 'fig_summary.pdf')
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
