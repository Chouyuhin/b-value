"""Foreshock traffic-light system (FTLS) of Gulia et al. (2024, SRL 95:3579, doi:10.1785/0220240163) applied to the
M >= 6 events of the six sub-regions since 2008, with the CENC catalogue.

Protocol (steps as numbered in the paper) and the deviations forced by our data:
1. T0 = 2008-01-01 (start of the high-quality phase), depth <= 40 km, M_c,overall = 1.0, Beijing time throughout.
2. Trigger = any M >= 6 in a sub-region box; analysis starts 1 h later.
3. Volume = fault box: nodal planes from GCMT, length/width from Wells & Coppersmith (1994) by Frohlich style,
   extended 3 km; the plane with more events in the first hour is the active one (ties: first 24 h).
   Deviation: the box is the map-view projection through the epicentre; CENC depths are not used.
4. b_ref = mean of b-positive over 150-event moving windows in the box from T0 to the mainshock (events after an
   earlier M >= 6 inside the box are excluded); with < 150 events, the 150 events nearest the epicentre before the
   mainshock give a single b_ref. Deviation (Yuxin, 2026-10-26): MAXC + 0.2 precut on every window, reference included.
5. Post series: first hour excluded; 300-event windows moved by one event; MAXC + 0.2 precut; b-positive with
   >= 50 positive differences; db = 100 b / b_ref.
6. Status every 12 h = mean db of the period: green >= 110, yellow 100-110, orange 90-100, red < 90; any single
   db >= 110 ends the analysis green.
7. Final status at 72 h = mean of the six 12-h means; valid 3 months within the Gardner-Knopoff radius.
Scoring: red + larger event within 3 months/GK radius = true positive; red + none = false alarm; green + larger =
missed; green + none = correct negative; yellow/orange neutral; no status when < 300 events by 72 h (then the time
of the first estimate and its status are reported as 'late').
Side row: pre-mainshock drop = last 150-event window before the mainshock with db < 90.
"""
import numpy as np
import pandas as pd
from seismostats.analysis import BPositiveBValueEstimator, estimate_mc_maxc

from claims_lib import (BOXES, DATA, EXP, REGIONS, DELTA_M, haversine_km, gardner_knopoff_km, targets)

T0 = pd.Timestamp('2008-01-01')
MC_OVERALL = 1.0
MAX_DEPTH = 40.0
N_PRE, N_POST, MIN_POS = 150, 300, 50
EXCL_HOURS = 1.0
STATUS_HOURS, FINAL_HOURS = 12.0, 72.0
VALID_DAYS = 91.0                         # three months
MARGIN_KM = 3.0
GREEN, YELLOW, ORANGE = 110.0, 100.0, 90.0
# Wells & Coppersmith (1994) Table 2A: log10 RLD = a + b M, log10 RW = c + d M (subsurface rupture length, width)
WC = {'R': (-2.42, 0.58, -1.61, 0.41), 'N': (-1.88, 0.50, -1.14, 0.35), 'SS': (-2.57, 0.62, -0.76, 0.27),
      'ALL': (-2.44, 0.59, -1.01, 0.32)}

_CENC = None


def load_cenc():
    """Full CENC catalogue (Mw-converted), Beijing time; cached."""
    global _CENC
    if _CENC is None:
        df = pd.read_csv(EXP.parent.parent / 'catalogs/CENCcat_Mw_final.txt', sep=r'\s+', header=None,
                         names=['y', 'mo', 'd', 'h', 'mi', 's', 'lat', 'lon', 'depth', 'type', 'mag'])
        df['time'] = pd.to_datetime(df[['y', 'mo', 'd', 'h', 'mi', 's']].rename(
            columns={'y': 'year', 'mo': 'month', 'd': 'day', 'h': 'hour', 'mi': 'minute', 's': 'second'}))
        df = df[['time', 'lon', 'lat', 'depth', 'mag']].sort_values('time', kind='stable').reset_index(drop=True)
        df['mag'] = np.round(df['mag'] / DELTA_M) * DELTA_M
        _CENC = df
    return _CENC


def load_gcmt_planes():
    """GCMT ndk -> DataFrame with both nodal planes and P/T/B plunges (UTC)."""
    import math
    rows = []
    for f in sorted((DATA / 'gcmt').glob('*.ndk')):
        L = open(f, errors='ignore').read().splitlines()
        for i in range(0, len(L) - 4, 5):
            l1 = L[i]; p = L[i + 4].split()
            try:
                t = pd.Timestamp(l1[5:15].replace('/', '-') + ' ' + l1[16:26].replace('60.0', '59.9'))
                lat, lon, dep = float(l1[27:33]), float(l1[34:41]), float(l1[42:47])
                ex = int(L[i + 3][0:2]); mw = (2 / 3) * (math.log10(float(p[10]) * 10 ** ex) - 16.1)
                s1, d1, r1, s2, d2, r2 = [float(x) for x in p[11:17]]
                pl_t, pl_b, pl_p = float(p[2]), float(p[5]), float(p[8])
            except (ValueError, IndexError):
                continue
            rows.append((t, lon, lat, dep, mw, s1, d1, r1, s2, d2, r2, pl_t, pl_b, pl_p))
    return pd.DataFrame(rows, columns=['time', 'lon', 'lat', 'depth', 'mw', 'strike1', 'dip1', 'rake1',
                                       'strike2', 'dip2', 'rake2', 'pl_t', 'pl_b', 'pl_p'])


def match_gcmt(ev, gcmt, max_s=300, max_deg=2.0):
    """GCMT entry of a CENC event (CENC time is UTC+8)."""
    t_utc = ev['time'] - pd.Timedelta(hours=8)
    dt = (gcmt['time'] - t_utc).abs().dt.total_seconds()
    c = gcmt[(dt < max_s) & ((gcmt['lon'] - ev['lon']).abs() < max_deg) & ((gcmt['lat'] - ev['lat']).abs() < max_deg)]
    if c.empty:
        return None
    return c.iloc[(c['time'] - t_utc).abs().argsort().iloc[0]]


def frohlich_style(pl_t, pl_b, pl_p):
    return 'R' if pl_t >= 50 else 'N' if pl_p >= 60 else 'SS' if pl_b >= 60 else 'ALL'


def fault_box(lon0, lat0, mag, strike, dip, style, margin_km=MARGIN_KM):
    """Map-view rectangle of a nodal plane through the epicentre: along-strike length RLD and the horizontal
    projection of the down-dip width RW cos(dip), each extended by margin_km. Returns a function (lon, lat) -> bool."""
    a, b, c, d = WC[style]
    L = 10 ** (a + b * mag); W = 10 ** (c + d * mag)
    half_l = L / 2 + margin_km; half_w = W * np.cos(np.radians(dip)) / 2 + margin_km
    phi = np.radians(strike)
    u = np.array([np.sin(phi), np.cos(phi)])            # along strike (x east, y north)
    v = np.array([np.cos(phi), -np.sin(phi)])           # dip direction (strike + 90)
    kx, ky = 111.32 * np.cos(np.radians(lat0)), 110.57

    def inside(lon, lat):
        dx = (np.asarray(lon, float) - lon0) * kx; dy = (np.asarray(lat, float) - lat0) * ky
        return (np.abs(dx * u[0] + dy * u[1]) <= half_l) & (np.abs(dx * v[0] + dy * v[1]) <= half_w)
    return inside, L, W


def bpos_precut(mags, times):
    """MAXC + 0.2 precut, then b-positive; NaN with fewer than MIN_POS positive differences."""
    mags = np.asarray(mags, float); times = np.asarray(times, dtype='datetime64[ns]')
    mc = float(estimate_mc_maxc(mags, fmd_bin=DELTA_M)[0])
    keep = mags >= mc - DELTA_M / 2
    if keep.sum() < 3:
        return np.nan, mc, 0
    est = BPositiveBValueEstimator()
    est.calculate(mags[keep], mc=mc, delta_m=DELTA_M, dmc=DELTA_M, times=times[keep])
    return (est.b_value if est.n >= MIN_POS else np.nan), mc, int(est.n)


GROW, N_CAP = 100, 3000


def moving_bpos(df, n, grow=GROW, n_cap=N_CAP):
    """b-positive (with precut) over moving windows ending at each event. The window holds n events; when the precut
    leaves fewer than MIN_POS positive differences it is widened backwards in steps of `grow` events (up to n_cap),
    because the CENC precut removes far more events than the local catalogues of Gulia et al. (2024).
    Value at the last event of the window; columns t, b, mc, n_pos, n_used."""
    mags = df['mag'].to_numpy(float); times = df['time'].to_numpy(dtype='datetime64[ns]')
    out = []
    for i in range(n - 1, len(df)):
        m = n
        while True:
            b, mc, npos = bpos_precut(mags[i - m + 1:i + 1], times[i - m + 1:i + 1])
            if npos >= MIN_POS or i - m + 1 == 0 or m >= n_cap:
                break
            m = min(m + grow, i + 1, n_cap)
        out.append((times[i], b, mc, npos, m))
    return pd.DataFrame(out, columns=['t', 'b', 'mc', 'n_pos', 'n_used'])


def ftls_event(ev, cat, gcmt_row, t0=T0, post_volume='box'):
    """Run the FTLS for one M >= 6 event. Returns (record dict, post-event series DataFrame)."""
    t_ev = ev['time']; M = float(ev['mag']); lon0, lat0 = float(ev['lon']), float(ev['lat'])
    rec = dict(region=ev['region'], time=t_ev, mag=M, lon=lon0, lat=lat0)
    near = cat[((cat['lon'] - lon0).abs() < 4) & ((cat['lat'] - lat0).abs() < 3) & (cat['depth'] <= MAX_DEPTH) &
               (cat['mag'] >= MC_OVERALL)]
    # ---- 3. fault box and active plane
    style = frohlich_style(gcmt_row['pl_t'], gcmt_row['pl_b'], gcmt_row['pl_p'])
    first_hour = near[(near['time'] > t_ev) & (near['time'] <= t_ev + pd.Timedelta(hours=1))]
    first_day = near[(near['time'] > t_ev) & (near['time'] <= t_ev + pd.Timedelta(hours=24))]
    boxes = []
    for k in (1, 2):
        inside, L, W = fault_box(lon0, lat0, M, gcmt_row[f'strike{k}'], gcmt_row[f'dip{k}'], style)
        boxes.append((k, inside, L, W, int(inside(first_hour['lon'], first_hour['lat']).sum()),
                      int(inside(first_day['lon'], first_day['lat']).sum())))
    boxes.sort(key=lambda x: (-x[4], -x[5]))
    k, inside, L, W, n1h, n24h = boxes[0]
    rec.update(style=style, plane=k, strike=gcmt_row[f'strike{k}'], dip=gcmt_row[f'dip{k}'], rupture_L_km=L,
               rupture_W_km=W, n_first_hour=n1h, n_first_day=n24h)
    inbox = near[inside(near['lon'], near['lat'])]
    # ---- 4. reference b
    pre = inbox[(inbox['time'] >= t0) & (inbox['time'] < t_ev)]
    prior_big = pre[pre['mag'] >= 6.0]
    if len(prior_big):
        pre = pre[pre['time'] < prior_big['time'].iloc[-1]]
    rec['n_ref_events'] = len(pre)
    if len(pre) >= N_PRE:
        ref = moving_bpos(pre, N_PRE)
        rec.update(ref_mode='box series', b_ref=float(ref['b'].mean()), n_ref_windows=int(ref['b'].notna().sum()))
        last = ref['b'].dropna()
        rec['db_pre_last'] = 100 * last.iloc[-1] / rec['b_ref'] if len(last) else np.nan
    else:
        cand = near[(near['time'] >= t0) & (near['time'] < t_ev)].copy()
        cand['dist'] = haversine_km(lon0, lat0, cand['lon'].to_numpy(), cand['lat'].to_numpy())
        m, b, npos, sel = N_PRE, np.nan, 0, cand.iloc[:0]
        while True:                                   # nearest 150, widened until the precut leaves >= 50 positive differences
            sel = cand.nsmallest(m, 'dist').sort_values('time')
            b, mc, npos = bpos_precut(sel['mag'], sel['time'])
            if npos >= MIN_POS or m >= len(cand) or m >= N_CAP:
                break
            m = min(m + GROW, len(cand), N_CAP)
        rec.update(ref_mode=f'nearest {len(sel)}', b_ref=b if npos >= MIN_POS else np.nan, n_ref_windows=1,
                   ref_radius_km=float(sel['dist'].max()) if len(sel) else np.nan, db_pre_last=np.nan)
    if np.isfinite(rec['b_ref']) and not (0.4 <= rec['b_ref'] <= 2.0):      # guard: implausible reference (tiny precut sample)
        rec.update(ref_note=f"b_ref {rec['b_ref']:.2f} implausible, reference rejected", b_ref=np.nan)
    rec['pre_drop'] = bool(rec['db_pre_last'] < ORANGE) if np.isfinite(rec.get('db_pre_last', np.nan)) else None
    # ---- 5. post-event series (volume: the fault box, or a circle of Gardner-Knopoff radius as a sensitivity variant)
    if post_volume == 'gk':
        d_all = haversine_km(lon0, lat0, near['lon'].to_numpy(), near['lat'].to_numpy())
        vol = near[d_all <= float(gardner_knopoff_km(M))]
    else:
        vol = inbox
    post = vol[vol['time'] > t_ev + pd.Timedelta(hours=EXCL_HOURS)]
    rec['n_post_72h'] = int((post['time'] <= t_ev + pd.Timedelta(hours=FINAL_HOURS)).sum())
    series = pd.DataFrame(columns=['t', 'b', 'mc', 'n_pos', 'hours', 'db'])
    if len(post) >= N_POST and np.isfinite(rec['b_ref']):
        series = moving_bpos(post.iloc[:min(len(post), 20000)], N_POST)
        series['hours'] = (series['t'] - t_ev.to_datetime64()) / np.timedelta64(1, 'h')
        series['db'] = 100 * series['b'] / rec['b_ref']
    # ---- 6-7. statuses
    def colour(x):
        return 'green' if x >= GREEN else 'yellow' if x >= YELLOW else 'orange' if x >= ORANGE else 'red'
    s72 = series[(series['hours'] <= FINAL_HOURS) & series['db'].notna()]
    means = []
    for q in range(int(FINAL_HOURS / STATUS_HOURS)):
        w = s72[(s72['hours'] > q * STATUS_HOURS) & (s72['hours'] <= (q + 1) * STATUS_HOURS)]
        means.append(float(w['db'].mean()) if len(w) else np.nan)
    rec['db_12h'] = means
    if len(s72) and (s72['db'] >= GREEN).any():
        rec.update(final_db=float(np.nanmean(means)) if np.isfinite(means).any() else np.nan, final_status='green',
                   status_time_h=float(s72.loc[s72['db'] >= GREEN, 'hours'].iloc[0]), status_mode='early green')
    elif np.isfinite(means).any():
        rec.update(final_db=float(np.nanmean(means)), final_status=colour(float(np.nanmean(means))),
                   status_time_h=FINAL_HOURS, status_mode='72 h')
    elif len(series) and series['db'].notna().any():
        first = series[series['db'].notna()].iloc[0]
        late = series[(series['db'].notna()) & (series['hours'] <= first['hours'] + FINAL_HOURS)]
        rec.update(final_db=float(late['db'].mean()), final_status=colour(float(late['db'].mean())),
                   status_time_h=float(first['hours']), status_mode='late')
    else:
        rec.update(final_db=np.nan, final_status='none', status_time_h=np.nan, status_mode='no estimate')
    # ---- outcome: larger event within 3 months and the G-K radius
    r_km = float(gardner_knopoff_km(M))
    after = cat[(cat['time'] > t_ev + pd.Timedelta(minutes=5)) & (cat['time'] <= t_ev + pd.Timedelta(days=VALID_DAYS)) &
                (cat['mag'] > M)]
    if len(after):
        d = haversine_km(lon0, lat0, after['lon'].to_numpy(), after['lat'].to_numpy())
        after = after[d <= r_km]
    rec.update(gk_radius_km=r_km, larger_followed=bool(len(after)),
               larger_event=f"{after['time'].iloc[0]:%Y-%m-%d} M{after['mag'].iloc[0]:.1f}" if len(after) else '')
    st, fol = rec['final_status'], rec['larger_followed']
    rec['class'] = ('true positive' if st == 'red' and fol else 'false alarm' if st == 'red' else
                    'missed' if st == 'green' and fol else 'correct negative' if st == 'green' else
                    'neutral' if st in ('yellow', 'orange') else 'no status')
    return rec, series


def run_ftls(min_mag=6.0, post_volume='box'):
    """FTLS for every M >= min_mag since T0 in the six regions (duplicate catalogue entries within 5 min and 50 km
    of a previous trigger are skipped). Returns (events DataFrame, dict of series)."""
    cat = load_cenc(); gcmt = load_gcmt_planes()
    recs, series = [], {}
    for region in REGIONS:
        ev = targets(region, min_mag, t0=T0, t1=pd.Timestamp('2023-08-14'))
        last = None
        for _, e in ev.iterrows():
            if last is not None and (e['time'] - last['time']).total_seconds() < 300 and \
                    haversine_km(e['lon'], e['lat'], last['lon'], last['lat']) < 50:
                continue
            last = e
            g = match_gcmt(e, gcmt)
            if g is None:
                recs.append(dict(region=region, time=e['time'], mag=e['mag'], lon=e['lon'], lat=e['lat'], final_status='none',
                                 class_='no GCMT')); continue
            rec, s = ftls_event(e, cat, g, post_volume=post_volume)
            rec['in_sequence'] = bool(e['in_sequence'])
            recs.append(rec); series[(region, e['time'])] = s
    return pd.DataFrame(recs), series


def ftls_summary(events):
    """Confusion counts per region and total."""
    rows = []
    for region in REGIONS + ['ALL']:
        e = events if region == 'ALL' else events[events.region == region]
        c = e['class'].value_counts()
        tp, fp, fn, tn = (int(c.get(k, 0)) for k in ('true positive', 'false alarm', 'missed', 'correct negative'))
        neu, non = int(c.get('neutral', 0)), int(c.get('no status', 0) + c.get('no GCMT', 0))
        pd_ = e['pre_drop'].dropna()
        rows.append(dict(region=region, n_events=len(e), with_status=tp + fp + fn + tn + neu, true_positive=tp,
                         false_alarm=fp, missed=fn, correct_negative=tn, neutral=neu, no_status=non,
                         accuracy=(tp + tn) / (tp + fp + fn + tn) if (tp + fp + fn + tn) else np.nan,
                         false_alarm_rate=fp / (tp + fp) if (tp + fp) else np.nan,
                         pre_drop_k=int(pd_.astype(bool).sum()), pre_drop_n=len(pd_),
                         larger_followed=int(e['larger_followed'].fillna(False).astype(bool).sum())))
    return pd.DataFrame(rows)
