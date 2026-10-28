"""Helpers for claims_test.ipynb: percentage test of four literature b-value claims per sub-region.

Constants and the functions `marked_events`, `marked_cutoff` and `spatial_b_value_map_2D` are copied
verbatim from bvalue4longmenshan.ipynb (cells 2 and 18) so the paper's settings are reused unchanged.
The sliding-window series `b_series` replicates seismostats.plots.plot_b_series_constant_nm
(seismostats 1.0.0rc2, plots/statistical.py), which returns only axes.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from seismostats.analysis import (BPositiveBValueEstimator, ClassicBValueEstimator,
                                  estimate_mc_maxc)
from seismostats.utils import bin_to_precision

# --------------------------------------------------------------------------- paths
EXP = Path(__file__).resolve().parents[1]                     # b_value/experiments
CLAIMS = EXP / 'claims'
DATA = CLAIMS / 'data'
OUT = CLAIMS / 'out'
HEATFLOW_XLSX = Path.home() / '🎓Academic/Research/🪺Dataset/热流数据/heatflow5_china.xlsx'
FAULTS_CFM = EXP.parent / 'plotting/faults/202409-3D_CFM_lonlat.txt'

REGION_CSV = {
    'LMS': EXP / '1chuandian/1.1龙门山/longmenshan.csv',
    'XSH': EXP / '1chuandian/1.2鲜水河/xianshuihe.csv',
    'WYN': EXP / '1chuandian/1.3左下/zuoxia.csv',
    'XJ': EXP / '1chuandian/1.4小江/xiaojiang.csv',
    'NCP': EXP / '2huabei/region2_huabei.csv',
    'NW': EXP / '3tibet/region3_tibet.csv',
}
REGIONS = list(REGION_CSV)
REGION_NAME = {'LMS': 'Longmenshan', 'XSH': 'Xianshuihe', 'WYN': 'Western Yunnan', 'XJ': 'Xiaojiang',
               'NCP': 'North China Plain', 'NW': 'Northwest China'}
# sub-region boxes (lon_min, lon_max, lat_min, lat_max), the extent of the regional catalogues
BOXES = {'LMS': (102.5, 106.0, 30.0, 33.5), 'XSH': (101.4, 102.75, 28.35, 30.6),
         'WYN': (97.0, 101.5, 24.0, 28.2), 'XJ': (101.8, 104.4, 26.0, 27.7),
         'NCP': (110.0, 117.0, 34.0, 42.0), 'NW': (73.0, 90.0, 35.0, 45.0)}

# ------------------------------------------------------------- verbatim from cell 2 (LMS notebook)
# N of the sliding window (panel d and b-significant), Table S1; not named N_M, which the MCSENS/MCCOMP cells rebind
WINDOW_N = {'LMS_2006-2010': 500, 'LMS_2011-2016': 500, 'LMS_2016-2023': 100, 'LMS_1970-2023': 1000,
            'XSH_2022-2023': 300, 'XSH_1970-2023': 800, 'XSH_2008-2023': 800,
            'WYN_2020-2023': 500, 'WYN_2000-2023': 1000, 'WYN_1970-2023': 1000,
            'XJ_2016-2023': 150, 'XJ_1970-2023': 1500,
            'NCP_2008-2023': 500, 'NCP_1970-2023': 1200,
            'NW_2006-2023': 3000, 'NW_1970-2023': 8000}
# main-text windows (Figs. 4-10): highlight threshold M >= 6; all other windows (SI long-term figures): M >= 7
MAIN_WINDOWS = {'LMS_2006-2010', 'LMS_2011-2016', 'XSH_2022-2023', 'WYN_2020-2023', 'XJ_2016-2023',
                'NCP_2008-2023', 'NW_2006-2023'}
N_LARGEST = 5              # panels (b), (d), (f) mark every event with M >= the fifth-largest magnitude
C_GR, C_BPOS = '#1f77b4', '#d62728'                   # panel (a): magnitudes (classic-b blue) and positive differences (b-positive red)
CLASS_COLOURS = [(7.0, '#d62728', 'M7'), (6.0, '#ff7f0e', 'M6'), (5.0, '#2ca02c', 'M5')]   # panel (b) circles
# events per grid node of the spatial b-map (panel e), per region notebook
N_PER_GRID = {'LMS': 150, 'XSH': 200, 'WYN': 300, 'XJ': 150, 'NCP': 300, 'NW': 500}
GRID_SIZE = 0.01
MAX_RADIUS = 0.5
DELTA_M = 0.1
# main-text analysis windows per region (start, end) with end exclusive, as in cell 2 ANALYSIS_WINDOWS
MAIN_WINDOW_SPAN = {'LMS_2006-2010': ('2006-01-01', '2011-01-01'), 'LMS_2011-2016': ('2011-01-01', '2017-01-01'),
                    'XSH_2022-2023': ('2022-01-01', '2024-01-01'), 'WYN_2020-2023': ('2020-01-01', '2024-01-01'),
                    'XJ_2016-2023': ('2016-01-01', '2024-01-01'), 'NCP_2008-2023': ('2008-01-01', '2024-01-01'),
                    'NW_2006-2023': ('2006-01-01', '2024-01-01')}


def marked_events(mags, times, thr, n_largest=N_LARGEST):
    """Marking rule. Returns (numbers, thr_set, top_set):
    numbers  dict {positional index: number}, numbering in time order over the union of the two sets;
    thr_set  indices of events with M >= thr (circles in panel b), or top_set when there is none;
    top_set  indices of the n_largest largest events (dashed lines in panel d, stars in panel f).
    Indices refer to `mags`/`times` in whatever order they are passed."""
    mags = np.asarray(mags, float)
    times = np.asarray(times, dtype='datetime64[ns]')
    m_r = np.round(mags, 1)
    cutoff = np.sort(m_r)[-n_largest] if len(m_r) >= n_largest else m_r.min()
    top_set = np.where(m_r >= cutoff)[0]                # five largest, ties at the fifth magnitude included
    thr_set = np.where(mags >= thr)[0]
    if len(thr_set) == 0:
        thr_set = top_set
    ordered = top_set[np.argsort(times[top_set], kind='stable')]   # numbered in time order (panels b, d, f)
    numbers = {int(i): k for k, i in enumerate(ordered, 1)}
    numbers.update({int(i): None for i in thr_set if int(i) not in numbers})
    by_time = lambda s: np.asarray(sorted(s, key=lambda i: times[int(i)]), dtype=int)
    return numbers, by_time(thr_set), by_time(top_set)


def marked_cutoff(mags, n_largest=N_LARGEST):
    """Magnitude of the fifth-largest event (rounded to 0.1): every event at or above it is marked."""
    m_r = np.round(np.asarray(mags, float), 1)
    return float(np.sort(m_r)[-n_largest]) if len(m_r) >= n_largest else float(m_r.min())


# ------------------------------------------------------------ verbatim from cell 18 (LMS notebook)
def spatial_b_value_map_2D(
    lon,
    lat,
    mags,
    grid_size=0.1,
    n_events=100,      # ZMAP逻辑的核心：目标样本数 N
    max_radius=None,   # 可选：最大搜索半径 (度)，防止连接到几百公里外的地震
    mc=None,
    delta_m=0.1,
    b_method=None,
    times=None,        # added: needed by BPositiveBValueEstimator (events sorted by time inside each node)
    verbose=True,
):
    """
    【ZMAP 模式】基于最近邻 (Nearest Neighbor) 的移动窗口扫描。
    逻辑：每个网格点强制取最近的 n_events 个地震进行计算。
    """
    # 1. 数据预清洗 (转为 numpy 数组)
    lons_arr = np.array(lon)
    lats_arr = np.array(lat)
    mags_arr = np.array(mags)
    times_arr = None if times is None else np.asarray(times, dtype='datetime64[ns]')

    valid_mask = ~np.isnan(mags_arr)
    lons_arr = lons_arr[valid_mask]
    lats_arr = lats_arr[valid_mask]
    mags_arr = mags_arr[valid_mask]
    if times_arr is not None:
        times_arr = times_arr[valid_mask]

    # 2. Mc 筛选
    # ZMAP 逻辑：只对大于 Mc 的“有效地震”进行最近邻搜索
    if mc is not None:
        idx = mags_arr >= (mc - delta_m / 2)
        lons_proc = lons_arr[idx]
        lats_proc = lats_arr[idx]
        mags_proc = mags_arr[idx]
        times_proc = None if times_arr is None else times_arr[idx]
    else:
        lons_proc = lons_arr
        lats_proc = lats_arr
        mags_proc = mags_arr
        times_proc = times_arr

    total_events = len(mags_proc)
    if verbose:
        print(f"筛选后用于扫描的有效事件库: {total_events} 个")

    if total_events < n_events:
        print(f"错误：总有效地震数 ({total_events}) 小于目标采样数 ({n_events})，无法计算。")
        return None, None, None, None, None

    # 3. 定义网格 (Grid Nodes)
    x_min, x_max = lons_proc.min(), lons_proc.max()
    y_min, y_max = lats_proc.min(), lats_proc.max()

    xi = np.arange(x_min, x_max, grid_size)
    yi = np.arange(y_min, y_max, grid_size)
    X, Y = np.meshgrid(xi, yi)

    # 初始化结果矩阵
    B_grid = np.full(X.shape, np.nan)
    Std_grid = np.full(X.shape, np.nan)
    Radius_grid = np.full(X.shape, np.nan) # 新增：记录每个点的搜索半径

    # 初始化估算器
    if b_method is None:
        raise ValueError("请传入 b_method 参数 (例如 ClassicBValueEstimator)")
    else:
        estimator = b_method()
    positive = b_method is BPositiveBValueEstimator
    if positive and times_proc is None:
        raise ValueError("BPositiveBValueEstimator needs `times`.")

    # 4. 空间扫描 (Moving Window)
    if verbose:
        print(f"开始 ZMAP 模式扫描... 网格: {X.shape}, 每个点取最近 {n_events} 个事件")

    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            node_lon = X[i, j]
            node_lat = Y[i, j]

            # 计算到所有有效地震的距离平方 (欧氏距离近似)
            dists_sq = (lons_proc - node_lon)**2 + (lats_proc - node_lat)**2

            # --- ZMAP 核心算法 ---
            nearest_idx_unsorted = np.argpartition(dists_sq, n_events)[:n_events]

            current_dists_sq = dists_sq[nearest_idx_unsorted]
            max_dist_sq = np.max(current_dists_sq)
            current_radius_deg = np.sqrt(max_dist_sq)

            local_mags = mags_proc[nearest_idx_unsorted]

            # --- 额外检查：最大半径限制 ---
            if max_radius is not None and current_radius_deg > max_radius:
                B_grid[i, j] = np.nan
                continue

            # 计算 b 值
            if positive:
                local_times = times_proc[nearest_idx_unsorted]
                order = np.argsort(local_times, kind='stable')
                estimator.calculate(local_mags[order], mc=mc, delta_m=delta_m, dmc=delta_m,
                                    times=local_times[order])
            else:
                estimator.calculate(local_mags, mc=mc, delta_m=delta_m)

            B_grid[i, j] = estimator.b_value
            Std_grid[i, j] = estimator.std
            Radius_grid[i, j] = current_radius_deg

    if verbose:
        print("计算完成。")
    return X, Y, B_grid, Std_grid, Radius_grid


# ------------------------------------------------------------------------ loading
def load_region(code, t0=None, t1=None):
    """Regional catalogue as DataFrame(time, lon, lat, depth, mag) sorted by time, magnitudes binned to 0.1.
    Same recipe as cell 4/6 of the regional notebooks: filter time >= t0 and time < t1."""
    df = pd.read_csv(REGION_CSV[code], usecols=['year', 'month', 'day', 'hour', 'minute', 'second',
                                                'latitude', 'longitude', 'depth', 'magnitude'])
    df['time'] = pd.to_datetime(df[['year', 'month', 'day', 'hour', 'minute', 'second']])
    if t0 is not None:
        df = df[df['time'] >= pd.Timestamp(t0)]
    if t1 is not None:
        df = df[df['time'] < pd.Timestamp(t1)]
    df = df.rename(columns={'latitude': 'lat', 'longitude': 'lon', 'magnitude': 'mag'})
    df = df[['time', 'lon', 'lat', 'depth', 'mag']].sort_values('time', kind='stable').reset_index(drop=True)
    df['mag'] = bin_to_precision(df['mag'].to_numpy(float), DELTA_M)
    return df


def load_figrun():
    """All per-window *_figrun.csv rows (Mc, N, whole-window b of both estimators), indexed by window label."""
    rows = []
    for code, p in REGION_CSV.items():
        for f in sorted(p.parent.glob(f'{code}_*_figrun.csv')):
            r = pd.read_csv(f)
            r['region'] = code
            rows.append(r)
    fr = pd.concat(rows, ignore_index=True).set_index('window')
    return fr


def window_span(label):
    """(t0, t1) of a window label: main-text windows from MAIN_WINDOW_SPAN, others from the years in the label."""
    if label in MAIN_WINDOW_SPAN:
        return MAIN_WINDOW_SPAN[label]
    y0, y1 = label.split('_')[1].split('-')
    return f'{y0}-01-01', f'{int(y1) + 1}-01-01'


# --------------------------------------------------------------------- b estimation
def b_window(mags, times, mc, estimator_cls, delta_m=DELTA_M):
    """One estimate on a sample (events must be sorted by time for b-positive). Returns the estimator."""
    est = estimator_cls()
    mags = np.asarray(mags, float)
    if estimator_cls is BPositiveBValueEstimator:
        est.calculate(mags, mc=mc, delta_m=delta_m, dmc=delta_m, times=np.asarray(times, dtype='datetime64[ns]'))
    else:
        est.calculate(mags, mc=mc, delta_m=delta_m)
    return est


def b_series(mags, times, mc, n_m, estimator_cls, delta_m=DELTA_M, step=1, min_num=2):
    """Sliding window of n_m consecutive events with M >= mc - delta_m/2, value at the last event of the window
    (plot_technique='right' of seismostats.plots.plot_b_series_constant_nm).
    Returns DataFrame(i_first, i_last, t_first, t_last, b, std, n); indices refer to the filtered, time-sorted sample."""
    mags = np.asarray(mags, float)
    times = np.asarray(times, dtype='datetime64[ns]')
    order = np.argsort(times, kind='stable')
    mags, times = mags[order], times[order]
    keep = mags >= mc - delta_m / 2
    mags, times = mags[keep], times[keep]
    positive = estimator_cls is BPositiveBValueEstimator
    est = estimator_cls()
    rec = []
    for ii in range(0, len(mags) - n_m + 1, step):
        w_m = mags[ii:ii + n_m]
        if positive:
            est.calculate(w_m, mc=mc, delta_m=delta_m, dmc=delta_m, times=times[ii:ii + n_m])
        else:
            est.calculate(w_m, mc=mc, delta_m=delta_m)
        ok = est.n >= min_num
        rec.append((ii, ii + n_m - 1, times[ii], times[ii + n_m - 1],
                    est.b_value if ok else np.nan, est.std if ok else np.nan, est.n))
    return pd.DataFrame(rec, columns=['i_first', 'i_last', 't_first', 't_last', 'b', 'std', 'n'])


def filtered_sample(df, mc, delta_m=DELTA_M):
    """The time-sorted sample a series is built on (M >= mc - delta_m/2), so series indices can be mapped back."""
    s = df[df['mag'] >= mc - delta_m / 2].sort_values('time', kind='stable').reset_index(drop=True)
    return s


# ----------------------------------------------------------------------- geometry
def haversine_km(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.radians, (lon1, lat1, lon2, lat2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def in_box(lon, lat, box):
    x0, x1, y0, y1 = box
    return (lon >= x0) & (lon <= x1) & (lat >= y0) & (lat <= y1)


def gardner_knopoff_days(m):
    """Aftershock time window of Gardner & Knopoff (1974), in days."""
    m = np.asarray(m, float)
    return np.where(m < 6.5, 10 ** (0.5409 * m - 0.547), 10 ** (0.032 * m + 2.7389))


def gardner_knopoff_km(m):
    """Aftershock distance window of Gardner & Knopoff (1974), in km."""
    return 10 ** (0.1238 * np.asarray(m, float) + 0.983)


# ================================================================== C4: temporal claims
# Units = the seven main-text windows (Figs 4-10), each with its own N and Mc (Table S1); reference = the window's
# whole-period b; targets = the events marked in panel (d) (fifth-largest rule); alarm duration 90 days.
FIG_OF = {'LMS_2006-2010': 'Fig. 4', 'LMS_2011-2016': 'Fig. 5', 'XSH_2022-2023': 'Fig. 6', 'WYN_2020-2023': 'Fig. 7',
          'XJ_2016-2023': 'Fig. 8', 'NCP_2008-2023': 'Fig. 9', 'NW_2006-2023': 'Fig. 10'}
MAIN_LABELS = list(FIG_OF)
T_START, T_END = pd.Timestamp('2008-01-01'), pd.Timestamp('2023-08-14')   # catalogue ends 2023-08-13 (SI epicentre row)
DROP = 0.10           # Gulia & Wiemer (2019): b at least 10% below the reference
RISE = 0.10           # b at least 10% above the pre-event level
POST_BLIND_DAYS = 7   # post-event window starts after this blind time (short-term aftershock incompleteness)
ALARM_DAYS = 7        # a drop is a hit when an M >= ALARM_TARGET_MAG follows within this many days (decided 2026-10-26)
MIN_REF = 5           # windows needed before a reference (expanding mean) exists
ALARM_TARGET_MAG = 4.0  # an alarm episode is a hit when an M >= this occurs in the window's region within ALARM_DAYS
EDGE_DEG = 0.25

_SERIES_CACHE = {}


def series_for_window(label, estimator_cls, figrun):
    """Sliding-window series of a paper window (its span, N and Mc), cached. Returns (series, sample)."""
    key = (label, estimator_cls.__name__)
    if key not in _SERIES_CACHE:
        region = figrun.loc[label, 'region']
        t0, t1 = window_span(label)
        df = load_region(region, t0, t1)
        mc, n_m = float(figrun.loc[label, 'mc']), int(figrun.loc[label, 'n_m'])
        s = b_series(df['mag'], df['time'], mc, n_m, estimator_cls)
        s['label'], s['mc'], s['n_m'] = label, mc, n_m
        _SERIES_CACHE[key] = (s, filtered_sample(df, mc), df)
    return _SERIES_CACHE[key][:2]


def window_catalogue(label, figrun):
    """Unfiltered catalogue of a paper window (cached by series_for_window)."""
    for key, v in _SERIES_CACHE.items():
        if key[0] == label:
            return v[2]
    series_for_window(label, ClassicBValueEstimator, figrun)
    return window_catalogue(label, figrun)


def flag_in_sequence(ev, cat, min_mag=5.0, cloud_days=7, cloud_mag=4.0):
    """True when ev lies inside the aftershock sequence of a larger earlier event: within that event's Gardner-Knopoff
    time window, and within its G-K distance of the mainshock or of any of its first-week M >= cloud_mag aftershocks
    (the early aftershock cloud stands in for the rupture extent, which the epicentre distance misses for an M7.9)."""
    big = cat[(cat['mag'] >= min_mag) & (cat['mag'] > ev['mag']) & (cat['time'] < ev['time'])]
    if big.empty:
        return False
    dt_days = (ev['time'] - big['time']).dt.total_seconds().to_numpy() / 86400
    big = big[dt_days <= gardner_knopoff_days(big['mag'])]
    for _, m in big.iterrows():
        cloud = cat[(cat['time'] >= m['time']) & (cat['time'] <= m['time'] + pd.Timedelta(days=cloud_days)) &
                    (cat['mag'] >= cloud_mag) & (cat['time'] < ev['time'])]
        d = haversine_km(ev['lon'], ev['lat'], cloud['lon'].to_numpy(), cloud['lat'].to_numpy())
        if np.any(d <= gardner_knopoff_km(m['mag'])):
            return True
    return False


def targets(region, min_mag=6.0, t0=T_START, t1=T_END, cat=None):
    """M >= min_mag events of a region in [t0, t1) with edge and in-sequence flags (used by the SI epicentre row)."""
    full = load_region(region) if cat is None else cat
    ev = full[(full['mag'] >= min_mag) & (full['time'] >= t0) & (full['time'] < t1)].copy()
    x0, x1, y0, y1 = BOXES[region]
    ev['near_edge'] = (np.minimum(ev['lon'] - x0, x1 - ev['lon']) < EDGE_DEG) | \
                      (np.minimum(ev['lat'] - y0, y1 - ev['lat']) < EDGE_DEG)
    ev['in_sequence'] = [flag_in_sequence(r, full) for _, r in ev.iterrows()]
    ev['region'] = region
    return ev.reset_index(drop=True)


def window_c4(label, estimator_cls, figrun, alarm_days=ALARM_DAYS):
    """C4 for one main-text window: marked events (drop before / rise after), alarm episodes (false alarms) and
    the alarm fraction, relative to the pre-event expanding-mean reference. Returns (events, episodes, summary)."""
    region = figrun.loc[label, 'region']
    est = 'classical' if estimator_cls is ClassicBValueEstimator else 'positive'
    b_ww = float(figrun.loc[label, 'b_classical' if est == 'classical' else 'b_positive'])   # Table S1 value (figure line)
    s, _ = series_for_window(label, estimator_cls, figrun)
    df = window_catalogue(label, figrun)
    full = load_region(region)
    s = s.copy()
    # reference of each window = mean of all earlier windows of the panel that share no event with it (pre-event
    # expanding mean: the Gulia & Wiemer 'years preceding' rule applied inside the panel); needs >= MIN_REF windows
    t_l = s['t_last'].to_numpy(dtype='datetime64[ns]'); t_f = s['t_first'].to_numpy(dtype='datetime64[ns]')
    bb = s['b'].to_numpy(); hi = np.searchsorted(t_l, t_f, side='left')
    csum = np.concatenate([[0.0], np.nancumsum(bb)]); cnt = np.concatenate([[0], np.cumsum(np.isfinite(bb))])
    s['n_ref'] = cnt[hi]
    s['ref'] = np.where(cnt[hi] >= MIN_REF, csum[hi] / np.maximum(cnt[hi], 1), np.nan)
    s['rel'] = (s['b'] - s['ref']) / s['ref']
    s['alarm'] = s['rel'] <= -DROP
    t_last = s['t_last'].to_numpy(dtype='datetime64[ns]'); t_first = s['t_first'].to_numpy(dtype='datetime64[ns]')
    b = s['b'].to_numpy(); ref = s['ref'].to_numpy()
    # ---- marked events of panel (d): every event with M >= the fifth-largest magnitude
    numbers, _, top = marked_events(df['mag'].to_numpy(), df['time'].to_numpy(), 6.0)
    ev_rows = []
    for i in top:
        r = df.iloc[i]
        t_ev = np.datetime64(r['time'], 'ns')
        rec = dict(window=label, fig=FIG_OF[label], region=region, estimator=est, number=numbers[int(i)],
                   time=r['time'], mag=r['mag'], lon=r['lon'], lat=r['lat'],
                   initiator=not flag_in_sequence(r, full), b_ww=b_ww, ref=np.nan,
                   b_pre=np.nan, pre_days=np.nan, drop_rel=np.nan, drop=None,
                   b_post=np.nan, post_start_days=np.nan, post_days=np.nan, rise_rel=np.nan, rise=None)
        k = np.searchsorted(t_last, t_ev, side='left') - 1          # last window ending before the event
        if k >= 0 and np.isfinite(b[k]):
            rec.update(b_pre=b[k], pre_days=(t_last[k] - t_first[k]) / np.timedelta64(1, 'D'), ref=ref[k])
            if np.isfinite(ref[k]):
                rec.update(drop_rel=(ref[k] - b[k]) / ref[k], drop=bool((ref[k] - b[k]) / ref[k] >= DROP))
        if True:
            j = np.searchsorted(t_first, t_ev + np.timedelta64(POST_BLIND_DAYS, 'D'), side='right')
            if j < len(s) and np.isfinite(b[j]) and np.isfinite(rec['b_pre']):
                rec.update(b_post=b[j], post_start_days=(t_first[j] - t_ev) / np.timedelta64(1, 'D'),
                           post_days=(t_last[j] - t_first[j]) / np.timedelta64(1, 'D'),
                           rise_rel=(b[j] - rec['b_pre']) / rec['b_pre'],
                           rise=bool((b[j] - rec['b_pre']) / rec['b_pre'] >= RISE))
        ev_rows.append(rec)
    events = pd.DataFrame(ev_rows)
    mk_t = events['time'].to_numpy(dtype='datetime64[ns]') if len(events) else np.array([], dtype='datetime64[ns]')
    mk_m = events['mag'].to_numpy(float) if len(events) else np.array([])
    tg_t = df.loc[df['mag'] >= ALARM_TARGET_MAG, 'time'].to_numpy(dtype='datetime64[ns]')   # alarm targets: M >= 4 in the window
    # ---- alarm episodes: runs of consecutive alarm windows
    alarm = s['alarm'].to_numpy(bool) & np.isfinite(b)
    edges = np.flatnonzero(np.diff(np.concatenate([[0], alarm.astype(int), [0]])))
    dur = np.timedelta64(alarm_days, 'D')
    ep_rows = []
    for i0, i1 in zip(edges[::2], edges[1::2] - 1):
        start, stop = t_last[i0], t_last[i1]
        hits = tg_t[(tg_t > start) & (tg_t <= start + dur)]
        in_seq = bool(np.any((mk_t < start) & (start <= mk_t + gardner_knopoff_days(mk_m).astype('timedelta64[D]'))))
        ep_rows.append(dict(window=label, fig=FIG_OF[label], region=region, estimator=est, t_start=start, t_last_alarm=stop,
                            n_windows=int(i1 - i0 + 1), min_rel=float(s['rel'].iloc[i0:i1 + 1].min()),
                            hit=len(hits) > 0, first_target=hits[0] if len(hits) else pd.NaT, in_sequence=in_seq))
    episodes = pd.DataFrame(ep_rows)
    # ---- alarm fraction of the window span
    t0 = np.datetime64(pd.Timestamp(window_span(label)[0]), 'ns'); t1 = t_last[-1]
    span = (t1 - t0) / np.timedelta64(1, 'D'); covered = 0.0
    if len(episodes):
        iv = sorted((a, min(a + dur, t1)) for a in episodes['t_start'].to_numpy())
        ca, cb = iv[0]
        for a, z in iv[1:]:
            if a <= cb:
                cb = max(cb, z)
            else:
                covered += (cb - ca) / np.timedelta64(1, 'D'); ca, cb = a, z
        covered += (cb - ca) / np.timedelta64(1, 'D')
    clean = episodes[~episodes['in_sequence']] if len(episodes) else episodes
    days = np.arange(t0, t1, np.timedelta64(1, 'D'))
    ref_share = float(np.mean([np.any((tg_t > d) & (tg_t <= d + dur)) for d in days])) if len(days) else np.nan
    summary = dict(window=label, fig=FIG_OF[label], region=region, estimator=est, n_m=int(s['n_m'].iloc[0]),
                   mc=float(s['mc'].iloc[0]), b_ww=b_ww, n_windows=len(s), n_alarm_windows=int(alarm.sum()),
                   n_marked=len(events),
                   drop_k=int(events['drop'].dropna().astype(bool).sum()) if len(events) else 0,
                   drop_n=int(events['drop'].notna().sum()) if len(events) else 0,
                   rise_k=int(events['rise'].dropna().astype(bool).sum()) if len(events) else 0,
                   rise_n=int(events['rise'].notna().sum()) if len(events) else 0,
                   n_episodes=len(episodes), n_episodes_clean=len(clean),
                   n_false_clean=int((~clean['hit']).sum()) if len(clean) else 0,
                   n_false_all=int((~episodes['hit']).sum()) if len(episodes) else 0,
                   n_hits=int(episodes['hit'].sum()) if len(episodes) else 0, hit_ref_share=ref_share,
                   alarm_fraction=float(alarm.mean()) if len(alarm) else np.nan,   # share of sliding windows under alarm
                   alarm_time_fraction=covered / span)
    return events, episodes, summary


# ================================================================ spatial claims (C1-C3)
from scipy.spatial import cKDTree, Delaunay

MAP_T0, MAP_T1 = '2008-01-01', '2024-01-01'
MAP_GRID = 0.05        # node spacing of the claims maps (display resolution only; N and radius are the paper's)
STRAIN_SIGMA_KM = 60.0
STRAIN_MAX_KM = 150.0
STRAIN_MIN_STA = 6


def bmap_nodes(df, mc, n_events, estimator_cls, box, grid=MAP_GRID, max_radius=MAX_RADIUS, delta_m=DELTA_M):
    """Nearest-n_events b-value at the nodes of a regular grid over `box` (same rule as spatial_b_value_map_2D:
    Euclidean distance in degrees, node invalid when the n-th neighbour is farther than max_radius).
    Catalogue coordinates are given to 0.01 deg, so several events often tie at the n-th distance; all events at
    that distance are included (deterministic; the verbatim function picks tie members by array order).
    Returns DataFrame(lon, lat, b, std, n, radius)."""
    d = df[df['mag'] >= mc - delta_m / 2]
    x0, x1, y0, y1 = box
    xi = np.arange(x0 + grid / 2, x1, grid)
    yi = np.arange(y0 + grid / 2, y1, grid)
    X, Y = np.meshgrid(xi, yi)
    nodes = np.column_stack([X.ravel(), Y.ravel()])
    tree = cKDTree(np.column_stack([d['lon'].to_numpy(), d['lat'].to_numpy()]))
    dist, idx = tree.query(nodes, k=n_events)
    radius = dist[:, -1]
    ok = radius <= max_radius
    mags = d['mag'].to_numpy(); times = d['time'].to_numpy(dtype='datetime64[ns]')
    b = np.full(len(nodes), np.nan); sd = np.full(len(nodes), np.nan); n = np.zeros(len(nodes), int)
    est = estimator_cls()
    positive = estimator_cls is BPositiveBValueEstimator
    for k in np.flatnonzero(ok):
        ii = np.asarray(tree.query_ball_point(nodes[k], radius[k] + 1e-9))
        if positive:
            o = np.argsort(times[ii], kind='stable')
            est.calculate(mags[ii][o], mc=mc, delta_m=delta_m, dmc=delta_m, times=times[ii][o])
        else:
            est.calculate(mags[ii], mc=mc, delta_m=delta_m)
        b[k], sd[k], n[k] = est.b_value, est.std, est.n
    return pd.DataFrame(dict(lon=nodes[:, 0], lat=nodes[:, 1], b=b, std=sd, n=n, radius=radius))


def load_gnss_ws():
    """Wang & Shen (2020) Table S4, Eurasia-fixed velocities. Columns lon, lat, ve, vn, se, sn, station."""
    rows = []
    for line in open(DATA / 'gnss/WangShen2020_TableS4_velocities.txt', errors='ignore'):
        p = line.split()
        if len(p) < 10 or line.startswith(('*', 'Table')):
            continue
        try:
            rows.append((float(p[1]), float(p[2]), float(p[4]), float(p[7]), float(p[5]), float(p[8]), p[0]))
        except ValueError:
            continue
    return pd.DataFrame(rows, columns=['lon', 'lat', 've', 'vn', 'se', 'sn', 'station'])


def load_gnss_ds01():
    """Sichuan-Yunnan field, South China Block-fixed (ProjectB ds01). Same columns."""
    g = pd.read_csv(DATA / 'gnss/ds01_gnss_velocity_field.txt', sep=r'\s+', comment='#', header=None,
                    names=['lon', 'lat', 've', 'vn', 'se', 'sn', 'corr', 'station'])
    return g[['lon', 'lat', 've', 'vn', 'se', 'sn', 'station']]


def strain_at(nodes, sta, sigma_km=STRAIN_SIGMA_KM, max_km=STRAIN_MAX_KM, min_sta=STRAIN_MIN_STA):
    """Strain-rate tensor at each node by weighted least squares v_i = v0 + L (x_i - x0), Gaussian weights
    exp(-d^2/2 sigma^2)/sigma_v^2 over stations within max_km; NaN with fewer than min_sta stations or outside
    the station hull. Returns DataFrame(exx, eyy, exy, e2inv, dilat, n_sta) in nanostrain/yr."""
    hull = Delaunay(np.column_stack([sta['lon'], sta['lat']]))
    inside = hull.find_simplex(np.column_stack([nodes['lon'], nodes['lat']])) >= 0
    out = np.full((len(nodes), 6), np.nan)
    slon, slat = sta['lon'].to_numpy(), sta['lat'].to_numpy()
    ve, vn = sta['ve'].to_numpy(), sta['vn'].to_numpy()
    sv2 = 0.5 * (sta['se'].to_numpy() ** 2 + sta['sn'].to_numpy() ** 2) + 0.01   # mm^2/yr^2, floor 0.1 mm/yr
    for k, (lon0, lat0) in enumerate(zip(nodes['lon'].to_numpy(), nodes['lat'].to_numpy())):
        if not inside[k]:
            continue
        dx = (slon - lon0) * 111.32 * np.cos(np.radians(lat0))
        dy = (slat - lat0) * 110.57
        d = np.hypot(dx, dy)
        sel = d <= max_km
        if sel.sum() < min_sta:
            continue
        w = np.exp(-d[sel] ** 2 / (2 * sigma_km ** 2)) / sv2[sel]
        A = np.column_stack([np.ones(sel.sum()), dx[sel], dy[sel]])
        sw = np.sqrt(w)
        try:
            pe = np.linalg.lstsq(A * sw[:, None], ve[sel] * sw, rcond=None)[0]
            pn = np.linalg.lstsq(A * sw[:, None], vn[sel] * sw, rcond=None)[0]
        except np.linalg.LinAlgError:
            continue
        exx, exy, eyx, eyy = pe[1], pe[2], pn[1], pn[2]          # mm/yr per km = 1e-6 /yr
        exy = 0.5 * (exy + eyx)
        out[k] = (exx * 1e3, eyy * 1e3, exy * 1e3, np.sqrt(exx ** 2 + eyy ** 2 + 2 * exy ** 2) * 1e3,
                  (exx + eyy) * 1e3, sel.sum())
    return pd.DataFrame(out, columns=['exx', 'eyy', 'exy', 'e2inv', 'dilat', 'n_sta'], index=nodes.index)


def parse_gcmt(files=None):
    """GCMT ndk files -> DataFrame(time, lon, lat, depth, mw, rake1, rake2, plunge_t/b/p, style); style by Frohlich (1992)."""
    import math
    files = sorted((DATA / 'gcmt').glob('*.ndk')) if files is None else files
    rows = []
    for f in files:
        L = open(f, errors='ignore').read().splitlines()
        for i in range(0, len(L) - 4, 5):
            l1, l4, l5 = L[i], L[i + 3], L[i + 4]
            try:
                date, hms = l1[5:15], l1[16:26]
                lat, lon, dep = float(l1[27:33]), float(l1[34:41]), float(l1[42:47])
                p = l5.split(); ex = int(l4[0:2]); m0 = float(p[10]) * 10 ** ex
                mw = (2 / 3) * (math.log10(m0) - 16.1)
                r1, r2 = float(p[13]), float(p[16])
                pl_t, pl_b, pl_p = float(p[2]), float(p[5]), float(p[8])   # eigenvector plunges: T, B (null), P
            except (ValueError, IndexError):
                continue
            t = pd.Timestamp(date.replace('/', '-') + ' ' + hms.replace('60.0', '59.9'), tz=None)
            rows.append((t, lon, lat, dep, mw, r1, r2, pl_t, pl_b, pl_p))
    g = pd.DataFrame(rows, columns=['time', 'lon', 'lat', 'depth', 'mw', 'rake1', 'rake2', 'plunge_t', 'plunge_b', 'plunge_p'])
    # Frohlich (1992) classification from the P/T/B-axis plunges
    g['style'] = np.select([g['plunge_t'] >= 50, g['plunge_p'] >= 60, g['plunge_b'] >= 60],
                           ['thrust', 'normal', 'strike-slip'], default='oblique')
    return g


def nearest_node_value(nodes, lon, lat, col='b', max_deg=0.1):
    """Value of `col` at the nearest valid node within max_deg of each (lon, lat); NaN otherwise."""
    v = nodes[np.isfinite(nodes[col])]
    tree = cKDTree(np.column_stack([v['lon'], v['lat']]))
    d, i = tree.query(np.column_stack([np.atleast_1d(lon), np.atleast_1d(lat)]))
    val = v[col].to_numpy()[i].astype(float)
    val[d > max_deg] = np.nan
    return val


# ================================================================ C3 covariates and SI rows
HF_RADIUS_KM = 50.0
HF_MIN_SITES = 2


def load_heatflow():
    """China heat-flow compilation (heatflow5_china.xlsx): DataFrame(lon, lat, hf) in mW/m^2."""
    hf = pd.read_excel(HEATFLOW_XLSX)
    hf.columns = ['id', 'lon', 'lat', 'hf']
    return hf.dropna(subset=['lon', 'lat', 'hf']).astype({'lon': float, 'lat': float, 'hf': float})


def heatflow_at(nodes, hf, radius_km=HF_RADIUS_KM, min_sites=HF_MIN_SITES):
    """Inverse-distance (power 2) mean of heat-flow sites within radius_km of each node; NaN below min_sites."""
    out = np.full(len(nodes), np.nan); cnt = np.zeros(len(nodes), int)
    hl, hb, hv = hf['lon'].to_numpy(), hf['lat'].to_numpy(), hf['hf'].to_numpy()
    for k, (lon0, lat0) in enumerate(zip(nodes['lon'].to_numpy(), nodes['lat'].to_numpy())):
        d = haversine_km(lon0, lat0, hl, hb)
        sel = d <= radius_km
        cnt[k] = sel.sum()
        if cnt[k] >= min_sites:
            w = 1.0 / np.maximum(d[sel], 1.0) ** 2
            out[k] = np.sum(w * hv[sel]) / np.sum(w)
    return pd.DataFrame(dict(hf=out, n_hf=cnt), index=nodes.index)




def epicentre_rows(region, events, figrun, years=10.0):
    """Pre-event b at each epicentre: map over [t_ev - years, t_ev) with Mc = MAXC+0.2 of that catalogue and the
    paper's N_PER_GRID; compared with the whole-window b of the same pre-event catalogue (both estimators)."""
    full = load_region(region)
    rows = []
    for _, ev in events.iterrows():
        t1 = ev['time']; t0 = max(t1 - pd.Timedelta(days=365.25 * years), pd.Timestamp('2008-01-01'))
        df = full[(full['time'] >= t0) & (full['time'] < t1)]
        if len(df) < 5 * N_PER_GRID[region]:
            rows.append(dict(region=region, time=t1, mag=ev['mag'], note='too few events')); continue
        mc = float(estimate_mc_maxc(df['mag'].to_numpy(), fmd_bin=0.1)[0])
        r = dict(region=region, time=t1, mag=ev['mag'], lon=ev['lon'], lat=ev['lat'], mc=mc, n=len(df), t0=t0)
        for est, cls in (('c', ClassicBValueEstimator), ('p', BPositiveBValueEstimator)):
            b0 = b_window(df['mag'], df['time'], mc, cls).b_value
            nodes = bmap_nodes(df, mc, N_PER_GRID[region], cls, BOXES[region])
            b_epi = nearest_node_value(nodes, ev['lon'], ev['lat'], col='b')[0]
            valid = nodes['b'].dropna()
            r[f'b0_{est}'] = b0; r[f'b_epi_{est}'] = b_epi
            r[f'low_{est}'] = bool(b_epi < b0) if np.isfinite(b_epi) else None
            r[f'pct_{est}'] = float((valid < b_epi).mean()) if np.isfinite(b_epi) else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def depth_rows(region, mc, bin_km=5.0, max_km=30.0, min_n=200, fixed_frac=0.05):
    """b per depth bin, post-2008, depth > 0; with and without 'fixed' depth values (holding > fixed_frac of events).
    Returns (bins DataFrame, summary dict with % of adjacent bin pairs where b decreases with depth)."""
    df = load_region(region, MAP_T0, MAP_T1)
    df = df[(df['depth'] > 0) & (df['mag'] >= mc - DELTA_M / 2)]
    vc = df['depth'].value_counts(normalize=True)
    fixed = vc[vc > fixed_frac].index.tolist()
    out = []; summ = dict(region=region, mc=mc, fixed_values=fixed, fixed_share=float(vc[vc > fixed_frac].sum()))
    for variant, d in (('all', df), ('no_fixed', df[~df['depth'].isin(fixed)])):
        edges = np.arange(0, max_km + bin_km, bin_km)
        res = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            s = d[(d['depth'] > lo) & (d['depth'] <= hi)]
            if len(s) < min_n:
                continue
            bc = b_window(s['mag'], s['time'], mc, ClassicBValueEstimator)
            bp = b_window(s['mag'], s['time'], mc, BPositiveBValueEstimator)
            res.append(dict(region=region, variant=variant, z_lo=lo, z_hi=hi, n=len(s), b_c=bc.b_value, std_c=bc.std,
                            b_p=bp.b_value, std_p=bp.std))
        res = pd.DataFrame(res); out.append(res)
        for est in ('c', 'p'):
            dif = np.diff(res[f'b_{est}'].to_numpy()) if len(res) > 1 else np.array([])
            summ[f'{variant}_dec_{est}'] = f'{int((dif < 0).sum())}/{len(dif)}' if len(dif) else 'n/a'
    return pd.concat(out, ignore_index=True), summ
