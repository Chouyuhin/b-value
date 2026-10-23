# Claims test: four literature b-value claims on the six sub-regions

Percentages of units that support each claim, computed with the paper's own settings. No hypothesis tests.

## Files

- `claims_lib.py` – constants and the paper's helpers (`marked_events`, `spatial_b_value_map_2D`, `WINDOW_N`, ...) copied verbatim from
  `1chuandian/1.1龙门山/bvalue4longmenshan.ipynb`; loaders; sliding-window series (replicates `plot_b_series_constant_nm`); C4 rules;
  b-maps; strain rate; heat flow; GCMT.
- `claims_report.py` – table (`build_table`, `table_markdown`) and figures (`fig_summary`, `figS_C12`, `figS_C3`, `figS_C4`).
- `build_notebook.py` – writes `claims_test.ipynb`; run the notebook with
  `KMP_DUPLICATE_LIB_OK=TRUE jupyter nbconvert --to notebook --execute --inplace claims_test.ipynb` (kernel `gmtenv`).
- `data/` – external inputs (see `data/README.md`). Heat flow is read from `🪺Dataset/热流数据/heatflow5_china.xlsx`.
- `out/` – `claims_table.csv` / `.md`, `nodes_<REG>.csv`, `events_C4.csv`, `alarms_C4.csv`, `alarms_C4_summary.csv`,
  `thresholds_2008-2023.csv`, `gcmt_summary.csv`, `epicentres_SI.csv`, `depth_SI.csv`, `fig_summary.pdf`, `figS_C12.pdf`,
  `figS_C3.pdf`, `figS_C4.pdf`.

## Rules (decided 2026-10-23)

Spatial claims (C1–C3), unit = grid node (0.05° spacing; N per node and 0.5° radius as in the paper; 2008–2023; both estimators):
- low / high b = below / above the whole-window b of the same region, period and estimator (b₀, Mc = MAXC + 0.2);
  the share of low-b (high-b) nodes is the reference for rows that condition on a covariate class.
- strain rate: weighted least squares of Wang & Shen (2020) velocities (Gaussian weight σ = 60 km, stations within 150 km,
  ≥ 6 stations, node inside the station hull); high = above the region's node median; dilatation terciles = contraction /
  shear / extension. The SCB-fixed field (ds01) is the check for the four SW regions.
- GCMT M ≥ 5 (1976–2023) classified by the Frohlich (1992) P/T/B plunges (thrust: T plunge ≥ 50°, normal: P plunge ≥ 60°,
  strike-slip: B plunge ≥ 60°, else oblique);
  rows reported where a region has ≥ 10 mechanisms of the class.
- heat flow at a node = inverse-distance mean of sites within 50 km (≥ 2 sites); high = above the median of the covered nodes.
  Rows need ≥ 20 nodes (XSH has no heat-flow site). Hot springs were dropped from the analysis (2026-10-23); the files stay in
  `data/hotsprings/` unused.

Temporal claims (C4), unit = M ≥ 6 event since 2008 (main) or M ≥ 5 in the main-text windows (secondary):
- series of the paper window containing the event (its N and Mc; shortest window where two overlap, `TILES`).
- pre-event b = last window ending before the event; reference = mean of windows ending in the year before the pre-event window
  starts; drop = (ref − b_pre)/ref ≥ 10 % (Gulia & Wiemer 2019 threshold; no std floor).
- post-event b = first N events after a 7-day blind time (short-term aftershock incompleteness); rise = ≥ 10 % above b_pre.
- alarms: windows ≥ 10 % below the running reference (same rule as for events, no std floor; decided 2026-10-23); contiguous
  windows = one episode; duration 1 yr; false-alarm target = M ≥ 5 in the box (the event rows keep M ≥ 6).
  False-alarm rate = episodes with no M ≥ 5 within 1 yr; alarm fraction = share of 2008–2023 covered by alarms. The row
  "targets preceded by an alarm" was dropped.
- events inside the Gardner–Knopoff time window of a larger earlier event AND within its G-K distance of that event or of its
  first-week M ≥ 4 aftershocks are `in_sequence`; alarm episodes starting inside the G-K time window of an M ≥ 6 in the box are
  `in_sequence` (the series is region-wide). Both are excluded from the clean counts.

Node percentages are area fractions (neighbourhoods overlap); node counts are printed with every percentage.
