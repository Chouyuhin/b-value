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

Temporal claims (C4), unit = the seven main-text windows (Figs 4–10), each with its own N and Mc (Table S1); the
1970–2023 SI windows are not used (decided 2026-10-23):
- reference of a window = mean of all earlier windows of the panel that share no event with it (pre-event expanding mean,
  \u2265 5 windows; the Gulia & Wiemer 'years preceding' rule inside the panel; decided 2026-10-26). The Table S1 b is only drawn.
- targets = the events marked in panel (d) (every event with M ≥ the fifth-largest magnitude); the initiator/aftershock
  distinction was removed on 2026-10-26 (the `initiator` flag stays in events_C4.csv for information).
- drop before = last window ending before the event ≥ 10 % below its reference (Gulia & Wiemer 2019 threshold).
- rise after = first N events after a 7-day blind time, ≥ 10 % above the pre-event b.
- alarm = window ≥ 10 % below the reference; contiguous windows = one episode; alarm duration 90 days; a drop (episode) is a
  hit when an M ≥ 4 occurs in the window region within the following 7 days (30 and 90 d examined, 7 d chosen 2026-10-26); false-alarm rate = 100 % − hit rate; all
  episodes are counted (the in-sequence flag stays in alarms_C4.csv); drops are identified on each estimator's own series, so
  the two estimators have different drop counts (kept deliberately, 2026-10-27);
  the 'random date' share (days of the window followed by an M ≥ 4 within 7 d) is printed as the chance reference; no pooled
  total and no alarm-fraction column (removed 2026-10-27).

Reference provenance: the per-window reference follows step 4 of Gulia et al. (2024, SRL 95:3579; PDF in
`data/srl-2024163.1.pdf`): mean of the moving-window b-values over all events before the mainshock, here applied inside
each panel with the paper's own N (not 150-event b-positive windows). Their full traffic-light procedure was tried on
2026-10-26 and dropped (`ftls.py` kept unused; its result: 0 larger events followed any M>=6, 5 reds = false alarms).

Node percentages are area fractions (neighbourhoods overlap); node counts are printed with every percentage.
