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
  `figS_C3.pdf`, `figS_C4.pdf`, `table_C4_events.csv` / `.md` (Form B: the C4 per-event table for the main text),
  `table_main.pdf` (`fig_table`: the two main-text tables, spatial rows (a) + C4 Form B (b), rendered as one figure; 2026-10-30),
  `tables_main.docx` (the same two tables as Word tables, written by `make_word_tables.py` from the CSVs with the base-miniconda
  python, which has python-docx; 2026-10-30).

## Rules (decided 2026-10-23)

Spatial claims (C1–C3), unit = grid node (0.05° spacing; N per node and 0.5° radius as in the paper; 2008–2023; both estimators):
- low / high b = below / above the whole-period b of the same region and estimator (b₀, Mc = MAXC + 0.2), i.e. classes
  relative to the regional average b as in Schorlemmer et al. (2005); the share of low-b (high-b) nodes is the reference
  for rows that condition on a covariate class. b₀ is event-weighted (one MLE over all events) while the node map is
  area-weighted, so the split is not 50/50 (LMS classical 34 % low / 66 % high). Sensitivity (2026-10-30, Yuxin asked for
  the literature rule with references, then kept the plain split because every cell keeps its direction): rows `*_utsu`
  in claims_table.csv restrict the classes to nodes whose b differs from b₀ at Utsu (1992) p < 0.05 (`utsu_p` in
  claims_report.py, formula as printed in Wiemer & Wyss 1997; n₂ = the regional fit's sample size, n_c / n_p in
  thresholds_2008-2023.csv), the anomaly definition of Wiemer & Wyss (1997) and Schorlemmer et al. (2004); 50–84 % of
  the nodes are then 'not distinguishable from b₀', b-positive loses more because its node std is about twice the
  classical one, and the direction of every cell is unchanged. A 2σ Shi & Bolt rule (92–97 % agreement with the Utsu
  classes) and the node median as threshold (2026-10-29 checks, scratch only) also gave the same directions, except NW
  contractional nodes with the median threshold (62 % vs 50 %).
- strain rate: weighted least squares of Wang & Shen (2020) velocities (Gaussian weight σ = 60 km, stations within 150 km,
  ≥ 6 stations, node inside the station hull); high = above the region's node median; dilatation terciles = contraction /
  shear / extension. The SCB-fixed field (ds01) is the check for the four SW regions.
- strain-rate row (merged 2026-10-29 from the former 'low-b nodes with high strain rate' and 'high-b nodes with low strain
  rate' rows): share of all nodes with strain data whose b and strain rate lie on opposite sides of their thresholds
  (b vs b₀, strain rate vs the node median); chance = 50 % because the strain split is at the median. The two former
  rows are still written to claims_table.csv (C1_rate, C2_rate). Their denominators differ (e.g. LMS classical 1,371
  low-b vs 2,640 high-b nodes, total 4,011 for both estimators) because b₀ is the event-weighted whole-window b, not
  the median of the node values: dense aftershock zones pull b₀ towards their own b, and the many peripheral nodes sit
  on the other side.
- GCMT mechanism rows REMOVED from the table 2026-10-29 (still in claims_table.csv as C1_mech, C2_mech, and drawn in
  figS_C12): mechanisms M ≥ 5 (1976–2023) classified by the Frohlich (1992) P/T/B plunges; a row needed ≥ 10 mechanisms
  of the class with a valid b node within 0.1°, which only LMS thrust (21), WYN normal (13) and NW thrust/normal (85/14)
  reached (box totals: LMS 36, XSH 8, WYN 58, XJ 11, NCP 11, NW 268; strike-slip dominates outside LMS and NW).
- heat flow at a node = inverse-distance mean of sites within 50 km (≥ 2 sites); high = above the median of the covered nodes.
  Rows need ≥ 20 nodes (XSH has no heat-flow site). Hot springs were dropped from the analysis (2026-10-23); the files stay in
  `data/hotsprings/` unused.

Temporal claims (C4), unit = the seven main-text windows (Figs 4–10), each with its own N and Mc (Table S1); the
1970–2023 SI windows are not used (decided 2026-10-23). Rules revised 2026-10-29 so that every element has a published
source (purpose: a descriptive count of drops before / rises after the marked events for both estimators, not a precursor test):
- series: fixed-N sliding windows, value at the last event, both estimators on the same windows (van der Elst 2021;
  Gulia et al. 2018). Classical windows are re-cut at max(panel Mc, MAXC + 0.2 of the window) and need ≥ 50 events above
  that Mc, else NaN (Gulia & Wiemer 2019, Methods steps 5–6a; Gulia et al. 2018; Dascher-Cousineau et al. 2020); the panel Mc
  is the 'pre-cut' floor. b-positive: panel Mc, dmc = 0.1, no re-cut and no exclusion (van der Elst 2021: "no windowing is
  required"; Lippiello & Petrillo 2024); dmc = 0.2 (van der Elst's default) is run as the SI sensitivity
  (`alarms_C4_summary_dmc02.csv`).
- reference of a window = median of all earlier windows of the panel that share no event with it (Gulia & Wiemer 2019, 6a
  "median of all individual b values"; Gulia et al. 2018), leaving out windows that end inside the Gardner–Knopoff time
  window of an earlier marked event (Gulia & Wiemer 2019, step 4: the pre-event catalogue avoids earlier sequences; G-K
  window as in Gulia et al. 2018). Fewer than 5 eligible windows → the region's 2008–2023 whole-period b (Gulia & Wiemer
  2019, 6b: regional background as fallback). The Table S1 b is only drawn.
- targets = the events marked in panel (d) (every event with M ≥ the fifth-largest magnitude); the initiator/aftershock
  distinction was removed on 2026-10-26 (the `initiator` flag stays in events_C4.csv for information). For an aftershock
  target the reference is still the pre-sequence background (Gulia & Wiemer 2019, step 8). An M ≥ 6-only variant
  (`TARGET_MAG = 6.0`, 20 events, 15 of them in NW) was computed on 2026-10-29 and not adopted: only 11 of the 44 marked
  events are M ≥ 6 (LMS 2+1, XSH 1, WYN 1, NW 6; XJ and NCP none), so `fig_summary` panel (b) now draws one glyph per
  marked event (size = magnitude, M ≥ 6 labelled, left half classical / right half positive, filled = criterion met,
  dashed outline = aftershock of an earlier marked event) with the k/n of both estimators at the end of each strip.
  On 2026-10-29 Yuxin chose a TABLE instead of that panel for the main text (Form B, `c4_event_table`): one row per
  window, one cell per marked event in time order (magnitude, * = aftershock of an earlier marked event, then two
  symbols classical | positive: filled = drop before, open = no drop, dot = no valid window), drop k/n per estimator in
  the last column; the rise results go into the footnote (1/6 | 1/6 in XSH, 0 elsewhere). The glyph panel stays in
  `fig_summary` for now (one of the two will be used, not both).
- drop before = last window ending before the event ≥ 10 % below its reference (Gulia & Wiemer 2019, step 7, red).
- rise after = first post-event window ≥ 10 % above the same reference (Gulia & Wiemer 2019, step 7, green). Post window:
  b-positive = the first window whose events all follow the event (van der Elst 2021); classical = the first such window
  whose per-window Mc is back at its pre-event level (median per-window Mc of the reference windows; Gulia et al. 2018:
  events removed until Mc(t) matches the pre-cut level). The 7-day blind time of the earlier version is gone.
- false alarms: REMOVED from the table on 2026-10-29. Alarm episodes (runs of windows ≥ 10 % below the reference) and
  their 7-day M ≥ 4 follow-up are still written to alarms_C4.csv / alarms_C4_summary.csv (n_episodes, n_hits,
  hit_ref_share) for the record. Reason: the question Yuxin wants answered ("is a large b drop followed by a larger
  event?") needs an a priori target magnitude, validity window and spatial tolerance; the only published b-value
  false-alarm counts (Gulia & Wiemer 2019 Fig. 6; Gulia et al. 2024: red status after a mainshock, target = a larger
  event within 3 months inside the G-K radius) answer a different question (what follows an M ≥ 6), and at the regional
  scale of our windows any target/duration choice is arbitrary (7 d / M ≥ 4 had been examined, with 30 and 90 d).
- differences from the sources to state in the text: regional windows instead of fault boxes; the paper's N (up to 3000)
  instead of 150–400; targets include in-sequence aftershocks; the per-window MAXC lands on the doubled converted-Mw
  bins (1.7 → Mc 1.9 in LMS/XSH), so the classical per-window Mc is often above the panel Mc even in quiet periods.

History: on 2026-10-26 the reference was the mean of the earlier windows (step 4 of Gulia et al. 2024, SRL 95:3579; PDF in
`data/srl-2024163.1.pdf`) and the rise was measured against the pre-event b after a 7-day blind time. Their full
traffic-light procedure was tried on 2026-10-26 and dropped (`ftls.py` kept unused; its result: 0 larger events followed
any M>=6, 5 reds = false alarms).

Node percentages are area fractions (neighbourhoods overlap); node counts are printed with every percentage.

## Manuscript status (2026-10-30)

Decided with Yuxin: no new section. Table 1 (= Table Y of `out/tables_main.docx`, C4 Form B) and two paragraphs went into
Section 4.3 after the "decreases before M 5 events are common" paragraph; a short paragraph on the spatial counts closes
Section 4.2; Text S5 (rules), Table S7 (= Table X) and Figure S18 (`figS_C4.pdf` rasterised at 300 dpi) were added to the SI;
the letter got one list item under AE comment 1 and one paragraph under Reviewer 1 comment 1. The Luding sentences in 4.3, 4.4
and Conclusions (iii) now name the window: the 2019–2022 decline is in the 2008–2023 window (Fig. S5d), while in the 2022–2023
window the last windows before the mainshock are above the reference of Table 1 (1.40 | 1.22 vs 1.10 | 1.09). Wang & Shen
(2020) was added to the reference list as plain text (Zotero entry still to be made). The heat-flow compilation is cited in
Text S5 as Jiang et al. (2019), to be confirmed by Yuxin. Edit scripts (`apply_manuscript.py`, `apply_si.py`,
`apply_letter.py`, `docxlib.py`) lived in the session scratchpad only.
