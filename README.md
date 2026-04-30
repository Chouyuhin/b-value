# Spatial and Temporal Frequency-Magnitude Distribution Analysis of Mainland China Seismicity

This repository contains the data, analysis notebooks, and mapping scripts for a study of Gutenberg–Richter *b*-value variations across three major seismic regions of mainland China: the **Sichuan–Yunnan block** (川滇), the **North China Plain** (华北平原), and the **Tibetan Plateau** (西北/青藏高原).

---

## Overview

The *b*-value of the Gutenberg–Richter relation log N = *a* − *b*M encodes the relative proportion of small to large earthquakes and is widely used as a proxy for crustal stress state, material heterogeneity, and seismic hazard potential. This study applies both classical and *b*-positive estimation to a ~1.5 million-event catalog spanning 1970–2023 (CENCcat), computes the magnitude of completeness M_c at each spatial grid cell via the MAXC method, and produces 2-D ZMAP-style *b*-value maps for each region and multiple time windows. Results are compared against major fault systems (Longmenshan, Xianshuihe, Xiaojiang, Anninghe) and significant earthquakes (M ≥ 6.6) including the 2008 Wenchuan M8.0, 2014 Ludian M6.6, and 2022 Luding M6.6 events.

---

## Repository Structure

```text
b_value/
├── bigEQs.csv                        9 major earthquakes M ≥ 6.6 (reference events)
├── experiments/
│   ├── bvalue4mainland.ipynb         Main analysis entry point (mainland-scale)
│   ├── 1chuandian/                   Sichuan–Yunnan region (710 k events)
│   │   ├── region1_chuandian.csv
│   │   ├── 1.1Longmenshan/                Longmenshan fault zone
│   │   ├── 1.2Xianshuihe/                Xianshuihe fault zone
│   │   ├── 1.3Zuoxia/                  SW sub-region
│   │   ├── 1.4Xiaojiang/                  Xiaojiang fault zone
│   │   └── 1.5Zhongxia/                  Central sub-region
│   ├── 2huabei/                      North China Plain (90 k events)
│   │   ├── region2_huabei.csv
│   │   └── bvalue4huabei.ipynb
│   └── 3tibet/                       Tibetan Plateau (238 k events)
│       ├── region3_tibet.csv
│       └── bvalue4xibei.ipynb
├── plotting/
│   ├── gmt.conf                      GMT 6.5 defaults
│   ├── nonmapfig.ipynb               M-T scatter, magnitude-type histogram
│   ├── fig1mainmap/                  Figure 1: mainland seismicity overview map
│   ├── 3chuandian/                   Figure 2: Sichuan–Yunnan b-value maps
│   ├── 4huabei/                      Figure 3: North China Plain maps
│   └── 5tibet/                       Figure 4: Tibetan Plateau maps
└── faults/                           Fault trace coordinates (lon/lat TXT)
    ├── Xiaojiang_Fault_all.txt
    ├── xianshuihe.txt
    ├── anninghe.txt
    └── ...
```

> **Catalog files** (CENCcat_Mw_final.txt, ~62 MB) are stored outside the repository under `../catalogs/` and are not tracked by git due to size.

---

## Methods

| Step | Tool / Method |
| ---- | ------------- |
| Catalog ingestion | `pandas`, custom CENCcat parser |
| Magnitude of completeness M_c | MAXC (Maximum Curvature), KS-distance, MBS |
| *b*-value estimation | Classical maximum-likelihood; *b*-positive (Bayesian) |
| Temporal analysis | Sliding-window *b*-value series with significance testing |
| 2-D spatial mapping | ZMAP-style moving window, fixed N = 200 events per cell |
| Grid export | ASCII `b_value_2D.txt` (lon, lat, *b*, σ) → GMT |
| Map production | GMT 6.5, Lambert Azimuthal Equal-Area projection |

---

## Key Results

- Mainland-scale b-value: **b = 0.59** (classical MLE), **b = 0.78** (*b*-positive), M_c ≈ 1.9, N ≈ 718 k events above M_c.
- Significant spatial heterogeneity: lower *b* (higher stress) along active thrust segments of the Longmenshan fault; elevated *b* in extensional back-arc regions.
- Temporal *b*-value series show statistically significant (p < 0.05) excursions preceding and following the 2008 Wenchuan, 2014 Ludian/Jinggu, and 2022 Luding earthquakes.
- Post-2022 *b*-value maps reveal stress reorganization in the central Sichuan–Yunnan block.

---

## Requirements

```text
Python  ≥ 3.10
jupyter
pandas
numpy
matplotlib
seismostats       # b-value and Mc estimation
GMT               ≥ 6.5   (for map figures)
```

Install Python dependencies:

```bash
pip install jupyter pandas numpy matplotlib seismostats
```

Install GMT via conda or the [official GMT installer](https://www.generic-mapping-tools.org/download/).

---

## Reproducing the Analysis

### 1. Regional *b*-value analysis (Python / Jupyter)

```bash
# Mainland overview
jupyter nbconvert --execute --inplace experiments/bvalue4mainland.ipynb

# Sichuan–Yunnan subregion example
jupyter nbconvert --execute --inplace "experiments/1chuandian/1.4小江/bvalue4xiaojiang.ipynb"
```

> **Note:** Several notebooks reference absolute paths under `/Users/chouyuhin/_Harvard/b-value/catalogs/`. Update these to your local catalog path before running.

### 2. GMT map figures

```bash
# Main seismicity overview map (Figure 1)
bash plotting/fig1mainmap/fig1_mainlandmap.sh

# Sichuan–Yunnan b-value maps (Figure 2)
bash plotting/3chuandian/fig2A_chuandian.sh

# Subregion map for a specific fault zone
bash "experiments/1chuandian/1.4小江/subregion_mainmap_gmt.sh"
```

All GMT scripts must be run from the repository root so that relative paths to `plotting/gmt.conf` and fault trace files resolve correctly.

### 3. Non-map publication figures

```bash
jupyter nbconvert --execute --inplace plotting/nonmapfig.ipynb
```

Produces `Fig1b_MT.png` (magnitude–time scatter) and `Fig1c_MagType.svg` (magnitude-type histogram).

---

## Data

| File | Events | Region | Period |
| ---- | ------ | ------ | ------ |
| CENCcat_Mw_final.txt* | ~1.5 M | Mainland China | 1970–2023 |
| region1_chuandian.csv | 710,692 | Sichuan–Yunnan | 1970–2023 |
| region2_huabei.csv | 90,172 | North China Plain | 1970–2023 |
| region3_tibet.csv | 238,126 | Tibetan Plateau | 1970–2023 |
| bigEQs.csv | 9 | All regions | 1970–2022 |

\* Not tracked in git; stored in `../catalogs/`.

The catalog is derived from the **China Earthquake Networks Center (CENC)** unified catalog, with magnitude converted to moment magnitude M_w using established regional conversion relations.

---

## Citation

If you use this code or data in your research, please cite:

> Zhou, Y. (2025). *Spatial and Temporal b-Value Analysis of Mainland China Seismicity*. GitHub repository, [github.com/chouyuhin/b-value](https://github.com/chouyuhin/b-value).

---

## License

MIT License © 2025 Yuxin Zhou
