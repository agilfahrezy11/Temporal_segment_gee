# Temporal Segmentation & Land Cover Classification — Timor-Leste

A Google Earth Engine (GEE) pipeline for multi-temporal land use/land cover (LULC) mapping of Timor-Leste using CCDC (Continuous Change Detection and Classification) temporal segmentation, Random Forest classification, and Bayesian post-processing.

Produces wall-to-wall LULC maps for **2020** and **2025** at 10 m (Sentinel-2) and 30 m (Landsat) resolution across an 18-class hierarchical land cover schema.

---

## Land Cover Classes

| ID | Class | ID | Class |
|---|---|---|---|
| 1 | Primary Forest | 10 | Mixed Garden |
| 2 | Secondary Forest | 11 | Paddy Field |
| 3 | Mangrove Forest | 12 | Other Cropland |
| 4 | Coastal Forest | 13 | Grassland |
| 5 | Teak Plantation | 14 | Shrubland |
| 6 | Eucalyptus Plantation | 15 | Cleared Land |
| 7 | Coffee Agroforestry | 16 | Built-up |
| 8 | Cocoa Agroforestry | 17 | Waterbody |
| 9 | Coconut Agroforestry | 18 | Dry Riverbed |

---

## Pipeline Overview

```
Raw Imagery (S2 / Landsat)
        │
        ▼
  Cloud Masking & Scaling
        │
        ▼
  CCDC Temporal Segmentation  ◄──  Google Global CCDC (1999–2019)
  (per-pixel harmonic fit +
   break detection 2020–2025)
        │
        ▼
  Feature Engineering
  ┌─────────────────────────────────────┐
  │  CCDC synthetic stacks (wet & dry)  │
  │  Harmonic coefficients (amp, phase) │
  │  Spectral indices (NDVI, EVI, MBI…) │
  │  SAR backscatter (S1 VV/VH)         │
  │  Topography (DEM, slope, TWI, TPI)  │
  │  Distance layers (road, urban)      │
  └─────────────────────────────────────┘
        │
        ▼
  Stability Filtering
  (remove training samples near CCDC breaks)
        │
        ▼
  Random Forest Classification
  (MULTIPROBABILITY mode, 400 trees)
        │
        ▼
  Post-Processing
  ┌────────────────────────────────┐
  │  Bayesian spatial smoothing    │
  │  Minimum Mapping Unit (MMU)    │
  │  Temporal consistency rules    │
  └────────────────────────────────┘
        │
        ▼
  Final LULC Maps (2020 & 2025)
```

---

## Repository Structure

```
Temporal_segment_gee/
├── src/                            # Reusable GEE modules
│   ├── ccdc/
│   │   ├── ccdc_landsat.js         # Landsat 8/9 CCDC pipeline
│   │   └── ccdc_sentinel2.js       # Sentinel-2 CCDC pipeline
│   └── classification/
│       ├── classify_model.js       # RF classification (hard / soft / multi-probability)
│       ├── extract_split_data.js   # Stratified train/test split
│       ├── tune_model.js           # RF hyperparameter grid search
│       ├── Post-process.js         # Bayesian smoothing, MMU, argmax
│       └── mask.js                 # Evidence mask layers
│
├── Sentinel-2/                     # S2 pipeline scripts (run in order)
│   ├── 0_test_ccdc_s2.js
│   ├── 1_temporal_segment_s2.js
│   ├── 2_1_feature_engineering_s2_2020.js
│   ├── 2_time_series_feature_engineering.js
│   ├── 3_time_series_model_tuning.js
│   └── 4_timer_series_classification_assembly.js
│
└── Landsat/                        # Landsat pipeline scripts (run in order)
    ├── 0_test_landsat_ccdc.js
    ├── 1_stable_data_filtering.js
    ├── 2_temporal_segment_landsat.js
    └── 3_feature_engineering_landsat_2025.js
```

Scripts are numbered to indicate execution order. Each step depends on GEE assets exported by the previous step.

---

## Key Dependencies

| Module | GEE Path | Description |
|---|---|---|
| `ccdc_sentinel2` | `users/rg2icraf/Luma:ccdc_sentinel2` | S2 time series prep + CCDC |
| `ccdc_landsat` | `users/rg2icraf/Luma:ccdc_landsat` | Landsat time series prep + CCDC |
| `classify_model` | `users/rg2icraf/Luma:classify_model` | RF classification modes |
| `extract_split_data` | `users/rg2icraf/Luma:extract_split_data` | Train/test splitting |
| `tune_model` | `users/rg2icraf/Luma:tune_model` | Hyperparameter tuning |
| `post_process` | `users/rg2icraf/Luma:post_process` | Spatial smoothing & consistency |
| `mask` | `users/rg2icraf/Luma:mask` | Ancillary evidence masks |
| `spectral` | `users/dmlmont/spectral:spectral` | Spectral index computation |

---

## Usage

All scripts run in the [GEE Code Editor](https://code.earthengine.google.com). There is no local build step.

1. Copy the script into the GEE Code Editor (or sync via GEE's Git integration)
2. Click **Run**
3. For export tasks, open the **Tasks** tab and click **Run** on each pending task
4. Outputs land in either a GEE Asset or Google Drive depending on the export call

> **Note:** Large-area CCDC exports are tiled into sub-regions to stay within GEE memory limits. Tile outputs are later mosaicked with `ee.ImageCollection([...]).mosaic()`.

---

## Output Specs

| Property | Value |
|---|---|
| CRS | `EPSG:32751` (UTM Zone 51S) |
| Resolution | 10 m (Sentinel-2), 30 m (Landsat) |
| Format | Cloud Optimized GeoTIFF |
| Bit depth | `Byte` (classification), `Float32` (probabilities) |
| Target years | 2020, 2025 |

---

## Study Area

Timor-Leste (East Timor), defined by `projects/ee-rg2v2/assets/AOI_Timor_Leste_10km`.

---

## References

- Zhu, Z. & Woodcock, C.E. (2014). Continuous change detection and classification of land cover using all available Landsat data. *Remote Sensing of Environment*, 144, 152–171.
- Camara, G. et al. (2024). Bayesian smoothing for land cover classification. [sits book](https://e-sensing.github.io/sitsbook/cl_smoothing.html)
- Montero, D. et al. (2023). A standardized catalogue of spectral indices. *Scientific Data*, 10, 197. [`spectral` library](https://github.com/davemlz/spectral)
