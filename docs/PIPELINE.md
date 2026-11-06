# Master Execution Pipeline & Architecture Guide

## 1. Overview & Architecture

The **Luma Geospatial Engine** (`Temporal_segment_gee`) provides an end-to-end, modular, and publication-ready Land Use / Land Cover (LULC) classification framework implemented on Google Earth Engine (GEE). It is designed to be completely project-agnostic, easily adaptable to any study region, sensor, and temporal resolution.

```
                      +---------------------------------------+
                      |       Raw Satellite Imagery           |
                      |   Sentinel-2 Harmonized (L2A) /       |
                      |   Landsat 5/7/8/9 Collection 2 SR     |
                      +-------------------+-------------------+
                                          |
                        [Select Temporal Mode in Config]
                                          |
                 +------------------------+------------------------+
                 |                                                 |
         (Mode: 'ccdc')                                   (Mode: 'composite')
                 v                                                 v
  +-------------------------------+                 +-------------------------------+
  | Stage 1: Cloud Score+ QA      |                 | Stage 1: QA Cloud Masking     |
  | Stage 2: Monthly Median Comp. |                 | Stage 2: Seasonal Percentile/ |
  | Stage 3: CCDC Harmonic Fit    |                 |          Median Aggregation   |
  |          (tStart, tEnd, coef) |                 |          (wet & dry seasons)  |
  +---------------+---------------+                 +---------------+---------------+
                  |                                                 |
                  +-----------------------+-------------------------+
                                          |
                                          v
                      +---------------------------------------+
                      | Stage 3: Covariate Feature Engine     |
                      | - Optical Wet/Dry Bands (B2-B12)      |
                      | - Spectral Indices (wet & dry)        |
                      | - Phenological Deltas (wet - dry)     |
                      | - GLCM Textures (dry-season NIR)      |
                      | - Radar Backscatter (S1 / PALSAR)     |
                      | - Topography (DEM, slope, TPI, TWI)   |
                      | - CCDC Harmonics (NIR, Red, SWIR1)    |
                      +-------------------+-------------------+
                                          |
                                          v
                      +---------------------------------------+
                      | Stage 4: Stratified Sample Split &    |
                      |          Random Forest Training       |
                      | - 70% Train / 30% Test Stratified     |
                      | - Optional Grid Search Tuning         |
                      | - MULTIPROBABILITY Soft Prediction    |
                      | - Log Loss, Overall Acc., Conf. Matrix|
                      +-------------------+-------------------+
                                          |
                                          v
                      +---------------------------------------+
                      | Stage 5: Spatial Regularization       |
                      | - Empirical Bayes Spatial Smoothing   |
                      | - Malinverni Confusion Index (CI)     |
                      | - False Hydrology / Roof Refinement   |
                      +-------------------+-------------------+
                                          |
                                          v
                      +---------------------------------------+
                      | Stage 6: Export Dispatcher            |
                      | - Probabilities: uint16 (mean pyr.)   |
                      | - Thematic LULC: uint8 (mode pyr.)    |
                      | - Superpixels:   int32 (sample pyr.)  |
                      | - Cloud-Optimized GeoTIFFs to Drive   |
                      +-------------------+-------------------+
                                          |
                                          v
                      +---------------------------------------+
                      | Stage 7: Local Python Post-Processing |
                      | (Improve_pixel_based_lc_classification)|
                      | - Windowed Tile Mosaicking            |
                      | - Seam Healing & Smoothing            |
                      | - Olofsson Good Practices Area Est.   |
                      +---------------------------------------+
```

---

## 2. Repository Structure

```
Temporal_segment_gee/
├── master_pipeline.js           # Single entry-point master execution orchestrator
├── src/
│   ├── config.js                # Centralized project configuration, schemas, exports
│   ├── features.js              # Sensor-agnostic feature engineering & index calculations
│   ├── inputs/
│   │   ├── s2.js                # Sentinel-2 collection preparation & Cloud Score+
│   │   ├── landsat.js           # Landsat 5/7/8/9 collection preparation & QA masking
│   │   └── composite.js         # Direct seasonal percentile/median composites
│   ├── ccdc/
│   │   ├── core.js              # Sensor-agnostic CCDC fitting, synthetic stacks, harmonics
│   │   ├── ccdc_sentinel2.js    # Backward-compatible S2 CCDC facade
│   │   └── ccdc_landsat.js      # Backward-compatible Landsat CCDC facade
│   └── classification/
│       ├── classify_model.js    # Multi-probability RF training & probability evaluation
│       ├── tune_model.js        # Random Forest grid search hyperparameter tuning
│       ├── Post-process.js      # Bayesian smoothing, SNIC OBIA, stability, hydrology rules
│       ├── extract_split_data.js# Stratified train/test pixel extraction
│       └── mask.js              # Ancillary environmental evidence layers
├── docs/
│   ├── PIPELINE.md              # This master workflow guide
│   ├── CCDC_TEMPORAL_SEGMENTATION.md # Detailed CCDC harmonic reconstruction math
│   ├── DESIGN_IMPROVEMENT_PLAN.md   # Architectural design specifications
│   └── HYBRID_PIXEL_OBJECT_WORKFLOW.md # Hybrid pixel-object assembly theory
├── Sentinel-2/                  # Sentinel-2 specific stage scripts (Timor-Leste testbed)
└── Landsat/                     # Landsat specific stage scripts (Timor-Leste testbed)
```

---

## 3. Quickstart: Running `master_pipeline.js`

To run the master pipeline for any project:

1. Open `master_pipeline.js` in Google Earth Engine Code Editor.
2. Edit Section 1 (`PROJECT CONFIGURATION & PARAMETERS`):
   ```javascript
   // Select Sensor & Mode
   var SENSOR        = 'S2';         // 'S2' or 'Landsat'
   var TEMPORAL_MODE = 'ccdc';       // 'ccdc' or 'composite'
   var REFERENCE_YEAR = 2025;        // Epoch for reference ground truth
   var TARGET_YEARS   = [2025, 2020];// Epochs to map

   // Area of Interest (pass null for default AOI or your own ee.Geometry)
   var USER_AOI = ee.FeatureCollection('your/aoi/asset').geometry();

   // Training Data
   var SAMPLES_ASSET = 'your/training/samples/asset';
   var CLASS_PROP    = 'class_id'; // property containing integer labels
   ```
3. Set execution toggles:
   - `RUN_TUNING: false` (set to `true` to run RF grid search on trees and variable split).
   - `TRAIN_MODEL: true` (trains multi-probability Random Forest).
   - `APPLY_BAYESIAN_SMOOTH: true` (applies empirical Bayes spatial regularization).
   - `EXPORT_PROBABILITIES: true` (exports 16-bit multi-class probability rasters).
   - `EXPORT_LULC_MAPS: true` (exports 8-bit discrete classification maps).
4. Click **Run**.
5. Check the **Tasks** tab in the GEE Code Editor to submit export jobs to Google Drive or Assets.

---

## 4. Key Engine Components

### A. Central Configuration (`src/config.js`)
All project constants are consolidated into `src/config.js`:
- `cfg.getAOI(customAoi)`: Dynamically returns custom AOI or falls back to default.
- `cfg.SCALE`: Canonical pixel resolutions (`S2`: 20m, `S2_fine`: 10m, `Landsat`: 30m).
- `cfg.T`: Calibrated seasonal fractions (`dry: 0.65`, `wet: 0.16`, `harm: 0.50`).
- `cfg.CLASS_SCHEMA`: 18-class canonical dictionary.
- `cfg.CLASS_BANDS` / `cfg.CLASS_VALUES`: Strictly ordered array sequences for multi-probability models.
- `cfg.CLASS_PALETTE`: Canonical 18-color hex palette.
- `cfg.SMOOTHNESS`: Empirically calibrated per-class prior variances for Bayesian smoothing.
- `cfg.toAsset()` & `cfg.toDrive()`: Automated export wrappers that enforce proper pyramiding policies and formats.

### B. Input Engines & Temporal Segmentation (`src/inputs/` & `src/ccdc/core.js`)
- **`src/inputs/s2.js`**: Filters Sentinel-2 L2A Harmonized, applies Cloud Score+ (`cs_cdf` or `cs`), standardizes band names to lowercase (`blue`, `green`, `red`, `nir`, `swir1`, `swir2`), and handles reflectance scaling ($0-1$ or $0-10,000$).
- **`src/inputs/landsat.js`**: Merges Landsat 5/7/8/9 Tier 1 Surface Reflectance with bitmask QA filtering (`QA_PIXEL`) and surface reflectance scaling.
- **`src/inputs/composite.js`**: Generates direct seasonal composites (median, mean, or percentile) without CCDC for rapid mapping or regions where CCDC assets are not precomputed.
- **`src/ccdc/core.js`**: Sensor-agnostic CCDC engine providing:
  - `runCCDC()`: Configurable CCDC fitting.
  - `nearestSegment()`: Vectorized nearest-segment lookup with temporal gap measurement.
  - `getSyntheticStack()`: Linear trend clamping to $[t_{\text{start}}, t_{\text{end}}]$ with unconstrained seasonal harmonics, unified reflectance bounds $[-0.05, 1.0]$, and coefficient case handling (`upper` vs `lower`).
  - `getHarmonics()`: Amplitude and phase extraction for NIR, Red, and SWIR1.
  - `stableInPeriod()`: Trajectory slope and breakpoint filtering for ground truth sample quality screening.

### C. Feature Engineering (`src/features.js`)
A single module that computes all candidate predictors:
- `computeSeasonalIndices(image, sensor, indexList, seasonSuffix)`: Evaluates Awesome Spectral Indices on wet and dry synthetic/composite stacks.
- `computePhenologicalDeltas(wetIndices, dryIndices, deltaList)`: Computes seasonal contrast $(\text{wet} - \text{dry})$ captures phenological dynamics (e.g., deciduous vs evergreen foliage, cropland cropping cycles).
- `getGlcmTexture(nirImage, options)`: Evaluates Gray-Level Co-occurrence Matrix texture metrics (contrast, entropy, ASM, correlation) on dry-season NIR.
- `getSeasonalS1(year, aoi)`: Generates dual-polarized SAR backscatter features (wet and dry VV, VH, and cross-ratio).
- `getTopography(aoi, options)`: Concatenates bare-earth DEM elevation, slope, aspect, Topographic Position Index (TPI) landforms, and Topographic Wetness Index (TWI).

### D. Classification & Spatial Post-Processing (`src/classification/`)
- **`classify_model.js`**: Trains Random Forest in `MULTIPROBABILITY` mode, returning soft class probability rasters across all 18 classes, and evaluates test samples computing Log Loss, Cohen's Kappa, and per-class accuracies.
- **`tune_model.js`**: Grid searches number of decision trees and split variables to prevent overfitting and identify optimal hyperparameters.
- **`Post-process.js`**:
  - `bayesianSmooth()`: Implements empirical Bayes spatial regularization (Camara et al., 2024 / *sits*), blending pixel probabilities with local neighborhood means based on per-class smoothness variance weights.
  - `computeStability()`: Computes winning class ($W$), runner-up class ($S$), Confusion Index ($CI = S / W$), and binary stability flag based on Malinverni & Sanità (2026).
  - `adaptiveBlend()`: Adaptively blends pixel-level probabilities with GEOBIA SNIC superpixel object-mean probabilities, dampening object weights where ambiguity or internal variance is high.
  - `suppressFalseHydrology()`: Intrinsic rule-based refinement resolving spectral overlap between urban asphalt/roofs, cleared land, and waterbodies/dry riverbeds.

---

## 5. EECU & Memory Management Rules

To prevent Earth Engine worker memory overflows (`User memory limit exceeded`) and computational timeouts:

1. **Single Synthetic Reconstruction Pass:**
   - Synthetic CCDC stacks are generated **once per season** (`synthWet`, `synthDry`). Spectral indices and textures are computed directly from these stacks rather than repeatedly querying CCDC array slices.
2. **Pyramiding Policies:**
   - **Continuous probabilities:** Must use `{ '.default': 'mean' }`. Never use `mode` on continuous probabilities.
   - **Thematic LULC classes:** Must use `{ '.default': 'mode' }` to prevent non-existent interpolated class IDs at coarse zoom levels.
   - **Segment ID images:** Must use `{ '.default': 'sample' }` to preserve integer cluster boundaries.
3. **Data Type Casting:**
   - Probability rasters are multiplied by $10,000$ and cast to `toUint16()`, cutting file size by 50% compared to `toFloat()`.
   - Thematic classifications are cast to `toUint8()` (classes $1-18$).
   - Cluster IDs are cast to `toInt32()`, with an `ID_STRIDE` (e.g. $300,000,000$) applied to prevent ID collisions across tile boundaries while remaining safely below the signed 32-bit integer ceiling ($2,147,483,647$).
4. **Cloud-Optimized GeoTIFF (COG) Exports:**
   - All Google Drive exports include `formatOptions: { cloudOptimized: true }`, ensuring compatibility with GDAL, QGIS, and Python rasterio windowed reads without full-raster decompression.

---

## 6. Handoff to Python Post-Processing

Downstream production assembly, windowed seam healing, and Olofsson area estimation reside in the external repository:
`C:\Users\AFahrezi\Documents\GitHub\Improve_pixel_based_lc_classification`

### Exported Artifacts Required by Python Pipeline:
1. `TL_S2_probs_{year}_{tile}_20m.tif`: 18-band uint16 probability rasters per regional tile or national extent.
2. `TL_SNIC_segments_{tile}_10m.tif` or `20m.tif`: int32 segment cluster rasters with unique global IDs.
3. `samples_extracted.csv`: Ground truth sample table with extracted features for RFECV feature selection.

Once exported to Google Drive, the Python scripts in `Improve_pixel_based_lc_classification` execute:
- Block-wise sliding-window hybrid blending.
- Edge seam mitigation along regional tile boundaries.
- Error matrix generation with stratified area-adjusted confidence intervals (Olofsson et al., 2014).
