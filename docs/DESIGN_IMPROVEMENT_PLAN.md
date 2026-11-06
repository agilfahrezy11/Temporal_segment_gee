# Codebase Design Improvement Plan
## Temporal Segment GEE — Timor-Leste LULC Pipeline

---

## Goals

The two non-negotiable properties of a research pipeline are:

- **Modular** — each piece of logic lives in exactly one place; changing it changes all callers.
- **Reproducible** — a new collaborator (or you, six months from now) can re-run the full pipeline
  from scratch to get the same outputs, without reverse-engineering asset IDs from script bodies
  or guessing which script to run in which order.

The codebase currently achieves neither well. The `src/` module layer exists and the CCDC and
Bayesian smoothing modules are genuinely good — but pipeline scripts bypass them, duplicate their
logic, and scatter configuration across 8+ files. And while scripts are numbered `0_` through `5_`,
the actual asset handoffs between steps are undocumented, several links are broken, and the export
parameters are hardcoded inconsistently in every script.

There are **six structural problems** to fix, grouped by goal.

---

## Executive Summary

The codebase has a solid module layer (`src/`) with genuinely deep modules for CCDC synthesis and
Bayesian smoothing. The problem is that pipeline scripts bypass or duplicate this module layer
instead of extending it. The result is the same logic appearing in 4–8 places simultaneously —
and because each copy drifts independently, it is now impossible to change a shared behaviour
(e.g. the Bayesian smoothness dictionary, the spectral index parameter set) without hunting down
every copy. Separately, the numbered pipeline sequence implies reproducibility but several
asset-generation links are missing or broken, making a fresh run from scratch impossible without
manual investigation.

There are **six structural problems**, listed from highest to lowest leverage for fixing.

---

## Problem 1 — The class schema and smoothness dictionary are scattered across 8+ scripts

**What it is.** The 18-class `classIdToName` mapping and the per-class `smoothnessTL` dictionary
are copied verbatim into every script that classifies or smooths. The values differ subtly between
copies (`Cleared_land: 0.5` in one script, `2.0` in another; `Built-up: 1` vs. `0.2`).

**Why it matters.** This is configuration masquerading as code. Changing a class name or smoothness
weight requires a manual search-and-edit across 8 files. The diverged values already cause quietly
inconsistent outputs between the Landsat and Sentinel-2 pipelines.

**The fix — create `src/config.js`**

A single module exporting everything that is true of the *project*, not of any particular script:

```js
// src/config.js
exports.AOI_ASSET    = 'projects/ee-rg2v2/assets/AOI_Timor_Leste_10km';
exports.CRS          = 'EPSG:32751';
exports.MAX_PIXELS   = 1e13;
exports.SEASONAL_T   = { dry: 0.65, wet: 0.16, harm: 0.50 };
exports.TOL_YEARS    = 1.0;

exports.CLASS_SCHEMA = {
  1: 'Primary_forest',      2: 'Secondary_forest',      3: 'Mangrove_forest',
  4: 'Coastal_forest',      5: 'Teak_plantation',        6: 'Eucalyptus_plantation',
  7: 'Coffe_agroforestry',  8: 'Cocoa_agroforestry',    9: 'Coconut_agroforestry',
  10: 'Mixed_garden',       11: 'Paddy_field',           12: 'Other_cropland',
  13: 'Grassland',          14: 'Shrubland',             15: 'Cleared_land',
  16: 'Built-up',           17: 'Waterbody',             18: 'Dry_riverbed'
};

// Canonical smoothness weights for Bayesian spatial smoothing.
// These are the single source of truth — do NOT copy this dict into pipeline scripts.
exports.SMOOTHNESS = {
  'Primary_forest':        4.5,  'Secondary_forest':      3.0,
  'Mangrove_forest':       3.0,  'Coastal_forest':        3.5,
  'Teak_plantation':       1.2,  'Eucalyptus_plantation': 1.2,
  'Coffe_agroforestry':    1.5,  'Cocoa_agroforestry':    1.5,
  'Coconut_agroforestry':  1.5,  'Mixed_garden':          1.2,
  'Paddy_field':           1.0,  'Other_cropland':        1.0,
  'Grassland':             1.5,  'Shrubland':             1.5,
  'Cleared_land':          0.5,  'Built-up':              0.2,
  'Waterbody':             4.0,  'Dry_riverbed':          0.1
};

// Landsat CCDC asset IDs (5-tile mosaic, 2000–2025)
exports.LANDSAT_CCDC_TILES = [
  'projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_central',
  'projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_east_timor',
  'projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_eastbox',
  'projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_off_dili',
  'projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_westbox'
];

// Sentinel-2 CCDC asset IDs (5-tile mosaic, 2020–2025)
exports.S2_CCDC_TILES = [
  'projects/ee-rg2v2/assets/CCDC_TL_2020_2025_20m_westbox',
  'projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_east-timor',
  'projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_off_dili',
  'projects/gee-v1-510606/assets/CCDC_TL_2020_2025_20m_central-timor',
  'projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_east_box'
];
```

Every pipeline script then becomes one line: `var cfg = require('users/rg2icraf/Luma:config');`
replacing 50+ lines of copied configuration per script.

**Impact:** High. Fixes 8 files, eliminates value drift, makes class schema a single source of truth.

---

## Problem 2 — `computeIndicesFromSynth` and `getSeasonalS1` are duplicated 4–5 times each

**What it is.** Three helper functions appear near-identically across pipeline scripts, never
extracted into a module:

| Function | Copies | Scripts |
|---|---|---|
| `computeIndicesFromSynth(synthImg, suffix)` | 4 | `L/3_integrated`, `L/5_end-to-end`, `S2/5_end-to-end`, `L/2_feature_engineering` |
| `getSeasonalS1(year, roi)` | 3 | `S2/2_time_series`, `S2/3_tuning`, `S2/5_end-to-end` |
| `getAllHarmonics()` / `nirHarmonics()` | 3 | `L/2_feature_engineering`, `L/2_1_feature_engineering`, `S2/2_time_series` |
| 5-tile CCDC mosaic loading | 5+ | nearly every script |

**Why it matters.** `computeIndicesFromSynth` renames bands to `B2/B3/...` notation before calling
the `spectral` library, then adds a season suffix. The Landsat and S2 copies have already diverged
in band mappings (e.g. `B8` vs `B11` for NIR). A bug in the index calculation would need to be
fixed in 4 places — and each fix risks introducing new divergence.

**The fix — create `src/features.js`**

```js
// src/features.js

/**
 * Computes a list of spectral indices from a CCDC synthetic stack.
 * Handles sensor-specific band renaming internally.
 *
 * @param {ee.Image} synthImg   - Synthetic stack with generic band names (blue, green, red, nir, swir1, swir2)
 * @param {string}   sensor     - 'S2' or 'Landsat' (controls B8 vs B11 naming for NIR)
 * @param {Array}    indexList  - Spectral index names from the `spectral` library
 * @param {string}   suffix     - Season suffix appended to each output band ('wet' or 'dry')
 * @returns {ee.Image} Multi-band image with bands named e.g. 'NDVI_wet', 'MNDWI_dry'
 */
exports.computeSeasonalIndices = function(synthImg, sensor, indexList, suffix) { ... };

/**
 * Builds Sentinel-1 GRD seasonal median composites (VV + VH, wet + dry).
 * Wet season = Dec–Apr, Dry season = May–Nov.
 *
 * @param {number}      year - Target year
 * @param {ee.Geometry} roi  - Region of interest
 * @returns {ee.Image} 4-band image: VH_dry, VH_wet, VV_dry, VV_wet
 */
exports.getSeasonalS1 = function(year, roi) { ... };

/**
 * Retrieves ALOS PALSAR / PALSAR-2 L-band SAR mosaic and converts to gamma0 dB.
 * Includes a 3-year fallback composite for years with missing orbital passes.
 *
 * @param {number}      year - Target year (proxy mapping: 2005→2007, 2025→2024)
 * @param {ee.Geometry} roi  - Region of interest
 * @returns {ee.Image} 3-band image: palsar_HH, palsar_HV, palsar_HH_HV_diff (all in dB)
 */
exports.getPalsarBackscatter = function(year, roi) { ... };

/**
 * Loads the 5-tile CCDC mosaic and clips to AOI.
 *
 * @param {Array}       tileIds - Array of GEE asset IDs (use cfg.LANDSAT_CCDC_TILES or cfg.S2_CCDC_TILES)
 * @param {ee.Geometry} aoi     - Clip geometry
 * @returns {ee.Image} Seamless CCDC mosaic
 */
exports.loadCCDCMosaic = function(tileIds, aoi) { ... };

/**
 * Builds the time-invariant topographic feature stack (elevation, slope, aspect, landform, TWI).
 *
 * @param {ee.Geometry} aoi - Clip geometry
 * @returns {ee.Image} 5-band image: elevation, slope, aspect, landform, TWI
 */
exports.getTopography = function(aoi) { ... };
```

**Impact:** High. `buildYearFeatures` in the end-to-end scripts shrinks to ~20 lines. Seasonal
calendar definitions (which months are "wet") become a single controlled value.

---

## Problem 3 — `tune_model.js` has a broken export and a redundant one

**What it is.**

- `trainMultiProbModel()` in `tune_model.js` calls `buildRfParams()`, which is a **private**
  function defined only inside `classify_model.js`. This throws a `ReferenceError` at runtime.
- `evaluateProbabilities()` in `tune_model.js` is a weaker, separately maintained version of the
  same-named function in `classify_model.js` — it's missing per-class breakdowns and uses
  `eps=1e-3` instead of `1e-7`.

**Why it matters.** Anyone who calls `tuning.trainMultiProbModel()` gets a silent runtime crash.
The duplicate `evaluateProbabilities` means callers using the wrong version get coarser metrics
without knowing it.

**The fix:**

1. Remove `trainMultiProbModel()` from `tune_model.js` — it duplicates
   `classify_model.multiProbabilityClassification` and is broken.
2. Remove `evaluateProbabilities()` from `tune_model.js` — redirect callers to
   `classify_model.evaluateProbabilities`.
3. Export `buildRfParams` from `classify_model.js` so it's available if genuinely needed elsewhere.

After the fix, `tune_model.js` has one clear purpose: **RF hyperparameter grid search**.

**Impact:** Medium. Prevents silent runtime bugs; clarifies module responsibilities.

---

## Problem 4 — `6_hybrid_pixel_object_assembly.js` bypasses the module system

**What it is.** This script imports no `src/` modules at all. It re-implements the full
pixel–object integration entirely inline — but with *more sophisticated* logic than what
`Post-process.js` currently provides:

| Logic | In `Post-process.js` | In `6_hybrid` |
|---|---|---|
| Object aggregation | `aggregateProbsByClusters` (simple mean) | `aggregate` + `integrate` (same) |
| Classification | `classifyFromProbs` (argmax) | `toClass` (same) |
| Stability assessment | ✗ | `computeStability` (W, S, CI, is_stable) |
| Adaptive object weights | ✗ | Dynamic decay by confusion index |
| Rule-based refinement | ✗ | `refineRules` (false hydrology suppression) |

**Why it matters.** The stability map and rule-based refinement in `6_hybrid` are the best
post-processing code in the project (Malinverni & Sanità 2026 methodology). But they exist only
as a script-local invention — unreachable from any other script, untestable via the module
interface, invisible to the Landsat pipeline.

**The fix — deepen `Post-process.js`**

Move the three core functions from `6_hybrid` into `Post-process.js` as exports:

```js
// Add to src/classification/Post-process.js

/**
 * Computes winning class (W), runner-up (S), Confusion Index (CI = S/W), and
 * binary stability map (CI <= tau => stable).
 * Based on Malinverni & Sanità (2026).
 *
 * @param {ee.Image}    probs        - Class probability image
 * @param {Array}       classBands   - Probability band names
 * @param {ee.List}     classIds     - Integer class IDs (matching band order)
 * @param {number}      [tau=0.65]   - CI threshold below which a pixel is 'stable'
 * @returns {{ W, S, winner, runner_up, CI, is_stable }}
 */
exports.computeStability = function(probs, classBands, classIds, tau) { ... };

/**
 * Hybrid pixel–object integration with adaptive object trust.
 * Aggregates probabilities over SNIC superpixel objects, decays object weight
 * proportionally to object-level confusion index, then blends with pixel probs.
 *
 * @param {ee.Image}    probs      - Class probability image
 * @param {ee.Image}    segments   - Integer segment label band ('segment')
 * @param {ee.Image}    weights    - Per-class object trust weights (ee.Image)
 * @param {Object}      [options]  - { maxSize, adaptiveW, adaptiveSensitivity, useGuard, guardSd }
 * @returns {ee.Image} LULC + conf + ci + runner_up + is_stable bands
 */
exports.hybridIntegrate = function(probs, segments, weights, options) { ... };

/**
 * Rule-based post-classification refinement.
 * Suppresses false hydrology (Waterbody/Dry_riverbed) in urban and bare contexts
 * when CI exceeds configurable thresholds.
 *
 * @param {ee.Image}    hybridResult  - Output of hybridIntegrate (needs lulc, runner_up, ci)
 * @param {Object}      [options]     - { applyRules, mergeWaterRiverbed, hydroCiThresh, riverbedCiThresh }
 * @returns {ee.Image} hybridResult with 'lulc_refined' band added
 */
exports.refineHydrology = function(hybridResult, options) { ... };
```

Then `6_hybrid_pixel_object_assembly.js` becomes ~50 lines of configuration and wiring, and the
Landsat pipeline can use the same GEOBIA logic without a second implementation.

**Impact:** Medium-high. Turns the best GEOBIA logic in the codebase into a reusable module.

---

## Problem 5 — Export parameters are hardcoded inconsistently across every script

**What it is.** Every `Export.image.toAsset` and `Export.image.toDrive` call independently
hardcodes `crs`, `scale`, `maxPixels`, and `pyramidingPolicy`. There is no shared export helper.
As a result:

- `pyramidingPolicy: {'.default': 'sample'}` is used in some scripts for probability images where
  `'mean'` is correct — a silent resampling error that degrades classification quality on export.
- `scale` values (`10`, `20`, `30`) are scattered; it's easy to export a Sentinel-2 product at
  30m or vice versa with no warning.
- Adding a new export (e.g. a second region or a new year) requires copying and manually editing
  the full boilerplate block.

**The fix — add `exports.toAsset` and `exports.toDrive` helpers to `src/config.js`**

```js
// In src/config.js

/**
 * Submits a tiled export-to-asset job with canonical project defaults.
 *
 * @param {ee.Image} image       - Image to export
 * @param {string}   assetId     - Full GEE asset path
 * @param {Object}   options     - { description, region, scale, pyramidingPolicy }
 */
exports.toAsset = function(image, assetId, options) {
  options = options || {};
  Export.image.toAsset({
    image:            image,
    description:      options.description || assetId.split('/').pop(),
    assetId:          assetId,
    region:           options.region,
    scale:            options.scale || 10,
    crs:              exports.CRS,
    maxPixels:        exports.MAX_PIXELS,
    pyramidingPolicy: options.pyramidingPolicy || {'.default': 'mean'}
  });
};

/**
 * Submits a Cloud-Optimized GeoTIFF export to Drive with canonical project defaults.
 */
exports.toDrive = function(image, fileName, options) {
  options = options || {};
  Export.image.toDrive({
    image:          image,
    description:    options.description || fileName,
    folder:         options.folder || 'GEE_Exports',
    fileNamePrefix: fileName,
    region:         options.region,
    scale:          options.scale || 10,
    crs:            exports.CRS,
    maxPixels:      exports.MAX_PIXELS,
    formatOptions:  { cloudOptimized: true }
  });
};
```

**Impact:** Medium-high. Prevents silent CRS/scale/pyramiding mismatches. Makes adding a new
tile export a one-liner. Reproducibility-critical — wrong pyramiding on a probability image
produces a subtly wrong final map with no error message.

---

## Problem 6 — Missing scripts and broken links make a fresh run impossible

**What it is.** The numbered pipeline sequence (`0_` → `5_`) implies you can follow the numbers
to reproduce the full output. You cannot:

| # | Problem | Location |
|---|---|---|
| 6a | No script in the repo generates the 5-tile Landsat CCDC 2000–2025 mosaic | `CCDC_2000_2025_central` etc. are loaded by scripts 3, 4, 5 but never created here |
| 6b | `3_time_series_model_tuning.js` (S2) ends without any `Export` calls | `Sentinel-2/3_time_series_model_tuning.js` — probability images are never written to an asset |
| 6c | `west_box` geometry is referenced before it is defined | `Sentinel-2/1_temporal_segment_s2.js` — the `Export` call will throw at runtime |
| 6d | `2_1_feature_engineering_s2_2020.js` reimplements cloud masking inline | Does not import `ccdc_sentinel2`; uses a different band selection than the module; superseded by `5_end-to-end_sentinel-2.js` |
| 6e | No pipeline README documenting asset handoffs | Asset IDs produced by each step and consumed by the next are only traceable by reading all scripts |

**The fixes:**

- **6a:** Restore or create `Landsat/0_run_ccdc_2000_2025.js` with the CCDC parameters and
  tiling geometry used to generate those five assets.  If the script was run interactively and
  not saved, document the parameters as a comment block at the top of
  `Landsat/1_temporal_segment_landsat.js` at minimum.
- **6b:** Add the missing `Export.image.toAsset` calls to `Sentinel-2/3_time_series_model_tuning.js`
  for the 2025 and 2020 probability images, or delete the script and route callers to
  `5_end-to-end_sentinel-2.js` which already exports them.
- **6c:** Define all five tile geometries (`west_box`, `east_box`, etc.) before the export block in
  `Sentinel-2/1_temporal_segment_s2.js`.
- **6d:** Delete `Sentinel-2/2_1_feature_engineering_s2_2020.js` — its use case is fully covered
  by `5_end-to-end_sentinel-2.js` — or replace the inline masking with
  `require('users/rg2icraf/Luma:ccdc_sentinel2')`.
- **6e:** Add `docs/PIPELINE.md` documenting each numbered step, which GEE assets it reads, which
  assets it writes, and the approximate EECU cost. Template below.

**Pipeline documentation template (`docs/PIPELINE.md`):**

```markdown
# Pipeline Execution Guide

## Sentinel-2 Pipeline

| Step | Script | Reads | Writes | Notes |
|------|--------|-------|--------|-------|
| 0 | 0_test_ccdc_s2.js | — | — | Sanity check only, no export |
| 1 | 1_temporal_segment_s2.js | S2 SR Harmonized, Cloud Score+ | CCDC_TL_2020_2025_20m_{tile} ×5 | Run one export per tile |
| 2 | 5_end-to-end_sentinel-2.js | CCDC mosaic, training samples | LC_prob_2025, LC_prob_2020 | Tuning optional (RUN_TUNING flag) |
| 3 | 6_hybrid_pixel_object_assembly.js | LC_prob_2025, SNIC segments | lulc_2025_hybrid | Requires SNIC from step 4_1 |
| 4 | 4_1_SNIC_segmentations.js | CCDC mosaic | snic/seg_{tile} ×5 | Run before step 3 |

## Landsat Pipeline
...
```

**Impact:** High for reproducibility. Without this, a fresh run requires reverse-engineering all
asset dependencies by reading every script body.

---

## Implementation Order

Each step is tagged **[M]** (modularity), **[R]** (reproducibility), or both.

| Step | Goal | Task | Files affected | Effort |
|---|---|---|---|---|
| **1** | M + R | Create `src/config.js` with export helpers, AOI, CRS, class schema, smoothness, tile IDs, seasonal fractions | New file + ~8 pipeline scripts | 1–2 days |
| **2** | M | Create `src/features.js` (`computeSeasonalIndices`, `getSeasonalS1`, `getPalsarBackscatter`, `loadCCDCMosaic`, `getTopography`) | New file + ~5 pipeline scripts | 1–2 days |
| **3** | M | Fix `tune_model.js`: remove broken `trainMultiProbModel` and duplicate `evaluateProbabilities`; export `buildRfParams` from `classify_model.js` | 2 src files | ~1 hour |
| **4** | M | Deepen `Post-process.js`: add `computeStability`, `hybridIntegrate`, `refineHydrology` from `6_hybrid` | `Post-process.js` + `6_hybrid_pixel_object_assembly.js` | 1 day |
| **5** | R | Fix broken pipeline links: missing tile geometries, missing S2 tuning exports, restore/document Landsat CCDC generation script | `1_temporal_segment_s2.js`, `3_time_series_model_tuning.js` | Half day |
| **6** | R | Write `docs/PIPELINE.md` with step-by-step execution guide, asset handoff table, and EECU cost notes | New file | Half day |

---

## What Not to Change

- The `buildYearFeatures` + `targetYears.map()` pattern in the end-to-end scripts is the right
  level of abstraction for a pipeline script. Keep it.
- `ccdc_landsat.js` and `ccdc_sentinel2.js` are genuinely deep modules. The only worthwhile
  internal change is extracting `nearestSegment` and `fillIsolated` (verbatim duplicates between
  the two files) into a shared `src/ccdc/ccdc_utils.js`.
- `classify_model.js` is well-designed. No structural changes needed.
- `Post-process.bayesianSmooth` is a good deep module. Don't simplify its interface — add to it
  (Problem 4).

---

## Module Depth Summary (current state)

| Module | Depth | Notes |
|---|---|---|
| `ccdc_sentinel2.getSyntheticStack` | **Deep** | Nearest-segment, harmonic eval, plausibility masking, hole fill |
| `ccdc_landsat.getSyntheticLandsatStack` | **Deep** | Same complexity; duplicates `nearestSegment` / `fillIsolated` from S2 module |
| `Post-process.bayesianSmooth` | **Deep** | Logit transform, robust local stats, Bayesian update |
| `classify_model.multiProbabilityClassification` | **Medium** | Encapsulates RF + array flatten; class-order sorting non-trivial |
| `classify_model.evaluateProbabilities` | **Medium** | Per-class per-sample log loss server-side |
| `ccdc_sentinel2.getS2TimeSeries` | **Medium** | Hides `linkCollection` + Cloud Score+ pattern |
| `ccdc_landsat.getLandsatTimeSeries` | **Medium** | Hides 4-sensor merge + QA_PIXEL bitmask |
| `tune_model.tuneRandomForest` | **Shallow** | Grid search plumbing only |
| `extract_split_data.stratifiedSplit` | **Shallow** | Thin boilerplate over `randomColumn` + `sampleRegions` |
| `mask.build` | **Shallow** | Dataset loading + threshold comparisons |
| `Post-process.snicSuperpixels` | **Shallow** | Pure pass-through to `ee.Algorithms.Image.Segmentation.SNIC` |
