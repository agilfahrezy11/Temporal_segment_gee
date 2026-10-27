/**
 * =========================================================================
 * Script: 1_stable_data_filtering_s2.js
 * Project: Temporal Segment GEE - Timor-Leste Land Cover Mapping (Sentinel-2)
 *
 * Description:
 * Dedicated CCDC temporal stability filtering for Sentinel-2 (20 m) across
 * target epochs (2025 reference year and 2020-2025 multi-temporal transfer),
 * aligned with 5_end-to-end_sentinel-2.js and src/ccdc/ccdc_sentinel2.js.
 *
 * Key Sentinel-2 Adaptations (vs. Landsat):
 *  1. Uses the 5-tile 20 m Sentinel-2 CCDC mosaic (2019-2025).
 *  2. Uses lowercase CCDC coefficient band names ('nir_coefs', 'swir1_coefs').
 *  3. Verifies continuous segment coverage across seasonal reconstruction dates
 *     (T.wet = year + 0.16 through T.dry = year + 0.65) so no sample is masked
 *     out by ccdc_sentinel2.getSyntheticStack()'s strict [tStart, tEnd] check.
 *  4. Evaluates gradual linear drift (|slope * deltaYears|) on 'nir_coefs' and
 *     'swir1_coefs' to filter out slow non-breakpoint transitions.
 * =========================================================================
 */

// 1. Area of Interest & Input Reference / Testing Samples
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry();

// Update this asset path to your 2025 ground-truth or testing FeatureCollection
var inputSamples = ee.FeatureCollection('projects/ee-rg2/assets/Ground_Reference_2025_combine');

// 2. Load 5-Tile Sentinel-2 CCDC Mosaic (20 m) — aligned with 5_end-to-end_sentinel-2.js
var tl_ccdc1 = ee.Image('projects/ee-rg2v2/assets/CCDC_TL_2020_2025_20m_westbox');
var tl_ccdc2 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_east-timor');
var tl_ccdc3 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_off_dili');
var tl_ccdc4 = ee.Image('projects/gee-v1-510606/assets/CCDC_TL_2020_2025_20m_central-timor');
var tl_ccdc5 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_east_box');

var s2_ccdc = ee.ImageCollection([tl_ccdc1, tl_ccdc2, tl_ccdc3, tl_ccdc4, tl_ccdc5])
  .mosaic()
  .clip(aoi);

// 3. Stability Hyperparameters (aligned with 5_end-to-end_sentinel-2.js)
var T = {
  wet:  0.16, // ~March (Peak wet season)
  harm: 0.50, // Mid-year harmonic evaluation
  dry:  0.65  // ~August (Peak dry season)
};

// In ccdc_sentinel2.js, getSyntheticStack requires strict [tStart <= t <= tEnd] (TOL = 0.0).
// Increase TOL (e.g. 0.25) only if you want lenient boundary tolerance.
var TOL                  = 0.0;
var CCDC_COEF_SCALE      = 10000; // Reflectance scaling factor in S2 CCDC run
var APPLY_SLOPE_FILTER   = true;  // Filter out slow non-breakpoint succession/drift
var MAX_CUMULATIVE_DRIFT = 0.10;  // Max allowable |slope * deltaYears| in surface reflectance (0..1)

/**
 * Finds the nearest CCDC segment index and temporal gap (in years) at fractional year t.
 * Prevents empty-array slicing errors when extracting segment coefficients.
 */
function nearestSegment(ccdcImage, t) {
  t = ee.Number(t);
  var tS = ccdcImage.select('tStart');
  var tE = ccdcImage.select('tEnd');
  // 0 if t is inside [tStart, tEnd]; >0 if t falls outside the segment
  var dist = tS.subtract(t).max(tE.multiply(-1).add(t)).max(0);
  return {
    idx: dist.multiply(-1).arrayArgmax().arrayGet([0]),
    gap: dist.arrayReduce(ee.Reducer.min(), [0]).arrayGet([0]).rename('gap_years'),
    tS:  tS,
    tE:  tE
  };
}

/**
 * Computes a binary stability mask (1 = stable, 0 = unstable) between startYear and endYear.
 *
 * @param {ee.Image} ccdcImage - Mosaicked Sentinel-2 CCDC image.
 * @param {number} startYear - Historical target start year (e.g., 2020).
 * @param {number} endYear - Reference end year (e.g., 2025).
 * @param {string} bandName - Output band name.
 * @returns {ee.Image} Binary stability band (1 = stable, 0 = changed/gap).
 */
function stableInPeriodS2(ccdcImage, startYear, endYear, bandName) {
  var windowStart = startYear;
  var windowEnd   = endYear + 0.99;

  // Condition 1: Zero confirmed CCDC breakpoints (changeProb == 1) inside [startYear, endYear + 0.99]
  // Note: 1_temporal_segment_s2.js uses cs_cdf >= 0.50 without TMask on monthly medians, which can
  // leave unconfirmed tail anomalies (0 < changeProb < 1) on the final segment. Checking changeProb.eq(1)
  // ensures only confirmed segment-splitting breaks are flagged.
  var tBreak     = ccdcImage.select('tBreak');
  var changeProb = ccdcImage.select('changeProb');
  var confirmedBreak = tBreak.multiply(changeProb.eq(1));
  var breaksInWindow = confirmedBreak.arrayMask(
    confirmedBreak.gte(windowStart).and(confirmedBreak.lte(windowEnd))
  );
  var noBreak = breaksInWindow.arrayLength(0).eq(0);

  // Condition 2: Continuous single-segment coverage from startYear wet season (startYear + 0.16)
  // through endYear dry season (endYear + 0.65), guaranteeing valid synthetic stacks in both years
  var segStartWet = nearestSegment(ccdcImage, startYear + T.wet);
  var segEndDry   = nearestSegment(ccdcImage, endYear + T.dry);

  var validCoverage = segStartWet.gap.lte(TOL)
    .and(segEndDry.gap.lte(TOL))
    .and(segStartWet.idx.eq(segEndDry.idx));

  var isStable = noBreak.and(validCoverage);

  // Condition 3 (Optional): Low gradual linear drift in lowercase 'nir_coefs' & 'swir1_coefs'
  if (APPLY_SLOPE_FILTER) {
    var deltaYears = Math.max(1.0, (endYear + T.dry) - (startYear + T.wet));

    var nirSlope = ccdcImage.select('nir_coefs')
      .arraySlice(0, segEndDry.idx, segEndDry.idx.add(1)).arrayProject([1])
      .arrayGet([1]).divide(CCDC_COEF_SCALE).abs();

    var swir1Slope = ccdcImage.select('swir1_coefs')
      .arraySlice(0, segEndDry.idx, segEndDry.idx.add(1)).arrayProject([1])
      .arrayGet([1]).divide(CCDC_COEF_SCALE).abs();

    var lowDrift = nirSlope.multiply(deltaYears).lte(MAX_CUMULATIVE_DRIFT)
      .and(swir1Slope.multiply(deltaYears).lte(MAX_CUMULATIVE_DRIFT));

    isStable = isStable.and(lowDrift);
  }

  return isStable.unmask(0).toByte().rename(bandName);
}

// 4. Build Sentinel-2 Multi-Epoch Stability Stack
// - stable_2025: Single-year stability in 2025 (for 2025 training/testing sample purity)
// - stable_2020: Multi-year stability from 2020 to 2025 (for 2020 accuracy assessment)
var stable2025 = stableInPeriodS2(s2_ccdc, 2025, 2025, 'stable_2025');
var stable2020 = stableInPeriodS2(s2_ccdc, 2020, 2025, 'stable_2020');

// Most recent confirmed breakpoint year (0 if no confirmed break ever occurred in 2019-2025)
var lastBreakYear = s2_ccdc.select('tBreak')
  .multiply(s2_ccdc.select('changeProb').eq(1))
  .arrayReduce(ee.Reducer.max(), [0])
  .arrayGet([0])
  .toFloat()
  .unmask(0)
  .rename('last_tBreak');

var stabilityStack = ee.Image.cat([
  stable2025,
  stable2020,
  lastBreakYear
]);

// 5. Extract Stability Flags onto Sample Features (at 20 m EPSG:32751 native CCDC grid)
var combinedReducer = ee.Reducer.min().forEachBand(
  stabilityStack.select(['stable_2025', 'stable_2020'])
).combine({
  reducer2: ee.Reducer.max().setOutputs(['last_tBreak']),
  sharedInputs: false
});

var samplesWithStability = stabilityStack.reduceRegions({
  collection: inputSamples,
  reducer:    combinedReducer,
  scale:      20,
  crs:        'EPSG:32751',
  tileScale:  4
});

// Filtered subsets
var stableSamples2025 = samplesWithStability.filter(ee.Filter.eq('stable_2025', 1));
var stableSamples2020 = samplesWithStability.filter(ee.Filter.eq('stable_2020', 1));

// 6. Console Diagnostics
print('Total Input Samples (2025):', inputSamples.size());
print('Stable Samples in 2025 (stable_2025 == 1):', stableSamples2025.size());
print('Stable Samples 2020-2025 (stable_2020 == 1):', stableSamples2020.size());
print('Stable 2020-2025 Samples by Class:', stableSamples2020.aggregate_histogram('Class_id'));

// 7. Map Visualization
Map.centerObject(aoi, 9);
var stabVis = {min: 0, max: 1, palette: ['d73027', '1a9850']};
Map.addLayer(stable2025, stabVis, 'S2 CCDC Stability: 2025 Only', false);
Map.addLayer(stable2020, stabVis, 'S2 CCDC Stability: 2020-2025', true);

// 8. Export Annotated Samples (contains stable_2025, stable_2020, last_tBreak)
Export.table.toAsset({
  collection:  samplesWithStability,
  description: 'Samples_TL_S2_2025_with_CCDC_Stability_Flags',
  assetId:     'projects/ee-rg2/assets/Samples/Samples_TL_S2_2025_with_CCDC_Stability_Flags'
});

Export.table.toDrive({
  collection:     samplesWithStability,
  description:    'TL_Sentinel2_Samples_with_Stability_GeoJSON',
  fileNamePrefix: 'TL_Sentinel2_Samples_with_Stability',
  fileFormat:     'GeoJSON'
});
