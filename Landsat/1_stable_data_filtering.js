/**
 * =========================================================================
 * Script: 1_stable_data_filtering.js
 * Project: Temporal Segment GEE - Timor-Leste Land Cover Mapping
 *
 * Description:
 * Evaluates CCDC temporal stability across target historical epochs
 * (2020, 2015, 2005, and full 2000-2025) using the unified 5-tile Landsat
 * CCDC mosaic from 5_end-to-end_Landsat.js.
 *
 * Stability Criteria per Epoch [startYear -> 2025]:
 *  1. Zero CCDC breakpoints (tBreak) inside [startYear, endYear].
 *  2. Valid segment coverage within TOL (1 year) at both startYear and 2025,
 *     belonging to the same continuous CCDC segment.
 *  3. Optional gradual-drift filter: cumulative linear trend drift in NIR and
 *     SWIR1 (|slope * deltaYears|) remains below MAX_CUMULATIVE_DRIFT.
 * =========================================================================
 */

var ccdc_module = require('users/rg2icraf/Luma:ccdc_landsat');

// 1. Area of Interest & Input Reference / Testing Samples
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry();

// Update this asset path to your uploaded testing / ground-truth FeatureCollection
var inputSamples = ee.FeatureCollection(
  'projects/ee-rg2/assets/Samples/Samples_TL_2025_augment_cocoa'
);

// 2. Load 5-Tile Unified Landsat CCDC Mosaic (2000-2025) — aligned with 5_end-to-end_Landsat.js
var img1 = ee.Image('projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_central');
var img2 = ee.Image('projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_east_timor');
var img3 = ee.Image('projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_eastbox');
var img4 = ee.Image('projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_off_dili');
var img5 = ee.Image('projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_westbox');

var landsat_ccdc = ee.ImageCollection([img1, img2, img3, img4, img5]).mosaic().clip(aoi);

// 3. Stability Hyperparameters (aligned with 5_end-to-end_Landsat.js)
var TOL                  = 1.0;   // Max segment gap tolerance in years
var CCDC_COEF_SCALE      = 10000; // Reflectance scaling factor in custom CCDC run
var APPLY_SLOPE_FILTER   = true;  // Filter out slow non-breakpoint succession/drift
var MAX_CUMULATIVE_DRIFT = 0.12;  // Max allowable |slope * deltaYears| in surface reflectance (0..1)

/**
 * Computes a binary stability mask (1 = stable, 0 = unstable) between startYear and endYear.
 *
 * @param {ee.Image} ccdcImage - Mosaicked CCDC image.
 * @param {number} startYear - Historical target start year (e.g., 2005.0).
 * @param {number} endYear - Reference end year (e.g., 2025.99).
 * @param {string} bandName - Output band name.
 * @returns {ee.Image} Binary stability band (1 = stable, 0 = changed/gap).
 */
function stableInPeriod(ccdcImage, startYear, endYear, bandName) {
  // Condition 1: No CCDC breakpoint inside [startYear, endYear]
  var tBreak = ccdcImage.select('tBreak');
  var changeProb = ccdcImage.select('changeProb');
  var breaksInWindow = tBreak.arrayMask(
    tBreak.gte(startYear).and(tBreak.lte(endYear)).and(changeProb.gt(0))
  );
  var noBreak = breaksInWindow.arrayLength(0).eq(0);

  // Condition 2: Valid segment coverage at both startYear and endYear within TOL,
  // and both epochs fall on the exact same continuous CCDC segment
  var segStart = ccdc_module.nearestSegment(ccdcImage, startYear + 0.50);
  var segEnd   = ccdc_module.nearestSegment(ccdcImage, Math.floor(endYear) + 0.50);

  var validCoverage = segStart.gap.lte(TOL)
    .and(segEnd.gap.lte(TOL))
    .and(segStart.idx.eq(segEnd.idx));

  var isStable = noBreak.and(validCoverage);

  // Condition 3 (Optional): Low gradual linear drift in NIR & SWIR1 across the period
  if (APPLY_SLOPE_FILTER) {
    var deltaYears = endYear - startYear;

    var nirSlope = ccdcImage.select('NIR_coefs')
      .arraySlice(0, segEnd.idx, segEnd.idx.add(1)).arrayProject([1])
      .arrayGet([1]).divide(CCDC_COEF_SCALE).abs();

    var swir1Slope = ccdcImage.select('SWIR1_coefs')
      .arraySlice(0, segEnd.idx, segEnd.idx.add(1)).arrayProject([1])
      .arrayGet([1]).divide(CCDC_COEF_SCALE).abs();

    var lowDrift = nirSlope.multiply(deltaYears).lte(MAX_CUMULATIVE_DRIFT)
      .and(swir1Slope.multiply(deltaYears).lte(MAX_CUMULATIVE_DRIFT));

    isStable = isStable.and(lowDrift);
  }

  return isStable.unmask(0).toByte().rename(bandName);
}

// 4. Build Multi-Epoch Stability Stack (2020, 2015, 2005, 2000) + Latest Break Year
var stable2020 = stableInPeriod(landsat_ccdc, 2020.0, 2025.99, 'stable_2020');
var stable2015 = stableInPeriod(landsat_ccdc, 2015.0, 2025.99, 'stable_2015');
var stable2005 = stableInPeriod(landsat_ccdc, 2005.0, 2025.99, 'stable_2005');
var stable2000 = stableInPeriod(landsat_ccdc, 2000.0, 2025.99, 'stable_2000');

// Most recent breakpoint year (0 if no break ever occurred)
var lastBreakYear = landsat_ccdc.select('tBreak')
  .arrayReduce(ee.Reducer.max(), [0])
  .arrayGet([0])
  .toFloat()
  .unmask(0)
  .rename('last_tBreak');

var stabilityStack = ee.Image.cat([
  stable2020,
  stable2015,
  stable2005,
  stable2000,
  lastBreakYear
]);

// 5. Extract Stability Flags onto Sample Features
// Uses ee.Reducer.min() for stability flags (so polygon geometries require all internal pixels to be stable)
// and ee.Reducer.max() for last_tBreak. For point geometries, min/max equals the exact pixel value.
var combinedReducer = ee.Reducer.min().forEachBand(
  stabilityStack.select(['stable_2020', 'stable_2015', 'stable_2005', 'stable_2000'])
).combine({
  reducer2: ee.Reducer.max().setOutputs(['last_tBreak']),
  sharedInputs: false
});

var samplesWithStability = stabilityStack.reduceRegions({
  collection: inputSamples,
  reducer: combinedReducer,
  scale: 30,
  tileScale: 4
});

// Filtered subsets per target historical epoch
var stableSamples2020 = samplesWithStability.filter(ee.Filter.eq('stable_2020', 1));
var stableSamples2015 = samplesWithStability.filter(ee.Filter.eq('stable_2015', 1));
var stableSamples2005 = samplesWithStability.filter(ee.Filter.eq('stable_2005', 1));
var stableSamples2000 = samplesWithStability.filter(ee.Filter.eq('stable_2000', 1));

// 6. Console Diagnostics
print('Total Input Samples (2025):', inputSamples.size());
print('Stable Samples (2020-2025):', stableSamples2020.size());
print('Stable Samples (2015-2025):', stableSamples2015.size());
print('Stable Samples (2005-2025):', stableSamples2005.size());
print('Stable Samples (2000-2025):', stableSamples2000.size());
print('Stable 2005-2025 Samples by Class:', stableSamples2005.aggregate_histogram('new_id'));

// 7. Map Visualization
Map.centerObject(aoi, 9);
var stabVis = {min: 0, max: 1, palette: ['d73027', '1a9850']};
Map.addLayer(stable2020, stabVis, 'CCDC Stability: 2020-2025', false);
Map.addLayer(stable2015, stabVis, 'CCDC Stability: 2015-2025', false);
Map.addLayer(stable2005, stabVis, 'CCDC Stability: 2005-2025', true);

// 8. Export Annotated Samples (contains stable_2020, stable_2015, stable_2005, stable_2000, last_tBreak)
// You can filter by `stable_2020 == 1`, `stable_2015 == 1`, or `stable_2005 == 1` downstream.
Export.table.toAsset({
  collection:  samplesWithStability,
  description: 'Samples_TL_2025_with_CCDC_Stability_Flags',
  assetId:     'projects/ee-rg2/assets/Samples/Samples_TL_2025_with_CCDC_Stability_Flags'
});

Export.table.toDrive({
  collection:     samplesWithStability,
  description:    'TL_Landsat_Samples_with_Stability_GeoJSON',
  fileNamePrefix: 'TL_Landsat_Samples_with_Stability',
  fileFormat:     'GeoJSON'
});

