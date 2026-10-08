//define the AOI
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry()
//Load the two CCDC Assets
var ccdc_old = ee.ImageCollection('GOOGLE/GLOBAL_CCDC/V1').filterBounds(aoi).mosaic(); // 1999-2019
var ccdc_new = ee.Image('projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_TL_2020_2025_Landsat'); // 2020-2025
var ccdc_module = require('users/rg2icraf/Luma:ccdc_landsat');
/**
 * Check for structural breaks within a specific time window
 * Returns an Image: 1 if ANY break occurred in window, 0 if clean/stable
 */
function checkBreakInPeriod(ccdcImg, startYear, endYear) {
  var tBreak = ccdcImg.select('tBreak'); // 1D Array Image of break dates
  
  // 1. Create a boolean mask: 1 where break date falls inside [startYear, endYear], 0 otherwise
  var maskInWindow = tBreak.gte(startYear).and(tBreak.lte(endYear));
  
  // 2. Filter the tBreak array to keep ONLY break dates that occurred inside the window
  var breaksInWindow = tBreak.arrayMask(maskInWindow);
  
  // 3. Count how many breaks happened inside the window
  var numBreaks = breaksInWindow.arrayLength(0);
  
  // 4. Return 1 if number of breaks > 0 (Unstable), else 0 (Stable)
  return numBreaks.gt(0).rename('has_break');
}
/**
 * Check 3: Spectral Boundary Continuity Check (2019.9 vs 2020.1)
 * Evaluates whether a hidden land change happened right at the split boundary.
 */
// Evaluate breaks in both assets
var break_2005_2019 = checkBreakInPeriod(ccdc_old, 2005.0, 2019.9);
var break_2020_2025 = checkBreakInPeriod(ccdc_new, 2020.0, 2025.9);

var synth_2019 = ccdc_module.getSyntheticLandsatStack(ccdc_old, 2019.9);
var synth_2020 = ccdc_module.getSyntheticLandsatStack(ccdc_new, 2020.1);

// Calculate spectral difference across key bands (e.g., NIR and SWIR1)
var diffNIR = synth_2019.select('nir').subtract(synth_2020.select('nir')).abs();
var diffSWIR = synth_2019.select('swir1').subtract(synth_2020.select('swir1')).abs();

// Threshold: If reflectance jumped by > 0.08 (~8%), treat as boundary break
var boundaryShift = diffNIR.gt(0.08).or(diffSWIR.gt(0.08));

/**
 * Combine all conditions into a single Stability Mask
 * 1 = Stable across 2005-2025
 * 0 = Unstable / Changed
 */
var isStable = break_2005_2019.eq(0)
  .and(break_2020_2025.eq(0))
  .and(boundaryShift.eq(0))
  .rename('is_stable');

//define the training data
var roi = ee.FeatureCollection('projects/ee-rg2/assets/Samples/Samples_TL_hierarchy_2025_v1');
// Filter training points to keep ONLY stable samples
var stableTrainingPoints = roi.map(function(feature) {
  var stabilityStatus = isStable.reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: feature.geometry(),
    scale: 30
  }).get('is_stable');
  
  return feature.set('is_stable', stabilityStatus);
}).filter(ee.Filter.eq('is_stable', 1));

print('Original 2025 Sample Count:', roi.size());
print('Filtered Stable (2005-2025) Sample Count:', stableTrainingPoints.size());