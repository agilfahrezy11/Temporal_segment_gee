/**
 * =========================================================================
 * Script: 3_Landsat_integrated_features_and_tuning.js
 * Project: Temporal Segment GEE - Timor-Leste Land Cover Mapping
 * 
 * Description:
 * End-to-end integrated script combining:
 *  1. On-the-fly CCDC feature engineering for 2025, 2015, and 2005
 *     (using RFECV optimal features: seasonal bands, indices, deltas, harmonics)
 *  2. EECU-optimized extraction (single synthetic reconstruction per season)
 *  3. Stratified sampling on 2025 ground truth samples
 *  4. Optional Random Forest hyperparameter tuning
 *  5. Multi-probability Random Forest model training & accuracy/log-loss evaluation
 *  6. Multi-year inference (2025, 2015, 2005) & Bayesian spatial smoothing
 *  7. Exporting multi-band probability assets for each target year
 * =========================================================================
 */

// =========================================================================
// 1. Setup & Modules
// =========================================================================

// Area of Interest (Timor-Leste + 10km buffer)
var cfg         = require('users/rg2icraf/Luma:config');
var aoi         = ee.FeatureCollection(cfg.AOI_ASSET).geometry();

// Import required modules
var spectral    = require('users/dmlmont/spectral:spectral');
var ccdc_module = require('users/rg2icraf/Luma:ccdc_landsat');
var feats       = require('users/rg2icraf/Luma:features');
var extraction  = require('users/rg2icraf/Luma:extract_split_data');
var clf         = require('users/rg2icraf/Luma:classify_model');
var tuning      = require('users/rg2icraf/Luma:tune_model');
var post        = require('users/rg2icraf/Luma:post_process');

// =========================================================================
// 2. Configuration & Hyperparameters
// =========================================================================

// Seasonal fraction settings — sourced from cfg
var T   = cfg.T;
var TOL = cfg.TOL_YEARS;

// RFECV optimal features selected
var landsat_wet_25 = ['B2_wet', 'B3_wet', 'B4_wet', 'B5_wet', 'B6_wet', 'B7_wet'];
var landsat_dry_25 = ['B2_dry', 'B3_dry', 'B4_dry', 'B5_dry', 'B6_dry', 'B7_dry'];
var index_ = [
  'BLFEI_delta', 'BLFEI_dry', 'BLFEI_wet',
  'DBSI_delta',  'DBSI_dry',  'DBSI_wet',
  'EVI_wet',
  'GBNDVI_dry',  'GBNDVI_wet',
  'MBI_delta',   'MBI_dry',
  'MNDWI_delta', 'MNDWI_dry', 'MNDWI_wet',
  'MSAVI_delta', 'MSAVI_wet',
  'NDMI_delta',  'NDMI_dry',
  'NDVI_dry',    'NDVI_wet'
];
var harm = ['NIR_amp', 'NIR_phase', 'RED_amp', 'SWIR1_amp'];

// Combine into final 36-feature optimal predictor list
var optimalFeatures = landsat_wet_25
  .concat(landsat_dry_25)
  .concat(index_)
  .concat(harm);

print('Total Selected Predictors (' + optimalFeatures.length + '):', optimalFeatures);

// Base indices to calculate from spectral library
var baseIndexList = ['BLFEI', 'DBSI', 'EVI', 'GBNDVI', 'MBI', 'MNDWI', 'MSAVI', 'NDMI', 'NDVI'];
// Indices requiring delta computation (wet - dry)
var deltaIndexList = ['BLFEI', 'DBSI', 'MBI', 'MNDWI', 'MSAVI', 'NDMI'];

// Target years to process
var targetYears = [2025, 2015, 2005];

// Class schema (18 classes) — sourced from cfg
var classIdToName = cfg.CLASS_SCHEMA;

// Bayesian smoothing parameters per class — sourced from cfg
var smoothnessTL = cfg.SMOOTHNESS;

// Execution toggles
var RUN_TUNING               = false; // Set to true to run RF grid search
var APPLY_BAYESIAN_SMOOTHING = true;  // Set to true to apply spatial Bayesian smoothing

// Random Forest default hyperparameters (calibrated from tuning)
var rfOptions = {
  nTrees:  500,
  vSplit:  6,    // ~sqrt(36 features), or specify e.g. 19
  minLeaf: 1,
  seed:    42
};

// Export configuration — sourced from cfg
var EXPORT_ASSET_PREFIX = 'projects/ee-rg2/assets/LC_prob_';
var EXPORT_SCALE        = cfg.SCALE.Landsat;

// =========================================================================
// 3. Load CCDC Dataset (Unified 2000-2025 Mosaic)
// =========================================================================

// Load tiled CCDC products — asset IDs sourced from cfg
var landsat_ccdc = cfg.loadCCDCTiles(cfg.LANDSAT_CCDC_TILES, aoi);

// Coefficient scale for custom Landsat CCDC run (0-10,000 range)
var CCDC_COEF_SCALE = cfg.CCDC_COEF_SCALE.own;

// =========================================================================
// 4. Feature Engineering Functions (EECU-Optimized)
// =========================================================================

/**
 * Builds the complete RFECV optimal feature stack for a given target year.
 * EECU Optimization: Reconstructs synthetic Landsat ONLY ONCE per season.
 */
function buildYearFeatures(ccdc, year, tol, coefScale) {
  // 1. Synthetic imagery for wet & dry seasons (single pass per season)
  var synthWet = ccdc_module.getSyntheticLandsatStack(ccdc, year + T.wet, false, tol, coefScale);
  var synthDry = ccdc_module.getSyntheticLandsatStack(ccdc, year + T.dry, false, tol, coefScale);

  // 2. Select and rename reflective bands
  var wetBands = synthWet.select(
    ['blue',   'green',  'red',    'nir',    'swir1',  'swir2'],
    ['B2_wet', 'B3_wet', 'B4_wet', 'B5_wet', 'B6_wet', 'B7_wet']
  );
  var dryBands = synthDry.select(
    ['blue',   'green',  'red',    'nir',    'swir1',  'swir2'],
    ['B2_dry', 'B3_dry', 'B4_dry', 'B5_dry', 'B6_dry', 'B7_dry']
  );

  // 3. Compute spectral indices directly from existing synthetic stacks
  var wetIndices = feats.computeSeasonalIndices(synthWet, 'Landsat', baseIndexList, 'wet');
  var dryIndices = feats.computeSeasonalIndices(synthDry, 'Landsat', baseIndexList, 'dry');

  // 4. Compute phenological seasonal deltas (wet - dry)
  var deltaIndices = feats.computePhenologicalDeltas(wetIndices, dryIndices, deltaIndexList);

  // 5. CCDC Harmonics (NIR, RED, SWIR1)
  var harmImgs = ccdc_module.getHarmonics(ccdc, year + T.harm, ['NIR', 'RED', 'SWIR1'], tol, coefScale);

  // 6. Concatenate all candidate bands and select the exact RFECV features
  var candidateStack = ee.Image.cat([
    wetBands,
    dryBands,
    wetIndices,
    dryIndices,
    deltaIndices,
    harmImgs
  ]);

  return candidateStack.select(optimalFeatures).toFloat().clip(aoi);
}

// Build feature stacks for each epoch from the unified Landsat CCDC mosaic
print('Building feature stacks for 2025, 2015, and 2005 from unified Landsat CCDC mosaic...');
var features_2025 = buildYearFeatures(landsat_ccdc, 2025, TOL, CCDC_COEF_SCALE);
var features_2015 = buildYearFeatures(landsat_ccdc, 2015, TOL, CCDC_COEF_SCALE);
var features_2005 = buildYearFeatures(landsat_ccdc, 2005, TOL, CCDC_COEF_SCALE);

print('2025 Feature Stack:', features_2025);
print('2015 Feature Stack:', features_2015);
print('2005 Feature Stack:', features_2005);

// =========================================================================
// 5. Sample Extraction & Model Training (2025 Ground Truth)
// =========================================================================

var roi = ee.FeatureCollection('projects/ee-rg2/assets/Samples/Samples_TL_hierarchy_2025_v1');

// Stratified train/test split (70% train, 30% test)
var samples = extraction.stratifiedSplit(roi, features_2025, 'new_id', 30, 0.7, 42);
var trainPix = samples.trainingPixels;
var testPix  = samples.testingPixels;

print('Training pixel count:', trainPix.size());
print('Testing pixel count:', testPix.size());

// Optional hyperparameter grid search
if (RUN_TUNING) {
  print('Running Random Forest hyperparameter tuning...');
  var num_trees     = ee.List([100, 200, 400, 500]);
  var variablesplit = ee.List([6, 9, 12, 19]);

  var tuning_result = tuning.tuneRandomForest(
    trainPix,
    testPix,
    'new_id',
    optimalFeatures,
    num_trees,
    variablesplit,
    ee.List([1]),
    42
  );
  var bestParams = tuning.getBestParams(tuning_result, 'accuracy');
  print('Best Hyperparameters found:', bestParams);
}

// Train Multi-Probability Random Forest Model
print('Training MultiProbability Random Forest Model...');
var rf_model_2025 = clf.multiProbabilityClassification(
  trainPix,
  'new_id',
  classIdToName,
  features_2025,
  rfOptions
);

var trainedModel = rf_model_2025.trainedModel;
var classBands   = rf_model_2025.classBands;
var classValues  = rf_model_2025.classValues;

// Evaluate probability performance on holdout test set (overall & per-class)
var eval_prob = clf.evaluateProbabilities(
  trainedModel,
  testPix,
  'new_id',
  classIdToName
);
print('--- 2025 Test Evaluation ---');
print('Overall Log Loss:', eval_prob.get('overallLogLoss'));
print('Overall Top-Class Accuracy:', eval_prob.get('overallAccuracy'));
print('Per-Class Log Loss:', eval_prob.get('perClassLogLoss'));
print('Per-Class Accuracy:', eval_prob.get('perClassAccuracy'));
print('Per-Class Sample Count:', eval_prob.get('perClassNSamples'));
print('Total Test Samples Evaluated:', eval_prob.get('nSamples'));

// =========================================================================
// 6. Multi-Year Multi-Probability Inference & Post-Processing
// =========================================================================

// 2025 Probability Image
var probs_2025 = rf_model_2025.probsImage;

// 2015 Probability Image
var probsArray_2015 = features_2015.classify(trainedModel);
var probs_2015 = probsArray_2015.rename('probs').arrayFlatten([classBands]);

// 2005 Probability Image
var probsArray_2005 = features_2005.classify(trainedModel);
var probs_2005 = probsArray_2005.rename('probs').arrayFlatten([classBands]);

// Apply Bayesian spatial smoothing
var finalProbs_2025, finalProbs_2015, finalProbs_2005;

if (APPLY_BAYESIAN_SMOOTHING) {
  print('Applying Bayesian spatial smoothing to probability surfaces...');
  finalProbs_2025 = post.bayesianSmooth(probs_2025, classBands, smoothnessTL, 5, 0.5);
  finalProbs_2015 = post.bayesianSmooth(probs_2015, classBands, smoothnessTL, 5, 0.5);
  finalProbs_2005 = post.bayesianSmooth(probs_2005, classBands, smoothnessTL, 5, 0.5);
} else {
  finalProbs_2025 = probs_2025;
  finalProbs_2015 = probs_2015;
  finalProbs_2005 = probs_2005;
}

// Generate quick preview land cover maps (argmax)
var lc_preview_2025 = post.classifyFromProbs(finalProbs_2025, classBands, classValues);
var lc_preview_2015 = post.classifyFromProbs(finalProbs_2015, classBands, classValues);
var lc_preview_2005 = post.classifyFromProbs(finalProbs_2005, classBands, classValues);

// =========================================================================
// 7. Map Visualization
// =========================================================================

Map.centerObject(aoi, 9);

var rgbVis = {bands: ['B4_dry', 'B3_dry', 'B2_dry'], min: 0.02, max: 0.22};
Map.addLayer(features_2025, rgbVis, '2025 Dry Season RGB (Synthetic)', false);
Map.addLayer(features_2015, rgbVis, '2015 Dry Season RGB (Synthetic)', false);
Map.addLayer(features_2005, rgbVis, '2005 Dry Season RGB (Synthetic)', false);

var classVis = {
  min: 1,
  max: 18,
  palette: [
    '#00441B', // 1  Primary_forest
    '#238B45', // 2  Secondary_forest
    '#006D2C', // 3  Mangrove_forest
    '#66C2A4', // 4  Coastal_forest
    '#8C6D31', // 5  Teak_plantation
    '#A6761D', // 6  Eucalyptus_plantation
    '#7FC97F', // 7  Coffe_agroforestry
    '#66A061', // 8  Cocoa_agroforestry
    '#41AB5D', // 9  Coconut_agroforestry
    '#A1D76A', // 10 Mixed_garden
    '#C7E9B4', // 11 Paddy_field
    '#FDD835', // 12 Other_cropland
    '#D9EF8B', // 13 Grassland
    '#C2A878', // 14 Shrubland
    '#d7d6d6', // 15 Cleared_land
    '#DB310D', // 16 Built-up
    '#0e8acc', // 17 Waterbody
    '#D2B48C'  // 18 Dry_riverbed
  ]
};

Map.addLayer(lc_preview_2025, classVis, 'Preview LC 2025 (Argmax)', false);
Map.addLayer(lc_preview_2015, classVis, 'Preview LC 2015 (Argmax)', false);
Map.addLayer(lc_preview_2005, classVis, 'Preview LC 2005 (Argmax)', false);

// =========================================================================
// 8. Export Multi-Band Probabilities to Assets
// =========================================================================

// Scale probabilities to Uint16 (0 - 10000 range) to save storage & bandwidth
var probYears = [
  {year: '2025', image: finalProbs_2025},
  {year: '2015', image: finalProbs_2015},
  {year: '2005', image: finalProbs_2005}
];

probYears.forEach(function(item) {
  var exportImg = item.image.multiply(10000).round().toUint16();
  cfg.toAsset(exportImg, EXPORT_ASSET_PREFIX + item.year + '_Landsat_TL', {
    description:      'smoothed_probabilities_' + item.year + '_TL',
    region:           aoi,
    scale:            EXPORT_SCALE,
    pyramidingPolicy: {'.default': 'mean'}
  });
});

print('Probability export tasks registered for 2025, 2015, and 2005.');
print('Open the Tasks tab in the Code Editor to run the exports.');
