/**
 * =========================================================================
 * Script: 5_End-to-End_Sentinel-2.js
 * Project: Temporal Segment GEE - Timor-Leste Land Cover Mapping
 * 
 * Description:
 * End-to-end integrated script combining:
 *  1. Candidate feature engineering for 2025 and 2020
 *     (seasonal bands, indices, deltas, texture, SAR, and topography)
 *  2. EECU-optimized extraction (single synthetic reconstruction per season)
 *  3. Stratified sampling on 2025 ground truth samples
 *  4. Optional Random Forest hyperparameter tuning
 *  5. Multi-probability Random Forest model training & accuracy/log-loss evaluation
 *  6. Multi-year inference (2025, 2020) & Bayesian spatial smoothing
 *  7. Exporting multi-band probability assets for both target years
 * =========================================================================
 */

/**
 * Setup 
 * 
*/
// Area of Interest (Timor-Leste + 10km buffer)
var cfg         = require('users/rg2icraf/Luma:config');
var aoi         = ee.FeatureCollection(cfg.AOI_ASSET).geometry();
//Import required modules
var spectral    = require('users/dmlmont/spectral:spectral');
var ccdc_module = require('users/rg2icraf/Luma:ccdc_sentinel2');
var feats       = require('users/rg2icraf/Luma:features');
var extraction  = require('users/rg2icraf/Luma:extract_split_data');
var clf         = require('users/rg2icraf/Luma:classify_model');
var tuning      = require('users/rg2icraf/Luma:tune_model');
var post        = require('users/rg2icraf/Luma:post_process');

/**
 * Configurations
 * 
*/

// Seasonal fraction settings — sourced from cfg
var T = cfg.T;

// Spectral indices to compute
var baseIndexList = ['BLFEI', 'DBSI', 'GBNDVI', 'MBI', 'MNDWI', 'MSAVI', 'NDMI', 'NDVI'];
// Indices requiring delta computation (wet - dry)
var deltaIndexList = ['BLFEI', 'DBSI', 'MBI', 'MNDWI', 'MSAVI', 'NDMI'];

// Target years to process
var targetYears = [2025, 2020];

// Class schema (18 classes) — sourced from cfg
var classIdToName = cfg.CLASS_SCHEMA;

// Bayesian smoothing parameters per class — sourced from cfg
var smoothnessTL = cfg.SMOOTHNESS;

// Execution toggles
var RUN_TUNING               = false; // Set to true to run RF grid search
var APPLY_BAYESIAN_SMOOTHING = true;  // Set to true to apply spatial Bayesian smoothing
/*
// Random Forest default hyperparameters (calibrated from tuning)
var rfOptions = {
  nTrees:  500,
  vSplit:  19,    // ~sqrt(36 features), or specify e.g. 19
  minLeaf: 1,
  seed:    42
};
*/
/**
 * Load CCDC Temporal Segmentation Data
 * 
*/

//define CCDC assets — asset IDs sourced from cfg
var ccdc_mosaic = cfg.loadCCDCTiles(cfg.S2_CCDC_TILES, aoi);
print(ccdc_mosaic);

/**
 * Build Feature set for the target year
 */

//function to build the feature stack for a given year using CCDC synthetic imagery
function buildYearFeatures(ccdc, year) {
  // 1. Synthetic imagery for wet & dry seasons (single pass per season)
  var synthWet = ccdc_module.getSyntheticStack(
    ccdc, year + T.wet, ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'], false
  );
  var synthDry = ccdc_module.getSyntheticStack(
    ccdc, year + T.dry, ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'], false
  );

  // 2. Select and rename reflective bands
  var wetBands = synthWet.select(
    ['blue',   'green',  'red',    'nir',    'swir1',  'swir2'],
    ['B2_wet', 'B3_wet', 'B4_wet', 'B8_wet', 'B11_wet', 'B12_wet']
  );
  var dryBands = synthDry.select(
    ['blue',   'green',  'red',    'nir',    'swir1',  'swir2'],
    ['B2_dry', 'B3_dry', 'B4_dry', 'B8_dry', 'B11_dry', 'B12_dry']
  );

  // 3. GLCM texture from dry-season NIR (sourced from feats)
  var glcmFeatures = feats.getGlcmTexture(synthDry.select('nir'));

  // 4. Compute spectral indices directly from existing synthetic stacks
  var wetIndices = feats.computeSeasonalIndices(synthWet, 'S2', baseIndexList, 'wet');
  var dryIndices = feats.computeSeasonalIndices(synthDry, 'S2', baseIndexList, 'dry');

  // 5. Compute phenological seasonal deltas (wet - dry)
  var deltaIndices = feats.computePhenologicalDeltas(wetIndices, dryIndices, deltaIndexList);

  // 6. CCDC Harmonics (NIR, RED, SWIR1)
  var harmImgs = ccdc_module.getHarmonics(
    ccdc, year + T.harm, ['nir', 'red', 'swir1']
  );

  // 7. Concatenate all candidate bands and select the exact RFECV features
  var candidateStack = ee.Image.cat([
    wetBands,
    dryBands,
    glcmFeatures,
    wetIndices,
    dryIndices,
    deltaIndices,
    harmImgs
  ]);

  return candidateStack.toFloat().clip(aoi);
}

// Annual PALSAR-2 ScanSAR & Seasonal S1 — sourced from feats module
var scansarFeatures2025 = feats.getScanSARFeatures(2025, aoi);
var scansarFeatures2020 = feats.getScanSARFeatures(2020, aoi);

// Sentinel-1 seasonal radar backscatter features
var s1_seasonal_2020 = feats.getSeasonalS1(2020, aoi);
var s1_seasonal_2025 = feats.getSeasonalS1(2025, aoi);

// Build feature stacks for each epoch from the unified S2 CCDC mosaic
var spectral_features_2025 = buildYearFeatures(ccdc_mosaic, 2025);
var spectral_features_2020 = buildYearFeatures(ccdc_mosaic, 2020);

// Topography — sourced from feats module
var landform = ee.Image('projects/ee-rg2v2/assets/TPI_landform_TL');
var twi = ee.Image('projects/ee-rg2/assets/Top_Wetness_index');
var topo = feats.getTopography(aoi, {
  landformImg: landform,
  twiImg:      twi
});

//Stack all features for each year
var final_stack_2025 = ee.Image.cat([
  spectral_features_2025,
  s1_seasonal_2025,
  scansarFeatures2025,
  topo
]).clip(aoi);
var final_stack_2020 = ee.Image.cat([
  spectral_features_2020,
  s1_seasonal_2020,
  scansarFeatures2020,
  topo
]).clip(aoi);
//define the training data
var roi = ee.FeatureCollection('projects/ee-rg2/assets/Samples/New_Samples_TL_2025');

//extract features
var samples = final_stack_2025.sampleRegions({
  collection: roi,
  properties: ['new_id'], 
  scale: 20,                 
  tileScale: 16,             
  geometries: false          //Table only 
});
print('Predictor bands:', final_stack_2025.bandNames());
print('Input sample count:', roi.size());
print('Extracted sample count:', samples.size());
print('Extracted samples by class:', samples.aggregate_histogram('new_id'));

//extract the table, perform RFECV in google collab
Export.table.toDrive({
  collection: samples,
  description: 'AFCLIM_Feature_extraction_S2_2025',
  fileFormat: 'CSV'
});

//optimize predictor stack
var selected_feature = [
  'B11_dry', 'B11_wet', 'B12_dry', 'B12_wet', 'B2_dry', 'B2_wet',
  'B3_dry', 'B3_wet', 'B4_dry', 'B4_wet', 'B8_dry', 'B8_wet',
  'BLFEI_dry', 'BLFEI_wet', 'DBSI_dry', 'GBNDVI_dry', 'GBNDVI_wet',
  'MBI_delta', 'MBI_dry', 'MBI_wet', 'MNDWI_dry', 'MNDWI_wet',
  'MSAVI_wet', 'NDMI_delta', 'NDMI_dry', 'NDVI_dry', 'NDVI_wet',
  'ScanSAR_HH_dB', 'ScanSAR_HV_dB', 'TWI', 'VH_dry', 'VH_wet',
  'VV_dry', 'VV_wet', 'elevation', 'red_amp', 'slope', 'swir1_amp'
];

var optimize_stack_25 = final_stack_2025.select(selected_feature)
var optimize_stack_20 = final_stack_2020.select(selected_feature)


//Model Training
var samples = extraction.stratifiedSplit(roi, optimize_stack_25, 'new_id', 20, 0.7, 42);
var trainPix = samples.trainingPixels;
var testPix  = samples.testingPixels;

print('Training pixel count:', trainPix.size());
print('Testing pixel count:', testPix.size());

// Optional hyperparameter grid search
if (RUN_TUNING) {
  var num_trees     = ee.List([100, 200, 400, 500]);
  var variablesplit = ee.List([6, 9, 12, 19]);

  var tuning_result = tuning.tuneRandomForest(
    trainPix,
    testPix,
    'new_id',
    optimize_stack_25.bandNames(),
    num_trees,
    variablesplit
  );
  var bestParams = tuning.getBestParams(tuning_result, 'accuracy');
  print('Best Hyperparameters found:', bestParams);
}
// Train Multi-Probability Random Forest Model
var rfOptions = { nTrees: 400, vSplit: 6, minLeaf: 1, seed: 42 };
var rf_model_2025 = clf.multiProbabilityClassification(
  trainPix,
  'new_id',
  classIdToName,
  optimize_stack_25,
  rfOptions
);

var trainedModel = rf_model_2025.trainedModel;
var classBands   = rf_model_2025.classBands;
var classValues  = rf_model_2025.classValues;

// Evaluate probability performance on holdout test set
var eval_prob = clf.evaluateProbabilities(
  trainedModel,
  testPix,
  'new_id',
  classIdToName
);
print('Overall Log Loss:',           eval_prob.get('overallLogLoss'));
print('Overall Accuracy:',        eval_prob.get('overallAccuracy'));
print('Per-Class Log Loss:',         eval_prob.get('perClassLogLoss'));
print('Per-Class Accuracy:',         eval_prob.get('perClassAccuracy'));
print('Per-Class Sample Count:',     eval_prob.get('perClassNSamples'));
print('Total Test Samples:',         eval_prob.get('nSamples'));

// =========================================================================
// 6. Multi-Year Multi-Probability Inference & Post-Processing
// =========================================================================

// 2025 Probability Image
var probs_2025 = rf_model_2025.probsImage;

// 2020 Probability Image
var probsArray_2020 = optimize_stack_20.classify(trainedModel);
var probs_2020 = probsArray_2020.rename('probs').arrayFlatten([classBands]);


// Apply Bayesian spatial smoothing
var finalProbs_2025, finalProbs_2020;

if (APPLY_BAYESIAN_SMOOTHING) {
  print('Applying Bayesian spatial smoothing to probability surfaces...');
  finalProbs_2025 = post.bayesianSmooth(probs_2025, classBands, smoothnessTL, 5, 0.5);
  finalProbs_2020 = post.bayesianSmooth(probs_2020, classBands, smoothnessTL, 5, 0.5);
} else {
  finalProbs_2025 = probs_2025;
  finalProbs_2020 = probs_2020;
}

// Export the (optionally smoothed) class probability bands as Earth Engine assets.
[
  {year: 2025, image: finalProbs_2025},
  {year: 2020, image: finalProbs_2020}
].forEach(function(result) {
  cfg.toAsset(result.image.clip(aoi),
    'projects/ee-rg2/assets/TL_S2_class_probabilities_' + result.year, {
      description:      'TL_S2_class_probabilities_' + result.year,
      region:           aoi,
      scale:            cfg.SCALE.S2,
      pyramidingPolicy: {'.default': 'mean'}
    });
});
