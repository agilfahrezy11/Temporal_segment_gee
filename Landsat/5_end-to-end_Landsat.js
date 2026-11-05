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

// Base indices to calculate from spectral library (Awesome Spectral Indices)
// 15 carefully selected non-thermal indices covering canopy, moisture, shadow, built, soil & residue
var baseIndexList = [
  // 1. Vegetation Canopy & Agroforestry
  'NDVI', 'EVI', 'MSAVI', 'GNDVI',
  // 2. Water, Moisture & Mountain Shadow
  'MNDWI', 'NDMI', 'AWEIsh',
  // 3. Built-Up & Impervious Surfaces
  'IBI', 'UI', 'NDBI', 'BLFEI',
  // 4. Bare Soil, Riverbeds & Dry Biomass Residue
  'BI', 'DBSI', 'MBI', 'NDTI'
];

// All 15 indices evaluated as wet-dry seasonal deltas (wet - dry)
var deltaIndexList = [
  'NDVI', 'EVI', 'MSAVI', 'GNDVI',
  'MNDWI', 'NDMI', 'AWEIsh',
  'IBI', 'UI', 'NDBI', 'BLFEI',
  'BI', 'DBSI', 'MBI', 'NDTI'
];

// Reflective optical bands
var landsat_wet_bands = ['B2_wet', 'B3_wet', 'B4_wet', 'B5_wet', 'B6_wet', 'B7_wet'];
var landsat_dry_bands = ['B2_dry', 'B3_dry', 'B4_dry', 'B5_dry', 'B6_dry', 'B7_dry'];

// Wet, Dry, and Delta index band names
var wetIndexBands   = baseIndexList.map(function(idx) { return idx + '_wet'; });
var dryIndexBands   = baseIndexList.map(function(idx) { return idx + '_dry'; });
var deltaIndexBands = deltaIndexList.map(function(idx) { return idx + '_delta'; });

// CCDC Harmonics
var harm = ['NIR_amp', 'NIR_phase', 'RED_amp', 'SWIR1_amp'];

// Topography & Geomorphology (Time-invariant across 2005-2025)
var elev = ee.Image('NASA/NASADEM_HGT/001')
            .select('elevation')
            .clip(aoi);
var terrain_product = ee.Terrain.products(elev).select('slope', 'aspect').toFloat();
var landform = ee.Image('projects/ee-rg2v2/assets/TPI_landform_TL')
                      .rename('landform')
                      .toByte()
                      .clip(aoi);
var twi = ee.Image('projects/ee-rg2/assets/Non_spectral_feature/Top_Wetness_index')
                  .rename('TWI')
                  .toFloat();

var topo_features = ee.Image.cat([elev, terrain_product, landform, twi]);
var topoBands = ['elevation', 'slope', 'aspect', 'landform', 'TWI'];

// GLCM Textural Feature Names (derived from dry-season NIR per year)
var glcmBandsList = ['nir_contrast', 'nir_ent', 'nir_asm', 'nir_corr'];

// ALOS PALSAR L-band SAR backscatter bands (dB)
var palsarBands = ['palsar_HH', 'palsar_HV', 'palsar_HH_HV_diff'];

// Full 73-feature candidate pool for RFECV
var allCandidateFeatures = landsat_wet_bands
  .concat(landsat_dry_bands)
  .concat(wetIndexBands)
  .concat(dryIndexBands)
  .concat(deltaIndexBands)
  .concat(harm)
  .concat(topoBands)
  .concat(glcmBandsList)
  .concat(palsarBands);

// Previous RFECV subset (for reference / deployment after running RFECV)
var optimalFeatures = ["AWEIsh_delta", "AWEIsh_dry", "AWEIsh_wet", "B2_dry", "B2_wet", "B3_dry", "B3_wet", 
                     "B4_dry", "B4_wet", "B5_dry", "B6_dry", "B6_wet", "B7_dry", "B7_wet", "BLFEI_dry", "BLFEI_wet", 
                     "DBSI_delta", "DBSI_dry", "DBSI_wet", "EVI_wet", "GNDVI_dry", "GNDVI_wet", "MBI_delta", "MBI_dry", "MNDWI_dry", 
                     "MNDWI_wet", "NDTI_delta", "NDTI_dry", "NDVI_wet", "NIR_phase", "RED_amp", "SWIR1_amp", 
                     "UI_dry", "UI_wet", "elevation", "nir_contrast", "palsar_HH", "palsar_HV", "slope"];

// Feature Pool Selector:
// Set to TRUE to extract and train on all 73 candidate features (e.g. for generating RFECV table)
// Set to FALSE to use your pruned optimalFeatures list after RFECV is complete
var USE_FULL_CANDIDATE_POOL = false;

var selectedFeatures = USE_FULL_CANDIDATE_POOL ? allCandidateFeatures : optimalFeatures;

print('Total Active Predictors (' + selectedFeatures.length + '):', selectedFeatures);

// Target years to process
var targetYears = [2025, 2020, 2015, 2005];

// Class schema (18 classes) — sourced from cfg
var classIdToName = cfg.CLASS_SCHEMA;

// Bayesian smoothing parameters per class — sourced from cfg
var smoothnessTL = cfg.SMOOTHNESS;

// Execution toggles
var RUN_TUNING               = false; // Set to true to run RF grid search
var APPLY_BAYESIAN_SMOOTHING = false;  // Set to true to apply spatial Bayesian smoothing

// Random Forest default hyperparameters (calibrated from tuning)
var rfOptions = {
  nTrees:  100,
  vSplit:  21,    // ~sqrt(36 features), or specify e.g. 19
  minLeaf: 1,
  seed:    42
};

// Export configuration — sourced from cfg
var EXPORT_ASSET_PREFIX = 'projects/gee-v1-510606/assets/LC_prob_';
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

  // 6. GLCM Textural Features on Dry-Season NIR (captures canopy volume roughness)
  var glcm = synthDry.select('nir').multiply(1000).toInt32().glcmTexture({size: 1});
  var glcmBands = glcm.select([
    'nir_contrast', 
    'nir_ent', 
    'nir_asm', 
    'nir_corr'
  ]);

  // 7. ALOS PALSAR L-band SAR Backscatter (HH, HV, HH-HV diff in dB) — sourced from feats
  var palsar = feats.getPalsarBackscatter(year, aoi);

  // 8. Concatenate all candidate bands including topography, GLCM texture & radar
  var candidateStack = ee.Image.cat([
    wetBands,
    dryBands,
    wetIndices,
    dryIndices,
    deltaIndices,
    harmImgs,
    topo_features,
    glcmBands,
    palsar
  ]);

  return candidateStack.select(selectedFeatures).toFloat().clip(aoi);
}

// Build feature stacks for each epoch from the unified Landsat CCDC mosaic
//print('Building feature stacks for 2025, 2015, and 2005 from unified Landsat CCDC mosaic...');
var features_2025 = buildYearFeatures(landsat_ccdc, 2025, TOL, CCDC_COEF_SCALE);
var features_2020 = buildYearFeatures(landsat_ccdc, 2020, TOL, CCDC_COEF_SCALE)
var features_2015 = buildYearFeatures(landsat_ccdc, 2015, TOL, CCDC_COEF_SCALE);
var features_2005 = buildYearFeatures(landsat_ccdc, 2005, TOL, CCDC_COEF_SCALE);
Map.addLayer(features_2020, {}, "Feature 2020")

//print('2025 Feature Stack:', features_2025);
//print('2015 Feature Stack:', features_2015);
//print('2005 Feature Stack:', features_2005);

// =========================================================================
// 5. Sample Extraction & Model Training (2025 Ground Truth)
// =========================================================================

var roi = ee.FeatureCollection('projects/ee-rg2/assets/Samples/Samples_TL_2025_augment_cocoa');

// Optional: Export 2025 Feature Matrix for RFECV in Python / Google Colab
var EXPORT_RFECV_SAMPLE_TABLE = true; // Set to true to export full CSV table to Google Drive

if (EXPORT_RFECV_SAMPLE_TABLE) {
  var sampled_rfecv = features_2025.sampleRegions({
    collection: roi,
    properties: ['new_id'],
    scale: 30,
    projection: 'EPSG:4326',
    tileScale: 4,
    geometries: false // tabular only for fast export and ML training
  });

  Export.table.toDrive({
    collection: sampled_rfecv,
    description: 'TL_Landsat_Predictor_Matrix_RFECV',
    fileFormat: 'CSV'
  });
  //print('RFECV Table Export registered: "TL_Landsat_Predictor_Matrix_RFECV" (Run via Tasks tab)');
}

// Stratified train/test split (75% train, 25% test)
var samples = extraction.stratifiedSplit(roi, features_2025, 'new_id', 30, 0.75, 42);
var trainPix = samples.trainingPixels;
var testPix  = samples.testingPixels;
var testFC   = samples.testFC; // Vector holdout set preserving original geometries & ground-truth attributes

//print('Training pixel count:', trainPix.size());
//print('Testing pixel count:', testPix.size());

// Export holdout testing vector data to Google Drive as GeoJSON
var EXPORT_TEST_DATA_GEOJSON = true;

if (EXPORT_TEST_DATA_GEOJSON) {
  Export.table.toDrive({
    collection:  testFC,
    description: 'TL_Landsat_2025_Testing_Data_GeoJSON',
    fileNamePrefix: 'TL_Landsat_2025_Testing_Data',
    fileFormat:  'GeoJSON'
  });
  print('Testing Data GeoJSON Export registered: "TL_Landsat_2025_Testing_Data_GeoJSON"');
}

// Optional hyperparameter grid search
if (RUN_TUNING) {
  print('Running Random Forest hyperparameter tuning...');
  var num_trees     = ee.List([100, 200, 400]);
  var variablesplit = ee.List([6, 9, 12, 19]);

  var tuning_result = tuning.tuneRandomForest(
    trainPix,
    testPix,
    'new_id',
    selectedFeatures,
    num_trees,
    variablesplit
  );
  var bestParams = tuning.getBestParams(tuning_result, 'accuracy');
  print('Best Hyperparameters found:', bestParams);
}
/*
accuracy: 0.7692307692307693
kappa: 0.7532857458068258
minLeafPopulation: 1
numberOfTrees: 100
variablesPerSplit: 12
*/
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
/*
// Evaluate probability performance on holdout test set
var eval_prob = clf.evaluateProbabilities(
  trainedModel,
  testPix,
  'new_id',
  classIdToName
);
print('--- 2025 Test Evaluation ---');
print('Overall Log Loss:',           eval_prob.get('overallLogLoss'));
print('Overall Accuracy:',        eval_prob.get('overallAccuracy'));
print('Per-Class Log Loss:',         eval_prob.get('perClassLogLoss'));
print('Per-Class Accuracy:',         eval_prob.get('perClassAccuracy'));
print('Per-Class Sample Count:',     eval_prob.get('perClassNSamples'));
print('Total Test Samples:',         eval_prob.get('nSamples'));
*/
// =========================================================================
// 6. Multi-Year Multi-Probability Inference & Post-Processing
// =========================================================================

// 2025 Probability Image
var probs_2025 = rf_model_2025.probsImage;

// 2015 Probability Image
var probsArray_2015 = features_2015.classify(trainedModel);
var probs_2015 = probsArray_2015.rename('probs').arrayFlatten([classBands]);

//2020
var probsArray_2020 = features_2020.classify(trainedModel)
var probs_2020 = probsArray_2020.rename('probs').arrayFlatten([classBands])

// 2005 Probability Image
var probsArray_2005 = features_2005.classify(trainedModel);
var probs_2005 = probsArray_2005.rename('probs').arrayFlatten([classBands]);

// Apply Bayesian spatial smoothing
var finalProbs_2025, finalProbs_2015, finalProbs_2020, finalProbs_2005;

if (APPLY_BAYESIAN_SMOOTHING) {
  print('Applying Bayesian spatial smoothing to probability surfaces...');
  finalProbs_2025 = post.bayesianSmooth(probs_2025, classBands, smoothnessTL, 5, 0.5);
  finalProbs_2020 = post.bayesianSmooth(probs_2020, classBands, smoothnessTL, 5, 0.5);
  finalProbs_2015 = post.bayesianSmooth(probs_2015, classBands, smoothnessTL, 5, 0.5);
  finalProbs_2005 = post.bayesianSmooth(probs_2005, classBands, smoothnessTL, 5, 0.5);
} else {
  finalProbs_2025 = probs_2025;
  finalProbs_2020 = probs_2020;
  finalProbs_2015 = probs_2015;
  finalProbs_2005 = probs_2005;
}
/*
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
*/
// =========================================================================
// 8. Export Multi-Band Probabilities to Assets
// =========================================================================

// Scale probabilities to Uint16 (0 - 10000 range) to save storage & bandwidth
var probYears = [
  {year: '2025', image: finalProbs_2025},
  {year: '2020', image: finalProbs_2020},
  {year: '2015', image: finalProbs_2015},
  {year: '2005', image: finalProbs_2005}
];

probYears.forEach(function(item) {
  var exportImg = item.image.multiply(10000).round().toUint16();
  cfg.toAsset(exportImg, EXPORT_ASSET_PREFIX + item.year + '_Landsat_TL', {
    description:      'probabilities_' + item.year + '_TL',
    region:           aoi,
    scale:            EXPORT_SCALE,
    pyramidingPolicy: {'.default': 'mean'}
  });
});

print('Probability export tasks registered for 2025, 2015, and 2005.');
print('Open the Tasks tab in the Code Editor to run the exports.');
