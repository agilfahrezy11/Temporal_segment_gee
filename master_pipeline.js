/**
 * =========================================================================
 * SCRIPT: master_pipeline.js
 * PROJECT: Luma Geospatial Engine - Land Use / Land Cover Mapping
 * 
 * DESCRIPTION:
 * Universal Master Execution Pipeline for Multi-Temporal LULC Classification.
 * Orchestrates the end-to-end workflow from raw observations or CCDC temporal
 * segments to publication-grade thematic maps and probability rasters.
 * 
 * CAPABILITIES:
 *  1. Multi-Sensor: Sentinel-2 (10-20m) or Landsat (30m)
 *  2. Dual Temporal Modes:
 *     - 'ccdc'      : Harmonic synthetic evaluation at wet/dry seasonal fractions
 *     - 'composite' : Direct cloud-masked seasonal percentile/median composites
 *  3. Integrated Covariates:
 *     - Optical seasonal bands (wet & dry)
 *     - Multi-spectral indices & phenological deltas (wet - dry)
 *     - GLCM textural metrics (dry-season NIR)
 *     - SAR backscatter (Sentinel-1 seasonal VV/VH, ALOS PALSAR/ScanSAR)
 *     - Multi-scale topography (Elevation, slope, aspect, TPI landform, TWI)
 *     - CCDC harmonics (amplitude & phase for NIR, Red, SWIR1)
 *  4. Machine Learning & Uncertainty:
 *     - Stratified sample train/validation splitting
 *     - Optional Random Forest hyperparameter grid search
 *     - Multi-probability Random Forest classifier (soft class probabilities)
 *     - Rigorous test set evaluation (Accuracy, Log Loss, Confusion Matrix)
 *  5. Spatial Regularization & Post-Processing:
 *     - Empirical Bayes spatial smoothing (Camara et al., 2024 / sits)
 *     - GEOBIA SNIC superpixel clustering & adaptive stability blending
 *     - Intrinsic rule-based hydrology/urban refinement
 *  6. Export Pipeline:
 *     - Pre-configured pyramiding policies ('mean' for probs, 'mode' for LULC, 'sample' for segments)
 *     - Cloud-Optimized GeoTIFF (COG) Drive exports and GEE Asset exports
 * =========================================================================
 */

// =========================================================================
// 1. PROJECT CONFIGURATION & PARAMETERS (CUSTOMIZE HERE)
// =========================================================================

// Module repository path prefix (adjust to your GEE account/repository)
var MODULE_ROOT = 'users/rg2icraf/Luma:';

// Import core modules
var cfg         = require(MODULE_ROOT + 'config');
var feats       = require(MODULE_ROOT + 'features');
var clf         = require(MODULE_ROOT + 'classify_model');
var tuning      = require(MODULE_ROOT + 'tune_model');
var post        = require(MODULE_ROOT + 'Post-process');
var extraction  = require(MODULE_ROOT + 'extract_split_data');

// Sensor and input modules
var ccdcCore    = require(MODULE_ROOT + 'ccdc_core');
var s2Inputs    = require(MODULE_ROOT + 'inputs_s2');
var landsatInputs = require(MODULE_ROOT + 'inputs_landsat');
var compInputs  = require(MODULE_ROOT + 'inputs_composite');

// -------------------------------------------------------------------------
// Pipeline Execution Mode Settings
// -------------------------------------------------------------------------
var SENSOR        = 'S2';         // 'S2' (Sentinel-2, 20m) or 'Landsat' (30m)
var TEMPORAL_MODE = 'ccdc';       // 'ccdc' (harmonic reconstruction) or 'composite' (direct medians)
var REFERENCE_YEAR = 2025;        // Epoch corresponding to training ground truth
var TARGET_YEARS   = [2025, 2020];// Epochs to classify and map

// Execution Toggles
var RUN_TUNING               = false; // Run RF hyperparameter grid search
var TRAIN_MODEL              = true;  // Train Random Forest classifier
var APPLY_BAYESIAN_SMOOTH    = true;  // Apply empirical Bayes spatial smoothing
var APPLY_HYDROLOGY_RULES    = true;  // Refine false hydrology & urban roof confusion
var MERGE_WATER_RIVERBED     = true;  // Merge Dry riverbed (18) into Waterbody (17)
var EXPORT_PROBABILITIES     = true;  // Export uint16 probability cubes
var EXPORT_LULC_MAPS         = true;  // Export uint8 thematic classification
var EXPORT_TO_DRIVE          = true;  // Export to Google Drive (COG format)
var EXPORT_TO_ASSET          = false; // Export to GEE Assets
var EXPORT_TABLE_SAMPLES     = false; // Export extracted sample CSV for offline RFECV

// -------------------------------------------------------------------------
// Geographic & Data Assets
// -------------------------------------------------------------------------
// Area of Interest (pass null to use project default AOI from cfg)
var USER_AOI = null;
var aoi = cfg.getAOI(USER_AOI);

// Reference training data FeatureCollection
var SAMPLES_ASSET = 'projects/ee-rg2/assets/Samples/New_Samples_TL_2025';
var CLASS_PROP    = 'new_id'; // Feature property with integer class labels

// Coordinate Reference System and Resolution
var CRS   = cfg.DEFAULT_CRS;
var SCALE = (SENSOR === 'S2') ? cfg.SCALE.S2 : cfg.SCALE.Landsat;

// Destination folder for Drive exports
var DRIVE_FOLDER = 'LULC_Master_Exports';

// Class Schema and Bayesian Smoothness Weights
var CLASS_SCHEMA = cfg.CLASS_SCHEMA;
var SMOOTHNESS   = cfg.SMOOTHNESS;

// Optional: Pre-computed CCDC asset(s) or null to compute / load via cfg
var CCDC_ASSET = null; 

// Spectral indices to compute
var BASE_INDICES = ['BLFEI', 'DBSI', 'GBNDVI', 'MBI', 'MNDWI', 'MSAVI', 'NDMI', 'NDVI'];
var DELTA_INDICES = ['BLFEI', 'DBSI', 'MBI', 'MNDWI', 'MSAVI', 'NDMI'];

// Optional: Optimized feature subset from RFECV (set null to use all candidate features)
var FEATURE_SUBSET = (SENSOR === 'S2') ? [
  'B11_dry', 'B11_wet', 'B12_dry', 'B12_wet', 'B2_dry', 'B2_wet',
  'B3_dry', 'B3_wet', 'B4_dry', 'B4_wet', 'B8_dry', 'B8_wet',
  'BLFEI_dry', 'BLFEI_wet', 'DBSI_dry', 'GBNDVI_dry', 'GBNDVI_wet',
  'MBI_delta', 'MBI_dry', 'MBI_wet', 'MNDWI_dry', 'MNDWI_wet',
  'MSAVI_wet', 'NDMI_delta', 'NDMI_dry', 'NDVI_dry', 'NDVI_wet',
  'ScanSAR_HH_dB', 'ScanSAR_HV_dB', 'TWI', 'VH_dry', 'VH_wet',
  'VV_dry', 'VV_wet', 'elevation', 'red_amp', 'slope', 'swir1_amp'
] : null;

// Random Forest Hyperparameters
var RF_OPTIONS = {
  nTrees:  400,
  vSplit:  (FEATURE_SUBSET ? Math.round(Math.sqrt(FEATURE_SUBSET.length)) : null),
  minLeaf: 1,
  seed:    42
};

// =========================================================================
// 2. TEMPORAL DATA PREPARATION (CCDC VS COMPOSITE)
// =========================================================================
print('--- Master Pipeline Initialized ---');
print('Sensor:', SENSOR, '| Temporal Mode:', TEMPORAL_MODE);
print('Target Years:', TARGET_YEARS);

var ccdcData = null;
if (TEMPORAL_MODE === 'ccdc') {
  if (CCDC_ASSET) {
    ccdcData = ee.Image(CCDC_ASSET);
  } else if (SENSOR === 'S2') {
    ccdcData = cfg.loadCCDCTiles(cfg.S2_CCDC_TILES, aoi);
  } else {
    ccdcData = cfg.loadCCDCTiles(cfg.LANDSAT_CCDC_TILES, aoi);
  }
}

// =========================================================================
// 3. FEATURE STACK GENERATION FUNCTION
// =========================================================================
function buildEpochFeatureStack(year) {
  var wetOptical, dryOptical;
  var harmImgs = ee.Image();

  if (TEMPORAL_MODE === 'ccdc') {
    // Stage 3 CCDC Synthetic Reconstruction at calibrated seasonal fractions
    var wetDate = year + cfg.T.wet;
    var dryDate = year + cfg.T.dry;
    var harmDate = year + cfg.T.harm;

    var coefCase = (SENSOR === 'Landsat') ? 'upper' : 'lower';
    var outCase  = (SENSOR === 'Landsat') ? 'upper' : 'lower';

    wetOptical = ccdcCore.getSyntheticStack(ccdcData, wetDate, {
      coefScale: 10000,
      coefCase: coefCase,
      outputCase: outCase,
      tolYears: cfg.TOL_YEARS
    });

    dryOptical = ccdcCore.getSyntheticStack(ccdcData, dryDate, {
      coefScale: 10000,
      coefCase: coefCase,
      outputCase: outCase,
      tolYears: cfg.TOL_YEARS
    });

    harmImgs = ccdcCore.getHarmonics(ccdcData, harmDate, ['nir', 'red', 'swir1'], {
      coefCase: coefCase,
      outputCase: outCase,
      tolYears: cfg.TOL_YEARS
    });
  } else {
    // Direct Cloud-Masked Seasonal Composite
    wetOptical = compInputs.seasonalComposite(SENSOR, aoi, year, 'wet', {
      scaleReflectance: true,
      reducer: 'median'
    });
    dryOptical = compInputs.seasonalComposite(SENSOR, aoi, year, 'dry', {
      scaleReflectance: true,
      reducer: 'median'
    });
  }

  // 1. Rename Reflective Bands for Wet & Dry
  var bandMap = (SENSOR === 'S2') ? {
    src: ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
    wet: ['B2_wet', 'B3_wet', 'B4_wet', 'B8_wet', 'B11_wet', 'B12_wet'],
    dry: ['B2_dry', 'B3_dry', 'B4_dry', 'B8_dry', 'B11_dry', 'B12_dry']
  } : {
    src: ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
    wet: ['B2_wet', 'B3_wet', 'B4_wet', 'B5_wet', 'B6_wet', 'B7_wet'],
    dry: ['B2_dry', 'B3_dry', 'B4_dry', 'B5_dry', 'B6_dry', 'B7_dry']
  };

  var wetBands = wetOptical.select(bandMap.src, bandMap.wet);
  var dryBands = dryOptical.select(bandMap.src, bandMap.dry);

  // 2. GLCM Textural Features from dry-season NIR
  var glcmFeatures = feats.getGlcmTexture(dryOptical.select('nir'));

  // 3. Spectral Indices (wet and dry)
  var wetIndices = feats.computeSeasonalIndices(wetOptical, SENSOR, BASE_INDICES, 'wet');
  var dryIndices = feats.computeSeasonalIndices(dryOptical, SENSOR, BASE_INDICES, 'dry');

  // 4. Phenological Seasonal Deltas (wet - dry)
  var deltaIndices = feats.computePhenologicalDeltas(wetIndices, dryIndices, DELTA_INDICES);

  // 5. Radar Features (Sentinel-1 and/or PALSAR/ScanSAR)
  var radarFeatures = ee.Image();
  if (year >= 2015) {
    var s1Seasonal = feats.getSeasonalS1(year, aoi);
    var scanSar = feats.getScanSARFeatures(year, aoi);
    radarFeatures = ee.Image.cat([s1Seasonal, scanSar]);
  } else {
    var palsar = feats.getPalsarBackscatter(year, aoi);
    radarFeatures = palsar;
  }

  // 6. Topography Features (Time-invariant)
  var topoFeatures = feats.getTopography(aoi);

  // 7. Concatenate all candidate feature layers
  var stackList = [
    wetBands,
    dryBands,
    glcmFeatures,
    wetIndices,
    dryIndices,
    deltaIndices,
    radarFeatures,
    topoFeatures
  ];
  if (TEMPORAL_MODE === 'ccdc') {
    stackList.push(harmImgs);
  }

  var fullStack = ee.Image.cat(stackList).toFloat().clip(aoi);

  // Apply feature selection if defined
  if (FEATURE_SUBSET && FEATURE_SUBSET.length > 0) {
    return fullStack.select(FEATURE_SUBSET);
  }
  return fullStack;
}

// Generate stacks for all target epochs
var featureStacks = {};
TARGET_YEARS.forEach(function(yr) {
  print('Building feature stack for epoch:', yr);
  featureStacks[yr] = buildEpochFeatureStack(yr);
});

var refStack = featureStacks[REFERENCE_YEAR];
print('Active Predictor Bands (' + refStack.bandNames().size().getInfo() + '):', refStack.bandNames());

// =========================================================================
// 4. SAMPLE EXTRACTION & MODEL TRAINING
// =========================================================================
var groundTruth = ee.FeatureCollection(SAMPLES_ASSET);

// Export candidate samples table if requested for offline RFECV
if (EXPORT_TABLE_SAMPLES) {
  var samplesTable = refStack.sampleRegions({
    collection: groundTruth,
    properties: [CLASS_PROP],
    scale: SCALE,
    tileScale: 16,
    geometries: false
  });
  Export.table.toDrive({
    collection: samplesTable,
    description: 'Luma_Samples_' + SENSOR + '_' + REFERENCE_YEAR,
    folder: DRIVE_FOLDER,
    fileFormat: 'CSV'
  });
  print('Sample extraction task submitted to Drive.');
}

// Stratified train/test partition (70% train, 30% test)
var sampleSplit = extraction.stratifiedSplit(groundTruth, refStack, CLASS_PROP, SCALE, 0.7, 42);
var trainPixels = sampleSplit.trainingPixels;
var testPixels  = sampleSplit.testingPixels;

print('Training Samples:', trainPixels.size());
print('Testing Samples:',  testPixels.size());

// Optional Hyperparameter Tuning
if (RUN_TUNING) {
  print('Running Random Forest Hyperparameter Grid Search...');
  var treeList  = ee.List([100, 200, 400, 500]);
  var splitList = ee.List([6, 9, 12, 19]);
  var tuningRes = tuning.tuneRandomForest(
    trainPixels,
    testPixels,
    CLASS_PROP,
    refStack.bandNames(),
    treeList,
    splitList
  );
  var bestOpt = tuning.getBestParams(tuningRes, 'accuracy');
  print('Optimal Hyperparameters Found:', bestOpt);
}

// Train Multi-Probability Random Forest Model
var trainedClfModel = null;
var classBands = null;
var classValues = null;

if (TRAIN_MODEL) {
  print('Training Multi-Probability Random Forest Model...');
  var rfResult = clf.multiProbabilityClassification(
    trainPixels,
    CLASS_PROP,
    CLASS_SCHEMA,
    refStack,
    RF_OPTIONS
  );

  trainedClfModel = rfResult.trainedModel;
  classBands      = rfResult.classBands;
  classValues     = rfResult.classValues;

  // Evaluate on holdout test set
  var evalProb = clf.evaluateProbabilities(
    trainedClfModel,
    testPixels,
    CLASS_PROP,
    CLASS_SCHEMA
  );
  print('--- Validation Performance ---');
  print('Overall Accuracy:',  evalProb.get('overallAccuracy'));
  print('Overall Log Loss:',  evalProb.get('overallLogLoss'));
  print('Per-Class Accuracy:', evalProb.get('perClassAccuracy'));
  print('Per-Class Log Loss:', evalProb.get('perClassLogLoss'));
}

// =========================================================================
// 5. MULTI-YEAR INFERENCE & SPATIAL POST-PROCESSING
// =========================================================================
var finalResults = {};

TARGET_YEARS.forEach(function(yr) {
  var yrStack = featureStacks[yr];
  
  // 1. Raw multi-probability inference
  var probsArray = yrStack.classify(trainedClfModel);
  var rawProbs = probsArray.rename('probs').arrayFlatten([classBands]);

  // 2. Empirical Bayes spatial smoothing
  var smoothedProbs = rawProbs;
  if (APPLY_BAYESIAN_SMOOTH) {
    smoothedProbs = post.bayesianSmooth(rawProbs, classBands, SMOOTHNESS, 5, 0.5);
  }

  // 3. Crisp Thematic Classification from probabilities
  var lulcMap = post.classifyFromProbs(smoothedProbs, classBands, classValues);

  // 4. Intrinsic Rule Refinement (Hydrology vs Urban Roofs / Cleared Land)
  var stabilityInfo = post.computeStability(smoothedProbs, classBands, classValues);
  var refinedLulc = lulcMap;

  if (APPLY_HYDROLOGY_RULES) {
    refinedLulc = post.suppressFalseHydrology(
      lulcMap,
      stabilityInfo.runner_up,
      stabilityInfo.CI.multiply(100),
      {
        ciThreshold: 65,
        riverbedCiThreshold: 75,
        waterClasses: [17, 18],
        terrestrialClasses: [15, 16],
        riverbedId: 18,
        waterbodyId: 17,
        clearedLandId: 15,
        mergeRiverbed: MERGE_WATER_RIVERBED
      }
    );
  } else if (MERGE_WATER_RIVERBED) {
    refinedLulc = lulcMap.where(lulcMap.eq(18), 17);
  }

  finalResults[yr] = {
    probs: smoothedProbs,
    lulc: refinedLulc.toUint8(),
    stability: stabilityInfo
  };

  // Add map layers for reference year
  if (yr === REFERENCE_YEAR) {
    var vizLulc = {min: 1, max: 18, palette: cfg.CLASS_PALETTE};
    Map.centerObject(aoi, 9);
    Map.addLayer(refinedLulc.clip(aoi), vizLulc, 'Classified LULC ' + yr);
    Map.addLayer(stabilityInfo.is_stable.clip(aoi), {min: 0, max: 1, palette: ['#D7191C', '#1A9641']}, 'Stability Map ' + yr, false);
    Map.addLayer(stabilityInfo.CI.clip(aoi), {min: 0, max: 1, palette: ['#1A9641', '#FFFFBF', '#D7191C']}, 'Confusion Index ' + yr, false);
  }
});

// =========================================================================
// 6. EXPORT DISPATCHER (DRIVE & ASSET)
// =========================================================================
TARGET_YEARS.forEach(function(yr) {
  var res = finalResults[yr];
  var prefix = 'Luma_' + SENSOR + '_' + yr;

  // A. Export Soft Probabilities (uint16, 0-10000)
  if (EXPORT_PROBABILITIES) {
    var scaledProbs = res.probs.multiply(10000).round().toUint16().clip(aoi);

    if (EXPORT_TO_DRIVE) {
      cfg.toDrive(scaledProbs, prefix + '_probabilities', {
        folder: DRIVE_FOLDER,
        region: aoi,
        scale: SCALE,
        crs: CRS
      });
    }
    if (EXPORT_TO_ASSET) {
      cfg.toAsset(scaledProbs, 'projects/gee-v1-510606/assets/' + prefix + '_probabilities', {
        region: aoi,
        scale: SCALE,
        crs: CRS,
        pyramidingPolicy: { '.default': 'mean' }
      });
    }
  }

  // B. Export Crisp Discrete LULC Map (uint8)
  if (EXPORT_LULC_MAPS) {
    var discreteMap = res.lulc.clip(aoi);

    if (EXPORT_TO_DRIVE) {
      cfg.toDrive(discreteMap, prefix + '_LULC_thematic', {
        folder: DRIVE_FOLDER,
        region: aoi,
        scale: SCALE,
        crs: CRS
      });
    }
    if (EXPORT_TO_ASSET) {
      cfg.toAsset(discreteMap, 'projects/gee-v1-510606/assets/' + prefix + '_LULC_thematic', {
        region: aoi,
        scale: SCALE,
        crs: CRS,
        pyramidingPolicy: { '.default': 'mode' }
      });
    }
  }
});

print('Master Pipeline execution complete. Check Tasks tab to run exports.');
