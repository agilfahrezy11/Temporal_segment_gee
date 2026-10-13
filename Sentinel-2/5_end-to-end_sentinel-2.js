/**
 * =========================================================================
 * Script: 5_End-to-End_Sentinel-2.js
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

/**
 * Setup 
 * 
*/
// Area of Interest (Timor-Leste + 10km buffer)
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry();
//Import required modules
var spectral    = require('users/dmlmont/spectral:spectral');
var ccdc_module = require('users/rg2icraf/Luma:ccdc_sentinel2');
var extraction  = require('users/rg2icraf/Luma:extract_split_data');
var clf         = require('users/rg2icraf/Luma:classify_model');
var tuning      = require('users/rg2icraf/Luma:tune_model');
var post        = require('users/rg2icraf/Luma:post_process');

/**
 * Configurations
 * 
*/

// Seasonal fraction settings
var T = {
  dry:  0.65,   // ~August (Peak dry season)
  wet:  0.16,   // ~March (Peak wet season)
  harm: 0.50    // Mid-year evaluation for nearest segment harmonics
};

// Base indices to calculate from spectral library
var baseIndexList = ['BLFEI', 'DBSI', 'GBNDVI', 'MBI', 'MNDWI', 'MSAVI', 'NDMI', 'NDVI', 'MCARI'];
// Indices requiring delta computation (wet - dry)
var deltaIndexList = ['BLFEI', 'DBSI', 'MBI', 'MNDWI', 'MSAVI', 'NDMI'];

// Target years to process
var targetYears = [2025, 2015, 2005];

// Class schema (18 classes)
var classIdToName = {
  1:  'Primary_forest',
  2:  'Secondary_forest',
  3:  'Mangrove_forest',
  4:  'Coastal_forest',
  5:  'Teak_plantation',
  6:  'Eucalyptus_plantation',
  7:  'Coffe_agroforestry',
  8:  'Cocoa_agroforestry',
  9:  'Coconut_agroforestry',
  10: 'Mixed_garden',
  11: 'Paddy_field',
  12: 'Other_cropland',
  13: 'Grassland',
  14: 'Shrubland',
  15: 'Cleared_land',
  16: 'Built-up',
  17: 'Waterbody',
  18: 'Dry_riverbed'
};

// Bayesian smoothing parameters per class
var smoothnessTL = {
  'Primary_forest':        4.5,
  'Secondary_forest':      3.0,
  'Mangrove_forest':       3.0,
  'Coastal_forest':        3.5,
  'Teak_plantation':       1.2,
  'Eucalyptus_plantation': 1.2,
  'Coffe_agroforestry':    1.5,
  'Cocoa_agroforestry':    1.5,
  'Coconut_agroforestry':  1.5,
  'Mixed_garden':          1.2,
  'Paddy_field':           1.0,
  'Other_cropland':        1.0,
  'Grassland':             1.5,
  'Shrubland':             1.5,
  'Cleared_land':          1.5,
  'Built-up':              1,
  'Waterbody':             4.0,
  'Dry_riverbed':          0.1
};

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

//define CCDC assets
var tl_ccdc1 = ee.Image('projects/ee-rg2v2/assets/CCDC_TL_2020_2025_20m_westbox')
var tl_ccdc2 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_east-timor')
var tl_ccdc3 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_off_dili')
var tl_ccdc4 = ee.Image('projects/gee-v1-510606/assets/CCDC_TL_2020_2025_20m_central-timor')
var tl_ccdc5 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_east_box')

//combine them into a single mosaic
var ccdc_mosaic = ee.ImageCollection([tl_ccdc1, tl_ccdc2, tl_ccdc3, tl_ccdc4, tl_ccdc5]).mosaic().clip(aoi);
print(ccdc_mosaic)

/**
 * Build Feature set for the target year
 * 
*/

/**
 * Computes required spectral indices on an already reconstructed synthetic image.
 */
function computeIndicesFromSynth(synthImg, suffix) {
  var imgForIndices = synthImg.select(
    ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
    ['B2',   'B3',    'B4',  'B8',  'B11',    'B12']
  );
  var param = {
    'B':  imgForIndices.select('B2'),
    'G':  imgForIndices.select('B3'),
    'R':  imgForIndices.select('B4'),
    'N':  imgForIndices.select('B8'),
    'S1': imgForIndices.select('B11'),
    'S2': imgForIndices.select('B12'),
    'C1': 6.0,
    'C2': 7.5,
    'g':  2.5,
    'L':  1.0
  };
  var computed = spectral.computeIndex(imgForIndices, baseIndexList, param);
  var renamed = computed.bandNames().map(function(b) {
    return ee.String(b).cat('_').cat(suffix);
  });
  return computed.rename(renamed);
}

function computeHarmonicsFromActiveSegment(ccdc, t, bands) {
  var time = ee.Number(t);
  var segmentMask = ccdc.select('tStart').lte(time)
    .and(ccdc.select('tEnd').gte(time))
    .toArray(0)
    .toArray(1);

  var harmonics = bands.map(function(band) {
    var coefficients = ccdc.select(band + '_coefs')
      .arrayMask(segmentMask)
      .arraySlice(0, 0, 1);
    var a1 = coefficients.arraySlice(1, 2, 3)
      .arrayProject([0])
      .arrayFlatten([['a1']])
      .divide(10000);
    var b1 = coefficients.arraySlice(1, 3, 4)
      .arrayProject([0])
      .arrayFlatten([['b1']])
      .divide(10000);

    return ee.Image.cat([
      a1.hypot(b1).rename(band + '_amp'),
      b1.atan2(a1).rename(band + '_phase')
    ]);
  });

  return ee.Image.cat(harmonics);
}

//function to build the feature stack for a given year using CCDC synthetic imagery
//use spectral index
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

  // 3. GLCM texture from dry-season NIR
  var glcm = synthDry.select('nir', ['B8'])
    .multiply(1000)
    .toInt32()
    .glcmTexture({size: 1});
  var glcmFeatures = glcm.select(['B8_contrast','B8_ent','B8_asm','B8_idm']);

  // 4. Compute spectral indices directly from existing synthetic stacks
  var wetIndices = computeIndicesFromSynth(synthWet, 'wet');
  var dryIndices = computeIndicesFromSynth(synthDry, 'dry');

  // 5. Compute phenological seasonal deltas (wet - dry)
  var deltaImages = deltaIndexList.map(function(idx) {
    return wetIndices.select(idx + '_wet')
      .subtract(dryIndices.select(idx + '_dry'))
      .rename(idx + '_delta');
  });
  var deltaIndices = ee.Image.cat(deltaImages);

  // 6. CCDC Harmonics (NIR, RED, SWIR1)
  var harmImgs = computeHarmonicsFromActiveSegment(
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
function getSeasonalS1(year, roi) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = ee.Date.fromYMD(year, 12, 31);

  var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .select(['VV', 'VH']);

  //Monsoonal Splitting
  var wet = s1.filter(ee.Filter.calendarRange(12, 4, 'month')).median();
  var dry = s1.filter(ee.Filter.calendarRange(5, 11, 'month')).median();
  //Rename base bands
  var vv_wet = wet.select('VV').rename('VV_wet');
  var vh_wet = wet.select('VH').rename('VH_wet');

  var vv_dry = dry.select('VV').rename('VV_dry');
  var vh_dry = dry.select('VH').rename('VH_dry');
  return ee.Image.cat([vh_dry, vh_wet, vv_dry, vv_wet]);
} 
//implement the code for radar backscatter features
var s1_seasonal_2020 = getSeasonalS1(2020, aoi).clip(aoi);
var s1_seasonal_2025 = getSeasonalS1(2025, aoi).clip(aoi);
// Build feature stacks for each epoch from the unified Landsat CCDC mosaic
var spectral_features_2025 = buildYearFeatures(ccdc_mosaic, 2025);
var spectral_features_2020 = buildYearFeatures(ccdc_mosaic, 2020);

//Topography
var elev = ee.Image("NASA/NASADEM_HGT/001")
            .select('elevation')
            .clip(aoi)
var terrain_product = ee.Terrain.products(elev).select('slope', 'aspect').toFloat();
//Landform 
var landform = ee.Image('projects/ee-rg2v2/assets/TPI_landform_TL')
                      .rename('landform')
                      .toByte()
                      .clip(aoi)
//Topographic Wetness Index
var twi = ee.Image('projects/ee-rg2/assets/Top_Wetness_index')
                  .rename('TWI')
                  .toFloat()

//Stack topographical data
var topo = elev.int16().addBands(terrain_product).addBands(landform).addBands(twi)//(2. Topographical Data Features)

//Stack all features for each year
var final_stack_2025 = ee.Image.cat([spectral_features_2025, s1_seasonal_2025, topo]).clip(aoi);
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
