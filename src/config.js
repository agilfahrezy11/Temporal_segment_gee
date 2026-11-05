/**
 * Module: Project Configuration
 * =========================================================================
 * Single source of truth for all project-level constants:
 *   - AOI, CRS, and export defaults
 *   - Seasonal fraction settings
 *   - 18-class land cover schema
 *   - Bayesian smoothing weights (canonical values)
 *   - CCDC tile asset IDs (Landsat and Sentinel-2)
 *   - Shared export helpers (toAsset / toDrive)
 *
 * Usage in any pipeline script:
 *   var cfg = require('users/rg2icraf/Luma:config');
 *   var aoi = ee.FeatureCollection(cfg.AOI_ASSET).geometry();
 *   cfg.toAsset(myImage, 'projects/.../assets/my_output', { scale: 20, region: aoi });
 * =========================================================================
 */

// ---------------------------------------------------------------------------
// 1. Study Area / AOI Configuration
// ---------------------------------------------------------------------------

/** Default/Example GEE asset path for AOI (Timor-Leste 10km buffer). */
exports.DEFAULT_AOI_ASSET = 'projects/ee-rg2v2/assets/AOI_Timor_Leste_10km';
exports.AOI_ASSET = exports.DEFAULT_AOI_ASSET;

/**
 * Helper to load an AOI geometry from an asset path or feature collection.
 * Defaults to exports.AOI_ASSET if no path or geometry is passed.
 *
 * @param {string|ee.Geometry|ee.FeatureCollection} [aoiInput]
 * @returns {ee.Geometry}
 */
exports.getAOI = function(aoiInput) {
  if (!aoiInput) {
    return ee.FeatureCollection(exports.AOI_ASSET).geometry();
  }
  if (typeof aoiInput === 'string') {
    return ee.FeatureCollection(aoiInput).geometry();
  }
  if (aoiInput instanceof ee.FeatureCollection) {
    return aoiInput.geometry();
  }
  return aoiInput; // Assumed ee.Geometry
};

// ---------------------------------------------------------------------------
// 2. Export Defaults
// ---------------------------------------------------------------------------

/** Default Coordinate reference system. Override via options.crs in toAsset/toDrive if needed. */
exports.DEFAULT_CRS = 'EPSG:32751';
exports.CRS = exports.DEFAULT_CRS;

/** Maximum pixel count passed to all Export calls. */
exports.MAX_PIXELS = 1e13;

/**
 * Default pixel scale in metres for each sensor.
 * Pass the appropriate value to toAsset() / toDrive() via options.scale.
 */
exports.SCALE = {
  S2:      20,   // Sentinel-2 (standard processing scale)
  S2_fine: 10,   // Sentinel-2 native fine resolution
  Landsat: 30    // Landsat 30m resolution
};

// ---------------------------------------------------------------------------
// 3. Seasonal Fractions (fractional year format: YYYY.fff)
// ---------------------------------------------------------------------------

/**
 * Fractional-year offsets added to a target year integer to produce the
 * evaluation date passed to getSyntheticStack() / getSyntheticLandsatStack().
 *
 *   e.g.  2025 + T.dry  ==  2025.65  ==  ~August (peak dry season)
 */
exports.T = {
  dry:  0.65,   // ~August  — peak dry season
  wet:  0.16,   // ~March   — peak wet season
  harm: 0.50    // ~July 1  — mid-year, used for harmonic extraction
};

/** Default maximum gap (years) to the nearest CCDC segment. */
exports.TOL_YEARS = 1.0;

// ---------------------------------------------------------------------------
// 4. Land Cover Class Schema (18 classes)
// ---------------------------------------------------------------------------

/**
 * Canonical mapping of integer class ID -> class name string.
 * Used by classify_model.multiProbabilityClassification() and
 * Post-process.classifyFromProbs().
 *
 * Do NOT redefine this in pipeline scripts — import it here.
 */
exports.CLASS_SCHEMA = {
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

/**
 * Ordered list of class names (ascending by ID), matching the band order
 * produced by ee.Classifier.smileRandomForest in MULTIPROBABILITY mode.
 * Use this as the classBands argument wherever an Array is needed.
 */
exports.CLASS_BANDS = (function() {
  var schema = exports.CLASS_SCHEMA;
  return Object.keys(schema)
    .map(Number)
    .sort(function(a, b) { return a - b; })
    .map(function(id) { return schema[id]; });
})();

/**
 * Ordered list of integer class IDs (ascending), matching CLASS_BANDS order.
 * Use this as the classValues argument for post.classifyFromProbs().
 */
exports.CLASS_VALUES = (function() {
  return Object.keys(exports.CLASS_SCHEMA)
    .map(Number)
    .sort(function(a, b) { return a - b; });
})();

/**
 * Canonical colour palette for map visualisation, indexed to CLASS_VALUES.
 * Use with: { min: 1, max: 18, palette: cfg.CLASS_PALETTE }
 */
exports.CLASS_PALETTE = [
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
  '#9ECAE1'  // 18 Dry_riverbed
];

// ---------------------------------------------------------------------------
// 5. Bayesian Smoothing Weights
// ---------------------------------------------------------------------------

/**
 * Canonical per-class smoothness values for Post-process.bayesianSmooth().
 *
 * Higher values = stronger pull toward the spatial neighbourhood mean
 * (good for large homogeneous classes: Primary_forest, Waterbody).
 * Lower values = trust the pixel-level RF probability more
 * (good for small or spectrally ambiguous classes: Dry_riverbed, Built-up).
 *
 * Reconciliation note — four versions of this dict existed across pipeline
 * scripts with diverged values. These are the calibrated values from
 * Sentinel-2/5_end-to-end_sentinel-2.js (most recent, highest-resolution run).
 * The main differences from older scripts:
 *   Cleared_land:     0.5  (was 2.0 in Landsat/4_Landsat_classification_assembly.js)
 *   Built-up:         0.2  (was 1.0–2.0 in older scripts)
 *   Secondary_forest: 3.0  (was 4.0 in Landsat/4_Landsat_classification_assembly.js)
 *   Waterbody:        4.0  (was 1.5–2.0 in older scripts)
 * If you are applying to Landsat outputs and need different weights,
 * override individual keys in the calling script rather than copying the dict.
 */
exports.SMOOTHNESS = {
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
  'Cleared_land':          0.5,
  'Built-up':              0.2,
  'Waterbody':             4.0,
  'Dry_riverbed':          0.1
};

// ---------------------------------------------------------------------------
// 6. CCDC Tile Asset IDs
// ---------------------------------------------------------------------------

/**
 * Five-tile Landsat CCDC mosaic (2000–2025, 30 m).
 * Load and merge with: cfg.loadCCDCTiles(cfg.LANDSAT_CCDC_TILES, aoi)
 * — or directly via features.loadCCDCMosaic() once src/features.js exists.
 */
exports.LANDSAT_CCDC_TILES = [
  'projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_central',
  'projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_east_timor',
  'projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_eastbox',
  'projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_off_dili',
  'projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_westbox'
];

/**
 * Five-tile Sentinel-2 CCDC mosaic (2020–2025, 20 m).
 * Spans three GEE project namespaces — do not edit unless re-exporting tiles.
 */
exports.S2_CCDC_TILES = [
  'projects/ee-rg2v2/assets/CCDC_TL_2020_2025_20m_westbox',
  'projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_east-timor',
  'projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_off_dili',
  'projects/gee-v1-510606/assets/CCDC_TL_2020_2025_20m_central-timor',
  'projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_east_box'
];

/**
 * Five-tile SNIC superpixel segment assets (Sentinel-2, 10 m).
 * Required by the hybrid pixel-object assembly script.
 */
exports.S2_SNIC_TILES = [
  'projects/gee-v1-510606/assets/snic/seg_central_timor',
  'projects/gee-v1-510606/assets/snic/seg_east_box',
  'projects/gee-v1-510606/assets/snic/seg_off_dili',
  'projects/gee-v1-510606/assets/snic/seg_west_box',
  'projects/gee-v1-510606/assets/snic/seg_east_timor'
];

// ---------------------------------------------------------------------------
// 7. CCDC Coefficient Scale
// ---------------------------------------------------------------------------

/**
 * Coefficient scale factor for CCDC outputs.
 *   10000 — own CCDC run (this project's Landsat and S2 pipelines)
 *   1     — GOOGLE/GLOBAL_CCDC/V1 (coefficients already in 0–1 reflectance)
 */
exports.CCDC_COEF_SCALE = {
  own:    10000,
  global: 1
};

// ---------------------------------------------------------------------------
// 8. Helper: Load and mosaic CCDC tiles
// ---------------------------------------------------------------------------

/**
 * Loads a list of GEE asset IDs as individual images, mosaics them, and
 * clips to the provided geometry.
 *
 * @param {Array<string>} tileIds - Array of GEE asset IDs
 * @param {ee.Geometry}   aoi     - Clip geometry
 * @returns {ee.Image} Seamless mosaicked image clipped to AOI
 *
 * Example:
 *   var ccdc = cfg.loadCCDCTiles(cfg.S2_CCDC_TILES, aoi);
 */
exports.loadCCDCTiles = function(tileIds, aoi) {
  var images = tileIds.map(function(id) { return ee.Image(id); });
  return ee.ImageCollection(images).mosaic().clip(aoi);
};

// ---------------------------------------------------------------------------
// 9. Export Helpers
// ---------------------------------------------------------------------------

/**
 * Registers an Export.image.toAsset task with canonical project defaults.
 * Reduces every Export block to ~5 lines instead of 8–10.
 *
 * @param {ee.Image}  image    - Image to export
 * @param {string}    assetId  - Full destination GEE asset path
 * @param {Object}    options  - Overrides: { description, region, scale,
 *                               pyramidingPolicy }
 *
 * pyramidingPolicy defaults:
 *   'mean'   — probability images and continuous rasters (DEFAULT)
 *   'sample' — CCDC array images
 *   'mode'   — categorical / label images
 *
 * Example:
 *   cfg.toAsset(probImage.multiply(10000).round().toUint16(),
 *     'projects/ee-rg2/assets/LC_prob_2025', {
 *       description:      'LC_prob_2025',
 *       region:           aoi,
 *       scale:            cfg.SCALE.S2,
 *       pyramidingPolicy: {'.default': 'mean'}
 *     });
 */
exports.toAsset = function(image, assetId, options) {
  options = options || {};
  Export.image.toAsset({
    image:            image,
    description:      options.description || assetId.split('/').pop(),
    assetId:          assetId,
    region:           options.region,
    scale:            options.scale  || exports.SCALE.S2,
    crs:              options.crs    || exports.CRS,
    maxPixels:        options.maxPixels || exports.MAX_PIXELS,
    pyramidingPolicy: options.pyramidingPolicy || {'.default': 'mean'}
  });
};

/**
 * Registers an Export.image.toDrive task as a Cloud-Optimized GeoTIFF
 * with canonical project defaults.
 *
 * @param {ee.Image}  image      - Image to export
 * @param {string}    fileName   - Output file name prefix (no extension)
 * @param {Object}    options    - Overrides: { description, folder, region,
 *                                 scale, crs, maxPixels }
 *
 * Example:
 *   cfg.toDrive(finalLulc, 'TL_LULC_2025_S2', {
 *     folder: 'TimorLeste_LULC',
 *     region: aoi,
 *     scale:  cfg.SCALE.S2
 *   });
 */
exports.toDrive = function(image, fileName, options) {
  options = options || {};
  Export.image.toDrive({
    image:          image,
    description:    options.description || fileName,
    folder:         options.folder      || 'GEE_Exports',
    fileNamePrefix: fileName,
    region:         options.region,
    scale:          options.scale || exports.SCALE.S2,
    crs:            options.crs   || exports.CRS,
    maxPixels:      options.maxPixels || exports.MAX_PIXELS,
    formatOptions:  { cloudOptimized: true }
  });
};
