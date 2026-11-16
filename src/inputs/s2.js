/**
 * Module: Sentinel-2 Input Data Preparation
 * =========================================================================
 * Provides cloud-masked, harmonized surface reflectance collections from
 * Sentinel-2 L2A (Harmonized) linked with Google Cloud Score+.
 *
 * Mapped to the unified 6-band lowercase schema:
 *   ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']
 *
 * Configurable reflectance scale:
 *   - reflectanceScale = 10000 (default, native S2 integer range for CCDC)
 *   - reflectanceScale = 1.0   (direct surface reflectance in 0..1 range)
 * =========================================================================
 */

/**
 * Prepares Sentinel-2 L2A Harmonized Surface Reflectance time series linked with Cloud Score+.
 *
 * @param {ee.Geometry} roi                      - Region of interest boundary
 * @param {string}      [startDate='2019-01-01'] - Start date (YYYY-MM-DD)
 * @param {string}      [endDate='2025-12-31']   - End date (YYYY-MM-DD)
 * @param {Object}      [options]                - Optional configurations:
 *   @param {number}      [options.clearThreshold=0.60]   - Cloud Score+ threshold (0.50 - 0.65)
 *   @param {string}      [options.qaBand='cs_cdf']       - Quality band ('cs_cdf' or 'cs')
 *   @param {number}      [options.reflectanceScale=10000]- Scale multiplier (10000 for CCDC, 1 for 0..1 float)
 * @returns {ee.ImageCollection} Harmonized 6-band collection with 'system:time_start'
 */
exports.getCollection = function(roi, startDate, endDate, options) {
  startDate = startDate || '2019-01-01';
  endDate   = endDate   || '2025-12-31';
  options   = options   || {};

  var clearThreshold = options.clearThreshold !== undefined ? options.clearThreshold : 0.60;
  var qaBand         = options.qaBand || 'cs_cdf';
  var scaleFactor    = options.reflectanceScale !== undefined ? options.reflectanceScale : 10000;

  var inputBands  = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'];
  var outputNames = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'];

  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(roi)
    .filterDate(startDate, endDate);

  var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED')
    .filterBounds(roi)
    .filterDate(startDate, endDate);

  var s2Linked = s2.linkCollection(csPlus, [qaBand]);

  // S2 SR Harmonized raw integers are in range 0..10000
  // When scaleFactor == 10000, keep native scale (divide by 1).
  // When scaleFactor == 1.0, divide by 10000 to get 0..1 float.
  var multiplier = scaleFactor / 10000.0;

  return s2Linked.map(function(img) {
    var cloudMask = img.select(qaBand).gte(clearThreshold);
    var scaled = img.select(inputBands, outputNames);

    if (multiplier !== 1.0) {
      scaled = scaled.multiply(multiplier);
    }

    return scaled.updateMask(cloudMask)
      .copyProperties(img, ['system:time_start']);
  });
};
