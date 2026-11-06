/**
 * Module: Landsat Input Data Preparation
 * =========================================================================
 * Provides cloud-masked, harmonized surface reflectance collections across
 * Landsat 5 (TM), Landsat 7 (ETM+), Landsat 8 (OLI), and Landsat 9 (OLI-2).
 *
 * All sensors are mapped to a unified 6-band lowercase schema:
 *   ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']
 *
 * Configurable reflectance scale:
 *   - reflectanceScale = 10000 (default, calibrated for CCDC LASSO lambda = 20)
 *   - reflectanceScale = 1.0   (direct surface reflectance in 0..1 range)
 * =========================================================================
 */

/**
 * Prepares a merged Landsat 5, 7, 8, and 9 Collection 2 Level-2 Surface Reflectance time series.
 *
 * @param {ee.Geometry} roi                    - Region of interest boundary
 * @param {string}      [startDate='2000-01-01'] - Start date (YYYY-MM-DD)
 * @param {string}      [endDate='2025-12-31']   - End date (YYYY-MM-DD)
 * @param {Object}      [options]              - Optional configurations:
 *   @param {number}      [options.reflectanceScale=10000] - Scale multiplier (10000 for CCDC, 1 for 0..1 float)
 *   @param {boolean}     [options.maskFill=true]          - Mask fill/SLC-off pixels (bit 0)
 *   @param {boolean}     [options.maskCloud=true]         - Mask cloud (bit 3) and dilated cloud (bit 1)
 *   @param {boolean}     [options.maskShadow=true]        - Mask cloud shadow (bit 4)
 *   @param {boolean}     [options.maskSnow=true]          - Mask snow (bit 5)
 * @returns {ee.ImageCollection} Harmonized 6-band collection with 'system:time_start'
 */
exports.getCollection = function(roi, startDate, endDate, options) {
  startDate = startDate || '2000-01-01';
  endDate   = endDate   || '2025-12-31';
  options   = options   || {};

  var scaleFactor = options.reflectanceScale !== undefined ? options.reflectanceScale : 10000;
  var maskFill    = options.maskFill !== undefined ? options.maskFill : true;
  var maskCloud   = options.maskCloud !== undefined ? options.maskCloud : true;
  var maskShadow  = options.maskShadow !== undefined ? options.maskShadow : true;
  var maskSnow    = options.maskSnow !== undefined ? options.maskSnow : true;

  var maskAndScale = function(img, sourceBands) {
    var qa = img.select('QA_PIXEL');
    var mask = ee.Image(1);

    if (maskFill) {
      mask = mask.and(qa.bitwiseAnd(1).eq(0)); // Bit 0: Fill / SLC-off
    }
    if (maskCloud) {
      mask = mask.and(qa.bitwiseAnd(1 << 1).eq(0))  // Bit 1: Dilated cloud
                 .and(qa.bitwiseAnd(1 << 3).eq(0)); // Bit 3: Cloud
    }
    if (maskShadow) {
      mask = mask.and(qa.bitwiseAnd(1 << 4).eq(0)); // Bit 4: Cloud shadow
    }
    if (maskSnow) {
      mask = mask.and(qa.bitwiseAnd(1 << 5).eq(0)); // Bit 5: Snow
    }

    // USGS Collection 2 Level-2 Surface Reflectance scaling:
    // Reflectance = (0.0000275 * DN - 0.2)
    // Scaled by scaleFactor (10000 -> integer range ~0..10000; 1.0 -> 0..1 float)
    var srBands = img.select(
      sourceBands,
      ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']
    ).multiply(0.0000275).subtract(0.2).multiply(scaleFactor);

    return srBands.updateMask(mask)
      .copyProperties(img, ['system:time_start']);
  };

  // Landsat 5/7 (TM/ETM+) vs Landsat 8/9 (OLI/OLI-2) band assignments
  var tmEtmBands = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'];
  var oliBands   = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'];

  var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .map(function(img) { return maskAndScale(img, tmEtmBands); });

  var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .map(function(img) { return maskAndScale(img, tmEtmBands); });

  var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .map(function(img) { return maskAndScale(img, oliBands); });

  var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .map(function(img) { return maskAndScale(img, oliBands); });

  return l5.merge(l7).merge(l8).merge(l9);
};
