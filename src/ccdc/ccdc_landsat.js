/**
 * Module: Landsat CCDC Pipeline (Compatibility Facade)
 * =========================================================================
 * Provides backward compatibility for all existing scripts requiring:
 *   require('users/rg2icraf/Luma:ccdc_landsat')
 *
 * Internally delegates to the stage-based architecture:
 *   - src/inputs/landsat.js
 *   - src/inputs/composite.js
 *   - src/ccdc/core.js
 * =========================================================================
 */

var inputsLandsat = require('users/rg2icraf/Luma:inputs/landsat');
var composite     = require('users/rg2icraf/Luma:inputs/composite');
var core          = require('users/rg2icraf/Luma:ccdc/core');

// 1. Data Ingestion
exports.getLandsatTimeSeries = function(roi, startDate, endDate, options) {
  return inputsLandsat.getCollection(roi, startDate, endDate, options);
};

// 2. Compositing
exports.toMonthlyMedian = function(col, startDate, endDate) {
  return composite.toMonthlyMedian(col, startDate, endDate);
};

// 3. Execution
exports.runLandsatCCDC = function(lsCollection, options) {
  options = options || {};
  options.breakpointBands = options.breakpointBands || ['green', 'red', 'nir', 'swir1', 'swir2'];
  options.tmaskBands = options.tmaskBands || ['green', 'swir1'];
  return core.runCCDC(lsCollection, options);
};

// 4. Nearest Segment Lookup
exports.nearestSegment = function(ccdc, t) {
  return core.nearestSegment(ccdc, t);
};

// 5. Hole Filling
exports.fillIsolated = function(img, minNeighbors) {
  return core.fillIsolated(img, minNeighbors);
};

// 6. Synthetic Reflectance Stack
exports.getSyntheticLandsatStack = function(ccdc, t, addIndices, tolYears, coefScale) {
  addIndices = addIndices !== undefined ? addIndices : true;
  tolYears  = tolYears !== undefined ? tolYears : 1.0;
  coefScale = coefScale !== undefined ? coefScale : 10000;

  var stack = core.getSyntheticStack(ccdc, t, {
    bands:          ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
    coefScale:      coefScale,
    tolYears:       tolYears,
    coefCase:       'upper', // Legacy Landsat assets use uppercase (e.g. BLUE_coefs)
    minReflectance: -0.05,
    maxReflectance: 1.0,
    fillHoles:      true
  });

  // Strip diagnostic bands for exact legacy band signature if desired, or keep
  var opticalBands = stack.select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']);

  if (addIndices) {
    var ndvi = opticalBands.normalizedDifference(['nir', 'red']).rename('NDVI');
    var ndwi = opticalBands.normalizedDifference(['green', 'nir']).rename('NDWI');
    var ndbi = opticalBands.normalizedDifference(['swir1', 'nir']).rename('NDBI');
    opticalBands = opticalBands.addBands([ndvi, ndwi, ndbi]);
  }

  return opticalBands.addBands(stack.select('gap_years'));
};

// 7. Harmonics
exports.getHarmonics = function(ccdc, t, bands, tolYears, coefScale) {
  tolYears  = tolYears !== undefined ? tolYears : 1.0;
  coefScale = coefScale !== undefined ? coefScale : 10000;

  return core.getHarmonics(ccdc, t, {
    bands:      bands || ['nir', 'red', 'swir1'],
    coefScale:  coefScale,
    tolYears:   tolYears,
    coefCase:   'upper', // Legacy Landsat assets use uppercase
    outputCase: 'upper'  // Legacy Landsat scripts expect 'NIR_amp', 'RED_amp'
  });
};
