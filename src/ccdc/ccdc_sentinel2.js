/**
 * Module: Sentinel-2 CCDC Pipeline (Compatibility Facade)
 * =========================================================================
 * Provides backward compatibility for all existing scripts requiring:
 *   require('users/rg2icraf/Luma:ccdc_sentinel2')
 *
 * Internally delegates to the stage-based architecture:
 *   - src/inputs/s2.js
 *   - src/inputs/composite.js
 *   - src/ccdc/core.js
 * =========================================================================
 */

var inputsS2  = require('users/rg2icraf/Luma:inputs/s2');
var composite = require('users/rg2icraf/Luma:inputs/composite');
var core      = require('users/rg2icraf/Luma:ccdc/core');

// 1. Cloud Masking & Scaling
exports.maskAndScaleS2 = function(image, clearThreshold, qaBand) {
  clearThreshold = clearThreshold !== undefined ? clearThreshold : 0.60;
  qaBand = qaBand || 'cs_cdf';

  var cloudMask = image.select(qaBand).gte(clearThreshold);
  var inputBands  = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'];
  var outputBands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'];

  return image.select(inputBands, outputBands)
    .divide(10000)
    .updateMask(cloudMask);
};

// 2. Data Ingestion Collection
exports.getS2TimeSeries = function(roi, startDate, endDate, clearThreshold, qaBand) {
  return inputsS2.getCollection(roi, startDate, endDate, {
    clearThreshold:   clearThreshold,
    qaBand:           qaBand,
    reflectanceScale: 10000 // Native 0..10000 integer scale for CCDC
  });
};

// 3. Monthly Compositing
exports.toMonthlyMedian = function(col, startDate, endDate) {
  return composite.toMonthlyMedian(col, startDate, endDate);
};

// 4. Execution
exports.runCCDC = function(s2Collection, options) {
  options = options || {};
  options.breakpointBands = options.breakpointBands || ['green', 'red', 'nir', 'swir1'];
  return core.runCCDC(s2Collection, options);
};

// 5. Nearest Segment
exports.nearestSegment = function(ccdc, t) {
  return core.nearestSegment(ccdc, t);
};

// 6. Hole Filling
exports.fillIsolated = function(img, minNeighbors) {
  return core.fillIsolated(img, minNeighbors);
};

// 7. Synthetic Stack
exports.getSyntheticStack = function(ccdcOutput, targetDateFractional, bands, addIndices, coefficientScale, tolYears) {
  // Support flexible argument order (both S2 style and Landsat style)
  var _bands, _addIndices, _coefScale, _tolYears;
  if (typeof bands === 'boolean') {
    _addIndices = bands;
    _tolYears   = addIndices !== undefined ? addIndices : 1.0;
    _coefScale  = coefficientScale !== undefined ? coefficientScale : 10000;
    _bands      = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'];
  } else {
    _bands      = bands || ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'];
    _addIndices = addIndices !== undefined ? addIndices : true;
    _coefScale  = coefficientScale !== undefined ? coefficientScale : 10000;
    _tolYears   = tolYears !== undefined ? tolYears : 1.0;
  }

  var stack = core.getSyntheticStack(ccdcOutput, targetDateFractional, {
    bands:          _bands,
    coefScale:      _coefScale,
    tolYears:       _tolYears,
    coefCase:       'lower', // Sentinel-2 assets use lowercase (blue_coefs)
    minReflectance: -0.05,
    maxReflectance: 1.0,     // Unified 1.0 ceiling
    fillHoles:      true
  });

  var opticalBands = stack.select(_bands);

  if (_addIndices) {
    var ndvi = opticalBands.normalizedDifference(['nir', 'red']).rename('NDVI');
    var ndwi = opticalBands.normalizedDifference(['green', 'nir']).rename('NDWI');
    var ndbi = opticalBands.normalizedDifference(['swir1', 'nir']).rename('NDBI');
    opticalBands = opticalBands.addBands([ndvi, ndwi, ndbi]);
  }

  return opticalBands;
};

// 8. Harmonics
exports.getHarmonics = function(ccdcOutput, targetDateFractional, bands, coefficientScale, tolYears) {
  // Support swapped argument order defensively
  var _bands     = bands || ['nir', 'red', 'swir1'];
  var _coefScale = 10000;
  var _tolYears  = 1.0;

  if (typeof coefficientScale === 'number' && coefficientScale <= 10.0) {
    // Caller passed tolYears in position 4
    _tolYears  = coefficientScale;
    _coefScale = tolYears !== undefined ? tolYears : 10000;
  } else {
    _coefScale = coefficientScale !== undefined ? coefficientScale : 10000;
    _tolYears  = tolYears !== undefined ? tolYears : 1.0;
  }

  return core.getHarmonics(ccdcOutput, targetDateFractional, {
    bands:      _bands,
    coefScale:  _coefScale,
    tolYears:   _tolYears,
    coefCase:   'lower', // Sentinel-2 assets use lowercase
    outputCase: 'lower'  // S2 pipelines expect 'nir_amp', 'red_amp'
  });
};

// 9. Multi-Year Synthetic Stacks
exports.getMultiYearSyntheticStacks = function(ccdcOutput, yearsList, fractionOfYear, bands, addIndices, coefficientScale, tolYears) {
  fractionOfYear = fractionOfYear !== undefined ? fractionOfYear : 0.5;
  bands = bands || ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'];
  var result = {};
  yearsList.forEach(function(yr) {
    var dateFrac = yr + fractionOfYear;
    result[yr] = exports.getSyntheticStack(ccdcOutput, dateFrac, bands, addIndices, coefficientScale, tolYears);
  });
  return result;
};
