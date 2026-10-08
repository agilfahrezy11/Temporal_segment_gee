/**
 * Landsat Continous Change Detection and Classification Pipeline
 */

// 1. Prepare Harmonized Landsat 8 & 9 Collection for CCDC (0-10000 range)
exports.getLandsatTimeSeries = function(roi, startDate, endDate) {
  startDate = startDate || '2020-01-01';
  endDate = endDate || '2025-12-31';

  var maskAndScale = function(img) {
    var qa = img.select('QA_PIXEL');
    // Mask clouds (bit 3) and cloud shadows (bit 4)
    var cloudShadowBitMask = 1 << 4;
    var cloudsBitMask = 1 << 3;
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));

    // Convert USGS Collection 2 SR (0.0000275 * DN - 0.2) to 0..10,000 scale
    var srBands = img.select(
      ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
      ['BLUE',  'GREEN', 'RED',   'NIR',   'SWIR1', 'SWIR2']
    ).multiply(0.0000275).subtract(0.2).multiply(10000);

    return srBands.updateMask(mask)
      .copyProperties(img, ['system:time_start']);
  };

  var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .map(maskAndScale);

  var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .map(maskAndScale);

  return l8.merge(l9);
};

// 2. Execute CCDC using Google Global CCDC Hyperparameters
exports.runLandsatCCDC = function(lsCollection) {
  return ee.Algorithms.TemporalSegmentation.Ccdc({
    collection: lsCollection,
    breakpointBands: ['GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'],
    tmaskBands: ['GREEN', 'SWIR1'],
    dateFormat: 1, // Fractional years
    minObservations: 6,
    chiSquareProbability: 0.99,
    minNumOfYearsScaler: 1.33,
    lambda: 20,
    maxIterations: 25000
  });
};

// 3. Synthesize Reflectance Stack from CCDC Output
exports.getSyntheticLandsatStack = function(ccdcOutput, targetDateFractional, addIndices) {
  addIndices = addIndices !== undefined ? addIndices : true;
  var bands = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
  var t = ee.Number(targetDateFractional);

  var tStart = ccdcOutput.select('tStart');
  var tEnd = ccdcOutput.select('tEnd');

  // Segment selection mask
  var segmentMask = tStart.lte(t).and(tEnd.gte(t));
  var segmentMask2D = segmentMask.toArray(0).toArray(1);

  // Harmonic terms matrix [8 x 1]
  var omega = 2.0 * Math.PI;
  var termList = [
    ee.Image.constant(1),
    ee.Image.constant(t),
    ee.Image.constant(t.multiply(omega).cos()),
    ee.Image.constant(t.multiply(omega).sin()),
    ee.Image.constant(t.multiply(omega * 2).cos()),
    ee.Image.constant(t.multiply(omega * 2).sin()),
    ee.Image.constant(t.multiply(omega * 3).cos()),
    ee.Image.constant(t.multiply(omega * 3).sin())
  ];
  var termsArray = ee.Image(termList).toArray().toArray(1);

  var syntheticBands = bands.map(function(band) {
    var coefs = ccdcOutput.select(band + '_coefs'); 
    var activeCoefs = coefs.arrayMask(segmentMask2D);
    var segCoefs = activeCoefs.arraySlice(0, 0, 1);
    var fitArray = segCoefs.matrixMultiply(termsArray);
    
    // Convert 0..10,000 integer range back to 0..1 reflectance
    return fitArray.arrayProject([0]).arrayFlatten([[band.toLowerCase()]]).divide(10000);
  });

  var stack = ee.Image(syntheticBands);

  if (addIndices) {
    var ndvi = stack.normalizedDifference(['nir', 'red']).rename('NDVI');
    var ndwi = stack.normalizedDifference(['green', 'nir']).rename('NDWI');
    var ndbi = stack.normalizedDifference(['swir1', 'nir']).rename('NDBI');
    stack = stack.addBands([ndvi, ndwi, ndbi]);
  }

  return stack;
};