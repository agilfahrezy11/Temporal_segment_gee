/**
 * Landsat Continuous Change Detection and Classification Pipeline
 *
 * Coefficient units differ by source:
 *   - own run (this module's getLandsatTimeSeries feeds CCDC in 0..10,000)  -> coefScale = 10000
 *   - GOOGLE/GLOBAL_CCDC/V1 (coefficients are already 0..1 reflectance)     -> coefScale = 1
 */

// 1. Prepare Landsat 8 & 9 Collection for CCDC (0-10000 range)
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

// 2. Execute CCDC using Google Global CCDC hyperparameters
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

// 3. Nearest segment in time, per pixel (index + distance in years)
exports.nearestSegment = function(ccdc, t) {
  t = ee.Number(t);
  var tS = ccdc.select('tStart'), tE = ccdc.select('tEnd');
  var dist = tS.subtract(t).max(tE.multiply(-1).add(t)).max(0);   // 0 if t is inside a segment
  return {
    idx: dist.multiply(-1).arrayArgmax().arrayGet([0]),
    gap: dist.arrayReduce(ee.Reducer.min(), [0]).arrayGet([0]).rename('gap_years'),
    tS: tS, tE: tE
  };
};

// 4. Fill isolated holes from a 3x3 median, only where enough neighbours are valid
exports.fillIsolated = function(img, minNeighbors) {              // e.g. 5 of 8
  var n = img.select(0).mask().reduceNeighborhood({
    reducer: ee.Reducer.sum(), kernel: ee.Kernel.square(1)});
  var med = img.focalMedian(1, 'square', 'pixels');
  return img.unmask(med.updateMask(n.gte(minNeighbors)));
};

// 5. Synthetic reflectance stack at fractional-year t.
//    tolYears : max distance (years) to the nearest segment
//    coefScale: 10000 for the own run, 1 for GOOGLE/GLOBAL_CCDC/V1
exports.getSyntheticLandsatStack = function(ccdc, t, addIndices, tolYears, coefScale) {
  addIndices = addIndices !== undefined ? addIndices : true;
  tolYears = tolYears !== undefined ? tolYears : 1;
  coefScale = coefScale !== undefined ? coefScale : 10000;
  t = ee.Number(t);
  var seg = exports.nearestSegment(ccdc, t);

  // clamp only the trend term; harmonics keep the target date's season
  var tc = ee.Image(t).max(seg.tS.arrayGet(seg.idx)).min(seg.tE.arrayGet(seg.idx));
  var w = 2 * Math.PI;
  var terms = ee.Image.cat([
    ee.Image(1), tc,
    ee.Image.constant(t.multiply(w).cos()),   ee.Image.constant(t.multiply(w).sin()),
    ee.Image.constant(t.multiply(2 * w).cos()), ee.Image.constant(t.multiply(2 * w).sin()),
    ee.Image.constant(t.multiply(3 * w).cos()), ee.Image.constant(t.multiply(3 * w).sin())
  ]).toArray();

  var names = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
  var bands = names.map(function(b) {
    var c = ccdc.select(b + '_coefs')
      .arraySlice(0, seg.idx, seg.idx.add(1)).arrayProject([1]);
    return c.multiply(terms).arrayReduce(ee.Reducer.sum(), [0])
      .arrayGet([0]).divide(coefScale).rename(b.toLowerCase());
  });
  var stack = ee.Image.cat(bands);

  // plausibility bounds are lenient: dark water can be slightly negative in SWIR
  var valid = seg.gap.lte(tolYears)
    .and(stack.reduce(ee.Reducer.min()).gt(-0.05))
    .and(stack.reduce(ee.Reducer.max()).lt(1));
  stack = exports.fillIsolated(stack.updateMask(valid), 5);

  if (addIndices) {
    stack = stack.addBands([
      stack.normalizedDifference(['nir', 'red']).rename('NDVI'),
      stack.normalizedDifference(['green', 'nir']).rename('NDWI'),
      stack.normalizedDifference(['swir1', 'nir']).rename('NDBI')
    ]);
  }

  return stack.addBands(seg.gap);
};

// 6. First-harmonic amplitude and phase of the nearest segment
exports.getHarmonics = function(ccdc, t, bands, tolYears, coefScale) {
  tolYears = tolYears !== undefined ? tolYears : 1;
  coefScale = coefScale !== undefined ? coefScale : 10000;
  var seg = exports.nearestSegment(ccdc, t);
  var imgs = bands.map(function(b) {
    var c = ccdc.select(b + '_coefs')
      .arraySlice(0, seg.idx, seg.idx.add(1)).arrayProject([1]);   // 8 coefficients
    var a1 = c.arrayGet([2]).divide(coefScale);                    // cos(wt)
    var b1 = c.arrayGet([3]).divide(coefScale);                    // sin(wt)
    return ee.Image.cat([
      a1.hypot(b1).rename(b + '_amp'),
      b1.atan2(a1).rename(b + '_phase')
    ]);
  });
  return ee.Image.cat(imgs).updateMask(seg.gap.lte(tolYears));
};