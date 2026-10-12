/**
 * Module: Post-Processing & Spatial Smoothing
 * Description: Empirical Bayes ("Bayesian") spatial smoothing for land cover classification 
 *              probability images. Based on Camara et al. (2024).
 * Reference: https://e-sensing.github.io/sitsbook/cl_smoothing.html
 */

/**
 * Module: Post-Processing & Spatial Smoothing (Fixed & Robust)
 */
var EPS = 1e-6;

function probsToLogit(probsImage, classBands) {
  var p = probsImage.select(classBands).clamp(EPS, 1 - EPS);
  var one = ee.Image.constant(1);
  return p.divide(one.subtract(p)).log().rename(classBands);
}

function logitToProbs(logitImage, classBands, normalize) {
  normalize = (normalize === undefined) ? true : normalize;
  var p = ee.Image.constant(1)
    .divide(ee.Image.constant(1).add(logitImage.multiply(-1).exp()))
    .rename(classBands);

  if (normalize) {
    var total = p.reduce(ee.Reducer.sum());
    p = p.divide(total).rename(classBands);
  }
  return p;
}

function buildSigma2Image(classBands, smoothness) {
  var bandsArray = Array.isArray(classBands) ? classBands : classBands.getInfo();
  
  var bands = bandsArray.map(function(b) {
    var val = (typeof smoothness === 'object' && smoothness !== null && smoothness[b] !== undefined)
      ? smoothness[b]
      : (typeof smoothness === 'number' ? smoothness : 1.0);
    return ee.Image.constant(ee.Number(val)).rename([b]);
  });
  return ee.Image.cat(bands);
}

function localMeanVariance(logitImage, classBands, windowSize, neighFraction) {
  windowSize = windowSize || 7;
  neighFraction = neighFraction || 0.5;

  var radius = Math.floor((windowSize - 1) / 2);
  var kernel = ee.Kernel.square(radius, 'pixels');

  var keepPct = Math.round((1 - neighFraction) * 100);
  var localThreshold = logitImage
    .reduceNeighborhood({reducer: ee.Reducer.percentile([keepPct]), kernel: kernel})
    .rename(classBands);

  var selected = logitImage.updateMask(logitImage.gte(localThreshold));

  var localMean = selected
    .reduceNeighborhood({reducer: ee.Reducer.mean(), kernel: kernel})
    .rename(classBands);

  var localVar = selected
    .reduceNeighborhood({reducer: ee.Reducer.variance(), kernel: kernel})
    .rename(classBands);

  localVar = localVar.unmask(1.0).where(localVar.eq(0), 1.0);
  localMean = localMean.unmask(logitImage);

  return { mean: localMean, variance: localVar };
}

exports.localVariance = function(probsImage, classBands, windowSize, neighFraction) {
  var logitImage = probsToLogit(probsImage, classBands);
  var stats = localMeanVariance(logitImage, classBands, windowSize, neighFraction);
  return stats.variance;
};

exports.bayesianSmooth = function(probsImage, classBands, smoothness, windowSize, neighFraction) {
  windowSize = windowSize || 7;
  neighFraction = neighFraction || 0.5;

  var logitImage = probsToLogit(probsImage, classBands);
  var stats = localMeanVariance(logitImage, classBands, windowSize, neighFraction);
  var sigma2 = buildSigma2Image(classBands, smoothness);

  var wX = stats.variance.divide(sigma2.add(stats.variance));
  var wM = sigma2.divide(sigma2.add(stats.variance));

  var posteriorLogit = logitImage.multiply(wX)
    .add(stats.mean.multiply(wM))
    .rename(classBands);

  return logitToProbs(posteriorLogit, classBands, true);
};

exports.classifyFromProbs = function(probsImage, classBands, classValues) {
  var arr = probsImage.select(classBands).toArray();
  var maxIdx = arr.arrayArgmax().arrayGet([0]).rename('classification');

  if (classValues) {
    var bandsArray = Array.isArray(classBands) ? classBands : classBands.getInfo();
    var seq = ee.List.sequence(0, bandsArray.length - 1);
    
    // Wrapped classValues in ee.List() to prevent type mismatch
    maxIdx = maxIdx.remap(seq, ee.List(classValues)).rename('classification');
  }

  return maxIdx.toByte();
};

/**
 * =========================================================================
 * GEOBIA: Geographic Object-Based Image Analysis (SNIC Superpixel Functions)
 * =========================================================================
 */

/**
 * Delineates spatial superpixel objects using SNIC (Simple Non-Iterative Clustering).
 * 
 * @param {ee.Image} featureImage - Multi-band image driving cluster boundaries (e.g., optical bands, NDVI)
 * @param {Object} [options] - SNIC parameters:
 *    - size: seed spacing in pixels (default: 10, ~1 ha at 10m resolution)
 *    - compactness: compactness factor (default: 0.5; lower = adheres more to spectral boundaries)
 *    - connectivity: pixel connectivity (4 or 8, default: 8)
 *    - neighborhoodSize: search neighborhood (default: 2 * size)
 *    - seeds: optional pre-defined seed grid image
 * @returns {ee.Image} SNIC result image containing the 'clusters' band and band cluster means
 */
exports.snicSuperpixels = function(featureImage, options) {
  options = options || {};
  var size = options.size || 10;
  var compactness = (options.compactness !== undefined) ? options.compactness : 0.5;
  var connectivity = options.connectivity || 8;
  var neighborhoodSize = options.neighborhoodSize || (size * 2);
  var seeds = options.seeds || null;

  var snic = ee.Algorithms.Image.Segmentation.SNIC({
    image: featureImage,
    size: size,
    compactness: compactness,
    connectivity: connectivity,
    neighborhoodSize: neighborhoodSize,
    seeds: seeds
  });

  return snic;
};

/**
 * Aggregates continuous class probabilities across SNIC superpixel objects by computing
 * the mean probability per object and re-normalizing.
 * 
 * @param {ee.Image} probsImage - Image containing class probability bands
 * @param {ee.Image} clusters - Single band integer image representing cluster IDs ('clusters')
 * @param {number} [maxSize=1024] - Maximum connected component size in pixels for reduceConnectedComponents
 * @returns {ee.Image} Object-mean probabilities stack
 */
exports.aggregateProbsByClusters = function(probsImage, clusters, maxSize) {
  maxSize = maxSize || 1024;
  var labelBand = clusters.select(['clusters']).toInt32();

  var clusterImage = probsImage.addBands(labelBand);
  var objectProbs = clusterImage.reduceConnectedComponents({
    reducer: ee.Reducer.mean(),
    labelBand: 'clusters',
    maxSize: maxSize
  });

  // Fallback to original pixel probabilities if cluster exceeds maxSize
  objectProbs = objectProbs.unmask(probsImage);

  // Normalize so probabilities sum to 1.0 across bands
  var sumProbs = objectProbs.reduce(ee.Reducer.sum());
  objectProbs = objectProbs.divide(sumProbs.where(sumProbs.eq(0), 1.0));

  return objectProbs;
};

/**
 * High-level GEOBIA post-processing pipeline:
 * 1. Clusters landscape into superpixels using SNIC on optical/spectral features.
 * 2. Averages class probabilities within each spatial object.
 * 3. Classifies objects via argmax into discrete LULC classes.
 * 
 * @param {ee.Image} probsImage - Class probabilities (optionally Bayesian-smoothed)
 * @param {ee.Image} featureImage - Optical bands / indices driving spatial segmentation
 * @param {Array|ee.List} classBands - Names of the class probability bands
 * @param {Array|ee.List} [classValues] - Target integer class IDs for remapping
 * @param {Object} [options] - Configuration for SNIC and connected components:
 *    - size: seed spacing (default: 10)
 *    - compactness: shape compactness factor (default: 0.5)
 *    - maxSize: max cluster size (default: 1024)
 * @returns {Object} { probs: ee.Image, classification: ee.Image, clusters: ee.Image }
 */
exports.applySnicOBIA = function(probsImage, featureImage, classBands, classValues, options) {
  options = options || {};
  var snic = exports.snicSuperpixels(featureImage, options);
  var clusters = snic.select('clusters');
  var objectProbs = exports.aggregateProbsByClusters(probsImage, clusters, options.maxSize);
  var objectClass = exports.classifyFromProbs(objectProbs, classBands, classValues);

  return {
    probs: objectProbs,
    classification: objectClass,
    clusters: clusters
  };
};