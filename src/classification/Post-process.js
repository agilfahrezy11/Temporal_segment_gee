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

/**
 * =========================================================================
 * HYBRID PIXEL-OBJECT & SPATIAL STABILITY HELPERS
 * =========================================================================
 * Reference: Malinverni & Sanità (2026), Section 3.2.6
 */

/**
 * Computes Winning class (W), Runner-up class (S), Confusion Index (CI = S / W),
 * and binary Stability Map based on Malinverni & Sanità (2026).
 *
 * @param {ee.Image} probs - Image containing class probability bands
 * @param {Array|ee.List} [classBands] - Band names of class probabilities (defaults to probs.bandNames())
 * @param {Array|ee.List} [classValues] - Integer class IDs corresponding to classBands
 * @param {Object} [options] - Options:
 *    - stabilityTau: threshold below which pixel is considered stable (default: 0.65)
 * @returns {Object} { W, S, winner, runner_up, CI, is_stable }
 */
exports.computeStability = function(probs, classBands, classValues, options) {
  options = options || {};
  var tau = (options.stabilityTau !== undefined) ? options.stabilityTau : 0.65;

  var bandsArray = classBands
    ? (Array.isArray(classBands) ? classBands : classBands.getInfo())
    : probs.bandNames().getInfo();

  var valuesList = classValues
    ? ee.List(classValues)
    : ee.List.sequence(0, bandsArray.length - 1);

  var W = probs.reduce(ee.Reducer.max()).rename('W');
  var top1Idx = probs.toArray().arrayArgmax().arrayGet([0]);
  var top1Class = top1Idx.remap(ee.List.sequence(0, bandsArray.length - 1), valuesList).rename('winner');

  var valuesArray = classValues
    ? (Array.isArray(classValues) ? classValues : classValues.getInfo())
    : bandsArray.map(function(_, i) { return i; });

  var maskedBands = bandsArray.map(function(bName, idx) {
    var bImg = probs.select([bName]);
    var bId = valuesArray[idx];
    return bImg.where(top1Class.eq(bId), -1.0).rename(bName);
  });
  var maskedStack = ee.Image.cat(maskedBands);

  var S = maskedStack.reduce(ee.Reducer.max()).max(0.0).rename('S');
  var top2Idx = maskedStack.toArray().arrayArgmax().arrayGet([0]);
  var top2Class = top2Idx.remap(ee.List.sequence(0, bandsArray.length - 1), valuesList).rename('runner_up');

  var CI = S.divide(W.max(1e-4)).min(1.0).rename('confusion_index');
  var isStable = CI.lte(tau).rename('is_stable').toUint8();

  return {
    W: W,
    S: S,
    winner: top1Class,
    runner_up: top2Class,
    CI: CI,
    is_stable: isStable
  };
};

/**
 * Adaptive Hybrid Pixel-Object convex blending with uncertainty damping.
 * Blends pixel-level RF probabilities with SNIC superpixel object-mean probabilities.
 *
 * @param {ee.Image} probsImage - Pixel-level class probability stack (e.g. 0-1)
 * @param {ee.Image} clustersImage - Integer cluster label image ('clusters' or 'segment')
 * @param {Object} [options] - Configuration options:
 *    - classBands: Array or ee.List of probability band names
 *    - classValues: Array or ee.List of integer class IDs
 *    - weights: Object mapping class name -> object weight (0..1), or constant number (default: 0.5)
 *    - maxSize: max cluster size for reduceConnectedComponents (default: 1024)
 *    - adaptive: whether to damp object weights where objects are ambiguous/unstable (default: true)
 *    - adaptiveSensitivity: strength of damping factor (default: 0.6)
 *    - stabilityTau: threshold for stability calculation (default: 0.65)
 *    - useGuard: whether to damp object weight where within-object SD is high (default: false)
 *    - guardSd: SD threshold for guard check (default: 0.25)
 * @returns {Object} { probs, classification, conf, ci, runner_up, is_stable, combined }
 */
exports.adaptiveBlend = function(probsImage, clustersImage, options) {
  options = options || {};
  var maxSize = options.maxSize || 1024;
  var adaptive = (options.adaptive !== undefined) ? options.adaptive : true;
  var adaptiveSensitivity = (options.adaptiveSensitivity !== undefined) ? options.adaptiveSensitivity : 0.6;
  var useGuard = options.useGuard || false;
  var guardSd = (options.guardSd !== undefined) ? options.guardSd : 0.25;

  var classBands = options.classBands || probsImage.bandNames();
  var bandsArray = Array.isArray(classBands) ? classBands : classBands.getInfo();
  var classValues = options.classValues || bandsArray.map(function(_, i) { return i; });

  // 1. Cluster label band
  var labelBandName = clustersImage.bandNames().get(0);
  var labelBand = clustersImage.select([labelBandName]).toInt32().rename('clusters');

  // 2. Aggregate object-level mean probabilities
  var clusterImage = probsImage.select(bandsArray).addBands(labelBand);
  var obj = clusterImage.reduceConnectedComponents({
    reducer: ee.Reducer.mean(),
    labelBand: 'clusters',
    maxSize: maxSize
  }).select(bandsArray).unmask(probsImage.select(bandsArray));

  // 3. Object-level stability
  var objStability = exports.computeStability(obj, bandsArray, classValues, options);

  // 4. Construct weight image
  var weightsOpt = options.weights || 0.5;
  var wBands = bandsArray.map(function(b) {
    var val = (typeof weightsOpt === 'object' && weightsOpt !== null && weightsOpt[b] !== undefined)
      ? weightsOpt[b]
      : (typeof weightsOpt === 'number' ? weightsOpt : 0.5);
    return ee.Image.constant(ee.Number(val)).rename([b]);
  });
  var w = ee.Image.cat(wBands);

  var wEff = w;
  if (adaptive) {
    var stabilityDamping = ee.Image.constant(1.0).subtract(objStability.CI.multiply(adaptiveSensitivity));
    wEff = w.multiply(stabilityDamping);
  }

  if (useGuard) {
    var sd = clusterImage.reduceConnectedComponents({
      reducer: ee.Reducer.stdDev(),
      labelBand: 'clusters',
      maxSize: maxSize
    }).select(bandsArray).reduce(ee.Reducer.max());
    wEff = wEff.multiply(sd.unmask(0).lte(guardSd));
  }

  // 5. Convex blending: blend = obj * wEff + pixel * (1 - wEff)
  var blend = obj.multiply(wEff).add(probsImage.select(bandsArray).multiply(ee.Image.constant(1.0).subtract(wEff)));
  var blendTotal = blend.reduce(ee.Reducer.sum());
  blend = blend.divide(blendTotal.where(blendTotal.eq(0), 1.0)).rename(bandsArray);

  // 6. Final thematic labeling and quality assessment layers
  var finalStability = exports.computeStability(blend, bandsArray, classValues, options);
  var lulc = finalStability.winner.rename('lulc').toUint8();
  var conf = finalStability.W.multiply(100).round().toUint8().rename('conf');
  var ciPct = finalStability.CI.multiply(100).round().toUint8().rename('ci');
  var runnerUp = finalStability.runner_up.toUint8().rename('runner_up');
  var isStable = finalStability.is_stable.rename('is_stable');

  var combined = lulc.addBands([conf, ciPct, runnerUp, isStable]);

  return {
    probs: blend,
    classification: lulc,
    conf: conf,
    ci: ciPct,
    runner_up: runnerUp,
    is_stable: isStable,
    combined: combined
  };
};

/**
 * Intrinsic Rule-Based Hydrology Refinement (Malinverni & Sanità 2026, Section 3.2.6)
 * Suppresses false hydrology classifications caused by spectral overlap between
 * urban/cleared surfaces and waterbody / dry riverbed features.
 *
 * @param {ee.Image} lulcImage - Thematic class image containing 'lulc' or discrete class IDs
 * @param {ee.Image} runnerUpImage - Runner-up class image
 * @param {ee.Image} ciImage - Confusion Index image (either 0-100 or 0-1)
 * @param {Object} [options] - Configuration:
 *    - ciThreshold: CI threshold above which runner-up replaces winner (default: 65 on 0-100 scale, or 0.65 on 0-1)
 *    - waterClasses: array of water/riverbed IDs (default: [17, 18])
 *    - terrestrialClasses: array of urban/cleared IDs (default: [15, 16])
 *    - riverbedId: ID for dry riverbed (default: 18)
 *    - waterbodyId: ID for waterbody (default: 17)
 *    - clearedLandId: ID for cleared land (default: 15)
 *    - mergeRiverbed: whether to merge riverbed into waterbody (default: true)
 * @returns {ee.Image} Refined thematic classification image
 */
exports.suppressFalseHydrology = function(lulcImage, runnerUpImage, ciImage, options) {
  options = options || {};
  var lulcBand = lulcImage.bandNames().contains('lulc')
    ? lulcImage.select('lulc')
    : lulcImage.select(0);

  var runnerUp = runnerUpImage.select(0);
  var ci = ciImage.select(0);

  var ciThresh = (options.ciThreshold !== undefined) ? options.ciThreshold : 65;
  var waterClasses = options.waterClasses || [17, 18];
  var terrestrialClasses = options.terrestrialClasses || [15, 16];
  var riverbedId = (options.riverbedId !== undefined) ? options.riverbedId : 18;
  var waterbodyId = (options.waterbodyId !== undefined) ? options.waterbodyId : 17;
  var clearedLandId = (options.clearedLandId !== undefined) ? options.clearedLandId : 15;
  var mergeRiverbed = (options.mergeRiverbed !== undefined) ? options.mergeRiverbed : true;

  // Build condition: winner is in waterClasses
  var isWater = ee.Image.constant(0);
  waterClasses.forEach(function(cId) {
    isWater = isWater.or(lulcBand.eq(cId));
  });

  // Runner-up is in terrestrialClasses
  var isTerrestrial = ee.Image.constant(0);
  terrestrialClasses.forEach(function(cId) {
    isTerrestrial = isTerrestrial.or(runnerUp.eq(cId));
  });

  // Rule 1: False Hydrology Suppression
  var falseHydrology = isWater.and(isTerrestrial).and(ci.gt(ciThresh));
  var refined = lulcBand.where(falseHydrology, runnerUp);

  // Rule 2: Cleared Land vs. Dry Riverbed Check
  var riverbedThresh = (options.riverbedCiThreshold !== undefined) ? options.riverbedCiThreshold : 75;
  var falseRiverbed = refined.eq(riverbedId).and(runnerUp.eq(clearedLandId)).and(ci.gt(riverbedThresh));
  refined = refined.where(falseRiverbed, runnerUp);

  // Rule 3: Merge Dry Riverbed into Waterbody post-classification
  if (mergeRiverbed) {
    refined = refined.where(refined.eq(riverbedId), waterbodyId);
  }

  return refined.rename('lulc_refined');
};