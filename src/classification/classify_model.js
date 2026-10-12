/**
 * Module: Classification & Model Evaluation
 * Description: Implements Hard Classification, One-vs-Rest (OVR) Soft Probability 
 *              Classification, Multi-Probability Classification , 
 *              Feature Importance extraction, and Model Accuracy Evaluation.
 */

/**
 * Helper to build Random Forest options dictionary for GEE constructor.
 * Omits variablesPerSplit if not provided so GEE defaults to sqrt(#bands).
 */
function buildRfParams(nTrees, vSplit, minLeaf, seed) {
  var params = {
    numberOfTrees: nTrees !== undefined ? nTrees : 100,
    minLeafPopulation: minLeaf !== undefined ? minLeaf : 1,
    seed: seed !== undefined ? seed : 0
  };
  if (vSplit !== undefined && vSplit !== null && vSplit > 0) {
    params.variablesPerSplit = vSplit;
  }
  return params;
}

/**
 * 1. Hard Classification
 * Performs standard multi-class Random Forest classification.
 */
exports.hardClassification = function(trainingPixels, classProperty, image, nTrees, vSplit, minLeaf, seed) {
  var bandNames = image.bandNames();
  var rfParams = buildRfParams(nTrees, vSplit, minLeaf, seed);

  var classifier = ee.Classifier.smileRandomForest(rfParams);
  var trainedModel = classifier.train({
    features: trainingPixels,
    classProperty: classProperty,
    inputProperties: bandNames
  });

  var classificationMap = image.classify(trainedModel).rename('classification');

  return {
    classificationMap: classificationMap,
    trainedModel: trainedModel
  };
};

/**
 * 2. Soft Classification (One-vs-Rest Binary Probability Stack & Argmax)
 * Trains individual binary Random Forest models in PROBABILITY mode for each class.
 */
exports.softClassification = function(
  trainingPixels,
  classProperty,
  image,
  includeFinalMap,
  nTrees,
  vSplit,
  minLeaf,
  probabilityScale,
  seed
) {
  includeFinalMap = includeFinalMap !== undefined ? includeFinalMap : true;
  probabilityScale = probabilityScale !== undefined ? probabilityScale : 100;

  var bandNames = image.bandNames();
  var rfParams = buildRfParams(nTrees, vSplit, minLeaf, seed);
  var classList = trainingPixels.aggregate_array(classProperty).distinct();

  var probImages = classList.map(function(classId) {
    classId = ee.Number(classId);

    var binaryTrain = trainingPixels.map(function(ft) {
      var isTarget = ee.Number(ft.get(classProperty)).eq(classId);
      return ft.set('binary', ee.Algorithms.If(isTarget, 1, 0));
    });

    var classifier = ee.Classifier.smileRandomForest(rfParams).setOutputMode('PROBABILITY');

    var trained = classifier.train({
      features: binaryTrain,
      classProperty: 'binary',
      inputProperties: bandNames
    });

    var probImg = image.classify(trained)
      .multiply(probabilityScale)
      .round()
      .toByte();

    var bandName = ee.String('prob_').cat(classId.int().format());
    return probImg.rename(bandName);
  });

  var probStack = ee.ImageCollection(probImages).toBands();
  
  // Clean up band names (e.g., "0_prob_1" -> "prob_1")
  var originalBandNames = probStack.bandNames();
  var cleanBandNames = originalBandNames.map(function(name) {
    var str = ee.String(name);
    var split = str.split('_');
    return ee.String('prob_').cat(split.get(2));
  });
  probStack = probStack.select(originalBandNames, cleanBandNames);

  if (!includeFinalMap) {
    return probStack;
  }

  var arrayImg = probStack.toArray();
  var maxProbIndex = arrayImg.arrayArgmax().arrayGet([0]);

  var classSequence = ee.List.sequence(0, classList.size().subtract(1));
  var finalMap = maxProbIndex.remap(classSequence, classList).rename('classification');

  var maxConfidence = arrayImg.arrayReduce(ee.Reducer.max(), [0])
    .arrayGet([0])
    .rename('confidence');

  return probStack.addBands([finalMap, maxConfidence]);
};

/**
 * 3. Multi-Probability Classification (Native RF Array Probabilities)
 * Returns a multi-band probability image [0,1] formatted for Bayesian smoothing.
 */
exports.multiProbabilityClassification = function(trainingPixels, classProperty, classIdToNameMap, image, options) {
  options = options || {};
  var rfParams = buildRfParams(options.nTrees, options.vSplit, options.minLeaf, options.seed);
  var bandNames = image.bandNames();

  var classifier = ee.Classifier.smileRandomForest(rfParams).setOutputMode('MULTIPROBABILITY');

  var trainedModel = classifier.train({
    features: trainingPixels,
    classProperty: classProperty,
    inputProperties: bandNames
  });

  var probsArrayImage = image.classify(trainedModel);

  // Parse class IDs in ascending numerical order (matching GEE array element order)
  var sortedIds = Object.keys(classIdToNameMap).map(Number).sort(function(a, b) { return a - b; });
  var classBands = sortedIds.map(function(id) { return classIdToNameMap[id]; });

  var probsImage = probsArrayImage.rename('probs').arrayFlatten([classBands]);

  return {
    probsImage: probsImage,
    classBands: classBands,
    classValues: sortedIds,
    trainedModel: trainedModel
  };
};

/**
 * 4. Unified Classification Selector
 * Supports 'hard', 'soft', and 'multi' / 'multiprobability' modes.
 */
exports.classify = function(mode, trainingPixels, classProperty, image, options) {
  options = options || {};
  mode = (mode || 'hard').toLowerCase();

  if (mode === 'hard') {
    return exports.hardClassification(
      trainingPixels,
      classProperty,
      image,
      options.nTrees,
      options.vSplit,
      options.minLeaf,
      options.seed
    );
  } else if (mode === 'soft') {
    return exports.softClassification(
      trainingPixels,
      classProperty,
      image,
      options.includeFinalMap,
      options.nTrees,
      options.vSplit,
      options.minLeaf,
      options.probabilityScale,
      options.seed
    );
  } else if (mode === 'multi' || mode === 'multiprobability') {
    if (!options.classIdToNameMap) {
      throw new Error("Mode 'multi' requires options.classIdToNameMap (e.g., {1: 'Water', 2: 'Forest'}).");
    }
    return exports.multiProbabilityClassification(
      trainingPixels,
      classProperty,
      options.classIdToNameMap,
      image,
      options
    );
  } else {
    throw new Error("Invalid mode: '" + mode + "'. Choose 'hard', 'soft', or 'multi'.");
  }
};

/**
 * 5. Feature Importance
 */
exports.getFeatureImportance = function(trainedModel) {
  var explanation = trainedModel.explain();
  return ee.Dictionary(explanation.get('importance'));
};

/**
 * 6. Model Evaluation Metrics
 *
 * @param {ee.Classifier} trainedModel - Trained hard classification model.
 * @param {ee.FeatureCollection} testPixels - Test/validation pixels.
 * @param {string} classProperty - Property name holding the ground truth class integer.
 * @param {Object} [classIdToNameMap] - Optional map of integer id -> class name, e.g. {1: 'Forest', 2: 'Water'}.
 *                                     If provided, adds labeled dictionaries for PA, CA, and F1-score.
 * @returns {Object} Evaluation metrics including error matrix, overall accuracy, kappa,
 *                   and per-class producer's accuracy, consumer's accuracy, and F1 score.
 */
exports.evaluateModel = function(trainedModel, testPixels, classProperty, classIdToNameMap) {
  var testClassified = testPixels.classify(trainedModel);
  var errorMatrix = testClassified.errorMatrix(classProperty, 'classification');

  var result = {
    errorMatrix: errorMatrix,
    overallAccuracy: errorMatrix.accuracy(),
    kappa: errorMatrix.kappa(),
    producersAccuracy: errorMatrix.producersAccuracy(),
    consumersAccuracy: errorMatrix.consumersAccuracy()
  };

  // If classIdToNameMap is provided, construct labeled dictionaries for each class
  if (classIdToNameMap) {
    var sortedIds = Object.keys(classIdToNameMap).map(Number).sort(function(a, b) { return a - b; });
    var classNames = sortedIds.map(function(id) { return classIdToNameMap[id]; });
    var eeClassNames = ee.List(classNames);

    // Flatten 2D matrices into 1D Lists
    var paList = ee.List(errorMatrix.producersAccuracy().flatten().toList());
    var caList = ee.List(errorMatrix.consumersAccuracy().flatten().toList());

    result.producersAccuracyByName = ee.Dictionary.fromLists(eeClassNames, paList);
    result.consumersAccuracyByName = ee.Dictionary.fromLists(eeClassNames, caList);

    // Compute per-class F1-score: 2 * (PA * CA) / (PA + CA)
    var f1List = ee.List.sequence(0, eeClassNames.size().subtract(1)).map(function(i) {
      i = ee.Number(i);
      var pa = ee.Number(paList.get(i));
      var ca = ee.Number(caList.get(i));
      var sum = pa.add(ca);
      return ee.Algorithms.If(sum.gt(0), pa.multiply(ca).multiply(2).divide(sum), ee.Number(0));
    });
    result.f1ScoreByName = ee.Dictionary.fromLists(eeClassNames, f1List);
  }

  return result;
};

/**
 * 7. Probability Evaluation — Log Loss & Top-Class Accuracy (overall and per-class)
 *
 * Computes log loss and top-class accuracy using predicted probabilities from a
 * MULTIPROBABILITY classifier against ground-truth labels in testPixels.
 *
 * Log loss per sample = -log(p_true_class + eps)
 * Overall log loss    = mean over all test samples
 * Per-class log loss  = mean over samples belonging to that class
 *
 * @param {ee.FeatureCollection|ee.Classifier} arg1 - testPixels or trainedModel (flexible argument order).
 * @param {ee.Classifier|ee.FeatureCollection} arg2 - trainedModel or testPixels.
 * @param {string}               classProperty - Property name holding the integer class label.
 * @param {Object}               classIdToNameMap - Mapping of integer class id -> band name string,
 *                               e.g. {1: 'Water', 2: 'Forest'}.
 * @param {number}               [eps=1e-7] - Small epsilon added to probabilities before log to
 *                               avoid log(0). Defaults to 1e-7.
 * @returns {ee.Dictionary} Dictionary with keys:
 *   - 'overallLogLoss'      {ee.Number}     — mean log loss across all samples
 *   - 'overallAccuracy'     {ee.Number}     — overall argmax top-class accuracy
 *   - 'perClassLogLoss'     {ee.Dictionary} — per-class name -> mean log loss for that class
 *   - 'perClassAccuracy'    {ee.Dictionary} — per-class name -> top-class accuracy (recall)
 *   - 'nSamples'            {ee.Number}     — total number of test samples evaluated
 *   - 'perClassNSamples'    {ee.Dictionary} — per-class name -> sample count
 */
exports.evaluateProbabilities = function(arg1, arg2, classProperty, classIdToNameMap, eps) {
  eps = (eps !== undefined && eps !== null) ? eps : 1e-7;

  // Flexible argument handling: supports (testPixels, trainedModel, ...) or (trainedModel, testPixels, ...)
  var testPixels, trainedModel;
  if (arg1 instanceof ee.Classifier) {
    trainedModel = arg1;
    testPixels   = arg2;
  } else {
    testPixels   = arg1;
    trainedModel = arg2;
  }

  // Parse sorted class ids and matching band names (must mirror multiProbabilityClassification)
  var sortedIds   = Object.keys(classIdToNameMap).map(Number).sort(function(a, b) { return a - b; });
  var classBands  = sortedIds.map(function(id) { return classIdToNameMap[id]; });
  var eeClassIds  = ee.List(sortedIds);
  var eeBandNames = ee.List(classBands);

  // Classify test pixels into a probability list property named 'probs'
  var scored = testPixels.classify(trainedModel, 'probs');

  // Compute per-sample metrics using server-side ee.List and ee.Array operations
  var withMetrics = scored.map(function(ft) {
    var probsArr    = ee.Array(ft.get('probs'));
    var probs       = probsArr.toList();
    var trueClassId = ee.Number(ft.get(classProperty));
    var trueIdx     = eeClassIds.indexOf(trueClassId);

    // Probability of the true ground-truth class
    var pTrue       = ee.Number(probs.get(trueIdx)).max(eps);
    var sampleLoss  = pTrue.log().multiply(-1);

    // Argmax predicted index (top-class); argmax() returns a list of coordinates
    var predIdx     = ee.Number(probsArr.argmax().get(0));
    var isCorrect   = predIdx.eq(trueIdx);

    return ft.set({
      '_log_loss': sampleLoss,
      '_correct':  ee.Number(isCorrect)
    });
  });

  // Overall metrics across all samples
  var overallLogLoss  = withMetrics.aggregate_mean('_log_loss');
  var overallAccuracy = withMetrics.aggregate_mean('_correct');
  var nSamples        = withMetrics.size();

  // Per-class metrics within each ground truth class subset
  var perClassLogLoss = ee.Dictionary.fromLists(
    eeBandNames,
    eeClassIds.map(function(classId) {
      classId = ee.Number(classId);
      var classSubset = withMetrics.filter(ee.Filter.eq(classProperty, classId));
      return ee.Algorithms.If(
        classSubset.size().gt(0),
        classSubset.aggregate_mean('_log_loss'),
        null
      );
    })
  );

  var perClassAccuracy = ee.Dictionary.fromLists(
    eeBandNames,
    eeClassIds.map(function(classId) {
      classId = ee.Number(classId);
      var classSubset = withMetrics.filter(ee.Filter.eq(classProperty, classId));
      return ee.Algorithms.If(
        classSubset.size().gt(0),
        classSubset.aggregate_mean('_correct'),
        null
      );
    })
  );

  var perClassNSamples = ee.Dictionary.fromLists(
    eeBandNames,
    eeClassIds.map(function(classId) {
      classId = ee.Number(classId);
      return withMetrics.filter(ee.Filter.eq(classProperty, classId)).size();
    })
  );

  return ee.Dictionary({
    overallLogLoss:    overallLogLoss,
    overallAccuracy:   overallAccuracy,
    perClassLogLoss:   perClassLogLoss,
    perClassAccuracy:  perClassAccuracy,
    nSamples:          nSamples,
    perClassNSamples:  perClassNSamples
  });
};