/**
 * Module: Unified CCDC Core
 * =========================================================================
 * Provides sensor-agnostic CCDC execution, nearest-segment lookup, trend-clamped
 * synthetic reflectance reconstruction, harmonic metric extraction, and
 * multi-epoch stability filtering.
 *
 * Applicable to both Landsat and Sentinel-2 CCDC array outputs.
 * =========================================================================
 */

/**
 * Executes the CCDC algorithm on a preprocessed optical collection.
 *
 * @param {ee.ImageCollection} collection - Harmonized collection (scaled 0..10,000)
 * @param {Object}             [options]  - CCDC hyperparameter overrides:
 *   @param {Array<string>}      [options.breakpointBands=['green', 'red', 'nir', 'swir1', 'swir2']]
 *   @param {Array<string>}      [options.tmaskBands] - Optional multi-temporal mask bands (e.g. ['green', 'swir1'])
 *   @param {number}             [options.minObservations=6]
 *   @param {number}             [options.chiSquareProbability=0.99]
 *   @param {number}             [options.minNumOfYearsScaler=1.33]
 *   @param {number}             [options.lambda=20]
 *   @param {number}             [options.maxIterations=25000]
 * @returns {ee.Image} CCDC array output image containing coefficients and break metadata
 */
exports.runCCDC = function(collection, options) {
  options = options || {};

  var params = {
    collection:           collection,
    breakpointBands:      options.breakpointBands || ['green', 'red', 'nir', 'swir1', 'swir2'],
    dateFormat:           1, // Fractional years (e.g., 2025.65)
    minObservations:      options.minObservations || 6,
    chiSquareProbability: options.chiSquareProbability || 0.99,
    minNumOfYearsScaler:  options.minNumOfYearsScaler || 1.33,
    lambda:               options.lambda || 20,
    maxIterations:        options.maxIterations || 25000
  };

  if (options.tmaskBands) {
    params.tmaskBands = options.tmaskBands;
  }

  return ee.Algorithms.TemporalSegmentation.Ccdc(params);
};

/**
 * Finds the nearest CCDC segment in time per pixel and calculates the temporal distance.
 * Returns gap = 0 if target date t is strictly inside an active segment.
 *
 * @param {ee.Image}         ccdc - CCDC array image containing 'tStart' and 'tEnd'
 * @param {number|ee.Number} t    - Target fractional year (e.g. 2025.65)
 * @returns {Object} { idx: ee.Image, gap: ee.Image, tS: ee.Image, tE: ee.Image }
 */
exports.nearestSegment = function(ccdc, t) {
  t = ee.Number(t);
  var tS = ccdc.select('tStart');
  var tE = ccdc.select('tEnd');

  // Distance is 0 if inside [tStart, tEnd], else distance to nearest segment edge
  var dist = tS.subtract(t).max(tE.multiply(-1).add(t)).max(0);

  return {
    idx: dist.multiply(-1).arrayArgmax().arrayGet([0]),
    gap: dist.arrayReduce(ee.Reducer.min(), [0]).arrayGet([0]).rename('gap_years'),
    tS:  tS,
    tE:  tE
  };
};

/**
 * Fills isolated nodata holes using 3x3 focal median only where enough neighbors are valid.
 *
 * @param {ee.Image} img                   - Image with potential isolated holes
 * @param {number}   [minNeighbors=5]      - Minimum valid neighbors required in 3x3 window
 * @returns {ee.Image} Gap-filled image
 */
exports.fillIsolated = function(img, minNeighbors) {
  minNeighbors = minNeighbors !== undefined ? minNeighbors : 5;
  var n = img.select(0).mask().reduceNeighborhood({
    reducer: ee.Reducer.sum(),
    kernel:  ee.Kernel.square(1)
  });
  var med = img.focalMedian(1, 'square', 'pixels');
  return img.unmask(med.updateMask(n.gte(minNeighbors)));
};

/**
 * Generates a synthetic surface reflectance stack evaluated at fractional year t.
 * Implements:
 *   - Nearest-segment lookup with configurable gap tolerance (tolYears)
 *   - Trend-clamping (tc) to segment bounds to prevent runaway linear extrapolation
 *   - Physical plausibility masking with unified bounds ([-0.05, 1.0])
 *   - 3x3 focal median hole repair for isolated pixels
 *   - Diagnostic bands: 'gap_years' and 'synth_valid'
 *
 * @param {ee.Image}         ccdcOutput            - CCDC array output image
 * @param {number|ee.Number} targetDateFractional  - Date in fractional years (e.g. 2025.65)
 * @param {Object}           [options]             - Configuration options:
 *   @param {Array<string>}    [options.bands=['blue', 'green', 'red', 'nir', 'swir1', 'swir2']]
 *   @param {number}           [options.coefScale=10000] - Scaling divisor for coefficients
 *   @param {number}           [options.tolYears=1.0]    - Max allowable gap to segment in years
 *   @param {string}           [options.coefCase='lower']- 'lower' (blue_coefs) or 'upper' (BLUE_coefs)
 *   @param {number}           [options.minReflectance=-0.05]
 *   @param {number}           [options.maxReflectance=1.0]
 *   @param {boolean}          [options.fillHoles=true]
 * @returns {ee.Image} Multi-band synthetic reflectance stack (0..1 scale) with diagnostic bands
 */
exports.getSyntheticStack = function(ccdcOutput, targetDateFractional, options) {
  options = options || {};

  var bands          = options.bands || ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'];
  var coefScale      = options.coefScale !== undefined ? options.coefScale : 10000;
  var tolYears       = options.tolYears !== undefined ? options.tolYears : 1.0;
  var coefCase       = options.coefCase || 'lower';
  var minReflectance = options.minReflectance !== undefined ? options.minReflectance : -0.05;
  var maxReflectance = options.maxReflectance !== undefined ? options.maxReflectance : 1.0;
  var fillHoles      = options.fillHoles !== undefined ? options.fillHoles : true;

  var t = ee.Number(targetDateFractional);
  var seg = exports.nearestSegment(ccdcOutput, t);

  // Clamp only linear trend term to segment start/end; harmonics use exact target date t
  var tc = ee.Image(t).max(seg.tS.arrayGet(seg.idx)).min(seg.tE.arrayGet(seg.idx));
  var omega = 2.0 * Math.PI;

  var terms = ee.Image.cat([
    ee.Image.constant(1),
    tc,
    ee.Image.constant(t.multiply(omega).cos()),
    ee.Image.constant(t.multiply(omega).sin()),
    ee.Image.constant(t.multiply(omega * 2).cos()),
    ee.Image.constant(t.multiply(omega * 2).sin()),
    ee.Image.constant(t.multiply(omega * 3).cos()),
    ee.Image.constant(t.multiply(omega * 3).sin())
  ]).toArray();

  var syntheticBands = bands.map(function(band) {
    var coefBandName = (coefCase === 'upper')
      ? band.toUpperCase() + '_coefs'
      : band.toLowerCase() + '_coefs';

    var coefs = ccdcOutput.select(coefBandName)
      .arraySlice(0, seg.idx, seg.idx.add(1))
      .arrayProject([1]);

    return coefs.multiply(terms).arrayReduce(ee.Reducer.sum(), [0])
      .arrayGet([0]).divide(coefScale).rename(band.toLowerCase());
  });

  var stack = ee.Image.cat(syntheticBands);

  // Plausibility bounds check: all bands must fall within [minReflectance, maxReflectance]
  // and pixel must be within temporal gap tolerance
  var plausible = stack.reduce(ee.Reducer.min()).gt(minReflectance)
    .and(stack.reduce(ee.Reducer.max()).lt(maxReflectance));
  var withinTol = seg.gap.lte(tolYears);
  var synthValid = plausible.and(withinTol).rename('synth_valid');

  var maskedStack = stack.updateMask(synthValid);

  if (fillHoles) {
    maskedStack = exports.fillIsolated(maskedStack, 5);
  }

  return maskedStack.addBands([seg.gap, synthValid]);
};

/**
 * Extracts first-order annual harmonic amplitude and phase for specified bands.
 *
 * @param {ee.Image}         ccdcOutput           - CCDC array output image
 * @param {number|ee.Number} targetDateFractional - Target fractional year
 * @param {Object}           [options]            - Configuration options:
 *   @param {Array<string>}    [options.bands=['nir', 'red', 'swir1']]
 *   @param {number}           [options.coefScale=10000]
 *   @param {number}           [options.tolYears=1.0]
 *   @param {string}           [options.coefCase='lower'] - 'lower' or 'upper'
 *   @param {string}           [options.outputCase='lower'] - Case for output band names ('lower' -> nir_amp, 'upper' -> NIR_amp)
 * @returns {ee.Image} Multi-band amplitude and phase image
 */
exports.getHarmonics = function(ccdcOutput, targetDateFractional, options) {
  options = options || {};

  var bands      = options.bands || ['nir', 'red', 'swir1'];
  var coefScale  = options.coefScale !== undefined ? options.coefScale : 10000;
  var tolYears   = options.tolYears !== undefined ? options.tolYears : 1.0;
  var coefCase   = options.coefCase || 'lower';
  var outputCase = options.outputCase || 'lower';

  var t = ee.Number(targetDateFractional);
  var seg = exports.nearestSegment(ccdcOutput, t);

  var harmonicImages = bands.map(function(band) {
    var coefBandName = (coefCase === 'upper')
      ? band.toUpperCase() + '_coefs'
      : band.toLowerCase() + '_coefs';

    var coefficients = ccdcOutput.select(coefBandName)
      .arraySlice(0, seg.idx, seg.idx.add(1))
      .arrayProject([1]);

    var a1 = coefficients.arrayGet([2]).divide(coefScale); // cos(wt)
    var b1 = coefficients.arrayGet([3]).divide(coefScale); // sin(wt)

    var ampName   = (outputCase === 'upper') ? band.toUpperCase() + '_amp'   : band.toLowerCase() + '_amp';
    var phaseName = (outputCase === 'upper') ? band.toUpperCase() + '_phase' : band.toLowerCase() + '_phase';

    return ee.Image.cat([
      a1.hypot(b1).rename(ampName),
      b1.atan2(a1).rename(phaseName)
    ]);
  });

  return ee.Image.cat(harmonicImages).updateMask(seg.gap.lte(tolYears));
};

/**
 * Evaluates CCDC stability between startYear and endYear for training sample transfer.
 * Criteria:
 *   1. Zero confirmed CCDC breakpoints (changeProb == 1) in [startYear, endYear].
 *   2. Valid segment coverage at both epochs within tolYears, mapping to the same segment.
 *   3. Optional gradual drift filter on linear slopes.
 *
 * @param {ee.Image} ccdcImage - Mosaicked CCDC array image
 * @param {number}   startYear - Start year (e.g. 2020.0)
 * @param {number}   endYear   - End year (e.g. 2025.99)
 * @param {Object}   [options] - { tolYears: 1.0, coefScale: 10000, coefCase: 'lower', applySlopeFilter: true, maxDrift: 0.12 }
 * @returns {ee.Image} Binary stability image (1 = stable, 0 = unstable)
 */
exports.stableInPeriod = function(ccdcImage, startYear, endYear, options) {
  options = options || {};

  var tolYears         = options.tolYears !== undefined ? options.tolYears : 1.0;
  var coefScale        = options.coefScale !== undefined ? options.coefScale : 10000;
  var coefCase         = options.coefCase || 'lower';
  var applySlopeFilter = options.applySlopeFilter !== undefined ? options.applySlopeFilter : true;
  var maxDrift         = options.maxDrift !== undefined ? options.maxDrift : 0.12;

  // Condition 1: No confirmed breakpoint in [startYear, endYear]
  var tBreak = ccdcImage.select('tBreak');
  var changeProb = ccdcImage.select('changeProb');
  var confirmedBreak = tBreak.multiply(changeProb.eq(1));
  var breaksInWindow = confirmedBreak.arrayMask(
    confirmedBreak.gte(startYear).and(confirmedBreak.lte(endYear))
  );
  var noBreak = breaksInWindow.arrayLength(0).eq(0);

  // Condition 2: Continuous single-segment coverage
  var segStart = exports.nearestSegment(ccdcImage, startYear + 0.50);
  var segEnd   = exports.nearestSegment(ccdcImage, Math.floor(endYear) + 0.50);

  var validCoverage = segStart.gap.lte(tolYears)
    .and(segEnd.gap.lte(tolYears))
    .and(segStart.idx.eq(segEnd.idx));

  var isStable = noBreak.and(validCoverage);

  // Condition 3: Gradual linear drift filter on NIR and SWIR1 slopes
  if (applySlopeFilter) {
    var deltaYears = endYear - startYear;

    var nirCoefBand = (coefCase === 'upper') ? 'NIR_coefs' : 'nir_coefs';
    var swir1CoefBand = (coefCase === 'upper') ? 'SWIR1_coefs' : 'swir1_coefs';

    var nirSlope = ccdcImage.select(nirCoefBand)
      .arraySlice(0, segEnd.idx, segEnd.idx.add(1)).arrayProject([1])
      .arrayGet([1]).divide(coefScale).abs();

    var swir1Slope = ccdcImage.select(swir1CoefBand)
      .arraySlice(0, segEnd.idx, segEnd.idx.add(1)).arrayProject([1])
      .arrayGet([1]).divide(coefScale).abs();

    var lowDrift = nirSlope.multiply(deltaYears).lte(maxDrift)
      .and(swir1Slope.multiply(deltaYears).lte(maxDrift));

    isStable = isStable.and(lowDrift);
  }

  return isStable.unmask(0).toByte();
};
