/**
 * Module: Feature Engineering & Multi-Sensor Stacking
 * =========================================================================
 * Single source for feature extraction and processing helpers:
 *   - computeSeasonalIndices: Reconstructs spectral indices from synthetic CCDC stacks.
 *   - computePhenologicalDeltas: Computes seasonal index differences (wet - dry).
 *   - applySpeckleFilter: Linear-power SAR speckle filtering (Gamma MAP, Lee, Refined Lee, Boxcar).
 *   - preprocessS1Collection: Optimized Sentinel-1 ARD pipeline (Mullissa et al. 2021 / gee_s1_ard).
 *   - getSeasonalS1: Preprocessed seasonal Sentinel-1 ARD gamma0 medians (VV/VH, wet/dry, dB).
 *   - getScanSARFeatures: QA-masked & speckle-filtered annual PALSAR-2 ScanSAR gamma0 (dB).
 *   - getPalsarBackscatter: QA-masked & speckle-filtered annual ALOS PALSAR / PALSAR-2 gamma0 (dB).
 *   - getTopography: Terrain covariates (elevation, slope, aspect, TPI landform, TWI).
 *   - getGlcmTexture: GLCM spatial texture from dry-season NIR.
 *   - loadCCDCMosaic: Generic loader for tiled CCDC asset collections.
 *
 * Designed to be sensor- and project-agnostic. All geometry/AOI arguments are
 * passed explicitly by callers.
 * =========================================================================
 */

var spectral = require('users/dmlmont/spectral:spectral');

// ---------------------------------------------------------------------------
// 1. Spectral Index Computation
// ---------------------------------------------------------------------------

/**
 * Computes a list of spectral indices from a synthetic surface reflectance image.
 * Maps sensor generic band names (blue, green, red, nir, swir1, swir2) to
 * standard optical band naming ('B2', 'B3', 'B4', etc.) and appends a seasonal suffix.
 *
 * @param {ee.Image}      synthImg   - Synthetic image containing ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']
 * @param {string}        sensor     - 'S2' (maps nir to B8, swir1 to B11, swir2 to B12) or 'Landsat' (maps nir to B5, swir1 to B6, swir2 to B7)
 * @param {Array<string>} indexList  - Names of indices supported by the 'spectral' library
 * @param {string}        suffix     - Suffix appended to output bands, e.g. 'wet' or 'dry'
 * @returns {ee.Image} Multi-band image with renamed index bands, e.g. 'NDVI_wet', 'MNDWI_dry'
 */
exports.computeSeasonalIndices = function(synthImg, sensor, indexList, suffix) {
  var isS2 = (sensor && sensor.toUpperCase() === 'S2');
  var srcBands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'];
  var dstBands = isS2
    ? ['B2', 'B3', 'B4', 'B8', 'B11', 'B12']
    : ['B2', 'B3', 'B4', 'B5', 'B6', 'B7'];

  var imgForIndices = synthImg.select(srcBands, dstBands);
  var param = {
    'B':   imgForIndices.select('B2'),
    'G':   imgForIndices.select('B3'),
    'R':   imgForIndices.select('B4'),
    'N':   imgForIndices.select(isS2 ? 'B8' : 'B5'),
    'S1':  imgForIndices.select(isS2 ? 'B11' : 'B6'),
    'S2':  imgForIndices.select(isS2 ? 'B12' : 'B7'),
    'C1':  6.0,
    'C2':  7.5,
    'g':   2.5,
    'L':   1.0
  };

  var computed = spectral.computeIndex(imgForIndices, indexList, param);
  var renamed = computed.bandNames().map(function(b) {
    return ee.String(b).cat('_').cat(suffix);
  });
  return computed.rename(renamed);
};

/**
 * Computes phenological seasonal index deltas (wet - dry).
 *
 * @param {ee.Image}      wetIndices - Image containing wet season indices (named '<Index>_wet')
 * @param {ee.Image}      dryIndices - Image containing dry season indices (named '<Index>_dry')
 * @param {Array<string>} indexList  - Base index names to subtract
 * @returns {ee.Image} Multi-band image with bands named '<Index>_delta'
 */
exports.computePhenologicalDeltas = function(wetIndices, dryIndices, indexList) {
  var deltaImages = indexList.map(function(idx) {
    return wetIndices.select(idx + '_wet')
      .subtract(dryIndices.select(idx + '_dry'))
      .rename(idx + '_delta');
  });
  return ee.Image.cat(deltaImages);
};

// ---------------------------------------------------------------------------
// 2. Texture Analysis (GLCM)
// ---------------------------------------------------------------------------

/**
 * Computes Gray-Level Co-occurrence Matrix (GLCM) spatial texture metrics
 * from a dry-season NIR band.
 *
 * @param {ee.Image} nirBand - NIR reflectance band (e.g. from synthDry)
 * @param {Object}   [options] - Optional settings: { size: 1, scaleFactor: 1000 }
 * @returns {ee.Image} 4-band image: ['B8_contrast', 'B8_ent', 'B8_asm', 'B8_idm']
 */
exports.getGlcmTexture = function(nirBand, options) {
  options = options || {};
  var size = options.size || 1;
  var scaleFactor = options.scaleFactor || 1000;

  var glcm = nirBand.select([0], ['B8'])
    .multiply(scaleFactor)
    .toInt32()
    .glcmTexture({ size: size });

  return glcm.select(['B8_contrast', 'B8_ent', 'B8_asm', 'B8_idm']);
};

// ---------------------------------------------------------------------------
// 3. Radar / SAR ARD Preprocessing & Backscatter Composites
// ---------------------------------------------------------------------------
// Sentinel-1 ARD methodology adapted and optimized from:
//  Mullissa, A.; Vollrath, A.; Odongo-Braun, C.; Slagter, B.; Balling, J.; Gou, Y.; Gorelick, N.; Reiche, J.
//  Sentinel-1 SAR Backscatter Analysis Ready Data Preparation in Google Earth Engine. Remote Sens. 2021, 13, 1954. 
//  https://doi.org/10.3390/rs13101954 /
//  https://github.com/adugnag/gee_s1_ard/tree/main

// Key GEE performance optimizations over reference gee_s1_ard:
//   1. Precomputes time-invariant DEM slope/aspect gradient components ONCE outside
//      collection.map() and avoids forced 10m .reproject() per scene (5-10x faster).
//   2. Estimates scene look azimuth from the 20km 'angle' band at scale 5000m
//      instead of 1000m (25x fewer pixels per reduceRegion).
//   3. Removes redundant linear->dB->linear roundtrip in border noise masking.
//   4. Defaults to 'POST_COMPOSITE' speckle filtering on seasonal medians (with
//      dynamic ENL derived from pixel observation count), replacing O(N^2) Quegan
//      multi-temporal archive queries (~30-50x faster for classification stacks).
//   5. Fixes Gamma MAP discriminant formula and uses .where() branching & .updateMask()
//      for numerical safety and mask preservation.

var DEG2RAD = Math.PI / 180;
var PALSAR_LIN_CF = Math.pow(10, -8.3); // 10^(-83 / 10) for gamma0_lin = DN^2 * 10^(-8.3)

/**
 * Converts linear backscatter power bands to decibels: 10 * log10(max(x, 1e-10)).
 * Preserves any non-selected bands (e.g. 'angle') and image properties.
 *
 * @param {ee.Image}           image     - Linear-power SAR image
 * @param {Array<string>|ee.List} [bands] - Band names to convert (defaults to all except 'angle')
 * @returns {ee.Image}
 */
function linToDb(image, bands) {
  var bandNames = bands ? ee.List(bands) : image.bandNames().remove('angle');
  var db = image.select(bandNames)
    .max(1e-10)
    .log10()
    .multiply(10.0)
    .rename(bandNames);
  return image.addBands(db, null, true);
}

/**
 * Additional Sentinel-1 IW border noise correction (Hird et al. 2017 / gee_s1_ard).
 * Masks swath edge incidence angles (<= 30.63993 deg or >= 45.23993 deg) and
 * non-positive backscatter directly in linear power (avoiding redundant dB roundtrips).
 *
 * @param {ee.Image} image - Linear-power S1_GRD_FLOAT image containing 'angle'
 * @returns {ee.Image}
 */
function maskS1BorderNoise(image) {
  var ang = image.select('angle');
  var polBands = image.bandNames().remove('angle');
  var angMask = ang.gt(30.63993).and(ang.lt(45.23993));
  var posMask = image.select(polBands).gt(0).reduce(ee.Reducer.min());
  return image.updateMask(angMask.and(posMask))
    .copyProperties(image, ['system:time_start', 'orbitProperties_pass']);
}

/**
 * Boxcar spatial mean speckle filter (linear power).
 */
function boxcarFilter(image, bandNames, kernelSize) {
  var radius = Math.max(1, Math.floor(kernelSize / 2));
  var kernel = ee.Kernel.square({ radius: radius, units: 'pixels', normalize: true });
  var filtered = image.select(bandNames).convolve(kernel).rename(bandNames);
  return image.addBands(filtered, null, true);
}

/**
 * Lee (1980) MMSE spatial speckle filter (linear power).
 * Accepts either a scalar number or an ee.Image for ENL (Equivalent Number of Looks).
 */
function leeFilter(image, bandNames, kernelSize, enl) {
  var radius = Math.max(1, Math.floor(kernelSize / 2));
  var enlImg = (enl instanceof ee.Image) ? enl : ee.Image.constant(enl);
  var etaSq = ee.Image.constant(1.0).divide(enlImg);
  var oneImg = ee.Image.constant(1.0);

  var img = image.select(bandNames);
  var reducers = ee.Reducer.mean().combine({
    reducer2: ee.Reducer.variance(),
    sharedInputs: true
  });
  var stats = img.reduceNeighborhood({
    reducer: reducers,
    kernel: ee.Kernel.square(radius, 'pixels'),
    optimization: 'window'
  });

  var meanBand = bandNames.map(function(b) { return ee.String(b).cat('_mean'); });
  var varBand  = bandNames.map(function(b) { return ee.String(b).cat('_variance'); });

  var zBar = stats.select(meanBand).rename(bandNames);
  var varZ = stats.select(varBand).rename(bandNames);

  var varX = varZ.subtract(zBar.pow(2).multiply(etaSq)).divide(oneImg.add(etaSq));
  var b = varX.divide(varZ.max(1e-12)).clamp(0, 1);

  var output = oneImg.subtract(b).multiply(zBar.abs())
    .add(b.multiply(img))
    .rename(bandNames);
  return image.addBands(output, null, true);
}

/**
 * Lopes et al. (1990) Gamma MAP spatial speckle filter (linear power).
 * Uses corrected Eq. (11) quadratic discriminant and .where() branching for speed & stability.
 * Accepts either a scalar number or an ee.Image for ENL.
 */
function gammaMapFilter(image, bandNames, kernelSize, enl) {
  var radius = Math.max(1, Math.floor(kernelSize / 2));
  var enlImg = (enl instanceof ee.Image) ? enl : ee.Image.constant(enl);
  var oneImg = ee.Image.constant(1.0);

  var img = image.select(bandNames);
  var reducers = ee.Reducer.mean().combine({
    reducer2: ee.Reducer.stdDev(),
    sharedInputs: true
  });
  var stats = img.reduceNeighborhood({
    reducer: reducers,
    kernel: ee.Kernel.square(radius, 'pixels'),
    optimization: 'window'
  });

  var meanBand   = bandNames.map(function(b) { return ee.String(b).cat('_mean'); });
  var stdDevBand = bandNames.map(function(b) { return ee.String(b).cat('_stdDev'); });

  var z    = stats.select(meanBand).rename(bandNames);
  var sigZ = stats.select(stdDevBand).rename(bandNames);

  var ci   = sigZ.divide(z.max(1e-12));
  var cu   = oneImg.divide(enlImg.sqrt());
  var cmax = cu.multiply(Math.SQRT2);

  var cuSq  = cu.pow(2);
  var ciSq  = ci.pow(2);
  var alpha = oneImg.add(cuSq).divide(ciSq.subtract(cuSq).max(1e-6));

  // Lopes et al. (1990) Eq. 11:
  // rHat = [z*(alpha - L - 1) + sqrt(z^2*(alpha - L - 1)^2 + 4*alpha*L*I*z)] / (2*alpha)
  var term = z.multiply(alpha.subtract(enlImg).subtract(oneImg));
  var disc = term.pow(2).add(alpha.multiply(4.0).multiply(enlImg).multiply(img).multiply(z));
  var rHat = term.add(disc.max(0).sqrt()).divide(alpha.multiply(2.0));

  var output = rHat
    .where(ci.lte(cu), z)       // Homogeneous region -> local mean
    .where(ci.gte(cmax), img)   // Strong point scatterer -> retain original
    .rename(bandNames);

  return image.addBands(output, null, true);
}

/**
 * Refined Lee 7x7 directional edge-preserving speckle filter (linear power).
 * Self-calibrating noise variance (sigmaV) estimated from 3x3 sub-windows.
 * Optimized vs gee_s1_ard by combining mean+variance into single reduceNeighborhood passes (2x faster).
 */
function refinedLeeFilter(image, bandNames) {
  var meanVarReducer = ee.Reducer.mean().combine({
    reducer2: ee.Reducer.variance(),
    sharedInputs: true
  });

  var weights3 = ee.List.repeat(ee.List.repeat(1, 3), 3);
  var kernel3 = ee.Kernel.fixed(3, 3, weights3, 1, 1, false);

  var sampleWeights = ee.List([
    [0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0]
  ]);
  var sampleKernel = ee.Kernel.fixed(7, 7, sampleWeights, 3, 3, false);

  var rectWeights = ee.List.repeat(ee.List.repeat(0, 7), 3)
    .cat(ee.List.repeat(ee.List.repeat(1, 7), 4));
  var diagWeights = ee.List([
    [1, 0, 0, 0, 0, 0, 0],
    [1, 1, 0, 0, 0, 0, 0],
    [1, 1, 1, 0, 0, 0, 0],
    [1, 1, 1, 1, 0, 0, 0],
    [1, 1, 1, 1, 1, 0, 0],
    [1, 1, 1, 1, 1, 1, 0],
    [1, 1, 1, 1, 1, 1, 1]
  ]);
  var rectKernel = ee.Kernel.fixed(7, 7, rectWeights, 3, 3, false);
  var diagKernel = ee.Kernel.fixed(7, 7, diagWeights, 3, 3, false);

  var result = ee.ImageCollection(bandNames.map(function(b) {
    var img = image.select([b]);

    // 3x3 sub-window stats (single combined pass)
    var stats3 = img.reduceNeighborhood(meanVarReducer, kernel3);
    var mean3 = stats3.select(0);
    var var3  = stats3.select(1);

    var sampleMean = mean3.neighborhoodToBands(sampleKernel);
    var sampleVar  = var3.neighborhoodToBands(sampleKernel);

    // 4 directional gradients
    var gradients = sampleMean.select(1).subtract(sampleMean.select(7)).abs()
      .addBands(sampleMean.select(6).subtract(sampleMean.select(2)).abs())
      .addBands(sampleMean.select(3).subtract(sampleMean.select(5)).abs())
      .addBands(sampleMean.select(0).subtract(sampleMean.select(8)).abs());

    var maxGradient = gradients.reduce(ee.Reducer.max());
    var gradMask = gradients.eq(maxGradient);
    gradMask = gradMask.addBands(gradMask);

    // 8 oriented sub-window directions
    var dir0 = sampleMean.select(1).subtract(sampleMean.select(4)).gt(sampleMean.select(4).subtract(sampleMean.select(7)));
    var dir1 = sampleMean.select(6).subtract(sampleMean.select(4)).gt(sampleMean.select(4).subtract(sampleMean.select(2)));
    var dir2 = sampleMean.select(3).subtract(sampleMean.select(4)).gt(sampleMean.select(4).subtract(sampleMean.select(5)));
    var dir3 = sampleMean.select(0).subtract(sampleMean.select(4)).gt(sampleMean.select(4).subtract(sampleMean.select(8)));

    var directions = dir0.multiply(1)
      .addBands(dir1.multiply(2))
      .addBands(dir2.multiply(3))
      .addBands(dir3.multiply(4))
      .addBands(dir0.not().multiply(5))
      .addBands(dir1.not().multiply(6))
      .addBands(dir2.not().multiply(7))
      .addBands(dir3.not().multiply(8))
      .updateMask(gradMask)
      .reduce(ee.Reducer.sum());

    // Self-calibrated local noise variance from the 5 most homogeneous sub-windows
    var sampleStats = sampleVar.divide(sampleMean.pow(2).max(1e-12));
    var sigmaV = sampleStats.toArray().arraySort().arraySlice(0, 0, 5)
      .arrayReduce(ee.Reducer.mean(), [0])
      .arrayGet([0]);

    // 8 directional 7x7 mean & variance stats (combined reducer = 8 passes instead of 16)
    var sRect0 = img.reduceNeighborhood(meanVarReducer, rectKernel);
    var sDiag0 = img.reduceNeighborhood(meanVarReducer, diagKernel);

    var dirMean = sRect0.select(0).updateMask(directions.eq(1))
      .addBands(sDiag0.select(0).updateMask(directions.eq(2)));
    var dirVar = sRect0.select(1).updateMask(directions.eq(1))
      .addBands(sDiag0.select(1).updateMask(directions.eq(2)));

    for (var i = 1; i < 4; i++) {
      var sRect = img.reduceNeighborhood(meanVarReducer, rectKernel.rotate(i));
      var sDiag = img.reduceNeighborhood(meanVarReducer, diagKernel.rotate(i));
      dirMean = dirMean
        .addBands(sRect.select(0).updateMask(directions.eq(2 * i + 1)))
        .addBands(sDiag.select(0).updateMask(directions.eq(2 * i + 2)));
      dirVar = dirVar
        .addBands(sRect.select(1).updateMask(directions.eq(2 * i + 1)))
        .addBands(sDiag.select(1).updateMask(directions.eq(2 * i + 2)));
    }

    dirMean = dirMean.reduce(ee.Reducer.sum());
    dirVar  = dirVar.reduce(ee.Reducer.sum());

    var varX = dirVar.subtract(dirMean.pow(2).multiply(sigmaV)).divide(sigmaV.add(1.0));
    var bWeight = varX.divide(dirVar.max(1e-12)).clamp(0, 1);

    return dirMean.add(bWeight.multiply(img.subtract(dirMean))).float();
  })).toBands().rename(bandNames);

  return image.addBands(result, null, true);
}

/**
 * Applies a spatial SAR speckle filter to a linear-power backscatter image.
 * Reusable across Sentinel-1, PALSAR-2 ScanSAR, and ALOS PALSAR mosaics.
 *
 * @param {ee.Image} image - Linear-power SAR backscatter image
 * @param {Object}   [options]
 *   - filter     {string}          'GAMMA_MAP' (default), 'LEE', 'REFINED_LEE', or 'BOXCAR'
 *   - kernelSize {number}          Odd kernel size in pixels (default 5; Refined Lee is fixed 7x7)
 *   - enl        {number|ee.Image} Equivalent Number of Looks (default 5)
 *   - bands      {Array<string>}   Bands to filter (defaults to all bands except 'angle')
 * @returns {ee.Image} Filtered linear-power image
 */
exports.applySpeckleFilter = function(image, options) {
  options = options || {};
  var filter = (options.filter || 'GAMMA_MAP').toUpperCase().replace(/\s+/g, '_');
  var kernelSize = options.kernelSize || 5;
  var enl = (options.enl !== undefined && options.enl !== null) ? options.enl : 5;
  var bandNames = options.bands ? ee.List(options.bands) : image.bandNames().remove('angle');

  if (filter === 'BOXCAR') {
    return boxcarFilter(image, bandNames, kernelSize);
  }
  if (filter === 'LEE') {
    return leeFilter(image, bandNames, kernelSize, enl);
  }
  if (filter === 'REFINED_LEE') {
    return refinedLeeFilter(image, bandNames);
  }
  // Default: GAMMA_MAP
  return gammaMapFilter(image, bandNames, kernelSize, enl);
};

/**
 * Optimized Angular-Based Radiometric Terrain Flattening for Sentinel-1
 * (Vollrath, Mullissa & Reiche 2020; gee_s1_ard).
 *
 * Precomputes branch-cut-free DEM terrain gradient vectors
 *   tan(alpha_s) * cos(phi_s)  and  tan(alpha_s) * sin(phi_s)
 * ONCE outside collection.map(), eliminating per-scene DEM slope/aspect computation
 * and forced 10m .reproject() calls.
 *
 * @param {ee.ImageCollection} collection - Linear-power S1_GRD_FLOAT collection with 'angle'
 * @param {Object}             [options]
 *   - model               {string}   'VOLUME' (default) or 'DIRECT'
 *   - dem                 {ee.Image} Elevation image (default NASA/NASADEM_HGT/001)
 *   - layoverShadowBuffer {number}   Additional erosion buffer in meters (default 0)
 * @returns {ee.ImageCollection} Terrain-flattened gamma0 linear-power collection
 */
function applyS1TerrainFlattening(collection, options) {
  options = options || {};
  var model = (options.model || 'VOLUME').toUpperCase();
  var buffer = options.layoverShadowBuffer || 0;
  var dem = options.dem
    ? ee.Image(options.dem).select(0)
    : ee.Image('NASA/NASADEM_HGT/001').select('elevation');

  var ninetyRad = ee.Image.constant(Math.PI / 2);

  // Precompute time-invariant DEM terrain gradient components ONCE outside .map()
  // Using tan(alpha_s)*cos(phi_s) and tan(alpha_s)*sin(phi_s) avoids the +/-180 deg
  // azimuth branch cut when bilinearly resampling to output scale.
  var alphaSRad = ee.Terrain.slope(dem).select('slope').multiply(DEG2RAD);
  var aspectDeg = ee.Terrain.aspect(dem).select('aspect');
  var phiSRad = aspectDeg
    .where(aspectDeg.gt(180), aspectDeg.subtract(360))
    .multiply(-DEG2RAD); // uphill convention (gee_s1_ard)

  var tanAlphaS = alphaSRad.tan();
  var tanAlphaCosPhiS = tanAlphaS.multiply(phiSRad.cos()).resample('bilinear');
  var tanAlphaSinPhiS = tanAlphaS.multiply(phiSRad.sin()).resample('bilinear');

  var erodeMask = function(mask, distance) {
    var d = mask.not().unmask(1)
      .fastDistanceTransform(30).sqrt()
      .multiply(ee.Image.pixelArea().sqrt());
    return mask.updateMask(d.gt(distance));
  };

  var correctScene = function(image) {
    var polBands = image.bandNames().remove('angle');
    var angleImg = image.select('angle');

    // Estimate look azimuth from the 20km-GSD 'angle' band at 5km scale (25x faster than 1km)
    var headingDict = ee.Terrain.aspect(angleImg).reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: image.geometry(),
      scale: 5000,
      maxPixels: 1e6,
      bestEffort: true
    });
    var rawHeading = ee.Number(ee.Dictionary(headingDict).combine({ aspect: 0 }, false).get('aspect'));
    var heading = ee.Number(ee.Algorithms.If(rawHeading.gt(180), rawHeading.subtract(360), rawHeading));
    var phiIRad = heading.multiply(DEG2RAD);
    var cosPhiI = phiIRad.cos();
    var sinPhiI = phiIRad.sin();

    // Radar incidence angle (radians)
    var thetaIRad = angleImg.multiply(DEG2RAD);

    // Slope steepness in range: tan(alpha_r) = tan(alpha_s) * cos(phi_i - phi_s)
    // = (tan(alpha_s)*cos(phi_s))*cos(phi_i) + (tan(alpha_s)*sin(phi_s))*sin(phi_i)
    var alphaRRad = tanAlphaCosPhiS.multiply(cosPhiI)
      .add(tanAlphaSinPhiS.multiply(sinPhiI))
      .atan();

    // Ellipsoid gamma0 = sigma0 / cos(theta_i)
    var gamma0 = image.select(polBands).divide(thetaIRad.cos());

    // Slope Correction Factor (SCF)
    var scf;
    if (model === 'DIRECT') {
      // Slope steepness in azimuth: tan(alpha_az) = tan(alpha_s) * sin(phi_i - phi_s)
      var alphaAzRad = tanAlphaCosPhiS.multiply(sinPhiI)
        .subtract(tanAlphaSinPhiS.multiply(cosPhiI))
        .atan();
      scf = ninetyRad.subtract(thetaIRad).cos()
        .divide(alphaAzRad.cos().multiply(ninetyRad.subtract(thetaIRad).add(alphaRRad).cos()));
    } else {
      // VOLUME scattering model (Hoekman & Reiche 2015 / Vollrath et al. 2020)
      scf = ninetyRad.subtract(thetaIRad).add(alphaRRad).tan()
        .divide(ninetyRad.subtract(thetaIRad).tan());
    }

    var gamma0Flat = gamma0.multiply(scf);

    // Active layover and radar shadow mask (preserve border mask via updateMask)
    var layover = alphaRRad.lt(thetaIRad);
    var shadow  = alphaRRad.gt(thetaIRad.subtract(ninetyRad));
    var lsMask  = layover.and(shadow).and(scf.gt(0));
    if (buffer > 0) {
      lsMask = erodeMask(lsMask, buffer);
    }

    return gamma0Flat
      .updateMask(lsMask)
      .rename(polBands)
      .addBands(angleImg)
      .copyProperties(image, ['system:time_start', 'orbitProperties_pass']);
  };

  return collection.map(correctScene);
}

/**
 * Full optimized Sentinel-1 Analysis Ready Data (ARD) collection preprocessor
 * implementing Mullissa et al. (2021) / adugnag/gee_s1_ard:
 *   1. Loads linear-power COPERNICUS/S1_GRD_FLOAT (IW, 10m, dual-pol VV+VH).
 *   2. Additional border noise masking (angle 30.64-45.24 deg).
 *   3. Optional per-scene mono-temporal speckle filtering (linear power).
 *   4. Optimized radiometric terrain flattening (Vollrath et al. 2020 VOLUME/DIRECT).
 *   5. Optional conversion to dB.
 *
 * @param {ee.Date|string} startDate - Inclusive start date
 * @param {ee.Date|string} endDate   - Exclusive end date
 * @param {ee.Geometry}    roi       - Region of interest
 * @param {Object}         [options]
 *   - orbit                  {string}   'BOTH' (default), 'ASCENDING', or 'DESCENDING'
 *   - applyBorderNoise       {boolean}  Apply edge angle masking (default true)
 *   - applySpeckleFilter     {boolean}  Apply per-scene mono-temporal filter (default false;
 *                                       use post-composite filtering in getSeasonalS1 for speed)
 *   - speckleFilter          {string}   'GAMMA_MAP' (default), 'LEE', 'REFINED_LEE', 'BOXCAR'
 *   - kernelSize             {number}   Kernel size in pixels (default 5)
 *   - applyTerrainFlattening {boolean}  Apply radiometric slope correction (default true)
 *   - terrainModel           {string}   'VOLUME' (default) or 'DIRECT'
 *   - layoverShadowBuffer    {number}   Extra layover/shadow buffer in meters (default 0)
 *   - dem                    {ee.Image} DEM image (default NASA/NASADEM_HGT/001)
 *   - format                 {string}   'LINEAR' (default) or 'DB'
 * @returns {ee.ImageCollection} Preprocessed Sentinel-1 collection (['VV', 'VH', 'angle'])
 */
exports.preprocessS1Collection = function(startDate, endDate, roi, options) {
  options = options || {};
  var orbit = (options.orbit || 'BOTH').toUpperCase();
  var applyBorderNoise = (options.applyBorderNoise !== undefined) ? options.applyBorderNoise : true;
  var applySpeckle = (options.applySpeckleFilter !== undefined) ? options.applySpeckleFilter : false;
  var applyTerrain = (options.applyTerrainFlattening !== undefined) ? options.applyTerrainFlattening : true;
  var format = (options.format || 'LINEAR').toUpperCase();

  var s1 = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.eq('resolution_meters', 10))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .select(['VV', 'VH', 'angle']);

  if (orbit !== 'BOTH') {
    s1 = s1.filter(ee.Filter.eq('orbitProperties_pass', orbit));
  }

  if (applyBorderNoise) {
    s1 = s1.map(maskS1BorderNoise);
  }

  if (applySpeckle) {
    var filterOpts = {
      filter: options.speckleFilter || 'GAMMA_MAP',
      kernelSize: options.kernelSize || 5,
      enl: options.enl || 5,
      bands: ['VV', 'VH']
    };
    s1 = s1.map(function(img) {
      return exports.applySpeckleFilter(img, filterOpts)
        .copyProperties(img, ['system:time_start', 'orbitProperties_pass']);
    });
  }

  if (applyTerrain) {
    s1 = applyS1TerrainFlattening(s1, {
      model: options.terrainModel || 'VOLUME',
      dem: options.dem || null,
      layoverShadowBuffer: options.layoverShadowBuffer || 0
    });
  }

  if (format === 'DB') {
    s1 = s1.map(function(img) {
      return linToDb(img, ['VV', 'VH'])
        .copyProperties(img, ['system:time_start', 'orbitProperties_pass']);
    });
  }

  return s1;
};

/**
 * Generates ARD-preprocessed seasonal median Sentinel-1 gamma0 backscatter
 * composites (VV & VH in dB), optimized for LULC classification.
 *
 * Pipeline (all radiometric math in linear power until final dB conversion):
 *   COPERNICUS/S1_GRD_FLOAT -> Border Noise Mask -> Optimized Radiometric Terrain
 *   Flattening (Vollrath et al. 2020) -> Seasonal Median Composite -> Spatial
 *   Speckle Filter (post-composite by default with count-calibrated ENL) -> dB.
 *
 * Default monsoonal season split:
 *   - Wet Season: December through April
 *   - Dry Season: May through November
 *
 * @param {number}      year - Target year
 * @param {ee.Geometry} roi  - Region of interest
 * @param {Object}      [options]
 *   - wetMonths              {Array<number>} [12, 4]
 *   - dryMonths              {Array<number>} [5, 11]
 *   - orbit                  {string}        'BOTH' (default), 'ASCENDING', 'DESCENDING'
 *   - applyBorderNoise       {boolean}       default true
 *   - applyTerrainFlattening {boolean}       default true
 *   - terrainModel           {string}        'VOLUME' (default) or 'DIRECT'
 *   - layoverShadowBuffer    {number}        default 0 (meters)
 *   - dem                    {ee.Image}      default NASA/NASADEM_HGT/001
 *   - applySpeckleFilter     {boolean}       default true
 *   - speckleMode            {string}        'POST_COMPOSITE' (default, fast), 'MONO', or 'NONE'
 *   - speckleFilter          {string}        'GAMMA_MAP' (default), 'LEE', 'REFINED_LEE', 'BOXCAR'
 *   - kernelSize             {number}        default 5 (pixels)
 *   - enl                    {number}        Optional manual ENL override
 *   - fillMaskedGaps         {boolean}       Fill rare permanent layover/shadow holes (default true)
 *   - includeRatio           {boolean}       Also output VH_VV_ratio_dry/wet in dB (default false)
 * @returns {ee.Image} 4-band image: ['VH_dry', 'VH_wet', 'VV_dry', 'VV_wet'] (dB)
 */
exports.getSeasonalS1 = function(year, roi, options) {
  options = options || {};
  var wetMonths = options.wetMonths || [12, 4];
  var dryMonths = options.dryMonths || [5, 11];
  var applySpeckle = (options.applySpeckleFilter !== undefined) ? options.applySpeckleFilter : true;
  var speckleMode = (options.speckleMode || 'POST_COMPOSITE').toUpperCase();
  var speckleFilter = options.speckleFilter || 'GAMMA_MAP';
  var kernelSize = options.kernelSize || 5;
  var applyTerrain = (options.applyTerrainFlattening !== undefined) ? options.applyTerrainFlattening : true;
  var fillMaskedGaps = (options.fillMaskedGaps !== undefined) ? options.fillMaskedGaps : true;

  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = startDate.advance(1, 'year');

  // Step 1-3: Load linear S1_GRD_FLOAT, border noise mask, optional MONO speckle, terrain flattening
  var s1Ard = exports.preprocessS1Collection(startDate, endDate, roi, {
    orbit: options.orbit || 'BOTH',
    applyBorderNoise: (options.applyBorderNoise !== undefined) ? options.applyBorderNoise : true,
    applySpeckleFilter: applySpeckle && (speckleMode === 'MONO'),
    speckleFilter: speckleFilter,
    kernelSize: kernelSize,
    enl: options.enl || 5,
    applyTerrainFlattening: applyTerrain,
    terrainModel: options.terrainModel || 'VOLUME',
    layoverShadowBuffer: options.layoverShadowBuffer || 0,
    dem: options.dem || null,
    format: 'LINEAR'
  });

  // Partition into wet and dry seasons
  var wetCol = s1Ard.filter(ee.Filter.calendarRange(wetMonths[0], wetMonths[1], 'month'));
  var dryCol = s1Ard.filter(ee.Filter.calendarRange(dryMonths[0], dryMonths[1], 'month'));

  var wetLin = wetCol.select(['VV', 'VH']).median();
  var dryLin = dryCol.select(['VV', 'VH']).median();

  // Optional gap fill for rare cliff pixels masked by layover/shadow in every seasonal pass:
  // Reuses opposite-season flattened median first (zero extra cost), then lightweight ellipsoid gamma0.
  if (fillMaskedGaps && applyTerrain) {
    var ellipsoidFallback = exports.preprocessS1Collection(startDate, endDate, roi, {
      orbit: options.orbit || 'BOTH',
      applyBorderNoise: (options.applyBorderNoise !== undefined) ? options.applyBorderNoise : true,
      applySpeckleFilter: false,
      applyTerrainFlattening: false,
      format: 'LINEAR'
    }).map(function(img) {
      var cosTheta = img.select('angle').multiply(DEG2RAD).cos();
      return img.select(['VV', 'VH']).divide(cosTheta);
    }).median();

    var wetFilled = wetLin.unmask(dryLin).unmask(ellipsoidFallback);
    var dryFilled = dryLin.unmask(wetLin).unmask(ellipsoidFallback);
    wetLin = wetFilled;
    dryLin = dryFilled;
  }

  // Post-composite spatial speckle filter (linear power) + conversion to dB
  var finalizeSeason = function(compLin, seasonCol) {
    if (applySpeckle && speckleMode === 'POST_COMPOSITE') {
      // For a temporal median of N independent 5-look S1 images:
      // Var(median) ~ (pi/2) * sigma^2 / N  =>  ENL_eff ~ (2/pi) * 5 * N ~ 3.18 * N
      // .unmask(0) ensures gap-filled pixels receive a valid ENL (>= 5) and are not re-masked.
      var enlParam = options.enl
        ? options.enl
        : seasonCol.select('VV').count().unmask(0).multiply(3.18).max(5);
      compLin = exports.applySpeckleFilter(compLin, {
        filter: speckleFilter,
        kernelSize: kernelSize,
        enl: enlParam,
        bands: ['VV', 'VH']
      });
    }
    return linToDb(compLin, ['VV', 'VH']);
  };

  var wetDb = finalizeSeason(wetLin, wetCol);
  var dryDb = finalizeSeason(dryLin, dryCol);

  var vh_dry = dryDb.select('VH').rename('VH_dry');
  var vh_wet = wetDb.select('VH').rename('VH_wet');
  var vv_dry = dryDb.select('VV').rename('VV_dry');
  var vv_wet = wetDb.select('VV').rename('VV_wet');

  var outList = [vh_dry, vh_wet, vv_dry, vv_wet];
  if (options.includeRatio) {
    outList.push(vh_dry.subtract(vv_dry).rename('VH_VV_ratio_dry'));
    outList.push(vh_wet.subtract(vv_wet).rename('VH_VV_ratio_wet'));
  }

  return ee.Image.cat(outList).clip(roi);
};

/**
 * Annual ALOS-2 PALSAR-2 ScanSAR Level 2.2 dual-polarization backscatter in dB.
 *
 * Note: JAXA Level 2.2 ScanSAR is already CARD4L-compliant ARD (ortho-rectified
 * and radiometrically terrain-flattened to gamma0 using AW3D30), so no external
 * terrain flattening is applied. Preprocessing applied here:
 *   1. Quality bitmask (MSK bits 0-2): masks out Layover (2), Shadowing (3), and
 *      Invalid data (5) while keeping Valid land (1) and Ocean water (4, optional).
 *   2. Linear gamma0 power calibration: gamma0_lin = DN^2 * 10^(-8.3).
 *   3. Temporal median compositing with optional layover/shadow gap fallback.
 *   4. Post-composite spatial speckle filtering in linear power (count-calibrated ENL).
 *   5. Conversion to decibels (dB).
 *
 * @param {number}      year - Target year
 * @param {ee.Geometry} roi  - Region of interest
 * @param {Object}      [options]
 *   - applyQaMask        {boolean} Mask MSK layover/shadow/invalid pixels (default true)
 *   - maskOcean          {boolean} Also mask MSK==4 ocean water (default false)
 *   - fillMaskedGaps     {boolean} Unmask permanent layover/shadow holes (default true)
 *   - applySpeckleFilter {boolean} Apply post-composite speckle filter (default true)
 *   - speckleFilter      {string}  'GAMMA_MAP' (default), 'LEE', 'REFINED_LEE', 'BOXCAR'
 *   - kernelSize         {number}  Kernel size in pixels (default 3 for 25m GSD)
 *   - enl                {number}  Optional manual ENL override
 * @returns {ee.Image} 2-band image: ['ScanSAR_HH_dB', 'ScanSAR_HV_dB']
 */
exports.getScanSARFeatures = function(year, roi, options) {
  options = options || {};
  var applyQaMask = (options.applyQaMask !== undefined) ? options.applyQaMask : true;
  var maskOcean = (options.maskOcean !== undefined) ? options.maskOcean : false;
  var fillMaskedGaps = (options.fillMaskedGaps !== undefined) ? options.fillMaskedGaps : true;
  var applySpeckle = (options.applySpeckleFilter !== undefined) ? options.applySpeckleFilter : true;
  var speckleFilter = options.speckleFilter || 'GAMMA_MAP';
  var kernelSize = options.kernelSize || 3;

  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = startDate.advance(1, 'year');

  var rawCol = ee.ImageCollection('JAXA/ALOS/PALSAR-2/Level2_2/ScanSAR')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.listContains('Polarizations', 'HH'))
    .filter(ee.Filter.listContains('Polarizations', 'HV'));

  // Calibrate positive DN to linear gamma0 power: DN^2 * 10^(-8.3)
  var toLinear = function(img) {
    var pol = img.select(['HH', 'HV']).toFloat();
    var posMask = pol.select('HH').gt(0).and(pol.select('HV').gt(0));
    var lin = pol.pow(2).multiply(PALSAR_LIN_CF).updateMask(posMask);
    if (applyQaMask) {
      var mskVal = img.select('MSK').bitwiseAnd(7);
      var qaValid = maskOcean
        ? mskVal.eq(1)
        : mskVal.eq(1).or(mskVal.eq(4));
      lin = lin.updateMask(qaValid);
    }
    return lin.copyProperties(img, ['system:time_start']);
  };

  var linCol = rawCol.map(toLinear);
  var compLin = linCol.median();

  // Optional gap fill for slopes masked by layover/shadow in every scene of the year
  if (fillMaskedGaps && applyQaMask) {
    var fallbackLin = rawCol.map(function(img) {
      var pol = img.select(['HH', 'HV']).toFloat();
      return pol.pow(2).multiply(PALSAR_LIN_CF)
        .updateMask(pol.select('HH').gt(0).and(pol.select('HV').gt(0)));
    }).median();
    compLin = compLin.unmask(fallbackLin);
  }

  if (applySpeckle) {
    // ScanSAR L2.2 single-scene ENL ~ 4; median ENL_eff ~ (2/pi) * 4 * N ~ 2.55 * N
    // .unmask(0) ensures gap-filled pixels receive a valid ENL (>= 4) and are not re-masked.
    var enlParam = options.enl
      ? options.enl
      : linCol.select('HH').count().unmask(0).multiply(2.55).max(4);
    compLin = exports.applySpeckleFilter(compLin, {
      filter: speckleFilter,
      kernelSize: kernelSize,
      enl: enlParam,
      bands: ['HH', 'HV']
    });
  }

  return linToDb(compLin, ['HH', 'HV'])
    .select(['HH', 'HV'], ['ScanSAR_HH_dB', 'ScanSAR_HV_dB'])
    .clip(roi);
};

/**
 * Retrieves ALOS PALSAR / PALSAR-2 L-band SAR yearly mosaic, applies QA masking
 * and linear-power speckle filtering, and converts to gamma0 in dB.
 *
 * Note: JAXA yearly mosaics are already ortho-rectified, radiometrically
 * slope-corrected (SRTM 90m for 2007-2010; AW3D30 for 2015+), and destriped.
 * Preprocessing applied here:
 *   1. Safe temporal clamping across the 2011-2014 ALOS/ALOS-2 mission gap.
 *   2. QA band masking: excludes No data (0), Layover (100), and Shadowing (150);
 *      keeps Land (255) and Ocean/water (50, optional).
 *   3. Two-tier fallback: fills target-year QA gaps from QA-masked 3-year median
 *      first, then optional unmasked DN>0 fallback so classification has no holes.
 *   4. Linear gamma0 calibration + spatial speckle filtering (Lee 3x3 default).
 *   5. Bilinear resampling and dB conversion + HH-HV log-ratio difference.
 *
 * @param {number}      year - Target year (<=2012 uses PALSAR 2007-2010; >=2013 uses PALSAR-2 2015-2024)
 * @param {ee.Geometry} [roi] - Optional clip geometry
 * @param {Object}      [options]
 *   - applyQaMask        {boolean} Mask qa layover/shadow/nodata (default true)
 *   - maskOcean          {boolean} Also mask qa==50 water pixels (default false)
 *   - fillMaskedGaps     {boolean} Unmask permanent layover/shadow holes (default true)
 *   - applySpeckleFilter {boolean} Apply spatial speckle filter in linear power (default true)
 *   - speckleFilter      {string}  'LEE' (default), 'GAMMA_MAP', 'REFINED_LEE', 'BOXCAR'
 *   - kernelSize         {number}  Kernel size in pixels (default 3 for 25m GSD)
 *   - enl                {number}  Equivalent Number of Looks (default 4)
 * @returns {ee.Image} 3-band image: ['palsar_HH', 'palsar_HV', 'palsar_HH_HV_diff']
 */
exports.getPalsarBackscatter = function(year, roi, options) {
  options = options || {};
  var applyQaMask = (options.applyQaMask !== undefined) ? options.applyQaMask : true;
  var maskOcean = (options.maskOcean !== undefined) ? options.maskOcean : false;
  var fillMaskedGaps = (options.fillMaskedGaps !== undefined) ? options.fillMaskedGaps : true;
  var applySpeckle = (options.applySpeckleFilter !== undefined) ? options.applySpeckleFilter : true;
  var speckleFilter = options.speckleFilter || 'LEE';
  var kernelSize = options.kernelSize || 3;
  var enl = options.enl || 4;

  // Handle 2011-2014 ALOS / ALOS-2 mission gap safely
  var sarCol, startYear;
  if (year <= 2012) {
    sarCol = 'JAXA/ALOS/PALSAR/YEARLY/SAR';
    startYear = (year < 2007) ? 2007 : (year > 2010 ? 2010 : year);
  } else {
    sarCol = 'JAXA/ALOS/PALSAR/YEARLY/SAR_EPOCH';
    startYear = (year < 2015) ? 2015 : (year > 2024 ? 2024 : year);
  }

  var startDate = startYear + '-01-01';
  var endDate   = (startYear + 1) + '-01-01';
  var winStart  = (startYear - 1) + '-01-01';
  var winEnd    = (startYear + 2) + '-01-01';

  var calibratePalsar = function(img, useQa) {
    var pol = img.select(['HH', 'HV']).toFloat();
    var posMask = pol.select('HH').gt(0).and(pol.select('HV').gt(0));
    var lin = pol.pow(2).multiply(PALSAR_LIN_CF).updateMask(posMask);
    if (useQa) {
      var qa = img.select('qa');
      var qaValid = maskOcean ? qa.eq(255) : qa.eq(255).or(qa.eq(50));
      lin = lin.updateMask(qaValid);
    }
    return lin;
  };

  var targetRaw = ee.ImageCollection(sarCol).filterDate(startDate, endDate).first();
  var windowCol = ee.ImageCollection(sarCol).filterDate(winStart, winEnd);

  var targetLin = calibratePalsar(targetRaw, applyQaMask);
  var fallbackQaLin = windowCol.map(function(img) {
    return calibratePalsar(img, applyQaMask);
  }).median();

  var compLin = targetLin.unmask(fallbackQaLin);

  if (fillMaskedGaps && applyQaMask) {
    var fallbackRawLin = windowCol.map(function(img) {
      return calibratePalsar(img, false);
    }).median();
    compLin = compLin.unmask(fallbackRawLin);
  }

  if (applySpeckle) {
    compLin = exports.applySpeckleFilter(compLin, {
      filter: speckleFilter,
      kernelSize: kernelSize,
      enl: enl,
      bands: ['HH', 'HV']
    });
  }

  var dbImg = linToDb(compLin.resample('bilinear'), ['HH', 'HV']);
  var hh_db   = dbImg.select('HH').rename('palsar_HH');
  var hv_db   = dbImg.select('HV').rename('palsar_HV');
  var diff_db = hh_db.subtract(hv_db).rename('palsar_HH_HV_diff');

  var result = ee.Image.cat([hh_db, hv_db, diff_db]);
  return roi ? result.clip(roi) : result;
};

// ---------------------------------------------------------------------------
// 4. Topography & Terrain Geomorphology
// ---------------------------------------------------------------------------

/**
 * Builds time-invariant topographic and geomorphological feature stacks.
 * Computes elevation, slope, and aspect from NASADEM. Optionally attaches
 * landform (TPI) and Topographic Wetness Index (TWI) images if provided.
 *
 * @param {ee.Geometry} roi            - Region of interest
 * @param {Object}      [options]      - Optional external layers: { landformImg, twiImg }
 * @returns {ee.Image} Multi-band terrain image
 */
exports.getTopography = function(roi, options) {
  options = options || {};
  var elev = ee.Image('NASA/NASADEM_HGT/001')
    .select('elevation')
    .clip(roi);

  var terrain = ee.Terrain.products(elev).select(['slope', 'aspect']).toFloat();
  var stack = elev.int16().addBands(terrain);

  if (options.landformImg) {
    stack = stack.addBands(options.landformImg.rename('landform').toByte().clip(roi));
  }
  if (options.twiImg) {
    stack = stack.addBands(options.twiImg.rename('TWI').toFloat().clip(roi));
  }

  return stack;
};

// ---------------------------------------------------------------------------
// 5. CCDC Tile Loader
// ---------------------------------------------------------------------------

/**
 * Loads a list of CCDC tiled images and combines them into a seamless mosaic.
 *
 * @param {Array<string>} tileIds - Array of GEE asset IDs
 * @param {ee.Geometry}   roi     - Clip geometry
 * @returns {ee.Image} Mosaicked image clipped to ROI
 */
exports.loadCCDCMosaic = function(tileIds, roi) {
  var images = tileIds.map(function(id) { return ee.Image(id); });
  return ee.ImageCollection(images).mosaic().clip(roi);
};
