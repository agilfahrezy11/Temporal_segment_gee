/**
 * Module: Feature Engineering & Multi-Sensor Stacking
 * =========================================================================
 * Single source of truth for all feature extraction and processing helpers:
 *   - computeSeasonalIndices: Reconstructs spectral indices from synthetic CCDC stacks.
 *   - computePhenologicalDeltas: Computes seasonal index differences (wet - dry).
 *   - getSeasonalS1: Multi-temporal Sentinel-1 GRD SAR medians (VV/VH, wet/dry).
 *   - getScanSARFeatures: Annual PALSAR-2 ScanSAR median backscatter in dB.
 *   - getPalsarBackscatter: Annual ALOS PALSAR / PALSAR-2 L-band backscatter in dB.
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
 * @param {string}        sensor     - 'S2' (maps nir to B8, swir1 to B11, swir2 to B12) or
 *                                     'Landsat' (maps nir to B5, swir1 to B6, swir2 to B7)
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
// 3. Radar / SAR Backscatter Composites
// ---------------------------------------------------------------------------

/**
 * Generates seasonal median Sentinel-1 GRD SAR backscatter composites (VV & VH).
 * Default monsoonal season split:
 *   - Wet Season: December through April
 *   - Dry Season: May through November
 *
 * @param {number}      year - Target year
 * @param {ee.Geometry} roi  - Region of interest
 * @param {Object}      [options] - { wetMonths: [12, 4], dryMonths: [5, 11] }
 * @returns {ee.Image} 4-band image: ['VH_dry', 'VH_wet', 'VV_dry', 'VV_wet']
 */
exports.getSeasonalS1 = function(year, roi, options) {
  options = options || {};
  var wetMonths = options.wetMonths || [12, 4];
  var dryMonths = options.dryMonths || [5, 11];

  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = ee.Date.fromYMD(year, 12, 31);

  var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .select(['VV', 'VH']);

  var wet = s1.filter(ee.Filter.calendarRange(wetMonths[0], wetMonths[1], 'month')).median();
  var dry = s1.filter(ee.Filter.calendarRange(dryMonths[0], dryMonths[1], 'month')).median();

  var vv_wet = wet.select('VV').rename('VV_wet');
  var vh_wet = wet.select('VH').rename('VH_wet');
  var vv_dry = dry.select('VV').rename('VV_dry');
  var vh_dry = dry.select('VH').rename('VH_dry');

  return ee.Image.cat([vh_dry, vh_wet, vv_dry, vv_wet]).clip(roi);
};

/**
 * Annual ALOS-2 PALSAR-2 ScanSAR dual-polarization backscatter in dB.
 * Converts positive digital numbers to gamma-nought backscatter in dB:
 *   gamma0_dB = 10 * log10(DN^2) - 83.0
 *
 * @param {number}      year - Target year
 * @param {ee.Geometry} roi  - Region of interest
 * @returns {ee.Image} 2-band image: ['ScanSAR_HH_dB', 'ScanSAR_HV_dB']
 */
exports.getScanSARFeatures = function(year, roi) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = startDate.advance(1, 'year');

  var collection = ee.ImageCollection('JAXA/ALOS/PALSAR-2/Level2_2/ScanSAR')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.listContains('Polarizations', 'HH'))
    .filter(ee.Filter.listContains('Polarizations', 'HV'));

  var dnMedian = collection.select(['HH', 'HV']).median();
  return dnMedian
    .updateMask(dnMedian.gt(0))
    .pow(2)
    .log10()
    .multiply(10)
    .subtract(83)
    .rename(['ScanSAR_HH_dB', 'ScanSAR_HV_dB'])
    .clip(roi);
};

/**
 * Retrieves ALOS PALSAR / PALSAR-2 L-band SAR yearly mosaic and converts to dB.
 * Uses a 3-year median window as fallback for orbital gaps.
 * Formula: gamma0_dB = 10 * log10(DN^2) - 83.0
 *
 * @param {number}      year - Target year (proxy mapping: <=2010 uses PALSAR, >2010 uses PALSAR-2)
 * @param {ee.Geometry} [roi] - Optional clip geometry
 * @returns {ee.Image} 3-band image: ['palsar_HH', 'palsar_HV', 'palsar_HH_HV_diff']
 */
exports.getPalsarBackscatter = function(year, roi) {
  var sarCol, startYear;
  if (year <= 2010) {
    sarCol = 'JAXA/ALOS/PALSAR/YEARLY/SAR';
    startYear = (year < 2007) ? 2007 : year;
  } else {
    sarCol = 'JAXA/ALOS/PALSAR/YEARLY/SAR_EPOCH';
    startYear = (year > 2024) ? 2024 : year;
  }

  var startDate = startYear + '-01-01';
  var endDate   = (startYear + 1) + '-01-01';

  var sarImg = ee.ImageCollection(sarCol)
    .filterDate(startDate, endDate)
    .first();

  var fallback = ee.ImageCollection(sarCol)
    .filterDate((startYear - 1) + '-01-01', (startYear + 2) + '-01-01')
    .median();

  var hh = sarImg.select('HH')
    .unmask(fallback.select('HH'))
    .toFloat()
    .resample('bilinear')
    .max(1.0);

  var hv = sarImg.select('HV')
    .unmask(fallback.select('HV'))
    .toFloat()
    .resample('bilinear')
    .max(1.0);

  var hh_db   = hh.pow(2).log10().multiply(10.0).subtract(83.0).rename('palsar_HH');
  var hv_db   = hv.pow(2).log10().multiply(10.0).subtract(83.0).rename('palsar_HV');
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
