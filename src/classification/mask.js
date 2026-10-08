/**
 * Module: mask.js
 * Description: Builds named mask / evidence layers used by the rule-based
 *              probability adjustment (Rules.js). Every layer is an ee.Image,
 *              so rules reference them by name (e.g. m.permanentWater).
 *
 * Quick start:
 *   var mask = require('users/YOUR_USER/YOUR_REPO:mask');
 *   var m = mask.build(aoi, {mapYear: 2023, startDate: '2023-01-01', endDate: '2023-12-31'});
 *
 * Binary layers (0/1):
 *   permanentWater, seasonalWater, everWater          (JRC Global Surface Water)
 *   wcWater, wcBuilt, wcBare, wcCrop, wcTrees, wcMangrove   (ESA WorldCover v200)
 *   hrslNear, ghslBuilt      (only if options.hrsl / options.ghsl are provided)
 *   noBuiltEvidence          (1 where NO built-up layer votes for settlement)
 *   mangroveNear             (dilated mangrove extent)
 *   forestLoss, recentLoss   (Hansen GFC, relative to options.mapYear)
 *   oilPalm                  (only if options.oilPalm is provided)
 * Count layer:
 *   builtVotes               (number of built-up evidence layers that agree, 0..N)
 * Continuous layers:
 *   dwBuilt, dwWater, dwBare (Dynamic World mean probability, 0..1)
 *   elevation (m), slope (deg), treecover2000 (%)
 *
 * NOTE: focalMax() dilation is in meters and is resolved at the scale you
 * export/inspect at. Use the same scale as your probability map when exporting.
 */

var DEFAULTS = {
  mapYear: 2023,
  startDate: '2023-01-01',
  endDate: '2023-12-31',
  waterPermanent: 75,          // JRC occurrence (%) treated as permanent water
  builtDilateMeters: 100,      // dilation of built-up evidence layers
  dwBuiltMin: 0.3,             // Dynamic World mean 'built' probability threshold
  ghslMin: 500,                // GHSL built surface value treated as built
  mangroveDilateMeters: 500,   // buffer around mapped mangrove extent
  recentYears: 5,              // window for "recent" forest loss
  treeCoverMin: 30,            // Hansen treecover2000 (%) counted as forest
  hansenId: 'UMD/hansen/global_forest_change_2024_v1_12', // update to latest version
  // Optional user-supplied layers (ee.Image or ee.ImageCollection)
  hrsl: null,       // HRSL settlement (values > 0 = settlement)
  ghsl: null,       // GHSL built surface image
  mangrove: null,   // e.g. Global Mangrove Watch; default is LANDSAT/MANGROVE_FORESTS
  oilPalm: null,    // oil palm layer (values > 0 = oil palm)
  dem: null,        // elevation image; default MERIT DEM (bare-earth corrected)
  slopeDem: null    // DEM used for slope; default SRTM 30 m
};

function merge(options) {
  var o = {};
  var k;
  for (k in DEFAULTS) { o[k] = DEFAULTS[k]; }
  options = options || {};
  for (k in options) { o[k] = options[k]; }
  return o;
}

function toImage(x) {
  return (x instanceof ee.ImageCollection) ? x.mosaic() : ee.Image(x);
}

function dilate(img, meters) {
  return img.focalMax({radius: meters, units: 'meters'});
}

/**
 * Build all mask layers.
 * @param {ee.Geometry} aoi
 * @param {Object} options  see DEFAULTS
 * @return {Object} dictionary name -> ee.Image (plus builtVoteLayers: JS array of voter names)
 */
exports.build = function(aoi, options) {
  var o = merge(options);
  var m = {};

  // ---- Water: JRC Global Surface Water ----
  var gsw = ee.Image('JRC/GSW1_4/GlobalSurfaceWater');
  var occ = gsw.select('occurrence').unmask(0);
  m.permanentWater = occ.gte(o.waterPermanent).rename('permanentWater');
  m.seasonalWater  = occ.gt(0).and(occ.lt(o.waterPermanent)).rename('seasonalWater');
  m.everWater      = gsw.select('max_extent').unmask(0).eq(1).rename('everWater');

  // ---- ESA WorldCover v200 ----
  var wc = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map');
  m.wcTrees    = wc.eq(10).unmask(0).rename('wcTrees');
  m.wcCrop     = wc.eq(40).unmask(0).rename('wcCrop');
  m.wcBuilt    = wc.eq(50).unmask(0).rename('wcBuilt');
  m.wcBare     = wc.eq(60).unmask(0).rename('wcBare');
  m.wcWater    = wc.eq(80).unmask(0).rename('wcWater');
  m.wcMangrove = wc.eq(95).unmask(0).rename('wcMangrove');

  // ---- Dynamic World (mean probability over the map period) ----
  var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
    .filterBounds(aoi)
    .filterDate(o.startDate, o.endDate)
    .select(['built', 'water', 'bare'])
    .mean();
  m.dwBuilt = dw.select('built').unmask(0).rename('dwBuilt');
  m.dwWater = dw.select('water').unmask(0).rename('dwWater');
  m.dwBare  = dw.select('bare').unmask(0).rename('dwBare');

  // ---- Built-up evidence vote ----
  // Each voter says "settlement nearby" (dilated). noBuiltEvidence = no voter agrees.
  // Note: WorldCover and Dynamic World are not fully independent of each other.
  var voters = [
    dilate(m.wcBuilt, o.builtDilateMeters),
    dilate(m.dwBuilt.gte(o.dwBuiltMin), o.builtDilateMeters)
  ];
  var voterNames = ['wcBuilt', 'dwBuilt'];

  if (o.hrsl) {
    m.hrslNear = dilate(toImage(o.hrsl).gt(0).unmask(0), o.builtDilateMeters).rename('hrslNear');
    voters.push(m.hrslNear);
    voterNames.push('hrsl');
  }
  if (o.ghsl) {
    m.ghslBuilt = dilate(toImage(o.ghsl).gt(o.ghslMin).unmask(0), o.builtDilateMeters).rename('ghslBuilt');
    voters.push(m.ghslBuilt);
    voterNames.push('ghsl');
  }
  m.builtVotes = ee.Image.cat(voters).reduce(ee.Reducer.sum()).rename('builtVotes');
  m.noBuiltEvidence = m.builtVotes.eq(0).rename('noBuiltEvidence');
  m.builtVoteLayers = voterNames;   // plain JS array, for logging

  // ---- Mangrove extent ----
  var mg = o.mangrove
    ? toImage(o.mangrove).gt(0)
    : ee.ImageCollection('LANDSAT/MANGROVE_FORESTS').first().select(0).gt(0);
  m.mangroveNear = dilate(mg.unmask(0).or(m.wcMangrove), o.mangroveDilateMeters)
    .rename('mangroveNear');

  // ---- Forest loss (Hansen GFC) ----
  var hansen = ee.Image(o.hansenId);
  var lossyear = hansen.select('lossyear');
  var t = o.mapYear - 2000;
  var forest2000 = hansen.select('treecover2000').gte(o.treeCoverMin);
  m.forestLoss = forest2000.and(lossyear.gt(0)).and(lossyear.lte(t))
    .unmask(0).rename('forestLoss');
  m.recentLoss = forest2000.and(lossyear.gt(t - o.recentYears)).and(lossyear.lte(t))
    .unmask(0).rename('recentLoss');
  m.treecover2000 = hansen.select('treecover2000').rename('treecover2000');

  // ---- Oil palm (optional) ----
  if (o.oilPalm) {
    m.oilPalm = toImage(o.oilPalm).gt(0).unmask(0).rename('oilPalm');
  }

  // ---- Terrain ----
  var demImg = o.dem ? ee.Image(o.dem) : ee.Image('MERIT/DEM/v1_0_3').select('dem');
  var slopeSrc = o.slopeDem ? ee.Image(o.slopeDem) : ee.Image('USGS/SRTMGL1_003');
  m.elevation = demImg.rename('elevation');
  m.slope = ee.Terrain.slope(slopeSrc).rename('slope');

  return m;
};

/**
 * Trapezoidal plausibility weight (0..1) for a continuous image.
 * spec = [a, b, c, d]:  0 below a, ramps up to 1 at b, stays 1 until c, ramps down to 0 at d.
 * Use null pairs for an open end, e.g. [null, null, 10, 25] = "1 up to 10, fades to 0 at 25".
 */
exports.trapezoid = function(img, spec) {
  var a = spec[0], b = spec[1], c = spec[2], d = spec[3];
  var w = ee.Image.constant(1);
  if (a !== null && b !== null) {
    w = (b > a)
      ? w.min(img.subtract(a).divide(b - a).clamp(0, 1))
      : w.min(img.gte(a));
  }
  if (c !== null && d !== null) {
    w = (d > c)
      ? w.min(ee.Image.constant(d).subtract(img).divide(d - c).clamp(0, 1))
      : w.min(img.lte(c));
  }
  return w.rename('weight');
};

/**
 * Terrain plausibility weight combining elevation and slope (fuzzy AND = min).
 * spec = {elev: [a,b,c,d], slope: [c,d]}  (either key optional)
 */
exports.terrainWeight = function(m, spec) {
  var w = ee.Image.constant(1);
  if (spec.elev)  { w = w.min(exports.trapezoid(m.elevation, spec.elev)); }
  if (spec.slope) { w = w.min(exports.trapezoid(m.slope, [null, null, spec.slope[0], spec.slope[1]])); }
  return w.rename('weight');
};

/**
 * Stack selected layers into one multiband float image (e.g. for export to asset).
 */
exports.stack = function(m, keys) {
  return ee.Image.cat(keys.map(function(k) {
    return m[k].rename(k).toFloat();
  }));
};

/**
 * Mean of each layer over the AOI. For binary layers this is the area fraction.
 */
exports.summarize = function(m, keys, aoi, scale) {
  return exports.stack(m, keys).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: scale,
    maxPixels: 1e10,
    bestEffort: true
  });
};

/**
 * Add binary layers to the map (hidden by default) for visual checking.
 */
exports.addToMap = function(m, keys, mapObj) {
  mapObj = mapObj || Map;
  keys.forEach(function(k) {
    mapObj.addLayer(m[k].selfMask(), {palette: ['#e31a1c']}, 'mask: ' + k, false);
  });
};
