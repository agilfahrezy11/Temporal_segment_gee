/**
 * SCRIPT: 8_test_topographic_settlement_constraints.js
 * PURPOSE: Interactive testing of Topographic Constraints (MERIT Hydro HAND + UPA)
 *          and High Resolution Settlement Layer (HRSL 30m) to resolve Built-up (16)
 *          vs. Dry Riverbed (18) / Waterbody (17) confusion.
 * 
 * FIX:
 *   Replaced all mapped client-side calls with pure Earth Engine server-side
 *   constructs to eliminate the "A mapped function's arguments cannot be used
 *   in client-side operations" error.
 */

// =========================================================================
// 1. CONFIGURATION & PRESET TEST REGIONS
// =========================================================================
var PRESETS = {
  'Comoro_Dili':     { lon: 125.534, lat: -8.560, zoom: 14, desc: 'Comoro River (Dili Airport/Town)' },
  'Laclo_Manatuto':  { lon: 126.012, lat: -8.514, zoom: 13, desc: 'Laclo Braided Riverbed (Manatuto)' },
  'Lois_Liquica':    { lon: 125.174, lat: -8.721, zoom: 13, desc: 'Lois Riverbed (Liquiçá/Bobonaro border)' },
  'Caraulun_Same':   { lon: 125.645, lat: -9.052, zoom: 13, desc: 'Caraulun River (South Coast)' }
};

// Choose test location: 'Comoro_Dili', 'Laclo_Manatuto', 'Lois_Liquica', 'Caraulun_Same'
var ACTIVE_PRESET = 'Comoro_Dili';

var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry();
var CRS = 'EPSG:32751';

// Topographic & Settlement thresholds
var HAND_MAX_RIVER   = 2.5;   // Pixels <= 2.5m above nearest drainage line
var UPA_MIN_RIVER    = 5.0;   // Upstream drainage area >= 5 km² (filters minor gullies)
var HAND_MIN_URBAN   = 5.0;   // Pixels > 5.0m above drainage are safe from flood channel
var MERGE_WATER_RIVERBED = true; // Final merge of 18 into 17

// Center Map on selected preset
if (PRESETS[ACTIVE_PRESET]) {
  var p = PRESETS[ACTIVE_PRESET];
  Map.setCenter(p.lon, p.lat, p.zoom);
  print('Active Test Region:', p.desc);
}

// =========================================================================
// 2. LOAD DATASETS
// =========================================================================

// A. 2025 S2 Probability Asset (scaled 0 to 1)
var lc_prob_2025 = ee.Image('projects/gee-v1-510606/assets/TL_S2_class_probabilities_2025').divide(10000);

// B. SNIC Segments
var pieces = [
  ee.Image('projects/gee-v1-510606/assets/snic/seg_central_timor'),
  ee.Image('projects/gee-v1-510606/assets/snic/seg_east_box'),
  ee.Image('projects/gee-v1-510606/assets/snic/seg_off_dili'),
  ee.Image('projects/gee-v1-510606/assets/snic/seg_west_box'),
  ee.Image('projects/gee-v1-510606/assets/snic/seg_east_timor')
];
var ID_STRIDE = 300000000;
var snic = ee.ImageCollection(pieces.map(function(img, i) {
  return img.select([0]).add(i * ID_STRIDE).toInt32();
})).mosaic().rename('segment');

// C. High Resolution Settlement Layer (HRSL 30m - Meta / CIESIN)
var hrslColl = ee.ImageCollection('projects/sat-io/open-datasets/hrsl/hrslpop');
var hrsl = hrslColl.filterBounds(aoi).mosaic();
var settlementMask = hrsl.gt(0).unmask(0).rename('settlement');

// D. Topography & Hydrology
// 1. MERIT Hydro v1.0.1 (HAND & Upstream Catchment)
var meritHydro = ee.Image('MERIT/Hydro/v1_0_1');
var meritHand = meritHydro.select('hnd').rename('hand'); // Height Above Nearest Drainage (m)
var meritUpa  = meritHydro.select('upa').rename('upa');  // Upstream drainage area (km²)

// 2. JRC Global Surface Water (Historical flood occurrence %)
var jrcWater = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0);

// 3. Sentinel-2 2025 Optical Composite for False Color reference
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate('2025-01-01', '2026-01-01')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .median();

// =========================================================================
// 3. DEFINE ACTIVE RIVERBED CORRIDOR MASK
// =========================================================================
// A pixel is in an active drainage corridor if:
// 1. HAND <= 2.5m (in the valley bottom)
// 2. AND (Upstream catchment >= 5 km² OR documented water occurrence > 0%)
var isRiverCorridor = meritHand.lte(HAND_MAX_RIVER)
  .and(meritUpa.gte(UPA_MIN_RIVER).or(jrcWater.gt(0)))
  .rename('river_corridor');

// =========================================================================
// 4. HYBRID CLASSIFICATION (CLIENT-SAFE SERVER-SIDE ASSEMBLY)
// =========================================================================
var classBandNames = [
  'Primary_forest',        'Secondary_forest',    'Mangrove_forest',
  'Coastal_forest',        'Teak_plantation',     'Eucalyptus_plantation',
  'Coffe_agroforestry',    'Cocoa_agroforestry',  'Coconut_agroforestry',
  'Mixed_garden',          'Paddy_field',         'Other_cropland',
  'Grassland',             'Shrubland',           'Cleared_land',
  'Built-up',              'Waterbody',           'Dry_riverbed'
];

var classIdsList = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18];

var smoothnessTL = {
  'Primary_forest': 4.5, 'Secondary_forest': 3.0, 'Mangrove_forest': 3.0, 'Coastal_forest': 3.5,
  'Teak_plantation': 1.2, 'Eucalyptus_plantation': 1.2, 'Coffe_agroforestry': 1.5, 'Cocoa_agroforestry': 1.5,
  'Coconut_agroforestry': 1.5, 'Mixed_garden': 1.2, 'Paddy_field': 1.0, 'Other_cropland': 1.0,
  'Grassland': 1.5, 'Shrubland': 1.5, 'Cleared_land': 1.5, 'Built-up': 1.0, 'Waterbody': 4.0, 'Dry_riverbed': 0.1
};

var wList = classBandNames.map(function(name) {
  var s = smoothnessTL[name] || 1.0;
  return Math.min(Math.max(s / 5.0, 0.2), 0.9);
});
var w = ee.Image.constant(wList).rename(classBandNames);

// Segment mean aggregation (10m)
var objMeans = lc_prob_2025.addBands(snic).reduceConnectedComponents({
  reducer: ee.Reducer.mean(),
  labelBand: 'segment',
  maxSize: 1024
}).select(classBandNames);

// Adaptive convex blending (Malinverni & Sanità 2026)
var diff = objMeans.subtract(lc_prob_2025).abs();
var wEff = w.multiply(ee.Image(1.0).subtract(diff.multiply(0.6))).clamp(0.05, 0.95);
var blended = objMeans.multiply(wEff).add(lc_prob_2025.multiply(ee.Image(1.0).subtract(wEff)));

// Winning class (W)
var W = blended.reduce(ee.Reducer.max()).rename('W');
var top1Idx = blended.toArray().arrayArgmax().arrayGet([0]);
var winner = top1Idx.remap(
  ee.List.sequence(0, 17),
  classIdsList
).rename('winner').toUint8();

// Runner-up class (S)
var maskedBands = classBandNames.map(function(bName, idx) {
  var bImg = blended.select([bName]);
  var bId = classIdsList[idx];
  return bImg.where(winner.eq(bId), -1.0).rename(bName);
});
var maskedStack = ee.Image.cat(maskedBands);
var S = maskedStack.reduce(ee.Reducer.max()).max(0.0).rename('S');
var top2Idx = maskedStack.toArray().arrayArgmax().arrayGet([0]);
var runner = top2Idx.remap(
  ee.List.sequence(0, 17),
  classIdsList
).rename('runner').toUint8();

var ci = S.divide(W.max(1e-4)).min(1.0).multiply(100).rename('ci');

// =========================================================================
// 5. BASELINE VS. REFINED WITH TOPOGRAPHY & HRSL
// =========================================================================

// Baseline classification (before new constraints)
var lulcBaseline = winner;
if (MERGE_WATER_RIVERBED) {
  lulcBaseline = lulcBaseline.where(lulcBaseline.eq(18), 17);
}

// -------------------------------------------------------------------------
// NEW CONSTRAINT LOGIC: Topography (HAND + UPA) + Settlement (HRSL)
// -------------------------------------------------------------------------

// Rule 1: False Built-up in River Corridor -> Dry Riverbed (18)
// Condition:
// - Winner is Built-up (16)
// - Pixel is inside active river corridor (HAND <= 2.5m, UPA >= 5km²)
// - NO confirmed human settlement in HRSL (settlementMask == 0)
var falseBuiltupInRiver = winner.eq(16)
  .and(isRiverCorridor)
  .and(settlementMask.eq(0));

var lulcRefined = winner.where(falseBuiltupInRiver, 18);

// Rule 2: False Riverbed in Real Towns -> Built-up (16)
// Condition:
// - Winner is Riverbed (18) or Waterbody (17)
// - HRSL confirms human settlement presence (settlementMask == 1)
// - Terrain is safely elevated above flood corridor (HAND > 5.0m)
var falseRiverbedInTown = lulcRefined.eq(17).or(lulcRefined.eq(18))
  .and(settlementMask.eq(1))
  .and(meritHand.gt(HAND_MIN_URBAN));

lulcRefined = lulcRefined.where(falseRiverbedInTown, 16);

// Rule 3: Merge Dry Riverbed (18) into Waterbody (17)
if (MERGE_WATER_RIVERBED) {
  lulcRefined = lulcRefined.where(lulcRefined.eq(18), 17);
}

// Difference Layer: Pixels corrected from Built-up to Riverbed/Water
var correctedBuiltupToWater = winner.eq(16).and(lulcRefined.eq(17));

// =========================================================================
// 6. VISUALIZATION & LAYERS
// =========================================================================
var classPalette = [
  '#00441B', // 1: Primary_forest
  '#238B45', // 2: Secondary_forest
  '#006D2C', // 3: Mangrove_forest
  '#66C2A4', // 4: Coastal_forest
  '#8C6D31', // 5: Teak_plantation
  '#A6761D', // 6: Eucalyptus_plantation
  '#7FC97F', // 7: Coffe_agroforestry
  '#66A061', // 8: Cocoa_agroforestry
  '#41AB5D', // 9: Coconut_agroforestry
  '#A1D76A', // 10: Mixed_garden
  '#C7E9B4', // 11: Paddy_field
  '#FDD835', // 12: Other_cropland
  '#D9EF8B', // 13: Grassland
  '#C2A878', // 14: Shrubland
  '#d7d6d6', // 15: Cleared_land
  '#DB310D', // 16: Built-up (RED)
  '#0e8acc', // 17: Waterbody (BLUE)
  '#9ECAE1'  // 18: Dry_riverbed
];
var vizLULC = {min: 1, max: 18, palette: classPalette};

// Add Layers to Map
Map.addLayer(s2, {bands: ['B8', 'B4', 'B3'], min: 500, max: 3500}, 'S2 False Color (NIR/Red/Green)', false);
Map.addLayer(s2, {bands: ['B4', 'B3', 'B2'], min: 200, max: 2000}, 'S2 True Color (RGB)', false);

// Constraint reference layers
Map.addLayer(meritHand, {min: 0, max: 20, palette: ['#08306b', '#2171b5', '#6baed6', '#c6dbef', '#ffffff']}, 'MERIT HAND (Dark Blue = Valley Floor <= 2m)', false);
Map.addLayer(isRiverCorridor.selfMask(), {palette: ['#00ffff']}, 'Active River Corridor Mask (Cyan)', false);
Map.addLayer(settlementMask.selfMask(), {palette: ['#ffff00']}, 'HRSL 30m Settlement Footprints (Yellow)', false);

// Classification layers
Map.addLayer(lulcBaseline, vizLULC, '1. BASELINE LULC (Has False Built-up in River)', false);
Map.addLayer(lulcRefined, vizLULC, '2. REFINED LULC (Corrected with Topo + HRSL)');
Map.addLayer(correctedBuiltupToWater.selfMask(), {palette: ['#ff00ff']}, '3. CORRECTION DIFF: Built-up -> Riverbed (Magenta)');

// =========================================================================
// 7. DIAGNOSTIC CONSOLE OUTPUT
// =========================================================================
print('--- Topographic & Settlement Refinement Diagnostics ---');
print('Active Preset:', ACTIVE_PRESET);
print('HAND Riverbed Threshold:', HAND_MAX_RIVER, 'm');
print('Upstream Catchment Threshold:', UPA_MIN_RIVER, 'km²');
print('Layers loaded: Turn on "3. CORRECTION DIFF (Magenta)" to inspect corrected river segments!');
