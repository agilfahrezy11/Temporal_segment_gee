///////// 1. Setup /////////
var testAOI = geometry;
var MAX_SIZE = 1024;        // Larger limits create oversized tile neighborhoods.
var USE_GUARD = false;      // true = extra within-object homogeneity guard
var GUARD_SD = 0.25;        // guard threshold
var W_DEFAULT = 0.6;        // default object weight
var classIdToName = {
  1:  'Primary_forest',
  2:  'Secondary_forest',
  3:  'Mangrove_forest',
  4:  'Coastal_forest',
  5:  'Teak_plantation',
  6:  'Eucalyptus_plantation',
  7:  'Coffe_agroforestry',
  8:  'Cocoa_agroforestry',
  9:  'Coconut_agroforestry',
  10: 'Mixed_garden',
  11: 'Paddy_field',
  12: 'Other_cropland',
  13: 'Grassland',
  14: 'Shrubland',
  15: 'Cleared_land',
  16: 'Built-up',
  17: 'Waterbody',
  18: 'Dry_riverbed'
};
// Bayesian smoothing / object trust weights per class
var smoothnessTL = {
  'Primary_forest':        4.5,
  'Secondary_forest':      3.0,
  'Mangrove_forest':       3.0,
  'Coastal_forest':        3.5,
  'Teak_plantation':       1.2,
  'Eucalyptus_plantation': 1.2,
  'Coffe_agroforestry':    1.5,
  'Cocoa_agroforestry':    1.5,
  'Coconut_agroforestry':  1.5,
  'Mixed_garden':          1.2,
  'Paddy_field':           1.0,
  'Other_cropland':        1.0,
  'Grassland':             1.5,
  'Shrubland':             1.5,
  'Cleared_land':          1.5,
  'Built-up':              1.0,
  'Waterbody':             4.0,
  'Dry_riverbed':          0.1
};
var lc_prob_2025 = ee.Image('projects/gee-v1-510606/assets/TL_S2_class_probabilities_2025').divide(10000);
var classBands = lc_prob_2025.bandNames();
// Build the class IDs and weights from the probability band's names.
var nameToId = {};
for (var id in classIdToName) {
  nameToId[classIdToName[id]] = parseInt(id, 10);
}
var classBandNames = classBands.getInfo();
var classIdsList = [];
var wClientDict = {};
classBandNames.forEach(function(bName) {
  if (nameToId[bName] === undefined) {
    throw new Error('No class ID configured for probability band: ' + bName);
  }
  classIdsList.push(nameToId[bName]);

  var smoothness = smoothnessTL[bName] || 1.0;
  wClientDict[bName] = Math.min(Math.max(smoothness / 5.0, 0.2), 0.9);
});
var classIds = ee.List(classIdsList);
print('Class bands:', classBands);
//define the segments
var pieces = [
  ee.Image('projects/gee-v1-510606/assets/snic/seg_central_timor'),
  ee.Image('projects/gee-v1-510606/assets/snic/seg_east_box'),
  ee.Image('projects/gee-v1-510606/assets/snic/seg_off_dili'),
  ee.Image('projects/gee-v1-510606/assets/snic/seg_west_box'),
  ee.Image('projects/gee-v1-510606/assets/snic/seg_east_timor')
];
// Check the ID range of each piece (approximate, sampled at 100 m).
pieces.forEach(function(s, i) {
  print('piece ' + (i + 1) + ' ID range', s.select([0]).reduceRegion({
    reducer: ee.Reducer.minMax(), geometry: s.geometry(), scale: 100,
    maxPixels: 1e10, bestEffort: true, tileScale: 4
  }));
});
// Give each piece its own ID range. ID_STRIDE must be larger than the biggest ID inside any
// piece, and 4 * ID_STRIDE + biggest ID must stay under 2,147,483,647 (int32 limit).
var ID_STRIDE = 300000000;
var snic = ee.ImageCollection(pieces.map(function(img, i) {
  return img.select([0]).add(i * ID_STRIDE).toInt32();
})).mosaic().rename('segment');
// ---------- 4. WEIGHTS AND CLASS IDS ----------
var wDict = ee.Dictionary(wClientDict);
var w = ee.Image.constant(classBands.map(function(b) {
  return wDict.get(b, W_DEFAULT);
})).rename(classBands);

// ---------- 5. FUNCTIONS ----------
function aggregate(img, reducer) {
  return img.addBands(snic).reduceConnectedComponents({
    reducer: reducer, labelBand: 'segment', maxSize: MAX_SIZE
  }).select(classBands);
}

function toClass(p) {
  var idx = p.toArray().arrayArgmax().arrayGet(0);
  return idx.remap(ee.List.sequence(0, classBands.size().subtract(1)), classIds)
    .rename('lulc').toUint8();
}

function integrate(probs) {
  var obj = aggregate(probs, ee.Reducer.mean()).unmask(probs);

  var wEff = w;
  if (USE_GUARD) {
    var sd = aggregate(probs, ee.Reducer.stdDev()).reduce(ee.Reducer.max());
    wEff = w.multiply(sd.unmask(0).lte(GUARD_SD));
  }

  var blend = obj.multiply(wEff).add(probs.multiply(ee.Image(1).subtract(wEff)));
  var conf = blend.reduce(ee.Reducer.max()).multiply(100).round().toUint8().rename('conf');
  return toClass(blend).addBands(conf);
}
// ---------- 6. TEST ON A SMALL AREA ----------
var pix2025 = toClass(lc_prob_2025);
var hyb2025 = integrate(lc_prob_2025);

var holes = snic.mask().and(aggregate(lc_prob_2025, ee.Reducer.mean()).select(0).mask().not()).selfMask();
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
  '#DB310D', // 16: Built-up
  '#0e8acc', // 17: Waterbody
  '#9ECAE1'  // 18: Dry_riverbed
];
var viz = {min: 1, max: 18, palette: classPalette};
Map.centerObject(testAOI, 12);
Map.addLayer(pix2025.clip(testAOI), viz, 'Pixel-only 2025');
Map.addLayer(hyb2025.select('lulc').clip(testAOI), viz, 'Hybrid 2025');
Map.addLayer(hyb2025.select('conf').clip(testAOI), {min: 0, max: 100}, 'Confidence', false);
Map.addLayer(holes.clip(testAOI), {palette: ['red']}, 'Holes (MAX_SIZE exceeded)', false);

print('Share of test area covered by segments:', snic.mask().reduceRegion({
  reducer: ee.Reducer.mean(), geometry: testAOI, scale: 10, maxPixels: 1e9, tileScale: 4
}));