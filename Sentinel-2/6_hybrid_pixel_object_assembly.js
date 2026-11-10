///////// 1. Setup /////////
var cfg = require('users/rg2icraf/Luma:config');
var post = require('users/rg2icraf/Luma:post_process');
var testAOI = geometry;
var MAX_SIZE = 1024;        // Larger limits create oversized tile neighborhoods.
var USE_GUARD = false;      // true = extra within-object homogeneity guard
var GUARD_SD = 0.25;        // guard threshold
var W_DEFAULT = 0.5;        // default object weight
var STABILITY_TAU = 0.65;   // Threshold from Malinverni & Sanità (2026): CI <= tau is stable
var ADAPTIVE_W = true;      // true = dynamically scale object weight by stability
var ADAPTIVE_SENSITIVITY = 0.6; // How strongly instability damps object weight (0.0 to 1.0)
var APPLY_RULE_REFINEMENT = true; // true = apply intrinsic rule-based refinement for common confusion
var MERGE_WATER_RIVERBED = true;  // true = merge Dry riverbed (18) into Waterbody (17) after classification
var classIdToName = cfg.CLASS_SCHEMA;
// Bayesian smoothing / object trust weights per class — sourced from cfg
var smoothnessTL = cfg.SMOOTHNESS;
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
//print('Class bands:', classBands);
//define the segments
var pieces = [
  ee.Image('projects/gee-v1-510606/assets/snic/seg_central_timor'),
  ee.Image('projects/gee-v1-510606/assets/snic/seg_east_box'),
  ee.Image('projects/gee-v1-510606/assets/snic/seg_off_dili'),
  ee.Image('projects/gee-v1-510606/assets/snic/seg_west_box'),
  ee.Image('projects/gee-v1-510606/assets/snic/seg_east_timor')
];
/*
// Check the ID range of each piece (approximate, sampled at 100 m).
pieces.forEach(function(s, i) {
  print('piece ' + (i + 1) + ' ID range', s.select([0]).reduceRegion({
    reducer: ee.Reducer.minMax(), geometry: s.geometry(), scale: 100,
    maxPixels: 1e10, bestEffort: true, tileScale: 4
  }));
});
*/
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

// ---------- 5. FUNCTIONS & HYBRID ASSEMBLY (via Post-process module) ----------
function integrate(probs) {
  return post.adaptiveBlend(probs, snic, {
    classBands: classBandNames,
    classValues: classIdsList,
    weights: wClientDict,
    maxSize: MAX_SIZE,
    adaptive: ADAPTIVE_W,
    adaptiveSensitivity: ADAPTIVE_SENSITIVITY,
    stabilityTau: STABILITY_TAU,
    useGuard: USE_GUARD,
    guardSd: GUARD_SD
  });
}

function refineRules(hybResult) {
  var refinedLulc = post.suppressFalseHydrology(
    hybResult.classification,
    hybResult.runner_up,
    hybResult.ci,
    {
      ciThreshold: 65,
      riverbedCiThreshold: 75,
      waterClasses: [17, 18],
      terrestrialClasses: [15, 16],
      riverbedId: 18,
      waterbodyId: 17,
      clearedLandId: 15,
      mergeRiverbed: MERGE_WATER_RIVERBED
    }
  );
  return hybResult.combined.addBands(refinedLulc.rename('lulc_refined'), null, true);
}

// ---------- 6. TEST ON A SMALL AREA ----------
var pix2025_raw = post.classifyFromProbs(lc_prob_2025, classBandNames, classIdsList);
var pix2025 = MERGE_WATER_RIVERBED ? pix2025_raw.where(pix2025_raw.eq(18), 17) : pix2025_raw;

var hyb2025_res = integrate(lc_prob_2025);
var hyb2025 = refineRules(hyb2025_res);

var holes = snic.mask().and(snic.reduceConnectedComponents({
  reducer: ee.Reducer.mean(), labelBand: 'segment', maxSize: MAX_SIZE
}).mask().not()).selfMask();

var viz = {min: 1, max: 18, palette: cfg.CLASS_PALETTE};
Map.centerObject(testAOI, 12);
Map.addLayer(pix2025.clip(testAOI), viz, 'Pixel-only 2025', false);
Map.addLayer(hyb2025.select('lulc').clip(testAOI), viz, 'Hybrid 2025 (Raw LULC)', false);
Map.addLayer(hyb2025.select('lulc_refined').clip(testAOI), viz, 'Hybrid 2025 (Refined LULC)');
Map.addLayer(hyb2025.select('is_stable').clip(testAOI), {min: 0, max: 1, palette: ['#D7191C', '#1A9641']}, 'Stability Map (Green=Stable, Red=Unstable)', false);
Map.addLayer(hyb2025.select('ci').clip(testAOI), {min: 0, max: 100, palette: ['#1A9641', '#FFFFBF', '#D7191C']}, 'Confusion Index % (S/W)', false);
Map.addLayer(hyb2025.select('conf').clip(testAOI), {min: 0, max: 100}, 'Confidence % (Winner)', false);
Map.addLayer(hyb2025.select('runner_up').clip(testAOI), viz, 'Runner-up Class', false);
Map.addLayer(holes.clip(testAOI), {palette: ['red']}, 'Holes (MAX_SIZE exceeded)', false);

print('Share of test area covered by segments:', snic.mask().reduceRegion({
  reducer: ee.Reducer.mean(), geometry: testAOI, scale: 10, maxPixels: 1e9, tileScale: 4
}));

print('Share of test area classified as Stable (CI <= ' + STABILITY_TAU + '):', hyb2025.select('is_stable').reduceRegion({
  reducer: ee.Reducer.mean(), geometry: testAOI, scale: 10, maxPixels: 1e9, tileScale: 4
}));

// Evaluate impact of rule refinement
var modifiedByRules = hyb2025.select('lulc').neq(hyb2025.select('lulc_refined'));
print('Share of test area adjusted by rule refinement:', modifiedByRules.reduceRegion({
  reducer: ee.Reducer.mean(), geometry: testAOI, scale: 10, maxPixels: 1e9, tileScale: 4
}));