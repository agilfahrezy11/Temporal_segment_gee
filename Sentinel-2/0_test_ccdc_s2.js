var ccdc = require("users/rg2icraf/Luma:ccdc_sentinel2");

var test = ee.Image('projects/ee-rg2v2/assets/CCDC_TL_2020_2025_20m_test');
var testArea = test.geometry();
Map.centerObject(testArea, 11);

// Synthetic stacks at the two target dates
var s2020 = ccdc.getSyntheticStack(test, 2020.5);
var s2025 = ccdc.getSyntheticStack(test, 2025.5);

// 1. Visual check
var rgbVis = {bands: ['swir1', 'nir', 'red'], min: 0.0, max: 0.4};
var ndviVis = {bands: ['NDVI'], min: -0.2, max: 0.9,
               palette: ['#8c510a', '#f6e8c3', '#c7eae5', '#5ab4ac', '#01665e']};

Map.addLayer(s2020, rgbVis, 'SWIR1/NIR/Red 2020.5');
Map.addLayer(s2025, rgbVis, 'SWIR1/NIR/Red 2025.5');
Map.addLayer(s2020, ndviVis, 'NDVI 2020.5');
Map.addLayer(s2025, ndviVis, 'NDVI 2025.5');

function validFraction(img, label) {
  var frac = img.select('NDVI').mask().reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: testArea,
    scale: 20,
    maxPixels: 1e10,
    bestEffort: true
  });
  print('Valid fraction ' + label, frac);
}
validFraction(s2020, '2020.5');
validFraction(s2025, '2025.5');

// 3. Value ranges (look for outliers)
function stats(img, label) {
  var s = img.select(['NDVI', 'NDWI', 'NDBI', 'nir', 'swir1']).reduceRegion({
    reducer: ee.Reducer.percentile([1, 50, 99]),
    geometry: testArea,
    scale: 100,          // coarse is enough for a sanity check
    maxPixels: 1e10,
    bestEffort: true
  });
  print('Percentiles ' + label, s);
}
stats(s2020, '2020.5');
stats(s2025, '2025.5');

// 4. Segments per pixel (how often CCDC found a break)
var nSeg = test.select('tStart').arrayLength(0);
Map.addLayer(nSeg, {min: 1, max: 4, palette: ['white', 'yellow', 'orange', 'red']}, 'Segments per pixel');