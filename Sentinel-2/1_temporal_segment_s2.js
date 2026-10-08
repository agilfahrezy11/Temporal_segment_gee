//define the AOI
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry()
Map.addLayer(aoi)
//////TESTING CCDC WORKFLOW////

//import the ccdc modules
var ccdc = require("users/rg2icraf/Luma:ccdc_sentinel2");
/*
//test point near DIli
var testPoint = ee.Geometry.Point([125.5736, -8.5568]); 

// Run CCDC on a small test geometry for fast feedback
var testRoi = testPoint.buffer(5000); 
var s2TS = ccdc.getS2TimeSeries(testRoi, '2018-12-01', '2025-12-31', 0.60, 'cs_cdf');
var ccdcTest = ccdc.runCCDC(s2TS);

// Generate synthetic stacks for 2019 and 2024 (Fraction 0.5 = July)
var synth2019 = ccdc.getSyntheticStack(ccdcTest, 2019.5);
var synth2024 = ccdc.getSyntheticStack(ccdcTest, 2024.5);

// Extract original surface reflectance and predicted NIR trajectory
var nirObs = s2TS.select('nir').map(function(img) {
  return img.select('nir').divide(10000)
    .copyProperties(img, ['system:time_start']);
});

// Plot time-series trajectory at test point
var chart = ui.Chart.image.series(nirObs, testPoint, ee.Reducer.mean(), 10)
  .setOptions({
    title: 'S2 NIR Raw Observations over Time',
    hAxis: {title: 'Date'},
    vAxis: {title: 'Surface Reflectance'},
    lineWidth: 1,
    pointSize: 3
  });
print(chart);

// Add to Map in True-Color (RGB)
var rgbVis = {bands: ['red', 'green', 'blue'], min: 0.02, max: 0.20};
Map.centerObject(testPoint, 13);
Map.addLayer(synth2019, rgbVis, 'Synthetic RGB 2019.5');
Map.addLayer(synth2024, rgbVis, 'Synthetic RGB 2024.5');
*/

//////RUN CCDC WORKFLOW////

//Prep collection
var raw_ts = ccdc.getS2TimeSeries(aoi, '2019-01-01', '2026-1-1', 0.50, 'cs_cdf');
var median_ts = ccdc.toMonthlyMedian(raw_ts, '2019-01-01', '2026-1-1');

//run CCDC algorithm
var ccdcResults = ccdc.runCCDC(median_ts, {
  breakpointBands: ['green', 'nir', 'swir1']
});
// drop what you don't need 
var selected_ccdc = ccdcResults.select(['tStart','tEnd','tBreak','changeProb','.*_coefs']);

// 3. Export each of your 3 custom bounding regions 
Export.image.toAsset({
  image: selected_ccdc,
  description: 'CCDC_TL_2020_2025_20m_westbox',
  assetId: 'projects/ee-rg2v2/assets/CCDC_TL_2020_2025_20m_westbox',
  region: west_box,
  scale: 20,
  crs: 'EPSG:32751',
  maxPixels: 1e13,
  pyramidingPolicy: { '.default': 'sample' }
});

