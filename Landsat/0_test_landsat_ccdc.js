var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km')
var landsat_ccdc = require('users/rg2icraf/Luma:ccdc_landsat')
var google_ccdc = ee.ImageCollection("GOOGLE/GLOBAL_CCDC/V1")


var test = ee.Image('projects/ee-rg2v2/assets/CCDC_TL_2020_2025_20m_test');
var testArea = test.geometry();
Map.centerObject(testArea, 11);

var testPoint = ee.Geometry.Point([125.5736, -8.5568]); 
var testRoi = testPoint.buffer(5000); 
var landsat_ts = landsat_ccdc.getLandsatTimeSeries(testRoi, '2020-12-01', '2025-12-31');
var ccdcTest = landsat_ccdc.runLandsatCCDC(landsat_ts);

// Generate synthetic stacks for 2019 and 2024 (Fraction 0.5 = July)
var synth2019 = landsat_ccdc.getSyntheticLandsatStack(ccdcTest, 2020.5);
var synth2024 = landsat_ccdc.getSyntheticLandsatStack(ccdcTest, 2024.5);

// Extract original surface reflectance and predicted NIR trajectory
var nirObs = landsat_ts.select('NIR').map(function(img) {
  return img.select('NIR').divide(10000)
    .copyProperties(img, ['system:time_start']);
});
// Plot time-series trajectory at test point
var chart = ui.Chart.image.series(nirObs, testPoint, ee.Reducer.mean(), 30)
  .setOptions({
    title: 'Landsat NIR Raw Observations over Time',
    hAxis: {title: 'Date'},
    vAxis: {title: 'Surface Reflectance'},
    lineWidth: 1,
    pointSize: 3
  });
print(chart);
// Add to Map in True-Color (RGB)
var rgbVis = {bands: ['RED', 'GREEN', 'BLUE'], min: 0.02, max: 0.20};

Map.addLayer(synth2019, rgbVis, 'Synthetic RGB 2019.5');
Map.addLayer(synth2024, rgbVis, 'Synthetic RGB 2024.5');

// Generate the 2020-2025 Landsat 8 & 9 time series (scaled to 0-10000)