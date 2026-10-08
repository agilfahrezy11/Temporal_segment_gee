//AOI
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry()
var start = '2020-01-01';
var end   = '2020-12-30'; 

/////////1. SENTINEL-2 IMAGERY/////////
//Threshold for cloud masking
var CLEAR_THRESHOLD = 0.60; 
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
//Cloud masking
function maskAndScaleS2(img) {
  var cloudMask = img.select('cs_cdf').gte(CLEAR_THRESHOLD);
  var scaled = img.select('B.*').divide(10000);
  return scaled
    .updateMask(cloudMask)
    .copyProperties(img, ['system:time_start', 'system:time_end']);
}
//Sentinel 2 Image Collection
var s2_coll = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(start, end)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
  .linkCollection(csPlus, ['cs_cdf'])
  .map(maskAndScaleS2)
  
var spectral_band = s2_coll.median()
                  .clip(aoi)
                  .select('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12')
Map.addLayer(spectral_band)
//Dynamic Spectral Index
var spectral = require("users/dmlmont/spectral:spectral");
function addIndices(img) {
  var param = {
    "B": img.select("B2"),
    "G": img.select("B3"),
    "R": img.select("B4"),
    "N": img.select("B8"),
    "S1": img.select("B11"),
    "S2": img.select("B12"), 
    "C1": 6,
    "C2": 7.5,
    "g": 2.5,
    "L":1
  };
  var indexList = ["NDVI", "MNDWI", "EVI", "NDMI"];
  var indices = spectral.computeIndex(img, indexList, param);
  return img.addBands(indices);
}
//Map index across image collection
var index_coll = s2_coll.map(addIndices)
//Retrieve index list
var indexList = ["NDVI", "MNDWI", "EVI", "NDMI"];
// Combine multiple statistical reducers into a single execution
var combinedReducer = ee.Reducer.mean()
  .combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true})
  .combine({reducer2: ee.Reducer.percentile([10, 90]), sharedInputs: true});
// Apply reducer over time-series
var index_summary = index_coll.select(indexList).reduce(combinedReducer);
var temporal_index = index_summary.clip(aoi);

//Static Spectral Index
var stat_param = {
  "B": spectral_band.select("B2"),
  "N": spectral_band.select("B8"), 
  "R": spectral_band.select("B4"), 
  "RE1": spectral_band.select("B5"), 
  "RE2": spectral_band.select("B6"),
  "S1": spectral_band.select("B11"), 
  "G": spectral_band.select("B3"),
  "S2": spectral_band.select("B12"), 
  "L": 1
}
//initiate calculation
//define selected index
var stat_list = ["RENDVI","SAVI","DBSI", "MBI", "BLFEI", "MCARI"];
var static_index = spectral.computeIndex(spectral_band, stat_list, stat_param).select("RENDVI","SAVI","DBSI", "MBI", "BLFEI", "MCARI");

var index_feature = temporal_index.addBands(static_index).toFloat() 

print('Spectral Index Features', index_feature) //(3. Spectral Index Feature)


var temporal_index = ["MNDWI_mean", "MNDWI_p90", "NDMI_p90", "NDVI_mean", "NDVI_p10", "NDVI_p90"]
var annual_index = ["BLFEI", "DBSI", "MBI","RENDVI", "SAVI"]
var ccdc_index = ["NDVI", "NDBI", "NDWI"]
//select index based on rfecv result
var final_temporal_index_2020 = index_feature.select(temporal_index)
var final_static_index_2020 = index_feature.select(annual_index)

print("Temporal Index 2020", final_temporal_index_2020)
print("Annual Index 2020", final_static_index_2020)

Export.image.toAsset({
  image: final_temporal_index_2020,
  description: 'temporal_index_2020',
  assetId: 'projects/ee-rg2/assets/temporal_indices_2020_TL',
  region: aoi,
  scale: 10,
  crs: 'EPSG:32751',
  maxPixels: 1e13,
  pyramidingPolicy: {'.default': 'mean'}
});
Export.image.toAsset({
  image: final_static_index_2020,
  description: 'annual_index_2020',
  assetId: 'projects/ee-rg2/assets/annual_indices_2020_TL',
  region: aoi,
  scale: 10,
  crs: 'EPSG:32751',
  maxPixels: 1e13,
  pyramidingPolicy: {'.default': 'mean'}
});