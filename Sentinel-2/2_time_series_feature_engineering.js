//Timor leste AOI
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry()

/////////1. CCDC spatiotemporal features/////////
//CCDC source code
var ccdcModule = require('users/rg2icraf/Luma:ccdc_sentinel2');

//define CCDC assets
var tl_ccdc1 = ee.Image('projects/ee-rg2v2/assets/CCDC_TL_2020_2025_20m_westbox')
var tl_ccdc2 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_east-timor')
var tl_ccdc3 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_off_dili')
var tl_ccdc4 = ee.Image('projects/gee-v1-510606/assets/CCDC_TL_2020_2025_20m_central-timor')
var tl_ccdc5 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_east_box')

//combine them into a single mosaic
var ccdc_mosaic = ee.ImageCollection([tl_ccdc1, tl_ccdc2, tl_ccdc3, tl_ccdc4, tl_ccdc5]).mosaic();
print(ccdc_mosaic)
//define target year
var targetYears = [2020, 2025];
var fractionOfYear = 0.5;

//build ccdc imagery
var syntheticStacks = ccdcModule.getMultiYearSyntheticStacks(
  ccdc_mosaic, 
  targetYears, 
  fractionOfYear
);

/////////2. Define predictors/////////
//A. Topography and distance
var elev = ee.Image("NASA/NASADEM_HGT/001")
            .select('elevation')
            .clip(aoi)
var slope = ee.Terrain.products(elev).select('slope').toFloat();
var landform = ee.Image('projects/ee-rg2v2/assets/TPI_landform_TL')
                      .rename('landform')
                      .toByte()
                      .clip(aoi)
var twi = ee.Image('projects/ee-rg2/assets/Top_Wetness_index')
                  .rename('TWI')
                  .toFloat()
var dist_urban = ee.Image('projects/ee-rg2/assets/dist_from_urban').rename('dist_urban').clip(aoi)
var dist_road = ee.Image('projects/ee-rg2/assets/dist_from_mainroad').rename('dist_road').clip(aoi)

//stack topography and distance
var topo_dist_stack = ee.Image.cat([elev, slope, landform, twi, dist_urban, dist_road]);

//B. Spectral Index
var index1 = ee.Image('projects/ee-rg2/assets/index_features_TL_2025')
var index2 = ee.Image('projects/ee-v2-508913/assets/MBI_s2_tl')
var index3 = ee.Image('projects/ee-rg2/assets/distance_TL_additional_index_TL').select('SAVI', 'DBSI')
var index_stack = ee.Image.cat([index1, index2, index3])

//C. Radar backscatter
//generate seasonal radar backscatter for 2020
function getSeasonalS1(year, roi) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = ee.Date.fromYMD(year, 12, 31);

  var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .select(['VV', 'VH']);

  //Monsoonal Splitting
  var wet = s1.filter(ee.Filter.calendarRange(12, 4, 'month')).median();
  var dry = s1.filter(ee.Filter.calendarRange(5, 11, 'month')).median();
  //Rename base bands
  var vv_wet = wet.select('VV').rename('VV_wet');
  var vh_wet = wet.select('VH').rename('VH_wet');

  var vv_dry = dry.select('VV').rename('VV_dry');
  var vh_dry = dry.select('VH').rename('VH_dry');
  return ee.Image.cat([vh_dry, vh_wet, vv_dry, vv_wet]);
} 
//since 2025 already in asset
var s1_seasonal_2025 = ee.Image('projects/ee-v2-508913/assets/s1_seasonal_tl_25')
var s1_seasonal_2020 = getSeasonalS1(2020, aoi).clip(aoi);

/////////3. NIR HARMONICS FROM THE ACTIVE SEGMENT, PER TARGET DATE/////////
function nirHarmonics(ccdcImg, t) {
  t = ee.Number(t);
  // same masking pattern as getSyntheticStack
  var mask2D = ccdcImg.select('tStart').lte(t)
                 .and(ccdcImg.select('tEnd').gte(t))
                 .toArray(0).toArray(1);
  var seg = ccdcImg.select('nir_coefs').arrayMask(mask2D).arraySlice(0, 0, 1); // [1 x 8]
  var a1 = seg.arraySlice(1, 2, 3).arrayProject([0]).arrayFlatten([['a1']]).divide(10000);
  var b1 = seg.arraySlice(1, 3, 4).arrayProject([0]).arrayFlatten([['b1']]).divide(10000);
  return ee.Image.cat([a1.hypot(b1).rename('nir_amplitude'),
                       b1.atan2(a1).rename('nir_phase')]);
}
var harm2020 = nirHarmonics(ccdc_mosaic, 2020.5);
var harm2025 = nirHarmonics(ccdc_mosaic, 2025.5);

/////////4. STABILITY BAND (distance to nearest break, in years)/////////
function breakDistance(ccdcImg, t) {
  return ccdcImg.select('tBreak').subtract(t).abs()
    .arrayReduce(ee.Reducer.min(), [0])
    .arrayProject([0])
    .arrayFlatten([['min_break_dist']]);
}
var bd2025 = breakDistance(ccdc_mosaic, 2025.5);
var bd2020 = breakDistance(ccdc_mosaic, 2020.5);

/////////5. MASTER STACKS/////////
var master_stack_2020 = syntheticStacks[2020]
  .addBands(s1_seasonal_2020).addBands(topo_dist_stack)
  .addBands(index_stack).addBands(harm2020).addBands(bd2020);

var master_stack_2025 = syntheticStacks[2025]
  .addBands(s1_seasonal_2025).addBands(topo_dist_stack)
  .addBands(index_stack).addBands(harm2025).addBands(bd2025);
print('Original Feature Stck',master_stack_2025)

/////////6. SAMPLE, THEN FILTER BY STABILITY/////////
var roi = ee.FeatureCollection('projects/ee-rg2/assets/Samples_TL_hierarchy_2025_v1');

var sampled_2025 = master_stack_2025.sampleRegions({
  collection: roi,
  properties: ['new_id', 'year'],
  scale: 10,
  projection: 'EPSG:32751',
  tileScale: 4,
  geometries: true //if for rfecv, set the geometries to false
});

// keep points with no break within 0.5 year of 2025.5
var stable_2025 = sampled_2025.filter(ee.Filter.gte('min_break_dist', 0.5));

print('Input points', roi.size());
print('Sampled (before stability filter)', sampled_2025.size());
print('Stable', stable_2025.size());

Export.table.toDrive({
  collection: stable_2025,
  description: 'TL_Predictor_Matrix_revised',
  fileFormat: 'CSV'
});
//selected bands (RFECV: 0.710)
var index = ["BLFEI", "DBSI", "MBI", "MNDWI_mean", "MNDWI_p90",  "NDMI_p90", "NDVI_mean", "NDVI_p10", "NDVI_p90", "RENDVI", "SAVI"]
var ccdc_index = ["NDBI", "NDVI", "NDWI"]
//var optical = ["blue","green","nir", "nir_amplitude", "nir_phase", "red","swir1", "swir2"]
var non_spectral = ["dist_road", "dist_urban", "elevation", "TWI","landform","slope",]
var radar = [ "VH_dry", "VH_wet", "VV_dry", "VV_wet",]

var optical = ['blue','green','nir','nir_amplitude','nir_phase','red','swir1','swir2'];

[2020, 2025].forEach(function(year) {
  var harm = (year === 2020) ? harm2020 : harm2025;
  var img = syntheticStacks[year].addBands(harm).select(ccdc_index).toFloat();

  Export.image.toAsset({
    image: img,
    description: 'TL_ccdc_indices_' + year,
    assetId: 'projects/ee-rg2/assets/TL_ccdc_indices_' + year,
    region: aoi,
    scale: 20,
    crs: 'EPSG:32751',
    maxPixels: 1e13,
    pyramidingPolicy: {'.default': 'mean', 'nir_phase': 'sample'}
  });
});
//check
var ccdc_feature = ee.Image('projects/ee-rg2/assets/TL_ccdc_features_2025').clip(aoi)
var vis = {
  min: 0.0,
  max: 0.3,
  bands: ['green', 'red', 'blue'],
};
Map.centerObject(aoi, 7)
Map.addLayer(ccdc_feature, vis, 'CCDC Feature 2025')
//export stable/filtered training data
var stable_ids = sampled_2025
  .filter(ee.Filter.gte('min_break_dist', 0.5))
  .select(['new_id', 'new_class']);   // geometry is kept

print('Stable points', stable_ids.size());

//export stable/filtered training data to asset
Export.table.toAsset({
  collection: stable_ids,
  description: 'Stable_training_data_TL_2025_revised',
  assetId: 'projects/ee-rg2/assets/filtered_sample_2025_revised'
});
