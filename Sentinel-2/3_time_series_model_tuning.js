//define the AOI
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry()

//import the required modules
var ccdc = require("users/rg2icraf/Luma:ccdc_sentinel2");
var extraction = require("users/rg2icraf/Luma:extract_split_data");
var clf = require("users/rg2icraf/Luma:classify_model")
var tuning = require("users/rg2icraf/Luma:tune_model")
var post = require("users/rg2icraf/Luma:post_process")

//rebuild the feature stacks for each year
//optical CCDC
var optical_2020 = ee.Image('projects/ee-rg2/assets/TL_ccdc_features_2020').clip(aoi)
var optical_2025 = ee.Image('projects/ee-rg2/assets/TL_ccdc_features_2025').clip(aoi)


//topography and distance: Use by both year
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
var topo_dist_features = ee.Image.cat([elev, slope, landform, twi, dist_urban, dist_road]);


//spectral indices 
//2025
var index1 = ee.Image('projects/ee-rg2/assets/index_features_TL_2025')
var index2 = ee.Image('projects/ee-v2-508913/assets/MBI_s2_tl')
var index3 = ee.Image('projects/ee-rg2/assets/distance_TL_additional_index_TL').select('SAVI', 'DBSI')
var ccdc_index_2025 = ee.Image('projects/ee-rg2/assets/TL_ccdc_indices_2025')
var index_features_2025 = ee.Image.cat([index1, index2, index3, ccdc_index_2025])

//2020
var temporal_index_2020 = ee.Image('projects/ee-rg2/assets/temporal_indices_2020_TL')
var annual_index_2020 = ee.Image('projects/ee-rg2/assets/annual_indices_2020_TL')
var ccdc_index_2020 = ee.Image('projects/ee-rg2/assets/TL_ccdc_indices_2020')
var index_features_2020 = ee.Image.cat([temporal_index_2020, annual_index_2020, ccdc_index_2020])

//radar backscatter
//2025
var s1_seasonal_2025 = ee.Image('projects/ee-v2-508913/assets/s1_seasonal_tl_25')

//2020
//not in asset, rebuild it
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

var s1_seasonal_2020 = getSeasonalS1(2020, aoi).clip(aoi);

var final_stack_2020 = ee.Image.cat([
                      optical_2020,         //CCDC optical band
                      index_features_2020, //spectral indices
                      s1_seasonal_2020,   //radar backscatter
                      topo_dist_features //topography and distance
                      ])
                    
print('stack for 2020', final_stack_2020)
var final_stack_2025 = ee.Image.cat([
                      optical_2025,         //CCDC optical band
                      index_features_2025, //spectral indices
                      s1_seasonal_2025,   //radar backscatter
                      topo_dist_features //topography and distance
                      ])
print('stack for 2025', final_stack_2025)
//print('Original feature stack', feature_stack_2025)

//Feature Split
//define filtered out samples
var roi = ee.FeatureCollection('projects/ee-rg2/assets/filtered_sample_2025_revised')
//Execute Stratified Split
//2025
var samples_2025 = extraction.stratifiedSplit(roi, final_stack_2025, 'new_id', 10, 0.7, 42);
var train_pix_2025 = samples_2025.trainingPixels;
var test_pix_2025  = samples_2025.testingPixels;

/*
//print('Training pixel count:', train_pix_2025.size());
//print('Testing pixel count:', test_pix_2025.size());

var num_trees = ee.List.sequence(200, 500, 100); //(min, max, sequence)
var variablesplit = ee.List([4, 5, 6, 8, 9]);
/*
//tuning code
var accuracies = numTreesList.map(function(numTrees){
  return variablesplit.map(function(variablesPerSplit){
    var classifier = ee.Classifier.smileRandomForest({
      numberOfTrees: numTrees, 
      variablesPerSplit: variablesPerSplit 
    })
    .train({
      features: trainingPixels, 
      classProperty: 'ID',
      inputProperties: input_features.bandNames() 
    });
    var accuracy = testingPixels 
    .classify(classifier) 
    .errorMatrix('ID', 'classification')
    .accuracy();
    return ee.Feature(null, {'accuracy': accuracy,
      'numberOfTrees': numTrees, 
      'variableSplit': variablesPerSplit 
    })
  })
}).flatten()
*/

/*
//use the module
var tuning_result = tuning.tuneRandomForest(
  train_pix_2025, //training
  test_pix_2025, //test
  'new_id', //attribute table
  final_stack_2025.bandNames(), //input
  num_trees, //number of trees list
  variablesplit, //number of varible 
  0 );


print('Tuning results', tuning_result);
var bestParams = tuning.getBestParams(tuning_result, 'accuracy');
print('Best parameters:', bestParams);
//best parameter: 

/Trees: 400, vsplit: 6
//oa: 0.814
*/

//perform multiprobability classification
//class name
var classIdToName = {
  1: 'Primary_forest',
  2: 'Secondary_forest',
  3: 'Mangrove_forest',
  4: 'Coastal_forest',
  5: 'Teak_plantation',
  6: 'Eucalyptus_plantation',
  7: 'Coffe_agroforestry',
  8: 'Cocoa_agroforestry',
  9: 'Coconut_agroforestry',
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
//perform classification on 2025 data
var rf_model_2025 = clf.multiProbabilityClassification(
  train_pix_2025,
  'new_id',
  classIdToName,
  final_stack_2025,
  { nTrees: 400, vSplit: 6, minLeaf: 1, seed: 42 });

//define result
var probsImage = rf_model_2025.probsImage;
var classBands = rf_model_2025.classBands;
var classValues = rf_model_2025.classValues;
/*
var selected = ["BLFEI", "DBSI", "MBI", "MNDWI_mean", "MNDWI_p90", "NDBI", "NDMI_p90", "NDVI", "NDVI_mean", "NDVI_p10", "NDVI_p90", "NDWI", "RENDVI",
"SAVI", "TWI", "VH_dry",
"VH_wet", "VV_dry", "VV_wet", "blue", "dist_road", "dist_urban", "elevation", "green", "landform", "nir", 
"nir_amplitude", "nir_phase", "red", "slope", "swir1", "swir2"];
var testArea = ee.Geometry.BBox(125.38, -8.76, 125.41, -8.73);
var model = rf_model_2025.trainedModel;

var s20_a = final_stack_2020.select(selected);                 // 2020 order
var s20_b = final_stack_2020.select(final_stack_2025.select(selected).bandNames());  // 2025 order

function lab(img) { return img.classify(model).arrayArgmax().arrayGet([0]); }
var diff = lab(s20_a).neq(lab(s20_b));

var pt = ee.Geometry.Point([125.57, -8.55]);   // any land point
print(final_stack_2020.select(selected).reduceRegion({reducer: ee.Reducer.first(), geometry: pt, scale: 10}));
print(final_stack_2025.select(selected).reduceRegion({reducer: ee.Reducer.first(), geometry: pt, scale: 10}));
print('Share of test-area pixels where label differs',
  diff.reduceRegion({reducer: ee.Reducer.mean(), geometry: testArea, scale: 10, maxPixels: 1e9}));
*/
//print('Prob band names', probsImage.bandNames());
var smoothnessTL = {
  'Primary_forest': 4.5,
  'Secondary_forest': 3.0,
  'Mangrove_forest': 3.0,
  'Coastal_forest': 3.5,
  'Teak_plantation': 1.2,
  'Eucalyptus_plantation': 1.2,
  'Coffe_agroforestry': 1.5,
  'Cocoa_agroforestry': 1.5,
  'Coconut_agroforestry': 1.5,
  'Mixed_garden': 1.2,
  'Paddy_field': 1.0,
  'Other_cropland': 1.0,
  'Grassland': 1.5,
  'Shrubland': 1.5,
  'Cleared_land': 0.5,
  'Built-up': 0.2,
  'Waterbody': 4.0,
  'Dry_riverbed': 0.1
}

//Testing
var testArea = ee.Geometry.BBox(125.38, -8.76, 125.41, -8.73);
var testProbs = probsImage.clip(testArea);
// Fast local variance check
var localVarTest = post.localVariance(testProbs, classBands, 7, 0.5);
print('Test Area Variance:', localVarTest.reduceRegion({
  reducer: ee.Reducer.percentile([50, 75, 90]),
  geometry: testArea,
  scale: 10,
  maxPixels: 1e7
}));
// Test smoothing locally
var testSmoothed = post.bayesianSmooth(testProbs, classBands, smoothnessTL, 7, 0.5);
var testLabel = post.classifyFromProbs(testSmoothed, classBands, classValues);
Map.centerObject(testArea, 10)
Map.addLayer(testLabel, {min: 1, max: 18}, 'Smoothed Test Area');

//Applied 2025 model on 2020 data
//get the trained model
var trainedModel2025 = rf_model_2025.trainedModel;
var probs_2020 = final_stack_2020.classify(trainedModel2025);
// Parse class IDs in ascending numerical order to maintain band order alignment
var sortedIds = Object.keys(classIdToName).map(Number).sort(function(a, b) { return a - b; });
var classBands = sortedIds.map(function(id) { return classIdToName[id]; });
// Flatten array image into separate class probability bands for 2020
var probsImage_2020 = probs_2020.rename('probs').arrayFlatten([classBands]);
var smoothedProbs = post.bayesianSmooth(
  probsImage,
  classBands,
  smoothnessTL,
  5,    //Window size (5x5)
  0.5   //Neighbor fraction
);