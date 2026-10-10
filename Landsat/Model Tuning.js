//define AOI
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry();

//define module
var extraction = require("users/rg2icraf/Luma:extract_split_data");
var clf = require("users/rg2icraf/Luma:classify_model")
var tuning = require("users/rg2icraf/Luma:tune_model")
var post = require("users/rg2icraf/Luma:post_process")

//Predefine asset relevant for spectral data
var index_dry = ee.Image('projects/gee-v1-510606/assets/Landsat_TL_Dry_IndexFeature_2025')
var index_wet = ee.Image('projects/gee-v1-510606/assets/Landsat_TL_Wet_IndexFeature_2025')
var landsat_dry = ee.Image('projects/gee-v1-510606/assets/Landsat_TL_Dry_Season_2025')
var landsat_wet = ee.Image('projects/gee-v1-510606/assets/Landsat_TL_Wet_Season_2025')

//topographical data
var elev = ee.Image("NASA/NASADEM_HGT/001")
            .select('elevation')
            .clip(aoi)
var terrain_product = ee.Terrain.products(elev).select('slope').toFloat();

//stack all data
var feature_stack = ee.Image.cat([
  landsat_wet,
  landsat_dry,
  index_dry,
  index_wet,
  elev,
  terrain_product
  ]);
//Sample data
print('feature_stack 2025', feature_stack)
var roi = ee.FeatureCollection('projects/ee-rg2/assets/Samples/Samples_TL_hierarchy_2025_v1');
var samples_2025 = extraction.stratifiedSplit(roi, feature_stack, 'new_id', 30, 0.7, 42);
var train_pix_2025 = samples_2025.trainingPixels;
var test_pix_2025  = samples_2025.testingPixels;

/*
//Define the parameter space
var num_trees = ee.List([100, 200, 400, 500]); 
var variablesplit = ee.List([ 9, 11, 15, 19]);

//use the module
var tuning_result = tuning.tuneRandomForest(
  train_pix_2025, //training
  test_pix_2025, //test
  'new_id', //attribute table
  feature_stack.bandNames(), //input
  num_trees, //number of trees list
  variablesplit, //number of varible per split
  0 );

//print('Tuning results', tuning_result);
var bestParams = tuning.getBestParams(tuning_result, 'accuracy');
print('Best parameters:', bestParams);

/*
//best parameter: 
accuracy: 0.7970687711386697
kappa: 0.780727776129734
minLeafPopulation: 1
numberOfTrees: 500
variablesPerSplit: 19
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
  feature_stack,
  { nTrees: 500, vSplit: 19, minLeaf: 1, seed: 42 });


//define result
var probsImage = rf_model_2025.probsImage;
var classBands = rf_model_2025.classBands;
var classValues = rf_model_2025.classValues;
//Probability evaluation
var eval_prob = clf.evaluateProbabilities(rf_model_2025.trainedModel, test_pix_2025, 'new_id', classValues)
print('Test log loss', eval_prob.logLoss);
print('Test top-class accuracy', eval_prob.accuracy);
print('Test rows used', eval_prob.n);
print('Random-guess log loss (ln 18)', Math.log(18));



//Build Feature stack for 2015 and 2005
var index_wet_05 = ee.Image('projects/gee-v1-510606/assets/Landsat_TL_Index_Wet_2005')
var index_dry_05 = ee.Image('projects/gee-v1-510606/assets/Landsat_TL_Index_Dry_2005')
var landsat_dry_05 = ee.Image('projects/gee-v1-510606/assets/Landsat_TL_Dry_Season_2005')
var landsat_wet_05 = ee.Image('projects/gee-v1-510606/assets/Landsat_TL_Wet_Season_2005')
//stack 2005
var feature_2005 = ee.Image.cat([
  landsat_dry_05,
  landsat_wet_05,
  index_dry_05,
  index_wet_05,
  elev,
  terrain_product])

//stack 2015
var landsat_dry_15 = ee.Image('projects/gee-v1-510606/assets/Landsat_TL_Dry_Season_2015')
var landsat_wet_15 = ee.Image('projects/gee-v1-510606/assets/Landsat_TL_Wet_Season_2015')
var index_wet_15 = ee.Image('projects/gee-v1-510606/assets/Landsat_TL_Index_Wet_2015')
var index_dry_15 = ee.Image('projects/gee-v1-510606/assets/Landsat_TL_Index_Dry_2015')

//stack 2015
var feature_2015 = ee.Image.cat([
  landsat_dry_15,
  landsat_wet_15,
  index_dry_15,
  index_wet_15,
  elev,
  terrain_product])

print(feature_2005)
print(feature_2015)


//Applied 2025 model on 2015 and 2005
//get the trained model
var trainedModel2025 = rf_model_2025.trainedModel;
//Applied trained model
var probs_2015 = feature_2015.classify(trainedModel2025);
var probs_2005 = feature_2005.classify(trainedModel2025)

// Parse class IDs in ascending numerical order to maintain band order alignment
var sortedIds = Object.keys(classIdToName).map(Number).sort(function(a, b) { return a - b; });
var classBands = sortedIds.map(function(id) { return classIdToName[id]; });
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

//Flatten array image into separate class probability bands for all year
var smoothedProbs_2025 = post.bayesianSmooth(probsImage, classBands,
  smoothnessTL,
  5,    
  0.5   
);
//print('Smoother Probability 2025', smoothedProbs_2025 )
var probsImage_2015 = probs_2015.rename('probs').arrayFlatten([classBands]);
var smoothedProbs_2015 = post.bayesianSmooth(probsImage_2015, classBands,
  smoothnessTL,
  5,    
  0.5   
);
//print('Smoother Probability 2015', smoothedProbs_2015 )
var probsImage_2005 = probs_2005.rename('probs').arrayFlatten([classBands]);
var smoothedProbs_2005 = post.bayesianSmooth(probsImage_2005, classBands,
  smoothnessTL,
  5,    
  0.5   
);
//print('Smoother Probability 2005', smoothedProbs_2005 )



// EXPORT SMOOTHED PROBABILITIES TO ASSETS
//.multiply(10000).round().toUint16() dont forget to rescale the prob

Export.image.toAsset({
  image: smoothedProbs_2025.multiply(10000).round().toUint16(),
  description: 'smoothed_probabilities_2025_TL',
  assetId: 'projects/ee-rg2/assets/LC_prob_2025_Landsat_TL',
  region: aoi,
  scale: 30,
  crs: 'EPSG:32751',
  maxPixels: 1e13,
  pyramidingPolicy: {'.default': 'mean'}
});

Export.image.toAsset({
  image: smoothedProbs_2015.multiply(10000).round().toUint16(),
  description: 'smoothed_probabilities_2015_TL',
  assetId: 'projects/ee-rg2/assets/LC_prob_2015_Landsat_TL',
  region: aoi,
  scale: 30,
  crs: 'EPSG:32751',
  maxPixels: 1e13,
  pyramidingPolicy: {'.default': 'mean'}
});

Export.image.toAsset({
  image: smoothedProbs_2005.multiply(10000).round().toUint16(),
  description: 'smoothed_probabilities_2005_TL',
  assetId: 'projects/ee-rg2/assets/LC_prob_2005_Landsat_TL',
  region: aoi,
  scale: 30,
  crs: 'EPSG:32751',
  maxPixels: 1e13,
  pyramidingPolicy: {'.default': 'mean'}
});
