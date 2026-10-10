/*
TIME SERIES CLASSIFICATION ASSEMBLY without masked layer, temporal consistency
*/
/////////1. Setup/////////
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry()
//post processing module
var post = require('users/rg2icraf/Luma:post_process')
var clf = require('users/rg2icraf/Luma:classify_model')
//Probability based land cover
//dont forget to modified
var lc_prob_2025 = ee.Image('projects/ee-rg2/assets/TL_Prob/LC_prob_2025_Landsat_TL').divide(10000)
var lc_prob_2015 = ee.Image('projects/ee-rg2/assets/TL_Prob/LC_prob_2015_Landsat_TL').divide(10000)
var lc_prob_2005 = ee.Image('projects/ee-rg2/assets/TL_Prob/LC_prob_2005_Landsat_TL').divide(10000)
/////////2. Define Land Cover Class and Smoothing Parameter/////////
//class name
var class_bands = ee.List(['Primary_forest', 'Secondary_forest', 'Mangrove_forest' ,'Coastal_forest',
                  'Teak_plantation', 'Eucalyptus_plantation', 'Coffe_agroforestry', 'Cocoa_agroforestry', 
                  'Coconut_agroforestry', 'Mixed_garden', 'Paddy_field', 'Other_cropland', 'Grassland', 'Shrubland', 
                  'Cleared_land', 'Built-up', 'Waterbody', 'Dry_riverbed']);

var class_value = ee.List([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18])
//Smoothing_parameter
var smoothness = {
  'Primary_forest': 4.5,       
  'Secondary_forest': 4,    
  'Mangrove_forest': 2.5,
  'Coastal_forest': 2.7,
  'Teak_plantation': 2,
  'Eucalyptus_plantation': 1.5,
  'Coffe_agroforestry': 2,   
  'Cocoa_agroforestry': 1.7,
  'Coconut_agroforestry': 2,
  'Mixed_garden': 1.5,
  'Paddy_field': 2.5,        
  'Other_cropland': 2,
  'Grassland': 2.5,
  'Shrubland': 2.5,
  'Cleared_land': 2,
  'Built-up': 2,            
  'Waterbody': 2,
  'Dry_riverbed': 1.3         
}

//Generate Land Cover 2025
var landcover_2025_raw = post.classifyFromProbs(lc_prob_2025, class_bands, class_value);
//smoothing applied to probability images
var smooth_probs_2025 = post.bayesianSmooth(lc_prob_2025,class_bands,smoothness, 7, 0.7  );
//generate land cover from smoothed prob using argmax
var lc_smooth_2025 = post.classifyFromProbs(smooth_probs_2025, class_bands, class_value);

//Generate Land Cover 2015
var landcover_2015_raw = post.classifyFromProbs(lc_prob_2015, class_bands, class_value);
//smoothing applied to probability images
var smooth_probs_2015 = post.bayesianSmooth(lc_prob_2015, class_bands, smoothness, 7, 0.7);
var lc_smooth_2015 = post.classifyFromProbs(smooth_probs_2015, class_bands, class_value);

//Generate Land Cover 2005
var landcover_2005 = post.classifyFromProbs(lc_prob_2005, class_bands, class_value);
//smoothing applied to probability images
var smooth_probs_2005 = post.bayesianSmooth(lc_prob_2005,class_bands, smoothness, 7, 0.7);
var lc_smooth_2005 = post.classifyFromProbs(smooth_probs_2005, class_bands, class_value);

//merge the waterbody and dryriver bed
var lc25 = lc_smooth_2025.where(lc_smooth_2025.eq(18), 17).toInt16();
var lc15 = lc_smooth_2015.where(lc_smooth_2015.eq(18), 17).toInt16();
var lc05 = lc_smooth_2005.where(lc_smooth_2005.eq(18), 17).toInt16();


//Identify which the existing transition
var names = {
  1: 'Primary_forest', 2: 'Secondary_forest', 3: 'Mangrove_forest', 4: 'Coastal_forest',
  5: 'Teak_plantation', 6: 'Eucalyptus_plantation', 7: 'Coffe_agroforestry',
  8: 'Cocoa_agroforestry', 9: 'Coconut_agroforestry', 10: 'Mixed_garden',
  11: 'Paddy_field', 12: 'Other_cropland', 13: 'Grassland', 14: 'Shrubland',
  15: 'Cleared_land', 16: 'Built-up', 17: 'Waterbody'
};

function topTransitions(lcA, lcB, region, mask, n) {
  // code = class in A * 100 + class in B, e.g. 201 = Secondary -> Primary
  var pair = lcA.toInt16().multiply(100).add(lcB.toInt16()).rename('pair');
  var changed = lcA.neq(lcB);
  var m = mask ? changed.and(mask) : changed;

  var hist = pair.updateMask(m).reduceRegion({
    reducer: ee.Reducer.frequencyHistogram(),
    geometry: region, scale: 100,
    maxPixels: 1e10, bestEffort: true, tileScale: 8
  }).get('pair');

  ee.Dictionary(hist).evaluate(function(d) {
    var rows = Object.keys(d).map(function(k) {
      var code = +k;
      return {code: code, from: names[Math.floor(code / 100)],
              to: names[code % 100], n: Math.round(d[k])};
    });
    rows.sort(function(a, b) { return b.n - a.n; });
    var total = rows.reduce(function(s, r) { return s + r.n; }, 0);
    rows.forEach(function(r) { r.share_pct = Math.round(1000 * r.n / total) / 10; });
    print('Top transitions (total sampled: ' + total + ')', rows.slice(0, n));
  });
}

var change_0515 = topTransitions(lc05, lc15, aoi, null, 30);
var change_15015 = topTransitions(lc15, lc25, aoi, null, 30);

var classVis = {
  min: 1,
  max: 17,
  palette: [
    '#00441B', // 1  Primary_forest
    '#238B45', // 2  Secondary_forest
    '#006D2C', // 3  Mangrove_forest
    '#66C2A4', // 4  Coastal_forest
    '#8C6D31', // 5  Teak_plantation
    '#A6761D', // 6  Eucalyptus_plantation
    '#7FC97F', // 7  Coffe_agroforestry
    '#66A061', // 8  Cocoa_agroforestry
    '#41AB5D', // 9  Coconut_agroforestry
    '#A1D76A', // 10 Mixed_garden
    '#C7E9B4', // 11 Paddy_field
    '#FDD835', // 12 Other_cropland
    '#D9EF8B', // 13 Grassland
    '#C2A878', // 14 Shrubland
    '#d7d6d6', // 15 Cleared_land
    '#DB310D', // 16 Built-up
    '#0e8acc'  // 17 Water
  ]};
  
Map.centerObject(aoi, 9);
Map.addLayer(lc25, classVis, 'Land Cover 2025');
Map.addLayer(lc15, classVis, 'Land Cover 2015');
Map.addLayer(lc05, classVis, 'Land Cover 2005');