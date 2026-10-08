/*
TIME SERIES CLASSIFICATION ASSEMBLY without masked layer, temporal consistency
*/
/////////1. Setup/////////
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry()
//post processing module
var post = require('users/rg2icraf/Luma:post_process')
var clf = require('users/rg2icraf/Luma:classify_model')
//Probability based land cover
var lc_prob_2025 = ee.Image('projects/ee-rg2/assets/LC_probs_2025_TL_CCDC').divide(10000)
var lc_prob_2020 = ee.Image('projects/ee-rg2/assets/LC_probs_2020_TL_CCDC').divide(10000)

/////////2. Define Land Cover Class and Smoothing Parameter/////////
//class name
var class_bands = ee.List(['Primary_forest', 'Secondary_forest', 'Mangrove_forest' ,'Coastal_forest',
                  'Teak_plantation', 'Eucalyptus_plantation', 'Coffe_agroforestry', 'Cocoa_agroforestry', 
                  'Coconut_agroforestry', 'Mixed_garden', 'Paddy_field', 'Other_cropland', 'Grassland', 'Shrubland', 
                  'Cleared_land', 'Built-up', 'Waterbody', 'Dry_riverbed']);
//class id
var class_value = ee.List([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18])
//Smoothing_parameter
var smoothness = {
  'Primary_forest': 4.5,       
  'Secondary_forest': 4,    
  'Mangrove_forest': 2.5,
  'Coastal_forest': 2.7,
  'Teak_plantation': 1.0,
  'Eucalyptus_plantation': 1.0,
  'Coffe_agroforestry': 1.0,   
  'Cocoa_agroforestry': 1.0,
  'Coconut_agroforestry': 1.0,
  'Mixed_garden': 0.8,
  'Paddy_field': 0.8,        
  'Other_cropland': 0.8,
  'Grassland': 2.5,
  'Shrubland': 2.5,
  'Cleared_land': 1.4,
  'Built-up': 1,            
  'Waterbody': 1.5,
  'Dry_riverbed': 0.5         
}

/////////2. Implement post-processing bayesian smoothing/////////

//Original landcover (unsmooth), from argmax
var landcover_2025 = post.classifyFromProbs(lc_prob_2025, class_bands, class_value);
//smoothing applied to probability images
var smoothedProbs_2025 = post.bayesianSmooth(lc_prob_2025,class_bands,
                      smoothness,
                      7,    //Window size
                      0.7   //Neighbor fraction
                    );
//generate land cover map using argmax
var lc_smooth_2025 = post.classifyFromProbs(smoothedProbs_2025, class_bands, class_value);

//2020 land cover
//2020 probabilitty
var landcover_2020 = post.classifyFromProbs(lc_prob_2020, class_bands, class_value);
//bayesian smooth for the probability
var smoothedProbs_2020 = post.bayesianSmooth(lc_prob_2020,class_bands,
                      smoothness,
                      7,    // Window size (5x5)
                      0.7   // Neighbor fraction
                    );
//generate the final mal using smoothed probability
var lc_smooth_2020 = post.classifyFromProbs(smoothedProbs_2020, class_bands, class_value);

/////////3. Adding Land Cover to layer/////////
var classVis = {
  min: 1,
  max: 18,
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
    '#A1D76A', // 10  Mixed_garden
    '#C7E9B4', // 11 Paddy_field
    '#FDD835', // 12 Other_cropland
    '#D9EF8B', // 13 Grassland
    '#C2A878', // 14 Shrubland
    '#d7d6d6', // 15 Cleared_land
    '#DB310D', // 16 Built-up
    '#0e8acc', // 17 Waterbody
    '#9ECAE1'  // 18 Dry_riverbed
  ]
};

//Map.addLayer(landcover_2025, classVis, 'Raw Land Cover (argmax')
Map.addLayer(lc_smooth_2020, classVis, 'Smooth Land Cover (argmax) 2020')
Map.addLayer(lc_smooth_2025, classVis, 'Smooth Land Cover (argmax) 2025')


/////// 4. CONSISTENCY ///////
/////// STABLE MASK ///////
/*
// CCDC assets (same ones used for the feature stack)
var tl_ccdc1 = ee.Image('projects/ee-rg2v2/assets/CCDC_TL_2020_2025_20m_westbox');
var tl_ccdc2 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_east-timor');
var tl_ccdc3 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_off_dili');
var tl_ccdc4 = ee.Image('projects/gee-v1-510606/assets/CCDC_TL_2020_2025_20m_central-timor');
var tl_ccdc5 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_east_box');
var ccdc_mosaic = ee.ImageCollection([tl_ccdc1, tl_ccdc2, tl_ccdc3, tl_ccdc4, tl_ccdc5]).mosaic();

// Window slightly wider than the map dates (break dates are only accurate to about a month)
var T0 = 2020.2, T1 = 2025.6;

var nBreaks = ccdc_mosaic.select('tBreak')
  .gte(T0).and(ccdc_mosaic.select('tBreak').lte(T1))
  .arrayReduce(ee.Reducer.sum(), [0])
  .arrayProject([0]).arrayFlatten([['n_breaks']])
  .unmask(1);                                  // no CCDC data = not stable

// Independent evidence of change: Hansen loss 2021-2025
var hansen = ee.Image('UMD/hansen/global_forest_change_2025_v1_13');
var ly = hansen.select('lossyear');
var hansenLoss = ly.gte(21).and(ly.lte(25)).unmask(0);

var stable = nBreaks.eq(0).and(hansenLoss.not())
  .rename('stable').toByte();

// Share of country flagged stable (should be close to what you saw before, slightly lower)
print('Share stable', stable.reduceRegion({
  reducer: ee.Reducer.mean(), geometry: aoi, scale: 100,
  maxPixels: 1e10, tileScale: 8, bestEffort: true}));

// Export once, then read it back to avoid recomputing the CCDC arrays
Export.image.toAsset({
  image: stable,
  description: 'TL_stable_mask_2020_2025',
  assetId: 'projects/ee-rg2/assets/TL_stable_mask_2020_2025',
  region: aoi, scale: 20, crs: 'EPSG:32751',
  maxPixels: 1e13,
  pyramidingPolicy: {'.default': 'mode'}
});
*/
//Labels from the smoothed probabilities, water merged first
var stable = ee.Image('projects/ee-rg2/assets/TL_stable_mask_2020_2025')
var lc20 = post.classifyFromProbs(smoothedProbs_2020, class_bands, class_value).toInt16();
var lc25 = post.classifyFromProbs(smoothedProbs_2025, class_bands, class_value).toInt16();
lc20 = lc20.where(lc20.eq(18), 17);
lc25 = lc25.where(lc25.eq(18), 17);
//MMU filter
//Set Minimum Mapping Unit (e.g., 0.5 ha = 50 pixels at 10m resolution)
//(0 = never remove)
var mmuTable = {     
  1: 50, 
  2: 50, 
  3: 50, 
  4: 50,
  5: 50, 
  6: 0, 
  7: 50, 
  8: 0, 
  9: 20, 
  10: 0,
  11: 20, 
  12: 20, 
  13: 50, 
  14: 50,
  15: 20, 
  16: 10, 
  17: 0
};

function applyMMU(lc, table, iterations) {
  var ids = Object.keys(table).map(Number);
  var vals = ids.map(function(k) { return table[k]; });
  var maxMmu = Math.max.apply(null, vals);

  var size = lc.connectedPixelCount(maxMmu + 1, true);
  var minSize = lc.remap(ids, vals, 0);          // each pixel's own-class threshold
  var small = size.lt(minSize);

  var kept = lc.updateMask(small.not());         // large patches only
  var filled = kept.focalMode({radius: 1, kernelType: 'square',
                               units: 'pixels', iterations: iterations || 6});
  // kept pixels on top, then filled values, then the original label as a fallback
  return lc.blend(filled).blend(kept).rename('classification').toByte();
}
//if apply MMU rule before enforcing consistency rule, uncomment
//lc20 = applyMMU(lc20, mmuTable);    
//lc25 = applyMMU(lc25, mmuTable);


// Rule groups: [class in 2020, class in 2025]
//ecologically impossible change
var impossible = [
  [2, 1], //secondary to primary forest
  [7, 1], //coffee to primary forest
  [14, 1],  // Shrubland -> Primary Forest
  [13, 1],  // Grassland -> Primary Forest
  [12, 1],  // Other Cropland -> Primary Forest
  [15, 1],  // Cleared Land -> Primary Forest
  [16, 1],  // Built-up -> Primary Forest
  [17, 1]   // Waterbody -> Primary Forest
  ]; 
//still possible change
var unlikely = [[1, 2], [2, 14], [14, 2], [1, 14], [13, 14],
                [12, 14], [11, 12], [12, 5], [9, 4], [6, 14], [11, 4]];

function enforceConsistency(lc2020, lc2025, invalid, stableMask) {
  var codes = invalid.map(function(p) { return p[0] * 100 + p[1]; });
  var pair = lc2020.toInt16().multiply(100).add(lc2025.toInt16()).rename('pair');
  var flag = pair.remap(codes, ee.List.repeat(1, codes.length), 0).rename('flag');
  if (stableMask) { flag = flag.and(stableMask); }
  return {
    label: lc2020.where(flag.eq(1), lc2025).toInt16(),
    flag: flag,
    pair: pair
  };
}

// 1. impossible transitions: everywhere
var step1 = enforceConsistency(lc20, lc25, impossible);
// 2. unlikely transitions: only where the stable mask says nothing changed
var step2 = enforceConsistency(step1.label, lc25, unlikely, stable);

var lc_final_2020 = step2.label.rename('classification').toByte();
var lc_final_2025 = lc25.rename('classification').toByte();

/////// C. REPORT ///////
function count(res, label) {
  print(label, res.pair.updateMask(res.flag.eq(1)).reduceRegion({
    reducer: ee.Reducer.frequencyHistogram(), geometry: aoi, scale: 100,
    maxPixels: 1e10, bestEffort: true, tileScale: 8}).get('pair'));
}
count(step1, 'Reclassified (impossible), by transition code');
count(step2, 'Reclassified (unlikely + stable), by transition code');

print('Share of pixels that changed class, before', lc20.neq(lc25).reduceRegion({
  reducer: ee.Reducer.mean(), geometry: aoi, scale: 100,
  maxPixels: 1e10, bestEffort: true, tileScale: 8}));
print('Share of pixels that changed class, after', lc_final_2020.neq(lc_final_2025).reduceRegion({
  reducer: ee.Reducer.mean(), geometry: aoi, scale: 100,
  maxPixels: 1e10, bestEffort: true, tileScale: 8}));

Map.addLayer(lc_final_2020, classVis, 'Final LC 2020');
Map.addLayer(lc_final_2025, classVis, 'Final LC 2025');

print(lc_final_2025)
//Export if already goood
Export.image.toDrive({
    image: lc_final_2020 , 
    description: 'LULC_S2_2020', 
    scale: 10, 
    region: aoi, 
    folder: 'AFCLIM_LC', 
    crs: 'EPSG:32751',
    fileFormat: 'GeoTIFF', 
    maxPixels: 1e13,
    formatOptions: {
    cloudOptimized: true
      },
      });
  
Export.image.toDrive({
    image: lc_final_2025 , 
    description: 'LULC_S2_2025_smoothing', 
    scale: 10, 
    region: aoi, 
    folder: 'AFCLIM_LC', 
    crs: 'EPSG:32751',
    fileFormat: 'GeoTIFF', 
    maxPixels: 1e13,
    formatOptions: {
    cloudOptimized: true
      },
      });
  