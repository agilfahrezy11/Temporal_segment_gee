//AOI
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry()
var start = '2025-01-01';
var end   = '2025-12-30'; 

//get the modules
var spectral = require("users/dmlmont/spectral:spectral");
var get_landsat = require('users/rg2icraf/Luma:get_landsat');
var ccdc_module = require('users/rg2icraf/Luma:ccdc_landsat');

//generate the synthetic image
//Get the CCDC image 
var landsat_ccdc = ee.Image('projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_TL_2020_2025_Landsat').clip(aoi);
//build ccdc imagery, breakdown the image perseasonality
//(1)
var landsat_wet = ccdc_module.getSyntheticLandsatStack(landsat_ccdc, 2025.65, false) //somewhere in august
var landsat_dry = ccdc_module.getSyntheticLandsatStack(landsat_ccdc, 2025.20, false) //somewhere in march

//Band Harmonics
function getAllHarmonics(ccdcImg, t, bands) {
  var seg = ccdc_module.nearestSegment(ccdcImg, t);
  var imgs = bands.map(function(b) {
    var c = ccdcImg.select(b + '_coefs')
      .arraySlice(0, seg.idx, seg.idx.add(1)).arrayProject([1]);
    var a1 = c.arrayGet([2]).divide(10000);
    var b1 = c.arrayGet([3]).divide(10000);
    return ee.Image.cat([a1.hypot(b1).rename(b + '_amp'),
                         b1.atan2(a1).rename(b + '_phase')]);
  });
  return ee.Image.cat(imgs).updateMask(seg.gap.lte(1));
}

// Extract harmonics for critical bands
//(2)
var multiHarmonics2025 = getAllHarmonics(landsat_ccdc, 2025.65, ['NIR', 'SWIR1', 'RED']); 
/**
 * Helper function: Generates synthetic image at target time and computes dynamic indices
 */
function getCCDCSeasonalIndices(ccdcAsset, fractionalYear, seasonSuffix, indexList) {
  //Generate CCDC synthetic Landsat stack for specific time t
  var synthImg = ccdc_module.getSyntheticLandsatStack(ccdcAsset, fractionalYear);
  //Standardize Landsat CCDC band names to match spectral library requirements
  var imgForIndices = synthImg.select(
    ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
    ['B2',   'B3',    'B4',  'B5',  'B6',    'B7']
  );
  //define parameter for index calculation
  var param = {
    "B": imgForIndices.select("B2"),
    "G": imgForIndices.select("B3"),
    "R": imgForIndices.select("B4"),
    "N": imgForIndices.select("B5"),
    "S1": imgForIndices.select("B6"),
    "S2": imgForIndices.select("B7"),
    "C1": 6.0,
    "C2": 7.5,
    "g": 2.5,
    "L": 1.0
  };
  //execute the computation
  var computedIndices = spectral.computeIndex(imgForIndices, indexList, param);
  //append seasonal suffix
  var renamedBands = computedIndices.bandNames().map(function(b) {
    return ee.String(b).cat('_').cat(seasonSuffix);
  });
  return computedIndices.rename(renamedBands);
}

//Execute the function
//define selected index
var myIndexList = ["NDVI", "MSAVI", "EVI", "GBNDVI",           //Vegetation
                  "MNDWI", "AWEIsh", "NDMI",                //Water
                   "DBSI", "MBI", "BLFEI"];              //Build
//define the target year
//Peak Dry Season: ~August (Julian Day 240 / 365 = 0.65)
var dryIndices2025 = getCCDCSeasonalIndices(landsat_ccdc, 2025.65, 'dry', myIndexList);
//Peak Wet Season: ~March (Julian Day 60 / 365 = 0.16)
var wetIndices2025 = getCCDCSeasonalIndices(landsat_ccdc, 2025.16, 'wet', myIndexList);
//Calculate Phenological Dynamics (Delta Layers)
//Quantifies maximum seasonal variation across monsoon cycles
var deltaIndicesList = myIndexList.map(function(idx) {
  var wetBand = wetIndices2025.select(idx + '_wet');
  var dryBand = dryIndices2025.select(idx + '_dry');
  return wetBand.subtract(dryBand).rename(idx + '_delta');
});

var deltaIndices2025 = ee.Image.cat(deltaIndicesList);

//Combine into Final CCDC Dynamic Index Feature Stack
//(3)
var temporal_indices = ee.Image.cat([
  dryIndices2025,
  wetIndices2025,
  deltaIndices2025
]).toFloat();

//Topography
var elev = ee.Image("NASA/NASADEM_HGT/001")
            .select('elevation')
            .clip(aoi)
var terrain_product = ee.Terrain.products(elev).select('slope', 'aspect').toFloat();
//Landform 
var landform = ee.Image('projects/ee-rg2v2/assets/TPI_landform_TL')
                      .rename('landform')
                      .toByte()
                      .clip(aoi)
//Topographic Wetness Index
var twi = ee.Image('projects/ee-rg2/assets/Non_spectral_feature/Top_Wetness_index')
                  .rename('TWI')
                  .toFloat()
//Stack topographical data
//(4)
var topo = elev.int16().addBands(terrain_product).addBands(landform).addBands(twi)
print('Topographical Features', topo)

//Grey level co-occurance matrix (GLCM) / Textural Feature
//(5)
var glcm = landsat_dry.select('nir').multiply(1000).toInt32().glcmTexture({size: 1});
var glcm_filter = glcm.select([
  'nir_contrast', 
  'nir_ent', 
  'nir_asm', 
  'nir_corr'
]);

//Stack all data
var feature_2025 = ee.Image.cat([
  multiHarmonics2025,
  temporal_indices,
  topo,
  glcm_filter
  ]).clip(aoi)

print('Original Feature Stck', feature_2025)

//Perform sampling to determine the optimal features
var roi = ee.FeatureCollection('projects/ee-rg2/assets/Samples/Samples_TL_hierarchy_2025_v1');
var sampled_2025 = feature_2025.sampleRegions({
  collection: roi,
  properties: ['new_id'],
  scale: 30,
  projection: 'EPSG:4326',
  tileScale: 4,
  geometries: false //if for rfecv, set the geometries to false
});
//export for RFECV in collab
Export.table.toDrive({
  collection: sampled_2025,
  description: 'TL_LandsatPredictor_Matrix_All',
  fileFormat: 'CSV'
});

//Selected bands 
var landsat_wet_25 = ['B2_wet', 'B3_wet', 'B4_wet', 'B6_wet', 'B7_wet']
var landsat_dry_2025 = ['B3_dry', 'B4_dry','B5_dry','B6_dry', 'B7_dry']
var index_2025_wet = ['AWEIsh_wet', 'BLFEI_wet', 'DBSI_wet', 'EVI_wet', 'GBNDVI_wet','MNDWI_wet','NDVI_wet']
var index_2025_dry = ['BLFEI_dry',  'DBSI_dry', 'GBNDVI_dry', 'MBI_dry', 'MNDWI_dry', 'NDMI_dry', 'NDVI_dry', 'SWIR1_amp',]
Map.addLayer(feature_2025.select(landsat_wet_25), {}, 'Landsat wet Season')
//export to asset
Export.image.toAsset({
  image:  feature_2025.select(landsat_wet_25),
  description: 'Export_TL_Feature_landsat_wet_2025',
  assetId: 'projects/ee-rg2/assets/Landsat_TL_Wet_Season_2025' ,
  pyramidingPolicy: {'.default': 'mean'},
  scale: 30,
  region: aoi,
  crs: 'EPSG:32751', 
  maxPixels: 1e13
});
