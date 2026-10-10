// 1.Define Area of Interest (AOI) - Timor-Leste
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry();

//define the required module
var spectral = require('users/dmlmont/spectral:spectral');
var ccdc_module = require('users/rg2icraf/Luma:ccdc_landsat');

//setup similar parameters for the data
var T = {
  dry: 0.65,    // ~August
  wet: 0.16,    // ~March
  harm: 0.5,    // 2025 script evaluated harmonics at 2025.5
  glcm: 0.20    // 2025 script built GLCM from NIR at 2025.20
};  // identical fractions in every epoch
var TOL = 1;                      // max gap in years for nearest-segment
var index_list = ["NDVI","MSAVI","EVI","GBNDVI","MNDWI","AWEIsh","NDMI","DBSI","MBI","BLFEI"];

// 2. Load Google Global Landsat CCDC Dataset
var global_ccdc = ee.ImageCollection('GOOGLE/GLOBAL_CCDC/V1')
  .filterBounds(aoi)
  .mosaic()
  .clip(aoi)
var ccdc2025 = ee.Image('projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_TL_2020_2025_Landsat').clip(aoi);  
//function to retrieve image harmonics
function getAllHarmonics(ccdcImg, t, bands) {          // revised version from before
  var seg = ccdc_module.nearestSegment(ccdcImg, t);
  var imgs = bands.map(function(b) {
    var c = ccdcImg.select(b + '_coefs')
      .arraySlice(0, seg.idx, seg.idx.add(1)).arrayProject([1]);
    var a1 = c.arrayGet([2]).divide(10000);
    var b1 = c.arrayGet([3]).divide(10000);
    return ee.Image.cat([a1.hypot(b1).rename(b + '_amp'),
                         b1.atan2(a1).rename(b + '_phase')]);
  });
  return ee.Image.cat(imgs).updateMask(seg.gap.lte(TOL));
}

//function to retrieve spectral index
function getCCDCSeasonalIndices(ccdc, t, suffix, tol, coefScale) {
  var synth = ccdc_module.getSyntheticLandsatStack(ccdc, t, false, tol, coefScale);
  var imgForIndices = synth.select(
    ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
    ['B2',   'B3',    'B4',  'B5',  'B6',    'B7']
  );
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
  var computed = spectral.computeIndex(imgForIndices, index_list, param);
  var renamed = computed.bandNames().map(function(b) {
    return ee.String(b).cat('_').cat(suffix);
  });
  return computed.rename(renamed);
}


function buildFeatures(ccdc, year, tol, coefScale) {
  var dryIdx = getCCDCSeasonalIndices(ccdc, year + T.dry, 'dry', tol, coefScale);
  var wetIdx = getCCDCSeasonalIndices(ccdc, year + T.wet, 'wet', tol, coefScale);
  var delta = ee.Image.cat(index_list.map(function(idx) {
    return wetIdx.select(idx + '_wet')
      .subtract(dryIdx.select(idx + '_dry')).rename(idx + '_delta');
  }));
  var harm = ccdc_module.getHarmonics(ccdc, year + T.harm, ['NIR', 'SWIR1', 'RED'], tol, coefScale);
  return ee.Image.cat([
    harm,
    dryIdx.toFloat(), wetIdx.toFloat(), delta.toFloat(),
  ]).clip(aoi);
}
 
function gapLayer(ccdc, year) {
  return ee.Image.cat([
    ccdc_module.nearestSegment(ccdc, year + T.dry).gap.rename('gap_dry'),
    ccdc_module.nearestSegment(ccdc, year + T.wet).gap.rename('gap_wet')
  ]).clip(aoi);
}

//Execute the function and build feature for 2015 and 2005
var f2015 = buildFeatures(global_ccdc, 2015, TOL, 1);
var f2005 = buildFeatures(global_ccdc, 2005, TOL, 1);
print('2015 features', f2015);
print('2005 features', f2005);

//select the optimal features based on 2025 RFECV
var landsat_wet = ['B2_wet', 'B3_wet', 'B4_wet', 'B6_wet', 'B7_wet']
var landsat_dry = ['B3_dry', 'B4_dry','B5_dry','B6_dry', 'B7_dry']
var index_wet = ['AWEIsh_wet', 'BLFEI_wet', 'DBSI_wet', 'EVI_wet', 'GBNDVI_wet','MNDWI_wet','NDVI_wet']
var index_dry = ['BLFEI_dry',  'DBSI_dry', 'GBNDVI_dry', 'MBI_dry', 'MNDWI_dry', 'NDMI_dry', 'NDVI_dry', 'SWIR1_amp',]

//2015 data
var landsat_wet_15 = f2015.select(landsat_wet)
var landsat_dry_15 = f2015.select(landsat_dry)
var index_wet_15 = f2015.select(index_wet)
var index_dry_15 = f2015.select(index_dry)

//2005 data
var landsat_wet_05 = f2005.select(landsat_wet)
var landsat_dry_05 = f2005.select(landsat_dry)
var index_wet_05 = f2005.select(index_wet)
var index_dry_05 = f2005.select(index_dry)


/*
//check
print('Seasonal wet for 2015', landsat_wet_15)
print('Seasonal dry for 2015', landsat_dry_15)
print('Seasonal wet for 2005', landsat_wet_05)
print('Seasonal dry for 2005', landsat_dry_05)


print('Index wet for 2015', index_wet_15)
print('Index dry for 2015', index_dry_15)
print('Index wet for 2005', index_wet_05)
print('Index dry for 2005', index_dry_05)
*/

// ---------------------------------------------------------------------
// Checks (replace the placeholders, then run)
// ---------------------------------------------------------------------
Map.centerObject(aoi, 9);
var rgbVis = {bands: ['red', 'green', 'blue'], min: 0.02, max: 0.25};
Map.addLayer(ccdc_module.getSyntheticLandsatStack(global_ccdc, 2015 + T.dry, false, TOL, 1),
  rgbVis, '2015 dry synthetic RGB');
Map.addLayer(ccdc_module.getSyntheticLandsatStack(global_ccdc, 2005 + T.dry, false, TOL, 1),
  rgbVis, '2005 dry synthetic RGB');
  
// (a) forest-pixel values: NIR should be ~0.25-0.4, red ~0.02-0.05
var forestPt = geometry;   
[['2015', 2015], ['2005', 2005]].forEach(function(p) {
  print('forest pixel ' + p[0],
    ccdc_module.getSyntheticLandsatStack(global_ccdc, p[1] + T.dry, false, TOL, 1)
      .select(['red', 'nir']).reduceRegion(ee.Reducer.first(), forestPt, 30));
});
 
// (b) valid share per band over land (verify the LSIB country code)
var land = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017')
  .filter(ee.Filter.eq('country_co', 'TT')).geometry();
[['2015', f2015], ['2005', f2005]].forEach(function(p) {
  print('valid share per band ' + p[0],
    p[1].mask().reduceRegion({reducer: ee.Reducer.mean(), geometry: land,
      scale: 30, maxPixels: 1e10, tileScale: 4}));
});
 
// (c) tolerance table for each target date
function tolTable(year, frac) {
  var gap = ccdc_module.nearestSegment(global_ccdc, year + frac).gap.unmask(99);
  var tols = [0, 0.25, 0.5, 1, 2, 3, 5];
  return ee.Dictionary.fromLists(
    tols.map(String),
    tols.map(function(tl) {
      return gap.lte(tl).reduceRegion({reducer: ee.Reducer.mean(), geometry: land,
        scale: 120, maxPixels: 1e10}).get('gap_years');
    }));
}
print('tolerance 2015 dry', tolTable(2015, T.dry));
print('tolerance 2015 wet', tolTable(2015, T.wet));
print('tolerance 2005 dry', tolTable(2005, T.dry));
print('tolerance 2005 wet', tolTable(2005, T.wet));  