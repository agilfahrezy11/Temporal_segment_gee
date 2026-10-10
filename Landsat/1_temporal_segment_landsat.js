var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km')
Map.addLayer(aoi)
var landsat_ccdc = require('users/rg2icraf/Luma:ccdc_landsat')
var google_ccdc = ee.ImageCollection("GOOGLE/GLOBAL_CCDC/V1")

//Perform CCDC for 2020-2025 (not covered by Global CCDC)
var time_series_landsat = landsat_ccdc.getLandsatTimeSeries(aoi, '2020-01-01', '2025-12-31');
//determine breakpoint bands
var tl_ccdc = landsat_ccdc.runLandsatCCDC(
  time_series_landsat,
  {breakpointBands: ['green', 'nir', 'swir1']}
  );
//select only a few of CCDC bands
var selected_ccdc = tl_ccdc.select(['tStart','tEnd','tBreak','changeProb','.*_coefs']);

//export to asset
Export.image.toAsset({
  image: selected_ccdc,
  description: 'CCDC_TL_2020_2025_Landsat',
  assetId: 'projects/ee-rg2v2/assets/Landsat_CCDC_TL/CCDC_TL_2020_2025_Landsat',
  region: aoi,
  scale: 30,
  crs: 'EPSG:32751',
  maxPixels: 1e13,
  pyramidingPolicy: { '.default': 'sample' }
});


//var stack2025 = ccdcModule.getSyntheticLandsatStack(customCCDC, 2025.5, true);
//Load Pre-computed CCDC Asset for 2005 & 2015
var globalCCDC = ee.ImageCollection("GOOGLE/GLOBAL_CCDC/V1").filterBounds(aoi).mosaic();
//Use CCCDC module to generate Synthetic Stacks
var stack2005 = ccdcModule.getSyntheticLandsatStack(globalCCDC, 2005.5, true);
var stack2015 = ccdcModule.getSyntheticLandsatStack(globalCCDC, 2015.5, true);

