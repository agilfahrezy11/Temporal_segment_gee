/**
 * =========================================================================
 * Script: Landsat/4_1_SNIC_segmentations_landsat.js
 * Project: Temporal Segment GEE - Timor-Leste Land Cover Mapping
 * 
 * Description:
 * Generates superpixel segments using Simple Non-Iterative Clustering (SNIC)
 * from the unified 2025 CCDC synthetic dry-season optical composite (Red, NIR, SWIR1, NDVI).
 * 
 * Design & Architecture:
 * - Master 2025 Segmentation: Run once at 30m scale to serve as the spatial 
 *   object framework across all mapping years (2025, 2020, 2015, 2005).
 * - Avoids "boundary churn" / false change slivers caused by re-segmenting across epochs.
 * - In downstream hybrid assembly, the Confusion Index (CI) and CCDC stability mask
 *   automatically gate object influence in areas that underwent historical change.
 * - Exports both:
 *     (A) A seamless national 30m asset (feasible given ~16.6M pixels at 30m), and
 *     (B) 5 tiled assets with globally unique integer offsets (identical to CCDC grid).
 * =========================================================================
 */

// =========================================================================
// 1. AOI & Projection Settings
// =========================================================================

var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry();
var CRS = 'EPSG:32751';
var SCALE = 30;

// Export target asset folder/prefix
var ASSET_FOLDER = 'projects/ee-agilfahrezy60/assets/snic_landsat/';

// Toggles for export
var EXPORT_SEAMLESS_NATIONAL = true; // Export single full-country segment asset
var EXPORT_TILED             = true; // Export 5 tiled assets matching CCDC grid

// =========================================================================
// 2. Load Unified CCDC Mosaic & Generate 2025 Synthetic Dry Composite
// =========================================================================

var ccdc_module = require('users/rg2icraf/Luma:ccdc_landsat');

// Load tiled CCDC products generated from Landsat CCDC pipeline (2000-2025)
var img1 = ee.Image('projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_central');
var img2 = ee.Image('projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_east_timor');
var img3 = ee.Image('projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_eastbox');
var img4 = ee.Image('projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_off_dili');
var img5 = ee.Image('projects/ee-agilfahrezy60/assets/Landsat_CCDC_TL/CCDC_2000_2025_westbox');

var landsat_ccdc = ee.ImageCollection([img1, img2, img3, img4, img5]).mosaic().clip(aoi);

// Reconstruct 2025 dry season synthetic reflectance (peak dry season ~ August / 2025.65)
var TOL = 0.005;
var CCDC_COEF_SCALE = 10000;
var synthDry = ccdc_module.getSyntheticLandsatStack(landsat_ccdc, 2025 + 0.65, false, TOL, CCDC_COEF_SCALE);

// Extract primary structural/vegetation bands
var red   = synthDry.select('red');
var nir   = synthDry.select('nir');
var swir1 = synthDry.select('swir1');
var ndvi  = nir.subtract(red).divide(nir.add(red)).rename('NDVI');

var feats = ee.Image.cat([red, nir, swir1, ndvi])
  .rename(['Red', 'NIR', 'SWIR1', 'NDVI']);

// =========================================================================
// 3. Feature Normalization & SNIC Superpixel Clustering
// =========================================================================

// Normalize bands to [0, 1] range to ensure balanced Euclidean distance in SNIC
// Baseline min/max values typical for Landsat surface reflectance across Timor-Leste
var LO = ee.Image.constant([0.01, 0.03, 0.02, -0.15]);
var HI = ee.Image.constant([0.25, 0.45, 0.40,  0.85]);

var segInput = feats.subtract(LO).divide(HI.subtract(LO)).clamp(0, 1)
  .rename(['Red', 'NIR', 'SWIR1', 'NDVI']);

// SNIC parameters calibrated for Landsat 30m:
// - size: 10 pixels -> 300m seed spacing (~9 ha mean object area, preserving smallholder parcels)
// - compactness: 0.8 -> balances natural landscape boundaries with compact polygon shapes
// - connectivity: 8 -> 8-neighbor connectivity
// - neighborhoodSize: 64 -> prevents boundary seams between computation tiles
var snic = ee.Algorithms.Image.Segmentation.SNIC({
  image: segInput,
  size: 10,
  compactness: 0.8,
  connectivity: 8,
  neighborhoodSize: 64
});

var clusters = snic.select('clusters');
var clustersInt = clusters.toInt32();

// =========================================================================
// 4. Map Visualization & Boundary Inspection
// =========================================================================

Map.centerObject(aoi, 9);

// False color composite (SWIR1, NIR, Red)
Map.addLayer(synthDry, {bands: ['swir1', 'nir', 'red'], min: 0.02, max: 0.35}, '2025 Dry Synthetic False Color', true);

// Compute 1-pixel boundary for visual inspection
var boundary = clustersInt.focalMax(1).neq(clustersInt.focalMin(1)).selfMask();
Map.addLayer(boundary, {palette: ['#ffff00']}, 'SNIC 30m Boundaries (Yellow)', false);

// =========================================================================
// 5. Asset Export Configuration
// =========================================================================

// Option A: Single Seamless National Asset (~16.6M pixels at 30m)
if (EXPORT_SEAMLESS_NATIONAL) {
  Export.image.toAsset({
    image:            clustersInt.rename('segment').clip(aoi),
    description:      'TL_Landsat_segments_2025_national',
    assetId:          ASSET_FOLDER + 'seg_landsat_2025_national',
    region:           aoi,
    crs:              CRS,
    scale:            SCALE,
    maxPixels:        1e13,
    pyramidingPolicy: {'.default': 'sample'}
  });
  print('Registered national seamless export: "TL_Landsat_segments_2025_national"');
}

// Option B: 5 Tiled Exports Matching CCDC Sub-Regions with Unique Integer Offsets
if (EXPORT_TILED) {
  var tileList = [
    { name: 'central',    geom: img1.geometry() },
    { name: 'east_timor', geom: img2.geometry() },
    { name: 'eastbox',    geom: img3.geometry() },
    { name: 'off_dili',   geom: img4.geometry() },
    { name: 'westbox',    geom: img5.geometry() }
  ];

  // Offset stride of 300,000,000 ensures unique cluster IDs across all 5 tiles
  // while staying well within the signed 32-bit integer limit (2,147,483,647)
  var OFFSET_STRIDE = 300000000;

  tileList.forEach(function(item, i) {
    var tileGeom = item.geom;
    var tileName = item.name;

    Export.image.toAsset({
      image:            clustersInt.add(i * OFFSET_STRIDE).toInt32().rename('segment').clip(tileGeom),
      description:      'TL_Landsat_segments_2025_' + tileName,
      assetId:          ASSET_FOLDER + 'seg_landsat_2025_' + tileName,
      region:           tileGeom,
      crs:              CRS,
      scale:            SCALE,
      maxPixels:        1e13,
      pyramidingPolicy: {'.default': 'sample'}
    });
  });
  print('Registered 5 tiled exports matching CCDC grid with unique offset stride.');
}

print('SNIC segmentation script initialized. Check Tasks tab to run exports.');
