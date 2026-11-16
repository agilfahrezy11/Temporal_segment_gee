/**
 * SCRIPT: 7_export_probabilities_and_segments.js
 * PURPOSE: Optimized export of 2020 & 2025 probability stacks and SNIC segments
 *          from GEE assets to Google Drive for local Python post-processing.
 * 
 * EECU & MEMORY OPTIMIZATIONS:
 * 1. Exports are divided by the 5 regional tile footprints (or full AOI) to prevent
 *    the 2 GB Google Drive file limit and avoid GEE worker memory throttling.
 * 2. Probabilities are exported as uint16 (0 - 10,000), reducing file size by 50%.
 * 3. Segments include ID_STRIDE (3e8) so cluster IDs are globally unique across tiles.
 * 4. Includes Cloud-Optimized GeoTIFF (COG) format options for fast windowed reads in Python.
 */

// =========================================================================
// 1. CONFIGURATION & ASSET PATHS
// =========================================================================
var cfg = require('users/rg2icraf/Luma:config');
var aoi = cfg.getAOI();
var CRS = cfg.DEFAULT_CRS;
var FOLDER_NAME = 'AFCLIM_EE'; // Destination folder in Google Drive

// Probability assets
var PROB_2025_ASSET = 'projects/gee-v1-510606/assets/TL_S2_class_probabilities_2025';
var PROB_2020_ASSET = 'projects/gee-v1-510606/assets/TL_S2_class_probabilities_2020'; 
// Landsat assets
var prob_2025_landsat = ee.Image('projects/gee-v1-510606/assets/LC_prob_2025_Landsat_TL');
var prob_2020_landsat = ee.Image('projects/gee-v1-510606/assets/LC_prob_2020_Landsat_TL');
var prob_2015_landsat = ee.Image('projects/gee-v1-510606/assets/LC_prob_2015_Landsat_TL');
var prob_2005_landsat = ee.Image('projects/gee-v1-510606/assets/LC_prob_2005_Landsat_TL');
// SNIC assets
var landsat_segments = ee.Image('projects/gee-v1-510606/assets/seg_landsat_2025_TL');
var segPieces = cfg.S2_SNIC_TILES.map(function(assetPath) {
  var name = assetPath.split('/').pop().replace('seg_', '');
  return { name: name, asset: ee.Image(assetPath) };
});

// =========================================================================
// 2. PREPARE LAYERS & TYPE CASTING
// =========================================================================

// Function to safely format probabilities into uint16 (0 to 10000)
function prepareProbabilities(imageAssetPath) {
  var img = ee.Image(imageAssetPath);
  return img.toUint16();
}

var probs2025 = prepareProbabilities(PROB_2025_ASSET);
var probs2020 = prepareProbabilities(PROB_2020_ASSET);

// Prepare unified segments with global ID_STRIDE (keeps segment IDs unique across tiles)
var ID_STRIDE = 300000000;
var segmentsMosaic = ee.ImageCollection(segPieces.map(function(item, i) {
  return item.asset.select([0]).add(i * ID_STRIDE).toInt32();
})).mosaic().rename('segment');

// =========================================================================
// 3. EXPORT STRATEGY: TILE-BY-TILE (RECOMMENDED FOR EECU & RAM SAFETY)
// =========================================================================
/*
segPieces.forEach(function(item, i) {
  var tileName = item.name;
  var tileGeom = item.asset.geometry();

  // 1. Export 2025 Probabilities for this tile (20m, uint16)
  cfg.toDrive(probs2025.clip(tileGeom), 'TL_S2_probs_2025_' + tileName + '_20m', {
    folder: FOLDER_NAME,
    region: tileGeom,
    scale: cfg.SCALE.S2,
    crs: CRS
  });

  // 2. Export 2020 Probabilities for this tile (20m, uint16)
  cfg.toDrive(probs2020.clip(tileGeom), 'TL_S2_probs_2020_' + tileName + '_20m', {
    folder: FOLDER_NAME,
    region: tileGeom,
    scale: cfg.SCALE.S2,
    crs: CRS
  });

  // 3. Export SNIC Segments for this tile (10m, int32, with global ID stride)
  var tileSegment = item.asset.select([0]).add(i * ID_STRIDE).toInt32().rename('segment');
  cfg.toDrive(tileSegment, 'TL_SNIC_segments_' + tileName + '_10m', {
    folder: FOLDER_NAME,
    region: tileGeom,
    scale: cfg.SCALE.S2_fine,
    crs: CRS
  });
});
*/

// =========================================================================
// 4. MONOLITHIC AOI-WIDE EXPORTS (TO ASSET & DRIVE)
// =========================================================================
// Pyramiding policies:
// - 'mean' for continuous probability rasters
// - 'sample' for SNIC superpixel ID rasters
// - 'mode' for discrete thematic LULC maps

// Landsat Probabilities & Segments (30m)
cfg.toDrive(prob_2005_landsat.toUint16(), 'TL_Landsat_probs_2005_30m', {
  folder: FOLDER_NAME,
  region: aoi,
  scale: cfg.SCALE.Landsat,
  crs: CRS
});

cfg.toDrive(landsat_segments.toInt32(), 'TL_SNIC_segments_Landsat_30m', {
  folder: FOLDER_NAME,
  region: aoi,
  scale: cfg.SCALE.Landsat,
  crs: CRS
});

// Sentinel-2 Probabilities (20m) & Full Segments Mosaic
cfg.toDrive(probs2020, 'TL_S2_probs_2020_national_20m', {
  folder: FOLDER_NAME,
  region: aoi,
  scale: cfg.SCALE.S2,
  crs: CRS
});

cfg.toDrive(probs2025, 'TL_S2_probs_2025_national_20m', {
  folder: FOLDER_NAME,
  region: aoi,
  scale: cfg.SCALE.S2,
  crs: CRS
});

// Also provide direct toAsset export option with explicit pyramiding policies
cfg.toAsset(probs2025, 'projects/gee-v1-510606/assets/TL_S2_probs_2025_export', {
  region: aoi,
  scale: cfg.SCALE.S2,
  crs: CRS,
  pyramidingPolicy: { '.default': 'mean' }
});

cfg.toAsset(segmentsMosaic, 'projects/gee-v1-510606/assets/TL_S2_SNIC_segments_mosaic', {
  region: aoi,
  scale: cfg.SCALE.S2_fine,
  crs: CRS,
  pyramidingPolicy: { '.default': 'sample' }
});