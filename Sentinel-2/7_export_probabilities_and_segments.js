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
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry();
var CRS = 'EPSG:32751'; // UTM Zone 51S
var FOLDER_NAME = 'AFCLIM_EE'; // Destination folder in Google Drive

// Probability assets
var PROB_2025_ASSET = 'projects/gee-v1-510606/assets/TL_S2_class_probabilities_2025';
var PROB_2020_ASSET = 'projects/gee-v1-510606/assets/TL_S2_class_probabilities_2020'; 
//asset for landsat
var prob_2025_landsat = ee.Image('projects/gee-v1-510606/assets/LC_prob_2025_Landsat_TL')
var prob_2020_landsat = ee.Image('projects/gee-v1-510606/assets/LC_prob_2020_Landsat_TL')
var prob_2015_landsat = ee.Image('projects/gee-v1-510606/assets/LC_prob_2015_Landsat_TL')
var prob_2005_landsat = ee.Image('projects/gee-v1-510606/assets/LC_prob_2005_Landsat_TL')
//SNIC asset
var landsat_segments = ee.Image('projects/gee-v1-510606/assets/seg_landsat_2025_TL')
// Existing 10m SNIC segment assets
var segPieces = [
  { name: 'central_timor', asset: ee.Image('projects/gee-v1-510606/assets/snic/seg_central_timor') },
  { name: 'east_box',      asset: ee.Image('projects/gee-v1-510606/assets/snic/seg_east_box') },
  { name: 'off_dili',      asset: ee.Image('projects/gee-v1-510606/assets/snic/seg_off_dili') },
  { name: 'west_box',      asset: ee.Image('projects/gee-v1-510606/assets/snic/seg_west_box') },
  { name: 'east_timor',    asset: ee.Image('projects/gee-v1-510606/assets/snic/seg_east_timor') }
];

// =========================================================================
// 2. PREPARE LAYERS & TYPE CASTING
// =========================================================================

// Function to safely format probabilities into uint16 (0 to 10000)
function prepareProbabilities(imageAssetPath) {
  var img = ee.Image(imageAssetPath);
  // The asset is already stored as scaled integers (0 to 10,000).
  // Casting to uint16 directly preserves exact values with zero computation overhead.
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
// Exporting by the 5 regional tiles splits the 1.5 GB stack into ~300 MB pieces.
// This guarantees:
// - Zero memory limits / timeouts in GEE.
// - Fast parallel downloads from Drive.
// - 1-to-1 matching bounding boxes between probabilities and segments.


/*
segPieces.forEach(function(item, i) {
  var tileName = item.name;
  var tileGeom = item.asset.geometry();

  // 1. Export 2025 Probabilities for this tile (20m, uint16)
  Export.image.toDrive({
    image: probs2025.clip(tileGeom),
    description: 'TL_S2_probs_2025_' + tileName,
    folder: FOLDER_NAME,
    fileNamePrefix: 'TL_S2_probs_2025_' + tileName + '_20m',
    region: tileGeom,
    scale: 20,
    crs: CRS,
    maxPixels: 1e13,
    formatOptions: {
      cloudOptimized: true
    }
  });

  // 2. Export 2020 Probabilities for this tile (20m, uint16)
  Export.image.toDrive({
    image: probs2020.clip(tileGeom),
    description: 'TL_S2_probs_2020_' + tileName,
    folder: FOLDER_NAME,
    fileNamePrefix: 'TL_S2_probs_2020_' + tileName + '_20m',
    region: tileGeom,
    scale: 20,
    crs: CRS,
    maxPixels: 1e13,
    formatOptions: {
      cloudOptimized: true
    }
  });

  // 3. Export SNIC Segments for this tile (10m, int32, with global ID stride)
  var tileSegment = item.asset.select([0]).add(i * ID_STRIDE).toInt32().rename('segment');
  Export.image.toDrive({
    image: tileSegment,
    description: 'TL_SNIC_segments_' + tileName,
    folder: FOLDER_NAME,
    fileNamePrefix: 'TL_SNIC_segments_' + tileName + '_10m',
    region: tileGeom,
    scale: 10,
    crs: CRS,
    maxPixels: 1e13,
    formatOptions: {
      cloudOptimized: true
    }
  });
});
*/
// =========================================================================
// 4. OPTIONAL: MONOLITHIC COUNTRY-WIDE EXPORT (IF PREFERRED)
Export.image.toDrive({
  image: prob_2005_landsat,
  description: 'TL_Landsat_probs_2005_30m',
  folder: FOLDER_NAME,
  fileNamePrefix: 'TL_Landsat_probs_2005_30m',
  region: aoi,
  scale: 30,
  crs: CRS,
  maxPixels: 1e13,
  fileDimensions: 4096, // Prevents Drive 2GB overflow
  formatOptions: { cloudOptimized: true }
});

Export.image.toDrive({
  image: landsat_segments,
  description: 'TL_SNIC_segments_Landsat_30m',
  folder: FOLDER_NAME,
  fileNamePrefix: 'TL_SNIC_segments_Landsat_30m',
  region: aoi,
  scale: 30,
  crs: CRS,
  maxPixels: 1e13,
  formatOptions: { cloudOptimized: true }
});

Export.image.toDrive({
  image: probs2020,
  description: 'TL_S2_probs_2020_national_20m',
  folder: FOLDER_NAME,
  fileNamePrefix: 'TL_S2_probs_2020_national_20m',
  region: aoi,
  scale: 20,
  crs: CRS,
  maxPixels: 1e13,
  fileDimensions: 4096, // Prevents Drive 2GB overflow
  formatOptions: { cloudOptimized: true }
});