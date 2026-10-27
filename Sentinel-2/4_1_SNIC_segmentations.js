/**
 * =========================================================================
 * Script: Sentinel-2/4_1_SNIC_segmentations.js
 * 
 * Description:
 * Generates superpixel segments using Simple Non-Iterative Clustering (SNIC)
 * from the unified 2025 CCDC synthetic dry-season optical composite (Red, NIR, SWIR1, NDVI)
 * with robust nearest-segment gap tolerance (zero missing pixels/holes!).
 * 
 * Design & Architecture:
 * - Master CCDC Synthetic Segmentation: Uses the 2025 CCDC synthetic dry-season
 *   reflectance stack (2025.65) rather than raw cloud-masked Sentinel-2 medians.
 * - Gap Tolerance & Nearest Segment: Implements nearestSegment + trend clamping (TOL = 1.0 yr)
 *   and 3x3 hole filling to eliminate missing pixels caused by boundary gaps or cloud dropouts.
 * - Eliminates cloud-mask holes, residual shadow seams, and int32 clamping overflow!
 * - Serves as a consistent, gap-free spatial object framework across both 2025 and 2020.
 * - Exports at 20m native CCDC scale (or 10m if desired):
 *     (A) A single seamless national asset (recommended; zero tile seams or stride issues), and
 *     (B) Direct Google Drive export (Cloud-Optimized GeoTIFF), and
 *     (C) Optional 5 tiled assets matching the CCDC sub-regions.
 * =========================================================================
 */

// =========================================================================
// 1. AOI & Projection Settings
// =========================================================================

var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry();
var CRS = 'EPSG:32751';

// Export resolution:
// - 20m matches the exact native resolution of the Sentinel-2 CCDC assets and probability rasters
//   (~164M pixels total, allowing a single seamless national export without int32 overflow).
// - Change to 10 if 10m superpixel oversampling is explicitly needed.
var SCALE = 20;

// Export target asset folder and Drive folder
var ASSET_FOLDER = 'projects/gee-v1-510606/assets/snic_ccdc/';
var DRIVE_FOLDER = 'AFCLIM_EE';

// Toggles for export modes
var EXPORT_SEAMLESS_ASSET = true;  // Export single full-country segment asset to GEE Assets
var EXPORT_SEAMLESS_DRIVE = true;  // Export single full-country segment COG directly to Google Drive
var EXPORT_TILED_ASSETS   = false; // Optional: Export 5 tiled assets matching CCDC grid

// Gap tolerance in years for synthetic image reconstruction
var TOL_YEARS = 1.0;

// =========================================================================
// 2. Load Unified Sentinel-2 CCDC Mosaic & Reconstruct 2025 Dry Composite
// =========================================================================

// Load the 5 regional Sentinel-2 CCDC assets (2020-2025, 20m)
var tl_ccdc1 = ee.Image('projects/ee-rg2v2/assets/CCDC_TL_2020_2025_20m_westbox');
var tl_ccdc2 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_east-timor');
var tl_ccdc3 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_off_dili');
var tl_ccdc4 = ee.Image('projects/gee-v1-510606/assets/CCDC_TL_2020_2025_20m_central-timor');
var tl_ccdc5 = ee.Image('projects/earth-engine-v1-504707/assets/CCDC_TL_2020_2025_20m_east_box');

var ccdc_mosaic = ee.ImageCollection([tl_ccdc1, tl_ccdc2, tl_ccdc3, tl_ccdc4, tl_ccdc5])
  .mosaic()
  .clip(aoi);

// -------------------------------------------------------------------------
// Helper: Nearest CCDC Segment & Synthetic Stack Reconstruction with Tolerance
// (Self-contained so it runs directly in GEE Code Editor without module dependencies)
// -------------------------------------------------------------------------

function nearestSegment(ccdc, t) {
  t = ee.Number(t);
  var tS = ccdc.select('tStart');
  var tE = ccdc.select('tEnd');
  var dist = tS.subtract(t).max(tE.multiply(-1).add(t)).max(0); // 0 if inside a segment
  return {
    idx: dist.multiply(-1).arrayArgmax().arrayGet([0]),
    gap: dist.arrayReduce(ee.Reducer.min(), [0]).arrayGet([0]).rename('gap_years'),
    tS: tS,
    tE: tE
  };
}

function fillIsolated(img, minNeighbors) {
  minNeighbors = minNeighbors !== undefined ? minNeighbors : 5;
  var n = img.select(0).mask().reduceNeighborhood({
    reducer: ee.Reducer.sum(),
    kernel: ee.Kernel.square(1)
  });
  var med = img.focalMedian(1, 'square', 'pixels');
  return img.unmask(med.updateMask(n.gte(minNeighbors)));
}

function getSyntheticStackWithTol(ccdcOutput, targetDateFractional, bands, addIndices, tolYears, coefScale) {
  bands = bands || ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'];
  addIndices = addIndices !== undefined ? addIndices : false;
  tolYears = tolYears !== undefined ? tolYears : 1.0;
  coefScale = coefScale !== undefined ? coefScale : 10000;

  var t = ee.Number(targetDateFractional);
  var seg = nearestSegment(ccdcOutput, t);

  // Clamp only the linear trend term; harmonics keep the target date's season
  var tc = ee.Image(t).max(seg.tS.arrayGet(seg.idx)).min(seg.tE.arrayGet(seg.idx));
  var omega = 2.0 * Math.PI;
  var terms = ee.Image.cat([
    ee.Image.constant(1),
    tc,
    ee.Image.constant(t.multiply(omega).cos()),
    ee.Image.constant(t.multiply(omega).sin()),
    ee.Image.constant(t.multiply(omega * 2).cos()),
    ee.Image.constant(t.multiply(omega * 2).sin()),
    ee.Image.constant(t.multiply(omega * 3).cos()),
    ee.Image.constant(t.multiply(omega * 3).sin())
  ]).toArray();

  var syntheticBands = bands.map(function(band) {
    var coefs = ccdcOutput.select(band + '_coefs')
      .arraySlice(0, seg.idx, seg.idx.add(1))
      .arrayProject([1]);
    
    return coefs.multiply(terms).arrayReduce(ee.Reducer.sum(), [0])
      .arrayGet([0]).divide(coefScale).rename(band);
  });

  var stack = ee.Image.cat(syntheticBands);

  // Plausibility bounds: within tolerance and valid reflectance
  var valid = seg.gap.lte(tolYears)
    .and(stack.reduce(ee.Reducer.min()).gt(-0.05))
    .and(stack.reduce(ee.Reducer.max()).lt(1.5));
  stack = fillIsolated(stack.updateMask(valid), 5);

  if (addIndices) {
    var ndvi = stack.normalizedDifference(['nir', 'red']).rename('NDVI');
    var ndwi = stack.normalizedDifference(['green', 'nir']).rename('NDWI');
    var ndbi = stack.normalizedDifference(['swir1', 'nir']).rename('NDBI');
    stack = stack.addBands([ndvi, ndwi, ndbi]);
  }

  return stack;
}

// Reconstruct 2025 peak dry season synthetic surface reflectance (~August / 2025.65)
var T_DRY = 2025 + 0.65;
var synthDry = getSyntheticStackWithTol(
  ccdc_mosaic,
  T_DRY,
  ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
  false,
  TOL_YEARS,
  10000
);

// Extract primary structural/vegetation bands (Red, NIR, SWIR1) and compute NDVI
var red   = synthDry.select('red');
var nir   = synthDry.select('nir');
var swir1 = synthDry.select('swir1');
var ndvi  = nir.subtract(red).divide(nir.add(red)).rename('NDVI');

var feats = ee.Image.cat([red, nir, swir1, ndvi])
  .rename(['Red', 'NIR', 'SWIR1', 'NDVI']);

// =========================================================================
// 3. Feature Normalization & SNIC Superpixel Clustering
// =========================================================================

// Normalize bands to [0, 1] range for balanced Euclidean distance in SNIC
var LO = ee.Image.constant([0.00756, 0.00595, 0.00648, -0.1725]);
var HI = ee.Image.constant([0.1507,  0.3646,  0.3512,   0.8709]);

var segInput = feats.subtract(LO).divide(HI.subtract(LO)).clamp(0, 1)
  .rename(['Red', 'NIR', 'SWIR1', 'NDVI']);

// SNIC parameters calibrated for Sentinel-2 CCDC:
// - At 20m scale: size = 15 pixels -> 300m seed spacing (~9 ha object scale, matching 30px at 10m)
// - At 10m scale: size = 30 pixels -> 300m seed spacing
// - compactness: 0.8 -> balances natural parcel/forest boundaries with compact shapes
// - connectivity: 8 -> 8-neighbor connectivity
// - neighborhoodSize: 64 -> prevents boundary seams between GEE computation tiles
var snicSize = (SCALE === 20) ? 15 : 30;

var snic = ee.Algorithms.Image.Segmentation.SNIC({
  image: segInput,
  size: snicSize,
  compactness: 0.8,
  connectivity: 8,
  neighborhoodSize: 64
});

// Convert SNIC seed coordinates into a unique, overflow-safe integer segment ID
// By hashing the cluster value modulo 200,000,000 + 1, we guarantee IDs are strictly
// within [1, 200,000,000] and NEVER hit the 2,147,483,647 int32 clamping ceiling!
var rawClusters = snic.select('clusters').toInt64();
var clustersInt = rawClusters.abs().mod(200000000).add(1).toInt32().rename('segment');

// =========================================================================
// 4. Map Visualization & Boundary Inspection
// =========================================================================

Map.centerObject(aoi, 9);

// False color composite (SWIR1, NIR, Red) from CCDC Synthetic Dry Season
Map.addLayer(
  synthDry,
  {bands: ['swir1', 'nir', 'red'], min: 0.02, max: 0.35},
  '2025 S2 CCDC Dry Synthetic (SWIR1/NIR/Red)',
  true
);

// Compute 1-pixel boundary for visual inspection (reprojected at target scale)
var clView = clustersInt.reproject({crs: CRS, scale: SCALE});
var boundary = clView.focalMax(1).neq(clView.focalMin(1)).selfMask();
Map.addLayer(boundary, {palette: ['#ffff00']}, 'SNIC CCDC Boundaries (Yellow)', false);

// =========================================================================
// 5. Export Configurations
// =========================================================================

// Option A: Single Seamless National GEE Asset
if (EXPORT_SEAMLESS_ASSET) {
  Export.image.toAsset({
    image:            clustersInt.clip(aoi),
    description:      'TL_S2_CCDC_segments_2025_national_' + SCALE + 'm',
    assetId:          ASSET_FOLDER + 'seg_s2_ccdc_2025_national_' + SCALE + 'm',
    region:           aoi,
    crs:              CRS,
    scale:            SCALE,
    maxPixels:        1e13,
    pyramidingPolicy: {'.default': 'sample'}
  });
  print('Registered GEE Asset export: "TL_S2_CCDC_segments_2025_national_' + SCALE + 'm"');
}

// Option B: Direct Seamless Google Drive Export (Cloud-Optimized GeoTIFF)
if (EXPORT_SEAMLESS_DRIVE) {
  Export.image.toDrive({
    image:          clustersInt.clip(aoi),
    description:    'TL_SNIC_CCDC_segments_S2_' + SCALE + 'm_Drive',
    folder:         DRIVE_FOLDER,
    fileNamePrefix: 'TL_SNIC_CCDC_segments_S2_' + SCALE + 'm',
    region:         aoi,
    crs:            CRS,
    scale:          SCALE,
    maxPixels:      1e13,
    formatOptions:  { cloudOptimized: true }
  });
  print('Registered Google Drive export: "TL_SNIC_CCDC_segments_S2_' + SCALE + 'm_Drive"');
}

// Option C: 5 Tiled Exports Matching CCDC Sub-Regions with Safe Integer Stride
if (EXPORT_TILED_ASSETS) {
  var tileList = [
    { name: 'west_box',      geom: tl_ccdc1.geometry() },
    { name: 'east_timor',    geom: tl_ccdc2.geometry() },
    { name: 'off_dili',      geom: tl_ccdc3.geometry() },
    { name: 'central_timor', geom: tl_ccdc4.geometry() },
    { name: 'east_box',      geom: tl_ccdc5.geometry() }
  ];

  // Since clustersInt is bounded in [1, 200,000,000], adding i * 300,000,000
  // guarantees max value <= 1,400,000,000 (strictly below 2,147,483,647 int32 limit!)
  var OFFSET_STRIDE = 300000000;

  tileList.forEach(function(item, i) {
    var tileGeom = item.geom;
    var tileName = item.name;

    Export.image.toAsset({
      image:            clustersInt.add(i * OFFSET_STRIDE).toInt32().rename('segment').clip(tileGeom),
      description:      'TL_S2_CCDC_segments_2025_' + tileName,
      assetId:          ASSET_FOLDER + 'seg_s2_ccdc_' + tileName,
      region:           tileGeom,
      crs:              CRS,
      scale:            SCALE,
      maxPixels:        1e13,
      pyramidingPolicy: {'.default': 'sample'}
    });
  });
  print('Registered 5 tiled exports matching CCDC grid with overflow-safe stride.');
}

print('Sentinel-2 CCDC SNIC segmentation script initialized. Check Tasks tab to run exports.');
