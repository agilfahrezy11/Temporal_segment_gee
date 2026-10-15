
// ---------- AOI and settings ----------
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry()

var CRS = 'EPSG:32751';

// ---------- 2025 composite (Cloud Score+ masked) ----------
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate('2025-01-01', '2026-01-01')
  .linkCollection(csPlus, ['cs_cdf'])
  .map(function(img) {
    return img.resample('bilinear')               // so 20 m B11 is bilinear at 10 m
      .updateMask(img.select('cs_cdf').gte(0.6))
      .select(['B4', 'B8', 'B11'])
      .divide(10000);
  });
var comp = s2.median().clip(aoi);

var red = comp.select('B4'), nir = comp.select('B8'), swir = comp.select('B11');
var ndvi = nir.subtract(red).divide(nir.add(red));
var feats = ee.Image.cat([red, nir, swir, ndvi])
  .rename(['Red', 'NIR', 'SWIR1', 'NDVI']);

// ---------- stretch to 0-1 (use the lo/hi printed by your Python run) ----------
var LO = ee.Image.constant([0.00756, 0.00595, 0.00648, -0.1725]);
var HI = ee.Image.constant([0.1507, 0.3646, 0.3512, 0.8709]);
var segInput = feats.subtract(LO).divide(HI.subtract(LO)).clamp(0, 1)
  .rename(['Red', 'NIR', 'SWIR1', 'NDVI']);

// ---------- SNIC ----------
var snic = ee.Algorithms.Image.Segmentation.SNIC({
  image: segInput,
  size: 30,               // seed spacing in px -> ~25 px per object (6 -> ~36)
  compactness: 1,      // depends on input scale; test 0.1 / 0.3 / 1
  connectivity: 8,
  neighborhoodSize: 64   // raise if you see tile-edge seams (more compute)
});
var clusters = snic.select('clusters');


// ---------- test view: small AOI only ----------
var testAOI = geometry;  // draw a small rectangle (a few km) over a mixed area
Map.centerObject(testAOI, 13);
Map.addLayer(comp, {bands: ['B8', 'B11', 'B4'], min: 0.02, max: 0.35}, 'False colour');

var cl10 = clusters.reproject({crs: CRS, scale: 10}).clip(testAOI);
var boundary = cl10.focalMax(1).neq(cl10.focalMin(1)).selfMask();
//Map.addLayer(boundary, {palette: ['yellow']}, 'SNIC boundaries');


var clustersInt = clusters.toInt32();   // 'clusters' from your SNIC step, using your tested parameters

//export loop
var tileList = [
  { name: 'east_box', geom: east_box },
  { name: 'east_timor', geom: east_timor },
  { name: 'west_box', geom: west_box },
  { name: 'off_dili', geom: off_dili },
  { name: 'central_timor', geom: central_timor }
];

var OFFSET = 1e6;   // set above the max ID found to keep segment IDs unique per tile

tileList.forEach(function(item, i) {
  var tile = item.geom;
  var tileName = item.name;

  Export.image.toAsset({
    image: clustersInt.add(i * OFFSET).toInt32().rename('segment').clip(tile),
    description: 'TL_segments_' + tileName,
    assetId: 'projects/gee-v1-510606/assets/snic/seg_' + tileName, // Make sure the parent folder exists
    region: tile,
    crs: CRS,
    scale: 10,
    maxPixels: 1e13,
    pyramidingPolicy: {'.default': 'sample'}
  });
});
