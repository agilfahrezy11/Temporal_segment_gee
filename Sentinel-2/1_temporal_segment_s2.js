// ---------------------------------------------------------------------------
// 1. Setup & Configuration
// ---------------------------------------------------------------------------
var cfg = require('users/rg2icraf/Luma:config');
var ccdc = require('users/rg2icraf/Luma:ccdc_sentinel2');

// Area of Interest (AOI)
var aoi = cfg.getAOI();
Map.centerObject(aoi, 9);
Map.addLayer(aoi, {color: 'FF0000'}, 'AOI Boundary', false);

// ---------------------------------------------------------------------------
// 2. Prep Sentinel-2 Harmonized Time Series & Monthly Medians
// ---------------------------------------------------------------------------
var raw_ts = ccdc.getS2TimeSeries(aoi, '2019-01-01', '2026-01-01', 0.50, 'cs_cdf');
var median_ts = ccdc.toMonthlyMedian(raw_ts, '2019-01-01', '2026-01-01');

// ---------------------------------------------------------------------------
// 3. Run CCDC Algorithm
// ---------------------------------------------------------------------------
var ccdcResults = ccdc.runCCDC(median_ts, {
  breakpointBands: ['green', 'nir', 'swir1']
});

// Retain essential segment temporal boundaries and harmonic coefficients
var selected_ccdc = ccdcResults.select(['tStart', 'tEnd', 'tBreak', 'changeProb', '.*_coefs']);

// ---------------------------------------------------------------------------
// 4. Export CCDC Image Asset
// ---------------------------------------------------------------------------
// Option A: Seamless AOI Export (Recommended for moderate-sized areas)
cfg.toAsset(selected_ccdc, 'CCDC_S2_2020_2025_20m', {
  region: aoi,
  scale: cfg.SCALE.S2,
  crs: cfg.DEFAULT_CRS,
  pyramidingPolicy: { '.default': 'sample' }
});

// Option B: Regional Tiling (for large jurisdictions or EECU memory management)
// Define regional bounding boxes if sub-tiling is required:
var tileBoxes = {
  west_box:      ee.Geometry.Rectangle([124.0, -9.5, 125.2, -8.8]),
  central_timor: ee.Geometry.Rectangle([125.0, -9.2, 126.2, -8.3]),
  east_box:      ee.Geometry.Rectangle([126.0, -8.9, 127.4, -8.2])
};

// To export a specific sub-tile:
// cfg.toAsset(selected_ccdc.clip(tileBoxes.west_box), 'CCDC_TL_2020_2025_20m_westbox', {
//   region: tileBoxes.west_box,
//   scale: cfg.SCALE.S2,
//   crs: cfg.DEFAULT_CRS,
//   pyramidingPolicy: { '.default': 'sample' }
// });

