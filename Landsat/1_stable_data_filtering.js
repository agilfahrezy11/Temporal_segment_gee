// Stable 2025 training samples from one continuous Landsat CCDC model.
// Update this to the asset exported by the unified 2000-2025 CCDC workflow.
var aoi = ee.FeatureCollection('projects/ee-rg2v2/assets/AOI_Timor_Leste_10km').geometry();
var unifiedCcdc = ee.Image(
  'projects/ee-rg2/assets/Landsat_CCDC_TL/CCDC_TL_2000_2025_Landsat'
).clip(aoi);

var trainingSamples = ee.FeatureCollection(
  'projects/ee-rg2/assets/Samples/Samples_TL_hierarchy_2025_v1'
);

/**
 * Returns 1 where no CCDC breakpoint occurs in the given time window and 0
 * otherwise. This is a spectral-stability flag, not proof that the ecological
 * class has remained unchanged.
 */
function stableInPeriod(ccdcImage, startYear, endYear) {
  var tBreak = ccdcImage.select('tBreak');
  var breaksInWindow = tBreak.arrayMask(
    tBreak.gte(startYear).and(tBreak.lte(endYear))
  );
  var hasBreak = breaksInWindow.arrayLength(0).gt(0);
  return hasBreak.not().rename('is_stable');
}

// No product-boundary comparison is needed: all years come from one CCDC run.
var isStable = stableInPeriod(unifiedCcdc, 2000.0, 2025.99);

// reduceRegions preserves the original training-feature properties and
// geometries while adding the is_stable value sampled from the CCDC mask.
var samplesWithStability = isStable.reduceRegions({
  collection: trainingSamples,
  reducer: ee.Reducer.first().setOutputs(['is_stable']),
  scale: 30,
  tileScale: 4
});

var stableTrainingPoints = samplesWithStability
  .filter(ee.Filter.eq('is_stable', 1));

print('Original 2025 sample count:', trainingSamples.size());
print('Stable 2000-2025 sample count:', stableTrainingPoints.size());
print('Stable samples by class:', stableTrainingPoints.aggregate_histogram('new_id'));
Map.centerObject(aoi, 9);
Map.addLayer(isStable, {min: 0, max: 1, palette: ['d73027', '1a9850']},
  'CCDC spectral stability: 2000-2025');

// Start this task in the Code Editor Tasks tab after verifying the counts.
Export.table.toAsset({
  collection: stableTrainingPoints,
  description: 'Samples_TL_hierarchy_2025_stable_2000_2025_v1',
  assetId: 'projects/ee-rg2/assets/Samples/Samples_TL_hierarchy_2025_stable_2000_2025_v1'
});
