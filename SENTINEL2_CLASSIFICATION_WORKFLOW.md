# Sentinel-2 Classification Workflow

This guide documents the workflow implemented in
[`Sentinel-2/5_end-to-end_sentinel-2.js`](Sentinel-2/5_end-to-end_sentinel-2.js).
The script uses precomputed CCDC assets, trains a multi-probability Random
Forest with 2025 reference samples, and produces class-probability images for
2020 and 2025. It exports probability bands rather than final hard-class maps.

## 1. Inputs and prerequisites

Run the script in the Google Earth Engine Code Editor. The account needs access
to the configured assets and imported GEE modules:

- Timor-Leste AOI with a 10 km buffer:
  `projects/ee-rg2v2/assets/AOI_Timor_Leste_10km`
- Five precomputed CCDC tiles, mosaicked and clipped to the AOI:
  `CCDC_TL_2020_2025_20m_westbox`,
  `CCDC_TL_2020_2025_20m_east-timor`,
  `CCDC_TL_2020_2025_20m_off_dili`,
  `CCDC_TL_2020_2025_20m_central-timor`, and
  `CCDC_TL_2020_2025_20m_east_box`
- 2025 reference samples:
  `projects/ee-rg2/assets/Samples/New_Samples_TL_2025`, with class labels in
  the `new_id` property.
- Imported modules for CCDC reconstruction, spectral indices, train/test
  splitting, Random Forest classification, optional tuning, and Bayesian
  post-processing.
- Topographic assets for the TPI landform and topographic wetness index (TWI).

The class IDs 1–18 are defined in the script and correspond to the class names
in the [README class table](README.md#land-cover-classes).

## 2. Build yearly predictors

The script creates predictor stacks for **2020** and **2025**. For each year,
it reconstructs wet- and dry-season synthetic optical imagery from the CCDC
mosaic:

- Wet season: fractional year `year + 0.16` (approximately March).
- Dry season: fractional year `year + 0.65` (approximately August).
- Harmonic features: evaluated at `year + 0.50`.

The optical stack includes the six reflective bands (blue, green, red, NIR,
SWIR1, and SWIR2) for both seasons; eight spectral indices
(BLFEI, DBSI, GBNDVI, MBI, MNDWI, MSAVI, NDMI, NDVI) for each season; and
wet-minus-dry deltas for BLFEI, DBSI, MBI, MNDWI, MSAVI, and NDMI. Four GLCM
texture metrics (contrast, entropy, ASM, and IDM) are derived from dry-season
NIR. Harmonic features for NIR, red, and SWIR1 are added from the CCDC module.

Additional predictors are:

- **Sentinel-1 GRD:** annual wet-season (December–April) and dry-season
  (May–November) median VV and VH backscatter.
- **ALOS PALSAR-2 ScanSAR:** annual HH and HV median backscatter, retaining
  dual-polarization acquisitions and converting positive digital numbers to
  gamma-nought dB.
- **Topography:** NASADEM elevation, derived slope and aspect, TPI landform,
  and TWI.

The script assembles all these features into a candidate stack for each year.
Before export and classification, it selects a fixed list of **38 predictors**
(`selected_feature`) derived from prior feature selection. The script itself
does not run RFECV.

## 3. Optional feature-selection export

For external feature selection, the script samples the full 2025 candidate
stack at the reference geometries, retaining `new_id`, and prints the input,
extracted, and per-class sample counts. It creates a Google Drive CSV export
named `AFCLIM_Feature_extraction_S2_2025`. Start this task from the GEE Tasks
tab if a new RFECV run is needed. Review the printed counts and predictor band
names before exporting; masked predictors can result in fewer extracted rows
than input samples.

The resulting feature choices must be reflected in the script's
`selected_feature` list before classification. The optional CSV export is
separate from the stratified train/test sampling used to fit and evaluate the
classifier.

## 4. Sampling, training, and evaluation

The classification sample is drawn from 2025 reference data and the selected
2025 predictor stack. The imported `stratifiedSplit` helper is called with a
20 m scale, a 70% training fraction, and seed 42. Training and holdout sample
counts are printed.

Hyperparameter grid search is disabled by default (`RUN_TUNING = false`). If
enabled, it tests tree counts of 100, 200, 400, and 500 against variables-per-
split values of 6, 9, 12, and 19, then reports the best parameters by accuracy.
The classifier currently uses these configured defaults:

| Parameter | Value |
|---|---:|
| Trees | 400 |
| Variables per split | 6 |
| Minimum leaf population | 1 |
| Seed | 42 |

The Random Forest is trained in multi-probability mode using the 2025 training
samples. The holdout set is evaluated for overall and per-class accuracy and
log loss, with per-class sample counts printed in the GEE Console.

## 5. Inference and post-processing

The trained 2025 model is applied to the selected predictor stacks for both
2025 and 2020, so both years use the classifier trained from 2025 reference
samples. Bayesian spatial smoothing is enabled by default
(`APPLY_BAYESIAN_SMOOTHING = true`) and is applied to each year's probability
bands using the configured per-class smoothness values, a 5-pixel neighborhood
parameter, and a 0.5 weighting parameter.

This script does **not** perform probability argmax conversion, minimum mapping
unit filtering, or temporal transition rules. Those steps, if required, belong
to a later classification/export workflow.

## 6. Outputs and checks

The script creates these Earth Engine image-asset export tasks:

- `projects/ee-rg2/assets/TL_S2_class_probabilities_2025`
- `projects/ee-rg2/assets/TL_S2_class_probabilities_2020`

Each output contains one probability band per class. Exports use a 20 m scale,
`EPSG:32751`, and mean pyramiding for the probability bands. Review the GEE
Console output for predictor band names, sample counts, holdout metrics, and
dual-polarization ScanSAR acquisition counts. Run the pending asset exports
from the Tasks tab and verify the resulting assets before using them downstream.
