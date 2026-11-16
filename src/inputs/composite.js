/**
 * Module: Temporal Composites (Direct Query / Non-CCDC Mode)
 * =========================================================================
 * Provides seasonal, annual, and monthly compositing directly from
 * harmonized optical collections (from src/inputs/landsat.js or src/inputs/s2.js).
 *
 * Used when temporalMode: 'composite' is selected in pipelines, eliminating the
 * need to run or export CCDC when multi-year harmonic modeling is not required.
 * =========================================================================
 */

/**
 * Creates monthly median composites from an optical time series collection.
 * Drops months with zero valid observations and sets timestamp to the 15th of the month.
 *
 * @param {ee.ImageCollection} col       - Harmonized collection
 * @param {string}             startDate - Start date (YYYY-MM-DD)
 * @param {string}             endDate   - End date (YYYY-MM-DD)
 * @returns {ee.ImageCollection} Monthly median collection
 */
exports.toMonthlyMedian = function(col, startDate, endDate) {
  var start = ee.Date(startDate);
  var nMonths = ee.Date(endDate).difference(start, 'month').floor();
  var composites = ee.List.sequence(0, nMonths.subtract(1)).map(function(month) {
    var monthStart = start.advance(month, 'month');
    var monthlyImages = col.filterDate(monthStart, monthStart.advance(1, 'month'));
    return monthlyImages.median()
      .set('system:time_start', monthStart.advance(15, 'day').millis())
      .set('n', monthlyImages.size());
  });

  return ee.ImageCollection.fromImages(composites)
    .filter(ee.Filter.gt('n', 0));
};

/**
 * Generates a seasonal composite (e.g. wet or dry season median/percentile) for a given year.
 *
 * @param {ee.ImageCollection} col       - Harmonized collection (should have reflectanceScale = 1.0)
 * @param {number}             year      - Target year
 * @param {number}             startMonth- Starting calendar month (1-12)
 * @param {number}             endMonth  - Ending calendar month (1-12, wrap-around supported e.g. 12 to 4)
 * @param {string}             [reducer='median'] - Reducer: 'median', 'mean', or percentile number (e.g. 25, 50)
 * @returns {ee.Image} Composite image with same bands as collection
 */
exports.seasonalComposite = function(col, year, startMonth, endMonth, reducer) {
  reducer = reducer || 'median';

  // If startMonth > endMonth (e.g. Dec to Apr), spans year-1 to year
  var filtered;
  if (startMonth > endMonth) {
    var startDate = ee.Date.fromYMD(year - 1, startMonth, 1);
    var endDate   = ee.Date.fromYMD(year, endMonth + 1, 1);
    filtered = col.filterDate(startDate, endDate)
      .filter(ee.Filter.calendarRange(startMonth, endMonth, 'month'));
  } else {
    var startDate = ee.Date.fromYMD(year, startMonth, 1);
    var endDate   = ee.Date.fromYMD(year, endMonth + 1, 1);
    filtered = col.filterDate(startDate, endDate);
  }

  if (reducer === 'mean') {
    return filtered.mean();
  }
  if (typeof reducer === 'number') {
    return filtered.reduce(ee.Reducer.percentile([reducer]));
  }
  return filtered.median();
};
