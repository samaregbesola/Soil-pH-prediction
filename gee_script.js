// Define the points dataset
var points = ee.FeatureCollection("projects/ee-samuelaregbesola/assets/LUCAS-biodiversity-updated");

// Load WorldClim Bioclimatic variables (Version 1.4)
var worldclim = ee.Image("WORLDCLIM/V1/BIO");

// Load Copernicus DEM (30m)
var dem = ee.ImageCollection("COPERNICUS/DEM/GLO30")
  .filterBounds(points)
  .select('DEM')
  .mosaic();

// Load Sentinel-2 Surface Reflectance
var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
           .filterBounds(points)
           .filterDate('2024-01-01', '2025-01-01')
           .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 40));          

// Load Terra Net NPP and GPP
var terranet = ee.ImageCollection("MODIS/061/MOD17A3HGF")
  .filterBounds(points)
  .filter(ee.Filter.calendarRange(2020, 2024, 'year')) // Keep only last 5 years
  .select(["Gpp", "Npp"]) // Select GPP and NPP bands
  .mean(); // Compute mean across the filtered years

// Function to mask clouds using the Sentinel-2 SCL band
function maskClouds(image) {
    var scl = image.select('SCL');
    var cloudMask = scl.neq(9).and(scl.neq(8)); // Remove clouds & cirrus
    return image.updateMask(cloudMask);
}

// Function to compute NDVI & EVI
function computeIndices(image) {
    var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
    var evi = image.expression(
        '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
        'NIR': image.select('B8'),
        'RED': image.select('B4'),
        'BLUE': image.select('B2')
    }).rename('EVI');

    return image.addBands([ndvi, evi]);
}

// Process Sentinel-2: Filter, mask clouds, compute NDVI & EVI
var s2Processed = s2.map(maskClouds).map(computeIndices);
 
// Function to extract all data (WorldClim, Sentinel-2) at each point
function extractData(feature) {
  try {
    // Sample WorldClim Bioclimatic data
    var bio = worldclim.sample(feature.geometry(), 30).first();

    // Sample Copernicus DEM
    var elevation = dem.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: feature.geometry(),
      scale: 30,
      maxPixels: 1e2
    });

    // Compute NDVI & EVI statistics for the point
    var ndviMean = s2Processed.select('NDVI')
      .reduce(ee.Reducer.mean())
      .sample(feature.geometry(), 10).first();

    var ndviSTD = s2Processed.select('NDVI')
      .reduce(ee.Reducer.stdDev())
      .sample(feature.geometry(), 10).first();

    var ndviPer = s2Processed.select('NDVI')
      .reduce(ee.Reducer.percentile([25, 50, 75]))
      .sample(feature.geometry(), 10).first();

    var eviMean = s2Processed.select('EVI')
      .reduce(ee.Reducer.mean())
      .sample(feature.geometry(), 10).first();

    var eviSTD = s2Processed.select('EVI')
      .reduce(ee.Reducer.stdDev())
      .sample(feature.geometry(), 10).first();

    var eviPer = s2Processed.select('EVI')
      .reduce(ee.Reducer.percentile([25, 50, 75]))
      .sample(feature.geometry(), 10).first();

     
    // GPP and NPP
    var gpp_npp = terranet.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: feature.geometry(),
      scale: 500,  // MODIS resolution is 500m
      maxPixels: 1e4

    });
 
    // Merge results into feature properties
    var extracted = feature.set({

      // WorldClim Bioclimatic Variables
      'Annual_Mean_Temperature': bio ? bio.get('bio01') : null,
      'Mean_Diurnal_Range': bio ? bio.get('bio02') : null,
      'Isothermality': bio ? bio.get('bio03') : null,
      'Temperature_Seasonality': bio ? bio.get('bio04') : null,
      'Max_Temperature_of_Warmest_Month': bio ? bio.get('bio05') : null,
      'Min_Temperature_of_Coldest_Month': bio ? bio.get('bio06') : null,
      'Temperature_Annual_Range': bio ? bio.get('bio07') : null,
      'Mean_Temperature_of_Wettest_Quarter': bio ? bio.get('bio08') : null,
      'Mean_Temperature_of_Driest_Quarter': bio ? bio.get('bio09') : null,
      'Mean_Temperature_of_Warmest_Quarter': bio ? bio.get('bio10') : null,
      'Mean_Temperature_of_Coldest_Quarter': bio ? bio.get('bio11') : null,
      'Annual_Precipitation': bio ? bio.get('bio12') : null,
      'Precipitation_of_Wettest_Month': bio ? bio.get('bio13') : null,
      'Precipitation_of_Driest_Month': bio ? bio.get('bio14') : null,
      'Precipitation_Seasonality': bio ? bio.get('bio15') : null,
      'Precipitation_of_Wettest_Quarter': bio ? bio.get('bio16') : null,
      'Precipitation_of_Driest_Quarter': bio ? bio.get('bio17') : null,
      'Precipitation_of_Warmest_Quarter': bio ? bio.get('bio18') : null,
      'Precipitation_of_Coldest_Quarter': bio ? bio.get('bio19') : null,
 
      // Copernicus DEM
      'Elevation': elevation ? elevation.get('DEM') : null,

      // GPP and NPP
      'GPP_mean_2020-2024': gpp_npp ? gpp_npp.get('Gpp') : null,
      'NPP_mean_2020-2024': gpp_npp ? gpp_npp.get('Npp') : null,

      // Sentinel-2 NDVI statistics
      'NDVI_Mean': ndviMean ? ndviMean.get('NDVI_mean') : null,
      'NDVI_Std': ndviSTD ? ndviSTD.get('NDVI_stdDev') : null,
      'NDVI_Q1': ndviPer ? ndviPer.get('NDVI_p25') : null,
      'NDVI_Q2': ndviPer ? ndviPer.get('NDVI_p50') : null,
      'NDVI_Q3': ndviPer ? ndviPer.get('NDVI_p75') : null,

      // Sentinel-2 EVI statistics
      'EVI_Mean': eviMean? eviMean.get('EVI_mean') : null,
      'EVI_Std': eviSTD ? eviSTD.get('EVI_stdDev') : null,
      'EVI_Q1': eviPer ? eviPer.get('EVI_p25') : null,
      'EVI_Q2': eviPer ? eviPer.get('EVI_p50') : null,
      'EVI_Q3': eviPer ? eviPer.get('EVI_p75') : null, 
    });

    return extracted;

  } catch (error) {
    print(error);

    return feature;
  }
}

// Get the total number of points first
var pointCount = points.size();

// Extract all points at once
var extracted = ee.FeatureCollection(points.toList(pointCount)).map(extractData);


// Export as CSV
Export.table.toDrive({
  collection: extracted,
  description: "Bioclimatic_WorldClim_Sentinel2_TerraNet_DEM_all",
  fileFormat: "CSV"
});