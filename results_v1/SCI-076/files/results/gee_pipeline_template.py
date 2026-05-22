
# ====================================================================
# Google Earth Engine Pipeline for Rice Yield Monitoring
# Target: Niigata Prefecture, Japan
# Imagery: Sentinel-2 MSI (10m) + Sentinel-1 SAR
# ====================================================================

import ee
ee.Initialize()

# --- 1. Define Area of Interest ---
niigata_aoi = ee.Geometry.Rectangle([138.8, 37.7, 139.3, 38.1])

# --- 2. Sentinel-2 NDVI Time Series ---
def compute_vegetation_indices(image):
    ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI')
    evi = image.expression(
        '2.5 * ((NIR - RED) / (NIR + 6*RED - 7.5*BLUE + 1))',
        {'NIR': image.select('B8'), 'RED': image.select('B4'),
         'BLUE': image.select('B2')}
    ).rename('EVI')
    ndre = image.normalizedDifference(['B8', 'B5']).rename('NDRE')
    lai = ndvi.expression('0.57 * exp(2.33 * NDVI)', {'NDVI': ndvi}).rename('LAI')
    return image.addBands([ndvi, evi, ndre, lai])

s2_collection = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(niigata_aoi)
    .filterDate('2025-05-01', '2025-10-31')
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    .map(compute_vegetation_indices))

# --- 3. Sentinel-1 SAR for paddy water detection ---
s1_collection = (ee.ImageCollection('COPERNICUS/S1_GRT_C')
    .filterBounds(niigata_aoi)
    .filterDate('2025-05-01', '2025-10-31')
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .select(['VV', 'VH']))

# --- 4. Weather data from ERA5 ---
era5 = (ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
    .filterDate('2025-05-01', '2025-10-31')
    .filterBounds(niigata_aoi)
    .select(['temperature_2m', 'total_precipitation_sum',
             'surface_solar_radiation_downwards_sum']))

# --- 5. Export NDVI time series to Drive ---
def export_ndvi_timeseries():
    ndvi_ts = s2_collection.select('NDVI')
    monthly = []
    for month in range(5, 11):
        start = f'2025-{month:02d}-01'
        end = f'2025-{month+1:02d}-01' if month < 10 else '2025-11-01'
        composite = ndvi_ts.filterDate(start, end).median()
        monthly.append(composite.set('month', month))
    
    for i, img in enumerate(monthly):
        task = ee.batch.Export.image.toDrive(
            image=img.clip(niigata_aoi),
            description=f'NDVI_month_{i+5}',
            scale=10,
            region=niigata_aoi,
            maxPixels=1e9
        )
        task.start()

# --- 6. Crop type classification using random forest ---
def classify_rice_paddies():
    training_points = ee.FeatureCollection('users/your_asset/rice_training_points')
    bands = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'NDVI', 'EVI']
    composite = s2_collection.filterDate('2025-07-01', '2025-08-31').median()
    training = composite.select(bands).sampleRegions(
        collection=training_points,
        properties=['class'],
        scale=10
    )
    classifier = ee.Classifier.smileRandomForest(100).train(
        features=training, classProperty='class', inputProperties=bands)
    classified = composite.select(bands).classify(classifier)
    return classified
