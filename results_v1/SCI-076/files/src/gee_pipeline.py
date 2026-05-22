"""
Module 6: GEE/GeoPandas解析パイプライン設計
Google Earth Engine & GeoPandas-based analysis pipeline architecture
"""
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, Polygon, box
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import json

FIGURES_DIR = Path(__file__).parent.parent / "figures"
RESULTS_DIR = Path(__file__).parent.parent / "results"
DATA_DIR = Path(__file__).parent.parent / "data"


# ============================================================
# GEE Pipeline Design (code templates — requires GEE auth)
# ============================================================

GEE_PIPELINE_CODE = """
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
"""


def create_geopandas_field_data():
    """Create GeoPandas-based field management data."""
    np.random.seed(42)
    
    # Create field polygons (rice paddies in Niigata)
    n_fields = 20
    base_lon, base_lat = 139.0, 37.9
    
    fields = []
    field_data = []
    
    for i in range(n_fields):
        # Random field location and size
        lon = base_lon + np.random.uniform(-0.1, 0.1)
        lat = base_lat + np.random.uniform(-0.05, 0.05)
        w = np.random.uniform(0.001, 0.003)
        h = np.random.uniform(0.001, 0.002)
        
        poly = box(lon, lat, lon + w, lat + h)
        fields.append(poly)
        
        area_ha = w * h * 111000 * 111000 * np.cos(np.radians(lat)) / 10000
        
        field_data.append({
            'field_id': f'F{i+1:03d}',
            'variety': np.random.choice(['Koshihikari', 'Akitakomachi', 'Hitomebore'], p=[0.5, 0.3, 0.2]),
            'area_ha': round(area_ha, 2),
            'transplant_date': pd.Timestamp(f'2025-05-{np.random.randint(15, 31):02d}'),
            'ndvi_peak': round(np.random.uniform(0.65, 0.85), 3),
            'yield_tha': round(np.random.uniform(4.5, 7.0), 2),
            'n_applied_kgha': round(np.random.uniform(60, 120), 0),
            'soil_ph': round(np.random.uniform(5.2, 6.5), 1),
        })
    
    gdf = gpd.GeoDataFrame(field_data, geometry=fields, crs='EPSG:4326')
    gdf.to_file(DATA_DIR / "rice_fields.geojson", driver='GeoJSON')
    
    return gdf


def generate_pipeline_architecture():
    """Generate pipeline architecture diagram."""
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 10)
    ax.axis('off')
    
    # Pipeline boxes
    boxes = [
        # Data Sources (top row)
        (1, 8.5, 2.5, 1, '#3498db', 'Sentinel-2 MSI\n(10m, 5 bands)'),
        (4, 8.5, 2.5, 1, '#2980b9', 'Sentinel-1 SAR\n(VV/VH)'),
        (7, 8.5, 2.5, 1, '#e74c3c', 'Weather Station\n(T, P, R)'),
        (10, 8.5, 2.5, 1, '#e67e22', 'Soil Sensors\n(VWC, EC, pH)'),
        (13, 8.5, 2.5, 1, '#27ae60', 'Drone UAV\n(Multispectral)'),
        
        # Processing (middle row)
        (1, 6, 3, 1.2, '#9b59b6', 'GEE Processing\nCloud Masking\nVI Calculation'),
        (5, 6, 3, 1.2, '#1abc9c', 'Crop Model\nDSSAT/APSIM\nGDD Tracking'),
        (9, 6, 3, 1.2, '#f39c12', 'Spatial Analysis\nKriging\nInterpolation'),
        (13, 6, 2.5, 1.2, '#16a085', 'GeoPandas\nField Mgmt'),
        
        # ML Pipeline (lower-middle)
        (4, 3.5, 8, 1.2, '#8e44ad', 'CNN-LSTM Deep Learning Pipeline\nSpatial Feature Extraction → Temporal Sequence Modeling → Yield Prediction'),
        
        # Output (bottom)
        (1, 1, 3.5, 1.5, '#2ecc71', 'Yield Map\n(10m resolution)'),
        (5.5, 1, 3.5, 1.5, '#e74c3c', 'VRA Prescription\n(N rate optimization)'),
        (10, 1, 3.5, 1.5, '#3498db', 'Decision Support\nDashboard'),
    ]
    
    for x, y, w, h, color, text in boxes:
        rect = plt.Rectangle((x, y), w, h, linewidth=2, edgecolor='black',
                              facecolor=color, alpha=0.8, zorder=2)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center',
                fontsize=8, fontweight='bold', color='white', zorder=3)
    
    # Arrows
    arrow_style = dict(arrowstyle='->', color='black', linewidth=1.5)
    connections = [
        ((2.25, 8.5), (2.5, 7.2)),
        ((5.25, 8.5), (6.5, 7.2)),
        ((8.25, 8.5), (6.5, 7.2)),
        ((11.25, 8.5), (10.5, 7.2)),
        ((14.25, 8.5), (14.25, 7.2)),
        ((2.5, 6.0), (8, 4.7)),
        ((6.5, 6.0), (8, 4.7)),
        ((10.5, 6.0), (8, 4.7)),
        ((14.25, 6.0), (12, 4.7)),
        ((6, 3.5), (2.75, 2.5)),
        ((8, 3.5), (7.25, 2.5)),
        ((10, 3.5), (11.75, 2.5)),
    ]
    
    for start, end in connections:
        ax.annotate('', xy=end, xytext=start,
                    arrowprops=dict(arrowstyle='->', color='#2c3e50', lw=1.5))
    
    ax.set_title('Multimodal Crop Yield Prediction System — Architecture Overview',
                 fontsize=15, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig11_pipeline_architecture.png", dpi=300, bbox_inches='tight')
    plt.close()


def run_geopandas_analysis():
    """Run GeoPandas field analysis and generate figures."""
    gdf = create_geopandas_field_data()
    generate_pipeline_architecture()
    
    # Save GEE pipeline code
    with open(RESULTS_DIR / "gee_pipeline_template.py", 'w') as f:
        f.write(GEE_PIPELINE_CODE)
    
    # --- Figure 12: Field Map and Statistics ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    ax = axes[0]
    gdf.plot(column='yield_tha', cmap='YlGn', legend=True, ax=ax, edgecolor='black', linewidth=0.5,
             legend_kwds={'label': 'Yield (t/ha)', 'shrink': 0.7})
    ax.set_title('Rice Paddy Fields — Yield Distribution', fontsize=11)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    
    ax = axes[1]
    gdf.plot(column='n_applied_kgha', cmap='RdYlGn_r', legend=True, ax=ax, edgecolor='black', linewidth=0.5,
             legend_kwds={'label': 'N Applied (kg/ha)', 'shrink': 0.7})
    ax.set_title('Nitrogen Application Rates', fontsize=11)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    
    ax = axes[2]
    variety_counts = gdf['variety'].value_counts()
    colors = ['#2ecc71', '#3498db', '#e74c3c']
    ax.pie(variety_counts.values, labels=variety_counts.index, autopct='%1.0f%%',
           colors=colors, startangle=90)
    ax.set_title('Rice Variety Distribution', fontsize=11)
    
    plt.suptitle('GeoPandas Field Management Analysis — Niigata Prefecture',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig12_field_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Summary
    summary = {
        'n_fields': len(gdf),
        'total_area_ha': round(gdf['area_ha'].sum(), 2),
        'mean_yield': round(gdf['yield_tha'].mean(), 2),
        'varieties': gdf['variety'].value_counts().to_dict(),
        'mean_n_rate': round(gdf['n_applied_kgha'].mean(), 1),
    }
    
    with open(RESULTS_DIR / "field_summary.json", 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    print(f"\n=== GeoPandas Field Analysis ===")
    print(f"Number of fields: {summary['n_fields']}")
    print(f"Total area: {summary['total_area_ha']:.2f} ha")
    print(f"Mean yield: {summary['mean_yield']:.2f} t/ha")
    print(f"Varieties: {summary['varieties']}")
    
    return gdf, summary


if __name__ == "__main__":
    run_geopandas_analysis()
