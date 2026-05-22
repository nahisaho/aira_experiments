"""
Module 2: 気象データと作物モデル（DSSAT/APSIM）の連携
Weather data integration and crop model simulation
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

FIGURES_DIR = Path(__file__).parent.parent / "figures"
RESULTS_DIR = Path(__file__).parent.parent / "results"
DATA_DIR = Path(__file__).parent.parent / "data"


def generate_weather_data(year=2025):
    """
    Generate synthetic weather data for a Japanese rice growing season.
    Location: Niigata Prefecture (37.9°N, 139.0°E) - major rice region.
    """
    np.random.seed(42)
    
    # Growing season: May 1 - October 31
    dates = pd.date_range(f'{year}-05-01', f'{year}-10-31')
    n_days = len(dates)
    doy = dates.dayofyear.values
    
    # Temperature: sinusoidal + noise (Niigata climate)
    tmax_base = 18 + 12 * np.sin(2 * np.pi * (doy - 100) / 365)
    tmin_base = tmax_base - 8 - np.random.uniform(1, 3, n_days)
    tmax = tmax_base + np.random.normal(0, 2, n_days)
    tmin = tmin_base + np.random.normal(0, 1.5, n_days)
    tavg = (tmax + tmin) / 2
    
    # Precipitation: monsoon pattern (Baiu June-July, typhoon Aug-Sep)
    precip_base = 3 + 5 * np.exp(-0.5 * ((doy - 185) / 20) ** 2)  # Baiu peak
    precip_base += 3 * np.exp(-0.5 * ((doy - 245) / 15) ** 2)  # typhoon
    precip = np.random.exponential(precip_base)
    precip = np.where(np.random.random(n_days) > 0.4, precip, 0)
    
    # Solar radiation (MJ/m²/day)
    srad_base = 12 + 8 * np.sin(2 * np.pi * (doy - 80) / 365)
    cloud_factor = np.where(precip > 0, np.random.uniform(0.3, 0.7, n_days), 
                            np.random.uniform(0.7, 1.0, n_days))
    srad = srad_base * cloud_factor
    
    weather_df = pd.DataFrame({
        'date': dates,
        'doy': doy,
        'tmax': np.round(tmax, 1),
        'tmin': np.round(tmin, 1),
        'tavg': np.round(tavg, 1),
        'precip_mm': np.round(precip, 1),
        'srad_MJ': np.round(srad, 1),
    })
    
    weather_df.to_csv(DATA_DIR / "weather_niigata_2025.csv", index=False)
    return weather_df


def simple_rice_model(weather_df, transplant_doy=155, base_temp=10.0,
                      optimal_temp=28.0, critical_gdd=2200):
    """
    Simplified DSSAT/APSIM-like rice growth model.
    Uses Growing Degree Days (GDD) with temperature response function.
    
    Phenological stages (GDD thresholds):
    - Transplanting: 0
    - Tillering: 300 GDD
    - Panicle Initiation: 900 GDD
    - Heading: 1400 GDD
    - Grain Filling: 1700 GDD
    - Maturity: 2200 GDD
    """
    stages = {
        'transplanting': 0, 'tillering': 300, 'panicle_init': 900,
        'heading': 1400, 'grain_filling': 1700, 'maturity': 2200
    }
    
    df = weather_df[weather_df['doy'] >= transplant_doy].copy()
    
    # GDD accumulation with cardinal temperatures
    def temp_response(tavg):
        if tavg <= base_temp:
            return 0
        elif tavg <= optimal_temp:
            return tavg - base_temp
        elif tavg <= 40:
            return max(0, (40 - tavg) / (40 - optimal_temp) * (optimal_temp - base_temp))
        return 0
    
    gdd_daily = df['tavg'].apply(temp_response).values
    gdd_cum = np.cumsum(gdd_daily)
    
    # Biomass accumulation (simplified radiation use efficiency model)
    # RUE = 1.2 g/MJ for rice (varies by growth stage)
    rue = np.where(gdd_cum < 900, 1.0, np.where(gdd_cum < 1700, 1.3, 0.8))
    
    # Light interception (Beer's Law with dynamic LAI)
    lai = np.minimum(6.0, 0.5 + 5.5 * (gdd_cum / 1400) ** 2 * np.exp(-((gdd_cum - 1400) / 800) ** 2 * 0.5))
    k_ext = 0.6  # extinction coefficient for rice
    fpar = 1 - np.exp(-k_ext * lai)
    
    biomass_daily = rue * df['srad_MJ'].values * fpar * 0.01  # t/ha/day
    biomass_cum = np.cumsum(biomass_daily)
    
    # Water stress factor
    precip_cum = np.cumsum(df['precip_mm'].values)
    et_demand = 0.6 * df['srad_MJ'].values  # simplified ET
    et_cum = np.cumsum(et_demand)
    water_ratio = np.minimum(1.0, precip_cum / (et_cum + 1))
    water_stress = np.where(water_ratio > 0.5, 1.0, water_ratio / 0.5)
    
    # Yield = Harvest Index × Final Biomass × Stress Factor
    harvest_index = np.where(gdd_cum > 1700, 0.45, 0.0)
    yield_estimate = biomass_cum * harvest_index * water_stress
    
    # Determine phenological dates
    stage_dates = {}
    for stage, threshold in stages.items():
        idx = np.where(gdd_cum >= threshold)[0]
        if len(idx) > 0:
            stage_dates[stage] = df.iloc[idx[0]]['date']
    
    results = pd.DataFrame({
        'date': df['date'].values,
        'doy': df['doy'].values,
        'gdd_cum': np.round(gdd_cum, 1),
        'lai': np.round(lai, 2),
        'biomass_tha': np.round(biomass_cum, 3),
        'water_stress': np.round(water_stress, 3),
        'yield_tha': np.round(yield_estimate, 3),
    })
    
    results.to_csv(RESULTS_DIR / "crop_model_output.csv", index=False)
    return results, stage_dates


def run_crop_model_analysis():
    """Run weather + crop model pipeline."""
    weather_df = generate_weather_data()
    model_results, stage_dates = simple_rice_model(weather_df)
    
    # --- Figure 3: Weather Data Overview ---
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    
    ax = axes[0]
    ax.fill_between(weather_df['date'], weather_df['tmin'], weather_df['tmax'], alpha=0.3, color='red')
    ax.plot(weather_df['date'], weather_df['tavg'], color='red', linewidth=1.5, label='Tavg')
    ax.set_ylabel('Temperature (°C)')
    ax.set_title('Weather Data — Niigata Prefecture, 2025 Growing Season', fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    ax = axes[1]
    ax.bar(weather_df['date'], weather_df['precip_mm'], color='steelblue', alpha=0.7, width=1)
    ax.set_ylabel('Precipitation (mm)')
    ax.grid(True, alpha=0.3)
    
    ax = axes[2]
    ax.plot(weather_df['date'], weather_df['srad_MJ'], color='orange', linewidth=1, alpha=0.7)
    ax.fill_between(weather_df['date'], 0, weather_df['srad_MJ'], color='orange', alpha=0.2)
    ax.set_ylabel('Solar Radiation (MJ/m²/day)')
    ax.set_xlabel('Date')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig03_weather_data.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # --- Figure 4: Crop Model Results ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    ax = axes[0, 0]
    ax.plot(model_results['date'], model_results['gdd_cum'], color='darkred', linewidth=2)
    for stage, date in stage_dates.items():
        gdd_val = model_results.loc[model_results['date'] == date, 'gdd_cum'].values[0]
        ax.axhline(y=gdd_val, color='gray', linestyle='--', alpha=0.5)
        ax.annotate(stage.replace('_', ' ').title(), xy=(date, gdd_val),
                    fontsize=8, ha='left', va='bottom')
    ax.set_ylabel('Cumulative GDD (°C·day)')
    ax.set_title('Growing Degree Days Accumulation')
    ax.grid(True, alpha=0.3)
    
    ax = axes[0, 1]
    ax.plot(model_results['date'], model_results['lai'], color='green', linewidth=2)
    ax.set_ylabel('LAI (m²/m²)')
    ax.set_title('Leaf Area Index')
    ax.grid(True, alpha=0.3)
    
    ax = axes[1, 0]
    ax.plot(model_results['date'], model_results['biomass_tha'], color='brown', linewidth=2, label='Biomass')
    ax.plot(model_results['date'], model_results['yield_tha'], color='goldenrod', linewidth=2, label='Grain Yield')
    ax.set_ylabel('Biomass / Yield (t/ha)')
    ax.set_title('Biomass and Yield Accumulation')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    ax = axes[1, 1]
    ax.plot(model_results['date'], model_results['water_stress'], color='blue', linewidth=2)
    ax.set_ylabel('Water Stress Factor')
    ax.set_title('Water Stress Index (1.0 = No Stress)')
    ax.set_ylim(0, 1.1)
    ax.grid(True, alpha=0.3)
    
    plt.suptitle('DSSAT/APSIM-like Rice Crop Model Output', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig04_crop_model_results.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    final_yield = model_results['yield_tha'].iloc[-1]
    max_lai = model_results['lai'].max()
    final_biomass = model_results['biomass_tha'].iloc[-1]
    
    print(f"=== Crop Model Results ===")
    print(f"Estimated grain yield: {final_yield:.2f} t/ha")
    print(f"Total aboveground biomass: {final_biomass:.2f} t/ha")
    print(f"Peak LAI: {max_lai:.2f}")
    print(f"Phenological stages:")
    for stage, date in stage_dates.items():
        print(f"  {stage}: {date.strftime('%Y-%m-%d')}")
    
    return weather_df, model_results, stage_dates


if __name__ == "__main__":
    run_crop_model_analysis()
