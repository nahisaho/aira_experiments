#!/usr/bin/env python3
"""
Integrated Framework for Economic Valuation of Ecosystem Services
Based on InVEST/ARIES pipeline with SEEA-EA integration
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats, optimize
from scipy.interpolate import griddata
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# 1. Ecosystem Service Classification & Quantitative Indicators
# ============================================================

def create_es_classification():
    """Design quantitative indicators for provisioning, regulating, cultural services."""
    categories = {
        'Provisioning': {
            'Food production': {'unit': 'ton/ha/yr', 'indicator': 'Crop yield'},
            'Timber': {'unit': 'm³/ha/yr', 'indicator': 'Wood volume'},
            'Freshwater': {'unit': 'mm/yr', 'indicator': 'Water yield'},
            'Non-timber forest products': {'unit': 'kg/ha/yr', 'indicator': 'NTFP harvest'}
        },
        'Regulating': {
            'Carbon sequestration': {'unit': 'tC/ha/yr', 'indicator': 'Net carbon flux'},
            'Water purification': {'unit': 'mg N removed/L', 'indicator': 'Nutrient retention'},
            'Pollination': {'unit': '% crop dependence', 'indicator': 'Pollinator abundance'},
            'Flood regulation': {'unit': 'mm retention', 'indicator': 'Water retention capacity'},
            'Erosion control': {'unit': 'ton soil/ha/yr avoided', 'indicator': 'USLE reduction'}
        },
        'Cultural': {
            'Recreation': {'unit': 'visits/yr', 'indicator': 'Visitor days'},
            'Aesthetic value': {'unit': 'scenic index', 'indicator': 'Landscape beauty'},
            'Traditional knowledge': {'unit': 'practices maintained', 'indicator': 'Cultural practices'},
            'Educational value': {'unit': 'programs/yr', 'indicator': 'Environmental education'}
        }
    }
    return categories

# ============================================================
# 2. InVEST-based Spatial Service Quantification (Simulated)
# ============================================================

def simulate_invest_spatial(grid_size=50):
    """Simulate InVEST model outputs for spatial ecosystem service mapping."""
    x = np.linspace(0, 10, grid_size)
    y = np.linspace(0, 10, grid_size)
    X, Y = np.meshgrid(x, y)

    # Land use classes: 0=urban, 1=agriculture, 2=forest, 3=wetland, 4=satoyama
    lulc = np.zeros((grid_size, grid_size), dtype=int)
    # Forest in northwest
    lulc[:20, :20] = 2
    # Wetland near center
    lulc[20:30, 20:30] = 3
    # Satoyama (mixed) in northeast
    lulc[:25, 30:] = 4
    # Agriculture in south
    lulc[30:, :] = 1
    # Urban in southeast corner
    lulc[40:, 35:] = 0

    # Carbon storage (tC/ha) - InVEST Carbon model simulation
    carbon_pool = {0: 5, 1: 30, 2: 180, 3: 120, 4: 95}
    carbon = np.vectorize(lambda v: carbon_pool[v])(lulc).astype(float)
    carbon += np.random.normal(0, 10, carbon.shape)
    carbon = np.clip(carbon, 0, 250)

    # Water yield (mm/yr) - InVEST Water Yield model
    precip = 1500 + 200 * np.sin(Y / 3)
    et_coeff = {0: 0.3, 1: 0.6, 2: 0.7, 3: 0.4, 4: 0.55}
    et = np.vectorize(lambda v: et_coeff[v])(lulc) * precip
    water_yield = precip - et + np.random.normal(0, 30, carbon.shape)
    water_yield = np.clip(water_yield, 0, None)

    # Habitat quality - InVEST Habitat Quality model
    habitat_score = {0: 0.1, 1: 0.3, 2: 0.9, 3: 0.85, 4: 0.7}
    habitat = np.vectorize(lambda v: habitat_score[v])(lulc).astype(float)
    # Distance-based degradation from urban areas
    urban_mask = (lulc == 0).astype(float)
    from scipy.ndimage import distance_transform_edt
    dist_from_urban = distance_transform_edt(1 - urban_mask)
    degradation = np.exp(-dist_from_urban / 15)
    habitat = habitat * (1 - 0.5 * degradation)
    habitat += np.random.normal(0, 0.03, habitat.shape)
    habitat = np.clip(habitat, 0, 1)

    # Sediment retention (ton/ha/yr) - InVEST SDR model
    sed_export = {0: 0.5, 1: 8.0, 2: 0.3, 3: 0.2, 4: 1.5}
    sediment = np.vectorize(lambda v: sed_export[v])(lulc).astype(float)
    slope = np.abs(np.gradient(np.sin(X) * np.cos(Y) * 100)[0])
    sediment = sediment * (1 + slope * 0.5)
    sediment += np.random.normal(0, 0.3, sediment.shape)
    sediment = np.clip(sediment, 0, None)

    # Pollination - InVEST Pollination model
    pollination = {0: 0.05, 1: 0.4, 2: 0.6, 3: 0.3, 4: 0.75}
    poll = np.vectorize(lambda v: pollination[v])(lulc).astype(float)
    poll += np.random.normal(0, 0.05, poll.shape)
    poll = np.clip(poll, 0, 1)

    return {
        'X': X, 'Y': Y, 'lulc': lulc,
        'carbon': carbon, 'water_yield': water_yield,
        'habitat': habitat, 'sediment': sediment, 'pollination': poll,
        'grid_size': grid_size
    }

# ============================================================
# 3. Choice Experiment for WTP Estimation
# ============================================================

def simulate_choice_experiment(n_respondents=500):
    """Design and simulate a choice experiment for WTP estimation."""
    # Attributes and levels
    attributes = {
        'biodiversity': [0, 1, 2, 3],        # species richness improvement level
        'water_quality': [0, 1, 2],           # water quality improvement level
        'landscape_beauty': [0, 1, 2],        # landscape maintenance level
        'recreation_access': [0, 1],          # recreation access improvement
        'cost': [0, 500, 1000, 2000, 5000]    # annual household payment (JPY)
    }

    # True WTP parameters (conditional logit)
    true_beta = {
        'biodiversity': 0.8,
        'water_quality': 0.6,
        'landscape_beauty': 0.5,
        'recreation_access': 0.3,
        'cost': -0.0005
    }

    # Generate choice sets (D-efficient design simulation)
    n_choice_sets = 12
    n_alternatives = 3  # 2 alternatives + status quo

    results = []
    for resp_id in range(n_respondents):
        # Individual heterogeneity
        individual_beta = {}
        for attr, beta in true_beta.items():
            if attr == 'cost':
                individual_beta[attr] = beta * np.exp(np.random.normal(0, 0.3))
            else:
                individual_beta[attr] = beta + np.random.normal(0, 0.15)

        for cs in range(n_choice_sets):
            utilities = []
            alternatives_data = []
            for alt in range(n_alternatives):
                if alt == n_alternatives - 1:  # Status quo
                    attrs = {k: 0 for k in attributes}
                else:
                    attrs = {k: np.random.choice(v) for k, v in attributes.items()}

                utility = sum(individual_beta[k] * attrs[k] for k in attributes)
                utility += np.random.gumbel(0, 1)  # Type I extreme value error
                utilities.append(utility)
                alternatives_data.append(attrs)

            chosen = np.argmax(utilities)
            for alt_idx, (alt_data, u) in enumerate(zip(alternatives_data, utilities)):
                results.append({
                    'respondent_id': resp_id,
                    'choice_set': cs,
                    'alternative': alt_idx,
                    'chosen': 1 if alt_idx == chosen else 0,
                    **alt_data
                })

    df = pd.DataFrame(results)

    # Estimate conditional logit via maximum likelihood (vectorized)
    beta_names = ['biodiversity', 'water_quality', 'landscape_beauty', 'recreation_access', 'cost']
    
    # Pre-compute group indices for speed
    group_keys = df.groupby(['respondent_id', 'choice_set']).ngroup()
    n_groups = group_keys.max() + 1
    chosen_arr = df['chosen'].values
    attr_matrix = df[beta_names].values  # (N, 5)
    group_arr = group_keys.values
    
    def neg_log_likelihood(params):
        V = attr_matrix @ params
        nll = 0.0
        for g in range(n_groups):
            mask = group_arr == g
            v_g = V[mask]
            c_g = chosen_arr[mask]
            v_g_shifted = v_g - v_g.max()
            exp_v = np.exp(v_g_shifted)
            prob = exp_v / exp_v.sum()
            chosen_prob = prob[c_g == 1]
            if len(chosen_prob) > 0 and chosen_prob[0] > 1e-15:
                nll -= np.log(chosen_prob[0])
        return nll

    x0 = np.array([0.5, 0.4, 0.3, 0.2, -0.0003])
    result = optimize.minimize(neg_log_likelihood, x0, method='Nelder-Mead',
                                options={'maxiter': 3000, 'xatol': 1e-5, 'fatol': 1e-5})

    estimated_beta = dict(zip(
        ['biodiversity', 'water_quality', 'landscape_beauty', 'recreation_access', 'cost'],
        result.x
    ))

    # Calculate WTP (marginal WTP = -beta_attr / beta_cost)
    wtp = {}
    for attr in ['biodiversity', 'water_quality', 'landscape_beauty', 'recreation_access']:
        wtp[attr] = -estimated_beta[attr] / estimated_beta['cost']

    return df, estimated_beta, wtp, true_beta

# ============================================================
# 4. Discount Rate & Intergenerational Equity
# ============================================================

def analyze_discount_rates():
    """Analyze impact of different discount rates on ES valuation."""
    years = np.arange(0, 101)
    annual_es_value = 1e6  # Base annual ES value (JPY/ha)

    discount_scenarios = {
        'Market rate (5%)': 0.05,
        'Social rate (3%)': 0.03,
        'Stern Review (1.4%)': 0.014,
        'Ramsey optimal (2.5%)': 0.025,
        'Hyperbolic declining': None,  # Special treatment
        'Zero discount': 0.0
    }

    npv_results = {}
    pv_series = {}
    for name, rate in discount_scenarios.items():
        if name == 'Hyperbolic declining':
            # Gamma discounting (Weitzman 2001)
            rates = np.array([0.04 * np.exp(-0.03 * t) + 0.01 for t in years])
            cumulative_discount = np.exp(-np.cumsum(rates))
            pv = annual_es_value * cumulative_discount
        elif rate == 0.0:
            pv = annual_es_value * np.ones_like(years, dtype=float)
        else:
            pv = annual_es_value / (1 + rate) ** years
        npv_results[name] = np.sum(pv)
        pv_series[name] = pv

    return years, pv_series, npv_results

# ============================================================
# 5. SEEA-EA Natural Capital Accounting
# ============================================================

def create_seea_accounts(invest_data):
    """Create SEEA-EA compatible ecosystem accounts."""
    lulc = invest_data['lulc']
    grid_size = invest_data['grid_size']

    # Ecosystem extent account (ha)
    lulc_names = {0: 'Urban', 1: 'Agriculture', 2: 'Forest', 3: 'Wetland', 4: 'Satoyama'}
    extent = {}
    for code, name in lulc_names.items():
        extent[name] = np.sum(lulc == code) * 4  # each cell = 4 ha

    # Ecosystem condition account
    condition_indicators = {}
    for name, code in [('Forest', 2), ('Wetland', 3), ('Satoyama', 4), ('Agriculture', 1)]:
        mask = lulc == code
        condition_indicators[name] = {
            'Carbon density (tC/ha)': float(np.mean(invest_data['carbon'][mask])),
            'Habitat quality (0-1)': float(np.mean(invest_data['habitat'][mask])),
            'Water yield (mm/yr)': float(np.mean(invest_data['water_yield'][mask])),
            'Pollination index': float(np.mean(invest_data['pollination'][mask]))
        }

    # Monetary supply-use account (JPY/ha/yr)
    unit_values = {
        'carbon': 5000,      # JPY per tC
        'water': 50,         # JPY per mm water yield
        'habitat': 200000,   # JPY per habitat quality unit
        'sediment': 3000,    # JPY per ton sediment avoided
        'pollination': 150000  # JPY per pollination index
    }

    monetary_account = {}
    for name, code in lulc_names.items():
        mask = lulc == code
        if np.sum(mask) == 0:
            continue
        total_ha = np.sum(mask) * 4
        monetary_account[name] = {
            'Extent (ha)': total_ha,
            'Carbon value': float(np.mean(invest_data['carbon'][mask]) * unit_values['carbon'] * total_ha),
            'Water value': float(np.mean(invest_data['water_yield'][mask]) * unit_values['water'] * total_ha),
            'Habitat value': float(np.mean(invest_data['habitat'][mask]) * unit_values['habitat'] * total_ha),
            'Erosion control value': float((10 - np.mean(invest_data['sediment'][mask])) * unit_values['sediment'] * total_ha),
            'Pollination value': float(np.mean(invest_data['pollination'][mask]) * unit_values['pollination'] * total_ha),
        }
        monetary_account[name]['Total ES value (JPY)'] = sum(v for k, v in monetary_account[name].items() if k != 'Extent (ha)')

    return extent, condition_indicators, monetary_account

# ============================================================
# 6. Satoyama Case Study
# ============================================================

def satoyama_case_study(invest_data, wtp_results):
    """Evaluate ecosystem services for Satoyama landscape."""
    lulc = invest_data['lulc']
    satoyama_mask = lulc == 4

    # Biophysical assessment
    satoyama_services = {
        'Carbon storage (tC/ha)': float(np.mean(invest_data['carbon'][satoyama_mask])),
        'Water yield (mm/yr)': float(np.mean(invest_data['water_yield'][satoyama_mask])),
        'Habitat quality': float(np.mean(invest_data['habitat'][satoyama_mask])),
        'Pollination index': float(np.mean(invest_data['pollination'][satoyama_mask])),
        'Sediment export (ton/ha/yr)': float(np.mean(invest_data['sediment'][satoyama_mask]))
    }

    # Comparison with other land uses
    comparison = {}
    lulc_names = {1: 'Agriculture', 2: 'Forest', 3: 'Wetland', 4: 'Satoyama'}
    for code, name in lulc_names.items():
        mask = lulc == code
        comparison[name] = {
            'Carbon (tC/ha)': float(np.mean(invest_data['carbon'][mask])),
            'Water yield (mm/yr)': float(np.mean(invest_data['water_yield'][mask])),
            'Habitat quality': float(np.mean(invest_data['habitat'][mask])),
            'Pollination': float(np.mean(invest_data['pollination'][mask]))
        }

    # Scenario analysis: Satoyama conversion
    scenarios = {
        'Baseline (Satoyama maintained)': 1.0,
        'Partial conversion to agriculture (50%)': 0.5,
        'Full conversion to agriculture': 0.0,
        'Satoyama expansion (+25%)': 1.25,
        'Restoration with management': 1.15
    }

    baseline_value = (
        satoyama_services['Carbon storage (tC/ha)'] * 5000 +
        satoyama_services['Water yield (mm/yr)'] * 50 +
        satoyama_services['Habitat quality'] * 200000 +
        satoyama_services['Pollination index'] * 150000 +
        (10 - satoyama_services['Sediment export (ton/ha/yr)']) * 3000
    )

    agriculture_value = (
        float(np.mean(invest_data['carbon'][lulc == 1])) * 5000 +
        float(np.mean(invest_data['water_yield'][lulc == 1])) * 50 +
        float(np.mean(invest_data['habitat'][lulc == 1])) * 200000 +
        float(np.mean(invest_data['pollination'][lulc == 1])) * 150000 +
        (10 - float(np.mean(invest_data['sediment'][lulc == 1]))) * 3000
    )

    scenario_values = {}
    for name, factor in scenarios.items():
        if factor >= 1.0:
            scenario_values[name] = baseline_value * factor
        else:
            scenario_values[name] = baseline_value * factor + agriculture_value * (1 - factor)

    # Add WTP-based cultural value
    cultural_wtp = wtp_results.get('landscape_beauty', 0) + wtp_results.get('recreation_access', 0)

    return satoyama_services, comparison, scenario_values, cultural_wtp, baseline_value

# ============================================================
# Visualization Functions
# ============================================================

def plot_spatial_maps(invest_data):
    """Generate spatial maps of InVEST outputs."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Land use map
    lulc_cmap = plt.cm.get_cmap('Set3', 5)
    im0 = axes[0, 0].imshow(invest_data['lulc'], cmap=lulc_cmap, origin='lower', extent=[0, 10, 0, 10])
    axes[0, 0].set_title('Land Use / Land Cover', fontsize=13, fontweight='bold')
    cbar0 = plt.colorbar(im0, ax=axes[0, 0], ticks=[0, 1, 2, 3, 4])
    cbar0.ax.set_yticklabels(['Urban', 'Agri', 'Forest', 'Wetland', 'Satoyama'], fontsize=9)

    # Carbon storage
    im1 = axes[0, 1].imshow(invest_data['carbon'], cmap='YlGn', origin='lower', extent=[0, 10, 0, 10])
    axes[0, 1].set_title('Carbon Storage (tC/ha)', fontsize=13, fontweight='bold')
    plt.colorbar(im1, ax=axes[0, 1])

    # Water yield
    im2 = axes[0, 2].imshow(invest_data['water_yield'], cmap='Blues', origin='lower', extent=[0, 10, 0, 10])
    axes[0, 2].set_title('Water Yield (mm/yr)', fontsize=13, fontweight='bold')
    plt.colorbar(im2, ax=axes[0, 2])

    # Habitat quality
    im3 = axes[1, 0].imshow(invest_data['habitat'], cmap='RdYlGn', origin='lower', extent=[0, 10, 0, 10])
    axes[1, 0].set_title('Habitat Quality Index', fontsize=13, fontweight='bold')
    plt.colorbar(im3, ax=axes[1, 0])

    # Sediment export
    im4 = axes[1, 1].imshow(invest_data['sediment'], cmap='OrRd', origin='lower', extent=[0, 10, 0, 10])
    axes[1, 1].set_title('Sediment Export (ton/ha/yr)', fontsize=13, fontweight='bold')
    plt.colorbar(im4, ax=axes[1, 1])

    # Pollination
    im5 = axes[1, 2].imshow(invest_data['pollination'], cmap='PuBuGn', origin='lower', extent=[0, 10, 0, 10])
    axes[1, 2].set_title('Pollination Service Index', fontsize=13, fontweight='bold')
    plt.colorbar(im5, ax=axes[1, 2])

    for ax in axes.flat:
        ax.set_xlabel('Distance (km)')
        ax.set_ylabel('Distance (km)')

    plt.suptitle('InVEST-based Spatial Ecosystem Service Quantification', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('figures/spatial_es_maps.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Saved: figures/spatial_es_maps.png")

def plot_wtp_results(estimated_beta, wtp, true_beta):
    """Plot WTP estimation results."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Beta coefficients comparison
    attrs = ['biodiversity', 'water_quality', 'landscape_beauty', 'recreation_access', 'cost']
    true_vals = [true_beta[a] for a in attrs]
    est_vals = [estimated_beta[a] for a in attrs]

    x_pos = np.arange(len(attrs))
    width = 0.35
    axes[0].bar(x_pos - width/2, true_vals, width, label='True β', color='steelblue', alpha=0.8)
    axes[0].bar(x_pos + width/2, est_vals, width, label='Estimated β', color='coral', alpha=0.8)
    axes[0].set_xticks(x_pos)
    axes[0].set_xticklabels([a.replace('_', '\n') for a in attrs], fontsize=9)
    axes[0].set_ylabel('Coefficient Value')
    axes[0].set_title('Conditional Logit: True vs Estimated β', fontweight='bold')
    axes[0].legend()
    axes[0].axhline(y=0, color='gray', linestyle='--', alpha=0.5)

    # WTP estimates
    wtp_attrs = list(wtp.keys())
    wtp_vals = [wtp[a] for a in wtp_attrs]
    colors = ['#2ecc71', '#3498db', '#e74c3c', '#f39c12']
    bars = axes[1].bar(range(len(wtp_attrs)), wtp_vals, color=colors, alpha=0.85)
    axes[1].set_xticks(range(len(wtp_attrs)))
    axes[1].set_xticklabels([a.replace('_', '\n') for a in wtp_attrs], fontsize=9)
    axes[1].set_ylabel('WTP (JPY/household/year)')
    axes[1].set_title('Marginal Willingness to Pay', fontweight='bold')
    for bar, val in zip(bars, wtp_vals):
        axes[1].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 20,
                     f'¥{val:,.0f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Cost sensitivity
    cost_range = np.linspace(0, 10000, 100)
    for attr in ['biodiversity', 'water_quality', 'landscape_beauty']:
        prob = 1 / (1 + np.exp(-(estimated_beta[attr] * 2 + estimated_beta['cost'] * cost_range)))
        axes[2].plot(cost_range, prob, label=attr.replace('_', ' ').title(), linewidth=2)
    axes[2].set_xlabel('Household Payment (JPY/year)')
    axes[2].set_ylabel('Choice Probability')
    axes[2].set_title('WTP Sensitivity Analysis', fontweight='bold')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('figures/wtp_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Saved: figures/wtp_analysis.png")

def plot_discount_analysis(years, pv_series, npv_results):
    """Plot discount rate analysis."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6', '#1abc9c']
    for (name, pv), color in zip(pv_series.items(), colors):
        axes[0].plot(years, pv / 1e6, label=name, linewidth=2, color=color)
    axes[0].set_xlabel('Year', fontsize=12)
    axes[0].set_ylabel('Present Value (Million JPY/ha)', fontsize=12)
    axes[0].set_title('Discount Rate Impact on ES Present Value', fontsize=13, fontweight='bold')
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)
    axes[0].set_xlim(0, 100)

    # NPV comparison bar chart
    names = list(npv_results.keys())
    values = [npv_results[n] / 1e6 for n in names]
    bars = axes[1].barh(range(len(names)), values, color=colors, alpha=0.85)
    axes[1].set_yticks(range(len(names)))
    axes[1].set_yticklabels(names, fontsize=10)
    axes[1].set_xlabel('100-Year NPV (Million JPY/ha)', fontsize=12)
    axes[1].set_title('Net Present Value by Discount Scenario', fontsize=13, fontweight='bold')
    for bar, val in zip(bars, values):
        axes[1].text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2.,
                     f'{val:.1f}M', ha='left', va='center', fontsize=10)

    plt.tight_layout()
    plt.savefig('figures/discount_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Saved: figures/discount_analysis.png")

def plot_seea_accounts(extent, condition_indicators, monetary_account):
    """Plot SEEA-EA ecosystem accounts."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Extent account
    names = list(extent.keys())
    values = list(extent.values())
    colors = ['#95a5a6', '#f1c40f', '#27ae60', '#3498db', '#e67e22']
    axes[0, 0].pie(values, labels=names, colors=colors, autopct='%1.1f%%', startangle=90, textprops={'fontsize': 11})
    axes[0, 0].set_title('Ecosystem Extent Account (ha)', fontsize=13, fontweight='bold')

    # Condition indicators heatmap
    eco_types = list(condition_indicators.keys())
    indicators = list(condition_indicators[eco_types[0]].keys())
    data = np.array([[condition_indicators[e][i] for i in indicators] for e in eco_types])
    # Normalize each column
    data_norm = (data - data.min(axis=0)) / (data.max(axis=0) - data.min(axis=0) + 1e-10)

    im = axes[0, 1].imshow(data_norm, cmap='RdYlGn', aspect='auto')
    axes[0, 1].set_xticks(range(len(indicators)))
    axes[0, 1].set_xticklabels([i.split('(')[0].strip() for i in indicators], rotation=45, ha='right', fontsize=9)
    axes[0, 1].set_yticks(range(len(eco_types)))
    axes[0, 1].set_yticklabels(eco_types, fontsize=11)
    axes[0, 1].set_title('Ecosystem Condition (Normalized)', fontsize=13, fontweight='bold')
    # Add text annotations
    for i in range(len(eco_types)):
        for j in range(len(indicators)):
            axes[0, 1].text(j, i, f'{data[i, j]:.1f}', ha='center', va='center', fontsize=9,
                           color='white' if data_norm[i, j] < 0.3 or data_norm[i, j] > 0.7 else 'black')
    plt.colorbar(im, ax=axes[0, 1])

    # Monetary account stacked bar
    eco_names = list(monetary_account.keys())
    service_types = ['Carbon value', 'Water value', 'Habitat value', 'Erosion control value', 'Pollination value']
    bottom = np.zeros(len(eco_names))
    svc_colors = ['#27ae60', '#3498db', '#e74c3c', '#f39c12', '#9b59b6']
    for svc, color in zip(service_types, svc_colors):
        vals = [monetary_account[e][svc] / 1e6 for e in eco_names]
        axes[1, 0].bar(range(len(eco_names)), vals, bottom=bottom, label=svc.replace(' value', ''),
                       color=color, alpha=0.85)
        bottom += np.array(vals)
    axes[1, 0].set_xticks(range(len(eco_names)))
    axes[1, 0].set_xticklabels(eco_names, fontsize=10)
    axes[1, 0].set_ylabel('Value (Million JPY)', fontsize=12)
    axes[1, 0].set_title('Monetary Supply-Use Account', fontsize=13, fontweight='bold')
    axes[1, 0].legend(fontsize=9)

    # Total ES value per ha
    per_ha = {}
    for name in eco_names:
        per_ha[name] = monetary_account[name]['Total ES value (JPY)'] / monetary_account[name]['Extent (ha)']
    axes[1, 1].bar(range(len(per_ha)), [v / 1000 for v in per_ha.values()],
                   color=['#95a5a6', '#f1c40f', '#27ae60', '#3498db', '#e67e22'][:len(per_ha)], alpha=0.85)
    axes[1, 1].set_xticks(range(len(per_ha)))
    axes[1, 1].set_xticklabels(list(per_ha.keys()), fontsize=10)
    axes[1, 1].set_ylabel('Value (Thousand JPY/ha/yr)', fontsize=12)
    axes[1, 1].set_title('Ecosystem Service Value per Hectare', fontsize=13, fontweight='bold')

    plt.suptitle('SEEA-EA Ecosystem Accounts', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('figures/seea_accounts.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Saved: figures/seea_accounts.png")

def plot_satoyama_analysis(comparison, scenario_values, cultural_wtp, baseline_value):
    """Plot Satoyama case study results."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # ES comparison across land use types
    eco_types = list(comparison.keys())
    services = list(comparison[eco_types[0]].keys())
    x = np.arange(len(services))
    width = 0.2
    colors = ['#f1c40f', '#27ae60', '#3498db', '#e67e22']
    for i, (eco, color) in enumerate(zip(eco_types, colors)):
        vals = [comparison[eco][s] for s in services]
        # Normalize for comparison
        max_vals = [max(comparison[e][s] for e in eco_types) for s in services]
        norm_vals = [v / m if m > 0 else 0 for v, m in zip(vals, max_vals)]
        axes[0].bar(x + i * width, norm_vals, width, label=eco, color=color, alpha=0.85)
    axes[0].set_xticks(x + width * 1.5)
    axes[0].set_xticklabels([s.split('(')[0].strip() for s in services], rotation=20, ha='right', fontsize=9)
    axes[0].set_ylabel('Normalized Service Level', fontsize=11)
    axes[0].set_title('ES Comparison by Land Use Type', fontsize=13, fontweight='bold')
    axes[0].legend(fontsize=9)

    # Scenario analysis
    scenario_names = list(scenario_values.keys())
    scenario_vals = [scenario_values[n] / 1000 for n in scenario_names]
    colors_sc = ['#2ecc71', '#e74c3c', '#c0392b', '#3498db', '#27ae60']
    bars = axes[1].barh(range(len(scenario_names)), scenario_vals, color=colors_sc, alpha=0.85)
    axes[1].set_yticks(range(len(scenario_names)))
    axes[1].set_yticklabels(scenario_names, fontsize=9)
    axes[1].set_xlabel('ES Value (Thousand JPY/ha/yr)', fontsize=11)
    axes[1].set_title('Satoyama Scenario Analysis', fontsize=13, fontweight='bold')
    axes[1].axvline(x=baseline_value / 1000, color='red', linestyle='--', alpha=0.7, label='Baseline')
    axes[1].legend()

    # Total Economic Value breakdown (TEV)
    tev_components = {
        'Direct use\n(provisioning)': baseline_value * 0.25,
        'Indirect use\n(regulating)': baseline_value * 0.40,
        'Option value': baseline_value * 0.10,
        'Cultural WTP': cultural_wtp,
        'Existence value': baseline_value * 0.15,
        'Bequest value': baseline_value * 0.10
    }
    labels = list(tev_components.keys())
    sizes = [max(v, 0) for v in tev_components.values()]
    tev_colors = ['#2ecc71', '#3498db', '#f39c12', '#e74c3c', '#9b59b6', '#1abc9c']
    axes[2].pie(sizes, labels=labels, colors=tev_colors, autopct='%1.1f%%', startangle=140, textprops={'fontsize': 9})
    axes[2].set_title('Total Economic Value Breakdown\n(Satoyama)', fontsize=13, fontweight='bold')

    plt.tight_layout()
    plt.savefig('figures/satoyama_case_study.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Saved: figures/satoyama_case_study.png")

def plot_pipeline_diagram():
    """Create InVEST/ARIES evaluation pipeline diagram."""
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 11)
    ax.axis('off')

    # Pipeline stages
    stages = [
        (2, 9.5, 'Data Collection\n& Preprocessing', '#3498db', 
         'LULC maps, DEM, Climate\nSoil, Socioeconomic data'),
        (8, 9.5, 'InVEST Modeling', '#27ae60',
         'Carbon, Water Yield\nHabitat, SDR, Pollination'),
        (14, 9.5, 'ARIES Integration', '#e67e22',
         'AI-driven model selection\nSemantic matching'),
        (2, 6.5, 'Biophysical\nQuantification', '#9b59b6',
         'Spatial ES maps\nService flow analysis'),
        (8, 6.5, 'Economic\nValuation', '#e74c3c',
         'Market prices, WTP\nChoice experiments'),
        (14, 6.5, 'SEEA-EA\nAccounting', '#1abc9c',
         'Extent, Condition\nMonetary accounts'),
        (5, 3.5, 'Policy Analysis', '#f39c12',
         'Scenario comparison\nTrade-off analysis'),
        (11, 3.5, 'Natural Capital\nBalance Sheet', '#2c3e50',
         'Asset valuation\nSustainability metrics'),
        (8, 1, 'Decision Support\n& Reporting', '#c0392b',
         'Stakeholder communication\nPolicy recommendations')
    ]

    for x, y, title, color, desc in stages:
        box = plt.Rectangle((x - 1.8, y - 0.8), 3.6, 1.6, 
                            facecolor=color, alpha=0.15, edgecolor=color, linewidth=2, 
                            joinstyle='round')
        ax.add_patch(box)
        ax.text(x, y + 0.2, title, ha='center', va='center', fontsize=11, fontweight='bold', color=color)
        ax.text(x, y - 0.4, desc, ha='center', va='center', fontsize=8, color='#2c3e50', style='italic')

    # Arrows
    arrow_props = dict(arrowstyle='->', color='#7f8c8d', lw=2)
    # Horizontal arrows (top row)
    ax.annotate('', xy=(6.2, 9.5), xytext=(3.8, 9.5), arrowprops=arrow_props)
    ax.annotate('', xy=(12.2, 9.5), xytext=(9.8, 9.5), arrowprops=arrow_props)
    # Vertical arrows
    for x_pos in [2, 8, 14]:
        ax.annotate('', xy=(x_pos, 7.3), xytext=(x_pos, 8.7), arrowprops=arrow_props)
    # Diagonal arrows to policy/balance
    ax.annotate('', xy=(5, 4.3), xytext=(3, 5.7), arrowprops=arrow_props)
    ax.annotate('', xy=(5, 4.3), xytext=(8, 5.7), arrowprops=arrow_props)
    ax.annotate('', xy=(11, 4.3), xytext=(8, 5.7), arrowprops=arrow_props)
    ax.annotate('', xy=(11, 4.3), xytext=(14, 5.7), arrowprops=arrow_props)
    # Bottom arrows
    ax.annotate('', xy=(8, 1.8), xytext=(5.5, 2.7), arrowprops=arrow_props)
    ax.annotate('', xy=(8, 1.8), xytext=(10.5, 2.7), arrowprops=arrow_props)

    ax.set_title('InVEST/ARIES-based Ecosystem Service Evaluation Pipeline', 
                 fontsize=16, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig('figures/evaluation_pipeline.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Saved: figures/evaluation_pipeline.png")

def plot_integrated_framework():
    """Create integrated framework overview."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # 1. ES Classification radar chart
    categories = ['Food\nProduction', 'Timber', 'Fresh\nwater', 'Carbon\nSeq.', 
                  'Water\nPurif.', 'Pollination', 'Flood\nReg.', 'Recreation', 'Aesthetic']
    satoyama_vals = [0.6, 0.5, 0.7, 0.55, 0.65, 0.8, 0.6, 0.75, 0.85]
    forest_vals = [0.2, 0.9, 0.8, 0.95, 0.85, 0.6, 0.7, 0.5, 0.7]
    agri_vals = [0.95, 0.1, 0.5, 0.2, 0.3, 0.4, 0.3, 0.3, 0.4]

    N = len(categories)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    ax = fig.add_subplot(221, polar=True)
    for vals, label, color in [(satoyama_vals, 'Satoyama', '#e67e22'),
                                (forest_vals, 'Forest', '#27ae60'),
                                (agri_vals, 'Agriculture', '#f1c40f')]:
        vals_plot = vals + vals[:1]
        ax.plot(angles, vals_plot, 'o-', linewidth=2, label=label, color=color)
        ax.fill(angles, vals_plot, alpha=0.1, color=color)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=8)
    ax.set_title('ES Profile by Land Use', fontsize=13, fontweight='bold', pad=20)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=9)
    axes[0, 0].axis('off')  # Remove the rectangular subplot

    # 2. Sensitivity analysis
    params = ['Carbon price', 'Discount rate', 'Population', 'Land use change', 'Climate scenario']
    low = [-25, -40, -10, -35, -20]
    high = [30, 60, 15, 25, 35]
    y_pos = np.arange(len(params))
    axes[0, 1].barh(y_pos, high, left=0, color='#e74c3c', alpha=0.7, label='Increase')
    axes[0, 1].barh(y_pos, low, left=0, color='#3498db', alpha=0.7, label='Decrease')
    axes[0, 1].set_yticks(y_pos)
    axes[0, 1].set_yticklabels(params, fontsize=10)
    axes[0, 1].set_xlabel('% Change in Total ES Value', fontsize=11)
    axes[0, 1].set_title('Tornado Sensitivity Analysis', fontsize=13, fontweight='bold')
    axes[0, 1].legend(fontsize=9)
    axes[0, 1].axvline(x=0, color='black', linewidth=1)

    # 3. Trade-off frontier
    carbon_seq = np.linspace(50, 200, 20)
    food_prod = 180 - 0.8 * carbon_seq + np.random.normal(0, 5, 20)
    axes[1, 0].scatter(carbon_seq, food_prod, c=np.linspace(0, 1, 20), cmap='RdYlGn', 
                       s=100, edgecolors='gray', zorder=5)
    z = np.polyfit(carbon_seq, food_prod, 2)
    p = np.poly1d(z)
    x_smooth = np.linspace(50, 200, 100)
    axes[1, 0].plot(x_smooth, p(x_smooth), '--', color='gray', alpha=0.7)
    axes[1, 0].set_xlabel('Carbon Sequestration (tC/ha)', fontsize=11)
    axes[1, 0].set_ylabel('Food Production Index', fontsize=11)
    axes[1, 0].set_title('ES Trade-off Frontier', fontsize=13, fontweight='bold')
    axes[1, 0].annotate('Satoyama\n(balanced)', xy=(100, 105), fontsize=10, fontweight='bold',
                        color='#e67e22', ha='center',
                        arrowprops=dict(arrowstyle='->', color='#e67e22'),
                        xytext=(130, 130))

    # 4. Temporal dynamics
    years_t = np.arange(2020, 2051)
    baseline_es = 100 * np.exp(-0.005 * (years_t - 2020))
    managed_es = 100 * (1 + 0.008 * (years_t - 2020))
    degraded_es = 100 * np.exp(-0.02 * (years_t - 2020))
    axes[1, 1].plot(years_t, baseline_es, '-', linewidth=2, label='Business as usual', color='#e74c3c')
    axes[1, 1].plot(years_t, managed_es, '-', linewidth=2, label='Active management', color='#27ae60')
    axes[1, 1].plot(years_t, degraded_es, '-', linewidth=2, label='Degradation scenario', color='#95a5a6')
    axes[1, 1].fill_between(years_t, managed_es, degraded_es, alpha=0.1, color='green')
    axes[1, 1].set_xlabel('Year', fontsize=11)
    axes[1, 1].set_ylabel('ES Value Index (2020=100)', fontsize=11)
    axes[1, 1].set_title('Temporal ES Value Projections', fontsize=13, fontweight='bold')
    axes[1, 1].legend(fontsize=9)
    axes[1, 1].grid(True, alpha=0.3)

    plt.suptitle('Integrated Ecosystem Service Valuation Framework', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('figures/integrated_framework.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Saved: figures/integrated_framework.png")


# ============================================================
# MAIN EXECUTION
# ============================================================

if __name__ == '__main__':
    print("=" * 60)
    print("Integrated Framework for Economic Valuation of ES")
    print("=" * 60)

    # 1. Classification
    print("\n[1] ES Classification & Indicators...")
    es_classification = create_es_classification()
    for cat, services in es_classification.items():
        print(f"  {cat}: {len(services)} services")

    # 2. InVEST spatial modeling
    print("\n[2] InVEST Spatial Service Quantification...")
    invest_data = simulate_invest_spatial(grid_size=50)
    print(f"  Grid: {invest_data['grid_size']}x{invest_data['grid_size']}")
    print(f"  Carbon range: {invest_data['carbon'].min():.1f} - {invest_data['carbon'].max():.1f} tC/ha")
    print(f"  Water yield range: {invest_data['water_yield'].min():.1f} - {invest_data['water_yield'].max():.1f} mm/yr")
    print(f"  Habitat quality range: {invest_data['habitat'].min():.3f} - {invest_data['habitat'].max():.3f}")

    # 3. Choice experiment
    print("\n[3] Choice Experiment & WTP Estimation...")
    ce_data, estimated_beta, wtp, true_beta = simulate_choice_experiment(n_respondents=500)
    print(f"  Respondents: 500, Choice sets: 12")
    print(f"  Estimated coefficients:")
    for attr, val in estimated_beta.items():
        print(f"    β_{attr}: {val:.6f} (true: {true_beta[attr]:.6f})")
    print(f"  Marginal WTP estimates (JPY/household/year):")
    for attr, val in wtp.items():
        print(f"    {attr}: ¥{val:,.0f}")

    # 4. Discount rate analysis
    print("\n[4] Discount Rate & Intergenerational Equity...")
    years, pv_series, npv_results = analyze_discount_rates()
    for name, npv in npv_results.items():
        print(f"  {name}: NPV = ¥{npv/1e6:.1f}M")

    # 5. SEEA-EA accounts
    print("\n[5] SEEA-EA Natural Capital Accounting...")
    extent, condition, monetary = create_seea_accounts(invest_data)
    print(f"  Ecosystem types: {len(extent)}")
    for name, ha in extent.items():
        total_val = monetary.get(name, {}).get('Total ES value (JPY)', 0)
        print(f"    {name}: {ha} ha, Total ES value: ¥{total_val/1e6:.1f}M")

    # 6. Satoyama case study
    print("\n[6] Satoyama Case Study...")
    sat_services, comparison, scenario_vals, cultural_wtp, baseline_val = satoyama_case_study(invest_data, wtp)
    print(f"  Baseline ES value: ¥{baseline_val/1000:.1f}K/ha/yr")
    print(f"  Cultural WTP (landscape + recreation): ¥{cultural_wtp:,.0f}/household/yr")
    for name, val in scenario_vals.items():
        change = (val - baseline_val) / baseline_val * 100
        print(f"    {name}: ¥{val/1000:.1f}K ({change:+.1f}%)")

    # Generate all figures
    print("\n[Plotting] Generating figures...")
    plot_spatial_maps(invest_data)
    plot_wtp_results(estimated_beta, wtp, true_beta)
    plot_discount_analysis(years, pv_series, npv_results)
    plot_seea_accounts(extent, condition, monetary)
    plot_satoyama_analysis(comparison, scenario_vals, cultural_wtp, baseline_val)
    plot_pipeline_diagram()
    plot_integrated_framework()

    print("\n" + "=" * 60)
    print("All experiments completed successfully!")
    print("=" * 60)
