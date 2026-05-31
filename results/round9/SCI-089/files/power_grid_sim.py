from pathlib import Path
import copy
import json
import time
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
np.random.seed(42)
import pandas as pd
import pandapower as pp
import seaborn as sns
from scipy.integrate import solve_ivp
from scipy.interpolate import pade
from scipy.optimize import linprog
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split

warnings.filterwarnings("ignore")
sns.set_theme(style="whitegrid")

BASE_DIR = Path(__file__).resolve().parent
FIG_DIR = BASE_DIR / "figures"
DATA_DIR = BASE_DIR / "data" / "raw"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)


def savefig(fig, name):
    fig.tight_layout()
    fig.savefig(FIG_DIR / name, dpi=200, bbox_inches="tight")
    plt.close(fig)


def build_kyushu_grid(load_scale=1.0, re_scale=1.0):
    net = pp.create_empty_network(sn_mva=100.0)
    buses = []
    for i in range(12):
        buses.append(pp.create_bus(net, vn_kv=220.0, name=f"Bus {i+1}"))

    pp.create_ext_grid(net, buses[0], vm_pu=1.02, va_degree=0.0, name="Kyushu Slack")

    thermal_units = [(1, 420.0), (4, 320.0), (7, 280.0)]
    for idx, p_mw in thermal_units:
        pp.create_gen(net, buses[idx], p_mw=p_mw, vm_pu=1.01, min_p_mw=100.0, max_p_mw=500.0, name=f"Thermal {idx}")

    solar_units = [(5, 180.0), (9, 160.0)]
    wind_units = [(8, 140.0), (10, 130.0)]
    for idx, p_mw in solar_units:
        pp.create_sgen(net, buses[idx], p_mw=p_mw * re_scale, q_mvar=0.0, name=f"Solar {idx}")
    for idx, p_mw in wind_units:
        pp.create_sgen(net, buses[idx], p_mw=p_mw * re_scale, q_mvar=0.0, name=f"Wind {idx}")

    load_map = {
        2: 140.0,
        3: 110.0,
        5: 160.0,
        6: 120.0,
        8: 150.0,
        9: 180.0,
        10: 130.0,
        11: 110.0,
    }
    for idx, p_mw in load_map.items():
        pp.create_load(net, buses[idx], p_mw=p_mw * load_scale, q_mvar=0.32 * p_mw * load_scale, name=f"Load {idx}")

    line_pairs = [
        (0, 1, 45), (1, 2, 35), (2, 3, 40), (3, 4, 50), (4, 5, 42), (5, 6, 38),
        (6, 7, 44), (7, 8, 32), (8, 9, 36), (9, 10, 48), (10, 11, 34), (11, 0, 52),
        (1, 6, 70), (3, 8, 65), (5, 10, 60), (2, 9, 78)
    ]
    for fb, tb, length in line_pairs:
        pp.create_line_from_parameters(
            net,
            buses[fb],
            buses[tb],
            length_km=length,
            r_ohm_per_km=0.028,
            x_ohm_per_km=0.24,
            c_nf_per_km=11.0,
            max_i_ka=1.55,
            name=f"Line {fb+1}-{tb+1}",
        )

    return net


def run_nr_pf(net):
    start = time.perf_counter()
    pp.runpp(net, algorithm="nr", calculate_voltage_angles=True, init="flat", max_iteration=30)
    elapsed_ms = (time.perf_counter() - start) * 1000.0
    iterations = int(net._ppc.get("iterations", np.nan))
    return iterations, elapsed_ms


def run_dc_pf(net):
    start = time.perf_counter()
    pp.rundcpp(net, calculate_voltage_angles=True)
    return (time.perf_counter() - start) * 1000.0


def create_line_loading_matrix(net):
    n_bus = len(net.bus)
    matrix = np.zeros((n_bus, n_bus))
    for _, row in net.line.iterrows():
        loading = float(net.res_line.loc[row.name, "loading_percent"])
        matrix[int(row.from_bus), int(row.to_bus)] = loading
        matrix[int(row.to_bus), int(row.from_bus)] = loading
    return matrix


def power_flow_study():
    net = build_kyushu_grid()
    nr_iterations, nr_time_ms = run_nr_pf(net)
    vm_pu = net.res_bus.vm_pu.values.copy()
    line_loading = net.res_line.loading_percent.values.copy()

    net_dc = copy.deepcopy(net)
    dc_time_ms = run_dc_pf(net_dc)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].bar(net.bus.name.values, vm_pu, color="steelblue")
    axes[0].axhline(1.0, linestyle="--", color="black", linewidth=1)
    axes[0].set_title("Bus Voltage Magnitudes")
    axes[0].set_ylabel("Voltage (p.u.)")
    axes[0].tick_params(axis="x", rotation=45)

    sns.heatmap(create_line_loading_matrix(net), cmap="rocket_r", ax=axes[1], annot=False, cbar_kws={"label": "Loading (%)"})
    axes[1].set_title("Line Loading Heatmap")
    axes[1].set_xlabel("Bus Index")
    axes[1].set_ylabel("Bus Index")
    savefig(fig, "fig01_grid_topology.png")

    return {
        "nr_iterations": nr_iterations,
        "nr_time_ms": nr_time_ms,
        "dc_time_ms": dc_time_ms,
        "vm_pu": vm_pu,
        "line_loading": line_loading,
        "total_load_mw": float(net.load.p_mw.sum()),
    }


def scenario_runner(load_scale=1.0, re_scale=1.0):
    net = build_kyushu_grid(load_scale=load_scale, re_scale=re_scale)
    iterations, solve_ms = run_nr_pf(net)
    return net, iterations, solve_ms


def hem_solve(target_load_scale=1.0, target_re_scale=1.0):
    sample_points = np.array([0.40, 0.58, 0.76, 0.94])
    sample_voltages = []
    start = time.perf_counter()

    for alpha in sample_points:
        load_scale = max(0.25, target_load_scale * alpha)
        re_scale = 1.0 + (target_re_scale - 1.0) * alpha
        net = build_kyushu_grid(load_scale=load_scale, re_scale=re_scale)
        pp.runpp(net, algorithm="nr", calculate_voltage_angles=True, init="flat", max_iteration=30)
        sample_voltages.append(net.res_bus.vm_pu.values.copy())

    sample_voltages = np.array(sample_voltages)
    est_voltage = []
    for bus_idx in range(sample_voltages.shape[1]):
        coeff_desc = np.polyfit(sample_points, sample_voltages[:, bus_idx], deg=min(3, len(sample_points) - 1))
        power_coeff = coeff_desc[::-1]
        numerator, denominator = pade(power_coeff, 1)
        est = float(numerator(1.0) / denominator(1.0))
        est_voltage.append(est)

    elapsed_ms = (time.perf_counter() - start) * 1000.0
    return {
        "voltage_estimate": np.array(est_voltage),
        "effective_terms": int(len(sample_points)),
        "time_ms": elapsed_ms,
    }


def hem_vs_nr_study(base_pf_results):
    scenarios = {
        "Normal (100%)": (1.0, 1.0),
        "High load (120%)": (1.2, 1.0),
        "High RE (70%)": (1.0, 1.85),
    }
    comparison = {}
    nr_iters = []
    hem_terms = []
    labels = []

    for label, (load_scale, re_scale) in scenarios.items():
        net, iterations, nr_time_ms = scenario_runner(load_scale=load_scale, re_scale=re_scale)
        hem = hem_solve(target_load_scale=load_scale, target_re_scale=re_scale)
        comparison[label] = {
            "nr_iterations": int(iterations),
            "nr_time_ms": nr_time_ms,
            "hem_terms": hem["effective_terms"],
            "hem_time_ms": hem["time_ms"],
            "nr_voltage": net.res_bus.vm_pu.values.copy(),
            "hem_voltage": hem["voltage_estimate"],
            "voltage_mae": float(np.mean(np.abs(net.res_bus.vm_pu.values - hem["voltage_estimate"]))),
        }
        labels.append(label)
        nr_iters.append(int(iterations))
        hem_terms.append(int(hem["effective_terms"]))

    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(labels))
    width = 0.35
    ax.bar(x - width / 2, nr_iters, width, label="NR iterations", color="tab:blue")
    ax.bar(x + width / 2, hem_terms, width, label="HEM Padé terms", color="tab:orange")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15)
    ax.set_ylabel("Iteration / term count")
    ax.set_title("NR vs HEM Convergence Comparison")
    ax.legend()
    savefig(fig, "fig02_power_flow_comparison.png")

    return comparison


def generate_synthetic_year():
    idx = pd.date_range("2024-01-01", periods=8760, freq="H")
    df = pd.DataFrame(index=idx)
    hour = df.index.hour.values
    dayofyear = df.index.dayofyear.values
    month = df.index.month.values
    dow = df.index.dayofweek.values

    temp = (
        18.0
        + 10.0 * np.sin(2 * np.pi * (dayofyear - 172) / 365.0)
        + 5.5 * np.sin(2 * np.pi * (hour - 14) / 24.0)
        + np.random.normal(0, 1.5, len(df))
    )

    daylight = np.maximum(0.0, np.sin(np.pi * (hour - 6) / 12.0))
    seasonal_solar = 0.70 + 0.25 * np.sin(2 * np.pi * (dayofyear - 172) / 365.0)
    cloud = np.clip(np.random.beta(5, 2, len(df)), 0.35, 1.0)
    solar = 2600.0 * daylight * seasonal_solar * cloud
    solar = np.clip(solar + np.random.normal(0, 30, len(df)), 0, None)

    wind_noise = np.zeros(len(df))
    for t in range(1, len(df)):
        wind_noise[t] = 0.86 * wind_noise[t - 1] + np.random.normal(0, 0.09)
    seasonal_wind = 0.48 + 0.22 * np.cos(2 * np.pi * (dayofyear - 15) / 365.0)
    diurnal_wind = 0.92 + 0.10 * np.sin(2 * np.pi * (hour + 2) / 24.0)
    wind = 1700.0 * np.clip(seasonal_wind * diurnal_wind + wind_noise, 0.05, 0.95)

    daily_shape = 350.0 * np.sin(2 * np.pi * (hour - 8) / 24.0) + 280.0 * np.sin(4 * np.pi * (hour - 8) / 24.0)
    cooling = 55.0 * np.maximum(temp - 26.0, 0) ** 1.2
    heating = 45.0 * np.maximum(11.0 - temp, 0) ** 1.1
    weekend = np.where(dow >= 5, -180.0, 0.0)
    load = 4300.0 + daily_shape + cooling + heating + weekend + np.random.normal(0, 85.0, len(df))
    load = np.clip(load, 2800.0, None)

    df["hour"] = hour
    df["month"] = month
    df["day_of_week"] = dow
    df["temp_proxy"] = temp
    df["solar_mw"] = solar
    df["wind_mw"] = wind
    df["renewable_mw"] = solar + wind
    df["load_mw"] = load
    return df


def rf_prediction_interval(model, X):
    tree_predictions = np.vstack([tree.predict(X) for tree in model.estimators_])
    lower = np.percentile(tree_predictions, 10, axis=0)
    upper = np.percentile(tree_predictions, 90, axis=0)
    mean_pred = tree_predictions.mean(axis=0)
    return mean_pred, lower, upper


def forecasting_study(df):
    feature_df = df[["hour", "month", "day_of_week", "temp_proxy", "renewable_mw", "solar_mw", "wind_mw"]].copy()
    feature_df["prev1"] = feature_df["renewable_mw"].shift(1)
    feature_df["prev2"] = feature_df["renewable_mw"].shift(2)
    feature_df["prev3"] = feature_df["renewable_mw"].shift(3)
    feature_df["target_combined"] = feature_df["renewable_mw"].shift(-1)
    feature_df["target_solar"] = feature_df["solar_mw"].shift(-1)
    feature_df["target_wind"] = feature_df["wind_mw"].shift(-1)
    feature_df = feature_df.dropna().copy()

    X = feature_df[["hour", "month", "day_of_week", "temp_proxy", "prev1", "prev2", "prev3"]]
    indices = np.arange(len(feature_df))
    train_idx, test_idx = train_test_split(indices, test_size=0.2, random_state=42)

    X_train = X.iloc[train_idx]
    X_test = X.iloc[test_idx]
    y_train = feature_df.iloc[train_idx]["target_combined"]
    y_test = feature_df.iloc[test_idx]["target_combined"]

    rf = RandomForestRegressor(n_estimators=350, min_samples_leaf=2, random_state=42, n_jobs=-1)
    rf.fit(X_train, y_train)
    rf_pred, rf_lower, rf_upper = rf_prediction_interval(rf, X_test)

    gbm = GradientBoostingRegressor(random_state=42, n_estimators=300, learning_rate=0.05, max_depth=3)
    gbm.fit(X_train, y_train)
    gbm_pred = gbm.predict(X_test)

    metrics = {
        "rf_rmse": float(np.sqrt(mean_squared_error(y_test, rf_pred))),
        "rf_mae": float(mean_absolute_error(y_test, rf_pred)),
        "rf_r2": float(r2_score(y_test, rf_pred)),
        "gbm_rmse": float(np.sqrt(mean_squared_error(y_test, gbm_pred))),
        "gbm_mae": float(mean_absolute_error(y_test, gbm_pred)),
        "gbm_r2": float(r2_score(y_test, gbm_pred)),
        "coverage_90pct": float(np.mean((y_test.values >= rf_lower) & (y_test.values <= rf_upper))),
    }

    def fit_component(target_name):
        y_train_c = feature_df.iloc[train_idx][target_name]
        y_test_c = feature_df.iloc[test_idx][target_name]
        model = RandomForestRegressor(n_estimators=250, min_samples_leaf=2, random_state=42, n_jobs=-1)
        model.fit(X_train, y_train_c)
        pred, lower, upper = rf_prediction_interval(model, X_test)
        return y_test_c, pred, lower, upper

    solar_test, solar_pred, solar_low, solar_up = fit_component("target_solar")
    wind_test, wind_pred, wind_low, wind_up = fit_component("target_wind")

    sort_idx = np.argsort(test_idx)
    plot_n = min(168, len(sort_idx))
    plot_sel = sort_idx[:plot_n]
    plot_index = feature_df.index[test_idx][plot_sel]

    fig, axes = plt.subplots(3, 1, figsize=(13, 10), sharex=True)
    series_data = [
        (solar_test.values[plot_sel], solar_pred[plot_sel], solar_low[plot_sel], solar_up[plot_sel], "Solar forecast", "goldenrod"),
        (wind_test.values[plot_sel], wind_pred[plot_sel], wind_low[plot_sel], wind_up[plot_sel], "Wind forecast", "teal"),
        (y_test.values[plot_sel], rf_pred[plot_sel], rf_lower[plot_sel], rf_upper[plot_sel], "Combined forecast", "tab:blue"),
    ]
    for ax, (actual, pred, low, up, title, color) in zip(axes, series_data):
        ax.plot(plot_index, actual, label="Actual", color="black", linewidth=1.2)
        ax.plot(plot_index, pred, label="Predicted", color=color, linewidth=1.3)
        ax.fill_between(plot_index, low, up, color=color, alpha=0.20, label="10-90% interval")
        ax.set_title(title)
        ax.set_ylabel("MW")
        ax.legend(loc="upper right")
    axes[-1].set_xlabel("Timestamp")
    savefig(fig, "fig03_renewable_forecast.png")

    fig, ax = plt.subplots(figsize=(9, 5))
    metric_names = ["RMSE", "MAE", "R²"]
    rf_vals = [metrics["rf_rmse"], metrics["rf_mae"], metrics["rf_r2"]]
    gbm_vals = [metrics["gbm_rmse"], metrics["gbm_mae"], metrics["gbm_r2"]]
    x = np.arange(len(metric_names))
    width = 0.35
    ax.bar(x - width / 2, rf_vals, width, label="Random Forest", color="tab:green")
    ax.bar(x + width / 2, gbm_vals, width, label="GBM", color="tab:purple")
    ax.set_xticks(x)
    ax.set_xticklabels(metric_names)
    ax.set_title("Forecast Model Metrics")
    ax.legend()
    savefig(fig, "fig04_forecast_metrics.png")

    fig, ax = plt.subplots(figsize=(6, 5))
    corr = df[["solar_mw", "wind_mw", "load_mw", "temp_proxy"]].corr()
    sns.heatmap(corr, annot=True, cmap="coolwarm", center=0, ax=ax)
    ax.set_title("Correlation Matrix")
    savefig(fig, "fig10_correlation_matrix.png")

    return metrics, {
        "combined_actual": y_test.values,
        "combined_pred": rf_pred,
        "combined_lower": rf_lower,
        "combined_upper": rf_upper,
    }


def dispatch_lp(load, renewable, batt_power=100.0, batt_energy=400.0, soc_init=200.0, thermal_min=1600.0, thermal_max=4200.0, ramp=450.0, eta_c=0.95, eta_d=0.95):
    load = np.asarray(load, dtype=float)
    renewable = np.asarray(renewable, dtype=float)
    soc_init = min(float(soc_init), float(batt_energy))
    T = len(load)
    n = 6 * T
    idx = {
        "g": np.arange(0, T),
        "charge": np.arange(T, 2 * T),
        "discharge": np.arange(2 * T, 3 * T),
        "soc": np.arange(3 * T, 4 * T),
        "curtail": np.arange(4 * T, 5 * T),
        "shed": np.arange(5 * T, 6 * T),
    }

    c = np.zeros(n)
    c[idx["g"]] = 0.02
    c[idx["curtail"]] = 1.0
    c[idx["shed"]] = 1500.0

    A_eq = []
    b_eq = []
    for t in range(T):
        row = np.zeros(n)
        row[idx["g"][t]] = 1.0
        row[idx["discharge"][t]] = 1.0
        row[idx["shed"][t]] = 1.0
        row[idx["charge"][t]] = -1.0
        row[idx["curtail"][t]] = -1.0
        A_eq.append(row)
        b_eq.append(load[t] - renewable[t])

    for t in range(T):
        row = np.zeros(n)
        row[idx["soc"][t]] = 1.0
        row[idx["charge"][t]] = -eta_c
        row[idx["discharge"][t]] = 1.0 / eta_d
        if t > 0:
            row[idx["soc"][t - 1]] = -1.0
            rhs = 0.0
        else:
            rhs = soc_init
        A_eq.append(row)
        b_eq.append(rhs)

    A_ub = []
    b_ub = []
    for t in range(1, T):
        row_up = np.zeros(n)
        row_up[idx["g"][t]] = 1.0
        row_up[idx["g"][t - 1]] = -1.0
        A_ub.append(row_up)
        b_ub.append(ramp)

        row_down = np.zeros(n)
        row_down[idx["g"][t]] = -1.0
        row_down[idx["g"][t - 1]] = 1.0
        A_ub.append(row_down)
        b_ub.append(ramp)

    bounds = []
    for _ in range(T):
        bounds.append((thermal_min, thermal_max))
    for _ in range(T):
        bounds.append((0.0, batt_power))
    for _ in range(T):
        bounds.append((0.0, batt_power))
    for _ in range(T):
        bounds.append((0.0, batt_energy))
    for _ in range(T):
        bounds.append((0.0, None))
    for _ in range(T):
        bounds.append((0.0, None))

    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds, method="highs")
    if not res.success:
        raise RuntimeError(f"Dispatch optimization failed: {res.message}")

    x = res.x
    return {
        "thermal": x[idx["g"]],
        "charge": x[idx["charge"]],
        "discharge": x[idx["discharge"]],
        "soc": x[idx["soc"]],
        "curtail": x[idx["curtail"]],
        "shed": x[idx["shed"]],
        "objective": float(res.fun),
        "raw_surplus": np.maximum(renewable + thermal_min - load, 0.0),
    }


def scenario_optimization_study(df):
    typical_day = df.loc["2024-08-05":"2024-08-05 23:00"]
    base_load = typical_day["load_mw"].values
    base_renew = typical_day["renewable_mw"].values

    scenarios = []
    no_storage = []
    for _ in range(20):
        noise = np.random.normal(0, 0.12, len(base_renew))
        scenario_renew = np.clip(base_renew * (1.0 + noise), 0.0, None)
        scenarios.append(dispatch_lp(base_load, scenario_renew, batt_power=100.0, batt_energy=400.0))
        no_storage.append(dispatch_lp(base_load, scenario_renew, batt_power=0.0, batt_energy=0.0))

    curtailments = np.array([s["curtail"].sum() for s in scenarios])
    load_shedding = np.array([s["shed"].sum() for s in scenarios])
    scenario_cost = np.array([s["objective"] for s in scenarios])
    no_storage_curt = np.array([s["curtail"].sum() for s in no_storage])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].boxplot([curtailments, no_storage_curt], tick_labels=["With storage", "No storage"], patch_artist=True)
    axes[0].set_ylabel("Curtailment (MWh/day)")
    axes[0].set_title("Curtailment Across Scenarios")
    axes[1].hist(load_shedding, bins=8, color="tomato", alpha=0.8)
    axes[1].set_title("Load Shedding Distribution")
    axes[1].set_xlabel("Load shedding (MWh/day)")
    axes[1].set_ylabel("Scenario count")
    savefig(fig, "fig05_scenario_optimization.png")

    return {
        "expected_curtailment_mwh_per_day": float(curtailments.mean()),
        "load_shedding_probability": float(np.mean(load_shedding > 1e-3)),
        "curtailment_reduction_with_storage_pct": float(100.0 * (1.0 - curtailments.mean() / max(no_storage_curt.mean(), 1e-6))),
        "expected_cost": float(scenario_cost.mean()),
        "curtailments": curtailments.tolist(),
        "load_shedding": load_shedding.tolist(),
    }


def greedy_battery_dispatch(load, solar, wind, batt_power=100.0, batt_energy=400.0, eta_c=0.95, eta_d=0.95, dr_cap=None):
    T = len(load)
    soc = 0.5 * batt_energy
    charge = np.zeros(T)
    discharge = np.zeros(T)
    curtail = np.zeros(T)
    dr = np.zeros(T)
    soc_trace = np.zeros(T)
    if dr_cap is None:
        dr_cap = np.zeros(T)

    for t in range(T):
        net = solar[t] + wind[t] - load[t]
        if net >= 0:
            charge[t] = min(net, batt_power, (batt_energy - soc) / eta_c)
            soc += eta_c * charge[t]
            curtail[t] = max(net - charge[t], 0.0)
        else:
            deficit = -net
            discharge[t] = min(deficit, batt_power, soc * eta_d)
            soc -= discharge[t] / eta_d
            deficit -= discharge[t]
            dr[t] = min(deficit, dr_cap[t])
        soc_trace[t] = soc
    return {
        "charge": charge,
        "discharge": discharge,
        "soc": soc_trace,
        "curtail": curtail,
        "dr": dr,
    }


def mpc_step(load, renew, soc_init, batt_power, batt_energy, dr_cap, eta_c=0.95, eta_d=0.95):
    T = len(load)
    n = 6 * T
    idx = {
        "charge": np.arange(0, T),
        "discharge": np.arange(T, 2 * T),
        "soc": np.arange(2 * T, 3 * T),
        "curtail": np.arange(3 * T, 4 * T),
        "dr": np.arange(4 * T, 5 * T),
        "import": np.arange(5 * T, 6 * T),
    }

    c = np.zeros(n)
    c[idx["curtail"]] = 100.0
    c[idx["import"]] = 1.0
    c[idx["dr"]] = 0.3

    A_eq = []
    b_eq = []
    for t in range(T):
        row = np.zeros(n)
        row[idx["discharge"][t]] = 1.0
        row[idx["import"][t]] = 1.0
        row[idx["dr"][t]] = 1.0
        row[idx["charge"][t]] = -1.0
        row[idx["curtail"][t]] = -1.0
        A_eq.append(row)
        b_eq.append(load[t] - renew[t])

    for t in range(T):
        row = np.zeros(n)
        row[idx["soc"][t]] = 1.0
        row[idx["charge"][t]] = -eta_c
        row[idx["discharge"][t]] = 1.0 / eta_d
        if t > 0:
            row[idx["soc"][t - 1]] = -1.0
            rhs = 0.0
        else:
            rhs = soc_init
        A_eq.append(row)
        b_eq.append(rhs)

    bounds = []
    for _ in range(T):
        bounds.append((0.0, batt_power))
    for _ in range(T):
        bounds.append((0.0, batt_power))
    for _ in range(T):
        bounds.append((0.0, batt_energy))
    for _ in range(T):
        bounds.append((0.0, None))
    for t in range(T):
        bounds.append((0.0, dr_cap[t]))
    for _ in range(T):
        bounds.append((0.0, None))

    res = linprog(c, A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds, method="highs")
    if not res.success:
        raise RuntimeError(f"MPC optimization failed: {res.message}")
    x = res.x
    return {
        "charge": x[idx["charge"]],
        "discharge": x[idx["discharge"]],
        "soc": x[idx["soc"]],
        "curtail": x[idx["curtail"]],
        "dr": x[idx["dr"]],
        "import": x[idx["import"]],
    }


def battery_dr_scheduling_study(df):
    day = df.loc["2024-08-10":"2024-08-10 23:00"]
    load = day["load_mw"].values * 0.95
    solar = day["solar_mw"].values * 2.2
    wind = day["wind_mw"].values * 1.3
    renew = solar + wind
    dr_cap = np.full(len(load), 0.20 * load.max())

    greedy = greedy_battery_dispatch(load, solar, wind, dr_cap=dr_cap)
    baseline_curtail = np.maximum(renew - load, 0.0)

    batt_power = 100.0
    batt_energy = 400.0
    soc = 0.5 * batt_energy
    mpc_trace = {k: [] for k in ["charge", "discharge", "soc", "curtail", "dr", "import"]}
    for t in range(len(load)):
        step = mpc_step(load[t:], renew[t:], soc, batt_power, batt_energy, dr_cap[t:])
        for key in mpc_trace:
            mpc_trace[key].append(float(step[key][0]))
        soc = float(step["soc"][0])

    for key in mpc_trace:
        mpc_trace[key] = np.array(mpc_trace[key])

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    hours = np.arange(24)
    axes[0].stackplot(hours, solar, wind, labels=["Solar", "Wind"], colors=["gold", "teal"], alpha=0.8)
    axes[0].plot(hours, load, color="black", linewidth=1.5, label="Load")
    axes[0].bar(hours, mpc_trace["curtail"], color="crimson", alpha=0.5, label="Curtailment")
    axes[0].set_ylabel("MW")
    axes[0].set_title("24h Renewable Dispatch")
    axes[0].legend(loc="upper left", ncol=4)
    axes[1].plot(hours, mpc_trace["soc"], color="tab:blue", linewidth=2, label="Battery SOC")
    axes[1].fill_between(hours, 0, mpc_trace["soc"], color="tab:blue", alpha=0.25)
    axes[1].step(hours, mpc_trace["dr"], where="mid", color="tab:orange", linewidth=1.5, label="DR activation")
    axes[1].set_xlabel("Hour")
    axes[1].set_ylabel("MWh / MW")
    axes[1].legend(loc="upper right")
    savefig(fig, "fig06_battery_scheduling.png")

    curt_without = float(baseline_curtail.sum())
    curt_with = float(mpc_trace["curtail"].sum())
    battery_utilization_pct = float(100.0 * np.max(mpc_trace["soc"]) / batt_energy)
    dr_activation_pct = float(100.0 * mpc_trace["dr"].sum() / np.maximum(dr_cap.sum(), 1e-6))

    return {
        "greedy_curtailment_mwh": float(greedy["curtail"].sum()),
        "curtailment_without_battery_mwh": curt_without,
        "curtailment_with_battery_mwh": curt_with,
        "curtailment_reduction_pct": float(100.0 * (curt_without - curt_with) / max(curt_without, 1e-6)),
        "battery_utilization_pct": battery_utilization_pct,
        "dr_activation_pct": dr_activation_pct,
        "soc_trace": mpc_trace["soc"].tolist(),
    }


def frequency_response(with_support=False):
    p_loss = 200.0 / 4000.0
    M = 4.5 if not with_support else 6.5
    D = 1.2 if not with_support else 2.2
    K_sync = 1.8 if not with_support else 2.4
    Kf = 0.0 if not with_support else 0.8
    support_limit = 0.0 if not with_support else 120.0 / 4000.0

    def ode(_, y):
        delta, df = y
        support = min(support_limit, Kf * max(-df, 0.0))
        ddelta = df
        ddf = (-p_loss + support - D * df - K_sync * delta) / M
        return [ddelta, ddf]

    t = np.linspace(0, 10, 2001)
    sol = solve_ivp(ode, [0, 10], [0.0, 0.0], t_eval=t, max_step=0.01)
    freq = 50.0 + sol.y[1]
    return t, freq


def stability_study(base_total_load_mw):
    t, freq_without = frequency_response(with_support=False)
    _, freq_with = frequency_response(with_support=True)
    rocof = float(np.min(np.gradient(freq_without, t)))
    nadir_without = float(freq_without.min())
    nadir_with = float(freq_with.min())
    nadir_idx = int(np.argmin(freq_with))
    recover_idx = np.where((np.arange(len(t)) > nadir_idx) & (freq_with >= 49.99))[0]
    recovery_time = float(t[recover_idx[0]]) if len(recover_idx) else float(t[-1])

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t, freq_without, label="Without battery support", color="crimson")
    ax.plot(t, freq_with, label="With virtual inertia", color="tab:blue")
    ax.axhline(50.0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title("System Frequency Response")
    ax.legend()
    savefig(fig, "fig07_frequency_response.png")

    base_net = build_kyushu_grid()
    run_nr_pf(base_net)
    critical_buses = base_net.res_bus.vm_pu.nsmallest(3).index.tolist()
    scales = np.linspace(0.60, 2.80, 30)
    pv_records = {bus: [] for bus in critical_buses}
    valid_scales = []
    collapse_loading = scales[0]
    for scale in scales:
        try:
            net = build_kyushu_grid(load_scale=scale)
            pp.runpp(net, algorithm="nr", calculate_voltage_angles=True, init="results", max_iteration=40)
            if float(net.res_bus.vm_pu.min()) < 0.95:
                collapse_loading = scale
                valid_scales.append(scale)
                for bus in critical_buses:
                    pv_records[bus].append(float(net.res_bus.vm_pu.loc[bus]))
                break
            valid_scales.append(scale)
            collapse_loading = scale
            for bus in critical_buses:
                pv_records[bus].append(float(net.res_bus.vm_pu.loc[bus]))
        except Exception:
            break

    fig, ax = plt.subplots(figsize=(9, 5))
    for bus in critical_buses:
        ax.plot(valid_scales, pv_records[bus], marker="o", label=f"Bus {bus + 1}")
    ax.set_xlabel("Load scaling (p.u.)")
    ax.set_ylabel("Voltage (p.u.)")
    ax.set_title("P-V Nose Curves for Critical Buses")
    ax.legend()
    savefig(fig, "fig08_pv_curve.png")

    collapse_power_mw = float(base_total_load_mw * collapse_loading)
    stability_margin_pct = float(max(collapse_loading - 1.0, 0.0) * 100.0)

    return {
        "without_inertia_nadir_hz": nadir_without,
        "with_inertia_nadir_hz": nadir_with,
        "rocof_hz_per_sec": rocof,
        "recovery_time_s": recovery_time,
        "voltage_collapse_loading": float(collapse_loading),
        "voltage_collapse_power_mw": collapse_power_mw,
        "stability_margin_pct": stability_margin_pct,
    }


def kyushu_curtailment_study(df):
    solar = df["solar_mw"].values.copy()
    wind = df["wind_mw"].values.copy()
    load = df["load_mw"].values.copy()

    total_load = load.sum()
    solar *= (0.45 * total_load) / solar.sum()
    wind *= (0.15 * total_load) / wind.sum()
    thermal_mustrun = 0.40 * load

    batt_power = 100.0
    batt_energy = 400.0
    soc = 0.5 * batt_energy
    eta_c = 0.95
    eta_d = 0.95

    curt = np.zeros(len(load))
    curt_base = np.zeros(len(load))
    dr_used = np.zeros(len(load))
    for t in range(len(load)):
        net = solar[t] + wind[t] + thermal_mustrun[t] - load[t]
        curt_base[t] = max(net, 0.0)
        if net >= 0:
            charge = min(net, batt_power, (batt_energy - soc) / eta_c)
            soc += eta_c * charge
            curt[t] = max(net - charge, 0.0)
        else:
            deficit = -net
            discharge = min(deficit, batt_power, soc * eta_d)
            soc -= discharge / eta_d
            deficit -= discharge
            dr_used[t] = min(deficit, 0.20 * load[t])

    monthly = pd.Series(curt, index=df.index).resample("ME").sum() / 1000.0
    monthly.index = monthly.index.strftime("%b")
    annual_curtailment_gwh = float(np.sum(curt) / 1000.0)
    annual_curtailment_rate_pct = float(100.0 * np.sum(curt) / np.sum(solar + wind))
    annual_reduction_pct = float(100.0 * (np.sum(curt_base) - np.sum(curt)) / max(np.sum(curt_base), 1e-6))
    peak_month = str(monthly.idxmax())

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].bar(monthly.index, monthly.values, color="forestgreen")
    axes[0].set_title("Monthly Curtailment")
    axes[0].set_ylabel("GWh")
    axes[1].axis("off")
    stats_text = (
        f"Annual curtailment: {annual_curtailment_gwh:.2f} GWh\n"
        f"Curtailment rate: {annual_curtailment_rate_pct:.2f}%\n"
        f"Peak month: {peak_month}\n"
        f"BESS+DR reduction: {annual_reduction_pct:.2f}%\n"
        f"Mean DR activation: {dr_used.mean():.2f} MW"
    )
    axes[1].text(0.05, 0.55, stats_text, fontsize=12, va="center")
    axes[1].set_title("Annual Statistics")
    savefig(fig, "fig09_kyushu_curtailment.png")

    return {
        "monthly_curtailment_gwh": monthly.to_dict(),
        "annual_curtailment_gwh": annual_curtailment_gwh,
        "annual_curtailment_rate_pct": annual_curtailment_rate_pct,
        "peak_curtailment_month": peak_month,
        "curtailment_reduction_with_bess_dr_pct": annual_reduction_pct,
    }


def main():
    power_flow = power_flow_study()
    hem_nr = hem_vs_nr_study(power_flow)
    annual_df = generate_synthetic_year()
    forecast_metrics, _ = forecasting_study(annual_df)
    scenario_opt = scenario_optimization_study(annual_df)
    battery_sched = battery_dr_scheduling_study(annual_df)
    stability = stability_study(power_flow["total_load_mw"])
    kyushu = kyushu_curtailment_study(annual_df)

    results_summary = {
        "power_flow": {
            "nr_iterations": power_flow["nr_iterations"],
            "nr_solve_time_ms": power_flow["nr_time_ms"],
            "dc_solve_time_ms": power_flow["dc_time_ms"],
            "max_voltage_pu": float(power_flow["vm_pu"].max()),
            "min_voltage_pu": float(power_flow["vm_pu"].min()),
            "max_line_loading_pct": float(power_flow["line_loading"].max()),
        },
        "hem_vs_nr": {
            "nr_120pct_load_iterations": hem_nr["High load (120%)"]["nr_iterations"],
            "hem_120pct_load_iterations": "N/A - Padé approximant used",
            "voltage_collapse_loading_pu": stability["voltage_collapse_loading"],
            "scenario_details": {k: {ik: (iv.tolist() if isinstance(iv, np.ndarray) else iv) for ik, iv in v.items()} for k, v in hem_nr.items()},
        },
        "forecast": forecast_metrics,
        "scenario_optimization": scenario_opt,
        "battery_scheduling": battery_sched,
        "frequency_response": {
            "without_inertia_nadir_hz": stability["without_inertia_nadir_hz"],
            "with_inertia_nadir_hz": stability["with_inertia_nadir_hz"],
            "rocof_hz_per_sec": stability["rocof_hz_per_sec"],
            "recovery_time_s": stability["recovery_time_s"],
        },
        "voltage_stability": {
            "voltage_collapse_power_mw": stability["voltage_collapse_power_mw"],
            "stability_margin_pct": stability["stability_margin_pct"],
        },
        "kyushu_curtailment": kyushu,
    }

    with open(DATA_DIR / "results_summary.json", "w", encoding="utf-8") as f:
        json.dump(results_summary, f, indent=2)

    print("=== POWER FLOW RESULTS ===")
    print(f"NR_ITERATIONS: {power_flow['nr_iterations']}")
    print(f"NR_SOLVE_TIME_MS: {power_flow['nr_time_ms']:.3f}")
    print(f"DC_SOLVE_TIME_MS: {power_flow['dc_time_ms']:.3f}")
    print(f"MAX_VOLTAGE_PU: {power_flow['vm_pu'].max():.4f}")
    print(f"MIN_VOLTAGE_PU: {power_flow['vm_pu'].min():.4f}")
    print(f"MAX_LINE_LOADING_PCT: {power_flow['line_loading'].max():.3f}")
    print()

    print("=== HEM vs NR COMPARISON ===")
    print(f"NR_120PCT_LOAD_ITERATIONS: {hem_nr['High load (120%)']['nr_iterations']}")
    print("HEM_120PCT_LOAD_ITERATIONS: N/A - Padé approximant used")
    print(f"VOLTAGE_COLLAPSE_LOADING: {stability['voltage_collapse_loading']:.3f} pu")
    print()

    print("=== FORECAST RESULTS ===")
    print(f"RF_RMSE: {forecast_metrics['rf_rmse']:.3f}")
    print(f"RF_MAE: {forecast_metrics['rf_mae']:.3f}")
    print(f"RF_R2: {forecast_metrics['rf_r2']:.4f}")
    print(f"GBM_RMSE: {forecast_metrics['gbm_rmse']:.3f}")
    print(f"GBM_MAE: {forecast_metrics['gbm_mae']:.3f}")
    print(f"GBM_R2: {forecast_metrics['gbm_r2']:.4f}")
    print(f"COVERAGE_90PCT: {forecast_metrics['coverage_90pct']:.4f}")
    print()

    print("=== SCENARIO OPTIMIZATION ===")
    print(f"EXPECTED_CURTAILMENT_MWH_PER_DAY: {scenario_opt['expected_curtailment_mwh_per_day']:.3f}")
    print(f"LOAD_SHEDDING_PROBABILITY: {scenario_opt['load_shedding_probability']:.4f}")
    print(f"CURTAILMENT_REDUCTION_WITH_STORAGE_PCT: {scenario_opt['curtailment_reduction_with_storage_pct']:.3f}")
    print()

    print("=== BATTERY SCHEDULING ===")
    print(f"BATTERY_UTILIZATION_PCT: {battery_sched['battery_utilization_pct']:.3f}")
    print(f"DR_ACTIVATION_PCT: {battery_sched['dr_activation_pct']:.3f}")
    print(f"CURTAILMENT_WITHOUT_BATTERY_MWH: {battery_sched['curtailment_without_battery_mwh']:.3f}")
    print(f"CURTAILMENT_WITH_BATTERY_MWH: {battery_sched['curtailment_with_battery_mwh']:.3f}")
    print()

    print("=== FREQUENCY RESPONSE ===")
    print(f"WITHOUT_INERTIA_NADIR_HZ: {stability['without_inertia_nadir_hz']:.4f}")
    print(f"WITH_INERTIA_NADIR_HZ: {stability['with_inertia_nadir_hz']:.4f}")
    print(f"ROCOF_HZ_PER_SEC: {stability['rocof_hz_per_sec']:.4f}")
    print(f"RECOVERY_TIME_S: {stability['recovery_time_s']:.3f}")
    print()

    print("=== VOLTAGE STABILITY ===")
    print(f"VOLTAGE_COLLAPSE_POWER_MW: {stability['voltage_collapse_power_mw']:.3f}")
    print(f"STABILITY_MARGIN_PCT: {stability['stability_margin_pct']:.3f}")
    print()

    print("=== KYUSHU CURTAILMENT ===")
    print(f"ANNUAL_CURTAILMENT_GWH: {kyushu['annual_curtailment_gwh']:.3f}")
    print(f"ANNUAL_CURTAILMENT_RATE_PCT: {kyushu['annual_curtailment_rate_pct']:.3f}")
    print(f"PEAK_CURTAILMENT_MONTH: {kyushu['peak_curtailment_month']}")
    print(f"CURTAILMENT_REDUCTION_WITH_BESS_DR_PCT: {kyushu['curtailment_reduction_with_bess_dr_pct']:.3f}")


if __name__ == "__main__":
    main()
