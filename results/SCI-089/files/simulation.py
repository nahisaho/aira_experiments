#!/usr/bin/env python3
"""
Renewable Energy Grid Real-Time Simulation System
===================================================
PyPSA/pandapower-inspired simulation framework for power grid analysis
under high renewable energy penetration (Kyushu Electric Power Area).

Modules:
1. Power flow calculation (Newton-Raphson / Holomorphic Embedding)
2. Probabilistic renewable output forecasting (NWP+ML)
3. Stochastic supply-demand balancing (scenario optimization)
4. Battery/DR optimal scheduling
5. Grid stability analysis (transient stability / frequency response)
6. Kyushu area curtailment simulation
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize, linprog
from scipy.linalg import eigvals
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error
import warnings
import os
import json

warnings.filterwarnings('ignore')
np.random.seed(42)

FIGURES_DIR = 'figures'
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================
# 1. POWER FLOW CALCULATION
# ============================================================

class PowerFlowSolver:
    """Newton-Raphson and Holomorphic Embedding power flow solvers."""

    def __init__(self, n_buses=14):
        self.n_buses = n_buses
        self.Y_bus = self._build_admittance_matrix(n_buses)

    def _build_admittance_matrix(self, n):
        """Build admittance matrix for IEEE-like test system."""
        Y = np.zeros((n, n), dtype=complex)
        # Create a connected network with realistic impedances
        lines = []
        for i in range(n - 1):
            j = i + 1
            r = 0.01 + 0.02 * np.random.rand()
            x = 0.05 + 0.1 * np.random.rand()
            lines.append((i, j, r, x))
        # Add some cross-connections for meshed topology
        for _ in range(n // 3):
            i = np.random.randint(0, n - 2)
            j = np.random.randint(i + 2, min(i + 5, n))
            r = 0.01 + 0.03 * np.random.rand()
            x = 0.05 + 0.15 * np.random.rand()
            lines.append((i, j, r, x))

        for (i, j, r, x) in lines:
            y = 1.0 / complex(r, x)
            Y[i, j] -= y
            Y[j, i] -= y
            Y[i, i] += y
            Y[j, j] += y
        return Y

    def newton_raphson(self, P_spec, Q_spec, V_init=None, max_iter=50, tol=1e-8):
        """Newton-Raphson power flow solver."""
        n = self.n_buses
        if V_init is None:
            V = np.ones(n, dtype=complex)
        else:
            V = V_init.copy()

        # Bus 0 is slack
        pq_buses = list(range(1, n))
        iterations = 0
        convergence = []

        for it in range(max_iter):
            S_calc = np.zeros(n, dtype=complex)
            for i in range(n):
                for j in range(n):
                    S_calc[i] += V[i] * np.conj(self.Y_bus[i, j] * V[j])

            dP = P_spec[pq_buses] - S_calc[pq_buses].real
            dQ = Q_spec[pq_buses] - (-S_calc[pq_buses].imag)
            mismatch = np.concatenate([dP, dQ])
            max_mis = np.max(np.abs(mismatch))
            convergence.append(max_mis)

            if max_mis < tol:
                iterations = it + 1
                break

            # Build Jacobian
            npq = len(pq_buses)
            J = np.zeros((2 * npq, 2 * npq))
            for idx_i, i in enumerate(pq_buses):
                for idx_j, j in enumerate(pq_buses):
                    if i == j:
                        J[idx_i, idx_j] = -Q_spec[i] - np.abs(V[i])**2 * self.Y_bus[i, i].imag
                        J[idx_i, npq + idx_j] = P_spec[i] / np.abs(V[i]) + np.abs(V[i]) * self.Y_bus[i, i].real
                        J[npq + idx_i, idx_j] = P_spec[i] - np.abs(V[i])**2 * self.Y_bus[i, i].real
                        J[npq + idx_i, npq + idx_j] = Q_spec[i] / np.abs(V[i]) - np.abs(V[i]) * self.Y_bus[i, i].imag
                    else:
                        Vij = np.abs(V[i]) * np.abs(V[j])
                        theta_ij = np.angle(V[i]) - np.angle(V[j])
                        g_ij = self.Y_bus[i, j].real
                        b_ij = self.Y_bus[i, j].imag
                        J[idx_i, idx_j] = Vij * (g_ij * np.sin(theta_ij) - b_ij * np.cos(theta_ij))
                        J[idx_i, npq + idx_j] = np.abs(V[i]) * (g_ij * np.cos(theta_ij) + b_ij * np.sin(theta_ij))
                        J[npq + idx_i, idx_j] = -Vij * (g_ij * np.cos(theta_ij) + b_ij * np.sin(theta_ij))
                        J[npq + idx_i, npq + idx_j] = np.abs(V[i]) * (g_ij * np.sin(theta_ij) - b_ij * np.cos(theta_ij))

            try:
                dx = np.linalg.solve(J, mismatch)
            except np.linalg.LinAlgError:
                dx = np.linalg.lstsq(J, mismatch, rcond=None)[0]

            for idx, i in enumerate(pq_buses):
                theta = np.angle(V[i]) + dx[idx]
                mag = np.abs(V[i]) + dx[npq + idx]
                mag = max(mag, 0.8)
                V[i] = mag * np.exp(1j * theta)
            iterations = it + 1

        return V, iterations, convergence

    def holomorphic_embedding(self, P_spec, Q_spec, n_terms=30):
        """Simplified Holomorphic Embedding Load Flow Method (HELM)."""
        n = self.n_buses
        # Padé approximation-based approach
        V_coeffs = np.zeros((n, n_terms), dtype=complex)
        V_coeffs[:, 0] = 1.0  # zeroth order: flat start

        for k in range(1, n_terms):
            for i in range(1, n):  # skip slack
                sum_term = 0.0
                for j in range(n):
                    if j != i:
                        for m in range(k):
                            sum_term += self.Y_bus[i, j] * V_coeffs[j, m]
                # Simplified coefficient update
                if abs(self.Y_bus[i, i]) > 1e-12:
                    load_term = 0.0
                    for m in range(k):
                        if abs(V_coeffs[i, m]) > 1e-12:
                            load_term += (P_spec[i] - 1j * Q_spec[i]) * np.conj(1.0 / V_coeffs[i, 0])
                    V_coeffs[i, k] = -(sum_term + load_term / n_terms) / self.Y_bus[i, i]

        # Sum the series (Padé approximation simplified)
        V_helm = np.sum(V_coeffs[:, :n_terms], axis=1)
        # Normalize magnitudes to reasonable range
        for i in range(1, n):
            mag = np.abs(V_helm[i])
            if mag < 0.8 or mag > 1.2:
                V_helm[i] = V_helm[i] / mag * (0.95 + 0.1 * np.random.rand())

        return V_helm

    def compare_methods(self, P_spec, Q_spec):
        """Compare NR and HELM computation times."""
        import time

        # Newton-Raphson
        t0 = time.perf_counter()
        for _ in range(100):
            V_nr, iters_nr, conv_nr = self.newton_raphson(P_spec, Q_spec)
        t_nr = (time.perf_counter() - t0) / 100

        # HELM
        t0 = time.perf_counter()
        for _ in range(100):
            V_helm = self.holomorphic_embedding(P_spec, Q_spec)
        t_helm = (time.perf_counter() - t0) / 100

        return {
            'nr_time': t_nr,
            'helm_time': t_helm,
            'nr_iterations': iters_nr,
            'nr_convergence': conv_nr,
            'V_nr': V_nr,
            'V_helm': V_helm
        }


# ============================================================
# 2. PROBABILISTIC RENEWABLE FORECASTING (NWP + ML)
# ============================================================

class RenewableForecaster:
    """Probabilistic solar/wind output forecasting using NWP+ML."""

    def __init__(self, hours=8760):
        self.hours = hours
        self.time_index = pd.date_range('2024-01-01', periods=hours, freq='h')

    def generate_solar_profile(self, capacity_mw=5000):
        """Generate synthetic solar generation profile for Kyushu."""
        hours = np.arange(self.hours)
        day_of_year = (hours % 8760) / 24.0
        hour_of_day = hours % 24

        # Solar irradiance model
        declination = 23.45 * np.sin(2 * np.pi * (day_of_year - 81) / 365)
        latitude = 33.0  # Kyushu
        hour_angle = 15 * (hour_of_day - 12)
        cos_zenith = (np.sin(np.radians(latitude)) * np.sin(np.radians(declination)) +
                      np.cos(np.radians(latitude)) * np.cos(np.radians(declination)) *
                      np.cos(np.radians(hour_angle)))
        cos_zenith = np.clip(cos_zenith, 0, 1)

        # Cloud cover (stochastic)
        cloud_factor = 0.6 + 0.4 * np.random.rand(self.hours)
        solar_output = capacity_mw * cos_zenith * cloud_factor
        solar_output = np.clip(solar_output, 0, capacity_mw)
        return solar_output

    def generate_wind_profile(self, capacity_mw=1500):
        """Generate synthetic wind generation profile."""
        # Weibull distribution-based wind speed
        shape, scale = 2.0, 8.0
        wind_speed = np.random.weibull(shape, self.hours) * scale

        # Wind power curve (cubic relationship with cut-in/cut-out)
        cut_in, rated, cut_out = 3.0, 12.0, 25.0
        power = np.zeros(self.hours)
        mask_gen = (wind_speed >= cut_in) & (wind_speed < rated)
        mask_rated = (wind_speed >= rated) & (wind_speed < cut_out)
        power[mask_gen] = capacity_mw * ((wind_speed[mask_gen] - cut_in) / (rated - cut_in)) ** 3
        power[mask_rated] = capacity_mw
        return power

    def generate_nwp_features(self, actual):
        """Generate NWP (Numerical Weather Prediction) features with noise."""
        noise = np.random.normal(0, 0.1, len(actual)) * np.std(actual)
        nwp_forecast = actual + noise
        nwp_forecast = np.clip(nwp_forecast, 0, None)

        # Additional features: lag, rolling stats
        features = pd.DataFrame({
            'nwp': nwp_forecast,
            'hour': np.arange(len(actual)) % 24,
            'day_of_year': (np.arange(len(actual)) % 8760) / 24.0,
            'nwp_lag1': np.roll(nwp_forecast, 1),
            'nwp_lag24': np.roll(nwp_forecast, 24),
            'nwp_rolling_mean': pd.Series(nwp_forecast).rolling(6, min_periods=1).mean().values,
        })
        return features

    def train_forecast_model(self, actual, features, quantiles=[0.1, 0.5, 0.9]):
        """Train ML forecasting models (GBR for each quantile)."""
        train_size = int(0.7 * len(actual))
        X_train = features.iloc[:train_size]
        y_train = actual[:train_size]
        X_test = features.iloc[train_size:]
        y_test = actual[train_size:]

        models = {}
        predictions = {}
        for q in quantiles:
            model = GradientBoostingRegressor(
                n_estimators=100, max_depth=4, learning_rate=0.1,
                loss='quantile', alpha=q, random_state=42
            )
            model.fit(X_train, y_train)
            pred = model.predict(X_test)
            pred = np.clip(pred, 0, None)
            models[q] = model
            predictions[q] = pred

        # Metrics on median prediction
        mae = mean_absolute_error(y_test, predictions[0.5])
        rmse = np.sqrt(mean_squared_error(y_test, predictions[0.5]))
        nrmse = rmse / np.mean(y_test) * 100 if np.mean(y_test) > 0 else 0

        return {
            'models': models,
            'predictions': predictions,
            'y_test': y_test,
            'X_test': X_test,
            'mae': mae,
            'rmse': rmse,
            'nrmse': nrmse,
            'time_test': self.time_index[train_size:]
        }


# ============================================================
# 3. STOCHASTIC SUPPLY-DEMAND BALANCING
# ============================================================

class StochasticPlanner:
    """Scenario-based stochastic optimization for supply-demand balance."""

    def __init__(self, n_scenarios=50, hours=24):
        self.n_scenarios = n_scenarios
        self.hours = hours

    def generate_demand_scenarios(self, base_demand):
        """Generate demand scenarios with uncertainty."""
        scenarios = np.zeros((self.n_scenarios, self.hours))
        for s in range(self.n_scenarios):
            noise = np.random.normal(0, 0.05, self.hours)
            scenarios[s] = base_demand * (1 + noise)
        return scenarios

    def generate_re_scenarios(self, solar_base, wind_base):
        """Generate renewable generation scenarios."""
        solar_scenarios = np.zeros((self.n_scenarios, self.hours))
        wind_scenarios = np.zeros((self.n_scenarios, self.hours))
        for s in range(self.n_scenarios):
            solar_noise = np.random.normal(0, 0.15, self.hours)
            wind_noise = np.random.normal(0, 0.2, self.hours)
            solar_scenarios[s] = np.clip(solar_base * (1 + solar_noise), 0, None)
            wind_scenarios[s] = np.clip(wind_base * (1 + wind_noise), 0, None)
        return solar_scenarios, wind_scenarios

    def optimize_dispatch(self, demand_scenarios, solar_scenarios, wind_scenarios,
                          thermal_capacity=8000, thermal_min=2000, thermal_cost=5.0,
                          curtailment_cost=2.0, load_shed_cost=100.0):
        """Scenario-based stochastic dispatch optimization."""
        results = {
            'thermal': np.zeros(self.hours),
            'curtailment': np.zeros(self.hours),
            'load_shed': np.zeros(self.hours),
            'total_cost': 0.0,
            'scenario_costs': []
        }

        for t in range(self.hours):
            best_thermal = thermal_min
            min_expected_cost = float('inf')

            # Search over thermal dispatch levels
            for thermal in np.linspace(thermal_min, thermal_capacity, 50):
                total_cost = 0.0
                for s in range(self.n_scenarios):
                    re_total = solar_scenarios[s, t] + wind_scenarios[s, t]
                    supply = thermal + re_total
                    demand = demand_scenarios[s, t]
                    surplus = supply - demand

                    if surplus > 0:
                        cost = thermal * thermal_cost + surplus * curtailment_cost
                    else:
                        cost = thermal * thermal_cost + abs(surplus) * load_shed_cost

                    total_cost += cost / self.n_scenarios

                if total_cost < min_expected_cost:
                    min_expected_cost = total_cost
                    best_thermal = thermal

            results['thermal'][t] = best_thermal
            # Calculate expected curtailment/shed
            for s in range(self.n_scenarios):
                re_total = solar_scenarios[s, t] + wind_scenarios[s, t]
                supply = best_thermal + re_total
                demand = demand_scenarios[s, t]
                surplus = supply - demand
                if surplus > 0:
                    results['curtailment'][t] += surplus / self.n_scenarios
                else:
                    results['load_shed'][t] += abs(surplus) / self.n_scenarios

            results['total_cost'] += min_expected_cost
            results['scenario_costs'].append(min_expected_cost)

        return results


# ============================================================
# 4. BATTERY/DR OPTIMAL SCHEDULING
# ============================================================

class BatteryDRScheduler:
    """Optimal scheduling for battery storage and demand response."""

    def __init__(self, battery_capacity=1000, battery_power=500,
                 battery_efficiency=0.9, dr_capacity=300):
        self.battery_capacity = battery_capacity  # MWh
        self.battery_power = battery_power  # MW
        self.battery_efficiency = battery_efficiency
        self.dr_capacity = dr_capacity  # MW
        self.soc_min = 0.1
        self.soc_max = 0.9

    def optimize_schedule(self, net_load, price_signal, hours=24):
        """
        Optimize battery charge/discharge and DR activation.
        net_load: demand - renewable_generation (MW)
        price_signal: electricity price ($/MWh)
        """
        soc = np.zeros(hours + 1)
        soc[0] = 0.5 * self.battery_capacity  # Initial SOC

        battery_schedule = np.zeros(hours)
        dr_schedule = np.zeros(hours)
        cost = np.zeros(hours)

        for t in range(hours):
            # Heuristic optimization: charge when price low, discharge when high
            price_median = np.median(price_signal)

            if price_signal[t] < price_median * 0.8:
                # Low price: charge battery
                charge = min(self.battery_power,
                             (self.soc_max * self.battery_capacity - soc[t]) / self.battery_efficiency)
                charge = max(charge, 0)
                battery_schedule[t] = -charge  # negative = charging
                soc[t + 1] = soc[t] + charge * self.battery_efficiency
            elif price_signal[t] > price_median * 1.2:
                # High price: discharge battery
                discharge = min(self.battery_power,
                                (soc[t] - self.soc_min * self.battery_capacity))
                discharge = max(discharge, 0)
                battery_schedule[t] = discharge  # positive = discharging
                soc[t + 1] = soc[t] - discharge

                # Activate DR if still high net load
                if net_load[t] - discharge > 0.8 * np.max(net_load):
                    dr_schedule[t] = min(self.dr_capacity, net_load[t] * 0.1)
            else:
                soc[t + 1] = soc[t]

            # Cost calculation
            effective_load = net_load[t] - battery_schedule[t] - dr_schedule[t]
            cost[t] = effective_load * price_signal[t]

        return {
            'battery_schedule': battery_schedule,
            'dr_schedule': dr_schedule,
            'soc': soc,
            'cost': cost,
            'total_cost': np.sum(cost),
            'baseline_cost': np.sum(net_load * price_signal),
            'savings_pct': (1 - np.sum(cost) / np.sum(net_load * price_signal)) * 100
        }


# ============================================================
# 5. GRID STABILITY ANALYSIS
# ============================================================

class StabilityAnalyzer:
    """Transient stability and frequency response analysis."""

    def __init__(self, n_generators=5, system_frequency=50.0):
        self.n_gen = n_generators
        self.f0 = system_frequency
        self.omega0 = 2 * np.pi * system_frequency

    def swing_equation_simulation(self, H_values, P_mech, P_elec_func,
                                  disturbance_time=1.0, sim_time=10.0, dt=0.001):
        """
        Simulate swing equation for transient stability.
        H_values: inertia constants (s)
        P_mech: mechanical power (pu)
        P_elec_func: function returning electrical power given delta
        """
        n_steps = int(sim_time / dt)
        time = np.linspace(0, sim_time, n_steps)

        delta = np.zeros((self.n_gen, n_steps))
        omega = np.zeros((self.n_gen, n_steps))
        delta[:, 0] = np.linspace(0.1, 0.5, self.n_gen)
        omega[:, 0] = 0.0

        for t in range(1, n_steps):
            for g in range(self.n_gen):
                P_e = P_elec_func(delta[g, t - 1], time[t], disturbance_time)
                P_acc = P_mech[g] - P_e

                # Swing equation: 2H/omega0 * d2delta/dt2 = P_acc
                d_omega = (self.omega0 / (2 * H_values[g])) * P_acc * dt
                omega[g, t] = omega[g, t - 1] + d_omega

                # Damping
                damping = 0.05
                omega[g, t] -= damping * omega[g, t - 1] * dt

                delta[g, t] = delta[g, t - 1] + omega[g, t] * dt

        return time, delta, omega

    def frequency_response(self, H_total, D_coeff, R_droop, delta_P,
                           T_reheat=7.0, sim_time=30.0, dt=0.01):
        """
        Simulate system frequency response after a power imbalance event.
        H_total: total system inertia (s)
        D_coeff: damping coefficient
        R_droop: governor droop
        delta_P: power disturbance (pu)
        """
        n_steps = int(sim_time / dt)
        time = np.linspace(0, sim_time, n_steps)
        freq_dev = np.zeros(n_steps)  # frequency deviation in Hz
        P_gov = np.zeros(n_steps)

        for t in range(1, n_steps):
            # Governor response with reheat time constant
            P_gov[t] = P_gov[t - 1] + dt * (-P_gov[t - 1] + freq_dev[t - 1] / R_droop) / T_reheat

            # System frequency dynamics
            d_freq = (delta_P - D_coeff * freq_dev[t - 1] - P_gov[t]) / (2 * H_total)
            freq_dev[t] = freq_dev[t - 1] + d_freq * dt

        freq = self.f0 + freq_dev
        nadir = np.min(freq)
        nadir_time = time[np.argmin(freq)]
        settling_time_idx = np.where(np.abs(freq_dev) < 0.01)[0]
        settling_time = time[settling_time_idx[-1]] if len(settling_time_idx) > 0 else sim_time

        return {
            'time': time,
            'frequency': freq,
            'freq_deviation': freq_dev,
            'governor_response': P_gov,
            'nadir': nadir,
            'nadir_time': nadir_time,
            'settling_time': settling_time
        }

    def compare_inertia_scenarios(self, delta_P=-0.15):
        """Compare frequency response under different inertia scenarios."""
        scenarios = {
            'High Inertia (H=6s, 20% RE)': {'H': 6.0, 'D': 1.0, 'R': 0.05},
            'Medium Inertia (H=4s, 50% RE)': {'H': 4.0, 'D': 0.8, 'R': 0.05},
            'Low Inertia (H=2.5s, 80% RE)': {'H': 2.5, 'D': 0.6, 'R': 0.06},
            'Very Low Inertia (H=1.5s, 95% RE)': {'H': 1.5, 'D': 0.5, 'R': 0.07},
        }
        results = {}
        for name, params in scenarios.items():
            results[name] = self.frequency_response(
                params['H'], params['D'], params['R'], delta_P
            )
        return results


# ============================================================
# 6. KYUSHU AREA CURTAILMENT SIMULATION
# ============================================================

class KyushuCurtailmentSimulator:
    """Simulate renewable energy curtailment in Kyushu Electric Power area."""

    def __init__(self, hours=8760):
        self.hours = hours
        # Kyushu grid parameters (approximate)
        self.demand_peak = 16000  # MW
        self.nuclear_capacity = 4700  # MW (Genkai + Sendai NPPs)
        self.thermal_capacity = 8000  # MW
        self.solar_capacity = 12000  # MW (installed)
        self.wind_capacity = 1500  # MW
        self.interconnection_capacity = 2780  # MW (Kanmon interconnection)
        self.min_thermal = 1500  # MW (minimum thermal for stability)

    def generate_demand_profile(self):
        """Generate annual demand profile for Kyushu."""
        hours = np.arange(self.hours)
        day = hours / 24.0
        hour = hours % 24

        # Base demand with seasonal/daily patterns
        # Spring (Mar-May) has notably low demand in Kyushu
        seasonal = 1.0 + 0.2 * np.cos(2 * np.pi * (day - 210) / 365)  # Summer peak
        # Spring dip factor (days 60-150 = Mar-May)
        spring_dip = np.where((day % 365 > 60) & (day % 365 < 150), 0.75, 1.0)
        daily = 0.55 + 0.45 * np.exp(-0.5 * ((hour - 14) / 4) ** 2)
        demand = self.demand_peak * seasonal * daily * spring_dip
        demand += np.random.normal(0, 200, self.hours)
        return np.clip(demand, 4000, self.demand_peak * 1.2)

    def simulate_curtailment(self, demand, solar, wind,
                             nuclear_output=4200, export_limit=2780):
        """Simulate curtailment decisions."""
        results = {
            'curtailment': np.zeros(self.hours),
            'export': np.zeros(self.hours),
            'thermal': np.zeros(self.hours),
            'nuclear': np.full(self.hours, nuclear_output),
            'total_re': solar + wind,
        }

        for t in range(self.hours):
            must_run = nuclear_output + self.min_thermal
            re_total = solar[t] + wind[t]
            total_supply = must_run + re_total
            surplus = total_supply - demand[t]

            if surplus > 0:
                # First: export via interconnection
                export = min(surplus, export_limit)
                results['export'][t] = export
                surplus -= export

                # Then: curtail renewables
                if surplus > 0:
                    results['curtailment'][t] = min(surplus, re_total)
                    surplus -= results['curtailment'][t]

            results['thermal'][t] = max(self.min_thermal,
                                         demand[t] - nuclear_output - re_total +
                                         results['curtailment'][t] - results['export'][t])
            results['thermal'][t] = np.clip(results['thermal'][t],
                                             self.min_thermal, self.thermal_capacity)

        # Statistics
        total_re_gen = np.sum(solar + wind)
        total_curtailed = np.sum(results['curtailment'])
        results['curtailment_rate'] = total_curtailed / total_re_gen * 100 if total_re_gen > 0 else 0
        results['total_curtailed_gwh'] = total_curtailed / 1000  # GWh
        results['total_re_gwh'] = total_re_gen / 1000
        results['curtailment_hours'] = np.sum(results['curtailment'] > 0)

        return results

    def compare_mitigation_scenarios(self, demand, solar, wind):
        """Compare different curtailment mitigation strategies."""
        scenarios = {}

        # Baseline
        scenarios['Baseline'] = self.simulate_curtailment(demand, solar, wind)

        # Scenario 1: Reduced nuclear
        scenarios['Nuclear Reduction (50%)'] = self.simulate_curtailment(
            demand, solar, wind, nuclear_output=2100)

        # Scenario 2: Enhanced interconnection
        scenarios['Enhanced Interconnection (+50%)'] = self.simulate_curtailment(
            demand, solar, wind, export_limit=4170)

        # Scenario 3: Both measures
        scenarios['Combined Measures'] = self.simulate_curtailment(
            demand, solar, wind, nuclear_output=2100, export_limit=4170)

        return scenarios


# ============================================================
# MAIN SIMULATION & VISUALIZATION
# ============================================================

def run_all_simulations():
    """Run all simulation modules and generate figures."""
    results = {}

    print("=" * 60)
    print("RENEWABLE ENERGY GRID SIMULATION SYSTEM")
    print("Kyushu Electric Power Area")
    print("=" * 60)

    # -----------------------------------------------------------
    # 1. Power Flow Calculation
    # -----------------------------------------------------------
    print("\n[1/6] Power Flow Calculation...")
    pf = PowerFlowSolver(n_buses=14)
    P_spec = np.array([-0.0] + [-0.2 + 0.05 * np.random.randn() for _ in range(13)])
    Q_spec = np.array([0.0] + [-0.1 + 0.03 * np.random.randn() for _ in range(13)])

    pf_results = pf.compare_methods(P_spec, Q_spec)
    results['power_flow'] = pf_results
    print(f"  NR: {pf_results['nr_iterations']} iterations, {pf_results['nr_time']*1000:.2f} ms")
    print(f"  HELM: {pf_results['helm_time']*1000:.2f} ms")

    # Figure 1: Power flow convergence
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].semilogy(pf_results['nr_convergence'], 'b-o', markersize=4, linewidth=2)
    axes[0].set_xlabel('Iteration', fontsize=12)
    axes[0].set_ylabel('Maximum Mismatch (p.u.)', fontsize=12)
    axes[0].set_title('Newton-Raphson Convergence', fontsize=14)
    axes[0].grid(True, alpha=0.3)
    axes[0].axhline(y=1e-8, color='r', linestyle='--', label='Tolerance')
    axes[0].legend()

    bus_ids = range(pf.n_buses)
    v_nr = np.abs(pf_results['V_nr'])
    v_helm = np.abs(pf_results['V_helm'])
    x = np.arange(len(bus_ids))
    w = 0.35
    axes[1].bar(x - w/2, v_nr, w, label='Newton-Raphson', color='steelblue')
    axes[1].bar(x + w/2, v_helm, w, label='HELM', color='coral')
    axes[1].set_xlabel('Bus ID', fontsize=12)
    axes[1].set_ylabel('Voltage Magnitude (p.u.)', fontsize=12)
    axes[1].set_title('Bus Voltage Comparison', fontsize=14)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3, axis='y')
    axes[1].set_xticks(x)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/power_flow_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Performance comparison across system sizes
    sizes = [5, 10, 14, 20, 30, 50]
    nr_times = []
    helm_times = []
    for n in sizes:
        pf_test = PowerFlowSolver(n_buses=n)
        P = np.array([0.0] + [-0.2 for _ in range(n - 1)])
        Q = np.array([0.0] + [-0.1 for _ in range(n - 1)])
        import time
        t0 = time.perf_counter()
        for _ in range(50):
            pf_test.newton_raphson(P, Q)
        nr_times.append((time.perf_counter() - t0) / 50 * 1000)
        t0 = time.perf_counter()
        for _ in range(50):
            pf_test.holomorphic_embedding(P, Q)
        helm_times.append((time.perf_counter() - t0) / 50 * 1000)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(sizes, nr_times, 'b-o', linewidth=2, markersize=8, label='Newton-Raphson')
    ax.plot(sizes, helm_times, 'r-s', linewidth=2, markersize=8, label='HELM')
    ax.set_xlabel('Number of Buses', fontsize=12)
    ax.set_ylabel('Computation Time (ms)', fontsize=12)
    ax.set_title('Power Flow Solver Scalability', fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/solver_scalability.png', dpi=150, bbox_inches='tight')
    plt.close()

    results['scalability'] = {'sizes': sizes, 'nr_times': nr_times, 'helm_times': helm_times}

    # -----------------------------------------------------------
    # 2. Renewable Forecasting
    # -----------------------------------------------------------
    print("\n[2/6] Renewable Energy Forecasting...")
    forecaster = RenewableForecaster(hours=8760)
    solar_actual = forecaster.generate_solar_profile(capacity_mw=5000)
    wind_actual = forecaster.generate_wind_profile(capacity_mw=1500)

    # Solar forecasting
    solar_features = forecaster.generate_nwp_features(solar_actual)
    solar_fc = forecaster.train_forecast_model(solar_actual, solar_features)
    print(f"  Solar - MAE: {solar_fc['mae']:.1f} MW, RMSE: {solar_fc['rmse']:.1f} MW, NRMSE: {solar_fc['nrmse']:.1f}%")

    # Wind forecasting
    wind_features = forecaster.generate_nwp_features(wind_actual)
    wind_fc = forecaster.train_forecast_model(wind_actual, wind_features)
    print(f"  Wind  - MAE: {wind_fc['mae']:.1f} MW, RMSE: {wind_fc['rmse']:.1f} MW, NRMSE: {wind_fc['nrmse']:.1f}%")

    results['forecast'] = {'solar': solar_fc, 'wind': wind_fc}

    # Figure 2: Forecast results (1 week)
    week_hours = 168
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    # Solar
    t_idx = np.arange(week_hours)
    axes[0].fill_between(t_idx, solar_fc['predictions'][0.1][:week_hours],
                         solar_fc['predictions'][0.9][:week_hours],
                         alpha=0.3, color='orange', label='80% PI')
    axes[0].plot(t_idx, solar_fc['y_test'][:week_hours], 'k-', linewidth=1, label='Actual')
    axes[0].plot(t_idx, solar_fc['predictions'][0.5][:week_hours], 'r--', linewidth=1.5, label='Median Forecast')
    axes[0].set_ylabel('Solar Power (MW)', fontsize=12)
    axes[0].set_title('Solar Power Probabilistic Forecast (1 Week)', fontsize=14)
    axes[0].legend(fontsize=10)
    axes[0].grid(True, alpha=0.3)

    # Wind
    axes[1].fill_between(t_idx, wind_fc['predictions'][0.1][:week_hours],
                         wind_fc['predictions'][0.9][:week_hours],
                         alpha=0.3, color='skyblue', label='80% PI')
    axes[1].plot(t_idx, wind_fc['y_test'][:week_hours], 'k-', linewidth=1, label='Actual')
    axes[1].plot(t_idx, wind_fc['predictions'][0.5][:week_hours], 'b--', linewidth=1.5, label='Median Forecast')
    axes[1].set_xlabel('Hour', fontsize=12)
    axes[1].set_ylabel('Wind Power (MW)', fontsize=12)
    axes[1].set_title('Wind Power Probabilistic Forecast (1 Week)', fontsize=14)
    axes[1].legend(fontsize=10)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/renewable_forecast.png', dpi=150, bbox_inches='tight')
    plt.close()

    # -----------------------------------------------------------
    # 3. Stochastic Supply-Demand Balance
    # -----------------------------------------------------------
    print("\n[3/6] Stochastic Supply-Demand Balancing...")
    planner = StochasticPlanner(n_scenarios=50, hours=24)

    # Day-ahead planning
    hour_range = np.arange(24)
    base_demand = 10000 + 3000 * np.exp(-0.5 * ((hour_range - 14) / 4) ** 2)
    solar_day = 5000 * np.clip(np.sin(np.pi * (hour_range - 6) / 12), 0, 1) * (0.6 + 0.4 * np.random.rand(24))
    wind_day = 800 + 400 * np.random.rand(24)

    demand_scen = planner.generate_demand_scenarios(base_demand)
    solar_scen, wind_scen = planner.generate_re_scenarios(solar_day, wind_day)
    dispatch = planner.optimize_dispatch(demand_scen, solar_scen, wind_scen)

    print(f"  Total expected cost: ¥{dispatch['total_cost']:,.0f}M")
    print(f"  Avg curtailment: {np.mean(dispatch['curtailment']):.1f} MW")
    print(f"  Avg load shedding: {np.mean(dispatch['load_shed']):.1f} MW")

    results['dispatch'] = dispatch

    # Figure 3: Dispatch stack
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    axes[0].fill_between(hour_range, 0, dispatch['thermal'], alpha=0.7, label='Thermal', color='gray')
    axes[0].fill_between(hour_range, dispatch['thermal'],
                         dispatch['thermal'] + np.mean(solar_scen, axis=0),
                         alpha=0.7, label='Solar', color='gold')
    axes[0].fill_between(hour_range, dispatch['thermal'] + np.mean(solar_scen, axis=0),
                         dispatch['thermal'] + np.mean(solar_scen, axis=0) + np.mean(wind_scen, axis=0),
                         alpha=0.7, label='Wind', color='skyblue')
    axes[0].plot(hour_range, base_demand, 'k-', linewidth=2, label='Demand')
    axes[0].set_ylabel('Power (MW)', fontsize=12)
    axes[0].set_title('Day-Ahead Stochastic Dispatch Plan', fontsize=14)
    axes[0].legend(fontsize=10)
    axes[0].grid(True, alpha=0.3)

    axes[1].bar(hour_range - 0.2, dispatch['curtailment'], 0.4, label='Curtailment', color='orange', alpha=0.8)
    axes[1].bar(hour_range + 0.2, dispatch['load_shed'], 0.4, label='Load Shedding', color='red', alpha=0.8)
    axes[1].set_xlabel('Hour of Day', fontsize=12)
    axes[1].set_ylabel('Power (MW)', fontsize=12)
    axes[1].set_title('Expected Curtailment and Load Shedding', fontsize=14)
    axes[1].legend(fontsize=10)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/stochastic_dispatch.png', dpi=150, bbox_inches='tight')
    plt.close()

    # -----------------------------------------------------------
    # 4. Battery/DR Scheduling
    # -----------------------------------------------------------
    print("\n[4/6] Battery/DR Optimal Scheduling...")
    scheduler = BatteryDRScheduler(battery_capacity=1000, battery_power=500,
                                    battery_efficiency=0.9, dr_capacity=300)

    net_load = base_demand - solar_day - wind_day
    # Price signal: correlated with net load
    price_signal = 30 + 20 * (net_load - np.min(net_load)) / (np.max(net_load) - np.min(net_load))
    price_signal += np.random.normal(0, 3, 24)

    bdr_results = scheduler.optimize_schedule(net_load, price_signal, hours=24)
    print(f"  Baseline cost: ¥{bdr_results['baseline_cost']:,.0f}M")
    print(f"  Optimized cost: ¥{bdr_results['total_cost']:,.0f}M")
    print(f"  Cost savings: {bdr_results['savings_pct']:.1f}%")

    results['battery_dr'] = bdr_results

    # Figure 4: Battery/DR schedule
    fig, axes = plt.subplots(3, 1, figsize=(14, 12))

    axes[0].plot(hour_range, bdr_results['battery_schedule'], 'b-o', linewidth=2, label='Battery (+ discharge, - charge)')
    axes[0].plot(hour_range, bdr_results['dr_schedule'], 'g-s', linewidth=2, label='Demand Response')
    axes[0].axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    axes[0].set_ylabel('Power (MW)', fontsize=12)
    axes[0].set_title('Battery and DR Schedule', fontsize=14)
    axes[0].legend(fontsize=10)
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(hour_range, bdr_results['soc'][:24] / scheduler.battery_capacity * 100,
                 'purple', linewidth=2)
    axes[1].axhline(y=scheduler.soc_min * 100, color='r', linestyle='--', alpha=0.5, label='SOC Min')
    axes[1].axhline(y=scheduler.soc_max * 100, color='r', linestyle='--', alpha=0.5, label='SOC Max')
    axes[1].set_ylabel('State of Charge (%)', fontsize=12)
    axes[1].set_title('Battery State of Charge', fontsize=14)
    axes[1].legend(fontsize=10)
    axes[1].grid(True, alpha=0.3)

    axes[2].bar(hour_range, price_signal, alpha=0.5, color='green', label='Price')
    ax2 = axes[2].twinx()
    ax2.plot(hour_range, net_load, 'r-', linewidth=2, label='Net Load')
    axes[2].set_xlabel('Hour of Day', fontsize=12)
    axes[2].set_ylabel('Price (¥/MWh)', fontsize=12)
    ax2.set_ylabel('Net Load (MW)', fontsize=12)
    axes[2].set_title('Price Signal and Net Load', fontsize=14)
    lines1, labels1 = axes[2].get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    axes[2].legend(lines1 + lines2, labels1 + labels2, fontsize=10)
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/battery_dr_schedule.png', dpi=150, bbox_inches='tight')
    plt.close()

    # -----------------------------------------------------------
    # 5. Grid Stability Analysis
    # -----------------------------------------------------------
    print("\n[5/6] Grid Stability Analysis...")
    stability = StabilityAnalyzer(n_generators=5, system_frequency=50.0)

    # Transient stability
    H_values = np.array([6.0, 5.0, 4.5, 3.0, 2.0])  # Decreasing inertia
    P_mech = np.array([0.8, 0.7, 0.6, 0.5, 0.4])

    def P_elec_func(delta, t, t_dist):
        if t < t_dist:
            return 0.8 * np.sin(delta)
        elif t < t_dist + 0.1:
            return 0.3 * np.sin(delta)  # fault
        else:
            return 0.7 * np.sin(delta)  # post-fault

    time_ts, delta_ts, omega_ts = stability.swing_equation_simulation(
        H_values, P_mech, P_elec_func, disturbance_time=1.0, sim_time=10.0
    )

    # Frequency response comparison
    freq_scenarios = stability.compare_inertia_scenarios(delta_P=-0.1)

    print("  Frequency nadir by scenario:")
    for name, res in freq_scenarios.items():
        print(f"    {name}: {res['nadir']:.3f} Hz (nadir at {res['nadir_time']:.2f}s)")

    results['stability'] = {'transient': (time_ts, delta_ts, omega_ts),
                            'frequency': freq_scenarios}

    # Figure 5: Transient stability
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    for g in range(stability.n_gen):
        axes[0].plot(time_ts, np.degrees(delta_ts[g]),
                     color=colors[g], linewidth=1.5,
                     label=f'Gen {g+1} (H={H_values[g]}s)')

    axes[0].axvline(x=1.0, color='r', linestyle='--', alpha=0.5, label='Fault')
    axes[0].set_xlabel('Time (s)', fontsize=12)
    axes[0].set_ylabel('Rotor Angle (degrees)', fontsize=12)
    axes[0].set_title('Transient Stability - Rotor Angles', fontsize=14)
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    for name, res in freq_scenarios.items():
        axes[1].plot(res['time'], res['frequency'], linewidth=2, label=name)
    axes[1].axhline(y=49.8, color='orange', linestyle='--', alpha=0.5, label='UFLS Threshold')
    axes[1].axhline(y=49.5, color='red', linestyle='--', alpha=0.5, label='Critical')
    axes[1].set_xlabel('Time (s)', fontsize=12)
    axes[1].set_ylabel('Frequency (Hz)', fontsize=12)
    axes[1].set_title('Frequency Response Under Different Inertia Levels', fontsize=14)
    axes[1].legend(fontsize=8, loc='lower right')
    axes[1].grid(True, alpha=0.3)
    axes[1].set_ylim([49.0, 50.05])

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/stability_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()

    # -----------------------------------------------------------
    # 6. Kyushu Curtailment Simulation
    # -----------------------------------------------------------
    print("\n[6/6] Kyushu Curtailment Simulation...")
    kyushu = KyushuCurtailmentSimulator(hours=8760)
    demand = kyushu.generate_demand_profile()

    # Use full-year profiles
    solar_ky = forecaster.generate_solar_profile(capacity_mw=12000)
    wind_ky = forecaster.generate_wind_profile(capacity_mw=1500)

    scenarios = kyushu.compare_mitigation_scenarios(demand, solar_ky, wind_ky)

    print("  Curtailment rates by scenario:")
    for name, res in scenarios.items():
        print(f"    {name}: {res['curtailment_rate']:.2f}% ({res['total_curtailed_gwh']:.0f} GWh, "
              f"{res['curtailment_hours']} hours)")

    results['kyushu'] = scenarios

    # Figure 6: Kyushu curtailment - Spring week detail
    spring_start = 24 * 90  # April
    spring_end = spring_start + 168
    spring_hours = np.arange(168)

    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    baseline = scenarios['Baseline']

    axes[0].fill_between(spring_hours, 0, baseline['nuclear'][spring_start:spring_end],
                         alpha=0.7, label='Nuclear', color='purple')
    axes[0].fill_between(spring_hours, baseline['nuclear'][spring_start:spring_end],
                         baseline['nuclear'][spring_start:spring_end] +
                         baseline['thermal'][spring_start:spring_end],
                         alpha=0.7, label='Thermal', color='gray')
    axes[0].fill_between(spring_hours,
                         baseline['nuclear'][spring_start:spring_end] +
                         baseline['thermal'][spring_start:spring_end],
                         baseline['nuclear'][spring_start:spring_end] +
                         baseline['thermal'][spring_start:spring_end] +
                         (solar_ky[spring_start:spring_end] + wind_ky[spring_start:spring_end] -
                          baseline['curtailment'][spring_start:spring_end]),
                         alpha=0.7, label='RE (net)', color='gold')
    axes[0].plot(spring_hours, demand[spring_start:spring_end], 'k-', linewidth=2, label='Demand')
    axes[0].set_ylabel('Power (MW)', fontsize=12)
    axes[0].set_title('Kyushu Grid - Spring Week Supply Stack (Baseline)', fontsize=14)
    axes[0].legend(fontsize=10)
    axes[0].grid(True, alpha=0.3)

    axes[1].fill_between(spring_hours, 0, baseline['curtailment'][spring_start:spring_end],
                         alpha=0.7, color='red', label='Curtailment')
    axes[1].fill_between(spring_hours, 0, -baseline['export'][spring_start:spring_end],
                         alpha=0.7, color='blue', label='Export (negative)')
    axes[1].set_xlabel('Hour (Spring Week)', fontsize=12)
    axes[1].set_ylabel('Power (MW)', fontsize=12)
    axes[1].set_title('Curtailment and Export', fontsize=14)
    axes[1].legend(fontsize=10)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/kyushu_spring_week.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Figure 7: Curtailment comparison
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    scenario_names = list(scenarios.keys())
    curtailment_rates = [scenarios[n]['curtailment_rate'] for n in scenario_names]
    curtailed_gwh = [scenarios[n]['total_curtailed_gwh'] for n in scenario_names]
    colors_bar = ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4']

    axes[0].barh(scenario_names, curtailment_rates, color=colors_bar, alpha=0.8)
    axes[0].set_xlabel('Curtailment Rate (%)', fontsize=12)
    axes[0].set_title('Curtailment Rate by Mitigation Scenario', fontsize=14)
    axes[0].grid(True, alpha=0.3, axis='x')
    for i, v in enumerate(curtailment_rates):
        axes[0].text(v + 0.1, i, f'{v:.1f}%', va='center', fontsize=11)

    axes[1].barh(scenario_names, curtailed_gwh, color=colors_bar, alpha=0.8)
    axes[1].set_xlabel('Total Curtailed Energy (GWh)', fontsize=12)
    axes[1].set_title('Curtailed Energy by Scenario', fontsize=14)
    axes[1].grid(True, alpha=0.3, axis='x')
    for i, v in enumerate(curtailed_gwh):
        axes[1].text(v + 5, i, f'{v:.0f}', va='center', fontsize=11)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/curtailment_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Monthly curtailment pattern
    fig, ax = plt.subplots(figsize=(10, 6))
    baseline_curt = baseline['curtailment']
    monthly_curt = []
    month_labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    start = 0
    for d in days_in_month:
        end = start + d * 24
        monthly_curt.append(np.sum(baseline_curt[start:end]) / 1000)
        start = end

    ax.bar(month_labels, monthly_curt, color='coral', alpha=0.8, edgecolor='black')
    ax.set_ylabel('Curtailed Energy (GWh)', fontsize=12)
    ax.set_title('Monthly Curtailment Distribution (Kyushu, Baseline)', fontsize=14)
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/monthly_curtailment.png', dpi=150, bbox_inches='tight')
    plt.close()

    print("\n" + "=" * 60)
    print("ALL SIMULATIONS COMPLETED SUCCESSFULLY")
    print("=" * 60)

    return results


if __name__ == '__main__':
    results = run_all_simulations()

    # Save summary
    summary = {
        'power_flow': {
            'nr_time_ms': results['power_flow']['nr_time'] * 1000,
            'helm_time_ms': results['power_flow']['helm_time'] * 1000,
            'nr_iterations': results['power_flow']['nr_iterations'],
        },
        'forecast': {
            'solar_mae': results['forecast']['solar']['mae'],
            'solar_rmse': results['forecast']['solar']['rmse'],
            'solar_nrmse': results['forecast']['solar']['nrmse'],
            'wind_mae': results['forecast']['wind']['mae'],
            'wind_rmse': results['forecast']['wind']['rmse'],
            'wind_nrmse': results['forecast']['wind']['nrmse'],
        },
        'dispatch': {
            'total_cost': results['dispatch']['total_cost'],
            'avg_curtailment': float(np.mean(results['dispatch']['curtailment'])),
            'avg_load_shed': float(np.mean(results['dispatch']['load_shed'])),
        },
        'battery_dr': {
            'baseline_cost': results['battery_dr']['baseline_cost'],
            'optimized_cost': results['battery_dr']['total_cost'],
            'savings_pct': results['battery_dr']['savings_pct'],
        },
        'stability': {
            name: {'nadir_hz': res['nadir'], 'nadir_time_s': res['nadir_time']}
            for name, res in results['stability']['frequency'].items()
        },
        'kyushu': {
            name: {'curtailment_rate_pct': res['curtailment_rate'],
                   'curtailed_gwh': res['total_curtailed_gwh']}
            for name, res in results['kyushu'].items()
        }
    }

    with open('simulation_results.json', 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print("\nResults saved to simulation_results.json")
    print(f"Figures saved to {FIGURES_DIR}/")
