"""
Online Analysis & Feedback Control System
===========================================
HPLC/IR integration with PID feedback control for continuous flow reactors.
Includes PAT (Process Analytical Technology) data processing.
"""

import numpy as np
import json, os

np.random.seed(42)

class PIDController:
    """Discrete PID controller with anti-windup."""

    def __init__(self, Kp, Ki, Kd, setpoint, output_limits=(0, 100), dt=1.0):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.setpoint = setpoint
        self.output_limits = output_limits
        self.dt = dt
        self.integral = 0.0
        self.prev_error = 0.0

    def update(self, measurement):
        error = self.setpoint - measurement
        self.integral += error * self.dt
        # Anti-windup clamp
        integral_term = self.Ki * self.integral
        integral_term = np.clip(integral_term, *self.output_limits)
        self.integral = integral_term / self.Ki if self.Ki != 0 else 0

        derivative = (error - self.prev_error) / self.dt
        self.prev_error = error

        output = self.Kp * error + integral_term + self.Kd * derivative
        output = np.clip(output, *self.output_limits)
        return output, error

class OnlineHPLC:
    """Simulated online HPLC analyzer."""

    def __init__(self, analysis_time=60, delay=30):
        self.analysis_time = analysis_time  # seconds per analysis
        self.delay = delay                   # transport delay
        self.measurement_history = []

    def analyze(self, true_composition, time_stamp):
        noise = np.random.normal(0, 0.5, len(true_composition))
        measured = {k: max(0, v + n) for (k, v), n in
                    zip(true_composition.items(), noise)}
        self.measurement_history.append({
            "time_s": time_stamp,
            "composition": {k: round(v, 2) for k, v in measured.items()},
        })
        return measured

class OnlineIR:
    """Simulated inline FTIR analyzer with chemometric model."""

    def __init__(self, calibration_r2=0.98, sampling_interval=5):
        self.calibration_r2 = calibration_r2
        self.sampling_interval = sampling_interval
        self.measurement_history = []

    def predict_conversion(self, true_conversion, time_stamp):
        noise = np.random.normal(0, 1.0)
        measured = max(0, min(100, true_conversion + noise))
        self.measurement_history.append({
            "time_s": time_stamp,
            "conversion_pct": round(measured, 2),
        })
        return measured

def simulate_reactor_with_control(duration=3600, dt=5, target_yield=85.0):
    """Simulate continuous flow reactor with PID feedback control."""

    # PID for temperature control
    pid_temp = PIDController(Kp=2.0, Ki=0.1, Kd=0.5, setpoint=target_yield,
                             output_limits=(50, 150), dt=dt)
    # PID for flow rate
    pid_flow = PIDController(Kp=0.05, Ki=0.005, Kd=0.01, setpoint=target_yield,
                             output_limits=(0.1, 5.0), dt=dt)

    hplc = OnlineHPLC(analysis_time=60, delay=30)
    ir = OnlineIR(sampling_interval=5)

    temperature = 80.0
    flow_rate = 1.0

    time_series = []
    n_steps = int(duration / dt)
    hplc_interval = 60  # HPLC every 60s

    for step in range(n_steps):
        t = step * dt

        # Disturbances
        temp_disturbance = 2 * np.sin(2 * np.pi * t / 600) + np.random.normal(0, 0.5)
        flow_disturbance = 0.05 * np.sin(2 * np.pi * t / 300) + np.random.normal(0, 0.02)

        actual_temp = temperature + temp_disturbance
        actual_flow = max(0.05, flow_rate + flow_disturbance)

        # Reaction model
        Ea = 60000
        R = 8.314
        k = 1e8 * np.exp(-Ea / (R * (actual_temp + 273.15)))
        tau = 30.0 / actual_flow
        Da = k * 1.0 * tau
        conversion = Da / (1 + Da) * 100

        # Step change at t=1800s (simulating process upset)
        if 1800 <= t < 2100:
            conversion *= 0.85

        # IR measurement (every dt)
        ir_reading = ir.predict_conversion(conversion, t)

        # HPLC measurement (every 60s)
        hplc_data = None
        if step % (hplc_interval // dt) == 0 and step > 0:
            composition = {
                "product": conversion * 0.95,
                "starting_material": (100 - conversion) * 0.9,
                "impurity_A": conversion * 0.03,
                "impurity_B": conversion * 0.02,
            }
            hplc_data = hplc.analyze(composition, t)

        # PID control updates (using IR feedback)
        temp_output, temp_error = pid_temp.update(ir_reading)
        flow_output, flow_error = pid_flow.update(ir_reading)

        # Apply control actions (with rate limiting)
        temperature += np.clip(temp_output - temperature, -2, 2) * 0.1
        flow_rate += np.clip(flow_output - flow_rate, -0.1, 0.1) * 0.1

        time_series.append({
            "time_s": t,
            "temperature_C": round(actual_temp, 2),
            "flow_rate_mL_min": round(actual_flow, 3),
            "conversion_pct": round(conversion, 2),
            "ir_reading_pct": round(ir_reading, 2),
            "temp_setpoint_C": round(temperature, 2),
            "flow_setpoint_mL_min": round(flow_rate, 3),
        })

    # Calculate control performance metrics
    conversions = [ts["ir_reading_pct"] for ts in time_series]
    steady_state = conversions[len(conversions)//3:]  # last 2/3

    results = {
        "control_config": {
            "target_yield_pct": target_yield,
            "pid_temperature": {"Kp": 2.0, "Ki": 0.1, "Kd": 0.5},
            "pid_flow_rate": {"Kp": 0.05, "Ki": 0.005, "Kd": 0.01},
            "hplc_interval_s": hplc_interval,
            "ir_interval_s": dt,
            "simulation_duration_s": duration,
        },
        "performance_metrics": {
            "mean_conversion_pct": round(np.mean(steady_state), 2),
            "std_conversion_pct": round(np.std(steady_state), 2),
            "max_deviation_pct": round(max(abs(c - target_yield) for c in steady_state), 2),
            "settling_time_s": _estimate_settling_time(conversions, target_yield, dt),
            "overshoot_pct": round(max(conversions[:200]) - target_yield, 2),
            "ise": round(sum((c - target_yield)**2 for c in steady_state) * dt, 1),
        },
        "pat_summary": {
            "n_ir_measurements": len(ir.measurement_history),
            "n_hplc_measurements": len(hplc.measurement_history),
            "ir_calibration_r2": ir.calibration_r2,
        },
    }

    return results, time_series

def _estimate_settling_time(conversions, target, dt, band=2.0):
    """Estimate 2% settling time."""
    for i in range(len(conversions) - 1, -1, -1):
        if abs(conversions[i] - target) > band:
            return (i + 1) * dt
    return 0

if __name__ == "__main__":
    results, time_series = simulate_reactor_with_control()
    os.makedirs("results", exist_ok=True)
    os.makedirs("data", exist_ok=True)

    with open("results/feedback_control_results.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    # Save time series (downsampled)
    downsampled = time_series[::12]  # every 60s
    with open("data/control_timeseries.json", "w") as f:
        json.dump(downsampled, f, indent=2)

    print("=== Feedback Control Results ===")
    print(f"Target Yield: {results['control_config']['target_yield_pct']}%")
    print(f"Mean Conversion: {results['performance_metrics']['mean_conversion_pct']}%")
    print(f"Std Deviation: {results['performance_metrics']['std_conversion_pct']}%")
    print(f"Settling Time: {results['performance_metrics']['settling_time_s']} s")
    print(f"ISE: {results['performance_metrics']['ise']}")
    print(f"\nPAT Measurements: {results['pat_summary']['n_ir_measurements']} IR, "
          f"{results['pat_summary']['n_hplc_measurements']} HPLC")
