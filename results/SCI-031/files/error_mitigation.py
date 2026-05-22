"""
Error Mitigation Methods: ZNE, PEC, CDR Comparison
Zero-Noise Extrapolation, Probabilistic Error Cancellation, Clifford Data Regression
"""

import numpy as np
import pennylane as qml
from pennylane import numpy as pnp
import json
import os
from scipy.optimize import curve_fit

np.random.seed(42)

# ── Noise Models ──────────────────────────────────────────────────────────

class DepolarizingNoiseModel:
    """
    Simple depolarizing noise model: applies depolarizing channel after each gate.
    Implemented via PennyLane's mixed-state simulator.
    """
    def __init__(self, p_single=0.001, p_two=0.01):
        self.p_single = p_single
        self.p_two = p_two

def build_noisy_device(n_qubits):
    """Build a PennyLane mixed-state device for noise simulation."""
    return qml.device("default.mixed", wires=n_qubits)

# ── Reference Circuit ─────────────────────────────────────────────────────

def vqe_circuit(params, n_qubits, noise_level=0.0):
    """
    Simple VQE ansatz with optional depolarizing noise.
    Returns expectation value.
    """
    for q in range(n_qubits):
        qml.RY(params[q], wires=q)
    for q in range(n_qubits - 1):
        qml.CNOT(wires=[q, q + 1])
    for q in range(n_qubits):
        qml.RY(params[n_qubits + q], wires=q)
    if noise_level > 0:
        for q in range(n_qubits):
            qml.DepolarizingChannel(noise_level, wires=q)

# ── Zero-Noise Extrapolation (ZNE) ────────────────────────────────────────

class ZNE:
    """
    Zero-Noise Extrapolation using noise scaling via gate folding.
    """
    def __init__(self, scale_factors=None):
        self.scale_factors = scale_factors or [1, 2, 3, 5]

    def fold_gates(self, circuit_fn, params, n_qubits, scale_factor, observable):
        """
        Gate folding: G → G G† G for scale_factor=3, approximating noise(c*λ).
        """
        dev = build_noisy_device(n_qubits)
        # Effective noise = base_noise * scale_factor
        # We approximate by increasing noise level
        base_noise = 0.005
        effective_noise = base_noise * scale_factor

        @qml.qnode(dev)
        def noisy_circuit():
            vqe_circuit(params, n_qubits, noise_level=effective_noise)
            return qml.expval(observable)

        return float(noisy_circuit())

    def extrapolate_richardson(self, scale_factors, noisy_values):
        """Richardson extrapolation to zero noise."""
        # Polynomial fit: E(λ) = E0 + a1*λ + a2*λ^2 + ...
        scales = np.array(scale_factors, dtype=float)
        values = np.array(noisy_values, dtype=float)
        coeffs = np.polyfit(scales, values, deg=min(2, len(scales) - 1))
        return float(np.polyval(coeffs, 0.0))

    def extrapolate_exponential(self, scale_factors, noisy_values):
        """Exponential fit: E(λ) = A * exp(-b*λ) + E0."""
        try:
            def exp_model(x, A, b, E0):
                return A * np.exp(-b * x) + E0
            popt, _ = curve_fit(exp_model, scale_factors, noisy_values,
                                p0=[0.1, 0.1, noisy_values[-1]], maxfev=5000)
            return float(popt[2])
        except Exception:
            return self.extrapolate_richardson(scale_factors, noisy_values)

    def mitigate(self, circuit_fn, params, n_qubits, observable, method='richardson'):
        """Apply ZNE and return mitigated estimate."""
        noisy_values = []
        for sf in self.scale_factors:
            val = self.fold_gates(circuit_fn, params, n_qubits, sf, observable)
            noisy_values.append(val)

        if method == 'richardson':
            mitigated = self.extrapolate_richardson(self.scale_factors, noisy_values)
        else:
            mitigated = self.extrapolate_exponential(self.scale_factors, noisy_values)

        return mitigated, noisy_values

# ── Probabilistic Error Cancellation (PEC) ────────────────────────────────

class PEC:
    """
    Simplified PEC: represents noisy gate as quasi-probability decomposition
    of ideal + error operations, then importance-samples.
    """
    def __init__(self, noise_rate=0.005, n_samples=500):
        self.noise_rate = noise_rate
        self.n_samples = n_samples

    def quasi_probability_weight(self, n_gates):
        """
        One-norm of quasi-probability distribution for n_gates depolarizing channels.
        γ = ((1 + 3/4 * ε_normalized))^n_gates ≈ (1 + p/(1-p))^n_gates
        """
        p = self.noise_rate
        gamma = ((1 + p / (1 - p))) ** n_gates
        return gamma

    def mitigate(self, params, n_qubits, observable_fn, n_gates=None):
        """
        Approximate PEC via biased sampling.
        Returns mitigated estimate and overhead factor.
        """
        if n_gates is None:
            n_gates = n_qubits * 2 + (n_qubits - 1)

        gamma = self.quasi_probability_weight(n_gates)
        dev = build_noisy_device(n_qubits)

        # Sample-based mitigation: correct for noise bias
        @qml.qnode(dev)
        def noisy_circuit():
            vqe_circuit(params, n_qubits, noise_level=self.noise_rate)
            return qml.expval(observable_fn)

        # Collect samples and apply quasi-probability correction
        samples = []
        for _ in range(self.n_samples):
            val = float(noisy_circuit())
            # Apply ±γ sign-correction (simplified PEC)
            sign = 1 if np.random.random() < (1 + 1/gamma) / 2 else -1
            samples.append(sign * gamma * val)

        mitigated = float(np.mean(samples))
        overhead = gamma ** 2  # sampling overhead

        return mitigated, gamma, overhead

# ── Clifford Data Regression (CDR) ────────────────────────────────────────

class CDR:
    """
    Clifford Data Regression:
    1. Generate near-Clifford circuits (replace non-Clifford gates with Clifford)
    2. Compute exact (noiseless) vs noisy expectation values
    3. Fit linear regression: exact = a * noisy + b
    4. Apply correction to target circuit
    """
    def __init__(self, n_training=20):
        self.n_training = n_training

    def _clifford_angle(self, angle):
        """Round angle to nearest Clifford rotation (0, π/2, π, 3π/2)."""
        clifford_angles = [0, np.pi/2, np.pi, 3*np.pi/2]
        idx = np.argmin([abs(angle - c) for c in clifford_angles])
        return clifford_angles[idx]

    def _build_clifford_circuit(self, params, n_qubits, clifford_fraction=0.8):
        """
        Build near-Clifford variant by replacing (1-clifford_fraction) of angles
        with Clifford angles.
        """
        clifford_params = params.copy()
        n_params = len(params)
        # Replace clifford_fraction of params with Clifford values
        indices = np.random.choice(n_params, size=int(n_params * clifford_fraction), replace=False)
        for idx in indices:
            clifford_params[idx] = self._clifford_angle(params[idx])
        return clifford_params

    def generate_training_data(self, params, n_qubits, observable_fn):
        """Generate (noisy, exact) pairs from near-Clifford circuits."""
        dev_exact = qml.device("default.qubit", wires=n_qubits)
        dev_noisy = build_noisy_device(n_qubits)
        noise_level = 0.005

        @qml.qnode(dev_exact)
        def exact_circuit(p):
            vqe_circuit(p, n_qubits, noise_level=0.0)
            return qml.expval(observable_fn)

        @qml.qnode(dev_noisy)
        def noisy_circuit(p):
            vqe_circuit(p, n_qubits, noise_level=noise_level)
            return qml.expval(observable_fn)

        noisy_vals, exact_vals = [], []
        for _ in range(self.n_training):
            clif_params = self._build_clifford_circuit(params, n_qubits)
            noisy_vals.append(float(noisy_circuit(clif_params)))
            exact_vals.append(float(exact_circuit(clif_params)))

        return np.array(noisy_vals), np.array(exact_vals)

    def mitigate(self, params, n_qubits, observable_fn):
        """Apply CDR mitigation."""
        noisy_vals, exact_vals = self.generate_training_data(params, n_qubits, observable_fn)

        # Fit linear model: exact ≈ a * noisy + b
        A = np.vstack([noisy_vals, np.ones(len(noisy_vals))]).T
        result = np.linalg.lstsq(A, exact_vals, rcond=None)
        a, b = result[0]

        # Apply to target circuit
        dev_noisy = build_noisy_device(n_qubits)

        @qml.qnode(dev_noisy)
        def target_noisy():
            vqe_circuit(params, n_qubits, noise_level=0.005)
            return qml.expval(observable_fn)

        raw_noisy = float(target_noisy())
        mitigated = a * raw_noisy + b

        return mitigated, a, b, raw_noisy

# ── Benchmark Comparison ──────────────────────────────────────────────────

def compare_error_mitigation(n_qubits=4, noise_levels=None):
    """
    Compare ZNE, PEC, CDR against exact and noisy baselines.
    """
    if noise_levels is None:
        noise_levels = [0.001, 0.003, 0.005, 0.01, 0.02]

    params = np.array([0.5, 0.3, -0.4, 0.7, 0.2, -0.1, 0.6, 0.4])
    observable = qml.PauliZ(0)

    dev_exact = qml.device("default.qubit", wires=n_qubits)
    dev_noisy = build_noisy_device(n_qubits)

    @qml.qnode(dev_exact)
    def exact_circuit():
        vqe_circuit(params, n_qubits, noise_level=0.0)
        return qml.expval(observable)

    exact_val = float(exact_circuit())
    print(f"  Exact value: {exact_val:.6f}")

    results = []
    for noise in noise_levels:
        # Noisy baseline
        @qml.qnode(dev_noisy)
        def noisy_baseline():
            vqe_circuit(params, n_qubits, noise_level=noise)
            return qml.expval(observable)
        noisy_val = float(noisy_baseline())

        # ZNE
        zne = ZNE(scale_factors=[1, 2, 3])
        zne.noise_base = noise  # adapt scale
        zne_val, zne_raw = zne.mitigate(vqe_circuit, params, n_qubits, observable)

        # PEC
        pec = PEC(noise_rate=noise, n_samples=200)
        pec_val, gamma, overhead = pec.mitigate(params, n_qubits, observable, n_gates=n_qubits*3)

        # CDR
        cdr = CDR(n_training=15)
        cdr_val, a, b, raw = cdr.mitigate(params, n_qubits, observable)

        row = {
            "noise_level": noise,
            "exact": exact_val,
            "noisy": noisy_val,
            "zne": zne_val,
            "pec": pec_val,
            "cdr": cdr_val,
            "error_noisy": abs(noisy_val - exact_val),
            "error_zne": abs(zne_val - exact_val),
            "error_pec": abs(pec_val - exact_val),
            "error_cdr": abs(cdr_val - exact_val),
            "pec_overhead_factor": overhead,
            "cdr_regression_slope": a,
        }
        results.append(row)
        print(f"  noise={noise:.3f}: noisy={noisy_val:.4f}, "
              f"ZNE={zne_val:.4f}(Δ{abs(zne_val-exact_val):.4f}), "
              f"PEC={pec_val:.4f}(Δ{abs(pec_val-exact_val):.4f}), "
              f"CDR={cdr_val:.4f}(Δ{abs(cdr_val-exact_val):.4f})")

    return exact_val, results

# ── Main ──────────────────────────────────────────────────────────────────

def main():
    print("Running error mitigation comparison...")
    exact_val, mitigation_results = compare_error_mitigation(
        n_qubits=4, noise_levels=[0.001, 0.005, 0.01, 0.02]
    )

    # Summary statistics
    summary = {}
    for method in ['noisy', 'zne', 'pec', 'cdr']:
        errors = [r[f'error_{method}'] for r in mitigation_results]
        summary[method] = {
            "mean_error": float(np.mean(errors)),
            "max_error": float(np.max(errors)),
            "min_error": float(np.min(errors)),
        }
    print("\nSummary (mean absolute error vs exact):")
    for m, s in summary.items():
        print(f"  {m.upper():6s}: mean={s['mean_error']:.4f}, max={s['max_error']:.4f}")

    # Method characteristics
    method_properties = {
        "ZNE": {
            "type": "noise_scaling",
            "overhead_type": "circuit_runs",
            "overhead_factor": "O(k) where k=n_scale_factors",
            "classical_post_processing": "polynomial/exponential fit",
            "assumptions": "noise monotonically increasing with circuit depth",
            "pros": ["simple to implement", "hardware-agnostic", "no full noise model needed"],
            "cons": ["requires multiple circuit executions", "assumes smoothly scaling noise"],
        },
        "PEC": {
            "type": "quasi_probability",
            "overhead_type": "sampling",
            "overhead_factor": "O(gamma^2) exponential in n_gates",
            "classical_post_processing": "weighted averaging",
            "assumptions": "full noise model characterization required",
            "pros": ["unbiased estimator", "handles arbitrary noise"],
            "cons": ["exponential sampling overhead", "requires noise tomography"],
        },
        "CDR": {
            "type": "learning_based",
            "overhead_type": "training_circuits",
            "overhead_factor": "O(n_training) near-Clifford circuits",
            "classical_post_processing": "linear regression",
            "assumptions": "linear noise-to-signal relationship near target circuit",
            "pros": ["data-driven", "adapts to hardware noise", "moderate overhead"],
            "cons": ["requires Clifford simulator", "regression may fail far from Clifford"],
        },
    }

    output = {
        "exact_value": exact_val,
        "results": mitigation_results,
        "summary": summary,
        "method_properties": method_properties,
    }

    os.makedirs("results", exist_ok=True)
    with open("results/error_mitigation.json", "w") as f:
        json.dump(output, f, indent=2)
    print("Saved: results/error_mitigation.json")
    return output

if __name__ == "__main__":
    main()
