"""
Biologically plausible neuron models: Hodgkin-Huxley, Izhikevich, AdEx.
All models implement the same interface for benchmarking.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Hodgkin-Huxley (HH) model
# ---------------------------------------------------------------------------

@dataclass
class HHParams:
    C_m: float = 1.0       # membrane capacitance [µF/cm²]
    g_Na: float = 120.0    # max Na conductance [mS/cm²]
    g_K: float = 36.0      # max K conductance [mS/cm²]
    g_L: float = 0.3       # leak conductance [mS/cm²]
    E_Na: float = 50.0     # Na reversal [mV]
    E_K: float = -77.0     # K reversal [mV]
    E_L: float = -54.387   # leak reversal [mV]
    V_thresh: float = 0.0  # spike detection threshold [mV]


def _hh_alpha_m(V): return 0.1*(V+40)/(1-np.exp(-(V+40)/10)+1e-7)
def _hh_beta_m(V):  return 4*np.exp(-(V+65)/18)
def _hh_alpha_h(V): return 0.07*np.exp(-(V+65)/20)
def _hh_beta_h(V):  return 1/(np.exp(-(V+35)/10)+1)
def _hh_alpha_n(V): return 0.01*(V+55)/(1-np.exp(-(V+55)/10)+1e-7)
def _hh_beta_n(V):  return 0.125*np.exp(-(V+65)/80)


def simulate_hh(I_ext: np.ndarray, dt: float = 0.01,
                params: Optional[HHParams] = None) -> dict:
    """Simulate single HH neuron.

    Args:
        I_ext: external current array [µA/cm²], shape (T,)
        dt: time step [ms]
        params: HHParams instance
    Returns:
        dict with V, m, h, n, spikes arrays
    """
    if params is None:
        params = HHParams()
    p = params
    T = len(I_ext)

    V = np.zeros(T); m = np.zeros(T); h = np.zeros(T); n = np.zeros(T)
    # initial conditions near resting state
    V[0] = -65.0
    m[0] = _hh_alpha_m(V[0]) / (_hh_alpha_m(V[0]) + _hh_beta_m(V[0]))
    h[0] = _hh_alpha_h(V[0]) / (_hh_alpha_h(V[0]) + _hh_beta_h(V[0]))
    n[0] = _hh_alpha_n(V[0]) / (_hh_alpha_n(V[0]) + _hh_beta_n(V[0]))

    spikes = []
    for t in range(T - 1):
        v, mv, hv, nv = V[t], m[t], h[t], n[t]
        I_Na = p.g_Na * mv**3 * hv * (v - p.E_Na)
        I_K  = p.g_K  * nv**4      * (v - p.E_K)
        I_L  = p.g_L               * (v - p.E_L)

        dV = (I_ext[t] - I_Na - I_K - I_L) / p.C_m
        dm = _hh_alpha_m(v)*(1-mv) - _hh_beta_m(v)*mv
        dh = _hh_alpha_h(v)*(1-hv) - _hh_beta_h(v)*hv
        dn = _hh_alpha_n(v)*(1-nv) - _hh_beta_n(v)*nv

        V[t+1] = v + dt*dV
        m[t+1] = mv + dt*dm
        h[t+1] = hv + dt*dh
        n[t+1] = nv + dt*dn

        if V[t] < p.V_thresh <= V[t+1]:
            spikes.append(t * dt)

    return {"V": V, "m": m, "h": h, "n": n,
            "spikes": np.array(spikes), "dt": dt, "model": "HH"}


# ---------------------------------------------------------------------------
# Izhikevich model
# ---------------------------------------------------------------------------

@dataclass
class IzhikevichParams:
    """Izhikevich neuron parameters. Default: regular spiking (RS)."""
    a: float = 0.02
    b: float = 0.2
    c: float = -65.0   # reset voltage [mV]
    d: float = 8.0     # reset recovery variable
    V_peak: float = 30.0
    V_thresh: float = 30.0


IZHIKEVICH_PRESETS = {
    "RS":   IzhikevichParams(0.02, 0.2,  -65.0, 8.0),   # Regular Spiking
    "IB":   IzhikevichParams(0.02, 0.2,  -55.0, 4.0),   # Intrinsically Bursting
    "CH":   IzhikevichParams(0.02, 0.2,  -50.0, 2.0),   # Chattering
    "FS":   IzhikevichParams(0.1,  0.2,  -65.0, 2.0),   # Fast Spiking
    "LTS":  IzhikevichParams(0.02, 0.25, -65.0, 2.0),   # Low-Threshold Spiking
    "TC":   IzhikevichParams(0.02, 0.25, -65.0, 0.05),  # Thalamocortical
}


def simulate_izhikevich(I_ext: np.ndarray, dt: float = 0.1,
                        params: Optional[IzhikevichParams] = None) -> dict:
    """Simulate single Izhikevich neuron."""
    if params is None:
        params = IzhikevichParams()
    p = params
    T = len(I_ext)

    V = np.zeros(T); u = np.zeros(T)
    V[0] = -65.0
    u[0] = p.b * V[0]
    spikes = []

    for t in range(T - 1):
        v, uv = V[t], u[t]
        if v >= p.V_peak:
            V[t] = p.V_peak   # mark spike peak
            v = p.c
            uv = uv + p.d
            spikes.append(t * dt)

        dV = 0.04*v**2 + 5*v + 140 - uv + I_ext[t]
        du = p.a * (p.b*v - uv)
        V[t+1] = v + dt*dV
        u[t+1] = uv + dt*du

    return {"V": V, "u": u, "spikes": np.array(spikes),
            "dt": dt, "model": "Izhikevich"}


# ---------------------------------------------------------------------------
# Adaptive Exponential Integrate-and-Fire (AdEx) model
# ---------------------------------------------------------------------------

@dataclass
class AdExParams:
    C:       float = 281.0    # membrane capacitance [pF]
    g_L:     float = 30.0     # leak conductance [nS]
    E_L:     float = -70.6    # leak reversal [mV]
    V_T:     float = -50.4    # threshold slope factor [mV]
    Delta_T: float = 2.0      # slope factor [mV]
    tau_w:   float = 144.0    # adaptation time constant [ms]
    a:       float = 4.0      # subthreshold coupling [nS]
    b:       float = 80.5     # spike-triggered adaptation increment [pA]
    V_peak:  float = 20.0     # spike peak [mV]
    V_reset: float = -70.6    # reset voltage [mV]
    # Brette & Gerstner (2005) params for regular spiking


def simulate_adex(I_ext: np.ndarray, dt: float = 0.1,
                  params: Optional[AdExParams] = None) -> dict:
    """Simulate single AdEx neuron."""
    if params is None:
        params = AdExParams()
    p = params
    T = len(I_ext)

    V = np.zeros(T); w = np.zeros(T)
    V[0] = p.E_L
    spikes = []

    for t in range(T - 1):
        v, wv = V[t], w[t]
        if v >= p.V_peak:
            V[t] = p.V_peak
            v = p.V_reset
            wv = wv + p.b
            spikes.append(t * dt)

        exp_term = p.Delta_T * np.exp((v - p.V_T) / p.Delta_T)
        dV = (-p.g_L*(v - p.E_L) + p.g_L*exp_term - wv + I_ext[t]) / p.C
        dw = (p.a*(v - p.E_L) - wv) / p.tau_w
        V[t+1] = v + dt*dV
        w[t+1] = wv + dt*dw

    return {"V": V, "w": w, "spikes": np.array(spikes),
            "dt": dt, "model": "AdEx"}


# ---------------------------------------------------------------------------
# Benchmark utility
# ---------------------------------------------------------------------------

def compute_firing_rate(spikes: np.ndarray, T_total_ms: float) -> float:
    """Mean firing rate in Hz."""
    return len(spikes) / (T_total_ms * 1e-3) if T_total_ms > 0 else 0.0


def compute_isi_cv(spikes: np.ndarray) -> float:
    """Coefficient of variation of inter-spike intervals."""
    if len(spikes) < 2:
        return float('nan')
    isi = np.diff(spikes)
    return float(np.std(isi) / (np.mean(isi) + 1e-12))


def benchmark_models(T_ms: float = 1000.0, dt_hh: float = 0.01,
                     dt_iz: float = 0.1, dt_adex: float = 0.1) -> dict:
    """Run all three models with a step current and return metrics."""
    import time

    results = {}

    # Hodgkin-Huxley
    N_hh  = int(T_ms / dt_hh)
    I_hh  = np.ones(N_hh) * 10.0   # [µA/cm²]
    t0 = time.perf_counter()
    res_hh = simulate_hh(I_hh, dt=dt_hh)
    elapsed_hh = time.perf_counter() - t0
    results["HH"] = {
        "V": res_hh["V"], "spikes": res_hh["spikes"],
        "firing_rate": compute_firing_rate(res_hh["spikes"], T_ms),
        "isi_cv": compute_isi_cv(res_hh["spikes"]),
        "elapsed_s": elapsed_hh, "dt": dt_hh
    }

    # Izhikevich (RS)
    N_iz = int(T_ms / dt_iz)
    I_iz = np.ones(N_iz) * 10.0
    t0 = time.perf_counter()
    res_iz = simulate_izhikevich(I_iz, dt=dt_iz)
    elapsed_iz = time.perf_counter() - t0
    results["Izhikevich-RS"] = {
        "V": res_iz["V"], "spikes": res_iz["spikes"],
        "firing_rate": compute_firing_rate(res_iz["spikes"], T_ms),
        "isi_cv": compute_isi_cv(res_iz["spikes"]),
        "elapsed_s": elapsed_iz, "dt": dt_iz
    }

    # AdEx
    N_adex = int(T_ms / dt_adex)
    I_adex = np.ones(N_adex) * 700.0   # [pA]
    t0 = time.perf_counter()
    res_adex = simulate_adex(I_adex, dt=dt_adex)
    elapsed_adex = time.perf_counter() - t0
    results["AdEx"] = {
        "V": res_adex["V"], "spikes": res_adex["spikes"],
        "firing_rate": compute_firing_rate(res_adex["spikes"], T_ms),
        "isi_cv": compute_isi_cv(res_adex["spikes"]),
        "elapsed_s": elapsed_adex, "dt": dt_adex
    }

    return results
