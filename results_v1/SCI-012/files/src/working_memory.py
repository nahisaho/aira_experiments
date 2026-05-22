"""
Working memory SNN model.

Implements a decision-making / working memory task using a network of
excitatory and inhibitory populations with selective persistent activity.

Based on Wang (2002) biophysical attractor model concept, adapted for
Izhikevich neurons.

Task: delayed match-to-sample (DMS)
  - Stimulus period (0–500 ms): stimulus A activates selective group A
  - Delay period (500–1500 ms): no input → test persistent activity
  - Test period (1500–2000 ms): probe stimulus → network signals match/mismatch
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, Dict, List, Tuple
import time

from .gpu_architecture import IzhikevichPopulation, CSRConnectivity
from .analysis_tools   import population_firing_rate, spikes_to_lfp


# ---------------------------------------------------------------------------
# Network configuration
# ---------------------------------------------------------------------------

@dataclass
class WMNetworkConfig:
    N_exc:          int   = 800    # total excitatory neurons
    N_inh:          int   = 200    # total inhibitory neurons
    n_selective:    int   = 120    # neurons in each selective assembly
    n_assemblies:   int   = 4      # number of stimulus-selective assemblies

    # Within-assembly (strong) connectivity
    w_plus:         float = 1.8    # potentiated weight (attractor)
    # Cross-assembly connectivity
    w_minus:        float = 0.7    # suppressed weight
    # E→I, I→E weights
    w_ei:           float = 0.9
    w_ie:           float = -1.5
    # Background weights
    w_bg:           float = 0.5

    # Poisson background rate [Hz]
    bg_rate:        float = 2.5
    p_conn:         float = 0.3


# ---------------------------------------------------------------------------
# Working memory network
# ---------------------------------------------------------------------------

class WorkingMemoryNetwork:

    def __init__(self, config: Optional[WMNetworkConfig] = None,
                 backend: str = "auto",
                 rng: Optional[np.random.Generator] = None):
        self.cfg = config or WMNetworkConfig()
        self.rng = rng or np.random.default_rng(42)
        self.backend = backend

        cfg = self.cfg
        self.N_exc = cfg.N_exc
        self.N_inh = cfg.N_inh

        # Excitatory (RS) and inhibitory (FS) populations
        self.exc = IzhikevichPopulation(cfg.N_exc, "RS",
                                        backend=backend, rng=self.rng)
        self.inh = IzhikevichPopulation(cfg.N_inh, "FS",
                                        backend=backend, rng=self.rng)

        # Build connectivity
        self._build_connectivity()

        # Spike records
        self.exc_records: List[Tuple[float, int]] = []
        self.inh_records: List[Tuple[float, int]] = []

    def _build_connectivity(self):
        """Build structured connectivity with selective assemblies."""
        cfg = self.cfg
        rng = self.rng
        N_e, N_i = cfg.N_exc, cfg.N_inh
        ns = cfg.n_selective
        na = cfg.n_assemblies

        # Assembly membership: neurons 0..na*ns-1 belong to assemblies
        # neurons na*ns..N_e-1 are non-selective

        # E→E connectivity matrix with structured weights
        W_ee = np.zeros((N_e, N_e), dtype=np.float32)
        # Baseline connectivity
        mask_ee = rng.random((N_e, N_e)) < cfg.p_conn
        W_ee[mask_ee] = cfg.w_bg
        # Within-assembly potentiation
        for a in range(na):
            idx = slice(a * ns, (a+1) * ns)
            W_ee[idx, idx] = 0.0
            sub = rng.random((ns, ns)) < cfg.p_conn
            W_ee[idx, idx][sub] = cfg.w_plus
        # Cross-assembly suppression
        for a in range(na):
            for b in range(na):
                if a != b:
                    ia = slice(a * ns, (a+1) * ns)
                    ib = slice(b * ns, (b+1) * ns)
                    W_ee[ia, ib][W_ee[ia, ib] > 0] = cfg.w_minus
        np.fill_diagonal(W_ee, 0.0)

        # E→I, I→E
        W_ei = np.zeros((N_i, N_e), dtype=np.float32)
        m_ei = rng.random((N_i, N_e)) < cfg.p_conn
        W_ei[m_ei] = cfg.w_ei

        W_ie = np.zeros((N_e, N_i), dtype=np.float32)
        m_ie = rng.random((N_e, N_i)) < cfg.p_conn
        W_ie[m_ie] = cfg.w_ie   # inhibitory (negative)

        def _dense_to_csr(W, n_pre, n_post):
            indptr  = np.zeros(n_post + 1, dtype=np.int32)
            indices_list = []
            weights_list = []
            for j in range(n_post):
                sources = np.where(W[j] != 0)[0]
                indptr[j+1] = indptr[j] + len(sources)
                indices_list.extend(sources.tolist())
                weights_list.extend(W[j, sources].tolist())
            return CSRConnectivity(
                n_pre=n_pre, n_post=n_post,
                indptr=indptr,
                indices=np.array(indices_list, dtype=np.int32),
                weights=np.array(weights_list, dtype=np.float32),
            )

        self.conn_ee = _dense_to_csr(W_ee, N_e, N_e)
        self.conn_ei = _dense_to_csr(W_ei, N_e, N_i)  # input to inh from exc
        self.conn_ie = _dense_to_csr(W_ie, N_i, N_e)  # input to exc from inh

    def _compute_isyn_exc(self) -> np.ndarray:
        e_spk = self.exc.spikes.astype(np.float32)
        i_spk = self.inh.spikes.astype(np.float32)
        I = np.zeros(self.N_exc, dtype=np.float32)
        conn = self.conn_ee
        for j in range(self.N_exc):
            s, e = conn.indptr[j], conn.indptr[j+1]
            if e > s:
                I[j] += np.dot(conn.weights[s:e], e_spk[conn.indices[s:e]])
        conn2 = self.conn_ie
        for j in range(self.N_exc):
            s, e = conn2.indptr[j], conn2.indptr[j+1]
            if e > s:
                I[j] += np.dot(conn2.weights[s:e], i_spk[conn2.indices[s:e]])
        return I

    def _compute_isyn_inh(self) -> np.ndarray:
        e_spk = self.exc.spikes.astype(np.float32)
        I = np.zeros(self.N_inh, dtype=np.float32)
        conn = self.conn_ei
        for j in range(self.N_inh):
            s, e = conn.indptr[j], conn.indptr[j+1]
            if e > s:
                I[j] += np.dot(conn.weights[s:e], e_spk[conn.indices[s:e]])
        return I

    def run(self, T_ms: float = 2000.0, dt: float = 0.5,
            stimulus_schedule: Optional[Dict] = None) -> Dict:
        """
        Run the working memory task.

        stimulus_schedule: dict mapping (t_start_ms, t_end_ms) → assembly_index
            e.g. {(0, 500): 0, (1500, 2000): 0}  → DMS match trial
        """
        if stimulus_schedule is None:
            # Default DMS match trial
            stimulus_schedule = {
                (0.0,   500.0): 0,   # encode stimulus to assembly 0
                (1500.0, 2000.0): 0,  # probe: match
            }

        T = int(T_ms / dt)
        p_bg = self.cfg.bg_rate * dt * 1e-3

        self.exc_records.clear()
        self.inh_records.clear()

        for step in range(T):
            t_ms = step * dt

            # Synaptic currents
            I_exc = self._compute_isyn_exc()
            I_inh = self._compute_isyn_inh()

            # Background Poisson
            bg_e = self.rng.random(self.N_exc) < p_bg
            bg_i = self.rng.random(self.N_inh) < p_bg
            I_exc += bg_e.astype(np.float32) * self.cfg.w_bg * 8.0
            I_inh += bg_i.astype(np.float32) * self.cfg.w_bg * 8.0
            # Add tonic baseline drive to maintain spontaneous activity
            I_exc += 4.0
            I_inh += 3.5

            # Task stimulus input
            for (t0, t1), assembly_idx in stimulus_schedule.items():
                if t0 <= t_ms < t1:
                    start = assembly_idx * self.cfg.n_selective
                    end   = start + self.cfg.n_selective
                    I_exc[start:end] += 8.0   # strong stimulus drive

            self.exc.step(dt, I_exc)
            self.inh.step(dt, I_inh)

            # Record spikes
            for idx in np.where(self.exc.spikes)[0]:
                self.exc_records.append((t_ms, int(idx)))
            for idx in np.where(self.inh.spikes)[0]:
                self.inh_records.append((t_ms, int(idx)))

        return self._analyse(T_ms)

    def _analyse(self, T_ms: float) -> Dict:
        cfg = self.cfg
        results = {}

        # Per-assembly firing rates in delay period (500–1500 ms)
        delay_start, delay_end = 500.0, min(1500.0, T_ms)
        encode_start, encode_end = 0.0, 500.0
        test_start = 1500.0

        for a in range(cfg.n_assemblies):
            start = a * cfg.n_selective
            end   = start + cfg.n_selective
            assembly_spikes = [(t, n) for (t, n) in self.exc_records
                               if start <= n < end]

            # Delay period rate
            delay_spikes = [t for (t, n) in assembly_spikes
                            if delay_start <= t < delay_end]
            delay_dur = (delay_end - delay_start) * 1e-3
            delay_rate = len(delay_spikes) / (cfg.n_selective * delay_dur + 1e-9)

            # Encoding rate
            enc_spikes = [t for (t, n) in assembly_spikes
                          if encode_start <= t < encode_end]
            enc_dur  = (encode_end - encode_start) * 1e-3
            enc_rate = len(enc_spikes) / (cfg.n_selective * enc_dur + 1e-9)

            results[f"assembly_{a}"] = {
                "encoding_rate_Hz": enc_rate,
                "delay_rate_Hz":    delay_rate,
                "persistent":       delay_rate > 5.0,   # threshold for WM
            }

        # LFP
        lfp_exc = spikes_to_lfp(self.exc_records, T_ms, dt=1.0)
        results["lfp"] = lfp_exc

        return results


# ---------------------------------------------------------------------------
# Experimental comparison (synthetic experimental data)
# ---------------------------------------------------------------------------

def compare_with_experiment(sim_results: Dict,
                             rng: Optional[np.random.Generator] = None) -> Dict:
    """
    Compare simulation results with synthetic experimental benchmarks.

    In a real study these would be loaded from electrophysiology recordings.
    Here we use published target values from Miller et al. (2018) and
    Funahashi et al. (1989) primate PFC working memory experiments.
    """
    if rng is None:
        rng = np.random.default_rng(99)

    # Published target values (mean ± SD from literature)
    exp_benchmarks = {
        "delay_rate_target_Hz":    8.0,    # PFC delay activity [Funahashi 1989]
        "delay_rate_sd":           3.0,
        "encoding_rate_target_Hz": 25.0,   # stimulus-evoked rate
        "encoding_rate_sd":        8.0,
    }

    comparison = {"benchmarks": exp_benchmarks, "assemblies": {}}

    for key, val in sim_results.items():
        if not key.startswith("assembly_"):
            continue
        delay  = val["delay_rate_Hz"]
        enc    = val["encoding_rate_Hz"]

        # z-scores vs experimental distributions
        z_delay = (delay - exp_benchmarks["delay_rate_target_Hz"]) / exp_benchmarks["delay_rate_sd"]
        z_enc   = (enc   - exp_benchmarks["encoding_rate_target_Hz"]) / exp_benchmarks["encoding_rate_sd"]

        comparison["assemblies"][key] = {
            "sim_delay_Hz":    delay,
            "sim_encoding_Hz": enc,
            "z_delay":         z_delay,
            "z_encoding":      z_enc,
            "match":           abs(z_delay) < 2.0 and abs(z_enc) < 2.0,
        }

    return comparison
