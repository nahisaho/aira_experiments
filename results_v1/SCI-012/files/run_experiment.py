#!/usr/bin/env python3
"""
Main experiment script: runs all SNN framework components and generates
figures + result files.
"""

import sys, os, json, time, warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator

# Colour-blind friendly palette (Wong, 2011)
CB = ["#000000","#E69F00","#56B4E9","#009E73",
      "#F0E442","#0072B2","#D55E00","#CC79A7"]

from src.neuron_models   import benchmark_models, simulate_izhikevich, IZHIKEVICH_PRESETS
from src.plasticity       import run_plasticity_demo
from src.gpu_architecture import benchmark_scale, IzhikevichPopulation, CSRConnectivity
from src.potjans_model    import build_potjans_model, run_potjans_simulation
from src.analysis_tools   import (analyse_simulation, oscillation_power_spectrum,
                                   mutual_information_binned, transfer_entropy,
                                   population_firing_rate)
from src.working_memory   import WorkingMemoryNetwork, compare_with_experiment, WMNetworkConfig


WORKSPACE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR    = os.path.join(WORKSPACE, "figures")
RES_DIR    = os.path.join(WORKSPACE, "results")
LOG_DIR    = os.path.join(WORKSPACE, "logs")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

LOG_PATH = os.path.join(LOG_DIR, "process-log.jsonl")

def log_event(phase, event_type, skill_or_tool, files_written=None,
              handoff_in=None, handoff_out=None, status="ok"):
    entry = {
        "timestamp":     time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "phase":         phase,
        "event_type":    event_type,
        "actor":         "co-scientist",
        "skill_or_tool": skill_or_tool,
        "handoff_in":    handoff_in  or {},
        "handoff_out":   handoff_out or {},
        "files_written": files_written or [],
        "status":        status,
    }
    with open(LOG_PATH, "a") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")


def save_fig(fig, name: str) -> str:
    path = os.path.join(FIG_DIR, name)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# ============================================================
# 1. Neuron Model Comparison
# ============================================================

def run_neuron_comparison():
    print("=== [1] Neuron model comparison ===")
    log_event("1", "run_started", "neuron_models")
    T_ms = 200.0

    results = benchmark_models(T_ms=T_ms, dt_hh=0.01, dt_iz=0.1, dt_adex=0.1)

    # --- Figure 1: voltage traces ---
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=False)
    model_names = ["HH", "Izhikevich-RS", "AdEx"]
    colors = [CB[1], CB[2], CB[3]]
    y_labels = ["Membrane potential (mV)"] * 3

    for ax, mname, col in zip(axes, model_names, colors):
        res = results[mname]
        V, dt = res["V"], res["dt"]
        t = np.arange(len(V)) * dt
        ax.plot(t, V, color=col, lw=0.8)
        ax.set_ylabel("V (mV)", fontsize=9)
        ax.set_title(f"{mname}  |  FR = {res['firing_rate']:.1f} Hz  |  "
                     f"ISI CV = {res['isi_cv']:.3f}  |  "
                     f"CPU = {res['elapsed_s']*1e3:.1f} ms", fontsize=9)
        ax.set_xlim(0, T_ms)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.spines[["top","right"]].set_visible(False)

    axes[-1].set_xlabel("Time (ms)", fontsize=10)
    fig.suptitle("Neuron Model Comparison: HH, Izhikevich, AdEx", fontsize=12, fontweight="bold")
    fig.tight_layout()
    p1 = save_fig(fig, "fig1_neuron_comparison.pdf")
    fig2, _ = plt.subplots(3,1,figsize=(10,8), sharex=False)
    # Also save PNG
    fig2.clf(); plt.close(fig2)

    # Reproduce as PNG
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=False)
    for ax, mname, col in zip(axes, model_names, colors):
        res = results[mname]; V, dt = res["V"], res["dt"]
        t = np.arange(len(V)) * dt
        ax.plot(t, V, color=col, lw=0.8)
        ax.set_ylabel("V (mV)", fontsize=9)
        ax.set_title(f"{mname}  |  FR = {res['firing_rate']:.1f} Hz  |  ISI CV = {res['isi_cv']:.3f}", fontsize=9)
        ax.set_xlim(0, T_ms); ax.spines[["top","right"]].set_visible(False)
    axes[-1].set_xlabel("Time (ms)", fontsize=10)
    fig.suptitle("Neuron Model Comparison: HH, Izhikevich, AdEx", fontsize=12, fontweight="bold")
    fig.tight_layout()
    p1png = save_fig(fig, "fig1_neuron_comparison.png")

    # --- Izhikevich firing pattern zoo ---
    fig, axes = plt.subplots(2, 3, figsize=(14, 7))
    from src.neuron_models import IZHIKEVICH_PRESETS, simulate_izhikevich
    for ax, (pname, pparams) in zip(axes.flat, IZHIKEVICH_PRESETS.items()):
        N = int(500 / 0.1)
        I = np.ones(N) * 10.0
        res = simulate_izhikevich(I, dt=0.1, params=pparams)
        t = np.arange(N) * 0.1
        ax.plot(t, res["V"], color=CB[5], lw=0.7)
        ax.set_title(f"Izhikevich {pname}", fontsize=9)
        ax.set_xlabel("t (ms)", fontsize=8); ax.set_ylabel("V (mV)", fontsize=8)
        ax.spines[["top","right"]].set_visible(False)
    fig.suptitle("Izhikevich Neuron Type Zoo", fontsize=12, fontweight="bold")
    fig.tight_layout()
    p2 = save_fig(fig, "fig2_izhikevich_zoo.png")

    # Save metrics
    metrics = {k: {m: v for m, v in r.items() if m not in ("V","spikes")}
               for k, r in results.items()}
    mpath = os.path.join(RES_DIR, "neuron_model_metrics.json")
    with open(mpath, "w") as f:
        json.dump(metrics, f, indent=2, default=float)

    log_event("1", "run_completed", "neuron_models",
              files_written=[p1png, p2, mpath])
    print(f"  HH:   {results['HH']['firing_rate']:.1f} Hz, "
          f"ISI CV={results['HH']['isi_cv']:.3f}, "
          f"t={results['HH']['elapsed_s']*1e3:.0f}ms")
    print(f"  Izh:  {results['Izhikevich-RS']['firing_rate']:.1f} Hz, "
          f"ISI CV={results['Izhikevich-RS']['isi_cv']:.3f}, "
          f"t={results['Izhikevich-RS']['elapsed_s']*1e3:.0f}ms")
    print(f"  AdEx: {results['AdEx']['firing_rate']:.1f} Hz, "
          f"ISI CV={results['AdEx']['isi_cv']:.3f}, "
          f"t={results['AdEx']['elapsed_s']*1e3:.0f}ms")
    return results, metrics


# ============================================================
# 2. Synaptic Plasticity
# ============================================================

def run_plasticity():
    print("=== [2] Synaptic plasticity ===")
    log_event("2", "run_started", "plasticity")

    result = run_plasticity_demo(n_pre=200, n_post=200,
                                  T_ms=10000.0, dt=1.0, input_rate=10.0)

    wh = result["weight_history"]
    rh = result["rate_history"]
    t_vals = [r["t_ms"] for r in wh]
    w_mean = [r["mean"] for r in wh]
    r_mean = [r["mean_rate"] for r in rh]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ax1.plot(t_vals, w_mean, color=CB[1], lw=1.5, label="Mean synaptic weight")
    ax1.set_ylabel("Mean weight (a.u.)", fontsize=10)
    ax1.legend(fontsize=9); ax1.spines[["top","right"]].set_visible(False)
    ax1.set_title("STDP + Homeostatic Synaptic Scaling", fontsize=11, fontweight="bold")

    ax2.plot(t_vals, r_mean, color=CB[3], lw=1.5, label="Mean firing rate")
    ax2.axhline(5.0, ls="--", color="gray", lw=1, label="Target rate (5 Hz)")
    ax2.set_ylabel("Firing rate (Hz)", fontsize=10)
    ax2.set_xlabel("Time (ms)", fontsize=10)
    ax2.legend(fontsize=9); ax2.spines[["top","right"]].set_visible(False)

    fig.tight_layout()
    p = save_fig(fig, "fig3_plasticity.png")

    # Weight distribution before/after
    final_W = result["final_W"]
    fig, ax = plt.subplots(figsize=(7,4))
    active = final_W[final_W > 0.01].flatten()
    ax.hist(active, bins=40, color=CB[2], edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Synaptic weight", fontsize=10)
    ax.set_ylabel("Count", fontsize=10)
    ax.set_title("Final Weight Distribution (STDP + Homeostatic)", fontsize=11, fontweight="bold")
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    pw = save_fig(fig, "fig4_weight_distribution.png")

    res_path = os.path.join(RES_DIR, "plasticity_results.json")
    with open(res_path, "w") as f:
        json.dump({"weight_history": wh, "rate_history": rh,
                   "final_rate_stats": {
                       "mean": float(result["final_rates"].mean()),
                       "std":  float(result["final_rates"].std())}},
                  f, indent=2, default=float)

    log_event("2", "run_completed", "plasticity",
              files_written=[p, pw, res_path])
    print(f"  Final mean weight: {w_mean[-1]:.4f}")
    print(f"  Final mean rate:   {r_mean[-1]:.2f} Hz  (target=5 Hz)")
    return result


# ============================================================
# 3. GPU/Parallel Architecture Scaling Benchmark
# ============================================================

def run_scale_benchmark():
    print("=== [3] GPU/parallel architecture scaling ===")
    log_event("3", "run_started", "gpu_architecture")

    sizes = [1000, 5000, 10000, 50000, 100000]
    bench = benchmark_scale(sizes, T_ms=50.0, dt=0.1,
                             p_conn=0.002, backend="auto")

    Ns     = [b["N"] for b in bench]
    times  = [b["elapsed_s"] for b in bench]
    tputs  = [b["neurons_per_s"] for b in bench]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.loglog(Ns, times, "o-", color=CB[1], lw=1.5, ms=6)
    ax1.set_xlabel("Network size N", fontsize=10)
    ax1.set_ylabel("Wall time (s)", fontsize=10)
    ax1.set_title("Simulation Time vs Network Size", fontsize=11, fontweight="bold")
    ax1.grid(True, which="both", alpha=0.3)
    ax1.spines[["top","right"]].set_visible(False)

    ax2.semilogx(Ns, [t/1e6 for t in tputs], "s-", color=CB[2], lw=1.5, ms=6)
    ax2.set_xlabel("Network size N", fontsize=10)
    ax2.set_ylabel("Throughput (M neurons·steps/s)", fontsize=10)
    ax2.set_title("Computational Throughput", fontsize=11, fontweight="bold")
    ax2.grid(True, alpha=0.3)
    ax2.spines[["top","right"]].set_visible(False)

    backend = bench[0]["backend"]
    fig.suptitle(f"Parallel Architecture Benchmark  [backend: {backend}]",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    p = save_fig(fig, "fig5_scaling_benchmark.png")

    rp = os.path.join(RES_DIR, "scaling_benchmark.json")
    with open(rp, "w") as f:
        json.dump(bench, f, indent=2)
    log_event("3", "run_completed", "gpu_architecture",
              files_written=[p, rp])
    for b in bench:
        print(f"  N={b['N']:>7,}  synapses={b['n_syn']:>8,}  "
              f"t={b['elapsed_s']:.2f}s  backend={b['backend']}")
    return bench


# ============================================================
# 4. Potjans-Diesmann Cortical Microcircuit
# ============================================================

def run_potjans():
    print("=== [4] Potjans-Diesmann microcircuit ===")
    log_event("4", "run_started", "potjans_model")

    model = build_potjans_model(scale=0.05, backend="auto")
    print(f"  Total neurons: {model['meta']['total_neurons']:,}  "
          f"synapses: {model['meta']['total_synapses']:,}")

    sim = run_potjans_simulation(model, T_ms=300.0, dt=0.2)

    # Raster plot
    pops = list(sim["spike_records"].keys())
    N_total = sum(sim["N_scaled"].values())
    offsets = {}; offset = 0
    for p in pops:
        offsets[p] = offset; offset += sim["N_scaled"][p]

    fig, axes = plt.subplots(2, 1, figsize=(12, 8),
                              gridspec_kw={"height_ratios": [3, 1]})
    ax_raster = axes[0]
    colors_pop = dict(zip(pops, CB[:8]))

    for pop_name in pops:
        records = sim["spike_records"][pop_name]
        if not records:
            continue
        ts  = [r[0] for r in records]
        nids = [r[1] + offsets[pop_name] for r in records]
        ax_raster.scatter(ts, nids, s=0.3, color=colors_pop[pop_name],
                          alpha=0.6, rasterized=True, label=pop_name)

    ax_raster.set_xlim(0, sim["T_ms"])
    ax_raster.set_ylabel("Neuron index", fontsize=10)
    ax_raster.set_title("Potjans-Diesmann Cortical Microcircuit — Raster Plot",
                         fontsize=11, fontweight="bold")
    ax_raster.legend(markerscale=6, fontsize=8, ncol=4,
                     loc="upper right")
    ax_raster.spines[["top","right"]].set_visible(False)

    # Population firing rates bar chart
    ax_bar = axes[1]
    rates = [sim["rates"].get(p, 0) for p in pops]
    bars = ax_bar.bar(pops, rates, color=[colors_pop[p] for p in pops], edgecolor="white")
    ax_bar.set_ylabel("Rate (Hz)", fontsize=10)
    ax_bar.set_xlabel("Population", fontsize=10)
    ax_bar.set_title("Mean Firing Rates per Population", fontsize=10)
    ax_bar.spines[["top","right"]].set_visible(False)

    fig.tight_layout()
    p_fig = save_fig(fig, "fig6_potjans_raster.png")

    # Spectral analysis
    analysis = analyse_simulation(sim["spike_records"], sim["N_scaled"], sim["T_ms"])
    fig2, axes2 = plt.subplots(2, 4, figsize=(14, 6))
    for ax, pop_name in zip(axes2.flat, pops):
        from src.analysis_tools import spikes_to_lfp, oscillation_power_spectrum
        lfp = spikes_to_lfp(sim["spike_records"][pop_name], sim["T_ms"], dt=1.0)
        freqs, psd = oscillation_power_spectrum(lfp, fs=1000.0)
        ax.semilogy(freqs, psd + 1e-12, color=colors_pop[pop_name], lw=1.2)
        ax.set_title(pop_name, fontsize=9)
        ax.set_xlabel("Freq (Hz)", fontsize=8)
        ax.set_ylabel("PSD", fontsize=8)
        ax.set_xlim(0, 150)
        ax.spines[["top","right"]].set_visible(False)
    fig2.suptitle("LFP Power Spectra per Population", fontsize=12, fontweight="bold")
    fig2.tight_layout()
    p_spec = save_fig(fig2, "fig7_potjans_spectra.png")

    # Save results
    rp = os.path.join(RES_DIR, "potjans_results.json")
    with open(rp, "w") as f:
        ser = {k: {m: v for m, v in r.items()
                   if not isinstance(v, np.ndarray)}
               for k, r in analysis.items()}
        json.dump({"rates": sim["rates"], "meta": sim["meta"],
                   "elapsed_s": sim["elapsed_s"],
                   "analysis": ser}, f, indent=2, default=float)

    log_event("4", "run_completed", "potjans_model",
              files_written=[p_fig, p_spec, rp])
    print(f"  Simulation: {sim['elapsed_s']:.1f}s for {sim['T_ms']}ms")
    for pop_name, r in sim["rates"].items():
        print(f"    {pop_name}: {r:.2f} Hz")
    return sim, analysis


# ============================================================
# 5. Analysis Tools Demo
# ============================================================

def run_analysis_demo(sim_data, analysis):
    print("=== [5] Analysis tools ===")
    log_event("5", "run_started", "analysis_tools")

    # Mutual information between two largest populations
    pop_a, pop_b = "L23E", "L4E"
    from src.analysis_tools import population_firing_rate, mutual_information_binned, transfer_entropy

    T_ms = sim_data["T_ms"]
    N_a = sim_data["N_scaled"][pop_a]
    N_b = sim_data["N_scaled"][pop_b]
    rate_a, _ = population_firing_rate(sim_data["spike_records"][pop_a], N_a, T_ms, bin_ms=5.0)
    rate_b, _ = population_firing_rate(sim_data["spike_records"][pop_b], N_b, T_ms, bin_ms=5.0)

    mi  = mutual_information_binned(rate_a, rate_b, n_bins=8)
    te  = transfer_entropy(rate_a, rate_b, lag=1, n_bins=4)

    print(f"  MI({pop_a}→{pop_b}): {mi:.4f} bits")
    print(f"  TE({pop_a}→{pop_b}): {te:.4f} bits")

    # Figure: rate comparison + MI
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    bin_ms = 5.0
    bins_a = np.arange(len(rate_a)) * bin_ms
    axes[0,0].plot(bins_a, rate_a, color=CB[1], lw=1.0)
    axes[0,0].set_title(f"Population Rate: {pop_a}", fontsize=10)
    axes[0,0].set_ylabel("Rate (Hz)", fontsize=9)
    axes[0,0].set_xlabel("Time (ms)", fontsize=9)
    axes[0,0].spines[["top","right"]].set_visible(False)

    axes[0,1].plot(bins_a, rate_b, color=CB[2], lw=1.0)
    axes[0,1].set_title(f"Population Rate: {pop_b}", fontsize=10)
    axes[0,1].set_ylabel("Rate (Hz)", fontsize=9)
    axes[0,1].set_xlabel("Time (ms)", fontsize=9)
    axes[0,1].spines[["top","right"]].set_visible(False)

    # Scatter: rates
    axes[1,0].scatter(rate_a, rate_b, s=8, alpha=0.5, color=CB[3])
    axes[1,0].set_xlabel(f"{pop_a} rate (Hz)", fontsize=9)
    axes[1,0].set_ylabel(f"{pop_b} rate (Hz)", fontsize=9)
    axes[1,0].set_title(f"Rate Correlation  MI={mi:.3f} bits", fontsize=10)
    axes[1,0].spines[["top","right"]].set_visible(False)

    # Bar chart of band powers
    pops_list = list(analysis.keys())
    gamma_pwr = [analysis[p]["gamma_power"] if isinstance(analysis[p], dict)
                 and "gamma_power" in analysis[p] else 0 for p in pops_list
                 if p != "cross"]
    beta_pwr  = [analysis[p]["beta_power"]  if isinstance(analysis[p], dict)
                 and "beta_power"  in analysis[p] else 0 for p in pops_list
                 if p != "cross"]
    pops_filt = [p for p in pops_list if p != "cross"]
    x = np.arange(len(pops_filt))
    w = 0.35
    axes[1,1].bar(x - w/2, gamma_pwr, w, label="Gamma (30-80 Hz)", color=CB[6])
    axes[1,1].bar(x + w/2, beta_pwr,  w, label="Beta (15-30 Hz)",  color=CB[7])
    axes[1,1].set_xticks(x); axes[1,1].set_xticklabels(pops_filt, rotation=45, fontsize=8)
    axes[1,1].set_ylabel("Band power (a.u.)", fontsize=9)
    axes[1,1].set_title("Oscillatory Band Power per Population", fontsize=10)
    axes[1,1].legend(fontsize=8)
    axes[1,1].spines[["top","right"]].set_visible(False)

    fig.suptitle("Information Analysis: Firing Rates & Oscillations", fontsize=12, fontweight="bold")
    fig.tight_layout()
    p = save_fig(fig, "fig8_analysis_tools.png")

    rp = os.path.join(RES_DIR, "information_analysis.json")
    with open(rp, "w") as f:
        json.dump({"MI_L23E_L4E_bits": float(mi),
                   "TE_L23E_L4E_bits": float(te)}, f, indent=2)

    log_event("5", "run_completed", "analysis_tools",
              files_written=[p, rp])
    return {"MI": mi, "TE": te}


# ============================================================
# 6. Working Memory Task
# ============================================================

def run_working_memory():
    print("=== [6] Working memory task ===")
    log_event("6", "run_started", "working_memory")

    cfg = WMNetworkConfig(N_exc=400, N_inh=100, n_selective=60, n_assemblies=3)
    wm  = WorkingMemoryNetwork(config=cfg, backend="auto",
                                rng=np.random.default_rng(7))

    print("  Running DMS match trial (T=2000ms)...")
    t0 = time.perf_counter()
    results = wm.run(T_ms=2000.0, dt=0.5)
    elapsed = time.perf_counter() - t0
    print(f"  Elapsed: {elapsed:.1f}s")

    comparison = compare_with_experiment(results)

    # --- Raster + assembly rate figure ---
    fig, axes = plt.subplots(3, 1, figsize=(12, 10),
                              gridspec_kw={"height_ratios": [2, 1, 1]})

    ax_r = axes[0]
    T_ms = 2000.0
    # Raster (sample first 400 neurons)
    sample_n = min(400, cfg.N_exc)
    xs = [t for (t, n) in wm.exc_records if n < sample_n]
    ys = [n for (t, n) in wm.exc_records if n < sample_n]
    ax_r.scatter(xs, ys, s=0.4, color=CB[0], alpha=0.4, rasterized=True)

    # Colour assembly bands
    for a in range(cfg.n_assemblies):
        start = a * cfg.n_selective
        end   = start + cfg.n_selective
        ax_r.axhspan(start, end, alpha=0.12, color=CB[a+1])
        ax_r.text(20, (start+end)/2, f"Asm {a}", fontsize=7, va="center")

    ax_r.axvspan(0,   500,  alpha=0.08, color="blue",  label="Encode")
    ax_r.axvspan(500, 1500, alpha=0.08, color="gray",  label="Delay")
    ax_r.axvspan(1500,2000, alpha=0.08, color="orange",label="Probe")
    ax_r.set_xlim(0, T_ms); ax_r.set_ylim(0, sample_n)
    ax_r.set_ylabel("Excitatory neuron", fontsize=10)
    ax_r.set_title("Working Memory DMS Task — Raster Plot", fontsize=11, fontweight="bold")
    ax_r.legend(loc="upper right", fontsize=8, markerscale=5)
    ax_r.spines[["top","right"]].set_visible(False)

    # Assembly-0 rate time course
    ax_rate = axes[1]
    from src.analysis_tools import population_firing_rate
    a0_records = [(t,n) for (t,n) in wm.exc_records
                  if 0 <= n < cfg.n_selective]
    rate_a0, bins_a0 = population_firing_rate(a0_records, cfg.n_selective,
                                               T_ms, bin_ms=20.0)
    ax_rate.plot(bins_a0, rate_a0, color=CB[1], lw=1.5, label="Assembly 0 (cued)")
    ax_rate.axvspan(0,   500,  alpha=0.1, color="blue")
    ax_rate.axvspan(500, 1500, alpha=0.1, color="gray")
    ax_rate.axvspan(1500,2000, alpha=0.1, color="orange")
    ax_rate.axhline(5.0, ls="--", color="gray", lw=0.8, label="Baseline (~5 Hz)")
    ax_rate.set_ylabel("Rate (Hz)", fontsize=10)
    ax_rate.set_title("Assembly 0 Firing Rate (Cued stimulus)", fontsize=10)
    ax_rate.legend(fontsize=8)
    ax_rate.spines[["top","right"]].set_visible(False)

    # Summary bar: delay rates per assembly
    ax_bar = axes[2]
    asm_labels = [f"Asm {a}" for a in range(cfg.n_assemblies)]
    delay_rates = [results[f"assembly_{a}"]["delay_rate_Hz"] for a in range(cfg.n_assemblies)]
    enc_rates   = [results[f"assembly_{a}"]["encoding_rate_Hz"] for a in range(cfg.n_assemblies)]
    x  = np.arange(cfg.n_assemblies); w = 0.35
    ax_bar.bar(x - w/2, enc_rates,   w, label="Encoding rate", color=CB[4])
    ax_bar.bar(x + w/2, delay_rates, w, label="Delay rate",    color=CB[5])
    ax_bar.set_xticks(x); ax_bar.set_xticklabels(asm_labels, fontsize=9)
    ax_bar.set_ylabel("Firing rate (Hz)", fontsize=10)
    ax_bar.set_title("Encoding vs Delay Period Rates", fontsize=10)
    ax_bar.legend(fontsize=8); ax_bar.spines[["top","right"]].set_visible(False)

    fig.tight_layout()
    p = save_fig(fig, "fig9_working_memory.png")

    # Save
    rp = os.path.join(RES_DIR, "working_memory_results.json")
    with open(rp, "w") as f:
        safe = {k: (v if not isinstance(v, np.ndarray) else v.tolist())
                for k, v in results.items()}
        json.dump({"simulation": safe, "comparison": comparison,
                   "elapsed_s": elapsed}, f, indent=2, default=float)

    log_event("6", "run_completed", "working_memory",
              files_written=[p, rp])

    for a in range(cfg.n_assemblies):
        k = f"assembly_{a}"
        print(f"  Assembly {a}: enc={results[k]['encoding_rate_Hz']:.1f}Hz "
              f"delay={results[k]['delay_rate_Hz']:.1f}Hz "
              f"persistent={results[k]['persistent']}")
    return results, comparison


# ============================================================
# Main
# ============================================================

def main():
    print("\n" + "="*60)
    print("  SNN Framework — Full Experiment Run")
    print("="*60 + "\n")
    log_event("0", "run_started", "main", status="ok")

    t_total = time.perf_counter()

    neuron_results, neuron_metrics = run_neuron_comparison()
    plasticity_result               = run_plasticity()
    scale_bench                     = run_scale_benchmark()
    potjans_sim, potjans_analysis   = run_potjans()
    info_metrics                    = run_analysis_demo(potjans_sim, potjans_analysis)
    wm_results, wm_comparison       = run_working_memory()

    elapsed_total = time.perf_counter() - t_total
    print(f"\nTotal elapsed: {elapsed_total:.1f}s")

    log_event("0", "run_completed", "main",
              handoff_out={"elapsed_total_s": elapsed_total},
              status="ok")

    # Return summary for report generation
    return {
        "neuron_metrics":   neuron_metrics,
        "plasticity":       plasticity_result,
        "scale_bench":      scale_bench,
        "potjans_rates":    potjans_sim["rates"],
        "potjans_elapsed":  potjans_sim["elapsed_s"],
        "info":             info_metrics,
        "wm":               wm_results,
        "wm_comparison":    wm_comparison,
        "elapsed_total_s":  elapsed_total,
    }


if __name__ == "__main__":
    summary = main()
