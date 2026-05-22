"""
Semiconductor Process Integration: 7nm and Below Patterning
BCP Self-Assembly for Advanced Lithography (EUV + DSA Hybrid)

Topics:
  1. Process flow for BCP DSA (chemoepitaxy + EUV guide)
  2. Critical dimension (CD) analysis and process window
  3. Line-width roughness (LWR) budget
  4. Defect density targets (ITRS/IRDS roadmap)
  5. Thermal budget optimization
  6. EUV + DSA hybrid patterning scheme

References:
  - IRDS 2023 Lithography Roadmap
  - Patel et al., Proc. SPIE 2022 (EUV+DSA)
  - Ji et al., ACS Nano 2019 (sub-7nm BCP)
  - Wan et al., Nature 2021 (high-chi BCP)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import json, os
from dataclasses import dataclass
from typing import List, Dict

# ─────────────────────────────────────────────────────────────────────────────
# Process Roadmap
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ProcessNode:
    name: str
    half_pitch_nm: float      # target half-pitch
    CD_nm: float              # critical dimension
    LWR_3sigma_nm: float      # line-width roughness spec
    defect_density_cm2: float # defects/cm^2 target
    BCP_N: int                # required chain length
    BCP_chi: float            # required χ parameter
    L0_nm: float              # BCP natural period
    DSA_n: int                # multiplication factor
    guide_pitch_nm: float     # EUV guide pitch
    note: str

PROCESS_NODES = [
    ProcessNode("28nm node",  28.0, 28.0, 6.0, 0.001,  80,  0.030, 28.0,  1, 28.0,  "PS-b-PMMA conventional"),
    ProcessNode("14nm node",  14.0, 14.0, 3.5, 0.0001, 200, 0.045, 14.0,  2, 28.0,  "PS-b-PMMA 2× DSA"),
    ProcessNode("10nm node",  10.0, 10.0, 2.5, 1e-5,   350, 0.060, 10.0,  3, 30.0,  "High-χ BCP 3× DSA"),
    ProcessNode("7nm node",    7.0,  7.0, 2.0, 1e-6,   500, 0.090,  7.0,  4, 28.0,  "High-χ BCP 4× DSA"),
    ProcessNode("5nm node",    5.0,  5.0, 1.5, 1e-7,   800, 0.120,  5.0,  5, 25.0,  "High-χ BCP 5× DSA"),
    ProcessNode("3nm node",    3.0,  3.0, 1.2, 1e-8,  1200, 0.180,  3.0,  8, 24.0,  "Ultra-high-χ BCP 8× DSA"),
]

# ─────────────────────────────────────────────────────────────────────────────
# High-Chi BCP Materials
# ─────────────────────────────────────────────────────────────────────────────

HIGH_CHI_BCPS = {
    "PS-b-PMMA":         {"chi_500K":0.040+4.9/500, "L0_min_nm":12, "Tg_K":373, "maturity":"production"},
    "PS-b-P2VP":         {"chi_500K":0.12,           "L0_min_nm":8,  "Tg_K":368, "maturity":"pilot"},
    "PLA-b-PS":          {"chi_500K":0.18,           "L0_min_nm":6,  "Tg_K":330, "maturity":"R&D"},
    "PDMS-b-PS":         {"chi_500K":0.15,           "L0_min_nm":7,  "Tg_K":423, "maturity":"pilot"},
    "PTMSS-b-PMOST":     {"chi_500K":0.25,           "L0_min_nm":5,  "Tg_K":420, "maturity":"R&D"},
    "Si-containing-BCP": {"chi_500K":0.30,           "L0_min_nm":4,  "Tg_K":450, "maturity":"research"},
    "PEO-b-PS":          {"chi_500K":0.10,           "L0_min_nm":9,  "Tg_K":340, "maturity":"R&D"},
}

def estimate_chiN_for_ordering(chi: float, N_target: float = 12.0) -> int:
    """Minimum N for χN > N_target (strong segregation)."""
    return int(np.ceil(N_target / chi))

# ─────────────────────────────────────────────────────────────────────────────
# Process Window Analysis
# ─────────────────────────────────────────────────────────────────────────────

class ProcessWindowAnalyzer:
    """
    Compute process window for BCP DSA lithography.
    Metrics: CD uniformity, LWR, defect density vs process conditions.
    """
    def __init__(self, node: ProcessNode):
        self.node = node

    def cd_bias(self, T_K: float, t_anneal_s: float) -> float:
        """
        CD bias from nominal: model as function of anneal conditions.
        CD_bias = CD_target * f(T, t)
        """
        T0  = 500.0  # reference temperature (K)
        t0  = 1e-3   # reference time (s)
        dt  = (T_K - T0) / T0
        tau = np.log(t_anneal_s / t0)
        return self.node.CD_nm * (0.02 * dt - 0.01 * tau + 0.005 * dt * tau)

    def lwr_model(self, T_K: float, chi: float) -> float:
        """
        LWR scaling: LWR ~ L0 / (chi*N)^0.5 * correction(T)
        """
        L0   = self.node.L0_nm
        chiN = chi * self.node.BCP_N
        base_lwr = L0 / np.sqrt(chiN) * 2.5  # nm (empirical coefficient)
        T_correction = 1.0 + 0.5 * np.exp(-(T_K - 400) / 100)
        return base_lwr * T_correction

    def defect_density_model(self, T_K: float, t_anneal_s: float,
                              rho0: float = 1e9) -> float:
        """
        Defect density after anneal: rho(T,t) using Arrhenius model.
        """
        Ea   = 50.0   # kJ/mol activation energy
        kB   = 8.314e-3
        k_ann = 1e12 * np.exp(-Ea / (kB * T_K))
        rho_eq = 1e6 * np.exp(-2*Ea / (kB * T_K))
        tau  = 1.0 / (k_ann * rho0)
        rho  = rho_eq + (rho0 - rho_eq) * np.exp(-t_anneal_s / tau)
        return max(rho, rho_eq)

    def compute_process_window(self, T_range=(420, 540, 13),
                               t_range=(1e-4, 0.1, 10)):
        """
        2D process window: (T, t_anneal) space satisfying all specs.
        """
        T_arr = np.linspace(*T_range)
        t_arr = np.logspace(np.log10(t_range[0]), np.log10(t_range[1]), t_range[2])
        TT, tt = np.meshgrid(T_arr, t_arr)

        CD_b   = np.vectorize(self.cd_bias)(TT, tt)
        LWR    = np.vectorize(lambda T,t: self.lwr_model(T, self.node.BCP_chi))(TT, tt)
        def_d  = np.vectorize(self.defect_density_model)(TT, tt)

        CD_pass  = np.abs(CD_b) < self.node.CD_nm * 0.05  # ±5% CD
        LWR_pass = LWR < self.node.LWR_3sigma_nm
        def_pass = def_d < self.node.defect_density_cm2 * 1e8  # cm^-2 -> m^-2

        process_window = CD_pass & LWR_pass & def_pass
        window_fraction = process_window.mean()

        return {
            "T_arr": T_arr.tolist(),
            "t_arr": t_arr.tolist(),
            "CD_bias": CD_b.tolist(),
            "LWR": LWR.tolist(),
            "defect_density": def_d.tolist(),
            "process_window": process_window.tolist(),
            "window_fraction": float(window_fraction),
        }


# ─────────────────────────────────────────────────────────────────────────────
# EUV + DSA Hybrid Patterning
# ─────────────────────────────────────────────────────────────────────────────

def euv_dsa_hybrid_flow() -> dict:
    """
    Process flow for EUV + DSA hybrid patterning (7nm node, lamellae).
    """
    steps = [
        {"step": 1, "process": "Substrate preparation",
         "description": "Si wafer + spin-on neutral brush (PS-r-PMMA, f=0.5), 250°C bake",
         "critical_params": {"brush_Mn": 5000, "bake_T_C": 250, "bake_t_min": 5}},
        {"step": 2, "process": "EUV guide patterning",
         "description": "EUV lithography: expose chemical stripe pattern at L_guide = 28nm",
         "critical_params": {"exposure_mJ_cm2": 35, "CD_guide_nm": 28, "LWR_euv_nm": 3.0}},
        {"step": 3, "process": "Guide pattern development",
         "description": "Develop EUV resist, rinse, descum O2 plasma",
         "critical_params": {"develop_time_s": 30, "plasma_W": 50, "plasma_t_s": 5}},
        {"step": 4, "process": "BCP coating",
         "description": "Spin-coat PS-b-PMMA (Mn=45k-45k, PDI<1.05), 1500 rpm, 30s",
         "critical_params": {"Mn_PS": 45000, "Mn_PMMA": 45000, "PDI": 1.04,
                              "thickness_nm": 28, "spin_rpm": 1500}},
        {"step": 5, "process": "Thermal annealing",
         "description": "RTP: 500°C/500ms or hotplate: 250°C/5min under N2",
         "critical_params": {"anneal_T_C": 250, "anneal_t_min": 5,
                              "atmosphere": "N2", "ramp_rate_C_s": 10}},
        {"step": 6, "process": "UV + acid treatment (PMMA removal)",
         "description": "254nm UV 10 min, acetic acid soak 1 min, N2 rinse",
         "critical_params": {"UV_nm": 254, "UV_min": 10, "acid_s": 60}},
        {"step": 7, "process": "Pattern transfer",
         "description": "Directional RIE: CF4/O2 to etch PMMA, then Si etch",
         "critical_params": {"etch_gas": "CF4:O2 = 9:1", "etch_power_W": 100,
                              "etch_rate_nm_min": 30, "selectivity": 8}},
        {"step": 8, "process": "Metrology",
         "description": "CD-SEM: CD uniformity, LWR, defect inspection (e-beam)",
         "critical_params": {"CD_target_nm": 7, "LWR_3sigma_limit_nm": 2.0,
                              "defect_limit_cm2": 1e-6}},
    ]
    return {"process_flow": steps, "target_node": "7nm", "L0_nm": 14.0,
            "multiplication": "2× (EUV 28nm guide → 14nm BCP period → 7nm CD)"}


# ─────────────────────────────────────────────────────────────────────────────
# Visualization
# ─────────────────────────────────────────────────────────────────────────────

def plot_semiconductor_analysis(output_path="figures/semiconductor_process.png"):
    fig, axes = plt.subplots(2, 3, figsize=(17, 11))

    # 1. Chi vs material / Temperature
    ax = axes[0][0]
    T_arr = np.linspace(300, 600, 100)
    for name, data in list(HIGH_CHI_BCPS.items())[:5]:
        chi_T = [data["chi_500K"] * 500/T for T in T_arr]
        ax.plot(T_arr, chi_T, lw=2, label=name)
    ax.axhline(0.04+4.9/500, color='k', ls=':', lw=1, alpha=0.5)
    ax.set_xlabel("Temperature (K)", fontsize=11)
    ax.set_ylabel("Flory-Huggins χ", fontsize=11)
    ax.set_title("High-χ BCP Materials Comparison", fontsize=12)
    ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=0.3)

    # 2. L0_min vs chi
    ax = axes[0][1]
    chi_arr = np.linspace(0.02, 0.35, 100)
    b = 0.68  # nm
    # L0_min from SSL: L0 = 1.03*b*N^(2/3)*chi^(1/6) at N=N_min(chi)
    # N_min: chiN=10.5 -> N_min=10.5/chi
    N_min = 10.5 / chi_arr
    L0_min = 1.03*b*N_min**(2/3)*chi_arr**(1/6)*(0.25)**(2/3)
    ax.plot(chi_arr, L0_min, 'b-', lw=2.5, label="L₀_min (SSL)")
    ax.axhline(7, color='red',   ls='--', lw=1.5, label="7nm half-pitch")
    ax.axhline(5, color='orange',ls='--', lw=1.5, label="5nm half-pitch")
    ax.axhline(3, color='purple',ls='--', lw=1.5, label="3nm half-pitch")
    for name, data in HIGH_CHI_BCPS.items():
        ax.scatter([data["chi_500K"]], [data["L0_min_nm"]],
                   s=80, zorder=5, label=f"{name.split('-b-')[0]}-b-*")
    ax.set_xlabel("Flory-Huggins χ (T=500K)", fontsize=11)
    ax.set_ylabel("Minimum L₀ (nm)", fontsize=11)
    ax.set_title("Accessible Period vs χ Parameter", fontsize=12)
    ax.legend(fontsize=7, loc="upper right"); ax.grid(alpha=0.3)
    ax.set_ylim(0, 30)

    # 3. IRDS Roadmap
    ax = axes[0][2]
    nodes  = [n.half_pitch_nm for n in PROCESS_NODES]
    dsaN   = [n.DSA_n         for n in PROCESS_NODES]
    bcpN   = [n.BCP_N         for n in PROCESS_NODES]
    ax2 = ax.twinx()
    ax.bar(range(len(nodes)), nodes, color='#3498DB', alpha=0.7, label="Half-pitch (nm)")
    ax2.plot(range(len(nodes)), dsaN, 'r-o', lw=2, ms=8, label="DSA multiplication n")
    ax2.plot(range(len(nodes)), [N/50 for N in bcpN], 'g-s', lw=2, ms=8,
             label="N/50 (chain length scale)")
    ax.set_xticks(range(len(nodes)))
    ax.set_xticklabels([n.name.replace(" node","") for n in PROCESS_NODES],
                       rotation=25, fontsize=9)
    ax.set_ylabel("Half-Pitch (nm)", fontsize=11)
    ax2.set_ylabel("DSA n / N/50", fontsize=11)
    ax.set_title("IRDS Process Roadmap", fontsize=12)
    lines1,lab1 = ax.get_legend_handles_labels()
    lines2,lab2 = ax2.get_legend_handles_labels()
    ax.legend(lines1+lines2, lab1+lab2, fontsize=8, loc="upper right")
    ax.grid(axis='y', alpha=0.3)

    # 4. LWR budget
    ax = axes[1][0]
    pitches = np.array([n.half_pitch_nm for n in PROCESS_NODES])
    lwr_specs = np.array([n.LWR_3sigma_nm for n in PROCESS_NODES])
    lwr_bcp  = 0.5*np.sqrt(pitches/25) * pitches**0.2
    lwr_guide = 0.7 * np.ones_like(pitches)  # EUV guide roughness
    lwr_process = np.sqrt(lwr_bcp**2 + lwr_guide**2)
    ax.bar(range(len(nodes)), lwr_specs, color='green', alpha=0.5, label="LWR spec (3σ)")
    ax.bar(range(len(nodes)), lwr_process, color='red', alpha=0.5, label="LWR predicted")
    ax.set_xticks(range(len(nodes)))
    ax.set_xticklabels([n.name.replace(" node","") for n in PROCESS_NODES],
                       rotation=25, fontsize=9)
    ax.set_ylabel("LWR 3σ (nm)", fontsize=11)
    ax.set_title("LWR Budget by Process Node", fontsize=12)
    ax.legend(fontsize=10); ax.grid(axis='y', alpha=0.3)

    # 5. Process window (7nm node)
    ax = axes[1][1]
    node7 = [n for n in PROCESS_NODES if n.name == "7nm node"][0]
    pwa = ProcessWindowAnalyzer(node7)
    pw_data = pwa.compute_process_window()
    T_arr = pw_data["T_arr"]
    t_arr = pw_data["t_arr"]
    PW = np.array(pw_data["process_window"])
    ax.contourf(T_arr, np.log10(t_arr), PW, levels=[-0.5, 0.5, 1.5],
                colors=["white", "#2ECC71"], alpha=0.6)
    ax.contour(T_arr, np.log10(t_arr), PW, levels=[0.5], colors="darkgreen", lw=2)
    ax.set_xlabel("Anneal Temperature (K)", fontsize=11)
    ax.set_ylabel("log₁₀(Anneal Time / s)", fontsize=11)
    ax.set_title(f"Process Window: 7nm Node\n(window fraction={pw_data['window_fraction']:.1%})",
                 fontsize=12)
    green_patch = mpatches.Patch(color="#2ECC71", alpha=0.6, label="Process Window")
    ax.legend(handles=[green_patch], fontsize=10)
    ax.grid(alpha=0.3)

    # 6. Defect density roadmap
    ax = axes[1][2]
    defects = [n.defect_density_cm2 for n in PROCESS_NODES]
    colors_d = plt.cm.RdYlGn(np.linspace(0.2, 0.8, len(defects)))[::-1]
    ax.bar(range(len(nodes)), np.log10(defects), color=colors_d, alpha=0.85)
    ax.set_xticks(range(len(nodes)))
    ax.set_xticklabels([n.name.replace(" node","") for n in PROCESS_NODES],
                       rotation=25, fontsize=9)
    ax.set_ylabel("log₁₀(Defect Density / cm⁻²)", fontsize=11)
    ax.set_title("Defect Density Roadmap (IRDS)", fontsize=12)
    for i, d in enumerate(defects):
        ax.text(i, np.log10(d)+0.1, f"{d:.0e}", ha='center', fontsize=8)
    ax.grid(axis='y', alpha=0.3)

    fig.suptitle("Semiconductor Process Integration: BCP DSA Lithography",
                 fontsize=14, fontweight='bold')
    fig.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_euv_dsa_flow(output_path="figures/euv_dsa_process_flow.png"):
    """Visualize EUV + DSA hybrid patterning flow."""
    flow = euv_dsa_hybrid_flow()
    steps = flow["process_flow"]

    fig, ax = plt.subplots(figsize=(14, 6))
    ax.set_xlim(0, 14); ax.set_ylim(0, 4)
    ax.axis('off')

    colors_step = ["#5DADE2","#E67E22","#2ECC71","#9B59B6",
                   "#E74C3C","#1ABC9C","#F39C12","#7F8C8D"]
    for i, (step, col) in enumerate(zip(steps, colors_step)):
        x = i * 1.7 + 0.5
        rect = mpatches.FancyBboxPatch((x, 0.5), 1.4, 3.0,
                                        boxstyle="round,pad=0.1",
                                        facecolor=col, alpha=0.8)
        ax.add_patch(rect)
        ax.text(x+0.7, 3.0, f"Step {step['step']}", ha='center',
                fontsize=8, fontweight='bold', color='white')
        process_short = step['process'][:12] + ("\n" + step['process'][12:20] if len(step['process'])>12 else "")
        ax.text(x+0.7, 2.2, process_short, ha='center', fontsize=7, color='white')
        if i < len(steps)-1:
            ax.annotate("", xy=(x+1.6, 2.0), xytext=(x+1.4, 2.0),
                        arrowprops=dict(arrowstyle="->", color="black", lw=2))

    ax.set_title(f"EUV + DSA Hybrid Patterning Flow (7nm Node, {flow['multiplication']})",
                 fontsize=13, fontweight='bold', pad=15)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=== Semiconductor Process Integration ===")

    print("\nProcess node summary:")
    process_summary = []
    for node in PROCESS_NODES:
        print(f"  {node.name:15s}  L0={node.L0_nm:.1f}nm  N={node.BCP_N}  "
              f"chi={node.BCP_chi:.3f}  n_DSA={node.DSA_n}  "
              f"defects={node.defect_density_cm2:.0e}/cm²")
        process_summary.append({
            "node": node.name, "half_pitch_nm": node.half_pitch_nm,
            "BCP_N": node.BCP_N, "chi": node.BCP_chi, "L0_nm": node.L0_nm,
            "DSA_n": node.DSA_n, "guide_pitch_nm": node.guide_pitch_nm,
            "LWR_spec_nm": node.LWR_3sigma_nm,
            "defect_density_cm2": node.defect_density_cm2, "note": node.note
        })

    print("\nHigh-χ BCP materials:")
    mat_data = []
    for name, data in HIGH_CHI_BCPS.items():
        N_min = estimate_chiN_for_ordering(data["chi_500K"])
        print(f"  {name:25s}  chi={data['chi_500K']:.3f}  "
              f"L0_min={data['L0_min_nm']:.1f}nm  N_min={N_min}  [{data['maturity']}]")
        mat_data.append({"material":name, **data, "N_min_ordering": N_min})

    # Process window for 7nm node
    node7 = [n for n in PROCESS_NODES if n.name == "7nm node"][0]
    pwa7  = ProcessWindowAnalyzer(node7)
    pw7   = pwa7.compute_process_window()
    print(f"\n7nm node process window fraction: {pw7['window_fraction']:.1%}")

    # EUV+DSA flow
    flow = euv_dsa_hybrid_flow()

    os.makedirs("results", exist_ok=True)
    with open("results/process_nodes.json", "w") as f:
        json.dump(process_summary, f, indent=2)
    with open("results/high_chi_materials.json", "w") as f:
        json.dump(mat_data, f, indent=2)
    with open("results/process_window_7nm.json", "w") as f:
        json.dump({k: v for k,v in pw7.items() if k != "defect_density"}, f, indent=2)
    with open("results/euv_dsa_flow.json", "w") as f:
        json.dump(flow, f, indent=2)

    plot_semiconductor_analysis()
    plot_euv_dsa_flow()

    print("\n✓ Semiconductor analysis complete.")
    print("  results/process_nodes.json, results/high_chi_materials.json")
    print("  results/process_window_7nm.json, results/euv_dsa_flow.json")
    print("  figures/semiconductor_process.png, figures/euv_dsa_process_flow.png")
