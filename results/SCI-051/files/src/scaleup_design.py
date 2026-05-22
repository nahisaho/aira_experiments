"""
Scale-Up Design: Numbering Up vs Scaling Up
=============================================
Comparative analysis of microreactor scale-up strategies.
"""

import numpy as np
import json, os

def numbering_up_analysis(single_reactor, n_parallel_options):
    """Analyze numbering-up (parallelization) strategy."""
    results = []
    base_throughput = single_reactor["throughput_g_h"]
    base_cost = single_reactor["reactor_cost_usd"]

    for n in n_parallel_options:
        throughput = base_throughput * n
        # Flow distribution challenge increases with n
        flow_uniformity = max(0.85, 1.0 - 0.005 * n)
        effective_throughput = throughput * flow_uniformity

        # Cost model: reactor + manifold + controls
        reactor_cost = base_cost * n * (1 - 0.05 * np.log2(max(n, 1)))  # volume discount
        manifold_cost = 500 * np.sqrt(n)
        control_cost = 2000 + 200 * n
        total_cost = reactor_cost + manifold_cost + control_cost

        # Heat transfer maintained (same channel dimensions)
        heat_transfer_coeff = single_reactor["heat_transfer_W_m2K"]

        results.append({
            "n_parallel": n,
            "throughput_g_h": round(effective_throughput, 2),
            "flow_uniformity": round(flow_uniformity, 3),
            "total_cost_usd": round(total_cost, 0),
            "cost_per_g_h": round(total_cost / effective_throughput, 2),
            "heat_transfer_W_m2K": heat_transfer_coeff,
            "mixing_time_ms": single_reactor["mixing_time_ms"],
            "pressure_drop_kPa": single_reactor["pressure_drop_kPa"],
        })

    return results

def scaling_up_analysis(single_reactor, scale_factors):
    """Analyze traditional scaling-up strategy."""
    results = []
    base_dh = single_reactor["hydraulic_diameter_um"]
    base_throughput = single_reactor["throughput_g_h"]
    base_cost = single_reactor["reactor_cost_usd"]

    for sf in scale_factors:
        dh_new = base_dh * sf
        throughput = base_throughput * sf**2

        # Heat transfer degrades with larger channels
        ht_ratio = 1.0 / sf  # surface-to-volume ratio decreases
        heat_transfer = single_reactor["heat_transfer_W_m2K"] * ht_ratio

        # Mixing time increases (diffusion limited)
        mixing_time = single_reactor["mixing_time_ms"] * sf**2

        # Reynolds number increases
        Re_new = single_reactor["reynolds_number"] * sf

        # Pressure drop changes
        dP = single_reactor["pressure_drop_kPa"] / sf

        # Cost
        reactor_cost = base_cost * sf**1.5
        total_cost = reactor_cost + 2000  # fixed control cost

        results.append({
            "scale_factor": sf,
            "hydraulic_diameter_um": round(dh_new, 0),
            "throughput_g_h": round(throughput, 2),
            "heat_transfer_W_m2K": round(heat_transfer, 0),
            "mixing_time_ms": round(mixing_time, 1),
            "reynolds_number": round(Re_new, 1),
            "pressure_drop_kPa": round(dP, 2),
            "total_cost_usd": round(total_cost, 0),
            "cost_per_g_h": round(total_cost / throughput, 2),
        })

    return results

def production_target_analysis(target_kg_per_day, single_reactor):
    """Determine optimal strategy to meet production targets."""
    target_g_h = target_kg_per_day * 1000 / 24
    base = single_reactor["throughput_g_h"]

    # Numbering up
    n_needed = int(np.ceil(target_g_h / base))
    nu_results = numbering_up_analysis(single_reactor, [n_needed])

    # Scaling up
    sf_needed = np.sqrt(target_g_h / base)
    su_results = scaling_up_analysis(single_reactor, [sf_needed])

    # Hybrid: moderate scale-up + numbering up
    sf_hybrid = min(3.0, np.sqrt(sf_needed))
    throughput_per_scaled = base * sf_hybrid**2
    n_hybrid = int(np.ceil(target_g_h / throughput_per_scaled))
    hybrid_reactor = single_reactor.copy()
    hybrid_reactor["throughput_g_h"] = throughput_per_scaled
    hybrid_reactor["hydraulic_diameter_um"] = single_reactor["hydraulic_diameter_um"] * sf_hybrid
    hybrid_reactor["heat_transfer_W_m2K"] = single_reactor["heat_transfer_W_m2K"] / sf_hybrid
    hybrid_reactor["mixing_time_ms"] = single_reactor["mixing_time_ms"] * sf_hybrid**2
    hy_results = numbering_up_analysis(hybrid_reactor, [n_hybrid])

    return {
        "production_target_kg_day": target_kg_per_day,
        "target_throughput_g_h": round(target_g_h, 1),
        "numbering_up": {
            "n_reactors": n_needed,
            "total_cost_usd": nu_results[0]["total_cost_usd"],
            "heat_transfer_preserved": True,
            "mixing_preserved": True,
        },
        "scaling_up": {
            "scale_factor": round(sf_needed, 2),
            "total_cost_usd": su_results[0]["total_cost_usd"],
            "heat_transfer_W_m2K": su_results[0]["heat_transfer_W_m2K"],
            "mixing_time_ms": su_results[0]["mixing_time_ms"],
        },
        "hybrid": {
            "scale_factor": round(sf_hybrid, 2),
            "n_reactors": n_hybrid,
            "total_cost_usd": hy_results[0]["total_cost_usd"],
            "heat_transfer_factor": round(1.0 / sf_hybrid, 3),
        },
        "recommendation": "hybrid" if sf_needed > 5 else "numbering_up",
    }

def run_scaleup_analysis():
    single_reactor = {
        "hydraulic_diameter_um": 285.7,
        "throughput_g_h": 2.5,
        "reactor_cost_usd": 5000,
        "heat_transfer_W_m2K": 5000,
        "mixing_time_ms": 50,
        "pressure_drop_kPa": 15.0,
        "reynolds_number": 4.76,
    }

    nu_results = numbering_up_analysis(single_reactor, [1, 2, 5, 10, 20, 50])
    su_results = scaling_up_analysis(single_reactor, [1, 2, 3, 5, 10])

    targets = [0.1, 1.0, 10.0, 100.0]  # kg/day
    target_analyses = [production_target_analysis(t, single_reactor) for t in targets]

    results = {
        "single_reactor_baseline": single_reactor,
        "numbering_up": nu_results,
        "scaling_up": su_results,
        "production_targets": target_analyses,
    }

    return results

if __name__ == "__main__":
    results = run_scaleup_analysis()
    os.makedirs("results", exist_ok=True)

    with open("results/scaleup_results.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print("=== Scale-Up Analysis ===")
    print("\n--- Numbering Up ---")
    print(f"{'N':<6} {'Throughput':<12} {'Cost/unit':<12} {'Uniformity':<10}")
    for r in results["numbering_up"]:
        print(f"{r['n_parallel']:<6} {r['throughput_g_h']:<12.1f} "
              f"{r['cost_per_g_h']:<12.1f} {r['flow_uniformity']:<10.3f}")

    print("\n--- Scaling Up ---")
    print(f"{'SF':<6} {'Dh (μm)':<10} {'Throughput':<12} {'HT (W/m²K)':<12} {'Mix (ms)':<10}")
    for r in results["scaling_up"]:
        print(f"{r['scale_factor']:<6} {r['hydraulic_diameter_um']:<10.0f} "
              f"{r['throughput_g_h']:<12.1f} {r['heat_transfer_W_m2K']:<12.0f} "
              f"{r['mixing_time_ms']:<10.1f}")

    print("\n--- Production Targets ---")
    for t in results["production_targets"]:
        print(f"\n  Target: {t['production_target_kg_day']} kg/day → "
              f"Recommended: {t['recommendation']}")
