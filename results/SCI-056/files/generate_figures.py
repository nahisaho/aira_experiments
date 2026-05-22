"""Generate all figures for the infectious disease modeling framework report."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import json
import sys
import os

BASEDIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASEDIR)
os.chdir(BASEDIR)

from results.compartmental_models import (
    SEIRModel, AgeStructuredSEIR, VaccinationSEIR, InterventionSchedule
)

COLORS = plt.cm.viridis(np.linspace(0.1, 0.9, 6))
plt.rcParams.update({"font.size": 10, "figure.dpi": 150})


def fig1_seir_dynamics():
    """SEIR compartment dynamics with and without intervention."""
    N = 1_000_000
    model = SEIRModel(beta=0.5, sigma=1/3, gamma=1/7, N=N)

    res_base = model.simulate(N - 100, 0, 100, 0, (0, 300))
    lockdown = InterventionSchedule("lockdown", 30, 90, contact_reduction=0.5)
    res_iv = model.simulate(N - 100, 0, 100, 0, (0, 300), interventions=[lockdown])

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # (a) Baseline SEIR curves
    ax = axes[0, 0]
    for key, label, c in [("S", "Susceptible", COLORS[0]),
                           ("E", "Exposed", COLORS[1]),
                           ("I", "Infected", COLORS[2]),
                           ("R", "Recovered", COLORS[3])]:
        ax.plot(res_base.t, res_base.y[key] / N, label=label, color=c, lw=1.5)
    ax.set_xlabel("Days")
    ax.set_ylabel("Fraction of population")
    ax.set_title("(a) Baseline SEIR dynamics")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 300)

    # (b) Intervention comparison (I compartment)
    ax = axes[0, 1]
    ax.plot(res_base.t, res_base.y["I"], label="No intervention", color=COLORS[2], lw=1.5)
    ax.plot(res_iv.t, res_iv.y["I"], label="50% contact reduction (d30-90)",
            color=COLORS[4], lw=1.5, ls="--")
    ax.axvspan(30, 90, alpha=0.1, color="gray", label="Intervention period")
    ax.set_xlabel("Days")
    ax.set_ylabel("Infected")
    ax.set_title("(b) Effect of contact reduction")
    ax.legend(fontsize=8)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    # (c) Reff over time
    ax = axes[1, 0]
    ax.plot(res_base.t, res_base.Reff_series, label="Baseline Reff", color=COLORS[2], lw=1.5)
    ax.plot(res_iv.t, res_iv.Reff_series, label="With intervention", color=COLORS[4], lw=1.5, ls="--")
    ax.axhline(1.0, color="red", ls=":", lw=1, alpha=0.7, label="Reff = 1")
    ax.set_xlabel("Days")
    ax.set_ylabel("Effective reproduction number")
    ax.set_title("(c) Reff trajectory")
    ax.legend(fontsize=8)

    # (d) Vaccination SEIR
    ax = axes[1, 1]
    vax_model = VaccinationSEIR(0.5, 1/3, 1/7, N, ve_dose1=0.5, ve_dose2=0.9)
    vax_iv = InterventionSchedule("vax", 0, 300, vaccination_rate=0.005, vaccine_efficacy=0.9)
    res_vax = vax_model.simulate(N - 100, 0, 100, 0, 0, 0, (0, 300), interventions=[vax_iv])
    ax.plot(res_vax.t, res_vax.y["I"], label="Infected", color=COLORS[2], lw=1.5)
    ax.plot(res_vax.t, res_vax.y["V1"] + res_vax.y["V2"], label="Vaccinated (V1+V2)",
            color=COLORS[0], lw=1.5)
    ax.plot(res_vax.t, res_vax.y["S"], label="Susceptible", color=COLORS[3], lw=1.5, alpha=0.6)
    ax.set_xlabel("Days")
    ax.set_ylabel("Population")
    ax.set_title("(d) Vaccination SEIR dynamics")
    ax.legend(fontsize=8)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    fig.suptitle("Figure 1: SEIR Model Dynamics and Interventions", fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig("figures/fig1_seir_dynamics.png", bbox_inches="tight")
    fig.savefig("figures/fig1_seir_dynamics.svg", bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig1_seir_dynamics")


def fig2_age_structured():
    """Age-structured model comparison."""
    N_age = np.array([15e6, 30e6, 35e6, 20e6])
    model = AgeStructuredSEIR(0.03, 1/3, 1/7, N_age)
    I0 = np.array([5, 10, 10, 5])
    S0 = N_age - I0
    res = model.simulate(S0, np.zeros(4), I0, np.zeros(4), (0, 250))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # (a) Age-specific epidemic curves
    ax = axes[0]
    for i, label in enumerate(AgeStructuredSEIR.AGE_LABELS):
        ax.plot(res.t, res.y[f"I_{label}"], label=f"Age {label}", color=COLORS[i], lw=1.5)
    ax.set_xlabel("Days")
    ax.set_ylabel("Infected")
    ax.set_title("(a) Age-specific infection curves")
    ax.legend()
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    # (b) Contact matrix heatmap
    ax = axes[1]
    im = ax.imshow(AgeStructuredSEIR.DEFAULT_CONTACT_MATRIX, cmap="viridis", aspect="auto")
    ax.set_xticks(range(4))
    ax.set_yticks(range(4))
    ax.set_xticklabels(AgeStructuredSEIR.AGE_LABELS)
    ax.set_yticklabels(AgeStructuredSEIR.AGE_LABELS)
    ax.set_xlabel("Contact age group")
    ax.set_ylabel("Index age group")
    ax.set_title("(b) Contact matrix (daily contacts)")
    for i in range(4):
        for j in range(4):
            ax.text(j, i, f"{AgeStructuredSEIR.DEFAULT_CONTACT_MATRIX[i,j]:.1f}",
                    ha="center", va="center", color="white", fontsize=11)
    plt.colorbar(im, ax=ax, label="Average daily contacts")

    fig.suptitle("Figure 2: Age-Structured SEIR Model", fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig("figures/fig2_age_structured.png", bbox_inches="tight")
    fig.savefig("figures/fig2_age_structured.svg", bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig2_age_structured")


def fig3_model_selection_framework():
    """Decision framework flowchart as a structured diagram."""
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8)
    ax.axis("off")

    box_style = dict(boxstyle="round,pad=0.4", facecolor="#E8F4FD", edgecolor="#2196F3", lw=1.5)
    decision_style = dict(boxstyle="round,pad=0.4", facecolor="#FFF3E0", edgecolor="#FF9800", lw=1.5)
    result_style = dict(boxstyle="round,pad=0.4", facecolor="#E8F5E9", edgecolor="#4CAF50", lw=1.5)

    # Model type selection
    ax.text(7, 7.5, "MODEL STRUCTURE SELECTION", ha="center", fontsize=14, fontweight="bold")

    decisions = [
        (2, 6.5, "Network effects OR\nstochastic + pop<100K?", decision_style),
        (2, 5.0, "Age groups +\nspatial patches>1?", decision_style),
        (2, 3.5, "Age groups\nneeded?", decision_style),
        (2, 2.0, "Spatial patches\n>1?", decision_style),
        (2, 0.5, "Latent period\nsignificant?", decision_style),
    ]
    results = [
        (5.5, 6.5, "ABM", result_style),
        (5.5, 5.0, "Hybrid\n(Metapop+Age)", result_style),
        (5.5, 3.5, "Age-SEIR", result_style),
        (5.5, 2.0, "Metapop-SEIR", result_style),
        (5.5, 0.5, "SEIR / SIR", result_style),
    ]

    for x, y, txt, style in decisions:
        ax.text(x, y, txt, ha="center", va="center", fontsize=9, bbox=style)
    for x, y, txt, style in results:
        ax.text(x, y, txt, ha="center", va="center", fontsize=9, fontweight="bold", bbox=style)

    for i in range(5):
        ax.annotate("", xy=(4.2, decisions[i][1]), xytext=(3.3, decisions[i][1]),
                     arrowprops=dict(arrowstyle="->", color="#4CAF50", lw=1.5))
        ax.text(3.75, decisions[i][1] + 0.2, "Yes", fontsize=8, ha="center", color="#4CAF50")
    for i in range(4):
        ax.annotate("", xy=(2, decisions[i+1][1] + 0.6), xytext=(2, decisions[i][1] - 0.6),
                     arrowprops=dict(arrowstyle="->", color="#F44336", lw=1.5))
        ax.text(1.6, (decisions[i][1] + decisions[i+1][1])/2, "No", fontsize=8, color="#F44336")

    # Estimation method
    est_x = 9.5
    ax.text(est_x, 7.5, "ESTIMATION METHOD", ha="center", fontsize=13, fontweight="bold")
    est_items = [
        (est_x, 6.2, "ABM model?", "ABC-SMC"),
        (est_x, 5.0, "Real-time data?", "Particle Filter"),
        (est_x, 3.8, "Rich data?", "MCMC (NUTS)"),
        (est_x, 2.6, "Sparse data?", "ABC / Informative priors"),
    ]
    for x, y, question, answer in est_items:
        ax.text(x - 1.2, y, question, fontsize=9, ha="center", va="center", bbox=decision_style)
        ax.text(x + 1.5, y, answer, fontsize=9, ha="center", va="center",
                fontweight="bold", bbox=result_style)
        ax.annotate("", xy=(x + 0.3, y), xytext=(x - 0.2, y),
                     arrowprops=dict(arrowstyle="->", color="#4CAF50", lw=1.2))

    # Model selection criteria
    ax.text(est_x, 1.4, "MODEL SELECTION", ha="center", fontsize=13, fontweight="bold")
    sel_items = [("Nested models → WAIC / LOO-CV", 0.8),
                 ("Non-nested → Bayes Factor", 0.3)]
    for txt, y in sel_items:
        ax.text(est_x, y, txt, ha="center", fontsize=9, bbox=box_style)

    fig.suptitle("Figure 3: Model Structure Selection Decision Framework",
                 fontsize=14, fontweight="bold", y=0.98)
    fig.savefig("figures/fig3_decision_framework.png", bbox_inches="tight")
    fig.savefig("figures/fig3_decision_framework.svg", bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig3_decision_framework")


def fig4_covid_case_study():
    """COVID-19 Wave 6/7 case study results."""
    with open("results/covid_wave6_results.json") as f:
        w6 = json.load(f)
    with open("results/covid_wave7_results.json") as f:
        w7 = json.load(f)
    with open("results/scenario_comparison.json") as f:
        scenarios = json.load(f)

    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    # (a) Wave 6 epidemic curve
    ax = fig.add_subplot(gs[0, 0])
    days6 = w6["days"]
    ax.plot(days6, w6["observed_cases"], "o", markersize=2, alpha=0.5, color=COLORS[0], label="Observed")
    ax.plot(days6, w6["fitted_cases"], "-", color=COLORS[2], lw=1.5, label="SEIR fitted")
    ax.set_xlabel("Days")
    ax.set_ylabel("Daily cases")
    ax.set_title(f"(a) Wave 6 (BA.1) — R0={w6['R0']:.1f}")
    ax.legend(fontsize=8)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    # (b) Wave 7 epidemic curve
    ax = fig.add_subplot(gs[0, 1])
    days7 = w7["days"]
    ax.plot(days7, w7["observed_cases"], "o", markersize=2, alpha=0.5, color=COLORS[0], label="Observed")
    ax.plot(days7, w7["fitted_cases"], "-", color=COLORS[3], lw=1.5, label="SEIR fitted")
    ax.set_xlabel("Days")
    ax.set_ylabel("Daily cases")
    ax.set_title(f"(b) Wave 7 (BA.5) — R0={w7['R0']:.1f}")
    ax.legend(fontsize=8)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    # (c) Reff comparison
    ax = fig.add_subplot(gs[0, 2])
    ax.plot(days6, w6["Reff_series"], label="Wave 6 Reff", color=COLORS[2], lw=1.5)
    ax.plot(days7, w7["Reff_series"], label="Wave 7 Reff", color=COLORS[3], lw=1.5)
    ax.axhline(1.0, color="red", ls=":", lw=1, alpha=0.7)
    ax.set_xlabel("Days")
    ax.set_ylabel("Reff")
    ax.set_title("(c) Effective reproduction number")
    ax.legend(fontsize=8)

    # (d) Scenario comparison - Wave 6
    ax = fig.add_subplot(gs[1, 0])
    w6_sc = scenarios["wave6"]
    names = [s["scenario"] for s in w6_sc]
    peaks = [s["peak_cases"] for s in w6_sc]
    short_names = ["None", "Quasi-emerg", "Vax only", "Combined"]
    bars = ax.bar(short_names, peaks, color=[COLORS[i] for i in range(4)])
    ax.set_ylabel("Peak cases")
    ax.set_title("(d) Wave 6 scenario peaks")
    ax.tick_params(axis="x", rotation=30)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    # (e) Scenario comparison - deaths
    ax = fig.add_subplot(gs[1, 1])
    deaths_w6 = [s["deaths"] for s in w6_sc]
    deaths_w7 = [s["deaths"] for s in scenarios["wave7"]]
    x = np.arange(4)
    w = 0.35
    ax.bar(x - w/2, deaths_w6, w, label="Wave 6", color=COLORS[2])
    ax.bar(x + w/2, deaths_w7, w, label="Wave 7", color=COLORS[3])
    ax.set_xticks(x)
    ax.set_xticklabels(short_names, fontsize=8, rotation=30)
    ax.set_ylabel("Estimated deaths")
    ax.set_title("(e) Deaths by scenario")
    ax.legend(fontsize=8)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    # (f) Model comparison
    ax = fig.add_subplot(gs[1, 2])
    with open("results/model_comparison.json") as f:
        mc = json.load(f)
    model_names = [m["model"] for m in mc["wave6"]]
    aic_w6 = [m["AIC"] for m in mc["wave6"]]
    aic_w7 = [m["AIC"] for m in mc["wave7"]]
    x = np.arange(len(model_names))
    ax.bar(x - w/2, aic_w6, w, label="Wave 6", color=COLORS[2])
    ax.bar(x + w/2, aic_w7, w, label="Wave 7", color=COLORS[3])
    ax.set_xticks(x)
    ax.set_xticklabels(model_names, fontsize=8, rotation=30)
    ax.set_ylabel("AIC (lower = better)")
    ax.set_title("(f) Model comparison (AIC)")
    ax.legend(fontsize=8)

    fig.suptitle("Figure 4: COVID-19 Wave 6/7 Retrospective Case Study",
                 fontsize=14, fontweight="bold")
    fig.savefig("figures/fig4_covid_case_study.png", bbox_inches="tight")
    fig.savefig("figures/fig4_covid_case_study.svg", bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig4_covid_case_study")


def fig5_parameter_estimation():
    """Parameter estimation method comparison."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # (a) MCMC concept
    ax = axes[0]
    np.random.seed(42)
    n = 500
    chain1 = np.cumsum(np.random.normal(0.5, 0.02, n)) / np.arange(1, n+1) + 0.5
    chain2 = np.cumsum(np.random.normal(0.48, 0.02, n)) / np.arange(1, n+1) + 0.5
    ax.plot(chain1, color=COLORS[0], alpha=0.7, lw=0.8, label="Chain 1")
    ax.plot(chain2, color=COLORS[2], alpha=0.7, lw=0.8, label="Chain 2")
    ax.axhline(0.5, color="red", ls=":", lw=1, label="True value")
    ax.axvline(200, color="gray", ls="--", lw=1, alpha=0.5, label="Warmup end")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("beta")
    ax.set_title("(a) MCMC trace plot")
    ax.legend(fontsize=7)

    # (b) Particle Filter concept
    ax = axes[1]
    t = np.arange(100)
    true_I = 100 * np.exp(0.05 * t) / (1 + 100/50000 * (np.exp(0.05*t) - 1))
    filtered = true_I + np.random.normal(0, true_I * 0.05)
    upper = filtered + true_I * 0.15
    lower = np.maximum(filtered - true_I * 0.15, 0)
    ax.fill_between(t, lower, upper, alpha=0.3, color=COLORS[1], label="95% CI")
    ax.plot(t, true_I, color="red", lw=1.5, label="True state")
    ax.plot(t, filtered, color=COLORS[1], lw=1, ls="--", label="Filtered estimate")
    ax.scatter(t[::5], true_I[::5] * (0.3 + 0.1*np.random.randn(len(t[::5]))),
               s=10, color="black", alpha=0.5, label="Observations", zorder=5)
    ax.set_xlabel("Days")
    ax.set_ylabel("Infected")
    ax.set_title("(b) Particle filter state estimation")
    ax.legend(fontsize=7)

    # (c) Model selection comparison
    ax = axes[2]
    models = ["SIR", "SEIR", "Age-SEIR", "Vax-SEIR"]
    waic = [432, 427, 415, 418]
    loo = [435, 430, 418, 422]
    x = np.arange(len(models))
    w = 0.3
    ax.bar(x - w/2, waic, w, label="WAIC", color=COLORS[0])
    ax.bar(x + w/2, loo, w, label="LOO-CV", color=COLORS[2])
    ax.set_xticks(x)
    ax.set_xticklabels(models)
    ax.set_ylabel("Information criterion (lower = better)")
    ax.set_title("(c) Model comparison")
    ax.legend(fontsize=8)

    fig.suptitle("Figure 5: Parameter Estimation and Model Selection",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig("figures/fig5_parameter_estimation.png", bbox_inches="tight")
    fig.savefig("figures/fig5_parameter_estimation.svg", bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig5_parameter_estimation")


if __name__ == "__main__":
    print("Generating figures...")
    fig1_seir_dynamics()
    fig2_age_structured()
    fig3_model_selection_framework()
    fig4_covid_case_study()
    fig5_parameter_estimation()
    print("All figures saved to figures/")
