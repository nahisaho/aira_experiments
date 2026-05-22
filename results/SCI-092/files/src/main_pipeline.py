#!/usr/bin/env python3
"""
Main pipeline: Social Acceptance Prediction Model for Emerging Technologies
Integrates all 6 components and generates report artifacts.
"""

import sys, os, json, datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from meta_analysis import run_meta_analysis
from sentiment_analysis import run_sentiment_analysis
from psychometric_risk import run_psychometric_analysis
from framing_effects import run_framing_analysis
from sem_model import run_sem_analysis
from japan_case_study import run_japan_case_study
from visualizations import generate_all_figures

def main():
    os.makedirs("data", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    os.makedirs("figures", exist_ok=True)
    os.makedirs("logs", exist_ok=True)

    log_entries = []
    def log(phase, event, **kwargs):
        entry = {"timestamp": datetime.datetime.now().isoformat(),
                 "phase": phase, "event_type": event,
                 "actor": "co-scientist", **kwargs}
        log_entries.append(entry)

    log("pipeline", "run_started")

    # Component 1: Meta-analysis
    print("=" * 60)
    print("Component 1: Meta-Analysis of Public Opinion Surveys")
    print("=" * 60)
    log("meta_analysis", "skill_selected")
    meta_df, meta_results = run_meta_analysis()
    for tech in ["gene_editing", "AI", "nuclear_fusion"]:
        r = meta_results[tech]
        print(f"  {tech}: pooled={r['pooled_effect']:.3f} "
              f"95%CI=[{r['ci_95'][0]:.3f}, {r['ci_95'][1]:.3f}] "
              f"I²={r['I2']:.1f}%")
    log("meta_analysis", "handoff_completed", files_written=["data/meta_analysis_studies.csv", "results/meta_analysis_results.json"])

    # Component 2: Sentiment Analysis
    print("\n" + "=" * 60)
    print("Component 2: Social Media Sentiment Analysis (BERT/Lexicon Hybrid)")
    print("=" * 60)
    log("sentiment", "skill_selected")
    sent_df, sent_results = run_sentiment_analysis()
    for tech in ["gene_editing", "AI", "nuclear_fusion"]:
        r = sent_results[tech]
        print(f"  {tech}: mean_hybrid={r['hybrid_score_stats']['mean']:.3f}, n={r['n_posts']}")
    log("sentiment", "handoff_completed", files_written=["data/social_media_sentiment.csv", "results/sentiment_analysis_results.json"])

    # Component 3: Psychometric Risk Analysis
    print("\n" + "=" * 60)
    print("Component 3: Psychometric Paradigm Risk Perception Model")
    print("=" * 60)
    log("psychometric", "skill_selected")
    psych_df, psych_loadings, psych_results = run_psychometric_analysis()
    for tech in ["gene_editing", "AI", "nuclear_fusion"]:
        r = psych_results[tech]
        print(f"  {tech}: dread={r['dread_factor_mean']:.3f}, unknown={r['unknown_factor_mean']:.3f}, "
              f"risk-accept r={r['risk_acceptance_correlation']:.3f}")
    log("psychometric", "handoff_completed", files_written=["data/psychometric_risk_data.csv", "results/psychometric_results.json"])

    # Component 4: Framing Effects
    print("\n" + "=" * 60)
    print("Component 4: Framing Effects Quantitative Evaluation")
    print("=" * 60)
    log("framing", "skill_selected")
    framing_df, framing_results = run_framing_analysis()
    for tech in ["gene_editing", "AI", "nuclear_fusion"]:
        r = framing_results[tech]["anova"]
        print(f"  {tech}: F={r['F_statistic']:.2f}, p={r['p_value']:.4f}, "
              f"η²={r['eta_squared']:.4f} ({r['effect_interpretation']})")
    log("framing", "handoff_completed", files_written=["data/framing_experiment_data.csv", "results/framing_effects_results.json"])

    # Component 5: SEM Path Analysis
    print("\n" + "=" * 60)
    print("Component 5: Trust-Acceptance SEM Path Analysis")
    print("=" * 60)
    log("sem", "skill_selected")
    sem_data, sem_model, sem_results = run_sem_analysis()
    for p in sem_results["structural_paths"]:
        sig = "***" if p["p_value"] < 0.001 else "**" if p["p_value"] < 0.01 else "*" if p["p_value"] < 0.05 else "ns"
        print(f"  {p['from']} → {p['to']}: β={p['estimate']:.3f} ({sig})")
    print(f"  Model Fit: CFI={sem_results['fit_indices'].get('CFI', 'N/A')}, "
          f"RMSEA={sem_results['fit_indices'].get('RMSEA', 'N/A')}")
    log("sem", "handoff_completed", files_written=["data/sem_data.csv", "results/sem_results.json"])

    # Component 6: Japan Case Study
    print("\n" + "=" * 60)
    print("Component 6: Japan Genome-Edited Food Case Study")
    print("=" * 60)
    log("japan_case", "skill_selected")
    japan_df, japan_results = run_japan_case_study()
    o = japan_results["overall"]
    print(f"  Acceptance: M={o['acceptance_mean']:.2f} (SD={o['acceptance_std']:.2f})")
    print(f"  % Accepting (≥4): {o['pct_accepting']:.1f}%")
    print(f"  % Rejecting (<3): {o['pct_rejecting']:.1f}%")
    print(f"  Purchase intention: M={o['purchase_intention_mean']:.2f}")
    r = japan_results["regression"]
    print(f"  Regression R²={r['r_squared']:.3f}")
    log("japan_case", "handoff_completed", files_written=["data/japan_genome_food_data.csv", "results/japan_case_study_results.json"])

    # Generate all figures
    print("\n" + "=" * 60)
    print("Generating Figures")
    print("=" * 60)
    log("visualization", "skill_selected")
    generate_all_figures(meta_df, meta_results, sent_df, sent_results,
                         psych_df, psych_results, framing_df, framing_results,
                         sem_results, japan_df, japan_results)
    log("visualization", "handoff_completed",
        files_written=["figures/fig1_meta_analysis_forest.png", "figures/fig2_sentiment_trends.png",
                        "figures/fig3_psychometric_space.png", "figures/fig4_framing_effects.png",
                        "figures/fig5_sem_path_diagram.png", "figures/fig6_japan_case_study.png",
                        "figures/fig7_correlation_heatmap.png"])

    log("pipeline", "run_completed", status="ok")

    with open("logs/process-log.jsonl", "w") as f:
        for entry in log_entries:
            f.write(json.dumps(entry, ensure_ascii=False) + "\n")

    print("\n" + "=" * 60)
    print("Pipeline completed successfully!")
    print("=" * 60)

    return {
        "meta": (meta_df, meta_results),
        "sentiment": (sent_df, sent_results),
        "psychometric": (psych_df, psych_results),
        "framing": (framing_df, framing_results),
        "sem": sem_results,
        "japan": (japan_df, japan_results)
    }


if __name__ == "__main__":
    main()
