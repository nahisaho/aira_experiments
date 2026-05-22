"""
Module 4: Agent-Based Model for Research Career Path Prediction
Pure-Python ABM (Mesa-compatible design) simulating researcher career
trajectories under different funding regimes.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import json
import os

np.random.seed(42)


class ResearcherAgent:
    """Agent representing a researcher in the funding ecosystem."""

    def __init__(self, agent_id, gender, region, field, career_stage):
        self.agent_id = agent_id
        self.gender = gender
        self.region = region
        self.field = field
        self.career_stage = career_stage
        self.career_years = {"early": 3, "mid": 12, "senior": 25}[career_stage]

        self.talent = np.random.beta(2, 5)
        self.productivity = max(1, int(np.random.lognormal(1.0, 0.5)))
        self.h_index = max(0, int(self.talent * self.career_years * 1.5))
        self.total_citations = max(0, int(self.h_index ** 2 * np.random.lognormal(0, 0.5)))
        self.reputation = self.talent * 0.5 + (self.h_index / 50) * 0.5

        self.funded = False
        self.cumulative_funding = 0
        self.funding_history = []
        self.papers_this_year = 0
        self.left_academia = False

    def step(self):
        if self.left_academia:
            return

        self.career_years += 1
        if self.career_years < 7:
            self.career_stage = "early"
        elif self.career_years < 20:
            self.career_stage = "mid"
        else:
            self.career_stage = "senior"

        base_output = self.talent * self.productivity * 0.5
        funding_boost = 1.5 if self.funded else 0.7
        noise = np.random.normal(1.0, 0.2)
        self.papers_this_year = max(0, int(base_output * funding_boost * noise))

        new_citations = int(self.papers_this_year * np.random.exponential(3))
        self.total_citations += new_citations
        self.h_index = min(self.h_index + max(0, self.papers_this_year // 3),
                           int(np.sqrt(self.total_citations) * 0.6))
        self.reputation = np.clip(
            self.reputation * 0.95 + 0.05 * (self.papers_this_year / 5), 0, 1
        )

        if self.career_stage == "early" and not self.funded:
            if self.career_years > 5 and np.random.random() < 0.08:
                self.left_academia = True

        self.funding_history.append(self.funded)
        self.funded = False


class FundingModel:
    """ABM model simulating research funding ecosystem."""

    def __init__(self, n_researchers=200, mechanism="peer_review",
                 budget=1_000_000, grant_size=50_000, n_steps=20):
        self.n_researchers = n_researchers
        self.mechanism = mechanism
        self.budget = budget
        self.grant_size = grant_size
        self.n_grants = int(budget / grant_size)
        self.n_steps = n_steps

        fields = ["Physics", "Biology", "CS", "Chemistry", "Medicine"]
        genders = ["M", "F", "Other"]
        regions = ["Asia", "Europe", "NorthAmerica", "LatinAmerica", "Africa"]
        stages = ["early", "mid", "senior"]

        self.agents = []
        for i in range(n_researchers):
            agent = ResearcherAgent(
                i,
                gender=np.random.choice(genders, p=[0.60, 0.35, 0.05]),
                region=np.random.choice(regions, p=[0.30, 0.30, 0.25, 0.10, 0.05]),
                field=np.random.choice(fields),
                career_stage=np.random.choice(stages, p=[0.40, 0.35, 0.25]),
            )
            self.agents.append(agent)

        self.history = []

    def step(self):
        active = [a for a in self.agents if not a.left_academia]
        self._allocate_funding(active)

        order = list(range(len(self.agents)))
        np.random.shuffle(order)
        for i in order:
            self.agents[i].step()

        self.history.append(self._collect_metrics())

    def _allocate_funding(self, active_agents):
        if self.mechanism == "peer_review":
            scored = sorted(active_agents,
                            key=lambda a: a.reputation + np.random.normal(0, 0.1),
                            reverse=True)
            for a in scored[:self.n_grants]:
                a.funded = True
                a.cumulative_funding += self.grant_size

        elif self.mechanism == "lottery":
            median_rep = np.median([a.reputation for a in active_agents])
            qualified = [a for a in active_agents if a.reputation >= median_rep * 0.6]
            if len(qualified) > self.n_grants:
                idx = np.random.choice(len(qualified), self.n_grants, replace=False)
                selected = [qualified[i] for i in idx]
            else:
                selected = qualified
            for a in selected:
                a.funded = True
                a.cumulative_funding += self.grant_size

        elif self.mechanism == "automated":
            for a in active_agents:
                a._auto_score = a.reputation * 0.6
                if a.career_stage == "early":
                    a._auto_score += 0.15
                if a.gender == "F":
                    a._auto_score += 0.05
                if a.region in ["LatinAmerica", "Africa"]:
                    a._auto_score += 0.1
            scored = sorted(active_agents, key=lambda a: a._auto_score, reverse=True)
            for a in scored[:self.n_grants]:
                a.funded = True
                a.cumulative_funding += self.grant_size

    def _collect_metrics(self):
        active = [a for a in self.agents if not a.left_academia]
        if not active:
            return {}
        m_rep = [a.reputation for a in active if a.gender == "M"]
        f_rep = [a.reputation for a in active if a.gender == "F"]
        funding = sorted([a.cumulative_funding for a in active])
        n = len(funding)
        total_f = sum(funding)
        if total_f > 0:
            idx = np.arange(1, n + 1)
            gini = (2 * np.sum(idx * np.array(funding)) - (n + 1) * total_f) / (n * total_f)
        else:
            gini = 0
        return {
            "TotalPapers": sum(a.papers_this_year for a in active),
            "AvgHIndex": np.mean([a.h_index for a in active]),
            "ActiveResearchers": len(active),
            "EarlyCareerPct": sum(1 for a in active if a.career_stage == "early") / len(active) * 100,
            "GenderGap": (np.mean(m_rep) - np.mean(f_rep)) if m_rep and f_rep else 0,
            "FundingGini": gini,
            "AttritionRate": sum(1 for a in self.agents if a.left_academia) / self.n_researchers,
        }

    def get_model_dataframe(self):
        return pd.DataFrame(self.history)


def run_abm_simulation(n_researchers=200, n_steps=20,
                       output_dir="figures", results_dir="results"):
    """Run ABM simulations for all mechanisms and compare."""
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    mechanisms = ["peer_review", "lottery", "automated"]
    all_data = {}
    final_stats = {}

    for mech in mechanisms:
        np.random.seed(42)
        model = FundingModel(n_researchers=n_researchers, mechanism=mech, n_steps=n_steps)
        for _ in range(n_steps):
            model.step()

        model_data = model.get_model_dataframe()
        all_data[mech] = model_data

        final_stats[mech] = {
            "final_total_papers": int(model_data["TotalPapers"].iloc[-1]),
            "final_avg_hindex": round(float(model_data["AvgHIndex"].iloc[-1]), 2),
            "final_active": int(model_data["ActiveResearchers"].iloc[-1]),
            "final_attrition": round(float(model_data["AttritionRate"].iloc[-1]), 3),
            "final_funding_gini": round(float(model_data["FundingGini"].iloc[-1]), 3),
            "final_gender_gap": round(float(model_data["GenderGap"].iloc[-1]), 4),
            "cumulative_papers": int(model_data["TotalPapers"].sum()),
        }

    _plot_abm_results(all_data, output_dir)

    with open(f"{results_dir}/abm_results.json", "w") as f:
        json.dump(final_stats, f, indent=2)

    print("ABM simulation complete.")
    for mech, stats in final_stats.items():
        print(f"  {mech}: papers={stats['cumulative_papers']}, "
              f"gini={stats['final_funding_gini']:.3f}, "
              f"attrition={stats['final_attrition']:.3f}")

    return final_stats


def _plot_abm_results(all_data, output_dir):
    """Plot ABM simulation results."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    colors = {"peer_review": "#2196F3", "lottery": "#4CAF50", "automated": "#FF9800"}
    labels = {"peer_review": "Peer Review", "lottery": "Lottery", "automated": "Automated"}

    metrics = ["TotalPapers", "AvgHIndex", "ActiveResearchers",
               "FundingGini", "GenderGap", "AttritionRate"]
    titles = ["Total Papers per Year", "Average h-index",
              "Active Researchers", "Funding Gini Coefficient",
              "Gender Gap (M-F Reputation)", "Attrition Rate"]

    for idx, (metric, title) in enumerate(zip(metrics, titles)):
        ax = axes[idx // 3][idx % 3]
        for mech, data in all_data.items():
            ax.plot(data.index + 1, data[metric], label=labels[mech],
                    color=colors[mech], linewidth=2)
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("Year", fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(f"{output_dir}/fig7_abm_career_simulation.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    run_abm_simulation()
