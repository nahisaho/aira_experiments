"""
Agent-Based Model for Research Funding Allocation Optimization.
Simulates peer review, lottery, and hybrid funding mechanisms,
evaluating efficiency, fairness, and diversity outcomes.
"""

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from dataclasses import dataclass, field
from typing import List, Dict, Tuple
from scipy import stats
from pathlib import Path
import json
import warnings
warnings.filterwarnings("ignore")

np.random.seed(42)
FIGURES_DIR = Path("figures")
FIGURES_DIR.mkdir(exist_ok=True)

# ── Configuration ──────────────────────────────────────────────────────────
@dataclass
class SimConfig:
    n_researchers: int = 200
    n_rounds: int = 30
    budget_per_round: float = 100.0
    grant_size: float = 10.0
    n_reviewers_per_proposal: int = 3
    reviewer_noise_std: float = 0.3
    reviewer_bias_std: float = 0.15
    gender_ratio_female: float = 0.35
    n_regions: int = 5
    n_fields: int = 6
    collaboration_prob: float = 0.05
    citation_decay: float = 0.9
    matthew_effect_strength: float = 0.1
    # KAKENHI-specific
    kakenhi_categories: List[str] = field(default_factory=lambda: [
        "S", "A", "B", "C", "Early_Career", "Challenging"
    ])


# ── Researcher Agent ───────────────────────────────────────────────────────
class Researcher:
    def __init__(self, uid: int, config: SimConfig):
        self.uid = uid
        self.intrinsic_quality = np.random.beta(2, 5)  # right-skewed
        self.productivity = np.random.gamma(2, 0.5)
        self.gender = "F" if np.random.random() < config.gender_ratio_female else "M"
        self.region = np.random.randint(0, config.n_regions)
        self.field = np.random.randint(0, config.n_fields)
        self.career_stage = np.random.choice(
            ["early", "mid", "senior"], p=[0.4, 0.35, 0.25]
        )
        self.cumulative_funding = 0.0
        self.publications = 0
        self.citations = 0
        self.h_index = 0
        self.funded_rounds = 0
        self.total_proposals = 0
        self.collaborators: set = set()
        self.publication_history: List[float] = []
        self.funding_history: List[bool] = []

    @property
    def current_quality(self):
        """Quality improves with funding (Matthew effect) but has diminishing returns."""
        funding_bonus = 0.05 * np.log1p(self.cumulative_funding)
        experience_bonus = 0.02 * np.log1p(self.publications)
        return min(1.0, self.intrinsic_quality + funding_bonus + experience_bonus)

    def produce_research(self, funded: bool, config: SimConfig):
        base_output = self.productivity * self.current_quality
        if funded:
            output = base_output * (1.5 + 0.5 * np.random.random())
            n_pubs = max(1, int(np.random.poisson(output * 3)))
        else:
            output = base_output * (0.3 + 0.2 * np.random.random())
            n_pubs = max(0, int(np.random.poisson(output * 1)))
        
        self.publications += n_pubs
        new_citations = sum(
            max(0, int(np.random.poisson(self.current_quality * 5 + 1)))
            for _ in range(n_pubs)
        )
        self.citations += new_citations
        self.publication_history.append(n_pubs)
        self._update_h_index()
        return output, n_pubs, new_citations

    def _update_h_index(self):
        """Approximate h-index from citation count and publication count."""
        if self.publications == 0:
            self.h_index = 0
            return
        avg_cit = self.citations / max(1, self.publications)
        self.h_index = int(min(self.publications, np.sqrt(self.citations)))

    def proposal_quality(self, config: SimConfig):
        """Generate a proposal with noise."""
        base = self.current_quality
        noise = np.random.normal(0, 0.1)
        return np.clip(base + noise, 0, 1)


# ── Funding Mechanisms ─────────────────────────────────────────────────────
def peer_review_allocation(researchers: List[Researcher], budget: float,
                           grant_size: float, config: SimConfig) -> List[int]:
    """Traditional peer review with reviewer noise and bias."""
    proposals = []
    for r in researchers:
        q = r.proposal_quality(config)
        # Reviewer scores with noise and potential bias
        scores = []
        for _ in range(config.n_reviewers_per_proposal):
            bias = np.random.normal(0, config.reviewer_bias_std)
            # Gender bias (small systematic)
            gender_bias = -0.03 if r.gender == "F" else 0.0
            # Matthew effect: known researchers get slight boost
            reputation_boost = config.matthew_effect_strength * np.log1p(r.h_index)
            noise = np.random.normal(0, config.reviewer_noise_std)
            score = q + bias + gender_bias + reputation_boost + noise
            scores.append(np.clip(score, 0, 1))
        avg_score = np.mean(scores)
        proposals.append((r.uid, avg_score))
    
    proposals.sort(key=lambda x: x[1], reverse=True)
    n_grants = int(budget / grant_size)
    return [p[0] for p in proposals[:n_grants]]


def lottery_allocation(researchers: List[Researcher], budget: float,
                       grant_size: float, config: SimConfig) -> List[int]:
    """Pure random lottery among eligible applicants."""
    n_grants = int(budget / grant_size)
    selected = np.random.choice(len(researchers), size=min(n_grants, len(researchers)),
                                replace=False)
    return [researchers[i].uid for i in selected]


def modified_lottery_allocation(researchers: List[Researcher], budget: float,
                                grant_size: float, config: SimConfig,
                                shortlist_ratio: float = 0.5) -> List[int]:
    """Two-stage: peer review shortlist, then lottery among shortlisted."""
    proposals = []
    for r in researchers:
        q = r.proposal_quality(config)
        scores = []
        for _ in range(config.n_reviewers_per_proposal):
            noise = np.random.normal(0, config.reviewer_noise_std)
            score = q + noise
            scores.append(np.clip(score, 0, 1))
        proposals.append((r.uid, np.mean(scores)))
    
    proposals.sort(key=lambda x: x[1], reverse=True)
    shortlist_size = int(len(proposals) * shortlist_ratio)
    shortlisted = proposals[:shortlist_size]
    
    n_grants = int(budget / grant_size)
    if len(shortlisted) <= n_grants:
        return [p[0] for p in shortlisted]
    
    selected_idx = np.random.choice(len(shortlisted), size=n_grants, replace=False)
    return [shortlisted[i][0] for i in selected_idx]


def diversity_constrained_allocation(researchers: List[Researcher], budget: float,
                                     grant_size: float, config: SimConfig,
                                     gender_target: float = 0.4,
                                     region_balance: float = 0.15) -> List[int]:
    """Peer review with diversity constraints (quotas)."""
    proposals = []
    for r in researchers:
        q = r.proposal_quality(config)
        scores = []
        for _ in range(config.n_reviewers_per_proposal):
            noise = np.random.normal(0, config.reviewer_noise_std)
            score = q + noise
            scores.append(np.clip(score, 0, 1))
        proposals.append((r.uid, np.mean(scores), r))
    
    proposals.sort(key=lambda x: x[1], reverse=True)
    n_grants = int(budget / grant_size)
    
    selected = []
    female_count = 0
    region_counts = {i: 0 for i in range(config.n_regions)}
    max_per_region = max(1, int(n_grants * region_balance * 2))
    
    for uid, score, r in proposals:
        if len(selected) >= n_grants:
            break
        if r.gender == "F":
            female_count += 1
        if region_counts.get(r.region, 0) < max_per_region:
            region_counts[r.region] = region_counts.get(r.region, 0) + 1
            selected.append(uid)
    
    # If under-target for female, do second pass
    if female_count / max(1, len(selected)) < gender_target:
        remaining = [p for p in proposals if p[0] not in selected and p[2].gender == "F"]
        for uid, score, r in remaining:
            if len(selected) >= n_grants:
                break
            selected.append(uid)
    
    return selected[:n_grants]


# ── Network Analysis ───────────────────────────────────────────────────────
def build_collaboration_network(researchers: List[Researcher],
                                config: SimConfig) -> nx.Graph:
    """Build co-authorship network with preferential attachment."""
    G = nx.Graph()
    for r in researchers:
        G.add_node(r.uid, gender=r.gender, region=r.region,
                   field=r.field, career_stage=r.career_stage)
    
    # Initial random connections
    for i, r1 in enumerate(researchers):
        for j, r2 in enumerate(researchers):
            if i >= j:
                continue
            # Same field increases probability
            field_factor = 3.0 if r1.field == r2.field else 0.5
            # Same region increases probability
            region_factor = 2.0 if r1.region == r2.region else 0.8
            prob = config.collaboration_prob * field_factor * region_factor
            if np.random.random() < prob:
                G.add_edge(r1.uid, r2.uid, weight=1)
                r1.collaborators.add(r2.uid)
                r2.collaborators.add(r1.uid)
    
    return G


def analyze_network(G: nx.Graph) -> Dict:
    """Compute network metrics."""
    metrics = {
        "n_nodes": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
        "density": nx.density(G),
        "avg_clustering": nx.average_clustering(G),
        "n_components": nx.number_connected_components(G),
    }
    
    degree_centrality = nx.degree_centrality(G)
    betweenness = nx.betweenness_centrality(G)
    
    metrics["avg_degree_centrality"] = np.mean(list(degree_centrality.values()))
    metrics["avg_betweenness"] = np.mean(list(betweenness.values()))
    metrics["max_degree_centrality"] = max(degree_centrality.values())
    
    # Power-law degree distribution check
    degrees = [d for _, d in G.degree()]
    metrics["degree_mean"] = np.mean(degrees)
    metrics["degree_std"] = np.std(degrees)
    metrics["degree_skewness"] = stats.skew(degrees) if len(degrees) > 2 else 0
    
    return metrics, degree_centrality, betweenness


# ── Research Impact Metrics ────────────────────────────────────────────────
def compute_impact_metrics(researchers: List[Researcher]) -> pd.DataFrame:
    """Compute multiple impact metrics for all researchers."""
    records = []
    for r in researchers:
        # Traditional metrics
        h = r.h_index
        total_cit = r.citations
        pubs = r.publications
        
        # Alternative metrics
        productivity_index = pubs / max(1, len(r.publication_history))
        citations_per_pub = total_cit / max(1, pubs)
        funding_efficiency = total_cit / max(1, r.cumulative_funding) if r.cumulative_funding > 0 else 0
        
        # Composite score (normalized)
        composite = 0.3 * min(h / 20, 1) + 0.3 * min(citations_per_pub / 10, 1) + \
                    0.2 * min(productivity_index / 5, 1) + 0.2 * min(funding_efficiency / 2, 1)
        
        records.append({
            "uid": r.uid,
            "gender": r.gender,
            "region": r.region,
            "field": r.field,
            "career_stage": r.career_stage,
            "h_index": h,
            "total_citations": total_cit,
            "publications": pubs,
            "cumulative_funding": r.cumulative_funding,
            "funded_rounds": r.funded_rounds,
            "productivity_index": productivity_index,
            "citations_per_pub": citations_per_pub,
            "funding_efficiency": funding_efficiency,
            "composite_score": composite,
            "intrinsic_quality": r.intrinsic_quality,
        })
    return pd.DataFrame(records)


# ── Gini Coefficient ───────────────────────────────────────────────────────
def gini_coefficient(values):
    """Compute Gini coefficient for inequality measurement."""
    values = np.array(values, dtype=float)
    values = values[values > 0]
    if len(values) == 0:
        return 0.0
    sorted_vals = np.sort(values)
    n = len(sorted_vals)
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * sorted_vals) / (n * np.sum(sorted_vals))) - (n + 1) / n


# ── KAKENHI Case Study ─────────────────────────────────────────────────────
def kakenhi_simulation(config: SimConfig) -> Dict:
    """Simulate KAKENHI-style multi-category funding."""
    categories = {
        "S": {"budget": 30, "grant_size": 15, "threshold": 0.85},
        "A": {"budget": 25, "grant_size": 8, "threshold": 0.7},
        "B": {"budget": 20, "grant_size": 5, "threshold": 0.5},
        "C": {"budget": 15, "grant_size": 3, "threshold": 0.3},
        "Early_Career": {"budget": 7, "grant_size": 2, "threshold": 0.0},
        "Challenging": {"budget": 3, "grant_size": 3, "threshold": 0.0},
    }
    
    researchers = [Researcher(i, config) for i in range(config.n_researchers)]
    results = {}
    
    for cat_name, cat_params in categories.items():
        if cat_name == "Early_Career":
            eligible = [r for r in researchers if r.career_stage == "early"]
        elif cat_name == "Challenging":
            eligible = researchers  # anyone can apply
        elif cat_name == "S":
            eligible = [r for r in researchers if r.career_stage == "senior"]
        else:
            eligible = researchers
        
        if not eligible:
            continue
            
        funded_ids = peer_review_allocation(
            eligible, cat_params["budget"], cat_params["grant_size"], config
        )
        
        n_funded = len(funded_ids)
        n_eligible = len(eligible)
        success_rate = n_funded / max(1, n_eligible)
        
        # Compute quality of funded
        funded_quality = [r.intrinsic_quality for r in eligible if r.uid in funded_ids]
        unfunded_quality = [r.intrinsic_quality for r in eligible if r.uid not in funded_ids]
        
        results[cat_name] = {
            "n_eligible": n_eligible,
            "n_funded": n_funded,
            "success_rate": success_rate,
            "avg_funded_quality": np.mean(funded_quality) if funded_quality else 0,
            "avg_unfunded_quality": np.mean(unfunded_quality) if unfunded_quality else 0,
            "quality_gap": (np.mean(funded_quality) - np.mean(unfunded_quality)) if funded_quality and unfunded_quality else 0,
        }
    
    return results


# ── Main Simulation ────────────────────────────────────────────────────────
def run_simulation(config: SimConfig = None):
    if config is None:
        config = SimConfig()
    
    researchers = [Researcher(i, config) for i in range(config.n_researchers)]
    
    # Build initial network
    G = build_collaboration_network(researchers, config)
    net_metrics, degree_cent, betweenness = analyze_network(G)
    
    mechanisms = {
        "Peer Review": peer_review_allocation,
        "Lottery": lottery_allocation,
        "Modified Lottery": modified_lottery_allocation,
        "Diversity-Constrained": diversity_constrained_allocation,
    }
    
    # Track results per mechanism
    mechanism_results = {name: {
        "total_output": [], "gini_funding": [], "gini_citations": [],
        "female_funding_rate": [], "field_entropy": [],
        "avg_quality_funded": [], "avg_quality_unfunded": [],
        "round_outputs": [],
    } for name in mechanisms}
    
    # Run separate simulation for each mechanism
    for mech_name, mech_func in mechanisms.items():
        # Fresh researchers for each mechanism
        mech_researchers = [Researcher(i, config) for i in range(config.n_researchers)]
        # Copy gender/region/field to match
        for i, r in enumerate(mech_researchers):
            r.gender = researchers[i].gender
            r.region = researchers[i].region
            r.field = researchers[i].field
            r.career_stage = researchers[i].career_stage
            r.intrinsic_quality = researchers[i].intrinsic_quality
            r.productivity = researchers[i].productivity
        
        for round_num in range(config.n_rounds):
            funded_ids = mech_func(mech_researchers, config.budget_per_round,
                                   config.grant_size, config)
            funded_set = set(funded_ids)
            
            round_output = 0
            for r in mech_researchers:
                r.total_proposals += 1
                is_funded = r.uid in funded_set
                if is_funded:
                    r.cumulative_funding += config.grant_size
                    r.funded_rounds += 1
                r.funding_history.append(is_funded)
                output, pubs, cits = r.produce_research(is_funded, config)
                round_output += output
            
            # Evolve network (new collaborations from co-funding)
            for i, r1 in enumerate(mech_researchers):
                if r1.uid in funded_set:
                    for r2 in mech_researchers:
                        if r2.uid != r1.uid and r2.uid in funded_set:
                            if r1.field == r2.field and np.random.random() < 0.1:
                                r1.collaborators.add(r2.uid)
                                r2.collaborators.add(r1.uid)
            
            # Compute round metrics
            fundings = [r.cumulative_funding for r in mech_researchers]
            citations = [r.citations for r in mech_researchers]
            
            mechanism_results[mech_name]["total_output"].append(round_output)
            mechanism_results[mech_name]["gini_funding"].append(gini_coefficient(fundings))
            mechanism_results[mech_name]["gini_citations"].append(gini_coefficient(citations))
            
            female_funded = sum(1 for r in mech_researchers if r.uid in funded_set and r.gender == "F")
            total_female = sum(1 for r in mech_researchers if r.gender == "F")
            mechanism_results[mech_name]["female_funding_rate"].append(
                female_funded / max(1, total_female)
            )
            
            # Field entropy
            field_counts = np.zeros(config.n_fields)
            for r in mech_researchers:
                if r.uid in funded_set:
                    field_counts[r.field] += 1
            field_probs = field_counts / max(1, field_counts.sum())
            field_probs = field_probs[field_probs > 0]
            mechanism_results[mech_name]["field_entropy"].append(
                -np.sum(field_probs * np.log2(field_probs)) if len(field_probs) > 0 else 0
            )
            
            funded_quality = [r.intrinsic_quality for r in mech_researchers if r.uid in funded_set]
            unfunded_quality = [r.intrinsic_quality for r in mech_researchers if r.uid not in funded_set]
            mechanism_results[mech_name]["avg_quality_funded"].append(
                np.mean(funded_quality) if funded_quality else 0
            )
            mechanism_results[mech_name]["avg_quality_unfunded"].append(
                np.mean(unfunded_quality) if unfunded_quality else 0
            )
        
        # Final impact metrics
        mechanism_results[mech_name]["final_metrics"] = compute_impact_metrics(mech_researchers)
    
    # KAKENHI simulation
    kakenhi_results = kakenhi_simulation(config)
    
    return mechanism_results, net_metrics, degree_cent, betweenness, G, kakenhi_results


# ── Visualization ──────────────────────────────────────────────────────────
def plot_all(mechanism_results, net_metrics, degree_cent, betweenness, G, kakenhi_results, config=None):
    if config is None:
        config = SimConfig()
    
    sns.set_style("whitegrid")
    colors = {"Peer Review": "#2196F3", "Lottery": "#FF9800",
              "Modified Lottery": "#4CAF50", "Diversity-Constrained": "#9C27B0"}
    
    # ── Fig 1: Network Structure ───────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Degree distribution
    degrees = [d for _, d in G.degree()]
    axes[0].hist(degrees, bins=30, color="#2196F3", alpha=0.7, edgecolor="white")
    axes[0].set_xlabel("Degree", fontsize=12)
    axes[0].set_ylabel("Frequency", fontsize=12)
    axes[0].set_title("Co-authorship Network Degree Distribution", fontsize=13)
    axes[0].axvline(np.mean(degrees), color="red", linestyle="--",
                    label=f"Mean={np.mean(degrees):.1f}")
    axes[0].legend()
    
    # Network visualization (subsample)
    sub_nodes = list(G.nodes())[:80]
    H = G.subgraph(sub_nodes)
    pos = nx.spring_layout(H, seed=42, k=0.3)
    node_colors = ["#E91E63" if G.nodes[n].get("gender") == "F" else "#2196F3"
                   for n in H.nodes()]
    node_sizes = [50 + 200 * degree_cent.get(n, 0) for n in H.nodes()]
    nx.draw_networkx(H, pos, ax=axes[1], node_color=node_colors,
                     node_size=node_sizes, with_labels=False, alpha=0.7,
                     edge_color="#CCCCCC", width=0.5)
    axes[1].set_title("Co-authorship Network (subset, n=80)", fontsize=13)
    
    # Centrality vs betweenness
    dc_vals = [degree_cent[n] for n in G.nodes()]
    bw_vals = [betweenness[n] for n in G.nodes()]
    gender_colors = ["#E91E63" if G.nodes[n].get("gender") == "F" else "#2196F3"
                     for n in G.nodes()]
    axes[2].scatter(dc_vals, bw_vals, c=gender_colors, alpha=0.5, s=20)
    axes[2].set_xlabel("Degree Centrality", fontsize=12)
    axes[2].set_ylabel("Betweenness Centrality", fontsize=12)
    axes[2].set_title("Centrality Measures by Gender", fontsize=13)
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0],[0],marker='o',color='w',markerfacecolor='#E91E63',markersize=8,label='Female'),
                       Line2D([0],[0],marker='o',color='w',markerfacecolor='#2196F3',markersize=8,label='Male')]
    axes[2].legend(handles=legend_elements)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "network_analysis.png", dpi=150, bbox_inches="tight")
    plt.close()
    
    # ── Fig 2: Efficiency Comparison ───────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    rounds = range(1, config.n_rounds + 1)
    
    for name in mechanism_results:
        axes[0, 0].plot(rounds, mechanism_results[name]["total_output"],
                        label=name, color=colors[name], linewidth=2)
    axes[0, 0].set_xlabel("Round")
    axes[0, 0].set_ylabel("Total Research Output")
    axes[0, 0].set_title("Research Output per Round")
    axes[0, 0].legend(fontsize=9)
    
    for name in mechanism_results:
        axes[0, 1].plot(rounds, mechanism_results[name]["gini_funding"],
                        label=name, color=colors[name], linewidth=2)
    axes[0, 1].set_xlabel("Round")
    axes[0, 1].set_ylabel("Gini Coefficient")
    axes[0, 1].set_title("Funding Inequality (Gini) Over Time")
    axes[0, 1].legend(fontsize=9)
    
    for name in mechanism_results:
        axes[1, 0].plot(rounds, mechanism_results[name]["avg_quality_funded"],
                        label=f"{name} (funded)", color=colors[name], linewidth=2)
        axes[1, 0].plot(rounds, mechanism_results[name]["avg_quality_unfunded"],
                        color=colors[name], linewidth=1, linestyle="--", alpha=0.5)
    axes[1, 0].set_xlabel("Round")
    axes[1, 0].set_ylabel("Avg Intrinsic Quality")
    axes[1, 0].set_title("Quality of Funded vs Unfunded (dashed)")
    axes[1, 0].legend(fontsize=8)
    
    for name in mechanism_results:
        axes[1, 1].plot(rounds, mechanism_results[name]["field_entropy"],
                        label=name, color=colors[name], linewidth=2)
    axes[1, 1].set_xlabel("Round")
    axes[1, 1].set_ylabel("Shannon Entropy")
    axes[1, 1].set_title("Field Diversity (Entropy) of Funded Projects")
    axes[1, 1].legend(fontsize=9)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "efficiency_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()
    
    # ── Fig 3: Fairness & Diversity ────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    for name in mechanism_results:
        axes[0, 0].plot(rounds, mechanism_results[name]["female_funding_rate"],
                        label=name, color=colors[name], linewidth=2)
    axes[0, 0].axhline(config.gender_ratio_female, color="gray", linestyle=":",
                        label=f"Population ratio ({config.gender_ratio_female})")
    axes[0, 0].set_xlabel("Round")
    axes[0, 0].set_ylabel("Female Funding Rate")
    axes[0, 0].set_title("Gender Equity in Funding")
    axes[0, 0].legend(fontsize=8)
    
    # Final citation Gini
    for name in mechanism_results:
        axes[0, 1].plot(rounds, mechanism_results[name]["gini_citations"],
                        label=name, color=colors[name], linewidth=2)
    axes[0, 1].set_xlabel("Round")
    axes[0, 1].set_ylabel("Gini Coefficient")
    axes[0, 1].set_title("Citation Inequality Over Time")
    axes[0, 1].legend(fontsize=9)
    
    # Box plot: Funding distribution by gender (final)
    box_data = []
    for name in mechanism_results:
        df = mechanism_results[name]["final_metrics"]
        for _, row in df.iterrows():
            box_data.append({"Mechanism": name, "Gender": row["gender"],
                             "Cumulative Funding": row["cumulative_funding"]})
    box_df = pd.DataFrame(box_data)
    sns.boxplot(data=box_df, x="Mechanism", y="Cumulative Funding", hue="Gender",
                ax=axes[1, 0], palette={"F": "#E91E63", "M": "#2196F3"})
    axes[1, 0].set_title("Funding Distribution by Gender (Final)")
    axes[1, 0].tick_params(axis='x', rotation=15)
    
    # Region balance
    region_data = []
    for name in mechanism_results:
        df = mechanism_results[name]["final_metrics"]
        for region in range(config.n_regions):
            region_funding = df[df["region"] == region]["cumulative_funding"].mean()
            region_data.append({"Mechanism": name, "Region": f"R{region}",
                                "Avg Funding": region_funding})
    region_df = pd.DataFrame(region_data)
    sns.barplot(data=region_df, x="Region", y="Avg Funding", hue="Mechanism",
                ax=axes[1, 1], palette=colors)
    axes[1, 1].set_title("Regional Funding Balance")
    axes[1, 1].legend(fontsize=8)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fairness_diversity.png", dpi=150, bbox_inches="tight")
    plt.close()
    
    # ── Fig 4: Career Paths ────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for name in mechanism_results:
        df = mechanism_results[name]["final_metrics"]
        axes[0].scatter(df["cumulative_funding"], df["total_citations"],
                        alpha=0.3, s=15, label=name, color=colors[name])
    axes[0].set_xlabel("Cumulative Funding")
    axes[0].set_ylabel("Total Citations")
    axes[0].set_title("Funding vs Citations")
    axes[0].legend(fontsize=9)
    
    # Career stage analysis
    career_data = []
    for name in mechanism_results:
        df = mechanism_results[name]["final_metrics"]
        for stage in ["early", "mid", "senior"]:
            stage_df = df[df["career_stage"] == stage]
            career_data.append({
                "Mechanism": name, "Career Stage": stage,
                "Avg Funding": stage_df["cumulative_funding"].mean(),
                "Avg Citations": stage_df["total_citations"].mean(),
            })
    career_df = pd.DataFrame(career_data)
    sns.barplot(data=career_df, x="Career Stage", y="Avg Funding", hue="Mechanism",
                ax=axes[1], palette=colors)
    axes[1].set_title("Funding by Career Stage")
    axes[1].legend(fontsize=8)
    
    # Composite score distribution
    for name in mechanism_results:
        df = mechanism_results[name]["final_metrics"]
        axes[2].hist(df["composite_score"], bins=20, alpha=0.4, label=name,
                     color=colors[name])
    axes[2].set_xlabel("Composite Impact Score")
    axes[2].set_ylabel("Count")
    axes[2].set_title("Distribution of Composite Impact Scores")
    axes[2].legend(fontsize=9)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "career_paths.png", dpi=150, bbox_inches="tight")
    plt.close()
    
    # ── Fig 5: Impact Metrics Comparison ───────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    metric_names = ["h_index", "citations_per_pub", "funding_efficiency"]
    metric_labels = ["h-index", "Citations per Publication", "Funding Efficiency (cit/$)"]
    
    for idx, (metric, label) in enumerate(zip(metric_names, metric_labels)):
        data_list = []
        for name in mechanism_results:
            df = mechanism_results[name]["final_metrics"]
            for val in df[metric]:
                data_list.append({"Mechanism": name, label: val})
        plot_df = pd.DataFrame(data_list)
        sns.violinplot(data=plot_df, x="Mechanism", y=label, ax=axes[idx],
                       palette=colors, inner="box", cut=0)
        axes[idx].set_title(f"{label} by Mechanism")
        axes[idx].tick_params(axis='x', rotation=15)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "impact_metrics.png", dpi=150, bbox_inches="tight")
    plt.close()
    
    # ── Fig 6: KAKENHI Case Study ──────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    cats = list(kakenhi_results.keys())
    success_rates = [kakenhi_results[c]["success_rate"] for c in cats]
    axes[0].bar(cats, success_rates, color=sns.color_palette("viridis", len(cats)))
    axes[0].set_ylabel("Success Rate")
    axes[0].set_title("KAKENHI Success Rate by Category")
    axes[0].tick_params(axis='x', rotation=30)
    
    funded_q = [kakenhi_results[c]["avg_funded_quality"] for c in cats]
    unfunded_q = [kakenhi_results[c]["avg_unfunded_quality"] for c in cats]
    x = np.arange(len(cats))
    axes[1].bar(x - 0.2, funded_q, 0.4, label="Funded", color="#4CAF50")
    axes[1].bar(x + 0.2, unfunded_q, 0.4, label="Unfunded", color="#FF5722")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(cats, rotation=30)
    axes[1].set_ylabel("Avg Intrinsic Quality")
    axes[1].set_title("Quality: Funded vs Unfunded")
    axes[1].legend()
    
    n_eligible = [kakenhi_results[c]["n_eligible"] for c in cats]
    n_funded = [kakenhi_results[c]["n_funded"] for c in cats]
    axes[2].bar(x - 0.2, n_eligible, 0.4, label="Eligible", color="#2196F3")
    axes[2].bar(x + 0.2, n_funded, 0.4, label="Funded", color="#FF9800")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(cats, rotation=30)
    axes[2].set_ylabel("Count")
    axes[2].set_title("Eligible vs Funded Researchers")
    axes[2].legend()
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "kakenhi_case_study.png", dpi=150, bbox_inches="tight")
    plt.close()
    
    # ── Fig 7: Pareto Front (Efficiency vs Fairness) ──────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))
    
    for name in mechanism_results:
        final_output = np.mean(mechanism_results[name]["total_output"][-5:])
        final_gini = np.mean(mechanism_results[name]["gini_funding"][-5:])
        ax.scatter(final_gini, final_output, s=200, color=colors[name],
                   label=name, zorder=5, edgecolors="black", linewidth=1.5)
    
    ax.set_xlabel("Funding Inequality (Gini) — lower is fairer", fontsize=12)
    ax.set_ylabel("Research Output — higher is more efficient", fontsize=12)
    ax.set_title("Efficiency–Fairness Trade-off (Pareto Front)", fontsize=14)
    ax.legend(fontsize=10)
    ax.annotate("Ideal →", xy=(0.1, 0.9), xycoords="axes fraction",
                fontsize=11, color="green", fontweight="bold")
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "pareto_front.png", dpi=150, bbox_inches="tight")
    plt.close()

    print("All figures saved to figures/")


# ── Summary Statistics ─────────────────────────────────────────────────────
def print_summary(mechanism_results, net_metrics, kakenhi_results):
    print("\n" + "="*70)
    print("SIMULATION RESULTS SUMMARY")
    print("="*70)
    
    print(f"\n── Network Analysis ──")
    for k, v in net_metrics.items():
        print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")
    
    print(f"\n── Mechanism Comparison (Final 5 Rounds Average) ──")
    summary_data = []
    for name in mechanism_results:
        r = mechanism_results[name]
        df = r["final_metrics"]
        row = {
            "Mechanism": name,
            "Avg Output": np.mean(r["total_output"][-5:]),
            "Funding Gini": np.mean(r["gini_funding"][-5:]),
            "Citation Gini": np.mean(r["gini_citations"][-5:]),
            "Female Rate": np.mean(r["female_funding_rate"][-5:]),
            "Field Entropy": np.mean(r["field_entropy"][-5:]),
            "Avg h-index": df["h_index"].mean(),
            "Avg Composite": df["composite_score"].mean(),
            "Funding Efficiency": df["funding_efficiency"].mean(),
        }
        summary_data.append(row)
    
    summary_df = pd.DataFrame(summary_data)
    print(summary_df.to_string(index=False))
    
    print(f"\n── KAKENHI Case Study ──")
    for cat, res in kakenhi_results.items():
        print(f"  {cat}: success_rate={res['success_rate']:.3f}, "
              f"quality_gap={res['quality_gap']:.3f}, "
              f"n_funded={res['n_funded']}/{res['n_eligible']}")
    
    return summary_df


# ── Main ───────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    config = SimConfig()
    print("Running simulation...")
    mechanism_results, net_metrics, degree_cent, betweenness, G, kakenhi_results = run_simulation(config)
    
    print("Generating plots...")
    plot_all(mechanism_results, net_metrics, degree_cent, betweenness, G, kakenhi_results, config)
    
    summary_df = print_summary(mechanism_results, net_metrics, kakenhi_results)
    
    # Save summary data
    summary_df.to_csv("summary_results.csv", index=False)
    
    # Save detailed metrics per mechanism
    for name in mechanism_results:
        df = mechanism_results[name]["final_metrics"]
        df.to_csv(f"metrics_{name.lower().replace(' ', '_').replace('-', '_')}.csv", index=False)
    
    print("\nSimulation complete. Files saved.")
