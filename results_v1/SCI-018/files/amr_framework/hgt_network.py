from __future__ import annotations

from collections import Counter
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from . import FIGURES_DIR, RESULTS_DIR, ROOT, log_event, save_json, seed_everything


class HGTNetwork:
    def __init__(self, seed: int = 42) -> None:
        seed_everything(seed)
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.graph = nx.DiGraph()
        self.plasmid_host_range = {
            "IncF": 0.45,
            "IncI": 0.55,
            "IncN": 0.60,
            "IncP": 0.78,
            "IncQ": 0.85,
        }
        self.last_transfer_history: pd.DataFrame | None = None
        self.last_infected_nodes: set[int] = set()

    def build_network(self, n_strains: int = 18, connectivity: float = 0.22) -> nx.DiGraph:
        species_pool = ["Escherichia", "Klebsiella", "Enterobacter", "Pseudomonas", "Acinetobacter"]
        self.graph.clear()
        for node in range(n_strains):
            species = species_pool[node % len(species_pool)]
            self.graph.add_node(node, label=f"{species}_{node}", species=species)
        mechanisms = ["conjugation", "transformation", "transduction"]
        for source in self.graph.nodes:
            for target in self.graph.nodes:
                if source == target:
                    continue
                if self.rng.random() < connectivity:
                    self.graph.add_edge(
                        source,
                        target,
                        weight=round(float(self.rng.uniform(0.15, 0.95)), 3),
                        mechanism=str(self.rng.choice(mechanisms, p=[0.55, 0.20, 0.25])),
                    )
        if self.graph.number_of_edges() == 0:
            self.graph.add_edge(0, 1, weight=0.5, mechanism="conjugation")
        return self.graph

    def _compatibility(self, source: int, target: int) -> float:
        source_species = self.graph.nodes[source]["species"]
        target_species = self.graph.nodes[target]["species"]
        if source_species == target_species:
            return 1.0
        if {source_species, target_species} <= {"Escherichia", "Klebsiella", "Enterobacter"}:
            return 0.85
        return 0.60

    def simulate_transfer(self, timesteps: int = 40, transfer_rate: float = 0.22, plasmid_type: str = "IncF") -> pd.DataFrame:
        if self.graph.number_of_nodes() == 0:
            self.build_network()
        host_range = self.plasmid_host_range[plasmid_type]
        infected = set(sorted(self.graph.nodes)[:2])
        first_generation = {node: 0 for node in infected}
        history = []
        for timestep in range(timesteps + 1):
            history.append({"time": timestep, "plasmid": plasmid_type, "infected": len(infected), "prevalence": len(infected) / self.graph.number_of_nodes()})
            newly_infected: set[int] = set()
            for source in infected:
                for target in self.graph.successors(source):
                    if target in infected or target in newly_infected:
                        continue
                    edge_weight = self.graph[source][target]["weight"]
                    probability = transfer_rate * edge_weight * host_range * self._compatibility(source, target)
                    if self.rng.random() < min(0.98, probability):
                        newly_infected.add(target)
                        if source in first_generation and first_generation[source] == 0:
                            first_generation[source] += 1
            retained = {node for node in infected if self.rng.random() > 0.025}
            infected = retained | newly_infected
            if not infected:
                infected = set(sorted(self.graph.nodes)[:1])
        self.last_transfer_history = pd.DataFrame(history)
        self.last_infected_nodes = infected
        self._last_R0 = float(np.mean(list(first_generation.values()))) if first_generation else 0.0
        return self.last_transfer_history.copy()

    def compute_network_metrics(self) -> dict[str, Any]:
        undirected = self.graph.to_undirected()
        betweenness = nx.betweenness_centrality(self.graph, weight="weight")
        degree_sequence = [degree for _node, degree in self.graph.degree()]
        return {
            "n_nodes": self.graph.number_of_nodes(),
            "n_edges": self.graph.number_of_edges(),
            "density": round(float(nx.density(self.graph)), 4),
            "average_degree": round(float(np.mean(degree_sequence)), 4),
            "clustering": round(float(nx.average_clustering(undirected)), 4),
            "betweenness": {str(node): round(value, 4) for node, value in betweenness.items()},
            "degree_distribution": degree_sequence,
        }

    def identify_superspreaders(self, top_n: int = 5) -> list[int]:
        betweenness = nx.betweenness_centrality(self.graph, weight="weight")
        ranked = sorted(betweenness.items(), key=lambda item: item[1], reverse=True)
        return [node for node, _score in ranked[:top_n]]

    def visualize_network(self) -> str:
        if self.graph.number_of_nodes() == 0:
            self.build_network()
        plt.figure(figsize=(10, 8))
        position = nx.spring_layout(self.graph, seed=self.seed, weight="weight")
        superspreaders = set(self.identify_superspreaders())
        node_colors = []
        for node in self.graph.nodes:
            if node in superspreaders:
                node_colors.append("tab:red")
            elif node in self.last_infected_nodes:
                node_colors.append("tab:orange")
            else:
                node_colors.append("tab:blue")
        edge_widths = [1.0 + 2.5 * self.graph[u][v]["weight"] for u, v in self.graph.edges]
        nx.draw_networkx(self.graph, pos=position, with_labels=False, node_color=node_colors, node_size=350, width=edge_widths, alpha=0.8, arrows=True)
        labels = {node: self.graph.nodes[node]["species"][0] + str(node) for node in self.graph.nodes}
        nx.draw_networkx_labels(self.graph, pos=position, labels=labels, font_size=7, font_color="black")
        plt.title("HGT network and superspreaders")
        plt.axis("off")
        output_path = FIGURES_DIR / "hgt_network.png"
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()
        log_event(
            phase="component_4",
            event_type="file_written",
            skill_or_tool="HGTNetwork",
            files_written=[str(output_path.relative_to(ROOT))],
        )
        return str(output_path)


class ARGSpreadSimulator:
    def __init__(self, network: HGTNetwork) -> None:
        self.network = network

    def simulate_arg_spread(self, n_timesteps: int = 40) -> pd.DataFrame:
        histories = []
        for plasmid_type, rate in {"IncF": 0.24, "IncN": 0.20, "IncP": 0.17}.items():
            history = self.network.simulate_transfer(timesteps=n_timesteps, transfer_rate=rate, plasmid_type=plasmid_type)
            histories.append(history)
        return pd.concat(histories, ignore_index=True)

    def compute_R0_plasmid(self) -> float:
        return round(float(getattr(self.network, "_last_R0", 0.0)), 4)

    def visualize_spread(self, spread_df: pd.DataFrame) -> str:
        plt.figure(figsize=(10, 5))
        for plasmid, subset in spread_df.groupby("plasmid"):
            plt.plot(subset["time"], subset["prevalence"], marker="o", label=plasmid)
        plt.xlabel("Time")
        plt.ylabel("ARG prevalence")
        plt.title("Plasmid-mediated ARG spread")
        plt.legend()
        plt.grid(alpha=0.3)
        output_path = FIGURES_DIR / "arg_spread_dynamics.png"
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()
        log_event(
            phase="component_4",
            event_type="file_written",
            skill_or_tool="ARGSpreadSimulator",
            files_written=[str(output_path.relative_to(ROOT))],
        )
        return str(output_path)


def run_component(seed: int = 42) -> dict[str, Any]:
    network = HGTNetwork(seed=seed)
    network.build_network(n_strains=20, connectivity=0.24)
    simulator = ARGSpreadSimulator(network)
    spread_df = simulator.simulate_arg_spread(n_timesteps=45)
    spread_path = RESULTS_DIR / "hgt_spread_timeseries.csv"
    spread_df.to_csv(spread_path, index=False)
    metrics = network.compute_network_metrics()
    superspreaders = network.identify_superspreaders()
    network.visualize_network()
    simulator.visualize_spread(spread_df)
    summary = {
        "metrics": metrics,
        "superspreaders": superspreaders,
        "R0_plasmid": simulator.compute_R0_plasmid(),
        "final_prevalence": spread_df.groupby("plasmid")["prevalence"].last().round(4).to_dict(),
    }
    save_json(RESULTS_DIR / "hgt_network_metrics.json", summary)
    log_event(
        phase="component_4",
        event_type="handoff_completed",
        skill_or_tool="hgt_network",
        handoff_out={"superspreaders": superspreaders, "R0_plasmid": summary["R0_plasmid"]},
        files_written=[str(spread_path.relative_to(ROOT))],
    )
    return summary
