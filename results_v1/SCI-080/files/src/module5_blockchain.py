from __future__ import annotations

import copy
import csv
import hashlib
import json
import math
import random
import statistics
import time
from dataclasses import dataclass, field
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

SEED = 42
DPI = 300
TOTAL_LOTS = 100
BATCH_SIZE = 25
BLOCK_DIFFICULTY = 2
SIMULATION_START = datetime(2024, 1, 1, 6, 0, 0)

BASE_DIR = Path(__file__).resolve().parents[1]
FIGURES_DIR = BASE_DIR / "figures"
RESULTS_DIR = BASE_DIR / "results"
DATA_DIR = BASE_DIR / "data"
LOGS_DIR = BASE_DIR / "logs"
REPORT_PATH = BASE_DIR / "report.md"
PROCESS_LOG_PATH = LOGS_DIR / "process-log.jsonl"
PREPROCESSING_LOG_PATH = DATA_DIR / "preprocessing-log.md"
CSV_PATH = DATA_DIR / "blockchain_transactions.csv"
METRICS_PATH = RESULTS_DIR / "module5_metrics.json"

STAGES = ["Farm", "Processing", "Distribution", "Retail", "Consumer"]
STAGE_ORDER = {stage: index for index, stage in enumerate(STAGES)}
STAGE_CONFIG = {
    "Farm": {
        "target_temp": 6.0,
        "spread": 1.4,
        "optimal_range": (2.0, 8.0),
        "hard_range": (0.0, 10.0),
        "locations": ["Farm North", "Farm East", "Farm South"],
        "certifications": ["GlobalG.A.P.", "Organic"],
        "deviation_rate": 0.10,
    },
    "Processing": {
        "target_temp": 4.0,
        "spread": 1.0,
        "optimal_range": (1.0, 5.0),
        "hard_range": (0.0, 8.0),
        "locations": ["Plant A", "Plant B"],
        "certifications": ["HACCP", "ISO 22000"],
        "deviation_rate": 0.12,
    },
    "Distribution": {
        "target_temp": 3.0,
        "spread": 0.9,
        "optimal_range": (1.0, 4.0),
        "hard_range": (0.0, 6.0),
        "locations": ["Hub West", "Hub Central", "Hub East"],
        "certifications": ["ColdChain", "GDP"],
        "deviation_rate": 0.10,
    },
    "Retail": {
        "target_temp": 4.0,
        "spread": 1.1,
        "optimal_range": (1.0, 5.0),
        "hard_range": (0.0, 7.0),
        "locations": ["Retail 01", "Retail 02", "Retail 03", "Retail 04"],
        "certifications": ["RetailQA", "FSMA"],
        "deviation_rate": 0.14,
    },
    "Consumer": {
        "target_temp": 5.0,
        "spread": 1.3,
        "optimal_range": (1.0, 6.0),
        "hard_range": (0.0, 8.0),
        "locations": ["Consumer Home North", "Consumer Home East", "Consumer Home South"],
        "certifications": ["ReceiptLogged"],
        "deviation_rate": 0.08,
    },
}
PRODUCT_TYPES = ["Leafy Greens", "Milk", "Chicken", "Berries", "Seafood"]
TRADITIONAL_TRACEBACK_HOURS = (6.0, 48.0)


@dataclass
class Block:
    index: int
    timestamp: str
    data: dict[str, Any]
    previous_hash: str
    nonce: int = 0
    hash: str = field(init=False)

    def __post_init__(self) -> None:
        self.hash = self.calculate_hash()

    def calculate_hash(self) -> str:
        payload = {
            "index": self.index,
            "timestamp": self.timestamp,
            "data": self.data,
            "previous_hash": self.previous_hash,
            "nonce": self.nonce,
        }
        serialized = json.dumps(payload, sort_keys=True, default=str)
        return hashlib.sha256(serialized.encode("utf-8")).hexdigest()

    def mine(self, difficulty: int) -> None:
        prefix = "0" * difficulty
        while not self.hash.startswith(prefix):
            self.nonce += 1
            self.hash = self.calculate_hash()


class MerkleTree:
    def __init__(self, transactions: list[dict[str, Any]]) -> None:
        self.transactions = [copy.deepcopy(tx) for tx in transactions]
        self.leaves = [self.hash_transaction(tx) for tx in self.transactions] or [self.hash_text("")]
        self.levels = self._build_levels(self.leaves)
        self.root = self.levels[-1][0]

    @staticmethod
    def hash_text(text: str) -> str:
        return hashlib.sha256(text.encode("utf-8")).hexdigest()

    @classmethod
    def hash_transaction(cls, transaction: dict[str, Any]) -> str:
        serialized = json.dumps(transaction, sort_keys=True, default=str)
        return cls.hash_text(serialized)

    def _build_levels(self, leaves: list[str]) -> list[list[str]]:
        levels = [leaves]
        current = leaves
        while len(current) > 1:
            if len(current) % 2 == 1:
                current = current + [current[-1]]
            next_level = []
            for i in range(0, len(current), 2):
                next_level.append(self.hash_text(current[i] + current[i + 1]))
            levels.append(next_level)
            current = next_level
        return levels

    def get_proof(self, index: int) -> list[dict[str, Any]]:
        proof: list[dict[str, Any]] = []
        current_index = index
        for level in self.levels[:-1]:
            padded_level = level if len(level) % 2 == 0 else level + [level[-1]]
            sibling_index = current_index ^ 1
            proof.append(
                {
                    "hash": padded_level[sibling_index],
                    "position": "left" if sibling_index < current_index else "right",
                }
            )
            current_index //= 2
        return proof

    @classmethod
    def verify_proof(
        cls, transaction: dict[str, Any], proof: list[dict[str, Any]], root: str
    ) -> bool:
        computed = cls.hash_transaction(transaction)
        for step in proof:
            if step["position"] == "left":
                computed = cls.hash_text(step["hash"] + computed)
            else:
                computed = cls.hash_text(computed + step["hash"])
        return computed == root


class Blockchain:
    def __init__(self, difficulty: int = 2) -> None:
        self.difficulty = difficulty
        self.chain = [self._create_genesis_block()]

    def _create_genesis_block(self) -> Block:
        data = {"transactions": [], "merkle_root": MerkleTree([]).root}
        block = Block(
            index=0,
            timestamp=SIMULATION_START.isoformat(),
            data=data,
            previous_hash="0",
        )
        block.mine(self.difficulty)
        return block

    def add_block(self, transactions: list[dict[str, Any]]) -> tuple[Block, MerkleTree]:
        copied_transactions = [copy.deepcopy(tx) for tx in transactions]
        merkle_tree = MerkleTree(copied_transactions)
        timestamp = (SIMULATION_START + timedelta(minutes=15 * len(self.chain))).isoformat()
        block = Block(
            index=len(self.chain),
            timestamp=timestamp,
            data={"transactions": copied_transactions, "merkle_root": merkle_tree.root},
            previous_hash=self.chain[-1].hash,
        )
        block.mine(self.difficulty)
        self.chain.append(block)
        return block, merkle_tree

    def is_valid(self) -> bool:
        prefix = "0" * self.difficulty
        for index in range(1, len(self.chain)):
            current = self.chain[index]
            previous = self.chain[index - 1]
            if current.previous_hash != previous.hash:
                return False
            if current.calculate_hash() != current.hash:
                return False
            if not current.hash.startswith(prefix):
                return False
            merkle_root = MerkleTree(current.data["transactions"]).root
            if merkle_root != current.data["merkle_root"]:
                return False
        return True

    def verify_transaction_inclusion(self, block_index: int, transaction_index: int) -> bool:
        block = self.chain[block_index]
        merkle_tree = MerkleTree(block.data["transactions"])
        proof = merkle_tree.get_proof(transaction_index)
        transaction = block.data["transactions"][transaction_index]
        return MerkleTree.verify_proof(transaction, proof, block.data["merkle_root"])

    def trace_lot(self, lot_id: str) -> list[dict[str, Any]]:
        matches: list[dict[str, Any]] = []
        for block in self.chain[1:]:
            for transaction in block.data["transactions"]:
                if transaction["lot_id"] == lot_id:
                    matches.append(transaction)
        return sorted(matches, key=lambda item: STAGE_ORDER[item["stage"]])


def ensure_directories() -> None:
    for directory in [FIGURES_DIR, RESULTS_DIR, DATA_DIR, LOGS_DIR]:
        directory.mkdir(parents=True, exist_ok=True)


def append_process_log(event_type: str, phase: str, handoff_in: dict[str, Any], handoff_out: dict[str, Any], files_written: list[str], status: str = "ok") -> None:
    record = {
        "timestamp": datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": "co-scientist-data-simulation",
        "handoff_in": handoff_in,
        "handoff_out": handoff_out,
        "files_written": files_written,
        "status": status,
    }
    with PROCESS_LOG_PATH.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record) + "\n")


def simulate_supply_chain() -> tuple[list[dict[str, Any]], dict[str, Any], Blockchain]:
    random.seed(SEED)
    np.random.seed(SEED)
    rng = np.random.default_rng(SEED)

    contaminated_numbers = sorted(rng.choice(np.arange(1, TOTAL_LOTS + 1), size=5, replace=False).tolist())
    contaminated_lots = {f"LOT-{number:03d}" for number in contaminated_numbers}
    contamination_stage_map = {
        lot_id: rng.choice(STAGES[:-1]).item() if hasattr(rng.choice(STAGES[:-1]), "item") else rng.choice(STAGES[:-1])
        for lot_id in contaminated_lots
    }

    transactions: list[dict[str, Any]] = []
    alert_count = 0
    hold_count = 0
    reject_count = 0
    contamination_edge_counts = {f"{src}->{dst}": 0 for src, dst in zip(STAGES[:-1], STAGES[1:])}

    for lot_number in range(1, TOTAL_LOTS + 1):
        lot_id = f"LOT-{lot_number:03d}"
        product_type = rng.choice(PRODUCT_TYPES).item() if hasattr(rng.choice(PRODUCT_TYPES), "item") else rng.choice(PRODUCT_TYPES)
        base_timestamp = SIMULATION_START + timedelta(minutes=25 * lot_number)
        contamination_stage = contamination_stage_map.get(lot_id)

        if contamination_stage is not None:
            for stage_index in range(STAGE_ORDER[contamination_stage], len(STAGES) - 1):
                edge_key = f"{STAGES[stage_index]}->{STAGES[stage_index + 1]}"
                contamination_edge_counts[edge_key] += 1

        for stage_index, stage in enumerate(STAGES):
            config = STAGE_CONFIG[stage]
            lower_opt, upper_opt = config["optimal_range"]
            lower_hard, upper_hard = config["hard_range"]
            temperature = float(rng.normal(config["target_temp"], config["spread"]))
            if float(rng.random()) < config["deviation_rate"]:
                direction = -1 if float(rng.random()) < 0.5 else 1
                temperature += direction * float(rng.uniform(2.5, 5.0))
            temperature = round(temperature, 2)
            timestamp = (base_timestamp + timedelta(hours=stage_index * 6) + timedelta(minutes=float(rng.uniform(0, 45)))).isoformat()
            location = rng.choice(config["locations"]).item() if hasattr(rng.choice(config["locations"]), "item") else rng.choice(config["locations"])
            certifications = list(config["certifications"])
            handler_id = f"{stage[:2].upper()}-{int(rng.integers(1000, 9999))}"
            alerts: list[str] = []
            action = "accept"
            is_contaminated = contamination_stage is not None and stage_index >= STAGE_ORDER[contamination_stage]
            recall_flag = is_contaminated and stage in {"Retail", "Consumer"}

            if temperature < lower_opt or temperature > upper_opt:
                alerts.append("Temperature deviation alert")
                action = "inspect"
            if temperature < lower_hard or temperature > upper_hard:
                alerts.append("Automatic hold threshold exceeded")
                action = "hold"
            if contamination_stage == stage:
                alerts.append("Contamination event registered")
            if is_contaminated:
                alerts.append("Contamination propagated")
                if action != "hold":
                    action = "hold"
            if recall_flag:
                alerts.append("Recall triggered")
                action = "reject"

            if alerts:
                alert_count += len(alerts)
            if action == "hold":
                hold_count += 1
            if action == "reject":
                reject_count += 1

            transactions.append(
                {
                    "lot_id": lot_id,
                    "product_type": product_type,
                    "stage": stage,
                    "temperature": temperature,
                    "timestamp": timestamp,
                    "location": location,
                    "handler_id": handler_id,
                    "certifications": certifications,
                    "alerts": alerts,
                    "action": action,
                    "is_contaminated": is_contaminated,
                    "contamination_stage": contamination_stage or "None",
                    "recall_flag": recall_flag,
                    "stage_index": stage_index,
                }
            )

    blockchain = Blockchain(difficulty=BLOCK_DIFFICULTY)
    merkle_checks: list[bool] = []
    linked_blocks = 0
    hashed_blocks = 0

    for start in range(0, len(transactions), BATCH_SIZE):
        batch = transactions[start : start + BATCH_SIZE]
        block_index = len(blockchain.chain)
        for transaction_index, transaction in enumerate(batch):
            transaction["block_index"] = block_index
            transaction["transaction_index"] = transaction_index
        block, _ = blockchain.add_block(batch)
        for transaction in batch:
            transaction["block_hash"] = block.hash
        merkle_checks.append(blockchain.verify_transaction_inclusion(block_index, 0))
        if blockchain.chain[block_index].previous_hash == blockchain.chain[block_index - 1].hash:
            linked_blocks += 1
        if blockchain.chain[block_index].calculate_hash() == blockchain.chain[block_index].hash:
            hashed_blocks += 1

    retail_contaminated = sorted(
        {
            transaction["lot_id"]
            for transaction in transactions
            if transaction["stage"] == "Retail" and transaction["is_contaminated"]
        }
    )
    consumer_exposed = sorted(
        {
            transaction["lot_id"]
            for transaction in transactions
            if transaction["stage"] == "Consumer" and transaction["is_contaminated"]
        }
    )

    trace_results: list[dict[str, Any]] = []
    traditional_seconds: list[float] = []
    blockchain_seconds: list[float] = []
    for lot_id in retail_contaminated:
        start_time = time.perf_counter()
        path = blockchain.trace_lot(lot_id)
        elapsed = time.perf_counter() - start_time
        blockchain_seconds.append(elapsed)
        traditional_time = float(rng.uniform(*TRADITIONAL_TRACEBACK_HOURS) * 3600.0)
        traditional_seconds.append(traditional_time)
        trace_results.append(
            {
                "lot_id": lot_id,
                "blockchain_seconds": elapsed,
                "traditional_seconds": traditional_time,
                "farm_origin": path[0]["location"] if path else "Unknown",
                "retail_location": next(item["location"] for item in path if item["stage"] == "Retail"),
                "contamination_stage": next(item["contamination_stage"] for item in path if item["contamination_stage"] != "None"),
                "path": [item["stage"] for item in path],
            }
        )

    speedup_samples = [traditional / blockchain for traditional, blockchain in zip(traditional_seconds, blockchain_seconds)]
    speedup_mean = float(np.mean(speedup_samples)) if speedup_samples else math.nan
    speedup_ci = np.percentile(speedup_samples, [2.5, 97.5]).tolist() if speedup_samples else [math.nan, math.nan]

    metrics = {
        "generated_at": datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "random_seed": SEED,
        "total_lots": TOTAL_LOTS,
        "contaminated_lots": len(contaminated_lots),
        "recalled_lots": len(retail_contaminated),
        "consumer_exposed_lots": len(consumer_exposed),
        "total_transactions": len(transactions),
        "transactions_per_lot": len(STAGES),
        "total_blocks": len(blockchain.chain),
        "chain_valid": blockchain.is_valid(),
        "merkle_verification_pass_rate": float(np.mean(merkle_checks)) if merkle_checks else 0.0,
        "linked_block_rate": linked_blocks / max(len(blockchain.chain) - 1, 1),
        "hashed_block_rate": hashed_blocks / max(len(blockchain.chain) - 1, 1),
        "alert_count": alert_count,
        "hold_count": hold_count,
        "reject_count": reject_count,
        "traceback_summary": {
            "blockchain_mean_seconds": float(np.mean(blockchain_seconds)) if blockchain_seconds else 0.0,
            "blockchain_max_seconds": float(np.max(blockchain_seconds)) if blockchain_seconds else 0.0,
            "blockchain_min_seconds": float(np.min(blockchain_seconds)) if blockchain_seconds else 0.0,
            "traditional_mean_seconds": float(np.mean(traditional_seconds)) if traditional_seconds else 0.0,
            "traditional_max_seconds": float(np.max(traditional_seconds)) if traditional_seconds else 0.0,
            "traditional_min_seconds": float(np.min(traditional_seconds)) if traditional_seconds else 0.0,
            "speedup_factor_mean": speedup_mean,
            "speedup_factor_95pct_interval": speedup_ci,
        },
        "recall_scope": {
            "retail_contaminated_lots": retail_contaminated,
            "consumer_exposed_lots": consumer_exposed,
            "affected_transactions": int(sum(tx["is_contaminated"] for tx in transactions)),
        },
        "trace_results": trace_results,
        "contamination_edge_counts": contamination_edge_counts,
    }
    return transactions, metrics, blockchain


def write_transactions_csv(transactions: list[dict[str, Any]]) -> None:
    dataframe = pd.DataFrame(transactions)
    dataframe["certifications"] = dataframe["certifications"].apply(lambda items: "; ".join(items))
    dataframe["alerts"] = dataframe["alerts"].apply(lambda items: "; ".join(items) if items else "None")
    dataframe.to_csv(CSV_PATH, index=False, quoting=csv.QUOTE_MINIMAL)


def save_metrics(metrics: dict[str, Any]) -> None:
    with METRICS_PATH.open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)


def plot_supply_chain_network(metrics: dict[str, Any]) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(10, 4.6), dpi=DPI)
    graph = nx.DiGraph()
    for stage in STAGES:
        graph.add_node(stage)
    for source, destination in zip(STAGES[:-1], STAGES[1:]):
        graph.add_edge(source, destination)

    positions = {stage: (index, 0) for index, stage in enumerate(STAGES)}
    node_colors = plt.cm.cividis(np.linspace(0.15, 0.85, len(STAGES)))
    nx.draw_networkx_nodes(graph, positions, node_size=2600, node_color=node_colors, ax=ax)
    nx.draw_networkx_labels(
        graph,
        positions,
        labels={stage: f"{stage}\n100 lots" for stage in STAGES},
        font_size=9,
        ax=ax,
    )
    nx.draw_networkx_edges(graph, positions, width=3.0, edge_color="#7f7f7f", arrows=True, ax=ax)

    contaminated_edges = []
    contaminated_widths = []
    edge_labels: dict[tuple[str, str], str] = {}
    for source, destination in zip(STAGES[:-1], STAGES[1:]):
        key = f"{source}->{destination}"
        contaminated_count = metrics["contamination_edge_counts"][key]
        edge_labels[(source, destination)] = f"Flow: 100\nContaminated: {contaminated_count}"
        if contaminated_count > 0:
            contaminated_edges.append((source, destination))
            contaminated_widths.append(2.0 + contaminated_count * 0.6)
    nx.draw_networkx_edges(
        graph,
        positions,
        edgelist=contaminated_edges,
        width=contaminated_widths,
        edge_color="#d55e00",
        arrows=True,
        ax=ax,
    )
    nx.draw_networkx_edge_labels(graph, positions, edge_labels=edge_labels, font_size=8, ax=ax)

    ax.set_title("Supply Chain Network and Contamination Paths")
    ax.text(
        1.5,
        -0.22,
        f"Recall scope: {metrics['recalled_lots']} retail lots, {metrics['consumer_exposed_lots']} consumer exposures",
        ha="center",
        va="center",
        fontsize=9,
    )
    ax.set_axis_off()
    fig.subplots_adjust(bottom=0.18, top=0.88)
    fig.savefig(FIGURES_DIR / "fig5_supply_chain_network.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def plot_traceback_time(metrics: dict[str, Any]) -> None:
    trace_results = metrics["trace_results"]
    blockchain_seconds = [item["blockchain_seconds"] for item in trace_results]
    traditional_seconds = [item["traditional_seconds"] for item in trace_results]
    means = [statistics.mean(blockchain_seconds), statistics.mean(traditional_seconds)]
    lowers = [
        np.percentile(blockchain_seconds, 2.5),
        np.percentile(traditional_seconds, 2.5),
    ]
    uppers = [
        np.percentile(blockchain_seconds, 97.5),
        np.percentile(traditional_seconds, 97.5),
    ]
    yerr = [
        [means[0] - lowers[0], means[1] - lowers[1]],
        [uppers[0] - means[0], uppers[1] - means[1]],
    ]

    fig, ax = plt.subplots(figsize=(6, 4), dpi=DPI)
    colors = [plt.cm.viridis(0.25), plt.cm.viridis(0.75)]
    ax.bar(["Blockchain", "Traditional"], means, yerr=yerr, color=colors, capsize=6)
    ax.set_yscale("log")
    ax.set_ylabel("Traceback time (seconds, log scale)")
    ax.set_title("Traceback Time Comparison")
    speedup = metrics["traceback_summary"]["speedup_factor_mean"]
    ax.text(0.5, max(means) * 1.1, f"Mean speedup: {speedup:,.0f}x", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig5b_traceback_time.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def plot_contamination_spread(transactions: list[dict[str, Any]], metrics: dict[str, Any]) -> None:
    contaminated_lots = metrics["recall_scope"]["retail_contaminated_lots"]
    fig, ax = plt.subplots(figsize=(10, 5), dpi=DPI)
    graph = nx.DiGraph()
    positions: dict[str, tuple[float, float]] = {}
    node_colors: list[Any] = []
    node_sizes: list[float] = []

    for lot_row, lot_id in enumerate(contaminated_lots):
        lot_records = sorted(
            [tx for tx in transactions if tx["lot_id"] == lot_id],
            key=lambda item: STAGE_ORDER[item["stage"]],
        )
        contamination_stage = next(record["contamination_stage"] for record in lot_records if record["contamination_stage"] != "None")
        for stage_col, record in enumerate(lot_records):
            node_id = f"{lot_id}-{record['stage']}"
            graph.add_node(node_id)
            positions[node_id] = (stage_col, -lot_row)
            if record["stage"] in {"Retail", "Consumer"} and record["recall_flag"]:
                node_colors.append("#e69f00")
                node_sizes.append(900)
            elif STAGE_ORDER[record["stage"]] >= STAGE_ORDER[contamination_stage]:
                node_colors.append("#d55e00")
                node_sizes.append(760)
            else:
                node_colors.append("#56b4e9")
                node_sizes.append(680)
        for left, right in zip(lot_records[:-1], lot_records[1:]):
            src = f"{lot_id}-{left['stage']}"
            dst = f"{lot_id}-{right['stage']}"
            graph.add_edge(src, dst)

    nx.draw_networkx_nodes(graph, positions, node_color=node_colors, node_size=node_sizes, ax=ax)
    edge_colors = []
    for src, dst in graph.edges():
        src_stage = src.split("-")[-1]
        dst_stage = dst.split("-")[-1]
        lot_id = src.rsplit("-", 1)[0]
        contamination_stage = next(
            item["contamination_stage"]
            for item in transactions
            if item["lot_id"] == lot_id and item["contamination_stage"] != "None"
        )
        if STAGE_ORDER[src_stage] >= STAGE_ORDER[contamination_stage] or STAGE_ORDER[dst_stage] >= STAGE_ORDER[contamination_stage]:
            edge_colors.append("#d55e00")
        else:
            edge_colors.append("#999999")
    nx.draw_networkx_edges(graph, positions, edge_color=edge_colors, width=2.2, arrows=True, ax=ax)
    nx.draw_networkx_labels(
        graph,
        positions,
        labels={node: node.replace("-", "\n", 1) for node in graph.nodes()},
        font_size=7,
        ax=ax,
    )
    ax.set_title("Contamination Propagation and Recall Scope")
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig5c_contamination_spread.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def plot_temperature_log(transactions: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(9, 5), dpi=DPI)
    cmap = plt.cm.viridis(np.linspace(0.15, 0.85, len(STAGES)))

    for index, stage in enumerate(STAGES):
        stage_temperatures = [tx["temperature"] for tx in transactions if tx["stage"] == stage]
        jitter = np.random.default_rng(SEED + index).normal(0, 0.04, len(stage_temperatures))
        x_values = np.full(len(stage_temperatures), index, dtype=float) + jitter
        lower_opt, upper_opt = STAGE_CONFIG[stage]["optimal_range"]
        lower_hard, upper_hard = STAGE_CONFIG[stage]["hard_range"]
        ax.fill_between([index - 0.32, index + 0.32], lower_hard, upper_hard, color="#d9d9d9", alpha=0.35)
        ax.fill_between([index - 0.24, index + 0.24], lower_opt, upper_opt, color="#009e73", alpha=0.18)
        ax.scatter(x_values, stage_temperatures, color=cmap[index], s=18, alpha=0.75)
        alert_points = [
            tx["temperature"]
            for tx in transactions
            if tx["stage"] == stage and any("Temperature deviation alert" == alert for alert in tx["alerts"])
        ]
        if alert_points:
            ax.scatter(
                np.full(len(alert_points), index),
                alert_points,
                color="#d55e00",
                s=28,
                marker="x",
                linewidths=1.0,
            )

    ax.set_xticks(range(len(STAGES)))
    ax.set_xticklabels(STAGES)
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("Temperature Monitoring Across the Supply Chain")
    ax.set_ylim(-2, 12)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig5d_temperature_log.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def write_preprocessing_log(transactions: list[dict[str, Any]], metrics: dict[str, Any]) -> None:
    content = f"""# Preprocessing Log

- Seed fixed to {SEED} for `random` and `numpy`.
- Simulated {metrics['total_lots']} lots across {len(STAGES)} supply-chain stages.
- Generated {len(transactions)} stage-level transaction records.
- Serialized certifications and alerts for CSV export.
- Batched transactions into {metrics['total_blocks'] - 1} mined blocks plus one genesis block.
- Verified blockchain link integrity and Merkle inclusion before saving outputs.
"""
    PREPROCESSING_LOG_PATH.write_text(content, encoding="utf-8")


def write_report(metrics: dict[str, Any]) -> None:
    traceback_summary = metrics["traceback_summary"]
    report = f"""# DRAFT — NOT FOR DISTRIBUTION

## Module 5: Blockchain-Based Supply Chain Traceability for Food Safety

**Timestamp:** {metrics['generated_at']}

## Methods

- Implemented a simplified blockchain with block fields `index`, `timestamp`, `data`, `previous_hash`, `hash`, and `nonce`.
- Used SHA-256 hashing, block mining, chain validation, and Merkle-tree verification.
- Simulated {metrics['total_lots']} product lots across Farm → Processing → Distribution → Retail → Consumer.
- Injected contamination events into 5% of lots and evaluated smart-contract style alert, hold, reject, and recall logic.
- Compared blockchain traceback times against simulated traditional traceback times.

## Results

- Total transactions: {metrics['total_transactions']}
- Total blocks: {metrics['total_blocks']}
- Chain valid: {metrics['chain_valid']}
- Merkle verification pass rate: {metrics['merkle_verification_pass_rate']:.2%}
- Alerts generated: {metrics['alert_count']}
- Holds generated: {metrics['hold_count']}
- Rejects generated: {metrics['reject_count']}
- Recalled lots at retail: {metrics['recalled_lots']}
- Mean blockchain traceback time: {traceback_summary['blockchain_mean_seconds']:.6f} seconds
- Mean traditional traceback time: {traceback_summary['traditional_mean_seconds'] / 3600:.2f} hours
- Mean speedup: {traceback_summary['speedup_factor_mean']:,.0f}x
- Speedup 95% interval: [{traceback_summary['speedup_factor_95pct_interval'][0]:,.0f}x, {traceback_summary['speedup_factor_95pct_interval'][1]:,.0f}x]

## Discussion

The simulation shows that tamper-evident lot records and indexed block traversal support near-instant traceback for contaminated retail products. Automated temperature-deviation alerts and recall triggers narrow the recall scope to affected lots while preserving chain-integrity checks. Traditional traceback remains orders of magnitude slower because it is modeled as manual, cross-organizational document retrieval.

## Limitations

- This is a synthetic simulation and not a production blockchain deployment.
- Traditional traceback time is modeled rather than observed from real-world operations.
- Temperature thresholds are generalized food-safety assumptions for demonstration.

## File Inventory

- `src/module5_blockchain.py`
- `data/blockchain_transactions.csv`
- `data/preprocessing-log.md`
- `results/module5_metrics.json`
- `figures/fig5_supply_chain_network.png`
- `figures/fig5b_traceback_time.png`
- `figures/fig5c_contamination_spread.png`
- `figures/fig5d_temperature_log.png`
- `logs/process-log.jsonl`
"""
    REPORT_PATH.write_text(report, encoding="utf-8")


def validate_outputs(metrics: dict[str, Any]) -> None:
    expected_files = [
        FIGURES_DIR / "fig5_supply_chain_network.png",
        FIGURES_DIR / "fig5b_traceback_time.png",
        FIGURES_DIR / "fig5c_contamination_spread.png",
        FIGURES_DIR / "fig5d_temperature_log.png",
        METRICS_PATH,
        CSV_PATH,
        REPORT_PATH,
        PREPROCESSING_LOG_PATH,
        PROCESS_LOG_PATH,
    ]
    missing = [str(path) for path in expected_files if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing expected outputs: {missing}")
    if not metrics["chain_valid"]:
        raise ValueError("Blockchain validation failed")


def main() -> None:
    ensure_directories()
    append_process_log(
        event_type="run_started",
        phase="PLAN",
        handoff_in={"seed": SEED},
        handoff_out={"objective": "Supply chain blockchain simulation"},
        files_written=[],
    )
    append_process_log(
        event_type="prompt_received",
        phase="PLAN",
        handoff_in={"lots": TOTAL_LOTS, "stages": STAGES},
        handoff_out={"required_outputs": [str(CSV_PATH), str(METRICS_PATH)]},
        files_written=[],
    )
    append_process_log(
        event_type="skill_selected",
        phase="PLAN",
        handoff_in={"skill": "co-scientist-data-simulation"},
        handoff_out={"method": "synthetic blockchain traceability simulation"},
        files_written=[],
    )
    append_process_log(
        event_type="handoff_started",
        phase="EXECUTE",
        handoff_in={"script": str(Path(__file__).name)},
        handoff_out={"status": "running"},
        files_written=[],
    )

    transactions, metrics, blockchain = simulate_supply_chain()
    write_transactions_csv(transactions)
    save_metrics(metrics)
    plot_supply_chain_network(metrics)
    plot_traceback_time(metrics)
    plot_contamination_spread(transactions, metrics)
    plot_temperature_log(transactions)
    write_preprocessing_log(transactions, metrics)
    write_report(metrics)

    append_process_log(
        event_type="handoff_completed",
        phase="EXECUTE",
        handoff_in={"transactions": len(transactions)},
        handoff_out={"blocks": len(blockchain.chain), "chain_valid": metrics["chain_valid"]},
        files_written=[str(CSV_PATH), str(METRICS_PATH)],
    )
    for file_path in [
        CSV_PATH,
        METRICS_PATH,
        FIGURES_DIR / "fig5_supply_chain_network.png",
        FIGURES_DIR / "fig5b_traceback_time.png",
        FIGURES_DIR / "fig5c_contamination_spread.png",
        FIGURES_DIR / "fig5d_temperature_log.png",
        PREPROCESSING_LOG_PATH,
    ]:
        append_process_log(
            event_type="file_written",
            phase="REPORT",
            handoff_in={"path": str(file_path)},
            handoff_out={"exists": file_path.exists()},
            files_written=[str(file_path)],
        )

    append_process_log(
        event_type="report_finalized",
        phase="REPORT",
        handoff_in={"report": str(REPORT_PATH)},
        handoff_out={"metrics_summary": metrics["traceback_summary"]},
        files_written=[str(REPORT_PATH)],
    )
    validate_outputs(metrics)
    append_process_log(
        event_type="run_completed",
        phase="LOG",
        handoff_in={"validation": "passed"},
        handoff_out={"status": "completed"},
        files_written=[str(PROCESS_LOG_PATH)],
    )
    print(json.dumps({
        "script": str(Path(__file__)),
        "chain_valid": metrics["chain_valid"],
        "recalled_lots": metrics["recalled_lots"],
        "mean_traceback_seconds": metrics["traceback_summary"]["blockchain_mean_seconds"],
        "speedup_factor_mean": metrics["traceback_summary"]["speedup_factor_mean"],
    }, indent=2))


if __name__ == "__main__":
    main()
