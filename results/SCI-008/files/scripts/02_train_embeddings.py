#!/usr/bin/env python3
"""
Step 2: Graph Embedding Training & Comparison (TransE, RotatE, ComplEx)
Uses PyKEEN for training and evaluation.
"""

import json
import os
import time
from datetime import datetime

import numpy as np
import pandas as pd
import torch
from pykeen.pipeline import pipeline
from pykeen.triples import TriplesFactory

np.random.seed(42)
torch.manual_seed(42)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
LOG_DIR = os.path.join(BASE_DIR, "logs")
os.makedirs(RESULTS_DIR, exist_ok=True)


def load_triples():
    df = pd.read_csv(os.path.join(DATA_DIR, "triples.tsv"), sep="\t")
    triples = df[["head", "relation", "tail"]].values
    tf = TriplesFactory.from_labeled_triples(triples)
    training, testing, validation = tf.split([0.8, 0.1, 0.1], random_state=42)
    return training, testing, validation, tf


def train_and_evaluate(model_name, training, testing, validation, embedding_dim=128, epochs=150):
    print(f"\n--- Training {model_name} (dim={embedding_dim}, epochs={epochs}) ---")
    start = time.time()

    result = pipeline(
        training=training,
        testing=testing,
        validation=validation,
        model=model_name,
        model_kwargs=dict(embedding_dim=embedding_dim),
        training_kwargs=dict(num_epochs=epochs, batch_size=64),
        optimizer="Adam",
        optimizer_kwargs=dict(lr=0.001),
        negative_sampler="basic",
        negative_sampler_kwargs=dict(num_negs_per_pos=10),
        evaluator_kwargs=dict(filtered=True),
        random_seed=42,
        device="cpu",
    )

    elapsed = time.time() - start
    metrics = result.metric_results.to_dict()

    # Extract key metrics
    summary = {
        "model": model_name,
        "embedding_dim": embedding_dim,
        "epochs": epochs,
        "training_time_sec": round(elapsed, 2),
        "hits_at_1": metrics.get("hits_at_1", {}).get("both", {}).get("realistic", None),
        "hits_at_3": metrics.get("hits_at_3", {}).get("both", {}).get("realistic", None),
        "hits_at_10": metrics.get("hits_at_10", {}).get("both", {}).get("realistic", None),
        "mean_rank": metrics.get("mean_rank", {}).get("both", {}).get("realistic", None),
        "mean_reciprocal_rank": metrics.get("mean_reciprocal_rank", {}).get("both", {}).get("realistic", None),
    }

    # Try alternate key structures
    if summary["hits_at_10"] is None:
        for key, val in metrics.items():
            if "hits_at_10" in key.lower() or "10" in key:
                if isinstance(val, dict):
                    for k2, v2 in val.items():
                        if isinstance(v2, dict):
                            for k3, v3 in v2.items():
                                if isinstance(v3, (int, float)):
                                    summary["hits_at_10"] = v3
                                    break
                        elif isinstance(v2, (int, float)):
                            summary["hits_at_10"] = v2
                            break
                elif isinstance(val, (int, float)):
                    summary["hits_at_10"] = val

    print(f"  Hits@1:  {summary['hits_at_1']}")
    print(f"  Hits@3:  {summary['hits_at_3']}")
    print(f"  Hits@10: {summary['hits_at_10']}")
    print(f"  MRR:     {summary['mean_reciprocal_rank']}")
    print(f"  MR:      {summary['mean_rank']}")
    print(f"  Time:    {elapsed:.1f}s")

    return result, summary, metrics


def main():
    print("=== Graph Embedding Training & Comparison ===")
    training, testing, validation, tf = load_triples()

    print(f"Training triples: {training.num_triples}")
    print(f"Testing triples:  {testing.num_triples}")
    print(f"Validation triples: {validation.num_triples}")
    print(f"Entities: {tf.num_entities}, Relations: {tf.num_relations}")

    models = ["TransE", "RotatE", "ComplEx"]
    all_results = {}
    all_summaries = []
    all_metrics = {}

    for model_name in models:
        result, summary, metrics = train_and_evaluate(
            model_name, training, testing, validation,
            embedding_dim=128, epochs=150,
        )
        all_results[model_name] = result
        all_summaries.append(summary)
        all_metrics[model_name] = metrics

        # Save model
        model_dir = os.path.join(RESULTS_DIR, f"model_{model_name.lower()}")
        os.makedirs(model_dir, exist_ok=True)
        result.save_to_directory(model_dir)

    # Save comparison summary
    comparison_df = pd.DataFrame(all_summaries)
    comparison_df.to_csv(os.path.join(RESULTS_DIR, "model_comparison.csv"), index=False)

    # Save full metrics
    with open(os.path.join(RESULTS_DIR, "full_metrics.json"), "w") as f:
        json.dump(all_metrics, f, indent=2, default=str)

    # Save training/test split info
    split_info = {
        "training_triples": training.num_triples,
        "testing_triples": testing.num_triples,
        "validation_triples": validation.num_triples,
        "num_entities": tf.num_entities,
        "num_relations": tf.num_relations,
    }
    with open(os.path.join(RESULTS_DIR, "split_info.json"), "w") as f:
        json.dump(split_info, f, indent=2)

    # Log
    log_entry = {
        "timestamp": datetime.now().isoformat(),
        "phase": "embedding_training",
        "event_type": "model_comparison",
        "actor": "co-scientist",
        "skill_or_tool": "02_train_embeddings.py",
        "handoff_out": {"models": models, "summaries": all_summaries},
        "files_written": [
            "results/model_comparison.csv",
            "results/full_metrics.json",
            "results/split_info.json",
        ],
        "status": "ok",
    }
    with open(os.path.join(LOG_DIR, "process-log.jsonl"), "a") as f:
        f.write(json.dumps(log_entry, default=str) + "\n")

    print("\n=== Embedding training complete ===")
    print(comparison_df.to_string())


if __name__ == "__main__":
    main()
