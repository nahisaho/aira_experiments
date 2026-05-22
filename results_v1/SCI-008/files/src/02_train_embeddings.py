"""
Graph Embedding Training: TransE, RotatE, ComplEx via PyKEEN
Compares models on link prediction (MRR, Hits@1, Hits@3, Hits@10)
"""

import json
import os
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch

warnings.filterwarnings("ignore")

BASE = Path(__file__).parent.parent
DATA_DIR = BASE / "data"
RESULTS_DIR = BASE / "results"
LOG_FILE = BASE / "logs" / "process-log.jsonl"
RESULTS_DIR.mkdir(exist_ok=True)

torch.manual_seed(42)
np.random.seed(42)


def log_event(event_type, details):
    entry = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "phase": "EXECUTE",
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": "co-scientist-drug-repurposing",
        **details,
    }
    with open(LOG_FILE, "a") as f:
        f.write(json.dumps(entry) + "\n")


log_event("handoff_started", {"script": "02_train_embeddings.py", "input": "data/kg_triples.tsv"})

# ─────────────────────────────────────────────
# Load triples
# ─────────────────────────────────────────────
df = pd.read_csv(DATA_DIR / "kg_triples.tsv", sep="\t")
print(f"Loaded {len(df)} triples")

from pykeen.triples import TriplesFactory
from pykeen.pipeline import pipeline
from pykeen.models import TransE, RotatE, ComplEx

# Build TriplesFactory
tf = TriplesFactory.from_labeled_triples(
    triples=df[["head", "relation", "tail"]].values,
    create_inverse_triples=True,
)

# Train/validation/test split: 80/10/10
training, testing, validation = tf.split([0.8, 0.1, 0.1], random_state=42)

print(f"Training triples: {training.num_triples}")
print(f"Validation triples: {validation.num_triples}")
print(f"Testing triples: {testing.num_triples}")

# ─────────────────────────────────────────────
# Train models
# ─────────────────────────────────────────────
MODELS_CONFIG = {
    "TransE": {
        "model": "TransE",
        "embedding_dim": 64,
        "epochs": 100,
    },
    "RotatE": {
        "model": "RotatE",
        "embedding_dim": 64,
        "epochs": 100,
    },
    "ComplEx": {
        "model": "ComplEx",
        "embedding_dim": 64,
        "epochs": 100,
    },
}

results_summary = {}

for model_name, config in MODELS_CONFIG.items():
    print(f"\n{'='*50}")
    print(f"Training {model_name}...")
    start = time.time()

    result = pipeline(
        training=training,
        testing=testing,
        validation=validation,
        model=config["model"],
        model_kwargs=dict(embedding_dim=config["embedding_dim"]),
        training_kwargs=dict(
            num_epochs=config["epochs"],
            batch_size=64,
        ),
        optimizer="Adam",
        optimizer_kwargs=dict(lr=0.01),
        negative_sampler="basic",
        negative_sampler_kwargs=dict(num_negs_per_pos=10),
        evaluator="RankBasedEvaluator",
        evaluator_kwargs=dict(filtered=True),
        random_seed=42,
        device="cpu",
    )

    elapsed = time.time() - start
    metrics = result.metric_results.to_dict()

    # Extract key metrics
    mrr = metrics.get("both.realistic.inverse_harmonic_mean_rank", 0)
    h1 = metrics.get("both.realistic.hits_at_1", 0)
    h3 = metrics.get("both.realistic.hits_at_3", 0)
    h10 = metrics.get("both.realistic.hits_at_10", 0)

    results_summary[model_name] = {
        "mrr": round(mrr, 4),
        "hits_at_1": round(h1, 4),
        "hits_at_3": round(h3, 4),
        "hits_at_10": round(h10, 4),
        "training_time_sec": round(elapsed, 1),
        "embedding_dim": config["embedding_dim"],
        "epochs": config["epochs"],
    }

    print(f"  MRR:      {mrr:.4f}")
    print(f"  Hits@1:   {h1:.4f}")
    print(f"  Hits@3:   {h3:.4f}")
    print(f"  Hits@10:  {h10:.4f}")
    print(f"  Time:     {elapsed:.1f}s")

    # Save model
    result.save_to_directory(RESULTS_DIR / f"model_{model_name.lower()}")

    log_event("model_trained", {
        "model": model_name,
        "metrics": results_summary[model_name],
        "status": "ok"
    })

# Save comparison table
df_results = pd.DataFrame(results_summary).T
df_results.index.name = "model"
df_results.to_csv(RESULTS_DIR / "embedding_comparison.csv")

print("\n=== Model Comparison ===")
print(df_results.to_string())

# Save best model info
best_model = df_results["mrr"].idxmax()
print(f"\nBest model by MRR: {best_model} (MRR={df_results.loc[best_model, 'mrr']:.4f})")

# Save factory for reuse
import pickle
with open(DATA_DIR / "triples_factory.pkl", "wb") as f:
    pickle.dump({"training": training, "testing": testing,
                 "validation": validation, "tf": tf}, f)

log_event("handoff_completed", {
    "files_written": [
        "results/embedding_comparison.csv",
        f"results/model_transe/",
        f"results/model_rotate/",
        f"results/model_complex/",
        "data/triples_factory.pkl",
    ],
    "best_model": best_model,
    "status": "ok"
})

print("\n[✓] Embedding training complete.")
