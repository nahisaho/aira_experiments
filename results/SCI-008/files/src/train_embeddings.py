"""
Train and evaluate KGE models (TransE, RotatE, ComplEx) using PyKEEN.
Performs link prediction and identifies drug repurposing candidates.
"""

import os
import json
import warnings
import numpy as np
import pandas as pd
import torch
from pykeen.pipeline import pipeline
from pykeen.triples import TriplesFactory
from pykeen.predict import predict_target

warnings.filterwarnings("ignore")

SEED = 42
torch.manual_seed(SEED)
np.random.seed(SEED)

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

MODELS = ["TransE", "RotatE", "ComplEx"]
EPOCHS = 200
EMBEDDING_DIM = 128
BATCH_SIZE = 64


def load_triples():
    """Load KG triples and create PyKEEN TriplesFactory."""
    df = pd.read_csv(os.path.join(DATA_DIR, "kg_triples.tsv"), sep="\t")
    triples = df[["head", "relation", "tail"]].values
    tf = TriplesFactory.from_labeled_triples(triples)
    training, testing, validation = tf.split([0.8, 0.1, 0.1], random_state=SEED)
    return tf, training, testing, validation


def train_model(model_name, training, testing, validation):
    """Train a single KGE model."""
    print(f"\n{'='*60}")
    print(f"Training {model_name}...")
    print(f"{'='*60}")

    model_kwargs = {"embedding_dim": EMBEDDING_DIM}

    result = pipeline(
        model=model_name,
        training=training,
        testing=testing,
        validation=validation,
        model_kwargs=model_kwargs,
        training_kwargs={
            "num_epochs": EPOCHS,
            "batch_size": BATCH_SIZE,
        },
        optimizer="Adam",
        optimizer_kwargs={"lr": 0.001},
        negative_sampler="basic",
        negative_sampler_kwargs={"num_negs_per_pos": 10},
        evaluator_kwargs={"filtered": True},
        random_seed=SEED,
        device="cpu",
    )

    return result


def extract_metrics(result):
    """Extract evaluation metrics from pipeline result."""
    metrics = result.metric_results.to_dict()

    # Extract key metrics
    extracted = {}
    for side in ["head", "tail", "both"]:
        prefix = side
        for metric_name in ["hits_at_1", "hits_at_3", "hits_at_5", "hits_at_10",
                            "mean_rank", "mean_reciprocal_rank",
                            "adjusted_mean_rank_index"]:
            key = f"{prefix}.{metric_name}"
            if "realistic" in metrics:
                val = metrics.get("realistic", {}).get(key)
            else:
                # try to find in flat structure
                for k, v in metrics.items():
                    if metric_name in k and side in k:
                        val = v
                        break
                else:
                    val = None
            if val is not None:
                extracted[f"{side}_{metric_name}"] = float(val) if val is not None else None

    # Flatten all available metrics
    flat = {}
    def flatten(d, prefix=""):
        if isinstance(d, dict):
            for k, v in d.items():
                flatten(v, f"{prefix}{k}.")
        else:
            flat[prefix.rstrip(".")] = d
    flatten(metrics)

    # Get key summary metrics
    summary = {}
    for k, v in flat.items():
        if any(m in k for m in ["hits_at_10", "hits_at_1", "mean_reciprocal_rank", "mean_rank"]):
            if "both" in k or "realistic" in k:
                clean_key = k.split(".")[-1] if "." in k else k
                summary[clean_key] = v

    return {**extracted, **summary, "all_metrics": flat}


def predict_drug_disease(result, tf, top_k=20):
    """Predict new drug-disease associations."""
    model = result.model
    model.eval()

    entity_to_id = tf.entity_to_id
    relation_to_id = tf.relation_to_id

    predictions = []

    # Load entity types
    with open(os.path.join(DATA_DIR, "entity_types.json")) as f:
        entity_types = json.load(f)

    drugs = [e for e, t in entity_types.items() if t == "Drug" and e in entity_to_id]
    diseases = [e for e, t in entity_types.items() if t == "Disease" and e in entity_to_id]

    treat_rel = "drug_treats_disease"
    if treat_rel not in relation_to_id:
        print("Warning: 'drug_treats_disease' relation not found")
        return pd.DataFrame()

    rel_id = relation_to_id[treat_rel]

    # Get existing drug-disease triples
    existing = set()
    triples_np = tf.mapped_triples.numpy()
    id_to_entity = {v: k for k, v in entity_to_id.items()}
    id_to_relation = {v: k for k, v in relation_to_id.items()}

    for row in triples_np:
        h, r, t = int(row[0]), int(row[1]), int(row[2])
        if r == rel_id:
            existing.add((id_to_entity[h], id_to_entity[t]))

    # Score all drug-disease pairs
    with torch.no_grad():
        for drug in drugs:
            for disease in diseases:
                if (drug, disease) in existing:
                    continue
                h_id = torch.tensor([entity_to_id[drug]], dtype=torch.long)
                r_id_t = torch.tensor([rel_id], dtype=torch.long)
                t_id = torch.tensor([entity_to_id[disease]], dtype=torch.long)
                triple = torch.stack([h_id, r_id_t, t_id], dim=1)
                score = model.score_hrt(triple).item()
                predictions.append({
                    "drug": drug,
                    "disease": disease,
                    "score": score,
                    "is_novel": True,
                })

    pred_df = pd.DataFrame(predictions)
    if len(pred_df) > 0:
        pred_df = pred_df.sort_values("score", ascending=False).head(top_k)
        pred_df["rank"] = range(1, len(pred_df) + 1)

    return pred_df


def predict_covid_candidates(result, tf, top_k=15):
    """Specifically predict COVID-19 drug candidates."""
    model = result.model
    model.eval()

    entity_to_id = tf.entity_to_id
    relation_to_id = tf.relation_to_id

    with open(os.path.join(DATA_DIR, "entity_types.json")) as f:
        entity_types = json.load(f)

    drugs = [e for e, t in entity_types.items() if t == "Drug" and e in entity_to_id]

    if "COVID-19" not in entity_to_id:
        return pd.DataFrame()

    treat_rel = "drug_treats_disease"
    if treat_rel not in relation_to_id:
        return pd.DataFrame()

    rel_id = relation_to_id[treat_rel]
    covid_id = entity_to_id["COVID-19"]

    # Known COVID treatments
    known = set()
    triples_np = tf.mapped_triples.numpy()
    id_to_entity = {v: k for k, v in entity_to_id.items()}
    for row in triples_np:
        h, r, t = int(row[0]), int(row[1]), int(row[2])
        if r == rel_id and t == covid_id:
            known.add(id_to_entity[h])

    predictions = []
    with torch.no_grad():
        for drug in drugs:
            h_id = torch.tensor([entity_to_id[drug]], dtype=torch.long)
            r_id_t = torch.tensor([rel_id], dtype=torch.long)
            t_id = torch.tensor([covid_id], dtype=torch.long)
            triple = torch.stack([h_id, r_id_t, t_id], dim=1)
            score = model.score_hrt(triple).item()
            predictions.append({
                "drug": drug,
                "score": score,
                "known_treatment": drug in known,
            })

    pred_df = pd.DataFrame(predictions)
    pred_df = pred_df.sort_values("score", ascending=False).head(top_k)
    pred_df["rank"] = range(1, len(pred_df) + 1)
    return pred_df


def main():
    print("Loading knowledge graph triples...")
    tf, training, testing, validation = load_triples()
    print(f"  Total triples: {tf.num_triples}")
    print(f"  Training: {training.num_triples}")
    print(f"  Testing: {testing.num_triples}")
    print(f"  Validation: {validation.num_triples}")
    print(f"  Entities: {tf.num_entities}")
    print(f"  Relations: {tf.num_relations}")

    all_results = {}
    all_metrics = {}
    all_covid_preds = {}
    all_drug_disease_preds = {}

    for model_name in MODELS:
        result = train_model(model_name, training, testing, validation)
        metrics = extract_metrics(result)
        all_results[model_name] = result
        all_metrics[model_name] = metrics

        # Predictions
        covid_preds = predict_covid_candidates(result, tf)
        drug_disease_preds = predict_drug_disease(result, tf)
        all_covid_preds[model_name] = covid_preds
        all_drug_disease_preds[model_name] = drug_disease_preds

        print(f"\n--- {model_name} COVID-19 Drug Candidates ---")
        if len(covid_preds) > 0:
            print(covid_preds[["rank", "drug", "score", "known_treatment"]].to_string(index=False))

    # Save metrics comparison
    metrics_summary = []
    for model_name, metrics in all_metrics.items():
        row = {"model": model_name}
        for k, v in metrics.items():
            if k != "all_metrics" and isinstance(v, (int, float)):
                row[k] = v
        metrics_summary.append(row)

    metrics_df = pd.DataFrame(metrics_summary)
    metrics_df.to_csv(os.path.join(RESULTS_DIR, "model_comparison.csv"), index=False)
    print(f"\nMetrics saved to {RESULTS_DIR}/model_comparison.csv")

    # Save COVID predictions
    for model_name, preds in all_covid_preds.items():
        if len(preds) > 0:
            preds.to_csv(os.path.join(RESULTS_DIR, f"covid_predictions_{model_name}.csv"), index=False)

    # Save drug-disease predictions
    for model_name, preds in all_drug_disease_preds.items():
        if len(preds) > 0:
            preds.to_csv(os.path.join(RESULTS_DIR, f"drug_disease_predictions_{model_name}.csv"), index=False)

    # Save all metrics as JSON
    serializable_metrics = {}
    for model_name, metrics in all_metrics.items():
        ser = {}
        for k, v in metrics.items():
            if isinstance(v, dict):
                ser[k] = {kk: float(vv) if isinstance(vv, (int, float, np.floating)) else str(vv) for kk, vv in v.items()}
            elif isinstance(v, (int, float, np.floating)):
                ser[k] = float(v)
            else:
                ser[k] = str(v)
        serializable_metrics[model_name] = ser

    with open(os.path.join(RESULTS_DIR, "all_metrics.json"), "w") as f:
        json.dump(serializable_metrics, f, indent=2)

    print("\nAll training and evaluation complete!")
    return all_results, all_metrics, all_covid_preds, all_drug_disease_preds


if __name__ == "__main__":
    main()
