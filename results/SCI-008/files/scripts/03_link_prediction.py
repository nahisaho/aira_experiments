#!/usr/bin/env python3
"""
Step 3: Link Prediction for Drug-Disease Associations & COVID-19 Case Study
"""

import json
import os
from datetime import datetime

import numpy as np
import pandas as pd
import torch
from pykeen.triples import TriplesFactory
from pykeen.models import Model

np.random.seed(42)
torch.manual_seed(42)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
LOG_DIR = os.path.join(BASE_DIR, "logs")


def load_data():
    df = pd.read_csv(os.path.join(DATA_DIR, "triples.tsv"), sep="\t")
    triples = df[["head", "relation", "tail"]].values
    tf = TriplesFactory.from_labeled_triples(triples)

    with open(os.path.join(DATA_DIR, "entities.json")) as f:
        entities = json.load(f)

    return df, tf, entities


def predict_drug_disease(tf, model, entities, relation="treats", target_disease=None):
    """Predict drug-disease associations using the trained model."""
    drugs = [eid for eid, info in entities.items() if info["type"] == "Drug"]
    diseases = [eid for eid, info in entities.items() if info["type"] == "Disease"]

    if target_disease:
        diseases = [target_disease] if target_disease in diseases else diseases

    predictions = []

    for drug in drugs:
        for disease in diseases:
            if drug not in tf.entity_to_id or disease not in tf.entity_to_id:
                continue
            if relation not in tf.relation_to_id:
                continue

            h_id = tf.entity_to_id[drug]
            r_id = tf.relation_to_id[relation]
            t_id = tf.entity_to_id[disease]

            h_tensor = torch.tensor([h_id], dtype=torch.long)
            r_tensor = torch.tensor([r_id], dtype=torch.long)
            t_tensor = torch.tensor([t_id], dtype=torch.long)

            with torch.no_grad():
                score = model.score_hrt(
                    torch.stack([h_tensor, r_tensor, t_tensor], dim=1)
                ).item()

            predictions.append({
                "drug_id": drug,
                "drug_name": entities[drug]["name"],
                "disease_id": disease,
                "disease_name": entities[disease]["name"],
                "score": score,
            })

    pred_df = pd.DataFrame(predictions)
    if len(pred_df) > 0:
        pred_df = pred_df.sort_values("score", ascending=False).reset_index(drop=True)
    return pred_df


def find_explanatory_paths(triples_df, entities, drug_id, disease_id, max_depth=3):
    """Find paths from drug to disease through the KG for explainability."""
    import networkx as nx

    G = nx.DiGraph()
    for _, row in triples_df.iterrows():
        G.add_edge(row["head"], row["tail"], relation=row["relation"])

    paths = []
    try:
        for path in nx.all_simple_paths(G, drug_id, disease_id, cutoff=max_depth):
            path_info = []
            for i in range(len(path) - 1):
                edge_data = G.edges[path[i], path[i + 1]]
                node_name = entities.get(path[i], {}).get("name", path[i])
                next_name = entities.get(path[i + 1], {}).get("name", path[i + 1])
                path_info.append({
                    "from": path[i],
                    "from_name": node_name,
                    "relation": edge_data["relation"],
                    "to": path[i + 1],
                    "to_name": next_name,
                })
            paths.append(path_info)
            if len(paths) >= 20:
                break
    except nx.NetworkXError:
        pass

    return paths


def covid_case_study(tf, model, entities, triples_df):
    """COVID-19 drug repurposing case study."""
    print("\n=== COVID-19 Drug Repurposing Case Study ===")

    covid_id = "DOID:0080600"

    # Predict drugs for COVID-19
    pred_df = predict_drug_disease(tf, model, entities, "treats", covid_id)

    # Get known COVID-19 drugs from the KG
    known_covid_drugs = set(
        triples_df[
            (triples_df["relation"] == "treats") & (triples_df["tail"] == covid_id)
        ]["head"].values
    )

    if len(pred_df) > 0:
        pred_df["known"] = pred_df["drug_id"].isin(known_covid_drugs)
        pred_df["rank"] = range(1, len(pred_df) + 1)

    print(f"\nKnown COVID-19 drugs in KG: {len(known_covid_drugs)}")
    print(f"Total predictions: {len(pred_df)}")
    print("\n--- Top 20 Predicted Drugs for COVID-19 ---")
    if len(pred_df) >= 20:
        top20 = pred_df.head(20)
    else:
        top20 = pred_df
    print(top20[["rank", "drug_name", "score", "known"]].to_string(index=False))

    # Novel predictions (not known)
    novel = pred_df[~pred_df["known"]].head(10) if len(pred_df) > 0 else pd.DataFrame()
    print("\n--- Top 10 Novel Predictions for COVID-19 ---")
    if len(novel) > 0:
        print(novel[["rank", "drug_name", "score"]].to_string(index=False))

    # Path explanations for top novel predictions
    explanations = {}
    if len(novel) > 0:
        for _, row in novel.head(5).iterrows():
            drug_id = row["drug_id"]
            paths = find_explanatory_paths(triples_df, entities, drug_id, covid_id)
            explanations[row["drug_name"]] = paths
            print(f"\n  Paths: {row['drug_name']} → COVID-19 ({len(paths)} paths found)")
            for i, path in enumerate(paths[:3]):
                path_str = " → ".join(
                    [f"{s['from_name']} --[{s['relation']}]--> {s['to_name']}" for s in path]
                )
                print(f"    Path {i+1}: {path_str}")

    return pred_df, novel, explanations


def main():
    print("=== Link Prediction & COVID-19 Case Study ===")
    triples_df, tf, entities = load_data()

    # Load best model (try each in priority order)
    best_model = None
    best_model_name = None
    for model_name in ["RotatE", "ComplEx", "TransE"]:
        model_dir = os.path.join(RESULTS_DIR, f"model_{model_name.lower()}")
        model_path = os.path.join(model_dir, "trained_model.pkl")
        if os.path.exists(model_path):
            print(f"Loading model: {model_name}")
            best_model = torch.load(model_path, map_location="cpu", weights_only=False)
            best_model_name = model_name
            break

    if best_model is None:
        print("ERROR: No trained model found. Run 02_train_embeddings.py first.")
        return

    print(f"Using model: {best_model_name}")

    # Full drug-disease predictions
    all_pred_df = predict_drug_disease(tf, best_model, entities, "treats")
    all_pred_df.to_csv(os.path.join(RESULTS_DIR, "all_drug_disease_predictions.csv"), index=False)
    print(f"\nTotal drug-disease predictions: {len(all_pred_df)}")

    # COVID-19 case study
    covid_pred_df, novel_df, explanations = covid_case_study(
        tf, best_model, entities, triples_df
    )

    # Save results
    if len(covid_pred_df) > 0:
        covid_pred_df.to_csv(os.path.join(RESULTS_DIR, "covid19_predictions.csv"), index=False)
    if len(novel_df) > 0:
        novel_df.to_csv(os.path.join(RESULTS_DIR, "covid19_novel_predictions.csv"), index=False)

    with open(os.path.join(RESULTS_DIR, "covid19_path_explanations.json"), "w") as f:
        json.dump(explanations, f, indent=2, default=str)

    # Log
    log_entry = {
        "timestamp": datetime.now().isoformat(),
        "phase": "link_prediction",
        "event_type": "covid19_case_study",
        "actor": "co-scientist",
        "skill_or_tool": "03_link_prediction.py",
        "handoff_out": {
            "model_used": best_model_name,
            "total_predictions": len(all_pred_df),
            "covid_predictions": len(covid_pred_df) if len(covid_pred_df) > 0 else 0,
            "novel_predictions": len(novel_df) if len(novel_df) > 0 else 0,
        },
        "files_written": [
            "results/all_drug_disease_predictions.csv",
            "results/covid19_predictions.csv",
            "results/covid19_novel_predictions.csv",
            "results/covid19_path_explanations.json",
        ],
        "status": "ok",
    }
    with open(os.path.join(LOG_DIR, "process-log.jsonl"), "a") as f:
        f.write(json.dumps(log_entry, default=str) + "\n")

    print("\n=== Link prediction complete ===")


if __name__ == "__main__":
    main()
