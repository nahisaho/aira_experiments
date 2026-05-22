from __future__ import annotations

import json
import random
import re
from collections import Counter
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
)
from sklearn.model_selection import train_test_split
from sklearn.svm import LinearSVC

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = ROOT / "figures"
RESULTS_DIR = ROOT / "results"
DATA_DIR = ROOT / "data"
LOGS_DIR = ROOT / "logs"

DATA_PATH = DATA_DIR / "recall_alerts.csv"
METRICS_PATH = RESULTS_DIR / "module2_metrics.json"
SUMMARY_PATH = RESULTS_DIR / "statistical-summary.md"
REPORT_PATH = ROOT / "report.md"
PREPROCESS_LOG_PATH = DATA_DIR / "preprocessing-log.md"
PROCESS_LOG_PATH = LOGS_DIR / "process-log.jsonl"

FIG_CLASSIFICATION = FIGURES_DIR / "fig2_nlp_classification.png"
FIG_TRENDS = FIGURES_DIR / "fig2b_alert_trends.png"
FIG_ENTITIES = FIGURES_DIR / "fig2c_entity_extraction.png"
FIG_TFIDF = FIGURES_DIR / "fig2d_tfidf_features.png"

CATEGORY_ORDER = ["biological", "chemical", "physical", "allergen"]
SEVERITY_ORDER = ["Class I", "Class II", "Class III"]

STOPWORDS = {
    "a",
    "an",
    "and",
    "are",
    "as",
    "at",
    "batch",
    "because",
    "by",
    "contains",
    "detected",
    "distribution",
    "distributed",
    "due",
    "during",
    "for",
    "found",
    "from",
    "in",
    "into",
    "is",
    "it",
    "its",
    "lot",
    "market",
    "may",
    "of",
    "on",
    "or",
    "potential",
    "presence",
    "product",
    "products",
    "recall",
    "risk",
    "sold",
    "that",
    "the",
    "this",
    "to",
    "undeclared",
    "was",
    "were",
    "with",
}

FDA_STATES = [
    "CA",
    "NY",
    "TX",
    "FL",
    "IL",
    "PA",
    "GA",
    "WA",
    "NJ",
    "MA",
    "OH",
    "MI",
]
RASFF_COUNTRIES = [
    "Spain",
    "Italy",
    "Germany",
    "France",
    "Netherlands",
    "Poland",
    "Belgium",
    "Ireland",
    "Denmark",
    "Sweden",
    "Greece",
    "Portugal",
]

BIOLOGICAL_PATHOGENS = [
    "Salmonella",
    "Listeria monocytogenes",
    "E. coli",
    "Campylobacter",
    "Norovirus",
    "Hepatitis A",
    "Clostridium botulinum",
    "Staphylococcus aureus",
]
CHEMICAL_HAZARDS = [
    "lead",
    "cadmium",
    "ethylene oxide",
    "aflatoxin",
    "histamine",
    "cleaning chemical residue",
    "pesticide residue",
    "PFAS",
]
PHYSICAL_HAZARDS = [
    "glass fragments",
    "metal shavings",
    "hard plastic pieces",
    "rubber fragments",
    "stones",
    "bone fragments",
    "wood splinters",
    "packaging film pieces",
]
ALLERGENS = [
    "milk",
    "peanut",
    "tree nuts",
    "soy",
    "sesame",
    "wheat",
    "egg",
    "shellfish",
    "almond",
]

PRODUCTS_BY_CATEGORY = {
    "biological": [
        "ground turkey",
        "romaine lettuce",
        "soft cheese",
        "frozen berries",
        "deli meat",
        "cantaloupe",
        "raw milk cheese",
        "chicken salad",
        "fresh spinach",
        "bagged salad mix",
    ],
    "chemical": [
        "canned tuna",
        "paprika powder",
        "sunflower oil",
        "herbal supplement",
        "rice noodles",
        "baby food puree",
        "dried chili powder",
        "smoked fish",
        "olive pomace oil",
        "protein powder",
    ],
    "physical": [
        "glass jar pasta sauce",
        "canned soup",
        "frozen pizza",
        "tortilla chips",
        "cereal bar",
        "bottled juice",
        "spice mix",
        "chocolate bar",
        "bread rolls",
        "frozen vegetables",
    ],
    "allergen": [
        "chocolate cookies",
        "ice cream sandwich",
        "granola bar",
        "pesto sauce",
        "curry sauce",
        "sandwich bread",
        "dark chocolate",
        "protein shake",
        "instant soup",
        "muffin mix",
    ],
}

ALL_PRODUCTS = sorted({item for values in PRODUCTS_BY_CATEGORY.values() for item in values}, key=len, reverse=True)
ALL_LOCATIONS = FDA_STATES + RASFF_COUNTRIES

PATHOGEN_PATTERN = re.compile(
    r"\\b(salmonella|listeria monocytogenes|listeria|e\\.?\\s?coli|escherichia coli|campylobacter|norovirus|hepatitis a|clostridium botulinum|staphylococcus aureus)\\b",
    re.IGNORECASE,
)
ALLERGEN_PATTERN = re.compile(
    r"\\b(milk|peanut|peanuts|tree nuts|soy|sesame|wheat|egg|eggs|shellfish|almond|walnut)\\b",
    re.IGNORECASE,
)
PRODUCT_PATTERN = re.compile(
    r"\\b(" + "|".join(re.escape(item.lower()) for item in ALL_PRODUCTS) + r")\\b",
    re.IGNORECASE,
)
LOCATION_PATTERN = re.compile(
    r"\\b(" + "|".join(re.escape(item) for item in sorted(ALL_LOCATIONS, key=len, reverse=True)) + r")\\b",
    re.IGNORECASE,
)

NORMALIZATION_MAP = {
    "listeria": "Listeria monocytogenes",
    "listeria monocytogenes": "Listeria monocytogenes",
    "e coli": "E. coli",
    "e. coli": "E. coli",
    "escherichia coli": "E. coli",
    "peanuts": "peanut",
    "eggs": "egg",
}


def ensure_dirs() -> None:
    for directory in (FIGURES_DIR, RESULTS_DIR, DATA_DIR, LOGS_DIR):
        directory.mkdir(parents=True, exist_ok=True)



def weighted_choice(options: List[str], weights: List[float]) -> str:
    return str(np.random.choice(options, p=np.array(weights) / np.sum(weights)))



def sample_date(index: int, n_samples: int) -> datetime:
    start = datetime(2023, 1, 1)
    end = datetime(2024, 12, 31)
    span_days = (end - start).days
    seasonal_anchor = int((index / max(n_samples - 1, 1)) * span_days)
    jitter = random.randint(-18, 18)
    day_offset = max(0, min(span_days, seasonal_anchor + jitter))
    return start + timedelta(days=day_offset)



def choose_source(index: int) -> str:
    return "FDA" if index % 2 == 0 else "RASFF"



def choose_category(index: int) -> str:
    return CATEGORY_ORDER[index % len(CATEGORY_ORDER)]



def choose_location(source: str) -> str:
    if source == "FDA":
        states = sorted(random.sample(FDA_STATES, k=2))
        return ", ".join(states)
    return random.choice(RASFF_COUNTRIES)



def choose_hazard(category: str, date_value: datetime) -> str:
    if category == "biological":
        weights = np.array([0.24, 0.18, 0.16, 0.09, 0.08, 0.06, 0.05, 0.14])
        if date_value >= datetime(2024, 7, 1):
            weights[1] += 0.18
        if date_value.month in {5, 6, 7, 8}:
            weights[0] += 0.12
        return weighted_choice(BIOLOGICAL_PATHOGENS, weights.tolist())
    if category == "chemical":
        weights = np.array([0.18, 0.1, 0.13, 0.16, 0.12, 0.09, 0.14, 0.08])
        if date_value >= datetime(2024, 9, 1):
            weights[2] += 0.22
        return weighted_choice(CHEMICAL_HAZARDS, weights.tolist())
    if category == "physical":
        weights = [0.23, 0.18, 0.14, 0.1, 0.11, 0.12, 0.06, 0.06]
        return weighted_choice(PHYSICAL_HAZARDS, weights)
    weights = np.array([0.24, 0.16, 0.12, 0.1, 0.14, 0.1, 0.07, 0.04, 0.03])
    if date_value >= datetime(2024, 4, 1):
        weights[4] += 0.2
    return weighted_choice(ALLERGENS, weights.tolist())



def choose_product(category: str) -> str:
    return random.choice(PRODUCTS_BY_CATEGORY[category])



def choose_severity(category: str, hazard: str) -> str:
    if category == "biological":
        if hazard in {"Clostridium botulinum", "Listeria monocytogenes"}:
            return weighted_choice(SEVERITY_ORDER, [0.78, 0.19, 0.03])
        return weighted_choice(SEVERITY_ORDER, [0.62, 0.31, 0.07])
    if category == "chemical":
        if hazard in {"lead", "ethylene oxide", "aflatoxin"}:
            return weighted_choice(SEVERITY_ORDER, [0.33, 0.53, 0.14])
        return weighted_choice(SEVERITY_ORDER, [0.22, 0.59, 0.19])
    if category == "physical":
        if hazard in {"glass fragments", "metal shavings"}:
            return weighted_choice(SEVERITY_ORDER, [0.26, 0.55, 0.19])
        return weighted_choice(SEVERITY_ORDER, [0.15, 0.55, 0.30])
    if hazard in {"milk", "peanut", "sesame", "tree nuts"}:
        return weighted_choice(SEVERITY_ORDER, [0.66, 0.28, 0.06])
    return weighted_choice(SEVERITY_ORDER, [0.52, 0.38, 0.10])



def generate_text(source: str, category: str, hazard: str, product: str, location: str) -> str:
    lot = random.randint(1000, 9999)
    batch = random.randint(10000, 99999)
    if source == "FDA":
        templates = {
            "biological": [
                f"{hazard} contamination detected in {product} product lot #{lot} distributed in {location}.",
                f"FDA recall initiated after routine sampling found {hazard} in refrigerated {product} shipped to {location}.",
                f"Consumer advisory: possible {hazard} in {product} batches #{batch} marketed across {location}.",
            ],
            "chemical": [
                f"Elevated {hazard} levels found in {product} distributed in {location}; voluntary FDA recall announced.",
                f"FDA alert reports {hazard} contamination in {product} lot #{lot} sold in {location}.",
                f"Residual {hazard} detected in imported {product} batch #{batch} with distribution in {location}.",
            ],
            "physical": [
                f"Possible presence of {hazard} in {product} packages shipped to {location}.",
                f"FDA warning for {product} after customer complaint identified {hazard} in lot #{lot} sold in {location}.",
                f"Recall notice: {product} batch #{batch} may contain {hazard} and was distributed in {location}.",
            ],
            "allergen": [
                f"Undeclared allergen ({hazard}) found in {product} distributed in {location}.",
                f"FDA recall issued because {product} label omitted {hazard}; units were sold in {location}.",
                f"Packaging error caused undeclared {hazard} in {product} lot #{lot} shipped to {location}.",
            ],
        }
    else:
        templates = {
            "biological": [
                f"RASFF notification: {hazard} detected in {product} originating from {location}; border alert issued.",
                f"Serious microbiological risk identified as {hazard} in {product} consignment #{batch} from {location}.",
                f"Follow-up RASFF alert reports {hazard} contamination in {product} on the EU market linked to {location}.",
            ],
            "chemical": [
                f"RASFF alert for {product} from {location} due to excessive {hazard} levels in batch #{lot}.",
                f"Chemical hazard notification: {hazard} found in {product} originating from {location}.",
                f"EU control authorities reported {hazard} in {product} imported from {location}.",
            ],
            "physical": [
                f"RASFF information notice: {product} from {location} may contain {hazard}.",
                f"Foreign body alert raised after {hazard} were identified in {product} lot #{lot} from {location}.",
                f"Market withdrawal of {product} linked to {hazard} contamination in shipments from {location}.",
            ],
            "allergen": [
                f"RASFF alert: undeclared {hazard} detected in {product} originating from {location}.",
                f"Incorrect labelling of {product} from {location} led to unlisted {hazard} on pack.",
                f"Border rejection issued for {product} because undeclared {hazard} was identified in batch #{batch} from {location}.",
            ],
        }
    return random.choice(templates[category])



def normalize_entity(value: str) -> str:
    cleaned = value.lower().replace("  ", " ").strip()
    cleaned = NORMALIZATION_MAP.get(cleaned, cleaned)
    if cleaned in FDA_STATES:
        return cleaned
    if cleaned.upper() in FDA_STATES:
        return cleaned.upper()
    if cleaned in {item.lower() for item in RASFF_COUNTRIES}:
        return cleaned.title()
    if cleaned in {item.lower() for item in ALL_PRODUCTS}:
        return cleaned
    if cleaned in {item.lower() for item in ALLERGENS}:
        return cleaned
    if cleaned in {item.lower() for item in BIOLOGICAL_PATHOGENS}:
        for item in BIOLOGICAL_PATHOGENS:
            if item.lower() == cleaned:
                return item
    return cleaned



def extract_entities(text: str) -> Dict[str, List[str]]:
    pathogens = sorted({normalize_entity(match) for match in PATHOGEN_PATTERN.findall(text)})
    allergens = sorted({normalize_entity(match) for match in ALLERGEN_PATTERN.findall(text)})
    products = sorted({normalize_entity(match) for match in PRODUCT_PATTERN.findall(text.lower())})
    locations = sorted({normalize_entity(match) for match in LOCATION_PATTERN.findall(text)})
    return {
        "pathogens": pathogens,
        "allergens": allergens,
        "products": products,
        "locations": locations,
    }



def preprocess_text(text: str) -> List[str]:
    tokens = re.findall(r"[a-zA-Z]+", text.lower())
    return [token for token in tokens if token not in STOPWORDS and len(token) > 2]



def build_vectorizer() -> TfidfVectorizer:
    return TfidfVectorizer(
        tokenizer=preprocess_text,
        preprocessor=None,
        token_pattern=None,
        ngram_range=(1, 2),
        min_df=2,
        max_df=0.95,
    )



def bootstrap_accuracy_ci(y_true: np.ndarray, y_pred: np.ndarray, n_bootstrap: int = 500) -> Tuple[float, float]:
    indices = np.arange(len(y_true))
    scores = []
    for _ in range(n_bootstrap):
        sample_idx = np.random.choice(indices, size=len(indices), replace=True)
        scores.append(accuracy_score(y_true[sample_idx], y_pred[sample_idx]))
    lower, upper = np.percentile(scores, [2.5, 97.5])
    return float(lower), float(upper)



def classification_models() -> Dict[str, object]:
    return {
        "Logistic Regression": LogisticRegression(max_iter=3000, random_state=SEED),
        "Linear SVM": LinearSVC(random_state=SEED),
        "Random Forest": RandomForestClassifier(
            n_estimators=300,
            min_samples_leaf=1,
            class_weight="balanced_subsample",
            random_state=SEED,
        ),
    }



def run_classification_task(df: pd.DataFrame, target_col: str, label_order: List[str]) -> Dict[str, Dict[str, object]]:
    x_train, x_test, y_train, y_test = train_test_split(
        df["text"],
        df[target_col],
        test_size=0.25,
        stratify=df[target_col],
        random_state=SEED,
    )
    vectorizer = build_vectorizer()
    x_train_tfidf = vectorizer.fit_transform(x_train)
    x_test_tfidf = vectorizer.transform(x_test)

    results: Dict[str, Dict[str, object]] = {}
    for name, model in classification_models().items():
        model.fit(x_train_tfidf, y_train)
        predictions = model.predict(x_test_tfidf)
        lower_ci, upper_ci = bootstrap_accuracy_ci(y_test.to_numpy(), np.asarray(predictions))
        results[name] = {
            "accuracy": float(accuracy_score(y_test, predictions)),
            "macro_f1": float(f1_score(y_test, predictions, average="macro")),
            "accuracy_ci95": [lower_ci, upper_ci],
            "confusion_matrix": confusion_matrix(y_test, predictions, labels=label_order).tolist(),
            "classification_report": classification_report(
                y_test,
                predictions,
                labels=label_order,
                output_dict=True,
                zero_division=0,
            ),
        }
    return results



def compute_keyword_zscores(df: pd.DataFrame, keywords: List[str], window: int = 3) -> Tuple[pd.DataFrame, pd.DataFrame, List[Dict[str, object]]]:
    working = df.copy()
    working["month"] = pd.to_datetime(working["date"]).dt.to_period("M").dt.to_timestamp()
    full_months = pd.date_range(working["month"].min(), working["month"].max(), freq="MS")

    monthly_counts = pd.DataFrame(index=full_months)
    for keyword in keywords:
        mask = working["text"].str.contains(re.escape(keyword), case=False, regex=True)
        monthly_counts[keyword] = (
            working.loc[mask].groupby("month").size().reindex(full_months, fill_value=0)
        )

    rolling_mean = monthly_counts.rolling(window=window, min_periods=2).mean().shift(1)
    rolling_std = monthly_counts.rolling(window=window, min_periods=2).std(ddof=0).shift(1).replace(0, np.nan)
    zscores = ((monthly_counts - rolling_mean) / rolling_std).replace([np.inf, -np.inf], np.nan).fillna(0.0)

    emerging = []
    for keyword in keywords:
        peak_month = zscores[keyword].idxmax()
        peak_score = float(zscores.loc[peak_month, keyword])
        emerging.append(
            {
                "keyword": keyword,
                "peak_month": peak_month.strftime("%Y-%m"),
                "max_zscore": peak_score,
                "peak_count": int(monthly_counts.loc[peak_month, keyword]),
            }
        )
    emerging.sort(key=lambda item: item["max_zscore"], reverse=True)
    return monthly_counts, zscores, emerging



def plot_confusion_matrices(results: Dict[str, Dict[str, object]], labels: List[str], output_path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    for ax, (name, metrics) in zip(axes, results.items()):
        disp = ConfusionMatrixDisplay(np.array(metrics["confusion_matrix"]), display_labels=labels)
        disp.plot(ax=ax, cmap="viridis", colorbar=False, values_format="d")
        ax.set_title(f"{name}\nMacro F1 = {metrics['macro_f1']:.2f}")
    fig.suptitle("Severity classification confusion matrices", fontsize=14)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)



def plot_alert_trends(monthly_category: pd.DataFrame, zscores: pd.DataFrame, output_path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True, constrained_layout=True)
    palette = plt.cm.cividis(np.linspace(0.15, 0.85, len(monthly_category.columns)))
    for color, category in zip(palette, monthly_category.columns):
        axes[0].plot(monthly_category.index, monthly_category[category], marker="o", linewidth=2, color=color, label=category.title())
    axes[0].set_title("Alert frequency over time by category")
    axes[0].set_ylabel("Monthly alerts")
    axes[0].legend(frameon=False, ncol=2)
    axes[0].grid(alpha=0.2)

    top_keywords = zscores.max().sort_values(ascending=False).head(5).index.tolist()
    palette_keywords = plt.cm.viridis(np.linspace(0.2, 0.9, len(top_keywords)))
    for color, keyword in zip(palette_keywords, top_keywords):
        axes[1].plot(zscores.index, zscores[keyword], marker="o", linewidth=2, color=color, label=keyword)
    axes[1].axhline(2.0, linestyle="--", color="grey", linewidth=1)
    axes[1].set_title("Emerging risk detection by rolling keyword z-score")
    axes[1].set_ylabel("Rolling z-score")
    axes[1].set_xlabel("Month")
    axes[1].legend(frameon=False, ncol=3)
    axes[1].grid(alpha=0.2)

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)



def plot_entity_counts(entity_counters: Dict[str, Counter], output_path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    plot_specs = [
        ("pathogens", "Top extracted pathogens"),
        ("allergens", "Top extracted allergens"),
        ("products", "Top extracted product types"),
    ]
    for idx, (ax, (entity_name, title)) in enumerate(zip(axes, plot_specs)):
        top_items = entity_counters[entity_name].most_common(8)
        labels = [item[0] for item in top_items][::-1]
        values = [item[1] for item in top_items][::-1]
        color = plt.cm.viridis(np.linspace(0.25, 0.85, len(values)))
        ax.barh(labels, values, color=color)
        ax.set_title(title)
        ax.set_xlabel("Count")
        ax.grid(axis="x", alpha=0.2)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)



def plot_top_tfidf_by_category(df: pd.DataFrame, output_path: Path) -> Dict[str, List[Dict[str, float]]]:
    vectorizer = build_vectorizer()
    matrix = vectorizer.fit_transform(df["text"])
    feature_names = np.array(vectorizer.get_feature_names_out())
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    feature_summary: Dict[str, List[Dict[str, float]]] = {}

    for ax, category in zip(axes.ravel(), CATEGORY_ORDER):
        category_mask = (df["category"] == category).to_numpy()
        mean_scores = np.asarray(matrix[category_mask].mean(axis=0)).ravel()
        top_indices = mean_scores.argsort()[-10:]
        words = feature_names[top_indices]
        scores = mean_scores[top_indices]
        order = np.argsort(scores)
        ax.barh(words[order], scores[order], color=plt.cm.cividis(np.linspace(0.2, 0.85, len(order))))
        ax.set_title(f"{category.title()} top TF-IDF features")
        ax.set_xlabel("Mean TF-IDF")
        ax.grid(axis="x", alpha=0.2)
        feature_summary[category] = [
            {"term": str(term), "mean_tfidf": float(score)} for term, score in zip(words[order][::-1], scores[order][::-1])
        ]

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return feature_summary



def write_markdown_report(
    df: pd.DataFrame,
    severity_results: Dict[str, Dict[str, object]],
    category_results: Dict[str, Dict[str, object]],
    emerging_keywords: List[Dict[str, object]],
) -> None:
    best_severity = max(severity_results.items(), key=lambda item: item[1]["macro_f1"])
    best_category = max(category_results.items(), key=lambda item: item[1]["macro_f1"])
    report = f"""# DRAFT — NOT FOR DISTRIBUTION

## NLP-based early detection of food recalls and alerts

**Timestamp:** {datetime.now(UTC).strftime('%Y-%m-%d %H:%M:%S UTC')}

## Objective
Create a reproducible NLP workflow for synthetic FDA and RASFF recall narratives, including entity extraction, multi-class classification, and keyword trend monitoring.

## Methods
- Generated {len(df)} synthetic recall/alert records with FDA and RASFF-style narratives.
- Applied tokenization, lowercasing, stopword removal, and TF-IDF vectorization.
- Extracted named entities with regex rules for pathogens, allergens, product types, and locations.
- Benchmarked Logistic Regression, Linear SVM, and Random Forest for severity and category prediction.
- Detected emerging risks with rolling 3-month z-scores on keyword frequency.

## Results
- Best severity model: **{best_severity[0]}** (macro F1 = {best_severity[1]['macro_f1']:.3f}, accuracy = {best_severity[1]['accuracy']:.3f}, 95% bootstrap CI = {best_severity[1]['accuracy_ci95'][0]:.3f}–{best_severity[1]['accuracy_ci95'][1]:.3f}).
- Best category model: **{best_category[0]}** (macro F1 = {best_category[1]['macro_f1']:.3f}, accuracy = {best_category[1]['accuracy']:.3f}, 95% bootstrap CI = {best_category[1]['accuracy_ci95'][0]:.3f}–{best_category[1]['accuracy_ci95'][1]:.3f}).
- Top emerging keywords by rolling z-score: {', '.join(f"{item['keyword']} ({item['peak_month']}, z={item['max_zscore']:.2f})" for item in emerging_keywords[:5])}.

## Figures
- `figures/fig2_nlp_classification.png` — severity confusion matrices across models.
- `figures/fig2b_alert_trends.png` — monthly alert trends and rolling keyword z-scores.
- `figures/fig2c_entity_extraction.png` — extracted entity counts.
- `figures/fig2d_tfidf_features.png` — top TF-IDF features by category.

## Discussion
The workflow shows that even simple TF-IDF features can separate alert categories strongly, while severity remains harder because multiple hazard types can map to overlapping classes. Rolling keyword z-scores provide a lightweight early-warning signal for shifting hazard patterns, but synthetic data should not be interpreted as regulatory evidence.

## File inventory
- `src/module2_nlp_detection.py`
- `data/recall_alerts.csv`
- `data/preprocessing-log.md`
- `results/module2_metrics.json`
- `results/statistical-summary.md`
- `figures/fig2_nlp_classification.png`
- `figures/fig2b_alert_trends.png`
- `figures/fig2c_entity_extraction.png`
- `figures/fig2d_tfidf_features.png`
- `logs/process-log.jsonl`
"""
    REPORT_PATH.write_text(report, encoding="utf-8")



def write_preprocessing_log() -> None:
    content = """# Preprocessing log

1. Set `random.seed(42)` and `numpy.random.seed(42)`.
2. Generated synthetic FDA and RASFF recall narratives with structured metadata.
3. Applied lowercase normalization and regex tokenization.
4. Removed custom English stopwords focused on administrative recall terms.
5. Built TF-IDF features with unigram and bigram support (`min_df=2`, `max_df=0.95`).
6. Extracted pathogens, allergens, product types, and locations using regex-based named entity rules.
7. Evaluated severity and category classifiers on stratified holdout sets.
8. Aggregated monthly counts and rolling z-scores for emerging keyword monitoring.
"""
    PREPROCESS_LOG_PATH.write_text(content, encoding="utf-8")



def write_summary_markdown(severity_results: Dict[str, Dict[str, object]], category_results: Dict[str, Dict[str, object]]) -> None:
    lines = ["# Statistical summary", "", "## Severity classification", ""]
    for name, metrics in severity_results.items():
        lines.append(
            f"- **{name}**: accuracy={metrics['accuracy']:.3f}, macro_f1={metrics['macro_f1']:.3f}, "
            f"accuracy 95% CI={metrics['accuracy_ci95'][0]:.3f}–{metrics['accuracy_ci95'][1]:.3f}"
        )
    lines.extend(["", "## Category classification", ""])
    for name, metrics in category_results.items():
        lines.append(
            f"- **{name}**: accuracy={metrics['accuracy']:.3f}, macro_f1={metrics['macro_f1']:.3f}, "
            f"accuracy 95% CI={metrics['accuracy_ci95'][0]:.3f}–{metrics['accuracy_ci95'][1]:.3f}"
        )
    SUMMARY_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")



def to_serializable(value):
    if isinstance(value, dict):
        return {str(key): to_serializable(val) for key, val in value.items()}
    if isinstance(value, list):
        return [to_serializable(item) for item in value]
    if isinstance(value, tuple):
        return [to_serializable(item) for item in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value



def main() -> None:
    ensure_dirs()

    process_events: List[Dict[str, object]] = []

    def log_event(event_type: str, skill_or_tool: str, handoff_in: Dict[str, object], handoff_out: Dict[str, object], files_written: List[str], status: str = "ok") -> None:
        process_events.append(
            {
                "timestamp": datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
                "phase": "PLAN" if event_type in {"run_started", "prompt_received", "skill_selected"} else "EXECUTE" if event_type in {"handoff_started", "handoff_completed", "file_written"} else "REPORT" if event_type == "report_finalized" else "LOG",
                "event_type": event_type,
                "actor": "co-scientist",
                "skill_or_tool": skill_or_tool,
                "handoff_in": handoff_in,
                "handoff_out": handoff_out,
                "files_written": files_written,
                "status": status,
            }
        )

    log_event("run_started", "module2_nlp_detection", {"seed": SEED}, {"cwd": str(ROOT)}, [])
    log_event(
        "prompt_received",
        "co-scientist-data-analysis",
        {"task": "Create NLP-based food recall detection script with synthetic FDA/RASFF data"},
        {"target_files": [str(DATA_PATH), str(METRICS_PATH)]},
        [],
    )
    log_event("skill_selected", "co-scientist-data-analysis", {"reason": "classification and trend analysis"}, {"backend": "matplotlib Agg + sklearn"}, [])

    n_samples = 320
    records = []
    log_event("handoff_started", "synthetic-data-generator", {"n_samples": n_samples}, {}, [])
    for index in range(n_samples):
        date_value = sample_date(index, n_samples)
        source = choose_source(index)
        category = choose_category(index)
        hazard = choose_hazard(category, date_value)
        product = choose_product(category)
        location = choose_location(source)
        severity = choose_severity(category, hazard)
        text = generate_text(source, category, hazard, product, location)
        year_fragment = date_value.strftime("%Y%m%d")
        alert_id = f"{source}-{year_fragment}-{index + 1:04d}"
        records.append(
            {
                "alert_id": alert_id,
                "date": date_value.strftime("%Y-%m-%d"),
                "source": source,
                "text": text,
                "category": category,
                "severity": severity,
                "product_type": product,
                "country_state": location,
            }
        )
    df = pd.DataFrame(records).sort_values("date").reset_index(drop=True)
    df.to_csv(DATA_PATH, index=False)
    log_event("handoff_completed", "synthetic-data-generator", {"n_samples": n_samples}, {"rows": len(df)}, [str(DATA_PATH)])
    log_event("file_written", "pandas.to_csv", {}, {"rows": len(df)}, [str(DATA_PATH)])

    extracted = df["text"].apply(extract_entities)
    entity_counters = {
        "pathogens": Counter(item for values in extracted for item in values["pathogens"]),
        "allergens": Counter(item for values in extracted for item in values["allergens"]),
        "products": Counter(item for values in extracted for item in values["products"]),
        "locations": Counter(item for values in extracted for item in values["locations"]),
    }

    severity_results = run_classification_task(df, "severity", SEVERITY_ORDER)
    category_results = run_classification_task(df, "category", CATEGORY_ORDER)

    monthly_category = (
        df.assign(month=pd.to_datetime(df["date"]).dt.to_period("M").dt.to_timestamp())
        .groupby(["month", "category"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=CATEGORY_ORDER)
    )
    keyword_list = ["salmonella", "listeria", "e. coli", "milk", "sesame", "lead", "ethylene oxide", "glass"]
    keyword_counts, keyword_zscores, emerging_keywords = compute_keyword_zscores(df, keyword_list)

    plot_confusion_matrices(severity_results, SEVERITY_ORDER, FIG_CLASSIFICATION)
    plot_alert_trends(monthly_category, keyword_zscores, FIG_TRENDS)
    plot_entity_counts(entity_counters, FIG_ENTITIES)
    top_tfidf_features = plot_top_tfidf_by_category(df, FIG_TFIDF)
    for path in (FIG_CLASSIFICATION, FIG_TRENDS, FIG_ENTITIES, FIG_TFIDF):
        log_event("file_written", "matplotlib", {}, {"file": str(path)}, [str(path)])

    write_preprocessing_log()
    write_summary_markdown(severity_results, category_results)
    log_event("file_written", "text-writer", {}, {"file": str(PREPROCESS_LOG_PATH)}, [str(PREPROCESS_LOG_PATH)])
    log_event("file_written", "text-writer", {}, {"file": str(SUMMARY_PATH)}, [str(SUMMARY_PATH)])

    metrics = {
        "seed": SEED,
        "sample_count": int(len(df)),
        "date_range": {"start": str(df['date'].min()), "end": str(df['date'].max())},
        "source_distribution": df["source"].value_counts().to_dict(),
        "category_distribution": df["category"].value_counts().to_dict(),
        "severity_distribution": df["severity"].value_counts().to_dict(),
        "severity_classification": severity_results,
        "category_classification": category_results,
        "entity_summary": {
            key: [{"entity": entity, "count": count} for entity, count in counter.most_common(10)]
            for key, counter in entity_counters.items()
        },
        "top_tfidf_features": top_tfidf_features,
        "monthly_category_counts": monthly_category.reset_index().assign(month=lambda frame: frame['month'].dt.strftime('%Y-%m')).to_dict(orient="records"),
        "keyword_monthly_counts": keyword_counts.reset_index().assign(index=lambda frame: frame['index'].dt.strftime('%Y-%m')).rename(columns={"index": "month"}).to_dict(orient="records"),
        "emerging_keywords": emerging_keywords,
    }
    METRICS_PATH.write_text(json.dumps(to_serializable(metrics), indent=2), encoding="utf-8")
    log_event("file_written", "json-writer", {}, {"file": str(METRICS_PATH)}, [str(METRICS_PATH)])

    write_markdown_report(df, severity_results, category_results, emerging_keywords)
    log_event("report_finalized", "report-writer", {}, {"file": str(REPORT_PATH)}, [str(REPORT_PATH)])

    log_event("run_completed", "module2_nlp_detection", {"artifacts": 9}, {"status": "success"}, [str(METRICS_PATH), str(REPORT_PATH)])
    PROCESS_LOG_PATH.write_text("\n".join(json.dumps(to_serializable(event)) for event in process_events) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
