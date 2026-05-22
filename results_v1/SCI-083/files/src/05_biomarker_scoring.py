#!/usr/bin/env python3
"""
Module 5: 疾患バイオマーカーの統合スコアリング
Integrated Biomarker Discovery & Scoring Pipeline

Components:
  1. Multi-omic feature selection (LASSO, Random Forest, Boruta)
  2. mixOmics DIABLO-style integration scoring
  3. Composite biomarker panel construction
  4. Cross-validated AUC evaluation
  5. Biomarker interpretation & clinical utility
"""

import os
import json
import logging

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Feature Selection Methods
# ---------------------------------------------------------------------------

class MultiOmicFeatureSelector:
    """Multi-omic feature selection combining multiple methods."""

    def __init__(self, n_top: int = 20):
        self.n_top = n_top
        self.selected_features = {}

    def univariate_selection(self, X: pd.DataFrame, y: np.ndarray,
                             omic_name: str) -> pd.DataFrame:
        """Mann-Whitney U test for case/control comparison"""
        logger.info(f"Univariate selection for {omic_name}")

        results = []
        mask_case = y == 1
        mask_ctrl = y == 0

        for col in X.columns:
            stat, pval = stats.mannwhitneyu(
                X.loc[mask_case, col], X.loc[mask_ctrl, col],
                alternative="two-sided"
            )
            fc = np.mean(X.loc[mask_case, col]) / (np.mean(X.loc[mask_ctrl, col]) + 1e-10)
            auc_single = stat / (mask_case.sum() * mask_ctrl.sum())

            results.append({
                "feature": col,
                "omic": omic_name,
                "statistic": round(stat, 2),
                "pvalue": pval,
                "fold_change": round(fc, 3),
                "auc_single": round(auc_single, 3),
            })

        df = pd.DataFrame(results)
        from statsmodels.stats.multitest import multipletests
        _, fdr, _, _ = multipletests(df["pvalue"], method="fdr_bh")
        df["fdr_qvalue"] = fdr
        df = df.sort_values("pvalue")

        self.selected_features[omic_name] = df.head(self.n_top)["feature"].tolist()
        return df

    def lasso_selection(self, X: np.ndarray, y: np.ndarray,
                        feature_names: list, alpha: float = 0.1) -> list:
        """
        LASSO-based feature selection (simplified coordinate descent).
        本番では sklearn.linear_model.LassoCV を使用。
        """
        logger.info("LASSO feature selection")
        np.random.seed(42)
        n, p = X.shape

        # Standardize
        X_std = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-10)

        # Coordinate descent (simplified)
        beta = np.zeros(p)
        for _ in range(100):
            for j in range(p):
                r = y - X_std @ beta + X_std[:, j] * beta[j]
                z = X_std[:, j] @ r / n
                beta[j] = np.sign(z) * max(abs(z) - alpha, 0)

        selected = [feature_names[i] for i in range(p) if abs(beta[i]) > 1e-4]
        logger.info(f"LASSO selected {len(selected)} features")
        return selected

    def random_forest_importance(self, X: np.ndarray, y: np.ndarray,
                                 feature_names: list, n_trees: int = 100) -> pd.DataFrame:
        """
        Random Forest feature importance (simplified permutation-based).
        本番では sklearn.ensemble.RandomForestClassifier を使用。
        """
        logger.info("Random Forest importance estimation")
        np.random.seed(42)
        n, p = X.shape

        # Simplified: correlation-based importance proxy
        importances = []
        for j in range(p):
            corr = abs(np.corrcoef(X[:, j], y)[0, 1])
            noise = np.random.uniform(0.8, 1.2)
            importances.append(corr * noise)

        df = pd.DataFrame({
            "feature": feature_names,
            "importance": np.round(importances, 4),
        }).sort_values("importance", ascending=False)

        return df


# ---------------------------------------------------------------------------
# DIABLO-style Integration
# ---------------------------------------------------------------------------

class DIABLOIntegration:
    """
    mixOmics DIABLO-inspired multi-omic integration.

    Design matrix specifies expected correlations between omic blocks:
      Taxa ↔ Metabolites: correlated (design = 1)
      Taxa ↔ Clinical: weakly correlated (design = 0.1)
      Metabolites ↔ Clinical: weakly correlated (design = 0.1)
    """

    def __init__(self, design_matrix: np.ndarray = None):
        self.design_matrix = design_matrix
        self.loadings = {}
        self.variates = {}

    def fit(self, blocks: dict, y: np.ndarray, n_components: int = 2,
            keepX: dict = None) -> dict:
        """
        Sparse PLS-DA integration across omic blocks.
        Simplified implementation; 本番では R mixOmics::block.splsda を使用。
        """
        logger.info("Fitting DIABLO multi-omic integration model")

        results = {}
        for block_name, X in blocks.items():
            if isinstance(X, pd.DataFrame):
                X_arr = X.values
                feat_names = X.columns.tolist()
            else:
                X_arr = X
                feat_names = [f"V{i}" for i in range(X_arr.shape[1])]

            # Standardize
            X_std = (X_arr - X_arr.mean(axis=0)) / (X_arr.std(axis=0) + 1e-10)

            # SVD-based component extraction
            U, S, Vt = np.linalg.svd(X_std, full_matrices=False)
            n_comp = min(n_components, len(S))

            loadings = Vt[:n_comp].T
            variates = X_std @ loadings

            # Sparse selection (keep top features per component)
            keep = keepX.get(block_name, 15) if keepX else 15
            for c in range(n_comp):
                abs_loadings = np.abs(loadings[:, c])
                threshold = np.sort(abs_loadings)[-min(keep, len(abs_loadings))]
                loadings[abs_loadings < threshold, c] = 0

            selected = [feat_names[i] for i in range(len(feat_names))
                        if np.any(np.abs(loadings[i]) > 0)]

            self.loadings[block_name] = loadings
            self.variates[block_name] = variates

            results[block_name] = {
                "n_features": X_arr.shape[1],
                "n_selected": len(selected),
                "selected_features": selected[:20],
                "variance_explained": [round(s**2 / np.sum(S**2) * 100, 2) for s in S[:n_comp]],
            }

        results["integration_quality"] = self._assess_integration(blocks, y)
        return results

    def _assess_integration(self, blocks: dict, y: np.ndarray) -> dict:
        """Integration quality metrics"""
        correlations = {}
        block_names = list(self.variates.keys())
        for i in range(len(block_names)):
            for j in range(i + 1, len(block_names)):
                name_i = block_names[i]
                name_j = block_names[j]
                if self.variates[name_i].shape[0] == self.variates[name_j].shape[0]:
                    r, p = stats.pearsonr(
                        self.variates[name_i][:, 0],
                        self.variates[name_j][:, 0]
                    )
                    correlations[f"{name_i} ↔ {name_j}"] = {
                        "correlation": round(r, 4),
                        "pvalue": round(p, 6),
                    }

        return {"cross_block_correlations": correlations}


# ---------------------------------------------------------------------------
# Composite Biomarker Panel
# ---------------------------------------------------------------------------

class BiomarkerPanel:
    """統合バイオマーカーパネルの構築と評価"""

    def __init__(self):
        self.panel_features = []
        self.weights = []

    def build_panel(self, taxa_scores: pd.DataFrame, met_scores: pd.DataFrame,
                    n_taxa: int = 5, n_met: int = 5) -> dict:
        """Top features from each omic to build a composite panel"""
        logger.info("Building composite biomarker panel")

        top_taxa = taxa_scores.head(n_taxa)
        top_met = met_scores.head(n_met)

        panel = []
        for _, row in top_taxa.iterrows():
            panel.append({
                "feature": row["feature"],
                "omic": "microbiome",
                "importance": round(row.get("importance", row.get("auc_single", 0)), 4),
                "direction": "↓" if row.get("fold_change", 1) < 1 else "↑",
            })
        for _, row in top_met.iterrows():
            panel.append({
                "feature": row["feature"],
                "omic": "metabolome",
                "importance": round(row.get("importance", row.get("auc_single", 0)), 4),
                "direction": "↓" if row.get("fold_change", 1) < 1 else "↑",
            })

        self.panel_features = panel
        return panel

    def compute_composite_score(self, X_taxa: np.ndarray, X_met: np.ndarray,
                                taxa_weights: np.ndarray, met_weights: np.ndarray) -> np.ndarray:
        """Weighted composite score"""
        score_taxa = X_taxa @ taxa_weights
        score_met = X_met @ met_weights
        composite = 0.5 * score_taxa + 0.5 * score_met
        return composite

    def evaluate_auc(self, scores: np.ndarray, labels: np.ndarray) -> dict:
        """
        AUC calculation (trapezoidal rule).
        本番では sklearn.metrics.roc_auc_score を使用。
        """
        # Sort by score descending
        sorted_idx = np.argsort(-scores)
        sorted_labels = labels[sorted_idx]

        n_pos = np.sum(labels == 1)
        n_neg = np.sum(labels == 0)
        if n_pos == 0 or n_neg == 0:
            return {"auc": 0.5, "ci_lower": 0.0, "ci_upper": 1.0}

        tpr_list = [0.0]
        fpr_list = [0.0]
        tp = 0
        fp = 0

        for lab in sorted_labels:
            if lab == 1:
                tp += 1
            else:
                fp += 1
            tpr_list.append(tp / n_pos)
            fpr_list.append(fp / n_neg)

        auc = np.trapz(tpr_list, fpr_list)

        # DeLong CI approximation
        se = np.sqrt(auc * (1 - auc) / min(n_pos, n_neg))
        ci_lower = max(0, auc - 1.96 * se)
        ci_upper = min(1, auc + 1.96 * se)

        return {
            "auc": round(auc, 4),
            "ci_lower": round(ci_lower, 4),
            "ci_upper": round(ci_upper, 4),
            "n_positive": int(n_pos),
            "n_negative": int(n_neg),
        }

    def cross_validate(self, X: np.ndarray, y: np.ndarray,
                       n_folds: int = 5) -> dict:
        """K-fold cross-validation AUC"""
        logger.info(f"Running {n_folds}-fold cross-validation")
        np.random.seed(42)
        n = len(y)
        indices = np.random.permutation(n)
        fold_size = n // n_folds

        aucs = []
        for k in range(n_folds):
            test_idx = indices[k * fold_size:(k + 1) * fold_size]
            train_idx = np.setdiff1d(indices, test_idx)

            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            # Simple logistic-like scoring
            w = np.corrcoef(X_train.T, y_train)[-1, :-1]
            scores = X_test @ w
            result = self.evaluate_auc(scores, y_test)
            aucs.append(result["auc"])

        return {
            "mean_auc": round(np.mean(aucs), 4),
            "std_auc": round(np.std(aucs), 4),
            "fold_aucs": [round(a, 4) for a in aucs],
            "n_folds": n_folds,
        }


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_biomarker_scoring_pipeline(output_dir: str = "results") -> dict:
    os.makedirs(output_dir, exist_ok=True)

    # Generate simulated multi-omic data
    np.random.seed(42)
    n_samples = 150
    n_taxa = 30
    n_met = 50

    y = np.array([0]*50 + [1]*50 + [1]*50)  # Control vs IBD
    taxa_names = [f"taxa_{i}" for i in range(n_taxa)]
    met_names = [f"met_{i}" for i in range(n_met)]

    X_taxa = np.random.lognormal(2, 1, (n_samples, n_taxa))
    X_met = np.random.lognormal(8, 2, (n_samples, n_met))

    # Introduce disease-related signals
    for i in range(50, 150):
        X_taxa[i, 0] *= 0.3; X_taxa[i, 1] *= 0.4  # depleted in IBD
        X_taxa[i, 5] *= 2.5; X_taxa[i, 6] *= 2.0  # enriched in IBD
        X_met[i, 0] *= 0.3; X_met[i, 2] *= 0.4    # depleted metabolites
        X_met[i, 5] *= 2.5; X_met[i, 8] *= 2.0    # enriched metabolites

    taxa_df = pd.DataFrame(X_taxa, columns=taxa_names)
    met_df = pd.DataFrame(X_met, columns=met_names)

    # 1. Feature Selection
    selector = MultiOmicFeatureSelector(n_top=15)
    taxa_univariate = selector.univariate_selection(taxa_df, y, "microbiome")
    met_univariate = selector.univariate_selection(met_df, y, "metabolome")

    rf = selector.random_forest_importance(X_taxa, y, taxa_names)

    # 2. DIABLO Integration
    diablo = DIABLOIntegration()
    diablo_results = diablo.fit(
        blocks={"microbiome": taxa_df, "metabolome": met_df},
        y=y,
        n_components=2,
        keepX={"microbiome": 10, "metabolome": 15}
    )

    # 3. Biomarker Panel
    panel = BiomarkerPanel()
    panel_features = panel.build_panel(taxa_univariate, met_univariate, n_taxa=5, n_met=5)

    # 4. Composite score and AUC
    X_combined = np.hstack([X_taxa[:, :5], X_met[:, :5]])
    cv_results = panel.cross_validate(X_combined, y, n_folds=5)

    # Single model AUC
    w = np.corrcoef(X_combined.T, y)[-1, :-1]
    composite_scores = X_combined @ w
    full_auc = panel.evaluate_auc(composite_scores, y)

    # Save results
    taxa_univariate.to_csv(os.path.join(output_dir, "taxa_biomarker_scores.csv"), index=False)
    met_univariate.to_csv(os.path.join(output_dir, "metabolite_biomarker_scores.csv"), index=False)

    all_results = {
        "feature_selection": {
            "taxa_significant": int((taxa_univariate["fdr_qvalue"] < 0.05).sum()),
            "met_significant": int((met_univariate["fdr_qvalue"] < 0.05).sum()),
        },
        "diablo_integration": diablo_results,
        "biomarker_panel": panel_features,
        "full_model_auc": full_auc,
        "cross_validation": cv_results,
    }

    with open(os.path.join(output_dir, "biomarker_scoring_results.json"), "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    summary = {
        "taxa_sig_features": int((taxa_univariate["fdr_qvalue"] < 0.05).sum()),
        "met_sig_features": int((met_univariate["fdr_qvalue"] < 0.05).sum()),
        "panel_size": len(panel_features),
        "full_auc": full_auc["auc"],
        "cv_mean_auc": cv_results["mean_auc"],
        "cv_std_auc": cv_results["std_auc"],
    }

    return summary


if __name__ == "__main__":
    summary = run_biomarker_scoring_pipeline(output_dir="../results")
    print(json.dumps(summary, indent=2))
