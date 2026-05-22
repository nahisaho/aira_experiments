import json
import random
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import cdist
from scipy.special import expit
from scipy.stats import t as student_t
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors


TRUE_EFFECT = 1.0
SEED = 42


def ensure_directories(root: Path) -> None:
    for name in ("figures", "results", "data", "logs"):
        (root / name).mkdir(parents=True, exist_ok=True)


def standardize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return (values - values.mean()) / values.std(ddof=0)


def fit_ols(X: np.ndarray, y: np.ndarray, target_index: int) -> dict:
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    n, p = X.shape
    xtx_inv = np.linalg.pinv(X.T @ X)
    beta = xtx_inv @ X.T @ y
    fitted = X @ beta
    resid = y - fitted
    dof = max(n - p, 1)
    sigma2 = float(resid @ resid / dof)
    cov = sigma2 * xtx_inv
    se = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    crit = float(student_t.ppf(0.975, dof))
    t_value = float(beta[target_index] / max(se[target_index], 1e-12))
    p_value = float(2 * student_t.sf(abs(t_value), dof))
    return {
        "beta": beta,
        "se": se,
        "ci": np.column_stack([beta - crit * se, beta + crit * se]),
        "fitted": fitted,
        "resid": resid,
        "sigma2": sigma2,
        "p_value": p_value,
        "dof": dof,
    }


def simulate_ecological_data(seed: int = SEED) -> dict:
    rng = np.random.default_rng(seed)
    random.seed(seed)

    n_side = 8
    coords = np.array([(i, j) for i in range(n_side) for j in range(n_side)], dtype=float)
    coords += rng.normal(0.0, 0.08, size=coords.shape)
    coords = (coords - coords.min(axis=0)) / (coords.max(axis=0) - coords.min(axis=0))
    n_areas = len(coords)

    distance = cdist(coords, coords)
    kernel = np.exp(-(distance ** 2) / (2 * 0.22 ** 2)) + 1e-6 * np.eye(n_areas)
    spatial_field = standardize(np.linalg.cholesky(kernel) @ rng.normal(size=n_areas))
    area_confounder = standardize(0.8 * spatial_field + rng.normal(scale=0.8, size=n_areas))
    area_sizes = rng.integers(90, 150, size=n_areas)

    area_id = []
    treatment = []
    outcome = []
    confounder = []
    x_coord = []
    y_coord = []

    for area in range(n_areas):
        size = int(area_sizes[area])
        z = rng.normal(
            loc=0.9 * area_confounder[area] + 0.35 * spatial_field[area],
            scale=1.0,
            size=size,
        )
        logits = -0.15 + 1.15 * z + 0.95 * area_confounder[area] + 0.55 * spatial_field[area]
        p_treat = expit(logits)
        t_ind = rng.binomial(1, p_treat, size=size)
        y = (
            TRUE_EFFECT * t_ind
            + 1.25 * z
            + 0.75 * area_confounder[area]
            + 0.85 * spatial_field[area]
            + rng.normal(scale=1.1, size=size)
        )
        area_id.append(np.full(size, area, dtype=int))
        treatment.append(t_ind)
        outcome.append(y)
        confounder.append(z)
        x_coord.append(np.full(size, coords[area, 0]))
        y_coord.append(np.full(size, coords[area, 1]))

    area_id = np.concatenate(area_id)
    treatment = np.concatenate(treatment)
    outcome = np.concatenate(outcome)
    confounder = np.concatenate(confounder)
    x_coord = np.concatenate(x_coord)
    y_coord = np.concatenate(y_coord)

    prevalence = np.array([treatment[area_id == a].mean() for a in range(n_areas)])
    mean_outcome = np.array([outcome[area_id == a].mean() for a in range(n_areas)])
    mean_confounder = np.array([confounder[area_id == a].mean() for a in range(n_areas)])

    return {
        "n_areas": n_areas,
        "coords": coords,
        "spatial_field": spatial_field,
        "area_confounder": area_confounder,
        "area_sizes": area_sizes,
        "area_id": area_id,
        "treatment": treatment,
        "outcome": outcome,
        "confounder": confounder,
        "x_coord": x_coord,
        "y_coord": y_coord,
        "prevalence": prevalence,
        "mean_outcome": mean_outcome,
        "mean_confounder": mean_confounder,
    }


def stratified_analysis(data: dict, seed: int) -> dict:
    treatment = data["treatment"]
    outcome = data["outcome"]
    area_id = data["area_id"]
    confounder = data["confounder"]
    coords = np.column_stack([data["x_coord"], data["y_coord"]])
    area_level = data["area_confounder"][area_id]
    spatial_level = data["spatial_field"][area_id]

    X_ps = np.column_stack([confounder, area_level, spatial_level, coords])
    ps_model = LogisticRegression(max_iter=500)
    ps_model.fit(X_ps, treatment)
    propensity = np.clip(ps_model.predict_proba(X_ps)[:, 1], 1e-4, 1 - 1e-4)

    overlap_mask = (propensity >= 0.05) & (propensity <= 0.95)
    trimmed_treatment = treatment[overlap_mask]
    trimmed_outcome = outcome[overlap_mask]
    trimmed_propensity = propensity[overlap_mask]
    quantiles = np.quantile(trimmed_propensity, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    strata = np.digitize(trimmed_propensity, quantiles[1:-1], right=True)

    weights = []
    estimates = []
    stratum_effects = []
    for stratum in range(5):
        mask = strata == stratum
        y1 = trimmed_outcome[mask & (trimmed_treatment == 1)]
        y0 = trimmed_outcome[mask & (trimmed_treatment == 0)]
        effect = float(y1.mean() - y0.mean())
        weight = float(mask.sum())
        weights.append(weight)
        estimates.append(effect)
        stratum_effects.append({
            "stratum": stratum + 1,
            "n": int(mask.sum()),
            "effect": effect,
        })

    weights = np.array(weights)
    estimates = np.array(estimates)
    weighted_effect = float(np.average(estimates, weights=weights))

    rng = np.random.default_rng(seed)
    boot = []
    indices = np.arange(len(trimmed_outcome))
    for _ in range(250):
        boot_values = []
        boot_weights = []
        for stratum in range(5):
            mask = indices[strata == stratum]
            draw = rng.choice(mask, size=len(mask), replace=True)
            y1 = trimmed_outcome[draw][trimmed_treatment[draw] == 1]
            y0 = trimmed_outcome[draw][trimmed_treatment[draw] == 0]
            if len(y1) == 0 or len(y0) == 0:
                continue
            boot_weights.append(len(draw))
            boot_values.append(float(y1.mean() - y0.mean()))
        if boot_values:
            boot.append(float(np.average(np.array(boot_values), weights=np.array(boot_weights))))
    ci_low, ci_high = np.quantile(boot, [0.025, 0.975])
    return {
        "estimate": weighted_effect,
        "ci": [float(ci_low), float(ci_high)],
        "stratum_effects": stratum_effects,
        "retained_fraction": float(overlap_mask.mean()),
    }


def fit_random_intercept_model(data: dict) -> dict:
    area_id = data["area_id"]
    treatment = data["treatment"]
    outcome = data["outcome"]
    confounder = data["confounder"]
    area_level = data["area_confounder"][area_id]

    X = np.column_stack([np.ones_like(treatment), treatment, confounder, area_level])
    y = outcome
    unique_areas = np.unique(area_id)
    grouped = [(X[area_id == area], y[area_id == area]) for area in unique_areas]

    def gls_terms(var_params: np.ndarray) -> tuple:
        sigma_u2, sigma_e2 = np.exp(var_params)
        xt_vinv_x = np.zeros((X.shape[1], X.shape[1]))
        xt_vinv_y = np.zeros(X.shape[1])
        log_det = 0.0
        for Xi, yi in grouped:
            n_i = len(yi)
            coeff = sigma_u2 / (sigma_e2 * (sigma_e2 + n_i * sigma_u2))
            sum_x = Xi.sum(axis=0)
            xt_vinv_x += (Xi.T @ Xi) / sigma_e2 - coeff * np.outer(sum_x, sum_x)
            xt_vinv_y += (Xi.T @ yi) / sigma_e2 - coeff * sum_x * yi.sum()
            log_det += (n_i - 1) * np.log(sigma_e2) + np.log(sigma_e2 + n_i * sigma_u2)
        beta = np.linalg.solve(xt_vinv_x, xt_vinv_y)
        quad = 0.0
        for Xi, yi in grouped:
            n_i = len(yi)
            coeff = sigma_u2 / (sigma_e2 * (sigma_e2 + n_i * sigma_u2))
            resid = yi - Xi @ beta
            quad += (resid @ resid) / sigma_e2 - coeff * resid.sum() ** 2
        return beta, xt_vinv_x, quad, log_det

    def objective(var_params: np.ndarray) -> float:
        _, _, quad, log_det = gls_terms(var_params)
        return 0.5 * (log_det + quad + len(y) * np.log(2 * np.pi))

    result = minimize(
        objective,
        x0=np.log([0.8, 1.2]),
        method="L-BFGS-B",
        bounds=[(-6.0, 4.0), (-6.0, 4.0)],
    )
    beta, xt_vinv_x, _, _ = gls_terms(result.x)
    cov = np.linalg.pinv(xt_vinv_x)
    se = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    ci = [float(beta[1] - 1.96 * se[1]), float(beta[1] + 1.96 * se[1])]
    return {
        "estimate": float(beta[1]),
        "ci": ci,
        "sigma_u2": float(np.exp(result.x[0])),
        "sigma_e2": float(np.exp(result.x[1])),
        "converged": bool(result.success),
    }


def spatial_basis(coords: np.ndarray, n_basis: int = 6) -> np.ndarray:
    distance = cdist(coords, coords)
    kernel = np.exp(-(distance ** 2) / (2 * 0.28 ** 2))
    n = len(coords)
    H = np.eye(n) - np.ones((n, n)) / n
    centered = H @ kernel @ H
    eigvals, eigvecs = np.linalg.eigh(centered)
    order = np.argsort(eigvals)[::-1]
    keep = order[eigvals[order] > 1e-8][:n_basis]
    return eigvecs[:, keep]


def restricted_spatial_regression(data: dict) -> dict:
    coords = data["coords"]
    y = data["mean_outcome"]
    prevalence = data["prevalence"]
    mean_confounder = data["mean_confounder"]
    basis = spatial_basis(coords, n_basis=6)

    X_base = np.column_stack([np.ones(len(y)), prevalence, mean_confounder])
    non_spatial = fit_ols(X_base, y, target_index=1)

    X_spatial = np.column_stack([X_base, basis])
    spatial_model = fit_ols(X_spatial, y, target_index=1)

    projection = X_base @ np.linalg.pinv(X_base.T @ X_base) @ X_base.T
    restricted_basis = basis - projection @ basis
    keep = np.linalg.norm(restricted_basis, axis=0) > 1e-6
    X_rsr = np.column_stack([X_base, restricted_basis[:, keep]])
    rsr_model = fit_ols(X_rsr, y, target_index=1)

    return {
        "non_spatial": {
            "estimate": float(non_spatial["beta"][1]),
            "ci": non_spatial["ci"][1].astype(float).tolist(),
        },
        "spatial": {
            "estimate": float(spatial_model["beta"][1]),
            "ci": spatial_model["ci"][1].astype(float).tolist(),
        },
        "rsr": {
            "estimate": float(rsr_model["beta"][1]),
            "ci": rsr_model["ci"][1].astype(float).tolist(),
        },
    }


def propensity_score_spatial_matching(data: dict, seed: int) -> dict:
    treatment = data["treatment"]
    outcome = data["outcome"]
    confounder = data["confounder"]
    area_id = data["area_id"]
    area_level = data["area_confounder"][area_id]
    spatial_level = data["spatial_field"][area_id]
    coords = np.column_stack([data["x_coord"], data["y_coord"]])

    X_ps = np.column_stack([confounder, area_level, spatial_level, coords])
    model = LogisticRegression(max_iter=500)
    model.fit(X_ps, treatment)
    propensity = np.clip(model.predict_proba(X_ps)[:, 1], 1e-4, 1 - 1e-4)
    logit_ps = np.log(propensity / (1 - propensity))

    treated_idx = np.where(treatment == 1)[0]
    control_idx = np.where(treatment == 0)[0]
    caliper = 0.25 * np.std(logit_ps)

    control_features = np.column_stack([logit_ps[control_idx], 0.75 * coords[control_idx]])
    treated_features = np.column_stack([logit_ps[treated_idx], 0.75 * coords[treated_idx]])
    neighbors = NearestNeighbors(n_neighbors=min(12, len(control_idx)))
    neighbors.fit(control_features)
    _, indices = neighbors.kneighbors(treated_features)

    matched_pairs = []
    for local_i, neighbor_row in enumerate(indices):
        t_idx = treated_idx[local_i]
        best_score = np.inf
        best_control = None
        for local_j in neighbor_row:
            c_idx = control_idx[local_j]
            ps_gap = abs(logit_ps[t_idx] - logit_ps[c_idx])
            if ps_gap > caliper:
                continue
            spatial_gap = np.linalg.norm(coords[t_idx] - coords[c_idx])
            score = ps_gap + 0.5 * spatial_gap
            if score < best_score:
                best_score = score
                best_control = c_idx
        if best_control is not None:
            matched_pairs.append((t_idx, best_control))

    pair_effects = np.array([outcome[t] - outcome[c] for t, c in matched_pairs], dtype=float)
    estimate = float(pair_effects.mean())
    rng = np.random.default_rng(seed)
    boot = []
    for _ in range(250):
        draw = rng.choice(pair_effects, size=len(pair_effects), replace=True)
        boot.append(float(draw.mean()))
    ci_low, ci_high = np.quantile(boot, [0.025, 0.975])
    return {
        "estimate": estimate,
        "ci": [float(ci_low), float(ci_high)],
        "matched_pairs": int(len(matched_pairs)),
        "propensity_mean": float(propensity.mean()),
    }


def save_fallacy_figure(data: dict, individual_fit: dict, ecological_fit: dict, output_path: Path) -> None:
    rng = np.random.default_rng(SEED)
    sample = rng.choice(len(data["outcome"]), size=min(1800, len(data["outcome"])), replace=False)
    jitter = rng.normal(0.0, 0.035, size=len(sample))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    axes[0].scatter(
        data["treatment"][sample] + jitter,
        data["outcome"][sample],
        s=10,
        alpha=0.25,
        color=plt.get_cmap("viridis")(0.55),
        edgecolor="none",
    )
    x_line = np.array([0.0, 1.0])
    y_line = individual_fit["beta"][0] + individual_fit["beta"][1] * x_line
    axes[0].plot(x_line, y_line, color=plt.get_cmap("cividis")(0.85), linewidth=2.5)
    axes[0].set_title("Individual-level association")
    axes[0].set_xlabel("Individual treatment")
    axes[0].set_ylabel("Outcome")
    axes[0].set_xticks([0, 1])

    scatter = axes[1].scatter(
        data["prevalence"],
        data["mean_outcome"],
        c=data["area_confounder"],
        cmap="cividis",
        s=65,
        edgecolor="black",
        linewidth=0.3,
    )
    x_grid = np.linspace(data["prevalence"].min(), data["prevalence"].max(), 200)
    y_grid = ecological_fit["beta"][0] + ecological_fit["beta"][1] * x_grid
    axes[1].plot(x_grid, y_grid, color=plt.get_cmap("viridis")(0.85), linewidth=2.5)
    axes[1].set_title("Area-level ecological regression")
    axes[1].set_xlabel("Area treatment prevalence")
    axes[1].set_ylabel("Mean outcome")
    fig.colorbar(scatter, ax=axes[1], label="Area confounder")

    fig.savefig(output_path, dpi=320)
    plt.close(fig)


def save_estimate_comparison(results: dict, output_path: Path) -> None:
    methods = [
        "Naive ecological",
        "Stratified",
        "Multilevel",
        "Spatial basis",
        "PS spatial match",
    ]
    estimates = [
        results["naive_ecological"]["estimate"],
        results["stratified"]["estimate"],
        results["multilevel"]["estimate"],
        results["spatially_adjusted"]["estimate"],
        results["propensity_score_spatial_matching"]["estimate"],
    ]
    cis = [
        results["naive_ecological"]["ci"],
        results["stratified"]["ci"],
        results["multilevel"]["ci"],
        results["spatially_adjusted"]["ci"],
        results["propensity_score_spatial_matching"]["ci"],
    ]

    lower = np.array([est - ci[0] for est, ci in zip(estimates, cis)])
    upper = np.array([ci[1] - est for est, ci in zip(estimates, cis)])
    y_pos = np.arange(len(methods))

    fig, ax = plt.subplots(figsize=(8.5, 5.5), constrained_layout=True)
    colors = plt.get_cmap("viridis")(np.linspace(0.2, 0.9, len(methods)))
    ax.errorbar(estimates, y_pos, xerr=np.vstack([lower, upper]), fmt="o", color="black", ecolor="black", capsize=4)
    ax.scatter(estimates, y_pos, s=95, c=colors, zorder=3)
    ax.axvline(TRUE_EFFECT, color=plt.get_cmap("cividis")(0.8), linestyle="--", linewidth=2, label="True effect")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(methods)
    ax.set_xlabel("Estimated treatment effect")
    ax.set_ylabel("Method")
    ax.set_title("Bias-corrected versus naive estimates")
    ax.legend(loc="lower right")
    fig.savefig(output_path, dpi=320)
    plt.close(fig)


def save_spatial_confounding_figure(spatial_results: dict, output_path: Path) -> None:
    methods = ["Non-spatial", "Spatial basis", "Restricted spatial"]
    estimates = [
        spatial_results["non_spatial"]["estimate"],
        spatial_results["spatial"]["estimate"],
        spatial_results["rsr"]["estimate"],
    ]
    cis = [
        spatial_results["non_spatial"]["ci"],
        spatial_results["spatial"]["ci"],
        spatial_results["rsr"]["ci"],
    ]
    lower = np.array([est - ci[0] for est, ci in zip(estimates, cis)])
    upper = np.array([ci[1] - est for est, ci in zip(estimates, cis)])

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    colors = plt.get_cmap("cividis")(np.linspace(0.25, 0.85, len(methods)))
    ax.bar(methods, estimates, color=colors, edgecolor="black")
    ax.errorbar(methods, estimates, yerr=np.vstack([lower, upper]), fmt="none", ecolor="black", capsize=5)
    ax.axhline(TRUE_EFFECT, color=plt.get_cmap("viridis")(0.8), linestyle="--", linewidth=2, label="True effect")
    ax.set_ylabel("Estimated treatment effect")
    ax.set_xlabel("Area-level model")
    ax.set_title("Spatial confounding adjustment")
    ax.legend(loc="upper right")
    fig.savefig(output_path, dpi=320)
    plt.close(fig)


def summarize_to_stdout(results: dict) -> None:
    print("Ecological bias analysis summary")
    print(f"True treatment effect: {TRUE_EFFECT:.3f}")
    print("Method                               Estimate   95% CI")
    print("-" * 62)
    ordered = [
        ("Naive ecological", results["naive_ecological"]),
        ("Stratified", results["stratified"]),
        ("Multilevel", results["multilevel"]),
        ("Spatial basis", results["spatially_adjusted"]),
        ("PS spatial matching", results["propensity_score_spatial_matching"]),
    ]
    for name, info in ordered:
        print(f"{name:30s} {info['estimate']:8.3f}   [{info['ci'][0]:6.3f}, {info['ci'][1]:6.3f}]")
    print(f"Ecological bias magnitude: {results['ecological_bias_magnitude']['absolute_bias']:.3f}")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    ensure_directories(root)

    data = simulate_ecological_data(seed=SEED)
    individual_fit = fit_ols(np.column_stack([np.ones(len(data["treatment"])), data["treatment"]]), data["outcome"], 1)
    ecological_fit = fit_ols(np.column_stack([np.ones(data["n_areas"]), data["prevalence"]]), data["mean_outcome"], 1)
    stratified = stratified_analysis(data, seed=SEED + 1)
    multilevel = fit_random_intercept_model(data)
    spatial_models = restricted_spatial_regression(data)
    psm = propensity_score_spatial_matching(data, seed=SEED + 2)

    results = {
        "true_effect": TRUE_EFFECT,
        "simulation": {
            "seed": SEED,
            "n_individuals": int(len(data["outcome"])),
            "n_areas": int(data["n_areas"]),
            "mean_area_size": float(np.mean(data["area_sizes"])),
            "treatment_prevalence_mean": float(np.mean(data["treatment"])),
        },
        "individual_naive": {
            "estimate": float(individual_fit["beta"][1]),
            "ci": individual_fit["ci"][1].astype(float).tolist(),
        },
        "naive_ecological": {
            "estimate": float(ecological_fit["beta"][1]),
            "ci": ecological_fit["ci"][1].astype(float).tolist(),
            "p_value": float(ecological_fit["p_value"]),
        },
        "stratified": stratified,
        "multilevel": multilevel,
        "spatially_adjusted": spatial_models["spatial"],
        "restricted_spatial": spatial_models["rsr"],
        "spatial_models": spatial_models,
        "propensity_score_spatial_matching": psm,
    }

    results["ecological_bias_magnitude"] = {
        "absolute_bias": float(abs(results["naive_ecological"]["estimate"] - TRUE_EFFECT)),
        "relative_bias_percent": float(100 * abs(results["naive_ecological"]["estimate"] - TRUE_EFFECT) / abs(TRUE_EFFECT)),
        "naive_minus_multilevel": float(results["naive_ecological"]["estimate"] - results["multilevel"]["estimate"]),
    }

    save_fallacy_figure(data, individual_fit, ecological_fit, root / "figures" / "ecological_fallacy_demo.png")
    save_estimate_comparison(results, root / "figures" / "ecological_bias_correction.png")
    save_spatial_confounding_figure(spatial_models, root / "figures" / "spatial_confounding.png")

    with open(root / "results" / "ecological_bias_results.json", "w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)

    summarize_to_stdout(results)


if __name__ == "__main__":
    main()
