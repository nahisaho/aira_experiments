from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import requests

from pipeline_common import FIGURES_DIR, RESULTS_DIR, append_log, save_json, set_global_seed

API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
UNIPROT_ID = "P0DTD1"
BINDING_RESIDUES = [41, 49, 140, 143, 145, 163, 166, 187]


class DockingSuitabilityAssessor:
    def compute_binding_site_plddt(self, plddt_scores: List[float], binding_residues: List[int]) -> float:
        selected = [plddt_scores[i - 1] for i in binding_residues if 1 <= i <= len(plddt_scores)]
        return float(np.mean(selected)) if selected else float(np.mean(plddt_scores))

    def classify_confidence(self, plddt: float) -> str:
        if plddt >= 90:
            return "high"
        if plddt >= 70:
            return "medium"
        if plddt >= 50:
            return "low"
        return "very_low"

    def compute_docking_suitability_score(self, plddt_array: List[float]) -> Dict[str, float | str]:
        arr = np.asarray(plddt_array, dtype=float)
        metrics = {
            "mean_plddt": float(arr.mean()),
            "median_plddt": float(np.median(arr)),
            "std_plddt": float(arr.std()),
            "min_plddt": float(arr.min()),
            "max_plddt": float(arr.max()),
            "fraction_high_confidence": float(np.mean(arr >= 90.0)),
            "fraction_medium_or_better": float(np.mean(arr >= 70.0)),
            "loop_penalty": float(np.mean(arr < 70.0)),
        }
        score = (
            0.45 * metrics["mean_plddt"]
            + 25.0 * metrics["fraction_high_confidence"]
            + 15.0 * metrics["fraction_medium_or_better"]
            - 10.0 * metrics["loop_penalty"]
            - 0.1 * metrics["std_plddt"]
        )
        metrics["docking_suitability_score"] = float(np.clip(score, 0, 100))
        metrics["overall_confidence"] = self.classify_confidence(metrics["mean_plddt"])
        metrics["docking_recommendation"] = (
            "suitable" if metrics["docking_suitability_score"] >= 75 else "caution" if metrics["docking_suitability_score"] >= 55 else "poor"
        )
        return metrics

    def generate_plddt_report(self, data: Dict) -> Dict:
        plddt_scores = data["plddt_scores"]
        binding_site_plddt = self.compute_binding_site_plddt(plddt_scores, BINDING_RESIDUES)
        summary = self.compute_docking_suitability_score(plddt_scores)
        summary.update(
            {
                "binding_site_plddt": binding_site_plddt,
                "binding_site_confidence": self.classify_confidence(binding_site_plddt),
                "residue_count": len(plddt_scores),
                "binding_residues": BINDING_RESIDUES,
                "source": data.get("source", "mock"),
            }
        )
        return summary


def mock_plddt_array(length: int = 306) -> List[float]:
    rng = np.random.default_rng(42)
    x = np.linspace(0, 1, length)
    baseline = 88 - 18 * np.exp(-((x - 0.12) / 0.08) ** 2) - 12 * np.exp(-((x - 0.82) / 0.1) ** 2)
    noise = rng.normal(0, 4, length)
    arr = np.clip(baseline + noise, 45, 98)
    return arr.round(2).tolist()


def extract_plddt_from_pdb(pdb_text: str) -> List[float]:
    residue_scores = {}
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            try:
                resseq = int(line[22:26].strip())
                bfactor = float(line[60:66].strip())
            except ValueError:
                continue
            residue_scores.setdefault(resseq, []).append(bfactor)
    if not residue_scores:
        return []
    return [float(np.mean(residue_scores[k])) for k in sorted(residue_scores)]


def fetch_alphafold_prediction(uniprot_id: str) -> Dict:
    metadata = {}
    plddt_scores = []
    try:
        response = requests.get(API_URL.format(uniprot_id=uniprot_id), timeout=20)
        response.raise_for_status()
        payload = response.json()
        metadata = payload[0] if isinstance(payload, list) and payload else payload
        pdb_url = metadata.get("pdbUrl") or metadata.get("cifUrl")
        if pdb_url and str(pdb_url).endswith(".pdb"):
            pdb_response = requests.get(pdb_url, timeout=20)
            pdb_response.raise_for_status()
            plddt_scores = extract_plddt_from_pdb(pdb_response.text)
        if not plddt_scores and metadata.get("averagePlddt"):
            avg = float(metadata["averagePlddt"])
            plddt_scores = [avg] * int(metadata.get("sequenceLength", 306))
    except Exception:
        metadata = {}
        plddt_scores = []

    if not plddt_scores:
        plddt_scores = mock_plddt_array()
        metadata = {"uniprot_id": uniprot_id, "averagePlddt": float(np.mean(plddt_scores)), "sequenceLength": len(plddt_scores)}
        source = "mock"
    else:
        source = "AlphaFoldDB"
    return {"metadata": metadata, "plddt_scores": plddt_scores, "source": source}


def make_figure(plddt_scores: List[float], output_path: Path) -> None:
    arr = np.asarray(plddt_scores, dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(14, 8), dpi=150)
    cmap = plt.get_cmap("viridis")
    axes[0].plot(np.arange(1, len(arr) + 1), arr, color=cmap(0.75), lw=1.8)
    axes[0].axhline(90, color=cmap(0.15), ls="--", label="High confidence")
    axes[0].axhline(70, color=cmap(0.45), ls=":", label="Medium confidence")
    axes[0].set_xlabel("Residue index")
    axes[0].set_ylabel("pLDDT")
    axes[0].set_title("AlphaFold pLDDT profile")
    axes[0].legend()
    axes[1].hist(arr, bins=20, color=cmap(0.6), edgecolor="black", alpha=0.9)
    axes[1].set_xlabel("pLDDT")
    axes[1].set_ylabel("Count")
    axes[1].set_title("pLDDT distribution")
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main() -> Dict:
    set_global_seed(42)
    append_log("alphafold_assessment", "run_started", "01_alphafold_assessment.py", {"uniprot_id": UNIPROT_ID})
    assessor = DockingSuitabilityAssessor()
    prediction = fetch_alphafold_prediction(UNIPROT_ID)
    report = assessor.generate_plddt_report(prediction)
    result = {"uniprot_id": UNIPROT_ID, "metadata": prediction["metadata"], "assessment": report, "plddt_scores": prediction["plddt_scores"]}
    json_path = RESULTS_DIR / "alphafold_assessment.json"
    fig_path = FIGURES_DIR / "plddt_distribution.png"
    save_json(json_path, result)
    make_figure(prediction["plddt_scores"], fig_path)
    append_log(
        "alphafold_assessment",
        "run_completed",
        "01_alphafold_assessment.py",
        {"uniprot_id": UNIPROT_ID},
        {"source": prediction["source"], "docking_score": report["docking_suitability_score"]},
        [str(json_path), str(fig_path)],
    )
    print(json.dumps({"status": "ok", "result_file": str(json_path), "figure_file": str(fig_path)}))
    return result


if __name__ == "__main__":
    main()
