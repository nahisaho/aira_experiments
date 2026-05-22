"""
Data Pipeline Orchestrator.
Coordinates all modules into a unified real-time pipeline with caching, validation, and scheduling.
"""

import json
import os
import sys
import time
import hashlib
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

warnings.filterwarnings("ignore")

# Add src to path for imports
_SRC_DIR = Path(__file__).parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

import numpy as np
import pandas as pd


def _load_config(config_path: str) -> Dict:
    """Load YAML config file."""
    if HAS_YAML and os.path.exists(config_path):
        with open(config_path) as f:
            return yaml.safe_load(f)
    return {}


def _log_event(log_path: str, phase: str, event_type: str,
               skill: str, handoff_in: Dict, handoff_out: Dict,
               files_written: list, status: str = "ok"):
    """Append a structured event to the process log."""
    entry = {
        "timestamp": datetime.now().isoformat(),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill,
        "handoff_in": handoff_in,
        "handoff_out": handoff_out,
        "files_written": files_written,
        "status": status,
    }
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "a") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")


class DataValidator:
    """Validates data quality for pipeline inputs."""

    def validate_dataframe(self, df: pd.DataFrame, required_cols: list,
                            min_rows: int = 1) -> Dict:
        """Check that a DataFrame has required columns and sufficient rows."""
        issues = []
        if df is None or len(df) == 0:
            issues.append("Empty dataframe")
            return {"valid": False, "issues": issues, "n_rows": 0}

        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            issues.append(f"Missing columns: {missing}")

        if len(df) < min_rows:
            issues.append(f"Too few rows: {len(df)} < {min_rows}")

        null_fractions = df.isnull().mean()
        high_null = null_fractions[null_fractions > 0.5].index.tolist()
        if high_null:
            issues.append(f"High null rate (>50%) in: {high_null}")

        return {
            "valid": len(issues) == 0,
            "issues": issues,
            "n_rows": len(df),
            "n_cols": len(df.columns),
        }

    def validate_numeric_range(self, value: float, name: str,
                                low: float, high: float) -> Dict:
        valid = low <= value <= high
        return {"valid": valid, "name": name, "value": value,
                "expected": f"[{low}, {high}]"}


class PandemicDataPipeline:
    """Orchestrates the full pandemic surveillance pipeline."""

    def __init__(self, config_path: str = "configs/pipeline_config.yaml"):
        self.config = _load_config(config_path)
        self.validator = DataValidator()
        self.workspace = Path(__file__).parent.parent
        self.log_path = str(self.workspace / "logs" / "process-log.jsonl")
        self.results_dir = self.workspace / "results"
        self.results_dir.mkdir(exist_ok=True)

    def run_full_pipeline(self, override_config: Optional[Dict] = None) -> Dict:
        """Run all pipeline modules and aggregate results."""
        cfg = {**self.config, **(override_config or {})}
        t_start = time.time()

        _log_event(self.log_path, "pipeline", "run_started",
                   "PandemicDataPipeline", {}, {}, [])

        print("[1/6] Running genomic surveillance...")
        genomic_result = self._run_genomic(cfg)

        print("[2/6] Running mutation hotspot analysis...")
        mutation_result = self._run_mutation(genomic_result, cfg)

        print("[3/6] Running epidemiology integration...")
        epi_result = self._run_epidemiology(cfg)

        print("[4/6] Running Rt estimation...")
        rt_result = self._run_rt(epi_result, cfg)

        print("[5/6] Running NLP alert analysis...")
        nlp_result = self._run_nlp(cfg)

        print("[6/6] Computing risk scores and generating alerts...")
        risk_result = self._run_risk(genomic_result, epi_result, rt_result, nlp_result, cfg)

        elapsed = round(time.time() - t_start, 2)
        results = {
            "genomic": genomic_result,
            "mutation": mutation_result,
            "epidemiology": epi_result,
            "rt_estimation": rt_result,
            "nlp": nlp_result,
            "risk_scoring": risk_result,
            "pipeline_metadata": {
                "elapsed_seconds": elapsed,
                "timestamp": datetime.now().isoformat(),
                "config": {k: v for k, v in cfg.items() if k != "sequences_df"},
            },
        }

        self._save_results(results)
        _log_event(self.log_path, "pipeline", "run_completed",
                   "PandemicDataPipeline", {}, {"elapsed_sec": elapsed}, [])
        print(f"\nPipeline completed in {elapsed}s")
        return results

    def _run_genomic(self, cfg: Dict) -> Dict:
        from genomic_surveillance import run_genomic_surveillance
        result = run_genomic_surveillance(cfg.get("genomics", {}))
        _log_event(self.log_path, "genomic_surveillance", "handoff_completed",
                   "genomic_surveillance",
                   {}, {"n_sequences": result["n_sequences_total"],
                        "n_emerging": result["n_novel_emerging"]}, [])
        return result

    def _run_mutation(self, genomic_result: Dict, cfg: Dict) -> Dict:
        from mutation_hotspot import run_mutation_analysis
        seqs_df = genomic_result.get("sequences_df")
        result = run_mutation_analysis(seqs_df)
        _log_event(self.log_path, "mutation_hotspot", "handoff_completed",
                   "mutation_hotspot",
                   {}, {"n_hotspots": result["n_hotspot_positions"],
                        "top_mutations": result["top_mutations"]}, [])
        return result

    def _run_epidemiology(self, cfg: Dict) -> Dict:
        from epidemiology_integration import run_epidemiology_pipeline
        result = run_epidemiology_pipeline(cfg.get("epidemiology", {}))
        _log_event(self.log_path, "epidemiology", "handoff_completed",
                   "epidemiology_integration",
                   {}, {"composite_score": result["integrated_signals"]["composite_epi_score"],
                        "n_anomalies": result["n_anomalies_detected"]}, [])
        return result

    def _run_rt(self, epi_result: Dict, cfg: Dict) -> Dict:
        from rt_estimation import run_rt_estimation
        cases_df = epi_result.get("cases_df")
        result = run_rt_estimation(cases_df)
        _log_event(self.log_path, "rt_estimation", "handoff_completed",
                   "rt_estimation",
                   {}, {"global_rt": result["global_rt_mean"],
                        "n_growing": result["n_countries_growing"]}, [])
        return result

    def _run_nlp(self, cfg: Dict) -> Dict:
        from nlp_alert_analysis import run_nlp_pipeline
        result = run_nlp_pipeline()
        _log_event(self.log_path, "nlp_analysis", "handoff_completed",
                   "nlp_alert_analysis",
                   {}, {"signal_score": result["signal_score"],
                        "n_alerts": result["n_alerts_processed"]}, [])
        return result

    def _run_risk(self, genomic: Dict, epi: Dict, rt: Dict,
                   nlp: Dict, cfg: Dict) -> Dict:
        from risk_scoring import run_risk_scoring
        weights = cfg.get("risk_scoring", {}).get("weights", {})
        result = run_risk_scoring(genomic, epi, rt, nlp, config={"weights": weights} if weights else None)
        _log_event(self.log_path, "risk_scoring", "handoff_completed",
                   "risk_scoring",
                   {}, {"composite_score": result["risk_result"]["composite_score"],
                        "risk_level": result["risk_result"]["risk_level"]}, [])
        return result

    def _save_results(self, results: Dict):
        """Persist key results to CSV/JSON in results/."""
        rd = self.results_dir

        # ── Rt estimates CSV ──
        rt_df = results["rt_estimation"].get("rt_summary_df", pd.DataFrame())
        if len(rt_df):
            path = str(rd / "rt_estimates.csv")
            rt_df.to_csv(path, index=False)
            _log_event(self.log_path, "io", "file_written", "pipeline",
                       {}, {}, [path])

        # ── Risk scores CSV ──
        risk_r = results["risk_scoring"]
        risk_path = str(rd / "risk_scores.csv")
        pd.DataFrame([{
            "timestamp": datetime.now().isoformat(),
            **risk_r["component_scores"],
            "composite_score": risk_r["risk_result"]["composite_score"],
            "risk_level": risk_r["risk_result"]["risk_level"],
        }]).to_csv(risk_path, index=False)
        _log_event(self.log_path, "io", "file_written", "pipeline", {}, {}, [risk_path])

        # ── Mutation hotspots CSV ──
        mut_df = results["mutation"].get("mutations_df", pd.DataFrame())
        if len(mut_df):
            mut_path = str(rd / "mutation_hotspots.csv")
            mut_df.head(20).to_csv(mut_path, index=False)
            _log_event(self.log_path, "io", "file_written", "pipeline", {}, {}, [mut_path])

        # ── Alert log CSV ──
        alert = results["risk_scoring"].get("alert", {})
        alert_path = str(rd / "alert_log.csv")
        pd.DataFrame([{
            "alert_id": alert.get("alert_id", ""),
            "timestamp": alert.get("timestamp", ""),
            "risk_score": alert.get("risk_score", 0),
            "risk_level": alert.get("risk_level", ""),
            "action_required": alert.get("action_required", False),
            "n_evidence": len(alert.get("evidence_summary", [])),
        }]).to_csv(alert_path, index=False)
        _log_event(self.log_path, "io", "file_written", "pipeline", {}, {}, [alert_path])

        # ── Component metrics JSON ──
        metrics = {
            "genomic": {
                "n_sequences": results["genomic"]["n_sequences_total"],
                "n_emerging_lineages": results["genomic"]["n_novel_emerging"],
                "evolutionary_rate": results["genomic"]["evolutionary_rate"],
                "diversity": results["genomic"]["diversity_metrics"],
            },
            "mutation": {
                "n_hotspot_positions": results["mutation"]["n_hotspot_positions"],
                "top_mutations": results["mutation"]["top_mutations"],
                "immune_escape_stats": results["mutation"]["escape_stats"],
            },
            "epidemiology": {
                "integrated_signals": results["epidemiology"]["integrated_signals"],
                "n_anomalies": results["epidemiology"]["n_anomalies_detected"],
                "ww_correlation": results["epidemiology"]["ww_correlation"],
            },
            "rt": {
                "global_rt_mean": results["rt_estimation"]["global_rt_mean"],
                "n_countries_growing": results["rt_estimation"]["n_countries_growing"],
                "latest_rts": results["rt_estimation"]["latest_rts"],
            },
            "nlp": {
                "signal_score": results["nlp"]["signal_score"],
                "severity_distribution": results["nlp"]["severity_distribution"],
                "n_alerts": results["nlp"]["n_alerts_processed"],
            },
            "risk": {
                "composite_score": results["risk_scoring"]["risk_result"]["composite_score"],
                "risk_level": results["risk_scoring"]["risk_result"]["risk_level"],
                "roc_auc": results["risk_scoring"]["threshold_optimization"]["roc_auc"],
                "optimal_threshold": results["risk_scoring"]["threshold_optimization"][
                    "optimal_threshold_95sens"],
                "uncertainty": results["risk_scoring"]["risk_result"].get("uncertainty", {}),
            },
        }
        metrics_path = str(rd / "component_metrics.json")
        with open(metrics_path, "w") as f:
            json.dump(metrics, f, indent=2, ensure_ascii=False)
        _log_event(self.log_path, "io", "file_written", "pipeline", {}, {}, [metrics_path])

        # ── Full pipeline JSON (serializable subset) ──
        pipeline_json_path = str(rd / "pipeline_results.json")
        serializable = {
            "metadata": results["pipeline_metadata"],
            "metrics": metrics,
            "alert": results["risk_scoring"].get("alert", {}),
        }
        with open(pipeline_json_path, "w") as f:
            json.dump(serializable, f, indent=2, ensure_ascii=False, default=str)
        _log_event(self.log_path, "io", "file_written", "pipeline",
                   {}, {}, [pipeline_json_path])

        print(f"  → Results saved to: {rd}")


if __name__ == "__main__":
    pipeline = PandemicDataPipeline()
    results = pipeline.run_full_pipeline({})
    risk = results["risk_scoring"]["risk_result"]
    print(f"\nFinal Risk Score: {risk['composite_score']:.1f} ({risk['risk_level']})")
    print(f"Alert: {results['risk_scoring']['alert']['alert_id']}")
