"""
Main pipeline orchestrator for MOF high-throughput screening.

Coordinates all modules:
1. Feature extraction (Zeo++)
2. GCMC simulation (RASPA)
3. Geometric analysis
4. ML prediction
5. Stability filtering
6. DAC ranking
"""
import json
import logging
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from src.config import PipelineConfig
from src.feature_extraction import (
    ChemicalAnalyzer, FeatureExtractionPipeline, MOFDatabaseLoader,
    MOFFeatures, ZeoppRunner,
)
from src.gcmc_simulation import GCMCSimulator, RASPAInputGenerator, RASPAOutputParser
from src.geometric_analysis import GeometricAdsorptionAnalyzer
from src.ml_prediction import AdsorptionPredictor, MultiTargetPredictor
from src.stability_filter import StabilityFilter
from src.dac_ranking import DACCandidate, DACRanker, ParetoFrontAnalysis

logger = logging.getLogger(__name__)


class MOFScreeningPipeline:
    """
    End-to-end MOF screening pipeline for DAC CO2 capture.

    Workflow:
    ┌──────────────────────────────────────────────────────────────┐
    │  Stage 1: Database Loading & Feature Extraction             │
    │  CoRE MOF / hMOF → CIF parsing → Zeo++ geometric analysis  │
    │                  → Chemical composition analysis            │
    ├──────────────────────────────────────────────────────────────┤
    │  Stage 2: Pre-screening (Geometric Filter)                  │
    │  LCD, PLD, ASA, porosity window filtering                   │
    ├──────────────────────────────────────────────────────────────┤
    │  Stage 3: GCMC Simulation (or ML surrogate)                 │
    │  RASPA GCMC → CO2/H2/N2 isotherms                          │
    │  OR: ML prediction from structural descriptors              │
    ├──────────────────────────────────────────────────────────────┤
    │  Stage 4: Stability & Synthesizability Filtering            │
    │  Water stability + thermal stability + synthesizability     │
    ├──────────────────────────────────────────────────────────────┤
    │  Stage 5: DAC Ranking                                       │
    │  Multi-criteria scoring → Pareto front → Top-N candidates   │
    └──────────────────────────────────────────────────────────────┘
    """

    def __init__(self, config: Optional[PipelineConfig] = None):
        self.config = config or PipelineConfig()
        self.start_time = None
        self.log_entries = []

        # Initialize components
        self.db_loader = MOFDatabaseLoader(self.config.data_dir)
        self.zeopp_runner = ZeoppRunner(
            probe_radius=self.config.zeopp.probe_radius_n2,
            n_sa_samples=self.config.zeopp.n_samples_sa,
            n_vol_samples=self.config.zeopp.n_samples_vol,
        )
        self.chem_analyzer = ChemicalAnalyzer()
        self.feature_pipeline = FeatureExtractionPipeline(
            self.zeopp_runner, self.chem_analyzer
        )
        self.gcmc_simulator = GCMCSimulator()
        self.geo_analyzer = GeometricAdsorptionAnalyzer()
        self.stability_filter = StabilityFilter(
            water_threshold=self.config.stability.water_stability_threshold,
            synth_threshold=self.config.stability.synthesizability_threshold,
        )
        self.dac_ranker = DACRanker(
            min_working_capacity=self.config.dac.min_working_capacity,
            min_selectivity=self.config.dac.min_selectivity_co2_n2,
        )
        self.ml_predictor = MultiTargetPredictor(
            model_type=self.config.ml.model_type,
            n_estimators=self.config.ml.n_estimators,
            max_depth=self.config.ml.max_depth,
            learning_rate=self.config.ml.learning_rate,
            random_state=self.config.ml.random_state,
        )

    def run(self, use_ml_surrogate: bool = True) -> Dict:
        """Execute full screening pipeline."""
        self.start_time = time.time()
        self._log("run_started", {"use_ml_surrogate": use_ml_surrogate})

        results = {}

        # Stage 1: Load database and extract features
        logger.info("=" * 60)
        logger.info("STAGE 1: Database Loading & Feature Extraction")
        logger.info("=" * 60)
        structures = self.db_loader.list_structures(
            max_n=self.config.max_structures
        )
        self._log("stage_1_db_loaded", {"n_structures": len(structures)})

        features = self.feature_pipeline.extract_batch(structures)
        self._save_features(features)
        results["n_total"] = len(features)

        # Stage 2: Geometric pre-screening
        logger.info("=" * 60)
        logger.info("STAGE 2: Geometric Pre-screening")
        logger.info("=" * 60)
        criteria = self.geo_analyzer.generate_screening_criteria("CO2")
        mof_dicts = self._features_to_dicts(features)
        filtered_geo = self.geo_analyzer.apply_geometric_filter(mof_dicts, criteria)
        results["n_after_geometric"] = len(filtered_geo)
        self._log("stage_2_geometric_filter", {
            "criteria": {k: list(v) for k, v in criteria.items()},
            "passed": len(filtered_geo),
        })

        # Stage 3: Adsorption prediction
        logger.info("=" * 60)
        logger.info("STAGE 3: Adsorption Prediction")
        logger.info("=" * 60)

        if use_ml_surrogate:
            predictions = self._ml_prediction_stage(features, filtered_geo)
        else:
            predictions = self._gcmc_simulation_stage(filtered_geo)

        results["n_with_predictions"] = len(predictions)

        # Stage 4: Stability filtering
        logger.info("=" * 60)
        logger.info("STAGE 4: Stability & Synthesizability Filtering")
        logger.info("=" * 60)
        passed, failed = self.stability_filter.filter_batch(predictions)
        results["n_after_stability"] = len(passed)
        self._log("stage_4_stability", {
            "passed": len(passed), "failed": len(failed)
        })

        # Stage 5: DAC ranking
        logger.info("=" * 60)
        logger.info("STAGE 5: DAC Ranking")
        logger.info("=" * 60)
        ranked = self.dac_ranker.rank_candidates(passed, top_n=50)
        ranked_filtered = self.dac_ranker.apply_hard_filters(ranked)

        # Save results
        self.dac_ranker.generate_report(
            ranked_filtered,
            self.config.results_dir / "dac_ranking.json"
        )

        # Pareto front analysis
        if len(ranked_filtered) >= 3:
            objectives = np.array([
                [c.working_capacity, c.co2_n2_selectivity, c.water_stability]
                for c in ranked_filtered
            ])
            pareto_mask = ParetoFrontAnalysis.compute_pareto_front(
                objectives, maximize=[True, True, True]
            )
            pareto_candidates = [c for c, m in zip(ranked_filtered, pareto_mask) if m]
            results["n_pareto_optimal"] = len(pareto_candidates)

        results["n_final_candidates"] = len(ranked_filtered)
        results["top_5"] = [c.to_dict() for c in ranked_filtered[:5]]

        elapsed = time.time() - self.start_time
        results["elapsed_seconds"] = round(elapsed, 1)
        self._log("run_completed", {"elapsed": elapsed, "n_final": len(ranked_filtered)})
        self._save_log()

        return results

    def _ml_prediction_stage(self, all_features: List[MOFFeatures],
                              filtered: List[Dict]) -> List[Dict]:
        """Use ML models to predict adsorption properties."""
        self._log("stage_3_ml_prediction", {"n_candidates": len(filtered)})

        # In production, this would load pre-trained models
        # Here we demonstrate the interface
        for mof in filtered:
            # Placeholder predictions based on geometric heuristics
            lcd = mof.get("LCD", 5.0)
            asa = mof.get("ASA", 1000)
            porosity = mof.get("porosity", 0.5)

            # Heuristic predictions (would be replaced by trained model)
            mof["co2_uptake_dac"] = max(0, 0.5 * porosity * np.log(asa / 100 + 1)
                                         * np.exp(-((lcd - 7.0) / 3.0) ** 2))
            mof["co2_uptake_1bar"] = mof["co2_uptake_dac"] * 10
            mof["working_capacity"] = mof["co2_uptake_1bar"] * 0.7
            mof["selectivity"] = max(10, 100 * porosity * (1 - porosity) * 4)
            mof["Qst"] = 25 + 15 * porosity

        return filtered

    def _gcmc_simulation_stage(self, filtered: List[Dict]) -> List[Dict]:
        """Run GCMC simulations for selected candidates."""
        self._log("stage_3_gcmc_simulation", {"n_candidates": len(filtered)})

        for mof in filtered:
            framework = mof.get("mof_id", "")
            pressures = self.config.gcmc.pressure_points

            iso_co2 = self.gcmc_simulator.run_isotherm(
                framework, "CO2", self.config.gcmc.temperature, pressures
            )

            if iso_co2.points:
                mof["co2_uptake_dac"] = float(np.interp(
                    0.0004, iso_co2.pressures(), iso_co2.loadings()
                ))
                mof["co2_uptake_1bar"] = float(np.interp(
                    1.0, iso_co2.pressures(), iso_co2.loadings()
                ))
                mof["working_capacity"] = self.gcmc_simulator.compute_working_capacity(
                    iso_co2, 1.0, 0.1
                )
                mof["Qst"] = iso_co2.heat_of_adsorption
                mof["selectivity"] = 100  # Would compute with N2 isotherm

        return filtered

    def _features_to_dicts(self, features: List[MOFFeatures]) -> List[Dict]:
        """Convert MOFFeatures to flat dictionaries for filtering."""
        dicts = []
        for f in features:
            d = {
                "mof_id": f.mof_id,
                "source_db": f.source_db,
                "LCD": f.geometric.lcd,
                "PLD": f.geometric.pld,
                "ASA": f.geometric.asa,
                "porosity": f.geometric.porosity,
                "AV": f.geometric.av,
                "density": f.geometric.density,
                "metal_type": f.chemical.metal_type,
                "linker_type": f.chemical.linker_type,
                "has_oms": f.chemical.has_open_metal_sites,
                "n_atom_types": f.chemical.n_atom_types,
                "n_atoms_per_uc": f.chemical.n_atoms_per_uc,
            }
            dicts.append(d)
        return dicts

    def _save_features(self, features: List[MOFFeatures]):
        """Save extracted features to CSV-like format."""
        output = self.config.results_dir / "extracted_features.json"
        data = []
        for f in features:
            data.append({
                "mof_id": f.mof_id,
                "source_db": f.source_db,
                **{name: float(val) for name, val in
                   zip(f.feature_names(), f.to_feature_vector())}
            })
        with open(output, "w") as fp:
            json.dump(data, fp, indent=2)

    def _log(self, event_type: str, data: Dict = None):
        entry = {
            "timestamp": datetime.now().isoformat(),
            "phase": "pipeline",
            "event_type": event_type,
            "actor": "co-scientist",
            "skill_or_tool": "mof-screening-pipeline",
            "data": data or {},
        }
        self.log_entries.append(entry)

    def _save_log(self):
        log_path = self.config.logs_dir / "process-log.jsonl"
        with open(log_path, "a") as f:
            for entry in self.log_entries:
                f.write(json.dumps(entry) + "\n")
        self.log_entries = []


def main():
    """Entry point for pipeline execution."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    )

    config = PipelineConfig()
    pipeline = MOFScreeningPipeline(config)
    results = pipeline.run(use_ml_surrogate=True)

    print("\n" + "=" * 60)
    print("MOF SCREENING PIPELINE — RESULTS SUMMARY")
    print("=" * 60)
    print(f"Total structures scanned:  {results.get('n_total', 0)}")
    print(f"After geometric filter:    {results.get('n_after_geometric', 0)}")
    print(f"After stability filter:    {results.get('n_after_stability', 0)}")
    print(f"Final DAC candidates:      {results.get('n_final_candidates', 0)}")
    print(f"Pareto-optimal:            {results.get('n_pareto_optimal', 0)}")
    print(f"Elapsed time:              {results.get('elapsed_seconds', 0)} s")

    if results.get("top_5"):
        print("\nTop 5 DAC Candidates:")
        for c in results["top_5"]:
            print(f"  #{c['rank']} {c['mof_id']} — "
                  f"WC={c['working_capacity_mmol_g']:.2f} mmol/g, "
                  f"Score={c['overall_score']:.4f}")


if __name__ == "__main__":
    main()
