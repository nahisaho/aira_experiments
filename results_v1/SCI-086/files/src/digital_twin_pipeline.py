"""
digital_twin_pipeline.py
========================
Main orchestrator for the cardiac digital twin framework.
Integrates all modules into a unified pipeline.
"""

import numpy as np
import json
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional

from src.preprocessing.cardiac_mri_segmentation import (
    MRIVolume, CardiacSegmentationPipeline, MeshGenerator,
    export_to_opencarp, export_to_febio
)
from src.electrophysiology.electrophysiology_models import (
    AlievPanfilovModel, TenTusscherModel, MonodomainSolver,
    TissueParams, IonicModelType, generate_opencarp_config
)
from src.mechanics.cardiac_mechanics import (
    HolzapfelOgdenParams, ActiveTensionParams, WindkesselParams,
    PassiveMechanicsModel, ActiveTensionModel, ElectroMechanicalCoupling,
    CardiacCycleSimulator, generate_febio_mechanics_config
)
from src.inverse.inverse_estimation import (
    ECGData, EchoData, ECGInverseSolver, MechanicsInverseSolver
)
from src.arrhythmia.arrhythmia_risk import (
    ReentryVulnerabilityAnalysis, APDDispersionAnalyzer,
    FibrosisMappingAnalyzer, ArrhythmiaRiskAssessor
)
from src.ablation.ablation_prediction import AblationCaseStudy

logger = logging.getLogger(__name__)


class CardiacDigitalTwinPipeline:
    """
    End-to-end cardiac digital twin construction and simulation pipeline.

    Workflow:
    1. Image Processing: MRI → Segmentation → 3D Mesh
    2. Electrophysiology: Ionic model → Tissue simulation
    3. Mechanics: Passive + Active → Cardiac cycle
    4. Personalization: Inverse problem (ECG + Echo)
    5. Risk Assessment: Arrhythmia vulnerability
    6. Intervention: Ablation strategy prediction
    """

    def __init__(self, patient_id: str, output_dir: str = "."):
        self.patient_id = patient_id
        self.output_dir = Path(output_dir)
        self.results = {}
        self.timestamp = datetime.now().isoformat()

        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "results").mkdir(exist_ok=True)
        (self.output_dir / "figures").mkdir(exist_ok=True)
        (self.output_dir / "data").mkdir(exist_ok=True)
        (self.output_dir / "configs").mkdir(exist_ok=True)

    def run_full_pipeline(self) -> Dict:
        """Execute the complete digital twin pipeline."""
        logger.info(f"{'='*60}")
        logger.info(f"Cardiac Digital Twin Pipeline - Patient {self.patient_id}")
        logger.info(f"{'='*60}")

        # Module 1: Image Processing & Mesh Generation
        mesh_results = self.run_image_processing()

        # Module 2: Electrophysiology Simulation
        ep_results = self.run_electrophysiology()

        # Module 3: Electro-Mechanical Coupling
        mech_results = self.run_mechanics()

        # Module 4: Patient-Specific Parameter Estimation
        inv_results = self.run_inverse_estimation()

        # Module 5: Arrhythmia Risk Assessment
        risk_results = self.run_arrhythmia_assessment()

        # Module 6: AF Ablation Case Study
        ablation_results = self.run_ablation_case_study()

        # Compile final results
        self.results = {
            "patient_id": self.patient_id,
            "timestamp": self.timestamp,
            "modules": {
                "image_processing": mesh_results,
                "electrophysiology": ep_results,
                "mechanics": mech_results,
                "inverse_estimation": inv_results,
                "arrhythmia_risk": risk_results,
                "ablation_prediction": ablation_results,
            }
        }

        # Save consolidated results
        self._save_results()

        return self.results

    def run_image_processing(self) -> Dict:
        """Module 1: Cardiac MRI segmentation and mesh generation."""
        logger.info("\n[Module 1] Image Processing & Mesh Generation")
        logger.info("-" * 50)

        # Create synthetic MRI volume
        volume = MRIVolume(
            data=np.random.default_rng(42).normal(500, 100, (128, 128, 64)).astype(np.float32),
            spacing=(1.25, 1.25, 2.5),
            patient_id=self.patient_id,
            sequence_type="cine_ssfp",
        )

        # Segmentation
        pipeline = CardiacSegmentationPipeline()
        seg_result = pipeline.segment(volume)
        volumes = seg_result.compute_volumes(volume.voxel_volume_mm3)

        # Mesh generation
        mesh_gen = MeshGenerator(target_edge_length=1.0)
        lv_mask = seg_result.get_region_mask(2)  # LV myocardium
        surface_mesh = mesh_gen.generate_surface_mesh(lv_mask, volume.spacing)
        volume_mesh = mesh_gen.generate_volume_mesh(surface_mesh)
        fibers = mesh_gen.assign_fiber_orientation(volume_mesh)

        # Export
        carp_files = export_to_opencarp(
            volume_mesh, fibers, str(self.output_dir / "data" / "opencarp")
        )
        feb_file = export_to_febio(
            volume_mesh, fibers, str(self.output_dir / "data" / "febio")
        )

        results = {
            "mri_shape": list(volume.shape),
            "spacing_mm": list(volume.spacing),
            "segmentation_labels": {k: v for k, v in volumes.items() if v > 0},
            "surface_mesh": {
                "n_vertices": len(surface_mesh["vertices"]),
                "n_faces": len(surface_mesh["faces"]),
                "quality": surface_mesh["quality"],
            },
            "volume_mesh": {
                "n_vertices": len(volume_mesh["vertices"]),
                "n_tetrahedra": len(volume_mesh["tetrahedra"]),
                "quality": volume_mesh["quality"],
            },
            "fiber_angle_range_deg": [-60, 60],
            "exported_files": {
                "opencarp": carp_files,
                "febio": feb_file,
            },
        }

        logger.info(f"  Segmented volumes (mL): {volumes}")
        logger.info(f"  Volume mesh: {results['volume_mesh']['n_tetrahedra']} tetrahedra")
        return results

    def run_electrophysiology(self) -> Dict:
        """Module 2: Cardiac electrophysiology simulation."""
        logger.info("\n[Module 2] Electrophysiology Simulation")
        logger.info("-" * 50)

        results = {}

        # 2a: Aliev-Panfilov (1D cable)
        logger.info("  Running Aliev-Panfilov model (1D cable)...")
        ap_model = AlievPanfilovModel()
        ap_solver = MonodomainSolver(ap_model, TissueParams(duration=200.0))
        ap_solver.setup_1d(n_cells=100, dx=0.2)
        ap_results = ap_solver.solve()

        results["aliev_panfilov"] = {
            "n_cells": 100,
            "duration_ms": 200.0,
            "cv_m_per_s": float(ap_results.get("conduction_velocity", 0)),
            "V_max": float(np.max(ap_results["V"])),
            "V_min": float(np.min(ap_results["V"])),
        }

        # 2b: ten Tusscher (single cell)
        logger.info("  Running ten Tusscher 2006 model (single cell)...")
        tt_model = TenTusscherModel()
        tt_solver = MonodomainSolver(tt_model, TissueParams(dt=0.02, duration=500.0))
        tt_solver.setup_1d(n_cells=1, dx=0.1)
        tt_results = tt_solver.solve([{
            "start": 10.0, "duration": 1.0,
            "amplitude": 52.0, "region": [0]
        }])

        V_trace = tt_results["V"][:, 0]
        results["ten_tusscher"] = {
            "n_states": tt_model.N_STATES,
            "V_rest_mV": float(V_trace[0]) if len(V_trace) > 0 else -86.2,
            "V_peak_mV": float(np.max(V_trace)),
            "duration_ms": 500.0,
        }

        # Generate OpenCARP config
        ep_config_path = generate_opencarp_config(
            TissueParams(), IonicModelType.TEN_TUSSCHER_2006,
            str(self.output_dir / "configs")
        )
        results["config_file"] = ep_config_path

        logger.info(f"  AP CV: {results['aliev_panfilov']['cv_m_per_s']:.3f} m/s")
        logger.info(f"  TT V_peak: {results['ten_tusscher']['V_peak_mV']:.1f} mV")
        return results

    def run_mechanics(self) -> Dict:
        """Module 3: Electro-mechanical coupling simulation."""
        logger.info("\n[Module 3] Cardiac Mechanics & EM Coupling")
        logger.info("-" * 50)

        # Cardiac cycle simulation
        simulator = CardiacCycleSimulator(n_elements=50)
        cycle_results = simulator.simulate_cycle(n_beats=1, bcl=800.0, dt=1.0)

        # Passive mechanics test
        passive = PassiveMechanicsModel()
        F = np.eye(3)
        F[0, 0] = 1.1  # 10% stretch in fiber direction
        f0 = np.array([1.0, 0.0, 0.0])
        s0 = np.array([0.0, 1.0, 0.0])
        W = passive.compute_strain_energy(F, f0, s0)
        S = passive.compute_pk2_stress(F, f0, s0)

        # Generate FEBio config
        feb_config = generate_febio_mechanics_config(
            HolzapfelOgdenParams(), ActiveTensionParams(),
            str(self.output_dir / "configs")
        )

        results = {
            "cardiac_cycle": {
                "EDV_mL": float(cycle_results["EDV"]),
                "ESV_mL": float(cycle_results["ESV"]),
                "SV_mL": float(cycle_results["SV"]),
                "EF_pct": float(cycle_results["EF"]),
                "peak_pressure_kPa": float(cycle_results["peak_pressure"]),
                "peak_pressure_mmHg": float(cycle_results["peak_pressure"] * 7.5),
            },
            "passive_mechanics": {
                "strain_energy_kPa": float(W),
                "pk2_stress_trace_kPa": float(np.trace(S)),
            },
            "constitutive_model": "Holzapfel-Ogden",
            "active_model": "Land et al. 2017",
            "febio_config": feb_config,
        }

        logger.info(f"  EDV={results['cardiac_cycle']['EDV_mL']:.1f} mL, "
                    f"ESV={results['cardiac_cycle']['ESV_mL']:.1f} mL, "
                    f"EF={results['cardiac_cycle']['EF_pct']:.1f}%")
        return results

    def run_inverse_estimation(self) -> Dict:
        """Module 4: Patient-specific parameter estimation."""
        logger.info("\n[Module 4] Inverse Parameter Estimation")
        logger.info("-" * 50)

        # Synthetic ECG data
        rng = np.random.default_rng(42)
        ecg = ECGData(
            signals=rng.normal(0, 0.5, (12, 5000)),
            sampling_rate=500.0,
        )

        # Synthetic echo data
        echo = EchoData(
            edv=130.0, esv=55.0, ef=57.7,
            gls=-19.5, e_prime=9.5, e_a_ratio=1.3,
            wall_thickness={"septal": 10.0, "lateral": 9.0, "anterior": 9.5, "inferior": 9.0},
        )

        # ECG-based inverse
        mesh_nodes = rng.normal(0, 20, (100, 3))
        ecg_solver = ECGInverseSolver(mesh_nodes)
        ecg_result = ecg_solver.estimate_conduction_params(ecg)

        # Echo-based inverse
        mech_solver = MechanicsInverseSolver()
        mech_result = mech_solver.estimate_from_echo(echo)

        results = {
            "ecg_inverse": {
                "estimated_params": ecg_result.estimated_params,
                "converged": ecg_result.converged,
                "final_residual": float(ecg_result.residual),
                "n_iterations": ecg_result.n_iterations,
                "sensitivity": ecg_result.sensitivity,
                "confidence_intervals": {
                    k: [float(v[0]), float(v[1])]
                    for k, v in ecg_result.confidence_intervals.items()
                },
            },
            "mechanics_inverse": {
                "estimated_params": {k: float(v) for k, v in mech_result.estimated_params.items()},
                "converged": mech_result.converged,
                "final_residual": float(mech_result.residual),
                "n_iterations": mech_result.n_iterations,
                "sensitivity": mech_result.sensitivity,
            },
            "clinical_data": {
                "ecg_qrs_ms": float(ecg.get_qrs_duration()),
                "ecg_qt_ms": float(ecg.get_qt_interval()),
                "echo_edv_mL": echo.edv,
                "echo_esv_mL": echo.esv,
                "echo_ef_pct": echo.ef,
                "echo_gls_pct": echo.gls,
            },
        }

        logger.info(f"  ECG inverse converged: {ecg_result.converged}")
        logger.info(f"  Mechanics inverse converged: {mech_result.converged}")
        return results

    def run_arrhythmia_assessment(self) -> Dict:
        """Module 5: Arrhythmia risk assessment."""
        logger.info("\n[Module 5] Arrhythmia Risk Assessment")
        logger.info("-" * 50)

        # Restitution analysis
        reentry = ReentryVulnerabilityAnalysis(cv_long=0.6, cv_trans=0.2, apd_90=280.0)
        restitution = reentry.compute_apd_restitution()
        vw = reentry.compute_vulnerability_window(erp=230.0, substrate_size=30.0)
        s1s2 = reentry.s1s2_protocol()

        # APD dispersion
        rng = np.random.default_rng(42)
        n_elements = 200
        apd_field = 280 + rng.normal(0, 15, n_elements)
        coords = rng.normal(0, 20, (n_elements, 3))

        dispersion_analyzer = APDDispersionAnalyzer()
        dispersion = dispersion_analyzer.compute_dispersion_map(apd_field, coords)

        # Fibrosis
        fibrosis_map = (rng.random(n_elements) < 0.12).astype(int)
        fibrosis_analyzer = FibrosisMappingAnalyzer()
        fibrosis = fibrosis_analyzer.analyze_fibrosis_pattern(fibrosis_map, coords)

        # Comprehensive risk
        assessor = ArrhythmiaRiskAssessor()
        risk = assessor.assess_risk(
            restitution_data=restitution,
            dispersion_data=dispersion,
            fibrosis_data=fibrosis,
            conduction_data=vw,
        )

        results = {
            "restitution": {
                "max_slope": float(restitution["max_slope"]),
                "alternans_prone": bool(restitution["alternans_prone"]),
                "critical_di_ms": float(restitution["critical_di"]),
            },
            "vulnerability": {
                "wavelength_mm": float(vw["wavelength_mm"]),
                "vw_width_ms": float(vw["vw_width_ms"]),
                "reentry_possible": bool(vw["reentry_possible"]),
            },
            "s1s2_protocol": {
                "erp_ms": float(s1s2["erp_ms"]),
                "vw_width_ms": float(s1s2["vw_width_ms"]),
            },
            "apd_dispersion": {
                "mean_apd_ms": float(dispersion["mean_apd_ms"]),
                "std_apd_ms": float(dispersion["std_apd_ms"]),
                "dispersion_index": float(dispersion["dispersion_index"]),
                "n_vulnerable_regions": int(dispersion["n_vulnerable_regions"]),
            },
            "fibrosis": {
                "burden_pct": float(fibrosis["fibrosis_burden_pct"]),
                "pattern": fibrosis["pattern"],
                "border_zone_pct": float(fibrosis["border_zone_fraction"]),
            },
            "overall_risk": {
                "score": float(risk.overall_risk),
                "category": risk.risk_category,
                "reentry_inducible": risk.reentry_inducible,
                "sub_scores": {k: float(v) for k, v in risk.sub_scores.items()},
            },
        }

        logger.info(f"  Risk score: {risk.overall_risk:.3f} ({risk.risk_category})")
        logger.info(f"  Re-entry inducible: {risk.reentry_inducible}")
        return results

    def run_ablation_case_study(self) -> Dict:
        """Module 6: AF ablation case study."""
        logger.info("\n[Module 6] AF Ablation Case Study")
        logger.info("-" * 50)

        case_study = AblationCaseStudy(patient_id=self.patient_id)
        results = case_study.run_case_study()

        logger.info(f"  Optimal strategy: {results['optimal_strategy']}")
        logger.info(f"  Optimal 1-yr recurrence: {results['optimal_recurrence']:.1%}")
        return results

    def _save_results(self):
        """Save all results to JSON."""
        output_file = self.output_dir / "results" / "pipeline_results.json"

        # Convert numpy types
        def convert(obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, np.bool_):
                return bool(obj)
            return obj

        serializable = json.loads(json.dumps(self.results, default=convert))

        with open(output_file, "w") as f:
            json.dump(serializable, f, indent=2, ensure_ascii=False)

        logger.info(f"\nResults saved to {output_file}")
