"""
Module 6: Automotive Parts Manufacturing Case Study
"""

import numpy as np
from typing import Dict


class AutomotiveCaseStudy:
    def __init__(self):
        self.part_spec = {
            'name': 'Bumper Bracket',
            'material': 'PA66-GF30',
            'dimensions_mm': [200, 100, 3],
            'weight_target_g': 45.0,
            'weight_tolerance_g': 0.5,
            'warpage_limit_mm': 0.5,
            'sink_mark_limit_mm': 0.02,
            'shrinkage_limit_pct': 1.0,
        }
        self.production_scenarios = {
            'nominal': dict(injection_pressure_MPa=80, packing_pressure_MPa=50,
                          cooling_time_s=20, melt_temp_C=280, mold_temp_C=80, injection_speed_mm_s=50),
            'high_speed': dict(injection_pressure_MPa=100, packing_pressure_MPa=60,
                             cooling_time_s=15, melt_temp_C=290, mold_temp_C=70, injection_speed_mm_s=75),
            'low_stress': dict(injection_pressure_MPa=65, packing_pressure_MPa=45,
                             cooling_time_s=30, melt_temp_C=275, mold_temp_C=90, injection_speed_mm_s=35),
            'optimized': dict(injection_pressure_MPa=85, packing_pressure_MPa=55,
                            cooling_time_s=22, melt_temp_C=282, mold_temp_C=85, injection_speed_mm_s=48),
        }

    def simulate_quality(self, params: Dict) -> Dict:
        P_inj = params['injection_pressure_MPa']
        P_pack = params['packing_pressure_MPa']
        t_cool = params['cooling_time_s']
        T_melt = params['melt_temp_C']
        T_mold = params['mold_temp_C']
        v_inj = params['injection_speed_mm_s']
        dT = T_melt - T_mold
        fill_time = 200 / v_inj
        max_shear_rate = 6 * v_inj * 0.1 * 0.003 / (0.1 * 0.003 ** 2)
        alpha_thermal = 0.30 / (1350 * 1700)
        fourier = alpha_thermal * t_cool / (0.003 / 2) ** 2
        T_center_final = T_mold + (T_melt - T_mold) * np.exp(-fourier * 0.5)
        crystallinity = 0.25 + 0.10 * (1 - np.exp(-t_cool / 15)) * (T_mold / 100)
        warpage = 0.15 * (dT / 200) ** 1.3 * (1 - P_pack / 100) * \
                  np.exp(-t_cool / 30) * (1 + 0.05 * (90 - T_mold) / 30)
        sink = 0.015 * (1 - P_pack / P_inj) * np.exp(-t_cool / 25) * (T_mold / 100)
        pvT_shrinkage = 0.005 * (T_melt / 280) * (80 / max(P_pack, 30))
        weight = 45.0 * (1 - pvT_shrinkage)
        shrinkage = 0.8 * (1 - P_pack / 100) * (dT / 200) * np.exp(-t_cool / 35) + 0.3
        residual_stress = 5.0 * (P_inj / 80) * (dT / 200) * (1 - t_cool / 50) + 3.0 * (v_inj / 50)
        cycle_time = fill_time + 8.0 + t_cool + 3.0
        quality_pass = (
            abs(warpage) < self.part_spec['warpage_limit_mm'] and
            abs(sink) < self.part_spec['sink_mark_limit_mm'] and
            abs(weight - self.part_spec['weight_target_g']) < self.part_spec['weight_tolerance_g'] and
            shrinkage < self.part_spec['shrinkage_limit_pct']
        )
        return {
            'quality_metrics': {
                'warpage_mm': round(float(abs(warpage)), 4),
                'sink_depth_mm': round(float(abs(sink)), 4),
                'weight_g': round(float(weight), 2),
                'shrinkage_pct': round(float(shrinkage), 3),
                'residual_stress_MPa': round(float(abs(residual_stress)), 2),
                'crystallinity_pct': round(float(crystallinity * 100), 1),
            },
            'process_metrics': {
                'fill_time_s': round(float(fill_time), 2),
                'cycle_time_s': round(float(cycle_time), 1),
                'max_shear_rate_1_s': round(float(max_shear_rate), 0),
                'center_temp_at_eject_C': round(float(T_center_final), 1),
            },
            'quality_pass': quality_pass,
            'defect_flags': {
                'warpage_fail': abs(warpage) >= self.part_spec['warpage_limit_mm'],
                'sink_fail': abs(sink) >= self.part_spec['sink_mark_limit_mm'],
                'weight_fail': abs(weight - self.part_spec['weight_target_g']) >= self.part_spec['weight_tolerance_g'],
                'shrinkage_fail': shrinkage >= self.part_spec['shrinkage_limit_pct'],
            }
        }

    def run_all_scenarios(self) -> Dict:
        return {name: self.simulate_quality(params) for name, params in self.production_scenarios.items()}

    def run_monte_carlo(self, scenario: str = 'nominal', n_samples: int = 500) -> Dict:
        base_params = self.production_scenarios[scenario].copy()
        np.random.seed(42)
        mc_results = {'warpage': [], 'sink': [], 'weight': [], 'shrinkage': [], 'stress': [], 'pass': []}
        for _ in range(n_samples):
            noisy_params = {key: val * (1.0 + 0.05 * np.random.randn()) for key, val in base_params.items()}
            result = self.simulate_quality(noisy_params)
            mc_results['warpage'].append(result['quality_metrics']['warpage_mm'])
            mc_results['sink'].append(result['quality_metrics']['sink_depth_mm'])
            mc_results['weight'].append(result['quality_metrics']['weight_g'])
            mc_results['shrinkage'].append(result['quality_metrics']['shrinkage_pct'])
            mc_results['stress'].append(result['quality_metrics']['residual_stress_MPa'])
            mc_results['pass'].append(result['quality_pass'])
        summary = {}
        for key in ['warpage', 'sink', 'weight', 'shrinkage', 'stress']:
            data = np.array(mc_results[key])
            summary[key] = {
                'mean': float(np.mean(data)), 'std': float(np.std(data)),
                'min': float(np.min(data)), 'max': float(np.max(data)),
                'p5': float(np.percentile(data, 5)), 'p95': float(np.percentile(data, 95)),
            }
        return {
            'n_samples': n_samples, 'scenario': scenario,
            'pass_rate_pct': float(sum(mc_results['pass']) / n_samples * 100),
            'statistics': summary, 'raw_data': mc_results,
        }
