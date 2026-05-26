"""
Mitigation Strategy Module
Quantifies cooling effects of green infrastructure and high-albedo materials.
"""
import numpy as np


class MitigationScenario:
    """Defines and evaluates UHI mitigation strategies."""

    def __init__(self, name, grid_size=(50, 50)):
        self.name = name
        self.nx, self.ny = grid_size
        self.green_fraction_delta = np.zeros((self.nx, self.ny))
        self.albedo_roof_delta = np.zeros((self.nx, self.ny))
        self.albedo_road_delta = np.zeros((self.nx, self.ny))
        self.cool_roof_fraction = np.zeros((self.nx, self.ny))
        self.tree_cover_delta = np.zeros((self.nx, self.ny))

    @staticmethod
    def create_baseline(grid_size=(50, 50)):
        return MitigationScenario("Baseline", grid_size)

    @staticmethod
    def create_green_scenario(grid_size=(50, 50)):
        """Scenario: 30% increase in urban greenery."""
        s = MitigationScenario("Green Infrastructure", grid_size)
        cx, cy = grid_size[0] // 2, grid_size[1] // 2
        for i in range(grid_size[0]):
            for j in range(grid_size[1]):
                dist = np.sqrt((i - cx)**2 + (j - cy)**2)
                if dist < 10:
                    s.green_fraction_delta[i, j] = 0.15
                    s.tree_cover_delta[i, j] = 0.10
                elif dist < 20:
                    s.green_fraction_delta[i, j] = 0.20
                    s.tree_cover_delta[i, j] = 0.15
                else:
                    s.green_fraction_delta[i, j] = 0.10
                    s.tree_cover_delta[i, j] = 0.05
        return s

    @staticmethod
    def create_cool_roof_scenario(grid_size=(50, 50)):
        """Scenario: Cool roofs with albedo 0.6 on 70% of buildings."""
        s = MitigationScenario("Cool Roofs", grid_size)
        s.cool_roof_fraction = np.full(grid_size, 0.70)
        s.albedo_roof_delta = np.full(grid_size, 0.35)  # from 0.20 to 0.55
        return s

    @staticmethod
    def create_combined_scenario(grid_size=(50, 50)):
        """Combined green + cool roof scenario."""
        s = MitigationScenario("Combined", grid_size)
        green = MitigationScenario.create_green_scenario(grid_size)
        cool = MitigationScenario.create_cool_roof_scenario(grid_size)
        s.green_fraction_delta = green.green_fraction_delta
        s.tree_cover_delta = green.tree_cover_delta
        s.cool_roof_fraction = cool.cool_roof_fraction
        s.albedo_roof_delta = cool.albedo_roof_delta
        return s

    def apply_to_ucm(self, ucm):
        """Apply mitigation measures to UCM parameters."""
        ucm.green_fraction = np.clip(ucm.green_fraction + self.green_fraction_delta, 0, 0.80)
        ucm.albedo_roof = np.clip(ucm.albedo_roof + self.albedo_roof_delta, 0.10, 0.85)
        ucm.albedo_road = np.clip(ucm.albedo_road + self.albedo_road_delta, 0.05, 0.60)

    def estimate_cooling(self, T_baseline, T_mitigated):
        """Compute cooling effectiveness metrics."""
        delta_T = T_baseline - T_mitigated
        return {
            'mean_cooling': np.mean(delta_T),
            'max_cooling': np.max(delta_T),
            'median_cooling': np.median(delta_T),
            'cooling_field': delta_T,
            'pct_area_cooled_1K': np.mean(delta_T > 1.0) * 100,
            'pct_area_cooled_2K': np.mean(delta_T > 2.0) * 100,
        }
