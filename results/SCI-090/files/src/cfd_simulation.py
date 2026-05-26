"""
CFD Natural Ventilation Simulation
Simplified CFD solver for cross-ventilation analysis using finite difference method.
"""
import numpy as np
from typing import Dict, Tuple, List


class CFDSimulation:
    """2D simplified CFD simulation for cross-ventilation analysis."""

    def __init__(self, params: Dict):
        self.params = params
        self.domain = params['domain']
        self.bc = params['boundary_conditions']
        self.nx = int(self.domain['length'] / self.domain['mesh_resolution'])
        self.ny = int(self.domain['width'] / self.domain['mesh_resolution'])
        self.dx = self.domain['mesh_resolution']
        self.dy = self.domain['mesh_resolution']

        # Flow fields
        self.u = np.zeros((self.ny, self.nx))  # x-velocity
        self.v = np.zeros((self.ny, self.nx))  # y-velocity
        self.p = np.zeros((self.ny, self.nx))  # pressure
        self.T = np.ones((self.ny, self.nx)) * self.bc['temperature_indoor']

        # Air properties
        self.rho = 1.2  # kg/m³
        self.mu = 1.8e-5  # Pa·s
        self.nu = self.mu / self.rho
        self.alpha = 2.2e-5  # thermal diffusivity m²/s

    def _setup_openings(self):
        """Define opening locations in the computational domain."""
        openings = []
        for opening in self.params['openings']:
            wall = opening['wall']
            count = opening['count']
            width_cells = max(1, int(opening['width'] / self.dx))

            if wall == 'south':
                j = 0
                spacing = self.nx // (count + 1)
                for k in range(count):
                    i_start = spacing * (k + 1) - width_cells // 2
                    openings.append(('inlet', j, i_start, i_start + width_cells))
            elif wall == 'north':
                j = self.ny - 1
                spacing = self.nx // (count + 1)
                for k in range(count):
                    i_start = spacing * (k + 1) - width_cells // 2
                    openings.append(('outlet', j, i_start, i_start + width_cells))

        return openings

    def _apply_boundary_conditions(self, openings):
        """Apply boundary conditions."""
        inlet_vel = self.bc['inlet_velocity']
        T_out = self.bc['temperature_outdoor']

        # Walls: no-slip
        self.u[0, :] = 0; self.u[-1, :] = 0
        self.u[:, 0] = 0; self.u[:, -1] = 0
        self.v[0, :] = 0; self.v[-1, :] = 0
        self.v[:, 0] = 0; self.v[:, -1] = 0

        # Apply openings
        for otype, j, i_start, i_end in openings:
            i_start = max(0, min(i_start, self.nx - 1))
            i_end = max(0, min(i_end, self.nx))
            if otype == 'inlet':
                cd = self.params['openings'][0]['discharge_coeff']
                self.v[j, i_start:i_end] = inlet_vel * cd
                self.T[j, i_start:i_end] = T_out
            elif otype == 'outlet':
                self.v[j, i_start:i_end] = self.v[j - 1, i_start:i_end]

    def run_simulation(self, max_iterations: int = 500) -> Dict:
        """Run simplified CFD simulation using iterative method."""
        openings = self._setup_openings()
        dt = 0.01
        residuals = []

        for iteration in range(max_iterations):
            u_old = self.u.copy()
            v_old = self.v.copy()

            # Simplified momentum equation (diffusion dominant)
            for j in range(1, self.ny - 1):
                for i in range(1, self.nx - 1):
                    # Laplacian of velocity
                    lap_u = (u_old[j, i+1] + u_old[j, i-1] + u_old[j+1, i] + u_old[j-1, i] - 4*u_old[j, i]) / (self.dx**2)
                    lap_v = (v_old[j, i+1] + v_old[j, i-1] + v_old[j+1, i] + v_old[j-1, i] - 4*v_old[j, i]) / (self.dx**2)

                    # Convection (upwind)
                    dudx = (u_old[j, i] - u_old[j, i-1]) / self.dx if u_old[j, i] > 0 else (u_old[j, i+1] - u_old[j, i]) / self.dx
                    dvdy = (v_old[j, i] - v_old[j-1, i]) / self.dy if v_old[j, i] > 0 else (v_old[j+1, i] - v_old[j, i]) / self.dy

                    self.u[j, i] = u_old[j, i] + dt * (self.nu * lap_u - u_old[j, i] * dudx)
                    self.v[j, i] = v_old[j, i] + dt * (self.nu * lap_v - v_old[j, i] * dvdy)

            # Temperature transport
            T_old = self.T.copy()
            for j in range(1, self.ny - 1):
                for i in range(1, self.nx - 1):
                    lap_T = (T_old[j, i+1] + T_old[j, i-1] + T_old[j+1, i] + T_old[j-1, i] - 4*T_old[j, i]) / (self.dx**2)
                    adv_T = u_old[j, i] * (T_old[j, i] - T_old[j, i-1]) / self.dx + v_old[j, i] * (T_old[j, i] - T_old[j-1, i]) / self.dy
                    self.T[j, i] = T_old[j, i] + dt * (self.alpha * lap_T - adv_T)

            self._apply_boundary_conditions(openings)

            # Residual
            res = np.sqrt(np.mean((self.u - u_old)**2 + (self.v - v_old)**2))
            residuals.append(res)

            if res < self.params['convergence_criterion'] and iteration > 50:
                break

        # Post-process
        velocity_magnitude = np.sqrt(self.u**2 + self.v**2)
        avg_velocity = np.mean(velocity_magnitude[1:-1, 1:-1])
        max_velocity = np.max(velocity_magnitude)
        avg_temp = np.mean(self.T[1:-1, 1:-1])

        # Ventilation effectiveness
        T_supply = self.bc['temperature_outdoor']
        T_exhaust = np.mean(self.T[-2, :])
        T_room = avg_temp
        if abs(T_exhaust - T_supply) > 0.01:
            vent_effectiveness = (T_exhaust - T_supply) / (T_room - T_supply) if abs(T_room - T_supply) > 0.01 else 1.0
        else:
            vent_effectiveness = 1.0

        # Air change rate estimation
        total_opening_area = sum(
            o['width'] * o['height'] * o['count'] * o['discharge_coeff']
            for o in self.params['openings']
        )
        volume = self.domain['length'] * self.domain['width'] * self.domain['height']
        ach = (avg_velocity * total_opening_area * 3600) / volume

        # Comfort assessment (ASHRAE 55 adaptive)
        comfort_velocity_range = (0.15, 0.8)  # m/s
        comfort_zone_fraction = np.mean(
            (velocity_magnitude[1:-1, 1:-1] >= comfort_velocity_range[0]) &
            (velocity_magnitude[1:-1, 1:-1] <= comfort_velocity_range[1])
        )

        return {
            'iterations': iteration + 1,
            'converged': res < self.params['convergence_criterion'],
            'final_residual': float(res),
            'velocity_field': velocity_magnitude.tolist(),
            'temperature_field': self.T.tolist(),
            'avg_indoor_velocity_ms': round(float(avg_velocity), 3),
            'max_velocity_ms': round(float(max_velocity), 3),
            'avg_indoor_temp_C': round(float(avg_temp), 2),
            'ventilation_effectiveness': round(float(vent_effectiveness), 3),
            'air_changes_per_hour': round(float(ach), 1),
            'comfort_zone_fraction': round(float(comfort_zone_fraction), 3),
            'residuals': residuals[-50:],  # last 50 residuals
            'domain_size': f"{self.domain['length']}m x {self.domain['width']}m",
            'mesh_cells': f"{self.nx} x {self.ny}",
        }

    def evaluate_cross_ventilation_scenarios(self) -> List[Dict]:
        """Evaluate multiple cross-ventilation scenarios."""
        scenarios = []
        wind_speeds = [1.0, 2.0, 3.0, 4.0, 5.0]
        wind_directions = [0, 45, 90, 135, 180]

        for speed in wind_speeds:
            self.bc['inlet_velocity'] = speed
            # Reset fields
            self.u = np.zeros((self.ny, self.nx))
            self.v = np.zeros((self.ny, self.nx))
            self.T = np.ones((self.ny, self.nx)) * self.bc['temperature_indoor']

            result = self.run_simulation(max_iterations=300)
            scenarios.append({
                'wind_speed': speed,
                'wind_direction': 180,
                'avg_velocity': result['avg_indoor_velocity_ms'],
                'ach': result['air_changes_per_hour'],
                'avg_temp': result['avg_indoor_temp_C'],
                'comfort_fraction': result['comfort_zone_fraction'],
                'vent_effectiveness': result['ventilation_effectiveness'],
            })

        return scenarios


if __name__ == "__main__":
    from ifc_converter import IFCConverter
    import json

    converter = IFCConverter()
    model = converter.create_reference_building()
    cfd_params = converter.generate_cfd_params()

    sim = CFDSimulation(cfd_params)
    result = sim.run_simulation()
    print("=== CFD Simulation Results ===")
    for k, v in result.items():
        if k not in ['velocity_field', 'temperature_field', 'residuals']:
            print(f"  {k}: {v}")
