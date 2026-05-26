"""
CFD-based Natural Ventilation and Cross-Ventilation Analysis
Simplified 2D CFD solver using finite difference method for indoor airflow.
"""
import numpy as np

class VentilationCFD:
    """2D steady-state CFD solver for cross-ventilation analysis."""
    
    def __init__(self, room_width=10.0, room_depth=8.0, resolution=50):
        self.W = room_width
        self.D = room_depth
        self.nx = resolution
        self.ny = int(resolution * room_depth / room_width)
        self.dx = room_width / self.nx
        self.dy = room_depth / self.ny
        
        self.u = np.zeros((self.ny, self.nx))  # x-velocity
        self.v = np.zeros((self.ny, self.nx))  # y-velocity
        self.p = np.zeros((self.ny, self.nx))  # pressure
        
    def setup_cross_ventilation(self, wind_speed=3.0, wind_angle=0, 
                                 inlet_pos=0.4, inlet_size=0.2,
                                 outlet_pos=0.4, outlet_size=0.2):
        """Configure cross-ventilation with inlet/outlet openings."""
        self.wind_speed = wind_speed
        self.wind_angle = wind_angle
        
        # Inlet on south wall
        inlet_start = int(inlet_pos * self.nx)
        inlet_end = int((inlet_pos + inlet_size) * self.nx)
        
        # Outlet on north wall
        outlet_start = int(outlet_pos * self.nx)
        outlet_end = int((outlet_pos + outlet_size) * self.nx)
        
        self.inlet = (inlet_start, inlet_end)
        self.outlet = (outlet_start, outlet_end)
        
        # Pressure coefficient
        Cp_inlet = 0.6
        Cp_outlet = -0.3
        
        # Effective velocity through openings
        Cd = 0.65
        delta_Cp = Cp_inlet - Cp_outlet
        v_eff = Cd * wind_speed * np.sqrt(abs(delta_Cp))
        
        self.v_inlet = v_eff
        self.delta_Cp = delta_Cp
        
        return v_eff
    
    def solve(self, iterations=500):
        """Iterative solver for velocity field (simplified potential flow)."""
        u = self.u.copy()
        v = self.v.copy()
        
        # Set inlet BC
        u[0, self.inlet[0]:self.inlet[1]] = 0
        v[0, self.inlet[0]:self.inlet[1]] = self.v_inlet
        
        # Set outlet BC
        v[-1, self.outlet[0]:self.outlet[1]] = self.v_inlet * 0.8
        
        for iteration in range(iterations):
            u_new = u.copy()
            v_new = v.copy()
            
            # Interior points - simplified diffusion
            for i in range(1, self.ny - 1):
                for j in range(1, self.nx - 1):
                    u_new[i, j] = 0.25 * (u[i+1, j] + u[i-1, j] + u[i, j+1] + u[i, j-1])
                    v_new[i, j] = 0.25 * (v[i+1, j] + v[i-1, j] + v[i, j+1] + v[i, j-1])
            
            # Wall BCs (no-slip)
            u_new[0, :] = 0
            u_new[-1, :] = 0
            u_new[:, 0] = 0
            u_new[:, -1] = 0
            v_new[:, 0] = 0
            v_new[:, -1] = 0
            
            # Inlet/outlet BCs
            v_new[0, self.inlet[0]:self.inlet[1]] = self.v_inlet
            v_new[-1, self.outlet[0]:self.outlet[1]] = self.v_inlet * 0.8
            
            u = u_new
            v = v_new
        
        self.u = u
        self.v = v
        self.velocity_mag = np.sqrt(u**2 + v**2)
        
        return self.velocity_mag
    
    def calculate_metrics(self):
        """Calculate ventilation performance metrics."""
        # Average indoor air velocity
        interior = self.velocity_mag[2:-2, 2:-2]
        avg_velocity = np.mean(interior)
        max_velocity = np.max(interior)
        
        # Air change rate estimation
        inlet_area = (self.inlet[1] - self.inlet[0]) * self.dx * 2.5  # assume 2.5m height
        Q = self.v_inlet * inlet_area
        room_volume = self.W * self.D * 3.0
        ACH = Q * 3600 / room_volume
        
        # Ventilation effectiveness
        # Based on age-of-air concept
        vel_std = np.std(interior)
        uniformity = 1 - vel_std / (avg_velocity + 1e-6)
        
        # Comfort zones (0.15 < v < 0.8 m/s per ASHRAE 55)
        comfort_zone = np.logical_and(interior > 0.15, interior < 0.80)
        comfort_ratio = np.sum(comfort_zone) / interior.size
        
        return {
            "avg_velocity_ms": avg_velocity,
            "max_velocity_ms": max_velocity,
            "ACH": ACH,
            "airflow_rate_m3s": Q,
            "uniformity_index": uniformity,
            "comfort_ratio": comfort_ratio,
            "Cd_effective": 0.65,
            "delta_Cp": self.delta_Cp,
        }
    
    def parametric_study(self, wind_speeds=[1, 2, 3, 4, 5],
                         opening_sizes=[0.10, 0.15, 0.20, 0.25, 0.30]):
        """Run parametric study on wind speed and opening size."""
        results = []
        for ws in wind_speeds:
            for os in opening_sizes:
                self.setup_cross_ventilation(wind_speed=ws, inlet_size=os, outlet_size=os)
                self.solve(iterations=300)
                metrics = self.calculate_metrics()
                metrics["wind_speed"] = ws
                metrics["opening_ratio"] = os
                results.append(metrics)
        return results


def run_multi_scenario():
    """Run multiple ventilation scenarios for the case study."""
    scenarios = [
        {"name": "Baseline", "wind": 3.0, "inlet_pos": 0.4, "inlet_size": 0.20, "outlet_pos": 0.4, "outlet_size": 0.20},
        {"name": "Large_Opening", "wind": 3.0, "inlet_pos": 0.3, "inlet_size": 0.40, "outlet_pos": 0.3, "outlet_size": 0.40},
        {"name": "Offset_Opening", "wind": 3.0, "inlet_pos": 0.2, "inlet_size": 0.20, "outlet_pos": 0.6, "outlet_size": 0.20},
        {"name": "Low_Wind", "wind": 1.5, "inlet_pos": 0.4, "inlet_size": 0.20, "outlet_pos": 0.4, "outlet_size": 0.20},
        {"name": "High_Wind", "wind": 5.0, "inlet_pos": 0.4, "inlet_size": 0.20, "outlet_pos": 0.4, "outlet_size": 0.20},
    ]
    
    all_results = {}
    for sc in scenarios:
        cfd = VentilationCFD()
        cfd.setup_cross_ventilation(
            wind_speed=sc["wind"],
            inlet_pos=sc["inlet_pos"],
            inlet_size=sc["inlet_size"],
            outlet_pos=sc["outlet_pos"],
            outlet_size=sc["outlet_size"],
        )
        vel = cfd.solve(iterations=400)
        metrics = cfd.calculate_metrics()
        all_results[sc["name"]] = {
            "metrics": metrics,
            "velocity_field": vel,
        }
    
    return all_results


if __name__ == "__main__":
    results = run_multi_scenario()
    for name, data in results.items():
        m = data["metrics"]
        print(f"{name}: ACH={m['ACH']:.1f}, Avg_v={m['avg_velocity_ms']:.3f} m/s, Comfort={m['comfort_ratio']:.1%}")
