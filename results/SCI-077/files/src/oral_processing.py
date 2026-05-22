"""
Module 4: 口腔内プロセシング（咀嚼、嚥下）のシミュレーション
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Optional


@dataclass
class FoodParticle:
    size: float
    mass: float
    moisture: float
    is_broken: bool = False


@dataclass
class JawDynamicsModel:
    max_bite_force: float = 500.0
    jaw_frequency: float = 1.5
    opening_angle: float = 25.0
    molar_position: float = 60.0

    def jaw_trajectory(self, t):
        omega = 2 * np.pi * self.jaw_frequency
        angle = self.opening_angle * np.abs(np.sin(omega * t / 2))
        y = self.molar_position * np.sin(np.radians(angle))
        x = 2.0 * np.cos(omega * t)
        return {'time': t, 'angle': angle, 'y_displacement': y, 'x_displacement': x}

    def bite_force_profile(self, t, food_height=10.0):
        omega = 2 * np.pi * self.jaw_frequency
        angle = self.opening_angle * np.abs(np.sin(omega * t / 2))
        gap = self.molar_position * np.sin(np.radians(angle))
        force = np.zeros_like(t)
        contact = gap < food_height
        if np.any(contact):
            penetration = (food_height - gap[contact]) / food_height
            v = np.gradient(gap[contact], t[contact])
            v_max = 100.0
            hill_factor = (v_max + np.abs(v)) / (v_max + 2*np.abs(v) + 1e-6)
            force[contact] = self.max_bite_force * penetration * hill_factor
        return force


@dataclass
class MasticationSimulator:
    jaw: JawDynamicsModel = field(default_factory=JawDynamicsModel)
    saliva_flow_rate: float = 0.5
    food_hardness: float = 100.0
    food_cohesiveness: float = 0.5
    critical_stress: float = 50.0

    def selection_function(self, particle_size, gap):
        if particle_size <= gap:
            return 0.0
        return 1 - np.exp(-2.0 * (particle_size / gap)**1.5)

    def breakage_function(self, parent_size, n_fragments=3):
        rng = np.random.default_rng()
        mean_size = parent_size / n_fragments**(1/3)
        fragments = rng.lognormal(np.log(mean_size), 0.3, n_fragments)
        total_vol = (parent_size/2)**3
        frag_vol = np.sum((fragments/2)**3)
        fragments *= (total_vol / frag_vol)**(1/3)
        return fragments.tolist()

    def simulate_mastication(self, initial_size=15.0, initial_mass=5.0, n_chews=30):
        particles = [FoodParticle(initial_size, initial_mass, 30.0)]
        results = {k: [] for k in ['chew_number','n_particles','mean_size',
                                    'median_size','d90','moisture','bolus_cohesion','swallowable']}
        for chew in range(n_chews):
            t_cycle = 1.0 / self.jaw.jaw_frequency
            t = np.linspace(0, t_cycle, 100)
            forces = self.jaw.bite_force_profile(t, food_height=initial_size)
            max_force = np.max(forces)
            jaw_traj = self.jaw.jaw_trajectory(t)
            min_gap = max(np.min(jaw_traj['y_displacement']), 0.5)
            new_particles = []
            for p in particles:
                S = self.selection_function(p.size, min_gap)
                if np.random.random() < S and max_force > self.critical_stress*0.1:
                    for fs in self.breakage_function(p.size):
                        new_particles.append(FoodParticle(fs, p.mass*(fs/p.size)**3, p.moisture))
                else:
                    new_particles.append(p)
            particles = new_particles
            saliva = self.saliva_flow_rate * t_cycle / 60
            total_mass = sum(p.mass for p in particles)
            for p in particles:
                p.moisture = min(p.moisture + saliva/total_mass*100*0.5, 90.0)
            sizes = np.array([p.size for p in particles])
            mean_moist = np.mean([p.moisture for p in particles])
            size_uni = 1 - np.std(sizes)/(np.mean(sizes)+1e-6)
            cohesion = np.clip((mean_moist/100)*max(size_uni,0)*self.food_cohesiveness, 0, 1)
            d90 = np.percentile(sizes, 90)
            swallowable = bool(d90 < 3.0 and mean_moist > 50.0 and cohesion > 0.15)
            results['chew_number'].append(chew+1)
            results['n_particles'].append(len(particles))
            results['mean_size'].append(np.mean(sizes))
            results['median_size'].append(np.median(sizes))
            results['d90'].append(d90)
            results['moisture'].append(mean_moist)
            results['bolus_cohesion'].append(cohesion)
            results['swallowable'].append(swallowable)
        return results


@dataclass
class SwallowingModel:
    pharynx_length: float = 120.0
    pharynx_diameter: float = 20.0
    peristaltic_pressure: float = 30.0

    def transport_velocity(self, viscosity):
        R = self.pharynx_diameter / 2 * 1e-3
        L = self.pharynx_length * 1e-3
        dP = self.peristaltic_pressure * 1e3
        return dP * R**2 / (8 * viscosity * L)

    def transit_time(self, viscosity):
        v = self.transport_velocity(viscosity)
        return (self.pharynx_length * 1e-3) / max(v, 1e-6)

    def residue_fraction(self, viscosity, yield_stress=0.0):
        if yield_stress <= 0:
            return 0.02
        tau_w = self.peristaltic_pressure * self.pharynx_diameter / (4*self.pharynx_length)
        if tau_w <= yield_stress:
            return 1.0
        return np.clip((yield_stress / tau_w)**2, 0.02, 1.0)

    def simulate_swallowing(self, bolus_volume=10.0, viscosity=0.1, yield_stress=0.0, n_steps=50):
        v = self.transport_velocity(viscosity)
        t_transit = self.transit_time(viscosity)
        residue = self.residue_fraction(viscosity, yield_stress)
        t = np.linspace(0, t_transit*1.5, n_steps)
        position = np.minimum(v*t*1e3, self.pharynx_length)
        vol_transported = np.minimum(bolus_volume*(1-residue)*t/t_transit, bolus_volume*(1-residue))
        return {'time': t, 'position': position, 'volume_transported': vol_transported,
                'transit_time': t_transit, 'transport_velocity': v,
                'residue_fraction': residue, 'residue_volume': bolus_volume*residue}


def simulate_oral_processing(food_hardness=100.0, food_cohesiveness=0.5,
                              bolus_viscosity=0.1, n_chews=30):
    masticator = MasticationSimulator(food_hardness=food_hardness,
                                      food_cohesiveness=food_cohesiveness)
    mast_result = masticator.simulate_mastication(n_chews=n_chews)
    swallow_chew = None
    for i, sw in enumerate(mast_result['swallowable']):
        if sw:
            swallow_chew = i + 1
            break
    swallower = SwallowingModel()
    swal_result = swallower.simulate_swallowing(viscosity=bolus_viscosity)
    return {'mastication': mast_result, 'swallowing': swal_result,
            'swallow_trigger_chew': swallow_chew if swallow_chew else n_chews}
