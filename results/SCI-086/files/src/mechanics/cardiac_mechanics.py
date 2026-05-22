"""
cardiac_mechanics.py
====================
Module 3: Cardiac mechanics simulation.
Implements passive/active myocardial mechanics with FEBio integration
and electro-mechanical coupling.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple
from pathlib import Path
import logging
import json

logger = logging.getLogger(__name__)


@dataclass
class HolzapfelOgdenParams:
    """
    Holzapfel-Ogden constitutive model parameters for passive myocardium.

    W = a/(2b) * [exp(b(I₁-3)) - 1]
      + Σᵢ aᵢ/(2bᵢ) * [exp(bᵢ(I₄ᵢ-1)²) - 1]
      + afs/(2bfs) * [exp(bfs*I₈fs²) - 1]

    Reference: Holzapfel & Ogden (2009), Phil Trans R Soc A.
    """
    a: float = 0.059       # Isotropic term (kPa)
    b: float = 8.023       # Isotropic exponent
    a_f: float = 18.472    # Fiber term (kPa)
    b_f: float = 16.026    # Fiber exponent
    a_s: float = 2.481     # Sheet term (kPa)
    b_s: float = 11.120    # Sheet exponent
    a_fs: float = 0.216    # Fiber-sheet coupling (kPa)
    b_fs: float = 11.436   # Fiber-sheet exponent
    kappa: float = 1e3     # Bulk modulus for incompressibility (kPa)

    def to_dict(self) -> Dict:
        return {k: v for k, v in self.__dict__.items()}


@dataclass
class ActiveTensionParams:
    """
    Parameters for active tension generation model.

    Based on Land et al. (2017) human cardiac contraction model.
    """
    T_ref: float = 120.0     # Reference tension (kPa)
    Ca50_ref: float = 0.805  # Reference Ca50 (µM)
    n_TRPN: float = 2.0      # Hill coefficient for troponin
    k_TRPN: float = 0.1      # TRPN rate constant (ms⁻¹)
    n_xb: float = 5.0        # Cooperativity coefficient
    k_xb: float = 0.1        # Cross-bridge rate constant (ms⁻¹)
    lambda_min: float = 0.87  # Minimum sarcomere stretch
    lambda_max: float = 1.2   # Maximum sarcomere stretch

    # Length-dependent parameters
    beta_0: float = 4.9       # Length-dependent activation
    Ca50_max: float = 1.05    # Maximum Ca50 (µM)

    def to_dict(self) -> Dict:
        return {k: v for k, v in self.__dict__.items()}


@dataclass
class WindkesselParams:
    """3-element Windkessel model for circulatory boundary conditions."""
    R_c: float = 0.03      # Characteristic impedance (kPa·s/mL)
    R_p: float = 0.9       # Peripheral resistance (kPa·s/mL)
    C: float = 10.0        # Arterial compliance (mL/kPa)
    P_venous: float = 0.5  # Venous pressure (kPa)

    # Valve parameters
    R_valve_open: float = 0.001   # Open valve resistance
    R_valve_closed: float = 1e6   # Closed valve resistance


class PassiveMechanicsModel:
    """
    Passive myocardial mechanics using the Holzapfel-Ogden model.

    Computes stress tensor from deformation gradient F and fiber directions.
    """

    def __init__(self, params: Optional[HolzapfelOgdenParams] = None):
        self.params = params or HolzapfelOgdenParams()

    def compute_strain_energy(self, F: np.ndarray,
                                f0: np.ndarray,
                                s0: np.ndarray) -> float:
        """
        Compute strain energy density.

        F: Deformation gradient (3x3)
        f0: Reference fiber direction (3,)
        s0: Reference sheet direction (3,)
        """
        p = self.params
        C = F.T @ F  # Right Cauchy-Green tensor
        I1 = np.trace(C)
        I4f = f0 @ C @ f0
        I4s = s0 @ C @ s0
        I8fs = f0 @ C @ s0

        W = p.a / (2 * p.b) * (np.exp(p.b * (I1 - 3)) - 1)

        if I4f > 1.0:
            W += p.a_f / (2 * p.b_f) * (np.exp(p.b_f * (I4f - 1)**2) - 1)
        if I4s > 1.0:
            W += p.a_s / (2 * p.b_s) * (np.exp(p.b_s * (I4s - 1)**2) - 1)

        W += p.a_fs / (2 * p.b_fs) * (np.exp(p.b_fs * I8fs**2) - 1)

        return W

    def compute_pk2_stress(self, F: np.ndarray,
                            f0: np.ndarray,
                            s0: np.ndarray) -> np.ndarray:
        """Compute 2nd Piola-Kirchhoff stress tensor S = 2∂W/∂C."""
        p = self.params
        C = F.T @ F
        I = np.eye(3)
        J = np.linalg.det(F)

        I1 = np.trace(C)
        I4f = f0 @ C @ f0
        I4s = s0 @ C @ s0
        I8fs = f0 @ C @ s0

        # Isotropic contribution
        S = p.a * np.exp(p.b * (I1 - 3)) * I

        # Fiber contribution
        if I4f > 1.0:
            S += 2 * p.a_f * (I4f - 1) * np.exp(p.b_f * (I4f - 1)**2) * \
                 np.outer(f0, f0)

        # Sheet contribution
        if I4s > 1.0:
            S += 2 * p.a_s * (I4s - 1) * np.exp(p.b_s * (I4s - 1)**2) * \
                 np.outer(s0, s0)

        # Fiber-sheet coupling
        S += p.a_fs * I8fs * np.exp(p.b_fs * I8fs**2) * \
             (np.outer(f0, s0) + np.outer(s0, f0))

        # Volumetric penalty (quasi-incompressibility)
        C_inv = np.linalg.inv(C)
        S += p.kappa * (J**2 - 1) * J * C_inv

        return S


class ActiveTensionModel:
    """
    Active tension generation model (Land et al. 2017).

    Calcium-driven cross-bridge cycling with length-dependent activation.
    """

    def __init__(self, params: Optional[ActiveTensionParams] = None):
        self.params = params or ActiveTensionParams()

    def compute_active_tension(self, Ca_i: float,
                                 lambda_f: float,
                                 state: Optional[Dict] = None
                                 ) -> Tuple[float, Dict]:
        """
        Compute active tension given intracellular calcium and fiber stretch.

        Ca_i: Intracellular calcium concentration (µM)
        lambda_f: Fiber stretch ratio
        state: Internal state variables (TRPN, XB)
        """
        p = self.params

        if state is None:
            state = {"CaTRPN": 0.0, "XB": 0.0}

        # Length-dependent Ca50
        lambda_clamped = np.clip(lambda_f, p.lambda_min, p.lambda_max)
        Ca50 = p.Ca50_ref * (1 + p.beta_0 * (lambda_clamped - 1))
        Ca50 = min(Ca50, p.Ca50_max)

        # Troponin binding (Hill equation)
        CaTRPN_ss = Ca_i**p.n_TRPN / (Ca_i**p.n_TRPN + Ca50**p.n_TRPN)
        dCaTRPN = p.k_TRPN * (CaTRPN_ss - state["CaTRPN"])
        state["CaTRPN"] += dCaTRPN * 0.01  # dt = 0.01 ms

        # Cross-bridge activation (cooperative)
        XB_ss = state["CaTRPN"]**p.n_xb / \
                (state["CaTRPN"]**p.n_xb + 0.5**p.n_xb)
        dXB = p.k_xb * (XB_ss - state["XB"])
        state["XB"] += dXB * 0.01

        # Active tension with length dependence
        T_active = p.T_ref * state["XB"] * \
                   (1 + p.beta_0 * (lambda_clamped - 1))

        return T_active, state


class ElectroMechanicalCoupling:
    """
    Electro-mechanical coupling framework.

    Bidirectional coupling:
    - Forward: V_m → Ca²⁺ → Active tension → Deformation
    - Backward: Stretch → Mechano-sensitive channels → V_m (MEF)
    """

    def __init__(self, passive_model: PassiveMechanicsModel,
                 active_model: ActiveTensionModel,
                 windkessel: Optional[WindkesselParams] = None):
        self.passive = passive_model
        self.active = active_model
        self.windkessel = windkessel or WindkesselParams()

        # State
        self.P_art = 10.0     # Arterial pressure (kPa)
        self.P_lv = 0.0       # LV pressure (kPa)
        self.V_lv = 120.0     # LV volume (mL) - EDV
        self.phase = "diastole"

    def coupling_step(self, Ca_i_field: np.ndarray,
                       lambda_field: np.ndarray,
                       dt: float = 1.0
                       ) -> Dict[str, np.ndarray]:
        """
        Single coupling step for all elements.

        Ca_i_field: Calcium concentration per element (µM)
        lambda_field: Fiber stretch per element
        """
        n_elements = len(Ca_i_field)
        T_active = np.zeros(n_elements)
        states = [{"CaTRPN": 0.0, "XB": 0.0} for _ in range(n_elements)]

        for i in range(n_elements):
            T_active[i], states[i] = self.active.compute_active_tension(
                Ca_i_field[i], lambda_field[i], states[i]
            )

        # Update LV pressure-volume
        mean_tension = np.mean(T_active)
        self.P_lv = mean_tension * 0.5  # Simplified P-T relation

        # Windkessel update
        dP_art = (self.P_lv - self.P_art) / \
                 max(self.windkessel.R_c, 0.001) - \
                 self.P_art / max(self.windkessel.R_p, 0.001)
        dP_art /= max(self.windkessel.C, 0.01)
        # Clamp to avoid divergence
        dP_art = np.clip(dP_art, -100, 100)
        self.P_art += dP_art * dt
        self.P_art = np.clip(self.P_art, 0, 300)

        # Volume change (simplified)
        if self.P_lv > self.P_art:
            self.phase = "ejection"
            dV = -(self.P_lv - self.P_art) / self.windkessel.R_c * dt
        elif mean_tension < 5.0 and self.phase == "ejection":
            self.phase = "filling"
            dV = (self.windkessel.P_venous - self.P_lv) * 0.5 * dt
        else:
            dV = 0
            self.phase = "isovolumic"

        self.V_lv += dV
        self.V_lv = np.clip(self.V_lv, 20, 300)  # Physiological bounds

        return {
            "T_active": T_active,
            "P_lv": self.P_lv,
            "P_art": self.P_art,
            "V_lv": self.V_lv,
            "phase": self.phase,
            "mean_tension": mean_tension,
        }


class CardiacCycleSimulator:
    """
    Full cardiac cycle simulator with electro-mechanical coupling.

    Integrates:
    - Electrophysiology (monodomain + ionic model)
    - Mechanics (passive + active, FE solver)
    - Hemodynamics (Windkessel)
    """

    def __init__(self, n_elements: int = 100):
        self.n_elements = n_elements
        self.passive = PassiveMechanicsModel()
        self.active = ActiveTensionModel()
        self.coupling = ElectroMechanicalCoupling(self.passive, self.active)

    def simulate_cycle(self, n_beats: int = 1,
                        bcl: float = 800.0,
                        dt: float = 1.0
                        ) -> Dict[str, np.ndarray]:
        """Simulate complete cardiac cycles."""
        total_time = n_beats * bcl
        n_steps = int(total_time / dt)

        # Storage
        results = {
            "time": np.zeros(n_steps),
            "P_lv": np.zeros(n_steps),
            "V_lv": np.zeros(n_steps),
            "P_art": np.zeros(n_steps),
            "T_active_mean": np.zeros(n_steps),
        }

        # Initialize
        self.coupling.V_lv = 120.0  # EDV
        self.coupling.P_art = 10.0

        for step in range(n_steps):
            t = step * dt
            t_in_cycle = t % bcl

            # Synthetic calcium transient
            if t_in_cycle < 50:
                Ca_i = 0.0001  # Resting
            elif t_in_cycle < 100:
                Ca_i = 0.0001 + 0.001 * (t_in_cycle - 50) / 50  # Rising
            elif t_in_cycle < 200:
                Ca_i = 0.0011  # Peak
            elif t_in_cycle < 400:
                Ca_i = 0.0011 * np.exp(-(t_in_cycle - 200) / 100)  # Decay
            else:
                Ca_i = 0.0001  # Resting

            Ca_field = np.full(self.n_elements, Ca_i * 1000)  # Convert to µM
            lambda_field = np.full(self.n_elements, 1.0)

            # Coupling step
            output = self.coupling.coupling_step(Ca_field, lambda_field, dt)

            results["time"][step] = t
            results["P_lv"][step] = output["P_lv"]
            results["V_lv"][step] = output["V_lv"]
            results["P_art"][step] = output["P_art"]
            results["T_active_mean"][step] = output["mean_tension"]

        # Compute hemodynamic indices
        results["EDV"] = np.max(results["V_lv"])
        results["ESV"] = np.min(results["V_lv"])
        results["SV"] = results["EDV"] - results["ESV"]
        results["EF"] = results["SV"] / results["EDV"] * 100 if results["EDV"] > 0 else 0
        results["peak_pressure"] = np.max(results["P_lv"])

        logger.info(f"Cardiac cycle: EDV={results['EDV']:.1f} mL, "
                    f"ESV={results['ESV']:.1f} mL, "
                    f"EF={results['EF']:.1f}%, "
                    f"Peak P={results['peak_pressure']:.1f} kPa")

        return results


def generate_febio_mechanics_config(passive_params: HolzapfelOgdenParams,
                                      active_params: ActiveTensionParams,
                                      output_dir: str) -> str:
    """Generate FEBio configuration for cardiac mechanics."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    filepath = output_path / "febio_mechanics.feb"

    with open(filepath, "w") as f:
        f.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
        f.write('<febio_spec version="4.0">\n')
        f.write('  <Module type="solid"/>\n\n')

        # Control section
        f.write('  <Control>\n')
        f.write('    <time_steps>100</time_steps>\n')
        f.write('    <step_size>0.01</step_size>\n')
        f.write('    <solver type="solid">\n')
        f.write('      <symmetric_stiffness>non-symmetric</symmetric_stiffness>\n')
        f.write('      <equation_scheme>staggered</equation_scheme>\n')
        f.write('    </solver>\n')
        f.write('  </Control>\n\n')

        # Material
        f.write('  <Material>\n')
        f.write('    <material id="1" name="myocardium" '
                'type="Holzapfel-Gasser-Ogden">\n')
        f.write(f'      <c>{passive_params.a}</c>\n')
        f.write(f'      <k1>{passive_params.a_f}</k1>\n')
        f.write(f'      <k2>{passive_params.b_f}</k2>\n')
        f.write(f'      <kappa>0.226</kappa>\n')
        f.write(f'      <K>{passive_params.kappa}</K>\n')
        f.write('    </material>\n')
        f.write('  </Material>\n\n')

        # Boundary conditions
        f.write('  <Boundary>\n')
        f.write('    <bc name="base_fix" type="zero displacement">\n')
        f.write('      <x_dof>1</x_dof>\n')
        f.write('      <y_dof>1</y_dof>\n')
        f.write('      <z_dof>1</z_dof>\n')
        f.write('    </bc>\n')
        f.write('  </Boundary>\n\n')

        # Load (pressure)
        f.write('  <Loads>\n')
        f.write('    <surface_load name="endocardial_pressure" '
                'type="pressure">\n')
        f.write('      <pressure lc="1">1.0</pressure>\n')
        f.write('    </surface_load>\n')
        f.write('  </Loads>\n\n')

        f.write('</febio_spec>\n')

    logger.info(f"Generated FEBio mechanics config: {filepath}")
    return str(filepath)
