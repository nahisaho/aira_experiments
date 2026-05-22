"""Transient stability simulation models and utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


@dataclass(slots=True)
class TransientStabilityResult:
    """Container for transient stability simulation outputs."""

    time: np.ndarray
    angles: dict[str, np.ndarray]
    speeds: dict[str, np.ndarray]
    voltages: dict[str, np.ndarray]
    is_stable: bool
    cct: float | None = None


class SynchronousMachineModel:
    """Classical synchronous machine model governed by the swing equation."""

    def __init__(
        self,
        name: str,
        inertia_constant: float,
        damping: float,
        transient_reactance: float,
        mechanical_power: float,
        internal_voltage: complex | float = 1.1,
        initial_angle: float = 0.0,
        initial_speed: float = 1.0,
        bus: str | int | None = None,
    ) -> None:
        self.name = name
        self.H = float(inertia_constant)
        self.D = float(damping)
        self.xd_prime = float(transient_reactance)
        self.mechanical_power = float(mechanical_power)
        self.internal_voltage = complex(internal_voltage)
        self.delta = float(initial_angle)
        self.omega = float(initial_speed)
        self.bus = bus if bus is not None else name

    @property
    def inertia_coefficient(self) -> float:
        """Return the lumped inertia term used in the classical model."""

        return max(2.0 * self.H, 1e-9)

    def compute_electrical_power(
        self,
        voltages: Mapping[str | int, complex],
        network: Any,
    ) -> float:
        """Compute electrical air-gap power from the current network state."""

        network_power = _call_network_method(
            network,
            "compute_electrical_power",
            self,
            voltages,
        )
        if network_power is not None:
            return float(network_power)

        terminal_voltage = complex(
            voltages.get(self.bus, _get_value(network, "infinite_bus_voltage", 1.0 + 0.0j))
        )
        transfer_reactance = _resolve_transfer_reactance(network, self)
        power_angle = self.delta - np.angle(terminal_voltage)
        e_internal = abs(self.internal_voltage)
        v_terminal = abs(terminal_voltage)
        return float((e_internal * v_terminal / max(transfer_reactance, 1e-9)) * np.sin(power_angle))

    def get_state(self) -> dict[str, float | complex | str | int | None]:
        """Return the current machine state and parameters."""

        return {
            "name": self.name,
            "bus": self.bus,
            "delta": self.delta,
            "omega": self.omega,
            "H": self.H,
            "D": self.D,
            "xd_prime": self.xd_prime,
            "mechanical_power": self.mechanical_power,
            "internal_voltage": self.internal_voltage,
        }


class TransientStabilitySimulator:
    """Time-domain simulator for transient stability studies."""

    def __init__(
        self,
        method: str = "rk4",
        angle_limit: float = np.pi,
        speed_limit: float = 0.2,
    ) -> None:
        supported_methods = {"rk4", "modified_euler"}
        if method not in supported_methods:
            raise ValueError(f"Unsupported integration method: {method}")
        self.method = method
        self.angle_limit = float(angle_limit)
        self.speed_limit = float(speed_limit)

    def simulate(
        self,
        network: Any,
        fault: Mapping[str, Any],
        t_end: float = 10.0,
        dt: float = 0.001,
    ) -> TransientStabilityResult:
        """Run a time-domain transient stability simulation."""

        machines = list(_get_machines(network))
        if not machines:
            raise ValueError("Network does not contain any synchronous machines.")

        time = np.arange(0.0, t_end + dt, dt, dtype=float)
        n_machines = len(machines)
        states = np.zeros((time.size, 2 * n_machines), dtype=float)
        states[0, :n_machines] = [machine.delta for machine in machines]
        states[0, n_machines:] = [machine.omega for machine in machines]

        initial_voltages = self._resolve_voltages(network, fault, 0.0, machines, states[0])
        voltage_history = {
            str(bus): np.zeros(time.size, dtype=complex) for bus in initial_voltages
        }
        for bus, value in initial_voltages.items():
            voltage_history[str(bus)][0] = complex(value)

        for index in range(time.size - 1):
            current_time = float(time[index])
            current_state = states[index]
            next_state = self._integrate_step(
                network=network,
                fault=fault,
                machines=machines,
                t=current_time,
                state=current_state,
                dt=dt,
            )
            states[index + 1] = next_state
            bus_voltages = self._resolve_voltages(network, fault, time[index + 1], machines, next_state)
            for bus, value in bus_voltages.items():
                key = str(bus)
                if key not in voltage_history:
                    voltage_history[key] = np.zeros(time.size, dtype=complex)
                voltage_history[key][index + 1] = complex(value)

        for machine, angle, speed in zip(
            machines,
            states[-1, :n_machines],
            states[-1, n_machines:],
        ):
            machine.delta = float(angle)
            machine.omega = float(speed)

        angles = {
            machine.name: np.unwrap(states[:, idx].copy())
            for idx, machine in enumerate(machines)
        }
        speeds = {
            machine.name: states[:, n_machines + idx].copy()
            for idx, machine in enumerate(machines)
        }
        is_stable = self._assess_stability(angles, speeds)

        return TransientStabilityResult(
            time=time,
            angles=angles,
            speeds=speeds,
            voltages=voltage_history,
            is_stable=is_stable,
            cct=fault.get("clearing_time"),
        )

    def compute_cct(
        self,
        network: Any,
        fault: Mapping[str, Any],
        tolerance: float = 1e-3,
        max_iterations: int = 25,
    ) -> float:
        """Compute critical clearing time by bisection on fault clearing time."""

        lower = float(fault.get("min_clearing_time", 0.0))
        upper = float(fault.get("max_clearing_time", 1.0))
        t_end = float(fault.get("cct_simulation_end", 10.0))
        dt = float(fault.get("dt", 0.002))

        lower_fault = dict(fault)
        lower_fault["clearing_time"] = lower
        if not self.simulate(network, lower_fault, t_end=t_end, dt=dt).is_stable:
            return lower

        upper_fault = dict(fault)
        upper_fault["clearing_time"] = upper
        if self.simulate(network, upper_fault, t_end=t_end, dt=dt).is_stable:
            return upper

        for _ in range(max_iterations):
            midpoint = 0.5 * (lower + upper)
            candidate_fault = dict(fault)
            candidate_fault["clearing_time"] = midpoint
            is_stable = self.simulate(network, candidate_fault, t_end=t_end, dt=dt).is_stable
            if is_stable:
                lower = midpoint
            else:
                upper = midpoint
            if upper - lower <= tolerance:
                break

        return lower

    def equal_area_criterion(
        self,
        machine: SynchronousMachineModel,
        fault: Mapping[str, Any],
    ) -> dict[str, float | bool]:
        """Evaluate the equal area criterion for an SMIB equivalent."""

        voltage = float(abs(fault.get("voltage", 1.0)))
        e_internal = float(abs(fault.get("internal_voltage", machine.internal_voltage)))
        p_mech = float(fault.get("mechanical_power", machine.mechanical_power))
        x_pre = float(fault.get("pre_fault_transfer", machine.xd_prime + 0.4))
        x_fault = float(fault.get("fault_transfer", machine.xd_prime + 5.0))
        x_post = float(fault.get("post_fault_transfer", machine.xd_prime + 0.6))

        pmax_pre = e_internal * voltage / max(x_pre, 1e-9)
        pmax_fault = e_internal * voltage / max(x_fault, 1e-9)
        pmax_post = e_internal * voltage / max(x_post, 1e-9)

        if abs(p_mech) >= min(pmax_pre, pmax_post):
            raise ValueError("Mechanical power exceeds transferable electrical power.")

        delta_0 = float(np.arcsin(np.clip(p_mech / pmax_pre, -1.0, 1.0)))
        delta_unstable = float(np.pi - np.arcsin(np.clip(p_mech / pmax_post, -1.0, 1.0)))

        def accelerating_area(delta_c: float) -> float:
            value, _ = quad(lambda delta: p_mech - pmax_fault * np.sin(delta), delta_0, delta_c)
            return float(value)

        def decelerating_area(delta_c: float) -> float:
            value, _ = quad(lambda delta: pmax_post * np.sin(delta) - p_mech, delta_c, delta_unstable)
            return float(value)

        def area_balance(delta_c: float) -> float:
            return accelerating_area(delta_c) - decelerating_area(delta_c)

        critical_angle = float(brentq(area_balance, delta_0 + 1e-6, delta_unstable - 1e-6))
        clearing_angle = float(fault.get("clearing_angle", critical_angle))
        accelerating = accelerating_area(clearing_angle)
        decelerating = decelerating_area(clearing_angle)

        return {
            "delta_0": delta_0,
            "critical_clearing_angle": critical_angle,
            "clearing_angle": clearing_angle,
            "accelerating_area": accelerating,
            "decelerating_area": decelerating,
            "margin": decelerating - accelerating,
            "is_stable": clearing_angle <= critical_angle,
        }

    def _integrate_step(
        self,
        network: Any,
        fault: Mapping[str, Any],
        machines: Sequence[SynchronousMachineModel],
        t: float,
        state: np.ndarray,
        dt: float,
    ) -> np.ndarray:
        if self.method == "modified_euler":
            k1 = self._derivatives(t, state, network, fault, machines)
            predictor = state + dt * k1
            k2 = self._derivatives(t + dt, predictor, network, fault, machines)
            return state + 0.5 * dt * (k1 + k2)

        k1 = self._derivatives(t, state, network, fault, machines)
        k2 = self._derivatives(t + 0.5 * dt, state + 0.5 * dt * k1, network, fault, machines)
        k3 = self._derivatives(t + 0.5 * dt, state + 0.5 * dt * k2, network, fault, machines)
        k4 = self._derivatives(t + dt, state + dt * k3, network, fault, machines)
        return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    def _derivatives(
        self,
        t: float,
        state: np.ndarray,
        network: Any,
        fault: Mapping[str, Any],
        machines: Sequence[SynchronousMachineModel],
    ) -> np.ndarray:
        n_machines = len(machines)
        angles = state[:n_machines]
        speeds = state[n_machines:]
        voltages = self._resolve_voltages(network, fault, t, machines, state)
        derivatives = np.zeros_like(state)

        for idx, machine in enumerate(machines):
            machine.delta = float(angles[idx])
            machine.omega = float(speeds[idx])
            electrical_power = machine.compute_electrical_power(voltages, self._phase_network(network, fault, t))
            derivatives[idx] = speeds[idx] - 1.0
            derivatives[n_machines + idx] = (
                machine.mechanical_power
                - electrical_power
                - machine.D * (speeds[idx] - 1.0)
            ) / machine.inertia_coefficient

        return derivatives

    def _resolve_voltages(
        self,
        network: Any,
        fault: Mapping[str, Any],
        t: float,
        machines: Sequence[SynchronousMachineModel],
        state: np.ndarray,
    ) -> dict[str | int, complex]:
        phase_network = self._phase_network(network, fault, t)
        machine_states = {
            machine.name: {
                "angle": float(state[idx]),
                "speed": float(state[len(machines) + idx]),
            }
            for idx, machine in enumerate(machines)
        }
        solved = _call_network_method(
            phase_network,
            "solve_bus_voltages",
            machine_states,
        )
        if solved is not None:
            return {key: complex(value) for key, value in dict(solved).items()}

        configured = _get_value(phase_network, "bus_voltages")
        if isinstance(configured, Mapping):
            return {key: complex(value) for key, value in configured.items()}

        default_voltage = complex(_get_value(phase_network, "infinite_bus_voltage", 1.0 + 0.0j))
        return {machine.bus: default_voltage for machine in machines}

    def _phase_network(self, network: Any, fault: Mapping[str, Any], t: float) -> Any:
        phase = _determine_fault_phase(fault, t)
        if isinstance(network, Mapping) and phase in network:
            return network[phase]
        candidate = _get_value(network, phase)
        return candidate if candidate is not None else network

    def _assess_stability(
        self,
        angles: Mapping[str, np.ndarray],
        speeds: Mapping[str, np.ndarray],
    ) -> bool:
        angle_matrix = np.vstack([np.asarray(values, dtype=float) for values in angles.values()])
        speed_matrix = np.vstack([np.asarray(values, dtype=float) for values in speeds.values()])
        angle_spread = np.max(angle_matrix.max(axis=0) - angle_matrix.min(axis=0))
        max_speed_deviation = np.max(np.abs(speed_matrix - 1.0))
        return bool(angle_spread <= self.angle_limit and max_speed_deviation <= self.speed_limit)


def _get_value(container: Any, key: str, default: Any = None) -> Any:
    if isinstance(container, Mapping):
        return container.get(key, default)
    return getattr(container, key, default)


def _call_network_method(container: Any, method_name: str, *args: Any) -> Any:
    method = getattr(container, method_name, None)
    if callable(method):
        return method(*args)
    return None


def _resolve_transfer_reactance(network: Any, machine: SynchronousMachineModel) -> float:
    transfer = _get_value(network, "transfer_reactance")
    if isinstance(transfer, Mapping):
        network_reactance = transfer.get(machine.name, transfer.get(machine.bus, 0.0))
    else:
        network_reactance = transfer if transfer is not None else _get_value(network, "x_transfer", 0.0)
    return float(machine.xd_prime + float(network_reactance or 0.0))


def _get_machines(network: Any) -> Sequence[SynchronousMachineModel]:
    machines = _get_value(network, "machines")
    if machines is None:
        raise ValueError("Network is missing a 'machines' collection.")
    return machines


def _determine_fault_phase(fault: Mapping[str, Any], t: float) -> str:
    start = float(fault.get("start_time", fault.get("start", 0.0)))
    clearing = float(fault.get("clearing_time", start))
    if t < start:
        return "pre_fault"
    if t < clearing:
        return "fault"
    return "post_fault"
