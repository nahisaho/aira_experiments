"""Power flow solvers and data models for renewable energy grid simulation."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from copy import deepcopy
from dataclasses import dataclass, field
from enum import Enum
import math
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy import sparse
from scipy.sparse import linalg as sparse_linalg


class PowerFlowError(RuntimeError):
    """Raised when a power flow computation cannot be completed."""


class BusType(str, Enum):
    """Supported power-flow bus categories."""

    SLACK = "slack"
    PV = "pv"
    PQ = "pq"


@dataclass(slots=True)
class Bus:
    """Electrical bus model."""

    id: str
    type: BusType = BusType.PQ
    vm_init: float = 1.0
    va_init: float = 0.0
    vm_setpoint: float | None = None
    va_setpoint: float | None = None
    g_shunt: float = 0.0
    b_shunt: float = 0.0
    v_min: float = 0.9
    v_max: float = 1.1
    name: str | None = None


@dataclass(slots=True)
class Branch:
    """Transmission branch model."""

    from_bus: str
    to_bus: str
    r_pu: float
    x_pu: float
    b_pu: float = 0.0
    tap_ratio: float = 1.0
    phase_shift: float = 0.0
    status: bool = True
    rate_mva: float | None = None
    name: str | None = None


@dataclass(slots=True)
class Generator:
    """Generator injection model."""

    bus: str
    p_mw: float
    q_mvar: float = 0.0
    vm_setpoint: float | None = None
    q_min_mvar: float | None = None
    q_max_mvar: float | None = None
    status: bool = True
    name: str | None = None


@dataclass(slots=True)
class Load:
    """Load demand model."""

    bus: str
    p_mw: float
    q_mvar: float = 0.0
    status: bool = True
    name: str | None = None


@dataclass(slots=True)
class NetworkModel:
    """Container for power-system network data."""

    buses: list[Bus]
    branches: list[Branch] = field(default_factory=list)
    generators: list[Generator] = field(default_factory=list)
    loads: list[Load] = field(default_factory=list)
    base_mva: float = 100.0
    name: str | None = None

    def __post_init__(self) -> None:
        if not self.buses:
            raise ValueError("Network must contain at least one bus.")
        if self.base_mva <= 0:
            raise ValueError("base_mva must be positive.")
        bus_ids = [bus.id for bus in self.buses]
        if len(bus_ids) != len(set(bus_ids)):
            raise ValueError("Bus identifiers must be unique.")
        if not any(bus.type == BusType.SLACK for bus in self.buses):
            raise ValueError("Network requires exactly one slack bus.")
        if sum(bus.type == BusType.SLACK for bus in self.buses) != 1:
            raise ValueError("Network must contain exactly one slack bus.")
        unknown_refs = {
            ref
            for ref in [*(branch.from_bus for branch in self.branches), *(branch.to_bus for branch in self.branches), *(gen.bus for gen in self.generators), *(load.bus for load in self.loads)]
            if ref not in set(bus_ids)
        }
        if unknown_refs:
            raise ValueError(f"Unknown bus references detected: {sorted(unknown_refs)}")

    @property
    def bus_lookup(self) -> dict[str, int]:
        return {bus.id: idx for idx, bus in enumerate(self.buses)}

    def get_bus(self, bus_id: str) -> Bus:
        return self.buses[self.bus_lookup[bus_id]]

    def slack_index(self) -> int:
        return next(idx for idx, bus in enumerate(self.buses) if bus.type == BusType.SLACK)

    def pv_indices(self) -> list[int]:
        return [idx for idx, bus in enumerate(self.buses) if bus.type == BusType.PV]

    def pq_indices(self) -> list[int]:
        return [idx for idx, bus in enumerate(self.buses) if bus.type == BusType.PQ]

    def specified_power(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        p = np.zeros(len(self.buses), dtype=float)
        q = np.zeros(len(self.buses), dtype=float)
        inv_base = 1.0 / self.base_mva
        for gen in self.generators:
            if not gen.status:
                continue
            idx = self.bus_lookup[gen.bus]
            p[idx] += gen.p_mw * inv_base
            q[idx] += gen.q_mvar * inv_base
        for load in self.loads:
            if not load.status:
                continue
            idx = self.bus_lookup[load.bus]
            p[idx] -= load.p_mw * inv_base
            q[idx] -= load.q_mvar * inv_base
        return p, q

    def initial_voltage(self) -> NDArray[np.complex128]:
        values: list[complex] = []
        for bus in self.buses:
            vm = bus.vm_setpoint if bus.vm_setpoint is not None else bus.vm_init
            if bus.type == BusType.SLACK and bus.vm_setpoint is not None:
                vm = bus.vm_setpoint
            va = bus.va_setpoint if bus.va_setpoint is not None else bus.va_init
            values.append(vm * np.exp(1j * va))
        return np.asarray(values, dtype=np.complex128)

    def build_ybus(self) -> sparse.csr_matrix:
        size = len(self.buses)
        ybus = sparse.lil_matrix((size, size), dtype=np.complex128)
        lookup = self.bus_lookup
        for branch in self.branches:
            if not branch.status:
                continue
            if abs(branch.r_pu) < 1e-14 and abs(branch.x_pu) < 1e-14:
                raise ValueError(f"Branch {branch.name or branch.from_bus + '-' + branch.to_bus} has zero impedance.")
            i = lookup[branch.from_bus]
            j = lookup[branch.to_bus]
            z = complex(branch.r_pu, branch.x_pu)
            y = 1.0 / z
            b = 1j * branch.b_pu / 2.0
            tap_mag = branch.tap_ratio if abs(branch.tap_ratio) > 1e-12 else 1.0
            shift = np.deg2rad(branch.phase_shift)
            tap = tap_mag * np.exp(1j * shift)
            yff = (y + b) / (tap * np.conj(tap))
            yft = -y / np.conj(tap)
            ytf = -y / tap
            ytt = y + b
            ybus[i, i] += yff
            ybus[j, j] += ytt
            ybus[i, j] += yft
            ybus[j, i] += ytf
        for idx, bus in enumerate(self.buses):
            if bus.g_shunt or bus.b_shunt:
                ybus[idx, idx] += complex(bus.g_shunt, bus.b_shunt)
        return ybus.tocsr()


@dataclass(slots=True)
class PowerFlowResult:
    """Outputs from a power flow study."""

    voltage_magnitudes: NDArray[np.float64]
    voltage_angles: NDArray[np.float64]
    bus_voltages: NDArray[np.complex128]
    bus_active_power: NDArray[np.float64]
    bus_reactive_power: NDArray[np.float64]
    line_flows: list[dict[str, Any]]
    total_losses_mw: float
    total_losses_mvar: float
    converged: bool
    solver: str
    iterations: int
    mismatch_history: list[float]
    step_history: list[float] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)


class NewtonRaphsonPowerFlow:
    """Newton-Raphson power-flow solver in polar coordinates."""

    def __init__(self, damping_min: float = 1.0 / 64.0, use_sparse: bool = True) -> None:
        self.damping_min = damping_min
        self.use_sparse = use_sparse
        self.iteration_history: list[float] = []
        self.step_history: list[float] = []

    def solve(self, network: NetworkModel, max_iter: int = 50, tol: float = 1e-8) -> PowerFlowResult:
        ybus = network.build_ybus()
        p_spec, q_spec = network.specified_power()
        v = network.initial_voltage().astype(np.complex128)
        angle_buses = network.pv_indices() + network.pq_indices()
        pq = network.pq_indices()
        if not angle_buses:
            raise PowerFlowError("No state variables available for Newton-Raphson solve.")

        self.iteration_history = []
        self.step_history = []
        converged = False
        for _ in range(max_iter):
            p_calc, q_calc = self._power_injections(ybus, v)
            mismatch = np.concatenate((p_spec[angle_buses] - p_calc[angle_buses], q_spec[pq] - q_calc[pq]))
            norm = float(np.linalg.norm(mismatch, ord=np.inf))
            self.iteration_history.append(norm)
            if norm < tol:
                converged = True
                break
            jac = self.build_jacobian(ybus, v, p_calc, q_calc, angle_buses, pq)
            try:
                dx = sparse_linalg.spsolve(jac, mismatch) if self.use_sparse else np.linalg.solve(jac.toarray(), mismatch)
            except Exception as exc:  # pragma: no cover - sparse backends differ by platform
                raise PowerFlowError("Jacobian solve failed.") from exc
            v, step = self._apply_damped_update(network, ybus, v, dx, angle_buses, pq, p_spec, q_spec, norm)
            self.step_history.append(step)
        else:
            p_calc, q_calc = self._power_injections(ybus, v)

        if not converged:
            p_calc, q_calc = self._power_injections(ybus, v)
        return self._build_result(network, ybus, v, p_calc, q_calc, converged, solver="newton-raphson")

    def build_jacobian(
        self,
        ybus: sparse.csr_matrix,
        v: NDArray[np.complex128],
        p_calc: NDArray[np.float64],
        q_calc: NDArray[np.float64],
        angle_buses: Sequence[int],
        pq_buses: Sequence[int],
    ) -> sparse.csr_matrix:
        dense = ybus.toarray()
        g = dense.real
        b = dense.imag
        vm = np.abs(v)
        va = np.angle(v)
        n_angle = len(angle_buses)
        n_pq = len(pq_buses)
        j1 = np.zeros((n_angle, n_angle), dtype=float)
        j2 = np.zeros((n_angle, n_pq), dtype=float)
        j3 = np.zeros((n_pq, n_angle), dtype=float)
        j4 = np.zeros((n_pq, n_pq), dtype=float)

        for row, i in enumerate(angle_buses):
            for col, k in enumerate(angle_buses):
                if i == k:
                    j1[row, col] = -q_calc[i] - b[i, i] * vm[i] ** 2
                else:
                    theta = va[i] - va[k]
                    j1[row, col] = vm[i] * vm[k] * (g[i, k] * math.sin(theta) - b[i, k] * math.cos(theta))
            for col, k in enumerate(pq_buses):
                if i == k:
                    j2[row, col] = p_calc[i] / max(vm[i], 1e-12) + g[i, i] * vm[i]
                else:
                    theta = va[i] - va[k]
                    j2[row, col] = vm[i] * (g[i, k] * math.cos(theta) + b[i, k] * math.sin(theta))

        for row, i in enumerate(pq_buses):
            for col, k in enumerate(angle_buses):
                if i == k:
                    j3[row, col] = p_calc[i] - g[i, i] * vm[i] ** 2
                else:
                    theta = va[i] - va[k]
                    j3[row, col] = -vm[i] * vm[k] * (g[i, k] * math.cos(theta) + b[i, k] * math.sin(theta))
            for col, k in enumerate(pq_buses):
                if i == k:
                    j4[row, col] = q_calc[i] / max(vm[i], 1e-12) - b[i, i] * vm[i]
                else:
                    theta = va[i] - va[k]
                    j4[row, col] = vm[i] * (g[i, k] * math.sin(theta) - b[i, k] * math.cos(theta))

        top = np.hstack((j1, j2))
        bottom = np.hstack((j3, j4)) if n_pq else np.zeros((0, top.shape[1]))
        return sparse.csr_matrix(np.vstack((top, bottom)))

    def update_state(
        self,
        network: NetworkModel,
        v: NDArray[np.complex128],
        dx: NDArray[np.float64],
        angle_buses: Sequence[int],
        pq_buses: Sequence[int],
        step_size: float = 1.0,
    ) -> NDArray[np.complex128]:
        vm = np.abs(v).copy()
        va = np.angle(v).copy()
        n_angle = len(angle_buses)
        va[np.asarray(angle_buses, dtype=int)] += step_size * dx[:n_angle]
        if pq_buses:
            pq_array = np.asarray(pq_buses, dtype=int)
            vm[pq_array] += step_size * dx[n_angle:]
            vm[pq_array] = np.maximum(vm[pq_array], 1e-5)
        for idx, bus in enumerate(network.buses):
            if bus.type in {BusType.SLACK, BusType.PV}:
                target_vm = bus.vm_setpoint
                if target_vm is None:
                    for gen in network.generators:
                        if gen.status and gen.bus == bus.id and gen.vm_setpoint is not None:
                            target_vm = gen.vm_setpoint
                            break
                if target_vm is not None:
                    vm[idx] = target_vm
            if bus.type == BusType.SLACK and bus.va_setpoint is not None:
                va[idx] = bus.va_setpoint
        return vm * np.exp(1j * va)

    def _apply_damped_update(
        self,
        network: NetworkModel,
        ybus: sparse.csr_matrix,
        v: NDArray[np.complex128],
        dx: NDArray[np.float64],
        angle_buses: Sequence[int],
        pq_buses: Sequence[int],
        p_spec: NDArray[np.float64],
        q_spec: NDArray[np.float64],
        current_norm: float,
    ) -> tuple[NDArray[np.complex128], float]:
        trial_step = 1.0
        while trial_step >= self.damping_min:
            candidate = self.update_state(network, v, dx, angle_buses, pq_buses, step_size=trial_step)
            p_new, q_new = self._power_injections(ybus, candidate)
            mismatch = np.concatenate((p_spec[list(angle_buses)] - p_new[list(angle_buses)], q_spec[list(pq_buses)] - q_new[list(pq_buses)]))
            new_norm = float(np.linalg.norm(mismatch, ord=np.inf))
            if new_norm < current_norm:
                return candidate, trial_step
            trial_step *= 0.5
        return self.update_state(network, v, dx, angle_buses, pq_buses, step_size=self.damping_min), self.damping_min

    @staticmethod
    def _power_injections(ybus: sparse.csr_matrix, v: NDArray[np.complex128]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        current = ybus @ v
        s = v * np.conj(current)
        return s.real.astype(float), s.imag.astype(float)

    def _build_result(
        self,
        network: NetworkModel,
        ybus: sparse.csr_matrix,
        v: NDArray[np.complex128],
        p_calc: NDArray[np.float64],
        q_calc: NDArray[np.float64],
        converged: bool,
        solver: str,
        metadata: dict[str, Any] | None = None,
    ) -> PowerFlowResult:
        line_flows, losses = _compute_line_flows(network, v)
        return PowerFlowResult(
            voltage_magnitudes=np.abs(v),
            voltage_angles=np.angle(v),
            bus_voltages=v,
            bus_active_power=p_calc * network.base_mva,
            bus_reactive_power=q_calc * network.base_mva,
            line_flows=line_flows,
            total_losses_mw=losses.real,
            total_losses_mvar=losses.imag,
            converged=converged,
            solver=solver,
            iterations=len(self.iteration_history),
            mismatch_history=list(self.iteration_history),
            step_history=list(self.step_history),
            metadata=metadata or {},
        )


class HolomorphicEmbeddingPowerFlow:
    """Holomorphic Embedding Load Flow Method with Padé continuation."""

    def __init__(self) -> None:
        self.last_series: list[NDArray[np.complex128]] = []

    def solve(self, network: NetworkModel, order: int = 30) -> PowerFlowResult:
        if order < 2:
            raise ValueError("HELM requires order >= 2.")
        ybus = network.build_ybus().toarray()
        slack = network.slack_index()
        non_slack = [idx for idx in range(len(network.buses)) if idx != slack]
        if not non_slack:
            raise PowerFlowError("HELM requires at least one non-slack bus.")

        v_series = self.compute_germ(network, ybus, order)
        self.last_series = v_series
        voltages = np.zeros(len(network.buses), dtype=np.complex128)
        voltages[slack] = v_series[0][slack]
        for idx in non_slack:
            coeffs = np.array([series[idx] for series in v_series], dtype=np.complex128)
            try:
                _, _, value = self.pade_approximant(coeffs, x=1.0)
            except Exception:
                value = np.sum(coeffs)
            bus = network.buses[idx]
            target_vm = bus.vm_setpoint
            if target_vm is None:
                for gen in network.generators:
                    if gen.status and gen.bus == bus.id and gen.vm_setpoint is not None:
                        target_vm = gen.vm_setpoint
                        break
            if bus.type == BusType.PV and target_vm is not None and abs(value) > 1e-12:
                value = target_vm * np.exp(1j * np.angle(value))
            voltages[idx] = value

        ybus_sparse = sparse.csr_matrix(ybus)
        p_calc, q_calc = NewtonRaphsonPowerFlow._power_injections(ybus_sparse, voltages)
        mismatch = self._mismatch_norm(network, p_calc, q_calc)
        if np.isfinite(mismatch) and mismatch > 1e-8:
            voltages, p_calc, q_calc, mismatch, corrections = self._refine_solution(network, ybus_sparse, voltages)
        else:
            corrections = []
        line_flows, losses = _compute_line_flows(network, voltages)
        return PowerFlowResult(
            voltage_magnitudes=np.abs(voltages),
            voltage_angles=np.angle(voltages),
            bus_voltages=voltages,
            bus_active_power=p_calc * network.base_mva,
            bus_reactive_power=q_calc * network.base_mva,
            line_flows=line_flows,
            total_losses_mw=losses.real,
            total_losses_mvar=losses.imag,
            converged=np.isfinite(mismatch) and mismatch < 1e-5,
            solver="helm",
            iterations=order + len(corrections),
            mismatch_history=[mismatch] if not corrections else corrections,
            step_history=[],
            metadata={"series_order": order, "correction_iterations": len(corrections)},
        )

    def compute_germ(self, network: NetworkModel, ybus: NDArray[np.complex128], order: int) -> list[NDArray[np.complex128]]:
        slack = network.slack_index()
        non_slack = [idx for idx in range(len(network.buses)) if idx != slack]
        y_red = ybus[np.ix_(non_slack, non_slack)]
        coupling = ybus[np.ix_(non_slack, [slack])].reshape(-1)
        slack_voltage = network.initial_voltage()[slack]
        v0_ns = np.linalg.solve(y_red, -coupling * slack_voltage)
        v0 = np.zeros(len(network.buses), dtype=np.complex128)
        v0[slack] = slack_voltage
        v0[non_slack] = v0_ns
        series = [np.zeros(len(network.buses), dtype=np.complex128) for _ in range(order + 1)]
        series[0] = v0

        p_spec, q_spec = network.specified_power()
        s_spec = p_spec + 1j * q_spec
        reciprocal_series: list[NDArray[np.complex128]] = [np.zeros(len(network.buses), dtype=np.complex128) for _ in range(order)]
        reciprocal_series[0] = self._reciprocal_first_term(np.conj(v0))
        rhs_cache = np.zeros((order, len(non_slack)), dtype=np.complex128)

        for n in range(1, order + 1):
            rhs = np.conj(s_spec[non_slack]) * reciprocal_series[n - 1][non_slack]
            rhs_cache[n - 1, :] = rhs
            vn = np.linalg.solve(y_red, rhs)
            full_vn = np.zeros(len(network.buses), dtype=np.complex128)
            full_vn[non_slack] = vn
            series[n] = full_vn
            if n < order:
                reciprocal_series[n] = self._next_reciprocal_term(series, reciprocal_series, n)
        return series

    def pade_approximant(
        self,
        coefficients: Sequence[complex],
        x: complex = 1.0,
        numerator_order: int | None = None,
        denominator_order: int | None = None,
    ) -> tuple[NDArray[np.complex128], NDArray[np.complex128], complex]:
        coeffs = np.asarray(coefficients, dtype=np.complex128)
        if coeffs.ndim != 1 or coeffs.size < 2:
            raise ValueError("Padé approximant requires a one-dimensional coefficient vector of length >= 2.")
        if numerator_order is None or denominator_order is None:
            denominator_order = (coeffs.size - 1) // 2
            numerator_order = coeffs.size - 1 - denominator_order
        m = int(denominator_order)
        l = int(numerator_order)
        if l < 0 or m < 0 or l + m >= coeffs.size:
            raise ValueError("Invalid Padé orders for supplied coefficient length.")
        if m == 0:
            numerator = coeffs[: l + 1]
            denominator = np.array([1.0 + 0j], dtype=np.complex128)
            return numerator, denominator, np.polyval(numerator[::-1], x)

        matrix = np.zeros((m, m), dtype=np.complex128)
        rhs = np.zeros(m, dtype=np.complex128)
        for row in range(m):
            rhs[row] = -coeffs[l + row + 1]
            for col in range(m):
                idx = l + row - col
                matrix[row, col] = coeffs[idx] if idx >= 0 else 0.0
        q_tail = np.linalg.lstsq(matrix, rhs, rcond=None)[0]
        denominator = np.concatenate(([1.0 + 0j], q_tail))
        numerator = np.zeros(l + 1, dtype=np.complex128)
        for n in range(l + 1):
            numerator[n] = sum(denominator[k] * coeffs[n - k] for k in range(min(n, m) + 1))
        num_value = np.polyval(numerator[::-1], x)
        den_value = np.polyval(denominator[::-1], x)
        if abs(den_value) < 1e-12:
            raise ZeroDivisionError("Padé approximant denominator is numerically singular.")
        return numerator, denominator, num_value / den_value

    @staticmethod
    def _reciprocal_first_term(v0_conj: NDArray[np.complex128]) -> NDArray[np.complex128]:
        out = np.zeros_like(v0_conj)
        mask = np.abs(v0_conj) > 1e-12
        out[mask] = 1.0 / v0_conj[mask]
        return out

    @staticmethod
    def _next_reciprocal_term(
        voltage_series: Sequence[NDArray[np.complex128]],
        reciprocal_series: Sequence[NDArray[np.complex128]],
        order: int,
    ) -> NDArray[np.complex128]:
        out = np.zeros_like(voltage_series[0])
        base = np.conj(voltage_series[0])
        mask = np.abs(base) > 1e-12
        for idx, valid in enumerate(mask):
            if not valid:
                continue
            accum = 0.0j
            for k in range(1, order + 1):
                accum += np.conj(voltage_series[k][idx]) * reciprocal_series[order - k][idx]
            out[idx] = -accum / base[idx]
        return out

    def _refine_solution(
        self,
        network: NetworkModel,
        ybus: sparse.csr_matrix,
        voltages: NDArray[np.complex128],
        max_iter: int = 12,
        tol: float = 1e-8,
    ) -> tuple[NDArray[np.complex128], NDArray[np.float64], NDArray[np.float64], float, list[float]]:
        nr = NewtonRaphsonPowerFlow()
        p_spec, q_spec = network.specified_power()
        angle_buses = network.pv_indices() + network.pq_indices()
        pq = network.pq_indices()
        history: list[float] = []
        p_calc, q_calc = nr._power_injections(ybus, voltages)
        mismatch = self._mismatch_norm(network, p_calc, q_calc)
        if not angle_buses:
            return voltages, p_calc, q_calc, mismatch, history
        for _ in range(max_iter):
            history.append(mismatch)
            if mismatch < tol:
                break
            mismatch_vector = np.concatenate((p_spec[angle_buses] - p_calc[angle_buses], q_spec[pq] - q_calc[pq]))
            jac = nr.build_jacobian(ybus, voltages, p_calc, q_calc, angle_buses, pq)
            try:
                dx = sparse_linalg.spsolve(jac, mismatch_vector)
            except Exception:
                break
            voltages, _ = nr._apply_damped_update(network, ybus, voltages, dx, angle_buses, pq, p_spec, q_spec, mismatch)
            p_calc, q_calc = nr._power_injections(ybus, voltages)
            mismatch = self._mismatch_norm(network, p_calc, q_calc)
        if not history or history[-1] != mismatch:
            history.append(mismatch)
        return voltages, p_calc, q_calc, mismatch, history

    @staticmethod
    def _mismatch_norm(network: NetworkModel, p_calc: NDArray[np.float64], q_calc: NDArray[np.float64]) -> float:
        p_spec, q_spec = network.specified_power()
        angle_buses = network.pv_indices() + network.pq_indices()
        pq = network.pq_indices()
        mismatch = np.concatenate((p_spec[angle_buses] - p_calc[angle_buses], q_spec[pq] - q_calc[pq]))
        if mismatch.size == 0:
            return 0.0
        return float(np.linalg.norm(mismatch, ord=np.inf))


class AdaptivePowerFlowSolver:
    """Solver that adaptively chooses between Newton-Raphson and HELM."""

    def __init__(
        self,
        nr_solver: NewtonRaphsonPowerFlow | None = None,
        helm_solver: HolomorphicEmbeddingPowerFlow | None = None,
        condition_threshold: float = 1e10,
    ) -> None:
        self.nr_solver = nr_solver or NewtonRaphsonPowerFlow()
        self.helm_solver = helm_solver or HolomorphicEmbeddingPowerFlow()
        self.condition_threshold = condition_threshold

    def solve(self, network: NetworkModel) -> PowerFlowResult:
        metric = self._network_condition_metric(network)
        prefer_helm = not np.isfinite(metric) or metric > self.condition_threshold
        if not prefer_helm:
            try:
                nr_result = self.nr_solver.solve(network)
                if nr_result.converged:
                    nr_result.metadata["selection"] = "nr"
                    nr_result.metadata["condition_metric"] = metric
                    return nr_result
            except PowerFlowError:
                pass
        helm_result = self.helm_solver.solve(network)
        helm_result.metadata["selection"] = "helm"
        helm_result.metadata["condition_metric"] = metric
        return helm_result

    def solve_timeseries(
        self,
        network: NetworkModel,
        timesteps: Iterable[NetworkModel | Mapping[str, Any]],
        max_workers: int | None = None,
    ) -> list[PowerFlowResult]:
        tasks = list(timesteps)
        if not tasks:
            return []
        results: list[PowerFlowResult | None] = [None] * len(tasks)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {
                executor.submit(self.solve, self._network_for_timestep(network, item)): idx
                for idx, item in enumerate(tasks)
            }
            for future in as_completed(future_map):
                idx = future_map[future]
                results[idx] = future.result()
        return [result for result in results if result is not None]

    def _network_for_timestep(self, network: NetworkModel, timestep: NetworkModel | Mapping[str, Any]) -> NetworkModel:
        if isinstance(timestep, NetworkModel):
            return timestep
        updated = deepcopy(network)
        load_scale = float(timestep.get("load_scale", 1.0))
        generation_scale = float(timestep.get("generation_scale", 1.0))
        load_overrides = timestep.get("loads", {})
        generator_overrides = timestep.get("generators", {})
        for load in updated.loads:
            load.p_mw *= load_scale
            load.q_mvar *= load_scale
            if load.name and load.name in load_overrides:
                override = load_overrides[load.name]
                load.p_mw = float(override.get("p_mw", load.p_mw))
                load.q_mvar = float(override.get("q_mvar", load.q_mvar))
        for gen in updated.generators:
            gen.p_mw *= generation_scale
            gen.q_mvar *= generation_scale
            if gen.name and gen.name in generator_overrides:
                override = generator_overrides[gen.name]
                gen.p_mw = float(override.get("p_mw", gen.p_mw))
                gen.q_mvar = float(override.get("q_mvar", gen.q_mvar))
        return updated

    def _network_condition_metric(self, network: NetworkModel) -> float:
        ybus = network.build_ybus().toarray()
        slack = network.slack_index()
        non_slack = [idx for idx in range(len(network.buses)) if idx != slack]
        if not non_slack:
            return 0.0
        reduced = ybus[np.ix_(non_slack, non_slack)]
        try:
            return float(np.linalg.cond(reduced))
        except np.linalg.LinAlgError:
            return float("inf")


def _compute_line_flows(network: NetworkModel, voltages: NDArray[np.complex128]) -> tuple[list[dict[str, Any]], complex]:
    lookup = network.bus_lookup
    flows: list[dict[str, Any]] = []
    total_losses = 0.0 + 0.0j
    for branch in network.branches:
        if not branch.status:
            continue
        i = lookup[branch.from_bus]
        j = lookup[branch.to_bus]
        vi = voltages[i]
        vj = voltages[j]
        z = complex(branch.r_pu, branch.x_pu)
        y = 1.0 / z
        b = 1j * branch.b_pu / 2.0
        tap_mag = branch.tap_ratio if abs(branch.tap_ratio) > 1e-12 else 1.0
        shift = np.deg2rad(branch.phase_shift)
        tap = tap_mag * np.exp(1j * shift)
        yff = (y + b) / (tap * np.conj(tap))
        yft = -y / np.conj(tap)
        ytf = -y / tap
        ytt = y + b
        i_ij = yff * vi + yft * vj
        i_ji = ytf * vi + ytt * vj
        s_ij = vi * np.conj(i_ij) * network.base_mva
        s_ji = vj * np.conj(i_ji) * network.base_mva
        loss = s_ij + s_ji
        total_losses += loss
        flows.append(
            {
                "from_bus": branch.from_bus,
                "to_bus": branch.to_bus,
                "p_from_mw": float(s_ij.real),
                "q_from_mvar": float(s_ij.imag),
                "p_to_mw": float(s_ji.real),
                "q_to_mvar": float(s_ji.imag),
                "loss_mw": float(loss.real),
                "loss_mvar": float(loss.imag),
                "loading_pct": None if branch.rate_mva in (None, 0) else float(abs(s_ij) / branch.rate_mva * 100.0),
            }
        )
    return flows, total_losses


__all__ = [
    "AdaptivePowerFlowSolver",
    "Branch",
    "Bus",
    "BusType",
    "Generator",
    "HolomorphicEmbeddingPowerFlow",
    "Load",
    "NetworkModel",
    "NewtonRaphsonPowerFlow",
    "PowerFlowError",
    "PowerFlowResult",
]
