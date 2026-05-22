"""Bridge utilities between the internal power-flow model and PyPSA."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

import numpy as np

from .power_flow import Branch, Bus, BusType, Generator, Load, NetworkModel, PowerFlowResult

try:  # pragma: no cover - optional dependency
    import pypsa
except ImportError:  # pragma: no cover - optional dependency
    pypsa = None


class PyPSABridge:
    """Conversion helpers for PyPSA network models and results."""

    @staticmethod
    def _require_pypsa() -> None:
        if pypsa is None:
            raise ImportError("PyPSA is required for PyPSA bridge operations.")

    @classmethod
    def from_pypsa(cls, network, snapshot: Any | None = None) -> NetworkModel:
        cls._require_pypsa()
        active_snapshot = snapshot if snapshot is not None else (network.snapshots[0] if len(network.snapshots) else None)
        buses: list[Bus] = []
        slack_buses = cls._slack_buses(network)
        pv_buses = set(network.generators.bus.astype(str).tolist()) if not network.generators.empty else set()
        for bus_id, row in network.buses.iterrows():
            bus_type = BusType.SLACK if str(bus_id) in slack_buses else BusType.PV if str(bus_id) in pv_buses else BusType.PQ
            buses.append(
                Bus(
                    id=str(bus_id),
                    type=bus_type,
                    vm_init=float(row.v_mag_pu_set) if "v_mag_pu_set" in row and not np.isnan(row.v_mag_pu_set) else 1.0,
                    va_init=0.0,
                    vm_setpoint=float(row.v_mag_pu_set) if "v_mag_pu_set" in row and not np.isnan(row.v_mag_pu_set) else None,
                    v_min=float(row.v_mag_pu_min) if "v_mag_pu_min" in row and not np.isnan(row.v_mag_pu_min) else 0.9,
                    v_max=float(row.v_mag_pu_max) if "v_mag_pu_max" in row and not np.isnan(row.v_mag_pu_max) else 1.1,
                    name=str(bus_id),
                )
            )
        generators = cls._extract_generators(network, active_snapshot)
        loads = cls._extract_loads(network, active_snapshot) + cls._extract_storage_as_loads(network, active_snapshot)
        branches = cls._extract_branches(network)
        return NetworkModel(buses=buses, branches=branches, generators=generators, loads=loads, base_mva=100.0, name=getattr(network, "name", None))

    @classmethod
    def apply_snapshot(cls, network, snapshot: Any) -> NetworkModel:
        return cls.from_pypsa(network, snapshot=snapshot)

    @classmethod
    def export_results(cls, pypsa_network, result: PowerFlowResult, snapshot: Any | None = None) -> None:
        cls._require_pypsa()
        snap = snapshot if snapshot is not None else (pypsa_network.snapshots[0] if len(pypsa_network.snapshots) else None)
        if snap is None:
            raise ValueError("PyPSA network has no snapshot for result export.")
        if pypsa_network.buses_t.v_mag_pu.empty:
            pypsa_network.buses_t.v_mag_pu = pypsa_network.buses_t.v_mag_pu.reindex(index=pypsa_network.snapshots, columns=pypsa_network.buses.index, fill_value=np.nan)
        if pypsa_network.buses_t.v_ang.empty:
            pypsa_network.buses_t.v_ang = pypsa_network.buses_t.v_ang.reindex(index=pypsa_network.snapshots, columns=pypsa_network.buses.index, fill_value=np.nan)
        for idx, bus_id in enumerate(pypsa_network.buses.index.astype(str)):
            pypsa_network.buses_t.v_mag_pu.loc[snap, bus_id] = result.voltage_magnitudes[idx]
            pypsa_network.buses_t.v_ang.loc[snap, bus_id] = result.voltage_angles[idx]

    @classmethod
    def snapshots(cls, pypsa_network) -> list[Any]:
        cls._require_pypsa()
        return list(pypsa_network.snapshots)

    @staticmethod
    def _slack_buses(network) -> set[str]:
        if hasattr(network, "generators") and not network.generators.empty and "control" in network.generators.columns:
            slack_rows = network.generators.loc[network.generators.control.astype(str).str.lower() == "slack"]
            return set(slack_rows.bus.astype(str).tolist())
        return set()

    @staticmethod
    def _extract_generators(network, snapshot: Any | None) -> list[Generator]:
        generators: list[Generator] = []
        if network.generators.empty:
            return generators
        p_series = getattr(network.generators_t, "p_set", None)
        for name, row in network.generators.iterrows():
            p_mw = float(row.p_nom)
            if snapshot is not None and p_series is not None and not p_series.empty and name in p_series.columns:
                value = p_series.loc[snapshot, name]
                if not np.isnan(value):
                    p_mw = float(value)
            generators.append(
                Generator(
                    bus=str(row.bus),
                    p_mw=p_mw,
                    q_mvar=0.0,
                    vm_setpoint=float(row.v_mag_pu_set) if "v_mag_pu_set" in row and not np.isnan(row.v_mag_pu_set) else None,
                    status=bool(row.active) if "active" in row else True,
                    name=str(name),
                )
            )
        return generators

    @staticmethod
    def _extract_loads(network, snapshot: Any | None) -> list[Load]:
        loads: list[Load] = []
        if network.loads.empty:
            return loads
        p_series = getattr(network.loads_t, "p_set", None)
        q_series = getattr(network.loads_t, "q_set", None)
        for name, row in network.loads.iterrows():
            p_mw = float(row.p_set) if "p_set" in row and not np.isnan(row.p_set) else 0.0
            q_mvar = float(row.q_set) if "q_set" in row and not np.isnan(row.q_set) else 0.0
            if snapshot is not None and p_series is not None and not p_series.empty and name in p_series.columns:
                val = p_series.loc[snapshot, name]
                if not np.isnan(val):
                    p_mw = float(val)
            if snapshot is not None and q_series is not None and not q_series.empty and name in q_series.columns:
                val = q_series.loc[snapshot, name]
                if not np.isnan(val):
                    q_mvar = float(val)
            loads.append(Load(bus=str(row.bus), p_mw=p_mw, q_mvar=q_mvar, status=bool(row.active) if "active" in row else True, name=str(name)))
        return loads

    @staticmethod
    def _extract_storage_as_loads(network, snapshot: Any | None) -> list[Load]:
        storage_loads: list[Load] = []
        if not hasattr(network, "storage_units") or network.storage_units.empty:
            return storage_loads
        p_series = getattr(network.storage_units_t, "p_set", None)
        for name, row in network.storage_units.iterrows():
            p_mw = float(row.p_set) if "p_set" in row and not np.isnan(row.p_set) else 0.0
            if snapshot is not None and p_series is not None and not p_series.empty and name in p_series.columns:
                val = p_series.loc[snapshot, name]
                if not np.isnan(val):
                    p_mw = float(val)
            if p_mw >= 0.0:
                storage_loads.append(Load(bus=str(row.bus), p_mw=p_mw, q_mvar=0.0, name=f"storage-load-{name}"))
            else:
                storage_loads.append(Load(bus=str(row.bus), p_mw=0.0, q_mvar=0.0, name=f"storage-load-{name}"))
        return storage_loads

    @staticmethod
    def _extract_branches(network) -> list[Branch]:
        branches: list[Branch] = []
        if hasattr(network, "lines") and not network.lines.empty:
            for name, row in network.lines.iterrows():
                x = float(row.x) if "x" in row else 0.0
                if abs(x) < 1e-12:
                    x = 1e-9
                branches.append(
                    Branch(
                        from_bus=str(row.bus0),
                        to_bus=str(row.bus1),
                        r_pu=float(row.r) if "r" in row else 0.0,
                        x_pu=x,
                        b_pu=float(row.b) if "b" in row and not np.isnan(row.b) else 0.0,
                        rate_mva=float(row.s_nom) if "s_nom" in row and not np.isnan(row.s_nom) else None,
                        status=bool(row.active) if "active" in row else True,
                        name=str(name),
                    )
                )
        if hasattr(network, "transformers") and not network.transformers.empty:
            for name, row in network.transformers.iterrows():
                x = float(row.x) if "x" in row else 0.0
                if abs(x) < 1e-12:
                    x = 1e-9
                tap_ratio = float(row.tap_ratio) if "tap_ratio" in row and not np.isnan(row.tap_ratio) else 1.0
                phase_shift = float(row.phase_shift) if "phase_shift" in row and not np.isnan(row.phase_shift) else 0.0
                branches.append(
                    Branch(
                        from_bus=str(row.bus0),
                        to_bus=str(row.bus1),
                        r_pu=float(row.r) if "r" in row else 0.0,
                        x_pu=x,
                        b_pu=float(row.b) if "b" in row and not np.isnan(row.b) else 0.0,
                        tap_ratio=tap_ratio,
                        phase_shift=phase_shift,
                        rate_mva=float(row.s_nom) if "s_nom" in row and not np.isnan(row.s_nom) else None,
                        status=bool(row.active) if "active" in row else True,
                        name=str(name),
                    )
                )
        return branches


__all__ = ["PyPSABridge"]
