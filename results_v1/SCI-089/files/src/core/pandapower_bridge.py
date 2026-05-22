"""Bridge utilities between the internal power-flow model and pandapower."""

from __future__ import annotations

from dataclasses import asdict
from typing import Any
import warnings

import numpy as np

from .power_flow import AdaptivePowerFlowSolver, Branch, Bus, BusType, Generator, Load, NetworkModel, PowerFlowResult

try:  # pragma: no cover - optional dependency
    import pandapower as pp
    import pandapower.networks as pp_networks
except ImportError:  # pragma: no cover - optional dependency
    pp = None
    pp_networks = None


class PandapowerBridge:
    """Conversion and validation helpers for pandapower interoperability."""

    @staticmethod
    def _require_pandapower() -> None:
        if pp is None:
            raise ImportError("pandapower is required for pandapower bridge operations.")

    @classmethod
    def to_pandapower(cls, network: NetworkModel):
        cls._require_pandapower()
        net = pp.create_empty_network(name=network.name or "internal-network", sn_mva=network.base_mva)
        bus_map: dict[str, int] = {}
        for bus in network.buses:
            bus_map[bus.id] = pp.create_bus(
                net,
                vn_kv=110.0,
                name=bus.name or bus.id,
                max_vm_pu=bus.v_max,
                min_vm_pu=bus.v_min,
            )
        for bus in network.buses:
            idx = bus_map[bus.id]
            if bus.type == BusType.SLACK:
                pp.create_ext_grid(
                    net,
                    bus=idx,
                    vm_pu=bus.vm_setpoint or bus.vm_init,
                    va_degree=np.rad2deg(bus.va_setpoint or bus.va_init),
                    name=bus.name or f"slack-{bus.id}",
                )
        for gen in network.generators:
            if not gen.status:
                continue
            bus_type = network.get_bus(gen.bus).type
            if bus_type == BusType.SLACK:
                continue
            pp.create_gen(
                net,
                bus=bus_map[gen.bus],
                p_mw=gen.p_mw,
                vm_pu=gen.vm_setpoint or network.get_bus(gen.bus).vm_setpoint or network.get_bus(gen.bus).vm_init,
                min_q_mvar=gen.q_min_mvar if gen.q_min_mvar is not None else -1e9,
                max_q_mvar=gen.q_max_mvar if gen.q_max_mvar is not None else 1e9,
                name=gen.name or f"gen-{gen.bus}",
                slack=False,
            )
        for load in network.loads:
            if not load.status:
                continue
            pp.create_load(net, bus=bus_map[load.bus], p_mw=load.p_mw, q_mvar=load.q_mvar, name=load.name or f"load-{load.bus}")
        for branch in network.branches:
            if not branch.status:
                continue
            if abs(branch.tap_ratio - 1.0) > 1e-9 or abs(branch.phase_shift) > 1e-9:
                warnings.warn(
                    f"Branch {branch.name or branch.from_bus + '-' + branch.to_bus} contains transformer settings; conversion uses equivalent impedance.",
                    RuntimeWarning,
                    stacklevel=2,
                )
            kwargs = {
                "from_bus": bus_map[branch.from_bus],
                "to_bus": bus_map[branch.to_bus],
                "rft_pu": branch.r_pu,
                "xft_pu": branch.x_pu,
                "rtf_pu": branch.r_pu,
                "xtf_pu": branch.x_pu,
                "sn_mva": branch.rate_mva or network.base_mva,
                "name": branch.name or f"{branch.from_bus}-{branch.to_bus}",
                "in_service": branch.status,
            }
            if "bf_pu" in pp.create_impedance.__code__.co_varnames:
                kwargs["bf_pu"] = branch.b_pu / 2.0
                kwargs["bt_pu"] = branch.b_pu / 2.0
            pp.create_impedance(net, **kwargs)
        return net

    @classmethod
    def from_pandapower(cls, net) -> NetworkModel:
        cls._require_pandapower()
        buses: list[Bus] = []
        for idx, row in net.bus.iterrows():
            bus_type = BusType.PQ
            vm_setpoint = None
            va_setpoint = None
            if hasattr(net, "ext_grid") and not net.ext_grid.empty and idx in set(net.ext_grid["bus"].tolist()):
                ext_row = net.ext_grid.loc[net.ext_grid["bus"] == idx].iloc[0]
                bus_type = BusType.SLACK
                vm_setpoint = float(ext_row.vm_pu)
                va_setpoint = np.deg2rad(float(ext_row.va_degree)) if "va_degree" in ext_row else 0.0
            elif hasattr(net, "gen") and not net.gen.empty and idx in set(net.gen["bus"].tolist()):
                gen_row = net.gen.loc[net.gen["bus"] == idx].iloc[0]
                bus_type = BusType.PV
                vm_setpoint = float(gen_row.vm_pu)
            buses.append(
                Bus(
                    id=str(idx),
                    type=bus_type,
                    vm_init=1.0,
                    va_init=0.0,
                    vm_setpoint=vm_setpoint,
                    va_setpoint=va_setpoint,
                    v_min=float(row.min_vm_pu) if "min_vm_pu" in row else 0.9,
                    v_max=float(row.max_vm_pu) if "max_vm_pu" in row else 1.1,
                    name=str(row.get("name", idx)),
                )
            )
        generators: list[Generator] = []
        if hasattr(net, "gen"):
            for _, row in net.gen.iterrows():
                generators.append(
                    Generator(
                        bus=str(int(row.bus)),
                        p_mw=float(row.p_mw),
                        vm_setpoint=float(row.vm_pu),
                        q_min_mvar=float(row.min_q_mvar) if "min_q_mvar" in row and not np.isnan(row.min_q_mvar) else None,
                        q_max_mvar=float(row.max_q_mvar) if "max_q_mvar" in row and not np.isnan(row.max_q_mvar) else None,
                        status=bool(row.in_service) if "in_service" in row else True,
                        name=str(row.get("name", f"gen-{row.bus}")),
                    )
                )
        loads: list[Load] = []
        if hasattr(net, "load"):
            for _, row in net.load.iterrows():
                loads.append(
                    Load(
                        bus=str(int(row.bus)),
                        p_mw=float(row.p_mw),
                        q_mvar=float(row.q_mvar),
                        status=bool(row.in_service) if "in_service" in row else True,
                        name=str(row.get("name", f"load-{row.bus}")),
                    )
                )
        branches: list[Branch] = []
        if hasattr(net, "impedance"):
            for _, row in net.impedance.iterrows():
                branches.append(
                    Branch(
                        from_bus=str(int(row.from_bus)),
                        to_bus=str(int(row.to_bus)),
                        r_pu=float(row.rft_pu),
                        x_pu=float(row.xft_pu),
                        b_pu=float(row.bf_pu + row.bt_pu) if "bf_pu" in row and "bt_pu" in row else 0.0,
                        rate_mva=float(row.sn_mva) if "sn_mva" in row else None,
                        status=bool(row.in_service) if "in_service" in row else True,
                        name=str(row.get("name", f"{row.from_bus}-{row.to_bus}")),
                    )
                )
        if hasattr(net, "line") and not net.line.empty:
            base_mva = float(getattr(net, "sn_mva", 100.0))
            for _, row in net.line.iterrows():
                z_base = (float(net.bus.loc[int(row.from_bus), "vn_kv"]) ** 2) / base_mva
                r_pu = float(row.r_ohm_per_km * row.length_km / z_base)
                x_pu = float(row.x_ohm_per_km * row.length_km / z_base)
                b_pu = float(2 * np.pi * 50.0 * row.c_nf_per_km * 1e-9 * row.length_km * z_base) if "c_nf_per_km" in row else 0.0
                branches.append(
                    Branch(
                        from_bus=str(int(row.from_bus)),
                        to_bus=str(int(row.to_bus)),
                        r_pu=r_pu,
                        x_pu=x_pu,
                        b_pu=b_pu,
                        rate_mva=float(row.max_i_ka) * float(net.bus.loc[int(row.from_bus), "vn_kv"]) * np.sqrt(3) if "max_i_ka" in row else None,
                        status=bool(row.in_service) if "in_service" in row else True,
                        name=str(row.get("name", f"line-{row.from_bus}-{row.to_bus}")),
                    )
                )
        return NetworkModel(buses=buses, branches=branches, generators=generators, loads=loads, base_mva=float(getattr(net, "sn_mva", 100.0)), name=getattr(net, "name", None))

    @classmethod
    def validate_against_pandapower(
        cls,
        network: NetworkModel,
        internal_result: PowerFlowResult | None = None,
        solver: AdaptivePowerFlowSolver | None = None,
        runpp_kwargs: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        cls._require_pandapower()
        runpp_kwargs = runpp_kwargs or {"init": "flat"}
        net = cls.to_pandapower(network)
        pp.runpp(net, **runpp_kwargs)
        internal_result = internal_result or (solver or AdaptivePowerFlowSolver()).solve(network)
        pp_vm = net.res_bus.vm_pu.to_numpy(dtype=float)
        pp_va = np.deg2rad(net.res_bus.va_degree.to_numpy(dtype=float))
        vm_error = np.abs(pp_vm - internal_result.voltage_magnitudes)
        va_error = np.abs(pp_va - internal_result.voltage_angles)
        return {
            "pandapower_converged": bool(net.converged),
            "internal_converged": internal_result.converged,
            "max_vm_error": float(vm_error.max(initial=0.0)),
            "max_va_error_rad": float(va_error.max(initial=0.0)),
            "mean_vm_error": float(vm_error.mean()) if vm_error.size else 0.0,
            "mean_va_error_rad": float(va_error.mean()) if va_error.size else 0.0,
            "internal_solver": internal_result.solver,
        }

    @classmethod
    def extract_topology(cls, net) -> dict[str, Any]:
        cls._require_pandapower()
        edges: list[dict[str, Any]] = []
        if hasattr(net, "line"):
            for idx, row in net.line.iterrows():
                edges.append({"component": "line", "index": int(idx), "from_bus": int(row.from_bus), "to_bus": int(row.to_bus), "name": row.get("name")})
        if hasattr(net, "impedance"):
            for idx, row in net.impedance.iterrows():
                edges.append({"component": "impedance", "index": int(idx), "from_bus": int(row.from_bus), "to_bus": int(row.to_bus), "name": row.get("name")})
        if hasattr(net, "trafo"):
            for idx, row in net.trafo.iterrows():
                edges.append({"component": "trafo", "index": int(idx), "from_bus": int(row.hv_bus), "to_bus": int(row.lv_bus), "name": row.get("name")})
        return {
            "bus_count": int(len(net.bus)),
            "edge_count": int(len(edges)),
            "slack_buses": [] if not hasattr(net, "ext_grid") else [int(bus) for bus in net.ext_grid.bus.tolist()],
            "edges": edges,
        }

    @classmethod
    def load_test_network(cls, name: str):
        cls._require_pandapower()
        if pp_networks is None or not hasattr(pp_networks, name):
            raise ValueError(f"Unknown pandapower test network: {name}")
        return getattr(pp_networks, name)()

    @classmethod
    def list_test_networks(cls) -> list[str]:
        cls._require_pandapower()
        return sorted(name for name in dir(pp_networks) if not name.startswith("_") and callable(getattr(pp_networks, name)))


__all__ = ["PandapowerBridge"]
