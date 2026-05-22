"""Formal language description and analysis utilities for synthetic gene circuits."""

from __future__ import annotations

from dataclasses import dataclass, field
import re
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np


DEFAULT_GATE_PARAMETERS: Dict[str, float] = {"Vmax": 1.0, "K": 1.0, "n": 2.0}
_GATE_ALIASES = {"BUF": "BUFFER"}
_GATE_TYPES = {"AND", "OR", "NOT", "NAND", "NOR", "XNOR", "BUFFER"}


def _normalize_gate_type(gate_type: str) -> str:
    normalized = _GATE_ALIASES.get(gate_type.upper(), gate_type.upper())
    if normalized not in _GATE_TYPES:
        raise ValueError(f"Unsupported gate type: {gate_type}")
    return normalized


def _split_top_level(text: str, delimiter: str = ",") -> List[str]:
    parts: List[str] = []
    current: List[str] = []
    depth = 0
    for char in text:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        if char == delimiter and depth == 0:
            item = "".join(current).strip()
            if item:
                parts.append(item)
            current = []
            continue
        current.append(char)
    tail = "".join(current).strip()
    if tail:
        parts.append(tail)
    return parts


def _coerce_value(value: str) -> Any:
    stripped = value.strip()
    try:
        if any(token in stripped.lower() for token in (".", "e")):
            return float(stripped)
        return int(stripped)
    except ValueError:
        return stripped


def _canonical_cycle(cycle: Sequence[str]) -> Tuple[str, ...]:
    sequence = list(cycle)
    if len(sequence) > 1 and sequence[0] == sequence[-1]:
        sequence = sequence[:-1]
    if not sequence:
        return tuple()
    return min(
        tuple(sequence[index:] + sequence[:index])
        for index in range(len(sequence))
    )


def _find_cycles(adjacency: Mapping[str, Iterable[str]]) -> List[List[str]]:
    nodes = sorted(adjacency)
    cycles: set[Tuple[str, ...]] = set()

    def dfs(start: str, node: str, path: List[str]) -> None:
        for neighbor in adjacency.get(node, ()):  # pragma: no branch - tiny graph helper
            if neighbor == start and len(path) > 1:
                cycles.add(_canonical_cycle(path))
            elif neighbor not in path:
                dfs(start, neighbor, path + [neighbor])

    for node in nodes:
        dfs(node, node, [node])
    return [list(cycle) for cycle in sorted(cycles)]


def _has_path(adjacency: Mapping[str, Iterable[str]], start: str, target: str) -> bool:
    if start == target:
        return True
    visited = {start}
    stack = [start]
    while stack:
        node = stack.pop()
        for neighbor in adjacency.get(node, ()):  # pragma: no branch - tiny graph helper
            if neighbor == target:
                return True
            if neighbor not in visited:
                visited.add(neighbor)
                stack.append(neighbor)
    return False


@dataclass
class Signal:
    """Represents a named signal in a gene circuit."""

    name: str
    signal_type: str = "internal"
    initial_value: float = 0.0
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "signal_type": self.signal_type,
            "initial_value": self.initial_value,
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Signal":
        return cls(
            name=str(data["name"]),
            signal_type=str(data.get("signal_type", "internal")),
            initial_value=float(data.get("initial_value", 0.0)),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass
class Gate:
    """Logic gate with Hill-function-based transfer behavior."""

    name: str
    gate_type: str
    inputs: List[str]
    output: str
    parameters: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.gate_type = _normalize_gate_type(self.gate_type)
        self.inputs = [str(item) for item in self.inputs]
        if not self.output:
            raise ValueError("Gate output cannot be empty")

    @property
    def resolved_parameters(self) -> Dict[str, float]:
        params = dict(DEFAULT_GATE_PARAMETERS)
        params.update(self.parameters)
        return {
            "Vmax": float(params["Vmax"]),
            "K": float(params["K"]),
            "n": float(params["n"]),
        }

    @staticmethod
    def _hill_activation(values: np.ndarray, k: float, n_value: float) -> np.ndarray:
        powered = np.power(values, n_value)
        denominator = np.power(k, n_value) + powered
        return np.divide(powered, denominator, out=np.zeros_like(powered), where=denominator != 0)

    def compute_output(self, input_values: Sequence[float]) -> float:
        """Compute the gate output using Hill-function response curves."""

        params = self.resolved_parameters
        values = np.asarray(list(input_values), dtype=float)
        if values.size == 0:
            raise ValueError(f"Gate {self.name} requires at least one input value")

        activation = self._hill_activation(values, params["K"], params["n"])

        if self.gate_type == "BUFFER":
            normalized = float(np.max(activation))
        elif self.gate_type == "NOT":
            normalized = 1.0 - float(np.max(activation))
        elif self.gate_type == "AND":
            normalized = float(np.prod(activation))
        elif self.gate_type == "OR":
            normalized = 1.0 - float(np.prod(1.0 - activation))
        elif self.gate_type == "NAND":
            normalized = 1.0 - float(np.prod(activation))
        elif self.gate_type == "NOR":
            normalized = float(np.prod(1.0 - activation))
        elif self.gate_type == "XNOR":
            normalized = float(np.prod(activation) + np.prod(1.0 - activation))
        else:  # pragma: no cover - safeguarded by validation
            raise ValueError(f"Unsupported gate type: {self.gate_type}")

        normalized = float(np.clip(normalized, 0.0, 1.0))
        return params["Vmax"] * normalized

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "gate_type": self.gate_type,
            "inputs": list(self.inputs),
            "output": self.output,
            "parameters": dict(self.parameters),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Gate":
        return cls(
            name=str(data["name"]),
            gate_type=str(data["gate_type"]),
            inputs=list(data.get("inputs", [])),
            output=str(data["output"]),
            parameters=dict(data.get("parameters", {})),
        )


@dataclass
class FeedbackLoop:
    """Explicit feedback connection from a signal to a target gate."""

    source_signal: str
    target_gate: str
    parameters: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "source_signal": self.source_signal,
            "target_gate": self.target_gate,
            "parameters": dict(self.parameters),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "FeedbackLoop":
        return cls(
            source_signal=str(data["source_signal"]),
            target_gate=str(data["target_gate"]),
            parameters=dict(data.get("parameters", {})),
        )


@dataclass
class CircuitSpec:
    """Complete circuit specification with topology and export helpers."""

    name: str
    signals: Dict[str, Signal] = field(default_factory=dict)
    gates: Dict[str, Gate] = field(default_factory=dict)
    feedback_loops: List[FeedbackLoop] = field(default_factory=list)

    def add_signal(
        self,
        signal: Signal | str,
        signal_type: str = "internal",
        initial_value: float = 0.0,
        metadata: Optional[Mapping[str, Any]] = None,
    ) -> Signal:
        if isinstance(signal, Signal):
            created = signal
        else:
            created = Signal(
                name=str(signal),
                signal_type=signal_type,
                initial_value=initial_value,
                metadata=dict(metadata or {}),
            )
        if created.name in self.signals:
            existing = self.signals[created.name]
            if existing.signal_type == "internal" and created.signal_type != "internal":
                existing.signal_type = created.signal_type
            if created.metadata:
                existing.metadata.update(created.metadata)
            return existing
        self.signals[created.name] = created
        return created

    def add_gate(
        self,
        gate: Gate | str,
        gate_type: Optional[str] = None,
        inputs: Optional[Sequence[str]] = None,
        output: Optional[str] = None,
        parameters: Optional[Mapping[str, Any]] = None,
    ) -> Gate:
        if isinstance(gate, Gate):
            created = gate
        else:
            if gate_type is None or output is None:
                raise ValueError("gate_type and output are required when constructing a gate")
            created = Gate(
                name=str(gate),
                gate_type=gate_type,
                inputs=list(inputs or []),
                output=output,
                parameters=dict(parameters or {}),
            )
        if created.name in self.gates:
            raise ValueError(f"Gate already exists: {created.name}")
        for signal_name in created.inputs:
            self.add_signal(signal_name)
        self.add_signal(created.output)
        self.gates[created.name] = created
        return created

    def add_feedback(
        self,
        source_signal: str,
        target_gate: str,
        parameters: Optional[Mapping[str, Any]] = None,
    ) -> FeedbackLoop:
        loop = FeedbackLoop(
            source_signal=source_signal,
            target_gate=target_gate,
            parameters=dict(parameters or {}),
        )
        self.add_signal(source_signal)
        self.feedback_loops.append(loop)
        return loop

    def define_feedback(
        self,
        source_signal: str,
        target_gate: str,
        parameters: Optional[Mapping[str, Any]] = None,
    ) -> FeedbackLoop:
        return self.add_feedback(source_signal, target_gate, parameters)

    def feedback_sources_for_gate(self, gate_name: str) -> List[str]:
        return [loop.source_signal for loop in self.feedback_loops if loop.target_gate == gate_name]

    def effective_inputs(self, gate_name: str) -> List[str]:
        gate = self.gates[gate_name]
        return list(gate.inputs) + self.feedback_sources_for_gate(gate_name)

    def _signal_drivers(self) -> Tuple[Dict[str, str], List[str]]:
        drivers: Dict[str, str] = {}
        errors: List[str] = []
        for gate in self.gates.values():
            if gate.output in drivers:
                errors.append(
                    f"Signal '{gate.output}' is driven by multiple gates: {drivers[gate.output]} and {gate.name}"
                )
            else:
                drivers[gate.output] = gate.name
        return drivers, errors

    def _dependency_graph(self) -> Dict[str, List[Tuple[str, str]]]:
        drivers, _ = self._signal_drivers()
        feedforward_edges: List[Tuple[str, str]] = []
        feedback_edges: List[Tuple[str, str]] = []

        for gate in self.gates.values():
            for signal_name in gate.inputs:
                source_gate = drivers.get(signal_name)
                if source_gate is not None:
                    feedforward_edges.append((source_gate, gate.name))

        for loop in self.feedback_loops:
            source_gate = drivers.get(loop.source_signal)
            if source_gate is not None:
                feedback_edges.append((source_gate, loop.target_gate))

        return {"feedforward": feedforward_edges, "feedback": feedback_edges}

    def validate_topology(self, raise_on_error: bool = False) -> Dict[str, Any]:
        """Validate signal/gate references and feedforward-versus-feedback topology."""

        errors: List[str] = []
        warnings: List[str] = []
        drivers, driver_errors = self._signal_drivers()
        errors.extend(driver_errors)

        for gate in self.gates.values():
            for signal_name in gate.inputs:
                if signal_name not in self.signals:
                    errors.append(f"Gate '{gate.name}' references unknown input signal '{signal_name}'")
            if gate.output not in self.signals:
                errors.append(f"Gate '{gate.name}' references unknown output signal '{gate.output}'")

        for loop in self.feedback_loops:
            if loop.target_gate not in self.gates:
                errors.append(f"Feedback target gate does not exist: {loop.target_gate}")
            if loop.source_signal not in self.signals:
                errors.append(f"Feedback source signal does not exist: {loop.source_signal}")
            elif loop.source_signal not in drivers:
                errors.append(
                    f"Feedback source signal '{loop.source_signal}' is not produced by any gate"
                )

        graph_parts = self._dependency_graph()
        gate_names = list(self.gates)
        feedforward_adj = {name: [] for name in gate_names}
        full_adj = {name: [] for name in gate_names}
        for source, target in graph_parts["feedforward"]:
            feedforward_adj[source].append(target)
            full_adj[source].append(target)
        for source, target in graph_parts["feedback"]:
            full_adj.setdefault(source, []).append(target)
            full_adj.setdefault(target, [])

        feedforward_cycles = _find_cycles(feedforward_adj)
        if feedforward_cycles:
            errors.append("Feedforward subgraph must be acyclic")

        full_cycles = _find_cycles(full_adj)
        if full_cycles and not self.feedback_loops:
            errors.append("Circuit contains cycles but no explicit feedback loops")

        for loop in self.feedback_loops:
            source_gate = drivers.get(loop.source_signal)
            if source_gate and not _has_path(full_adj, loop.target_gate, source_gate):
                warnings.append(
                    f"Feedback edge {loop.source_signal} -> {loop.target_gate} does not close a detected cycle"
                )

        result = {
            "is_valid": not errors,
            "errors": errors,
            "warnings": warnings,
            "feedforward_is_dag": not feedforward_cycles,
            "feedforward_cycles": feedforward_cycles,
            "feedback_cycles": full_cycles,
            "edges": graph_parts,
        }
        if raise_on_error and errors:
            raise ValueError("; ".join(errors))
        return result

    def transfer_functions(self, signal_values: Mapping[str, float]) -> Dict[str, float]:
        """Compute the instantaneous transfer response of each gate."""

        outputs: Dict[str, float] = {}
        for gate in self.gates.values():
            effective_inputs = self.effective_inputs(gate.name)
            values = [signal_values.get(name, self.signals[name].initial_value) for name in effective_inputs]
            outputs[gate.name] = gate.compute_output(values)
        return outputs

    def evaluate(
        self,
        signal_values: Optional[Mapping[str, float]] = None,
        iterations: int = 32,
        atol: float = 1e-6,
    ) -> Dict[str, float]:
        """Iteratively evaluate the circuit, including explicit feedback loops."""

        state = {name: signal.initial_value for name, signal in self.signals.items()}
        if signal_values:
            state.update({name: float(value) for name, value in signal_values.items()})

        for _ in range(max(1, iterations)):
            updates = dict(state)
            for gate in self.gates.values():
                input_names = self.effective_inputs(gate.name)
                gate_inputs = [state.get(name, self.signals[name].initial_value) for name in input_names]
                updates[gate.output] = gate.compute_output(gate_inputs)
            max_delta = max(abs(updates[name] - state.get(name, 0.0)) for name in updates)
            state = updates
            if max_delta <= atol:
                break
        return state

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "signals": {name: signal.to_dict() for name, signal in self.signals.items()},
            "gates": {name: gate.to_dict() for name, gate in self.gates.items()},
            "feedback_loops": [loop.to_dict() for loop in self.feedback_loops],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "CircuitSpec":
        spec = cls(name=str(data["name"]))
        for signal_data in data.get("signals", {}).values():
            spec.add_signal(Signal.from_dict(signal_data))
        for gate_data in data.get("gates", {}).values():
            spec.add_gate(Gate.from_dict(gate_data))
        for loop_data in data.get("feedback_loops", []):
            spec.feedback_loops.append(FeedbackLoop.from_dict(loop_data))
        return spec

    @classmethod
    def from_verilog_like(cls, text: str) -> "CircuitSpec":
        """Parse a compact Verilog-like circuit description."""

        cleaned = re.sub(r"//.*", "", text)
        match = re.search(r"module\s+(\w+)\s*\((.*?)\)\s*;(?P<body>.*)endmodule", cleaned, flags=re.S)
        if not match:
            raise ValueError("Could not parse module declaration")

        spec = cls(name=match.group(1))
        header = match.group(2)
        for part in _split_top_level(header):
            port_match = re.fullmatch(r"(input|output)\s+(\w+)", part.strip())
            if not port_match:
                raise ValueError(f"Unsupported port declaration: {part}")
            spec.add_signal(port_match.group(2), signal_type=port_match.group(1))

        body = match.group("body")
        statements = [statement.strip() for statement in body.split(";") if statement.strip()]

        for statement in statements:
            if statement.startswith("assign "):
                assign_match = re.fullmatch(r"assign\s+(\w+)\s*=\s*(\w+)", statement)
                if not assign_match:
                    raise ValueError(f"Unsupported assign statement: {statement}")
                target, source = assign_match.groups()
                gate_name = f"assign_{target}"
                suffix = 1
                while gate_name in spec.gates:
                    suffix += 1
                    gate_name = f"assign_{target}_{suffix}"
                spec.add_gate(gate_name, "BUFFER", [source], target)
                continue

            if statement.startswith("feedback"):
                feedback_match = re.fullmatch(r"feedback\s*\((.*)\)", statement, flags=re.S)
                if not feedback_match:
                    raise ValueError(f"Unsupported feedback declaration: {statement}")
                for edge in _split_top_level(feedback_match.group(1)):
                    edge_match = re.fullmatch(r"(\w+)\s*->\s*(\w+)", edge)
                    if not edge_match:
                        raise ValueError(f"Unsupported feedback edge: {edge}")
                    spec.add_feedback(edge_match.group(1), edge_match.group(2))
                continue

            gate_match = re.fullmatch(r"(\w+)\s+(\w+)\s*\((.*)\)", statement, flags=re.S)
            if not gate_match:
                raise ValueError(f"Unsupported statement: {statement}")
            gate_type, gate_name, port_blob = gate_match.groups()
            port_entries = dict(re.findall(r"\.(\w+)\(([^()]*)\)", port_blob))
            if not port_entries:
                raise ValueError(f"Could not parse gate ports: {statement}")

            inputs: List[str] = []
            output: Optional[str] = None
            parameters: Dict[str, Any] = {}
            for port_name, raw_value in port_entries.items():
                values = [item.strip() for item in _split_top_level(raw_value) if item.strip()]
                port_lower = port_name.lower()
                if port_lower == "out":
                    if len(values) != 1:
                        raise ValueError(f"Gate output must have exactly one signal: {statement}")
                    output = values[0]
                elif port_lower.startswith("in"):
                    inputs.extend(values)
                else:
                    parameters[port_name] = _coerce_value(values[0]) if len(values) == 1 else [_coerce_value(v) for v in values]
            if output is None:
                raise ValueError(f"Gate is missing an .out(...) port: {statement}")
            spec.add_gate(gate_name, gate_type, inputs, output, parameters)

        spec.validate_topology(raise_on_error=True)
        return spec

    def to_sbol_like_dict(self) -> Dict[str, Any]:
        """Export a lightweight SBOL-compatible dictionary representation."""

        component_definitions = [
            {
                "id": signal.name,
                "type": "signal",
                "roles": [signal.signal_type],
                "initial_value": signal.initial_value,
            }
            for signal in self.signals.values()
        ]
        component_definitions.extend(
            {
                "id": gate.name,
                "type": "logic-gate",
                "roles": [gate.gate_type.lower()],
                "parameters": dict(gate.parameters),
            }
            for gate in self.gates.values()
        )

        functional_components = [
            {
                "id": signal.name,
                "definition": signal.name,
                "access": "public",
                "direction": signal.signal_type,
            }
            for signal in self.signals.values()
        ]
        functional_components.extend(
            {
                "id": gate.name,
                "definition": gate.name,
                "access": "private",
                "direction": "none",
            }
            for gate in self.gates.values()
        )

        interaction_types = {
            "BUFFER": "stimulation",
            "AND": "stimulation",
            "OR": "stimulation",
            "NOT": "inhibition",
            "NAND": "inhibition",
            "NOR": "inhibition",
            "XNOR": "logic",
        }
        interactions = []
        for gate in self.gates.values():
            participants = [
                {"participant": signal_name, "role": "input"}
                for signal_name in gate.inputs
            ]
            participants.extend(
                {"participant": signal_name, "role": "feedback-input"}
                for signal_name in self.feedback_sources_for_gate(gate.name)
            )
            participants.append({"participant": gate.output, "role": "product"})
            interactions.append(
                {
                    "id": gate.name,
                    "type": interaction_types[gate.gate_type],
                    "gate_type": gate.gate_type,
                    "participants": participants,
                }
            )

        return {
            "displayId": self.name,
            "componentDefinitions": component_definitions,
            "functionalComponents": functional_components,
            "interactions": interactions,
        }

    @classmethod
    def toggle_switch(cls) -> "CircuitSpec":
        """Construct a canonical toggle-switch circuit."""

        spec = cls(name="toggle_switch")
        spec.add_signal("IPTG", signal_type="input")
        spec.add_signal("aTc", signal_type="input")
        spec.add_signal("LacI")
        spec.add_signal("TetR")
        spec.add_signal("GFP", signal_type="output")
        spec.add_signal("RFP", signal_type="output")
        spec.add_gate("g_lacI", "NOT", ["IPTG"], "LacI")
        spec.add_gate("g_tetR", "NOT", ["aTc"], "TetR")
        spec.add_gate("g_gfp", "BUFFER", ["LacI"], "GFP")
        spec.add_gate("g_rfp", "BUFFER", ["TetR"], "RFP")
        spec.add_feedback("LacI", "g_tetR")
        spec.add_feedback("TetR", "g_lacI")
        return spec

    @classmethod
    def repressilator(cls) -> "CircuitSpec":
        """Construct a three-node repressilator oscillator."""

        spec = cls(name="repressilator")
        for signal_name in ("A", "B", "C"):
            spec.add_signal(signal_name)
        spec.add_signal("GFP", signal_type="output")
        spec.add_gate("gA", "NOT", [], "A")
        spec.add_gate("gB", "NOT", [], "B")
        spec.add_gate("gC", "NOT", [], "C")
        spec.add_gate("g_reporter", "BUFFER", ["C"], "GFP")
        spec.add_feedback("C", "gA")
        spec.add_feedback("A", "gB")
        spec.add_feedback("B", "gC")
        return spec


__all__ = ["Signal", "Gate", "FeedbackLoop", "CircuitSpec"]
