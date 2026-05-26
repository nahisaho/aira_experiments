"""
Circuit Specification Module — Formal language for describing
synthetic gene circuits with logic gates and feedback loops.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
from enum import Enum


class GateType(Enum):
    NOT = "NOT"
    AND = "AND"
    OR = "OR"
    NAND = "NAND"
    NOR = "NOR"
    BUFFER = "BUFFER"


class FeedbackType(Enum):
    POSITIVE = "positive"
    NEGATIVE = "negative"
    NONE = "none"


@dataclass
class Signal:
    name: str
    is_input: bool = False
    is_output: bool = False


@dataclass
class LogicGate:
    gate_id: str
    gate_type: GateType
    inputs: List[str]
    output: str
    promoter: str = ""
    rbs: str = ""
    cds: str = ""
    terminator: str = ""


@dataclass
class FeedbackLoop:
    source: str
    target: str
    feedback_type: FeedbackType
    delay: float = 0.0  # minutes


@dataclass
class CircuitSpec:
    """Formal specification of a genetic circuit."""
    name: str
    inputs: List[Signal] = field(default_factory=list)
    outputs: List[Signal] = field(default_factory=list)
    gates: List[LogicGate] = field(default_factory=list)
    feedbacks: List[FeedbackLoop] = field(default_factory=list)

    def add_gate(self, gate: LogicGate):
        self.gates.append(gate)

    def add_feedback(self, fb: FeedbackLoop):
        self.feedbacks.append(fb)

    def get_truth_table(self) -> Dict[Tuple[int, ...], Dict[str, int]]:
        """Compute expected truth table for combinational logic."""
        import itertools
        n = len(self.inputs)
        table = {}
        for bits in itertools.product([0, 1], repeat=n):
            vals = {self.inputs[i].name: bits[i] for i in range(n)}
            for gate in self.gates:
                ins = [vals.get(x, 0) for x in gate.inputs]
                if gate.gate_type == GateType.NOT:
                    vals[gate.output] = 1 - ins[0]
                elif gate.gate_type == GateType.AND:
                    vals[gate.output] = int(all(ins))
                elif gate.gate_type == GateType.OR:
                    vals[gate.output] = int(any(ins))
                elif gate.gate_type == GateType.NAND:
                    vals[gate.output] = 1 - int(all(ins))
                elif gate.gate_type == GateType.NOR:
                    vals[gate.output] = 1 - int(any(ins))
                elif gate.gate_type == GateType.BUFFER:
                    vals[gate.output] = ins[0]
            out = {o.name: vals.get(o.name, 0) for o in self.outputs}
            table[bits] = out
        return table

    def to_verilog(self) -> str:
        """Export circuit spec in Verilog-like format (Cello-compatible)."""
        lines = [f"module {self.name}("]
        ins = ", ".join(s.name for s in self.inputs)
        outs = ", ".join(s.name for s in self.outputs)
        lines.append(f"  input {ins},")
        lines.append(f"  output {outs}")
        lines.append(");")
        for g in self.gates:
            if g.gate_type == GateType.NOT:
                lines.append(f"  assign {g.output} = ~{g.inputs[0]};")
            elif g.gate_type == GateType.AND:
                expr = " & ".join(g.inputs)
                lines.append(f"  assign {g.output} = {expr};")
            elif g.gate_type == GateType.OR:
                expr = " | ".join(g.inputs)
                lines.append(f"  assign {g.output} = {expr};")
            elif g.gate_type == GateType.NAND:
                expr = " & ".join(g.inputs)
                lines.append(f"  assign {g.output} = ~({expr});")
            elif g.gate_type == GateType.NOR:
                expr = " | ".join(g.inputs)
                lines.append(f"  assign {g.output} = ~({expr});")
        lines.append("endmodule")
        return "\n".join(lines)


def make_toggle_switch() -> CircuitSpec:
    """Create the classic Gardner-Collins toggle switch."""
    spec = CircuitSpec(name="toggle_switch")
    spec.inputs = [
        Signal("IPTG", is_input=True),
        Signal("aTc", is_input=True),
    ]
    spec.outputs = [
        Signal("GFP", is_output=True),
    ]
    spec.gates = [
        LogicGate("G1", GateType.NOT, ["TetR_protein"], "LacI_protein",
                  promoter="pTet", rbs="B0034", cds="LacI", terminator="B0015"),
        LogicGate("G2", GateType.NOT, ["LacI_protein"], "TetR_protein",
                  promoter="pLac", rbs="B0034", cds="TetR", terminator="B0015"),
    ]
    spec.feedbacks = [
        FeedbackLoop("LacI_protein", "G2", FeedbackType.NEGATIVE),
        FeedbackLoop("TetR_protein", "G1", FeedbackType.NEGATIVE),
    ]
    return spec


def make_repressilator() -> CircuitSpec:
    """Create the Elowitz-Leibler repressilator."""
    spec = CircuitSpec(name="repressilator")
    spec.outputs = [
        Signal("GFP", is_output=True),
    ]
    spec.gates = [
        LogicGate("G1", GateType.NOT, ["cI_protein"], "LacI_protein",
                  promoter="pLambda", rbs="B0034", cds="LacI", terminator="B0015"),
        LogicGate("G2", GateType.NOT, ["LacI_protein"], "TetR_protein",
                  promoter="pLac", rbs="B0034", cds="TetR", terminator="B0015"),
        LogicGate("G3", GateType.NOT, ["TetR_protein"], "cI_protein",
                  promoter="pTet", rbs="B0034", cds="cI", terminator="B0015"),
    ]
    spec.feedbacks = [
        FeedbackLoop("LacI_protein", "G2", FeedbackType.NEGATIVE),
        FeedbackLoop("TetR_protein", "G3", FeedbackType.NEGATIVE),
        FeedbackLoop("cI_protein", "G1", FeedbackType.NEGATIVE),
    ]
    return spec
