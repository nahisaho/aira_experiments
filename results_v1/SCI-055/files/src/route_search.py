"""Multi-step retrosynthetic route search using MCTS and A* algorithms."""

import math
import heapq
import numpy as np
from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass, field
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import Descriptors

from .sa_score import improved_sa_score
from .template_based import TemplateBasedRetroSynth


# Simple set of commercially available building blocks (SMILES)
BUILDING_BLOCKS = {
    "c1ccccc1",         # benzene
    "CC",               # ethane
    "CCO",              # ethanol
    "CC=O",             # acetaldehyde
    "CC(=O)O",          # acetic acid
    "CC(=O)Cl",         # acetyl chloride
    "CN",               # methylamine
    "CCN",              # ethylamine
    "c1ccc(N)cc1",      # aniline
    "c1ccc(O)cc1",      # phenol
    "c1ccc(Br)cc1",     # bromobenzene
    "c1ccc(Cl)cc1",     # chlorobenzene
    "OB(O)c1ccccc1",    # phenylboronic acid
    "C=C",              # ethylene
    "O=C(O)c1ccccc1",   # benzoic acid
    "OC(=O)CC(=O)O",    # malonic acid
    "CC(=O)C",          # acetone
    "c1ccncc1",         # pyridine
    "C1CCNCC1",         # piperidine
    "C1CCOC1",          # THF ring
    "O",                # water
    "N",                # ammonia
    "Cl",               # HCl
    "O=C=O",            # CO2
    "CC(C)=O",          # acetone
    "CC(=O)OCC",        # ethyl acetate
    "CCCC",             # butane
    "CCCCO",            # 1-butanol
    "c1ccc2ccccc2c1",   # naphthalene
    "c1ccc(-c2ccccc2)cc1",  # biphenyl
}


def is_building_block(smiles: str) -> bool:
    """Check if a molecule is a known building block."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    canonical = Chem.MolToSmiles(mol)
    if canonical in BUILDING_BLOCKS:
        return True
    if mol.GetNumHeavyAtoms() <= 4:
        return True
    return False


def estimate_cost(smiles: str) -> float:
    """Heuristic cost estimate for reaching building blocks from given molecule."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return float("inf")
    if is_building_block(smiles):
        return 0.0
    sa = improved_sa_score(mol)
    return sa["sa_score"]


@dataclass
class SynthNode:
    """Node in the retrosynthetic search tree."""
    smiles: str
    depth: int = 0
    parent: Optional["SynthNode"] = None
    children: List["SynthNode"] = field(default_factory=list)
    reaction_name: str = ""
    is_solved: bool = False
    cost: float = 0.0

    def __hash__(self):
        return hash(self.smiles)

    def __eq__(self, other):
        return self.smiles == other.smiles


# ============= MCTS-based Route Search ================

class MCTSNode:
    """Monte Carlo Tree Search node for retrosynthesis."""

    def __init__(self, smiles: str, parent=None, depth: int = 0):
        self.smiles = smiles
        self.parent = parent
        self.depth = depth
        self.children: List["MCTSNode"] = []
        self.visits = 0
        self.value = 0.0
        self.is_terminal = is_building_block(smiles)
        self.is_expanded = False
        self.reaction_used = ""

    def ucb1(self, exploration_weight: float = 1.414) -> float:
        if self.visits == 0:
            return float("inf")
        exploitation = self.value / self.visits
        exploration = exploration_weight * math.sqrt(math.log(self.parent.visits) / self.visits)
        return exploitation + exploration


class MCTSRetroSynthesis:
    """MCTS-based multi-step retrosynthetic route search."""

    def __init__(
        self,
        retro_model: TemplateBasedRetroSynth,
        max_depth: int = 6,
        n_iterations: int = 200,
        exploration_weight: float = 1.414,
    ):
        self.retro_model = retro_model
        self.max_depth = max_depth
        self.n_iterations = n_iterations
        self.exploration_weight = exploration_weight

    def search(self, target_smiles: str) -> Dict:
        """Run MCTS to find retrosynthetic routes."""
        root = MCTSNode(target_smiles)
        best_route = None
        best_score = -float("inf")

        for iteration in range(self.n_iterations):
            # Selection
            node = self._select(root)

            # Expansion
            if not node.is_terminal and not node.is_expanded and node.depth < self.max_depth:
                self._expand(node)

            # Simulation
            reward = self._simulate(node)

            # Backpropagation
            self._backpropagate(node, reward)

            # Track best route
            route = self._extract_best_route(root)
            if route and route["score"] > best_score:
                best_score = route["score"]
                best_route = route

        stats = {
            "target": target_smiles,
            "iterations": self.n_iterations,
            "root_visits": root.visits,
            "root_value": root.value / max(root.visits, 1),
        }

        return {
            "best_route": best_route,
            "stats": stats,
        }

    def _select(self, node: MCTSNode) -> MCTSNode:
        while node.is_expanded and node.children and not node.is_terminal:
            node = max(node.children, key=lambda c: c.ucb1(self.exploration_weight))
        return node

    def _expand(self, node: MCTSNode):
        predictions = self.retro_model.predict(node.smiles, top_k=5)
        for pred in predictions:
            reactant_smiles_list = pred["reactants"].split(".")
            for r_smi in reactant_smiles_list:
                child = MCTSNode(r_smi, parent=node, depth=node.depth + 1)
                child.reaction_used = pred["template_name"]
                node.children.append(child)
        node.is_expanded = True

    def _simulate(self, node: MCTSNode) -> float:
        if node.is_terminal:
            return 1.0
        current_smiles = node.smiles
        depth = node.depth
        while depth < self.max_depth:
            if is_building_block(current_smiles):
                return 1.0
            predictions = self.retro_model.predict(current_smiles, top_k=3)
            if not predictions:
                break
            chosen = predictions[0]
            reactants = chosen["reactants"].split(".")
            current_smiles = reactants[0]
            depth += 1

        mol = Chem.MolFromSmiles(current_smiles)
        if mol is None:
            return 0.0
        sa = improved_sa_score(mol)
        return max(0, 1.0 - sa["sa_score"] / 10.0)

    def _backpropagate(self, node: MCTSNode, reward: float):
        while node is not None:
            node.visits += 1
            node.value += reward
            node = node.parent

    def _extract_best_route(self, root: MCTSNode) -> Optional[Dict]:
        route_steps = []
        self._dfs_best_route(root, route_steps)
        if not route_steps:
            return None

        total_score = sum(s.get("confidence", 0.5) for s in route_steps) / max(len(route_steps), 1)
        all_solved = all(
            is_building_block(r)
            for step in route_steps
            for r in step.get("reactants", "").split(".")
        )

        return {
            "steps": route_steps,
            "num_steps": len(route_steps),
            "score": total_score,
            "all_building_blocks": all_solved,
        }

    def _dfs_best_route(self, node: MCTSNode, steps: List[Dict]):
        if node.is_terminal or not node.children:
            return
        best_child = max(node.children, key=lambda c: c.value / max(c.visits, 1))
        if best_child.reaction_used:
            steps.append({
                "product": node.smiles,
                "reactants": best_child.smiles,
                "reaction": best_child.reaction_used,
                "confidence": best_child.value / max(best_child.visits, 1),
            })
        self._dfs_best_route(best_child, steps)


# ============= A* Route Search ================

@dataclass(order=True)
class AStarState:
    priority: float
    smiles: str = field(compare=False)
    path: List[Dict] = field(default_factory=list, compare=False)
    g_cost: float = field(default=0.0, compare=False)
    depth: int = field(default=0, compare=False)


class AStarRetroSynthesis:
    """A* search-based multi-step retrosynthetic route planning."""

    def __init__(
        self,
        retro_model: TemplateBasedRetroSynth,
        max_depth: int = 6,
        max_iterations: int = 500,
    ):
        self.retro_model = retro_model
        self.max_depth = max_depth
        self.max_iterations = max_iterations

    def search(self, target_smiles: str) -> Dict:
        """Run A* search to find optimal retrosynthetic route."""
        h0 = estimate_cost(target_smiles)
        initial_state = AStarState(
            priority=h0,
            smiles=target_smiles,
            path=[],
            g_cost=0.0,
            depth=0,
        )

        open_set = [initial_state]
        closed_set: Set[str] = set()
        best_route = None
        best_cost = float("inf")
        iterations = 0

        while open_set and iterations < self.max_iterations:
            iterations += 1
            current = heapq.heappop(open_set)

            if current.smiles in closed_set:
                continue
            closed_set.add(current.smiles)

            if is_building_block(current.smiles):
                if current.g_cost < best_cost:
                    best_cost = current.g_cost
                    best_route = {
                        "steps": current.path,
                        "num_steps": len(current.path),
                        "total_cost": current.g_cost,
                        "all_building_blocks": True,
                    }
                continue

            if current.depth >= self.max_depth:
                continue

            predictions = self.retro_model.predict(current.smiles, top_k=5)
            for pred in predictions:
                reactants = pred["reactants"].split(".")
                for r_smi in reactants:
                    if r_smi in closed_set:
                        continue
                    step_cost = 1.0 - pred.get("confidence", 0.5)
                    new_g = current.g_cost + step_cost
                    h = estimate_cost(r_smi)
                    f = new_g + h

                    new_path = current.path + [{
                        "product": current.smiles,
                        "reactants": r_smi,
                        "reaction": pred["template_name"],
                        "confidence": pred["confidence"],
                        "step_cost": step_cost,
                    }]

                    new_state = AStarState(
                        priority=f,
                        smiles=r_smi,
                        path=new_path,
                        g_cost=new_g,
                        depth=current.depth + 1,
                    )
                    heapq.heappush(open_set, new_state)

        return {
            "best_route": best_route,
            "stats": {
                "target": target_smiles,
                "iterations": iterations,
                "explored_nodes": len(closed_set),
                "best_cost": best_cost if best_route else None,
            },
        }
