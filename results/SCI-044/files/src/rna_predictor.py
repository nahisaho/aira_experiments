"""
RNA Secondary Structure Prediction with Integrated Approaches
- Turner nearest-neighbor thermodynamic model
- Pseudoknot detection via hierarchical decomposition
- SHAPE/DMS chemical probing constraint integration
- MSA-based covariation scoring
- Riboswitch structure-function prediction
"""

import numpy as np
import math
from typing import List, Tuple, Dict, Optional

# ============================================================
# Turner Nearest-Neighbor Energy Parameters (simplified, kcal/mol)
# Based on Mathews et al., 2004 (doi:10.1073/pnas.0401799101)
# ============================================================

# Stacking energies for canonical base pairs (5'->3' / 3'->5')
STACKING_ENERGIES = {
    ('AU', 'AU'): -0.9, ('AU', 'CG'): -2.2, ('AU', 'GC'): -2.1, ('AU', 'UA'): -1.1,
    ('AU', 'GU'): -1.4, ('AU', 'UG'): -0.6,
    ('CG', 'AU'): -2.1, ('CG', 'CG'): -3.3, ('CG', 'GC'): -2.4, ('CG', 'UA'): -2.1,
    ('CG', 'GU'): -2.1, ('CG', 'UG'): -1.4,
    ('GC', 'AU'): -2.4, ('GC', 'CG'): -3.4, ('GC', 'GC'): -3.3, ('GC', 'UA'): -2.2,
    ('GC', 'GU'): -2.5, ('GC', 'UG'): -1.5,
    ('UA', 'AU'): -1.3, ('UA', 'CG'): -2.4, ('UA', 'GC'): -2.1, ('UA', 'UA'): -0.9,
    ('UA', 'GU'): -1.3, ('UA', 'UG'): -0.5,
    ('GU', 'AU'): -1.3, ('GU', 'CG'): -2.5, ('GU', 'GC'): -2.1, ('GU', 'UA'): -1.4,
    ('GU', 'GU'): -0.5, ('GU', 'UG'): 1.3,
    ('UG', 'AU'): -0.6, ('UG', 'CG'): -1.4, ('UG', 'GC'): -1.5, ('UG', 'UA'): -0.6,
    ('UG', 'GU'): -0.5, ('UG', 'UG'): -0.3,
}

# Loop initiation energies (kcal/mol)
HAIRPIN_INIT = {3: 5.4, 4: 5.6, 5: 5.7, 6: 5.4, 7: 6.0, 8: 5.5, 9: 6.4}
INTERNAL_INIT = {1: 1.0, 2: 0.5, 3: 1.6, 4: 1.7, 5: 1.8, 6: 2.0}
BULGE_INIT = {1: 3.8, 2: 2.8, 3: 3.2, 4: 3.6, 5: 4.0, 6: 4.4}
MULTILOOP_PENALTY = 3.4
MULTILOOP_BRANCH = 0.4
MULTILOOP_UNPAIRED = 0.0

# Terminal AU/GU penalty
TERMINAL_AU_PENALTY = 0.5

CANONICAL_PAIRS = {('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G'), ('G', 'U'), ('U', 'G')}


def can_pair(b1: str, b2: str) -> bool:
    return (b1.upper(), b2.upper()) in CANONICAL_PAIRS


def get_pair_type(b1: str, b2: str) -> str:
    return b1.upper() + b2.upper()


def stacking_energy(bp1: str, bp2: str) -> float:
    return STACKING_ENERGIES.get((bp1, bp2), 0.0)


def hairpin_energy(loop_size: int) -> float:
    if loop_size in HAIRPIN_INIT:
        return HAIRPIN_INIT[loop_size]
    return HAIRPIN_INIT[9] + 1.75 * 0.616 * math.log(loop_size / 9.0)


def internal_loop_energy(size: int) -> float:
    if size in INTERNAL_INIT:
        return INTERNAL_INIT[size]
    return INTERNAL_INIT[6] + 1.75 * 0.616 * math.log(size / 6.0)


def bulge_energy(size: int) -> float:
    if size in BULGE_INIT:
        return BULGE_INIT[size]
    return BULGE_INIT[6] + 1.75 * 0.616 * math.log(size / 6.0)


# ============================================================
# Core Zuker-style DP Algorithm for MFE Structure Prediction
# ============================================================

class RNAStructurePredictor:
    """
    Zuker-style dynamic programming for RNA secondary structure prediction
    with Turner energy model and optional SHAPE/DMS constraints.
    """

    def __init__(self, sequence: str, shape_data: Optional[List[float]] = None,
                 dms_data: Optional[List[float]] = None,
                 shape_slope: float = 1.8, shape_intercept: float = -0.6,
                 covariation_scores: Optional[np.ndarray] = None):
        self.seq = sequence.upper().replace('T', 'U')
        self.n = len(self.seq)
        self.shape_data = shape_data
        self.dms_data = dms_data
        self.shape_slope = shape_slope
        self.shape_intercept = shape_intercept
        self.covariation_scores = covariation_scores

        self.W = np.full(self.n, np.inf)       # minimum energy for [0..i]
        self.V = np.full((self.n, self.n), np.inf)  # energy of closed structure (i,j)
        self.traceback_V = {}
        self.traceback_W = {}
        self.min_hairpin = 3

    def shape_constraint(self, i: int, j: int) -> float:
        """Convert SHAPE reactivity to pseudo-energy (Deigan et al. method)."""
        penalty = 0.0
        if self.shape_data is not None:
            ri = self.shape_data[i] if i < len(self.shape_data) else 0.0
            rj = self.shape_data[j] if j < len(self.shape_data) else 0.0
            penalty += self.shape_slope * ri + self.shape_intercept
            penalty += self.shape_slope * rj + self.shape_intercept
        return penalty

    def dms_constraint(self, i: int, j: int) -> float:
        """Convert DMS reactivity to pseudo-energy penalty."""
        penalty = 0.0
        if self.dms_data is not None:
            di = self.dms_data[i] if i < len(self.dms_data) else 0.0
            dj = self.dms_data[j] if j < len(self.dms_data) else 0.0
            # High DMS reactivity at A/C = unpaired → penalize pairing
            if self.seq[i] in ('A', 'C'):
                penalty += 2.0 * di
            if self.seq[j] in ('A', 'C'):
                penalty += 2.0 * dj
        return penalty

    def covariation_bonus(self, i: int, j: int) -> float:
        """Apply covariation-based bonus from MSA analysis."""
        if self.covariation_scores is not None and i < self.covariation_scores.shape[0] and j < self.covariation_scores.shape[1]:
            return -self.covariation_scores[i, j]  # negative = bonus
        return 0.0

    def compute_V(self, i: int, j: int) -> float:
        """Compute energy of a closed structure between paired positions i and j."""
        if j - i - 1 < self.min_hairpin:
            return np.inf
        if not can_pair(self.seq[i], self.seq[j]):
            return np.inf

        best = np.inf
        best_trace = None
        bp_ij = get_pair_type(self.seq[i], self.seq[j])

        # Hairpin loop
        loop_size = j - i - 1
        e_hairpin = hairpin_energy(loop_size)
        if bp_ij.startswith(('A', 'U')) or bp_ij.startswith(('G', 'U')) or bp_ij.endswith(('A', 'U', 'G')):
            if bp_ij[0] in ('A', 'U') or bp_ij[1] in ('A', 'U'):
                pass  # simplified penalty already in params
        e_hairpin += self.shape_constraint(i, j) + self.dms_constraint(i, j) + self.covariation_bonus(i, j)
        if e_hairpin < best:
            best = e_hairpin
            best_trace = ('H', i, j)

        # Stacking / Internal / Bulge loops
        for p in range(i + 1, min(i + 31, j)):
            for q in range(max(p + self.min_hairpin + 1, j - 30), j):
                if not can_pair(self.seq[p], self.seq[q]):
                    continue
                if self.V[p][q] == np.inf:
                    continue

                left_size = p - i - 1
                right_size = j - q - 1

                if left_size == 0 and right_size == 0:
                    # Stacking
                    bp_pq = get_pair_type(self.seq[p], self.seq[q])
                    e_stack = stacking_energy(bp_ij, bp_pq) + self.V[p][q]
                    e_stack += self.shape_constraint(i, j) + self.covariation_bonus(i, j)
                    if e_stack < best:
                        best = e_stack
                        best_trace = ('S', p, q)
                elif left_size == 0 or right_size == 0:
                    # Bulge loop
                    bsize = max(left_size, right_size)
                    e_bulge = bulge_energy(bsize) + self.V[p][q]
                    e_bulge += self.shape_constraint(i, j) + self.covariation_bonus(i, j)
                    if e_bulge < best:
                        best = e_bulge
                        best_trace = ('B', p, q)
                else:
                    # Internal loop
                    isize = left_size + right_size
                    e_int = internal_loop_energy(isize) + self.V[p][q]
                    e_int += self.shape_constraint(i, j) + self.covariation_bonus(i, j)
                    if e_int < best:
                        best = e_int
                        best_trace = ('I', p, q)

        # Multiloop
        for k in range(i + self.min_hairpin + 2, j - self.min_hairpin - 1):
            if self.V[i + 1][k] < np.inf and self.V[k + 1][j - 1] < np.inf:
                e_multi = MULTILOOP_PENALTY + 2 * MULTILOOP_BRANCH + self.V[i + 1][k] + self.V[k + 1][j - 1]
                e_multi += self.covariation_bonus(i, j)
                if e_multi < best:
                    best = e_multi
                    best_trace = ('M', k)

        self.V[i][j] = best
        if best_trace is not None:
            self.traceback_V[(i, j)] = best_trace
        return best

    def fold(self) -> Tuple[float, str]:
        """Run the DP algorithm to find MFE structure."""
        # Fill V table (bottom-up by span length)
        for span in range(self.min_hairpin + 2, self.n):
            for i in range(0, self.n - span):
                j = i + span
                self.compute_V(i, j)

        # Fill W array
        self.W[0] = 0.0
        for j in range(1, self.n):
            self.W[j] = self.W[j - 1]  # j unpaired
            self.traceback_W[j] = ('U', j - 1)
            for i in range(0, j - self.min_hairpin):
                if self.V[i][j] < np.inf:
                    e = (self.W[i - 1] if i > 0 else 0.0) + self.V[i][j]
                    if e < self.W[j]:
                        self.W[j] = e
                        self.traceback_W[j] = ('P', i)

        # Traceback
        structure = ['.'] * self.n
        self._traceback_W(self.n - 1, structure)
        return self.W[self.n - 1], ''.join(structure)

    def _traceback_W(self, j: int, structure: List[str]):
        if j < 0:
            return
        if j not in self.traceback_W:
            return
        trace = self.traceback_W[j]
        if trace[0] == 'U':
            self._traceback_W(j - 1, structure)
        elif trace[0] == 'P':
            i = trace[1]
            structure[i] = '('
            structure[j] = ')'
            self._traceback_V(i, j, structure)
            if i > 0:
                self._traceback_W(i - 1, structure)

    def _traceback_V(self, i: int, j: int, structure: List[str]):
        if (i, j) not in self.traceback_V:
            return
        trace = self.traceback_V[(i, j)]
        if trace[0] == 'H':
            pass  # hairpin, nothing inside
        elif trace[0] in ('S', 'B', 'I'):
            p, q = trace[1], trace[2]
            structure[p] = '('
            structure[q] = ')'
            self._traceback_V(p, q, structure)
        elif trace[0] == 'M':
            k = trace[1]
            structure[i + 1] = '('
            structure[k] = ')'
            structure[k + 1] = '('
            structure[j - 1] = ')'
            self._traceback_V(i + 1, k, structure)
            self._traceback_V(k + 1, j - 1, structure)


# ============================================================
# Pseudoknot Detection via Hierarchical Decomposition
# ============================================================

class PseudoknotDetector:
    """
    Hierarchical approach to detect pseudoknots after initial prediction.
    Uses O(n²l²) algorithm on sections between stems.
    """

    def __init__(self, sequence: str, base_structure: str):
        self.seq = sequence.upper().replace('T', 'U')
        self.n = len(self.seq)
        self.structure = list(base_structure)
        self.pairs = self._parse_pairs()

    def _parse_pairs(self) -> Dict[int, int]:
        pairs = {}
        stack = []
        for i, c in enumerate(self.structure):
            if c == '(':
                stack.append(i)
            elif c == ')':
                if stack:
                    j = stack.pop()
                    pairs[j] = i
                    pairs[i] = j
        return pairs

    def find_unpaired_regions(self) -> List[Tuple[int, int]]:
        """Identify contiguous unpaired regions (sections)."""
        regions = []
        start = None
        for i in range(self.n):
            if self.structure[i] == '.':
                if start is None:
                    start = i
            else:
                if start is not None:
                    if i - start >= 4:
                        regions.append((start, i - 1))
                    start = None
        if start is not None and self.n - start >= 4:
            regions.append((start, self.n - 1))
        return regions

    def detect_pseudoknots(self) -> List[Tuple[int, int, int, int]]:
        """
        Scan for potential pseudoknot base pairs between unpaired regions.
        Returns list of (i, j, k, l) describing crossing pairs.
        """
        regions = self.find_unpaired_regions()
        pseudoknots = []

        for r1_idx in range(len(regions)):
            for r2_idx in range(r1_idx + 1, len(regions)):
                r1_start, r1_end = regions[r1_idx]
                r2_start, r2_end = regions[r2_idx]

                # Check if regions can form base pairs
                pk_pairs = []
                max_len = min(r1_end - r1_start + 1, r2_end - r2_start + 1)
                for offset in range(max_len):
                    i = r1_start + offset
                    j = r2_end - offset
                    if j <= i:
                        break
                    if can_pair(self.seq[i], self.seq[j]):
                        pk_pairs.append((i, j))
                    else:
                        break

                if len(pk_pairs) >= 2:
                    # Verify this would create a pseudoknot crossing
                    crosses = False
                    for pi, pj in pk_pairs:
                        for existing_i, existing_j in self.pairs.items():
                            if existing_i < existing_j:
                                if (existing_i < pi < existing_j < pj) or (pi < existing_i < pj < existing_j):
                                    crosses = True
                                    break
                        if crosses:
                            break
                    if crosses and len(pk_pairs) >= 2:
                        pseudoknots.append((pk_pairs[0][0], pk_pairs[0][1],
                                           pk_pairs[-1][0], pk_pairs[-1][1]))
        return pseudoknots

    def add_pseudoknots(self) -> str:
        """Add detected pseudoknots to the structure using [] notation."""
        pks = self.detect_pseudoknots()
        result = list(self.structure)
        for pk in pks:
            i_start, j_start, i_end, j_end = pk
            max_len = min(abs(i_end - i_start) + 1, abs(j_start - j_end) + 1)
            for offset in range(max_len):
                i = i_start + offset
                j = j_start - offset
                if j <= i:
                    break
                if can_pair(self.seq[i], self.seq[j]) and result[i] == '.' and result[j] == '.':
                    result[i] = '['
                    result[j] = ']'
        return ''.join(result)


# ============================================================
# SHAPE/DMS Data Simulator
# ============================================================

def simulate_shape_data(sequence: str, structure: str, noise: float = 0.1) -> List[float]:
    """Simulate SHAPE reactivity based on known structure."""
    reactivities = []
    paired = set()
    stack = []
    for i, c in enumerate(structure):
        if c == '(':
            stack.append(i)
        elif c == ')' and stack:
            j = stack.pop()
            paired.add(i)
            paired.add(j)

    for i in range(len(sequence)):
        if i in paired:
            r = np.random.exponential(0.2)  # low reactivity when paired
        else:
            r = 0.5 + np.random.exponential(0.5)  # high reactivity when unpaired
        r = max(0, r + np.random.normal(0, noise))
        reactivities.append(r)
    return reactivities


def simulate_dms_data(sequence: str, structure: str, noise: float = 0.1) -> List[float]:
    """Simulate DMS reactivity (reactive at unpaired A and C)."""
    reactivities = []
    paired = set()
    stack = []
    for i, c in enumerate(structure):
        if c == '(':
            stack.append(i)
        elif c == ')' and stack:
            j = stack.pop()
            paired.add(i)
            paired.add(j)

    for i in range(len(sequence)):
        base = sequence[i].upper()
        if base in ('A', 'C') and i not in paired:
            r = 0.6 + np.random.exponential(0.4)
        else:
            r = np.random.exponential(0.15)
        r = max(0, r + np.random.normal(0, noise))
        reactivities.append(r)
    return reactivities


# ============================================================
# MSA Covariation Analysis
# ============================================================

def compute_mutual_information(msa: List[str]) -> np.ndarray:
    """
    Compute mutual information matrix from multiple sequence alignment.
    Used to detect covarying positions (potential base pairs).
    """
    n_seq = len(msa)
    seq_len = len(msa[0])
    mi_matrix = np.zeros((seq_len, seq_len))

    bases = ['A', 'C', 'G', 'U', '-']

    for i in range(seq_len):
        for j in range(i + 1, seq_len):
            # Joint and marginal frequencies
            joint = np.zeros((len(bases), len(bases)))
            for s in range(n_seq):
                bi = bases.index(msa[s][i].upper()) if msa[s][i].upper() in bases else 4
                bj = bases.index(msa[s][j].upper()) if msa[s][j].upper() in bases else 4
                joint[bi][bj] += 1

            joint /= n_seq
            marginal_i = joint.sum(axis=1)
            marginal_j = joint.sum(axis=0)

            mi = 0.0
            for a in range(len(bases)):
                for b in range(len(bases)):
                    if joint[a][b] > 0 and marginal_i[a] > 0 and marginal_j[b] > 0:
                        mi += joint[a][b] * math.log2(joint[a][b] / (marginal_i[a] * marginal_j[b]))

            mi_matrix[i][j] = mi
            mi_matrix[j][i] = mi

    return mi_matrix


def generate_synthetic_msa(sequence: str, structure: str, n_sequences: int = 50,
                           mutation_rate: float = 0.15) -> List[str]:
    """Generate synthetic MSA with covarying positions at base pairs."""
    paired = {}
    stack = []
    for i, c in enumerate(structure):
        if c == '(':
            stack.append(i)
        elif c == ')' and stack:
            j = stack.pop()
            paired[j] = i
            paired[i] = j

    complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    wobble = {'G': 'U', 'U': 'G'}

    msa = [sequence]
    seq_len = len(sequence)
    for _ in range(n_sequences - 1):
        new_seq = list(sequence)
        for i in range(seq_len):
            if np.random.random() < mutation_rate:
                if i in paired:
                    j = paired[i]
                    if i < j and j < seq_len:
                        new_base = np.random.choice(['A', 'U', 'G', 'C'])
                        new_seq[i] = new_base
                        if np.random.random() < 0.8:
                            new_seq[j] = complement.get(new_base, 'A')
                        else:
                            new_seq[j] = wobble.get(new_base, complement.get(new_base, 'A'))
                elif i not in paired:
                    new_seq[i] = np.random.choice(['A', 'U', 'G', 'C'])
        msa.append(''.join(new_seq))
    return msa


# ============================================================
# Structure Comparison Metrics
# ============================================================

def parse_pairs_from_structure(structure: str) -> set:
    """Parse base pairs from dot-bracket notation."""
    pairs = set()
    stack = []
    for i, c in enumerate(structure):
        if c in ('(', '['):
            stack.append(i)
        elif c in (')', ']') and stack:
            j = stack.pop()
            pairs.add((j, i))
    return pairs


def compute_metrics(predicted: str, reference: str) -> Dict[str, float]:
    """Compute sensitivity, PPV, and F1 score."""
    pred_pairs = parse_pairs_from_structure(predicted)
    ref_pairs = parse_pairs_from_structure(reference)

    tp = len(pred_pairs & ref_pairs)
    fp = len(pred_pairs - ref_pairs)
    fn = len(ref_pairs - pred_pairs)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * sensitivity * ppv / (sensitivity + ppv) if (sensitivity + ppv) > 0 else 0.0
    mcc_denom = math.sqrt((tp + fp) * (tp + fn) * (1 + fp) * (1 + fn)) if (tp + fp) * (tp + fn) > 0 else 1.0

    return {
        'sensitivity': sensitivity,
        'ppv': ppv,
        'f1': f1,
        'tp': tp, 'fp': fp, 'fn': fn,
        'total_ref_pairs': len(ref_pairs),
        'total_pred_pairs': len(pred_pairs)
    }


# ============================================================
# Riboswitch Structure-Function Predictor
# ============================================================

class RiboswitchPredictor:
    """Predict functional states of riboswitches based on structure."""

    KNOWN_APTAMER_MOTIFS = {
        'TPP': 'UGAGA',
        'SAM': 'GAUC',
        'FMN': 'AGGU',
        'purine': 'UACCU',
    }

    def __init__(self, sequence: str, structure: str):
        self.seq = sequence.upper()
        self.structure = structure

    def identify_functional_elements(self) -> Dict[str, any]:
        """Identify potential riboswitch elements."""
        results = {
            'stem_loops': self._count_stem_loops(),
            'aptamer_candidates': self._scan_aptamer_motifs(),
            'expression_platform': self._identify_expression_platform(),
            'structural_switch_potential': self._assess_switch_potential()
        }
        return results

    def _count_stem_loops(self) -> int:
        count = 0
        in_stem = False
        for c in self.structure:
            if c == '(' and not in_stem:
                in_stem = True
            elif c == '.' and in_stem:
                count += 1
                in_stem = False
        return max(count, 1)

    def _scan_aptamer_motifs(self) -> List[str]:
        found = []
        for name, motif in self.KNOWN_APTAMER_MOTIFS.items():
            if motif in self.seq:
                found.append(name)
        return found

    def _identify_expression_platform(self) -> str:
        # Look for anti-terminator / terminator patterns
        paired_count = self.structure.count('(') + self.structure.count(')')
        unpaired_count = self.structure.count('.')
        ratio = paired_count / max(unpaired_count, 1)
        if ratio > 1.5:
            return "terminator-like (highly structured)"
        elif ratio > 0.8:
            return "balanced (potential switching)"
        else:
            return "anti-terminator-like (mostly unstructured)"

    def _assess_switch_potential(self) -> float:
        """Score 0-1 for potential to undergo conformational switch."""
        n = len(self.structure)
        if n == 0:
            return 0.0
        # Count regions that could refold
        flexible_regions = 0
        window = 10
        for i in range(0, n - window):
            window_struct = self.structure[i:i + window]
            dot_ratio = window_struct.count('.') / window
            if 0.3 < dot_ratio < 0.7:
                flexible_regions += 1
        return min(1.0, flexible_regions / max(n / window, 1))
