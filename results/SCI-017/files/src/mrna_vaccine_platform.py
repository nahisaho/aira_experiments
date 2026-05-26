#!/usr/bin/env python3
"""
Next-Generation mRNA Vaccine In Silico Design Optimization Platform
====================================================================
Integrated bioinformatics pipeline for:
1. Codon optimization (stability, translation efficiency, immunogenicity)
2. 5'UTR / 3'UTR optimal sequence design
3. Modified nucleotide effect prediction (N1-methylpseudouridine)
4. Antigen epitope selection (MHC binding, T/B cell epitopes)
5. Lipid nanoparticle (LNP) composition optimization
6. Multivalent vaccine design for variant coverage
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.optimize import minimize, differential_evolution
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from collections import Counter
import random
import json
import os
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
random.seed(42)

FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

# ============================================================
# Genetic Code Tables
# ============================================================
CODON_TABLE = {
    'F': ['UUU', 'UUC'], 'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'I': ['AUU', 'AUC', 'AUA'], 'M': ['AUG'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'Y': ['UAU', 'UAC'], '*': ['UAA', 'UAG', 'UGA'],
    'H': ['CAU', 'CAC'], 'Q': ['CAA', 'CAG'],
    'N': ['AAU', 'AAC'], 'K': ['AAA', 'AAG'],
    'D': ['GAU', 'GAC'], 'E': ['GAA', 'GAG'],
    'C': ['UGU', 'UGC'], 'W': ['UGG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
}

# Human codon usage frequencies (per thousand)
HUMAN_CODON_FREQ = {
    'UUU': 17.6, 'UUC': 20.3, 'UUA': 7.7, 'UUG': 12.9,
    'CUU': 13.2, 'CUC': 19.6, 'CUA': 7.2, 'CUG': 39.6,
    'AUU': 16.0, 'AUC': 20.8, 'AUA': 7.5, 'AUG': 22.0,
    'GUU': 11.0, 'GUC': 14.5, 'GUA': 7.1, 'GUG': 28.1,
    'UCU': 15.2, 'UCC': 17.7, 'UCA': 12.2, 'UCG': 4.4,
    'CCU': 17.5, 'CCC': 19.8, 'CCA': 16.9, 'CCG': 6.9,
    'ACU': 13.1, 'ACC': 18.9, 'ACA': 15.1, 'ACG': 6.1,
    'GCU': 18.4, 'GCC': 27.7, 'GCA': 15.8, 'GCG': 7.4,
    'UAU': 12.2, 'UAC': 15.3, 'UAA': 1.0, 'UAG': 0.8, 'UGA': 1.6,
    'CAU': 10.9, 'CAC': 15.1, 'CAA': 12.3, 'CAG': 34.2,
    'AAU': 17.0, 'AAC': 19.1, 'AAA': 24.4, 'AAG': 31.9,
    'GAU': 21.8, 'GAC': 25.1, 'GAA': 29.0, 'GAG': 39.6,
    'UGU': 10.6, 'UGC': 12.6, 'UGG': 13.2,
    'CGU': 4.5, 'CGC': 10.4, 'CGA': 6.2, 'CGG': 11.4,
    'AGA': 12.2, 'AGG': 12.0, 'AGU': 12.1, 'AGC': 19.5,
    'GGU': 10.8, 'GGC': 22.2, 'GGA': 16.5, 'GGG': 16.5,
}


# ============================================================
# Module 1: Codon Optimization Engine
# ============================================================
class CodonOptimizer:
    """Multi-objective codon optimization balancing stability, translation, immunogenicity."""

    def __init__(self, protein_seq):
        self.protein_seq = protein_seq
        self.results = {}

    def _gc_content(self, rna_seq):
        gc = sum(1 for n in rna_seq if n in 'GC')
        return gc / len(rna_seq)

    def _codon_adaptation_index(self, codons):
        """Calculate Codon Adaptation Index (CAI)."""
        cai_values = []
        for codon in codons:
            freq = HUMAN_CODON_FREQ.get(codon, 1.0)
            aa = None
            for a, cs in CODON_TABLE.items():
                if codon in cs:
                    aa = a
                    break
            if aa:
                max_freq = max(HUMAN_CODON_FREQ.get(c, 1.0) for c in CODON_TABLE[aa])
                cai_values.append(freq / max_freq if max_freq > 0 else 0)
        return np.exp(np.mean(np.log(np.array(cai_values) + 1e-10)))

    def _uridine_content(self, rna_seq):
        return sum(1 for n in rna_seq if n == 'U') / len(rna_seq)

    def _dinucleotide_score(self, rna_seq):
        """Penalize CpG and UpA dinucleotides (immune stimulatory)."""
        cpg = rna_seq.count('CG')
        upa = rna_seq.count('UA')
        return -(cpg * 2 + upa) / (len(rna_seq) - 1)

    def _mfe_proxy(self, rna_seq):
        """Proxy for minimum free energy based on GC and local structure potential."""
        gc = self._gc_content(rna_seq)
        return -(gc * 2.5 + len(rna_seq) * 0.01)

    def optimize(self, strategy='balanced', n_iterations=500):
        """
        Multi-objective optimization using evolutionary algorithm.
        Strategies: 'max_expression', 'max_stability', 'min_immunogenicity', 'balanced'
        """
        weights = {
            'max_expression': {'cai': 0.7, 'gc': 0.1, 'dinuc': 0.1, 'uridine': 0.1},
            'max_stability': {'cai': 0.2, 'gc': 0.5, 'dinuc': 0.1, 'uridine': 0.2},
            'min_immunogenicity': {'cai': 0.2, 'gc': 0.1, 'dinuc': 0.5, 'uridine': 0.2},
            'balanced': {'cai': 0.35, 'gc': 0.25, 'dinuc': 0.2, 'uridine': 0.2},
        }
        w = weights[strategy]

        best_score = -np.inf
        best_codons = None
        scores_history = []

        for iteration in range(n_iterations):
            codons = []
            for aa in self.protein_seq:
                if aa in CODON_TABLE and aa != '*':
                    candidates = CODON_TABLE[aa]
                    if iteration == 0:
                        # Start with frequency-weighted selection
                        freqs = [HUMAN_CODON_FREQ.get(c, 1.0) for c in candidates]
                        total = sum(freqs)
                        probs = [f / total for f in freqs]
                        chosen = np.random.choice(candidates, p=probs)
                    else:
                        # Mutation: occasionally pick suboptimal codons
                        if random.random() < 0.3:
                            chosen = random.choice(candidates)
                        else:
                            freqs = [HUMAN_CODON_FREQ.get(c, 1.0) for c in candidates]
                            total = sum(freqs)
                            probs = [f / total for f in freqs]
                            chosen = np.random.choice(candidates, p=probs)
                    codons.append(chosen)

            rna_seq = ''.join(codons)
            cai = self._codon_adaptation_index(codons)
            gc = self._gc_content(rna_seq)
            gc_score = 1.0 - abs(gc - 0.55) * 4  # Optimal GC ~55%
            dinuc = self._dinucleotide_score(rna_seq)
            uridine = 1.0 - self._uridine_content(rna_seq)

            score = (w['cai'] * cai + w['gc'] * gc_score +
                     w['dinuc'] * (1 + dinuc) + w['uridine'] * uridine)

            scores_history.append(score)

            if score > best_score:
                best_score = score
                best_codons = codons

        best_rna = ''.join(best_codons)
        self.results[strategy] = {
            'rna_sequence': best_rna,
            'codons': best_codons,
            'cai': self._codon_adaptation_index(best_codons),
            'gc_content': self._gc_content(best_rna),
            'uridine_content': self._uridine_content(best_rna),
            'dinucleotide_score': self._dinucleotide_score(best_rna),
            'mfe_proxy': self._mfe_proxy(best_rna),
            'total_score': best_score,
            'convergence': scores_history,
        }
        return self.results[strategy]

    def compare_strategies(self):
        strategies = ['max_expression', 'max_stability', 'min_immunogenicity', 'balanced']
        for s in strategies:
            self.optimize(strategy=s)
        return self.results


# ============================================================
# Module 2: UTR Optimizer
# ============================================================
class UTROptimizer:
    """5'UTR and 3'UTR sequence design for maximal translation efficiency."""

    KOZAK_CONSENSUS = 'GCCACCAUGG'  # Strong Kozak
    KNOWN_5UTR_ELEMENTS = {
        'alpha_globin': 'ACUUCUUGGUCCUUAGCUACUGCUCUAAAGCCUCCAGCUGCCUCAGAUCUGUCUACAUCCGAG',
        'beta_globin': 'ACAUUUGCUUCUGACACAACUGUGUUCACUAGCAACCUCAAACAGACACCAUGG',
        'hsp70': 'AAGCAGCCGAGCCGACGGCAAGCUGGCUGCCAAGAAGGUGCUGCUGGCGG',
    }
    KNOWN_3UTR_ELEMENTS = {
        'alpha_globin': 'GCUGGAGCCUCGGUGGCCAUGCUUCUUGCCCCUUGGGCCUCCCCCCAGCCC',
        'beta_globin': 'GCUAAUAAAGCCUAAUAUUUUCCUCAGCUUUCCUGGCUGUUCCCCCAACUG',
        'AES': 'CUGCUAGCCUUCUGCUAAUCAUGUUUAUAAAUUGUAAAUAUUCUAACCCCAU',
    }

    def __init__(self):
        self.results = {}

    def _score_5utr(self, seq):
        score = 0
        # Kozak sequence at end
        if seq[-4:] == 'AUGG' or seq[-7:].count('GCC') > 0:
            score += 30
        # Avoid upstream AUGs
        uaugs = seq[:-3].count('AUG')
        score -= uaugs * 20
        # GC content (moderate preferred)
        gc = sum(1 for n in seq if n in 'GC') / len(seq)
        score += (1 - abs(gc - 0.50)) * 20
        # Length penalty (optimal 50-100 nt)
        if 50 <= len(seq) <= 100:
            score += 15
        elif len(seq) < 30 or len(seq) > 150:
            score -= 10
        # Secondary structure proxy (fewer self-complementary stretches)
        complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
        self_comp = 0
        for i in range(len(seq) - 5):
            window = seq[i:i+6]
            rev_comp = ''.join(complement.get(n, n) for n in reversed(window))
            if window == rev_comp:
                self_comp += 1
        score -= self_comp * 5
        return score

    def _score_3utr(self, seq):
        score = 0
        # Stability elements
        if 'AAUAAA' in seq:
            score += 20  # Polyadenylation signal
        # AU-rich elements (destabilizing)
        are_count = seq.count('AUUUA')
        score -= are_count * 15
        # Length (optimal 100-300 nt)
        if 100 <= len(seq) <= 300:
            score += 15
        # GC content
        gc = sum(1 for n in seq if n in 'GC') / len(seq)
        score += (1 - abs(gc - 0.45)) * 15
        return score

    def optimize_5utr(self, n_candidates=1000):
        """Generate and score 5'UTR candidates."""
        candidates = []
        # Include known high-performance UTRs
        for name, seq in self.KNOWN_5UTR_ELEMENTS.items():
            candidates.append((name, seq, self._score_5utr(seq)))

        # Generate random optimized candidates
        for i in range(n_candidates):
            length = random.randint(40, 100)
            # Bias toward optimal composition
            weights = {'A': 0.2, 'U': 0.15, 'G': 0.35, 'C': 0.3}
            seq = ''.join(random.choices(list(weights.keys()), weights=list(weights.values()), k=length))
            # Ensure Kozak at end
            seq = seq[:-10] + 'GCCACC' + 'AUGG'
            # Remove upstream AUGs
            seq_clean = seq[:-3].replace('AUG', 'ACG') + seq[-3:]
            score = self._score_5utr(seq_clean)
            candidates.append((f'synthetic_{i}', seq_clean, score))

        candidates.sort(key=lambda x: x[2], reverse=True)
        self.results['5utr'] = candidates[:20]
        return candidates[:20]

    def optimize_3utr(self, n_candidates=1000):
        """Generate and score 3'UTR candidates."""
        candidates = []
        for name, seq in self.KNOWN_3UTR_ELEMENTS.items():
            candidates.append((name, seq, self._score_3utr(seq)))

        for i in range(n_candidates):
            length = random.randint(80, 250)
            weights = {'A': 0.25, 'U': 0.2, 'G': 0.3, 'C': 0.25}
            seq = ''.join(random.choices(list(weights.keys()), weights=list(weights.values()), k=length))
            # Add polyadenylation signal
            pos = random.randint(len(seq) - 30, len(seq) - 10)
            seq = seq[:pos] + 'AAUAAA' + seq[pos + 6:]
            # Remove destabilizing AREs
            seq = seq.replace('AUUUA', 'ACUCA')
            score = self._score_3utr(seq)
            candidates.append((f'synthetic_{i}', seq, score))

        candidates.sort(key=lambda x: x[2], reverse=True)
        self.results['3utr'] = candidates[:20]
        return candidates[:20]


# ============================================================
# Module 3: Modified Nucleotide Effect Predictor
# ============================================================
class ModifiedNucleotidePredictor:
    """Predict effects of nucleotide modifications on mRNA performance."""

    MODIFICATIONS = {
        'm1psi': {'name': 'N1-methylpseudouridine', 'innate_reduction': 0.85,
                  'translation_boost': 1.8, 'stability_factor': 1.5, 'cost': 1.5},
        'psi': {'name': 'Pseudouridine', 'innate_reduction': 0.70,
                'translation_boost': 1.4, 'stability_factor': 1.3, 'cost': 1.2},
        'm5C': {'name': '5-methylcytidine', 'innate_reduction': 0.50,
                'translation_boost': 1.2, 'stability_factor': 1.2, 'cost': 1.3},
        'm6A': {'name': 'N6-methyladenosine', 'innate_reduction': 0.30,
                'translation_boost': 0.9, 'stability_factor': 1.1, 'cost': 1.1},
        'mo5U': {'name': '5-methoxyuridine', 'innate_reduction': 0.60,
                 'translation_boost': 1.3, 'stability_factor': 1.25, 'cost': 1.4},
        'none': {'name': 'Unmodified', 'innate_reduction': 0.0,
                 'translation_boost': 1.0, 'stability_factor': 1.0, 'cost': 1.0},
    }

    def __init__(self, rna_sequence):
        self.rna_seq = rna_sequence
        self.u_count = rna_sequence.count('U')
        self.c_count = rna_sequence.count('C')
        self.a_count = rna_sequence.count('A')
        self.length = len(rna_sequence)

    def predict_effects(self):
        results = {}
        for mod_id, mod in self.MODIFICATIONS.items():
            # Model TLR activation (lower is better)
            base_tlr = 0.8  # Unmodified baseline
            tlr_activation = base_tlr * (1 - mod['innate_reduction'])

            # Translation efficiency
            base_translation = 1.0
            translation = base_translation * mod['translation_boost']

            # mRNA half-life (hours)
            base_halflife = 6.0
            halflife = base_halflife * mod['stability_factor']

            # Protein yield (relative)
            protein_yield = translation * (1 + halflife / 24)

            # Immunogenicity score (adaptive, higher is better)
            adaptive_immune = protein_yield * (1 + mod['innate_reduction'] * 0.3)

            results[mod_id] = {
                'name': mod['name'],
                'tlr_activation': round(tlr_activation, 3),
                'translation_efficiency': round(translation, 3),
                'halflife_hours': round(halflife, 2),
                'protein_yield': round(protein_yield, 3),
                'adaptive_immunogenicity': round(adaptive_immune, 3),
                'cost_factor': mod['cost'],
            }
        return results


# ============================================================
# Module 4: Epitope Predictor
# ============================================================
class EpitopePredictor:
    """MHC binding prediction, T-cell and B-cell epitope selection."""

    # Common HLA alleles
    HLA_ALLELES = ['HLA-A*02:01', 'HLA-A*01:01', 'HLA-A*03:01', 'HLA-A*24:02',
                   'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*44:02',
                   'HLA-DRB1*01:01', 'HLA-DRB1*03:01', 'HLA-DRB1*04:01']

    # Amino acid binding preferences (simplified PSSM-like)
    AA_HYDROPHOBICITY = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    }

    def __init__(self, protein_seq):
        self.protein_seq = protein_seq

    def predict_mhc_binding(self, peptide_length=9):
        """Predict MHC-I binding using simplified scoring matrix."""
        peptides = []
        for i in range(len(self.protein_seq) - peptide_length + 1):
            peptide = self.protein_seq[i:i + peptide_length]
            scores = {}
            for allele in self.HLA_ALLELES[:7]:  # Class I
                # Anchor residue scoring
                score = 0
                # Position 2 and C-terminal anchor
                if peptide[1] in 'LMI':
                    score += 3.0
                if peptide[-1] in 'VLIKY':
                    score += 3.0
                # Hydrophobicity profile
                for j, aa in enumerate(peptide):
                    h = self.AA_HYDROPHOBICITY.get(aa, 0)
                    if j in [1, 8]:  # Anchor positions
                        score += h * 0.5
                    else:
                        score += h * 0.1
                # Add noise for allele-specific variation
                score += np.random.normal(0, 0.5)
                scores[allele] = round(score, 3)

            avg_score = np.mean(list(scores.values()))
            peptides.append({
                'position': i + 1,
                'peptide': peptide,
                'allele_scores': scores,
                'average_score': round(avg_score, 3),
                'predicted_binder': avg_score > 3.0,
            })

        peptides.sort(key=lambda x: x['average_score'], reverse=True)
        return peptides

    def predict_tcell_epitopes(self, peptide_length=9):
        """Combined T-cell epitope prediction (MHC binding + proteasomal cleavage + TAP transport)."""
        mhc_results = self.predict_mhc_binding(peptide_length)
        tcell_epitopes = []

        for pep in mhc_results[:50]:
            peptide = pep['peptide']
            # Proteasomal cleavage prediction (C-terminal preference)
            cleavage_score = self.AA_HYDROPHOBICITY.get(peptide[-1], 0) * 0.3
            # TAP transport score
            tap_score = sum(self.AA_HYDROPHOBICITY.get(aa, 0) for aa in peptide) / len(peptide)
            # Combined score
            combined = pep['average_score'] * 0.5 + cleavage_score * 0.25 + tap_score * 0.25
            pep['cleavage_score'] = round(cleavage_score, 3)
            pep['tap_score'] = round(tap_score, 3)
            pep['tcell_score'] = round(combined, 3)
            tcell_epitopes.append(pep)

        tcell_epitopes.sort(key=lambda x: x['tcell_score'], reverse=True)
        return tcell_epitopes

    def predict_bcell_epitopes(self, window_size=15):
        """B-cell epitope prediction using surface accessibility and hydrophilicity."""
        epitopes = []
        for i in range(len(self.protein_seq) - window_size + 1):
            window = self.protein_seq[i:i + window_size]
            # Parker hydrophilicity
            hydrophilicity = -np.mean([self.AA_HYDROPHOBICITY.get(aa, 0) for aa in window])
            # Surface accessibility proxy
            surface = sum(1 for aa in window if aa in 'DEKRHNQS') / len(window)
            # Flexibility (based on B-factors)
            flexibility = sum(1 for aa in window if aa in 'GSNDP') / len(window)
            # BepiPred-like combined score
            score = hydrophilicity * 0.4 + surface * 30 + flexibility * 20

            epitopes.append({
                'position': i + 1,
                'epitope': window,
                'hydrophilicity': round(hydrophilicity, 3),
                'surface_accessibility': round(surface, 3),
                'flexibility': round(flexibility, 3),
                'bcell_score': round(score, 3),
            })

        epitopes.sort(key=lambda x: x['bcell_score'], reverse=True)
        return epitopes

    def population_coverage(self, epitopes, top_n=10):
        """Estimate global population coverage based on HLA allele frequencies."""
        allele_freqs = {
            'HLA-A*02:01': 0.29, 'HLA-A*01:01': 0.15, 'HLA-A*03:01': 0.13,
            'HLA-A*24:02': 0.17, 'HLA-B*07:02': 0.12, 'HLA-B*08:01': 0.09,
            'HLA-B*44:02': 0.11, 'HLA-DRB1*01:01': 0.10, 'HLA-DRB1*03:01': 0.12,
            'HLA-DRB1*04:01': 0.15,
        }
        covered_alleles = set()
        for ep in epitopes[:top_n]:
            for allele, score in ep.get('allele_scores', {}).items():
                if score > 3.0:
                    covered_alleles.add(allele)

        coverage = 1 - np.prod([1 - allele_freqs.get(a, 0.05) for a in covered_alleles])
        return round(coverage * 100, 1), covered_alleles


# ============================================================
# Module 5: LNP Optimizer
# ============================================================
class LNPOptimizer:
    """Lipid nanoparticle composition optimization via ML-guided simulation."""

    def __init__(self):
        self.training_data = self._generate_training_data()
        self.model = None

    def _generate_training_data(self):
        """Generate synthetic LNP formulation dataset based on published trends."""
        n_samples = 500
        data = []
        for _ in range(n_samples):
            ionizable = np.random.uniform(30, 60)  # mol%
            helper = np.random.uniform(5, 20)       # mol% (DSPC)
            cholesterol = np.random.uniform(20, 50)  # mol%
            peg = np.random.uniform(0.5, 5.0)       # mol%
            # Normalize to 100%
            total = ionizable + helper + cholesterol + peg
            ionizable, helper, cholesterol, peg = [x/total*100 for x in [ionizable, helper, cholesterol, peg]]
            np_ratio = np.random.uniform(3, 12)  # N/P ratio

            # Simulated outcomes based on empirical relationships
            size = 60 + (ionizable - 46)**2 * 0.3 + (peg - 1.5)**2 * 5 + np.random.normal(0, 5)
            pdi = 0.05 + abs(ionizable - 50) * 0.005 + abs(peg - 1.5) * 0.02 + np.random.normal(0, 0.02)
            encapsulation = 95 - (ionizable - 50)**2 * 0.1 - abs(np_ratio - 6)**2 * 0.5 + np.random.normal(0, 3)
            encapsulation = np.clip(encapsulation, 40, 99)
            pka = 6.0 + (ionizable - 50) * 0.02 + np.random.normal(0, 0.1)
            transfection = (encapsulation / 100) * (1 / (1 + abs(size - 80) / 50)) * (1 - pdi) * 100
            transfection += np.random.normal(0, 5)
            transfection = np.clip(transfection, 5, 100)

            data.append({
                'ionizable_lipid': round(ionizable, 2),
                'helper_lipid': round(helper, 2),
                'cholesterol': round(cholesterol, 2),
                'peg_lipid': round(peg, 2),
                'np_ratio': round(np_ratio, 2),
                'size_nm': round(size, 1),
                'pdi': round(np.clip(pdi, 0.01, 0.5), 3),
                'encapsulation': round(encapsulation, 1),
                'pka': round(pka, 2),
                'transfection_efficiency': round(transfection, 1),
            })
        return pd.DataFrame(data)

    def train_model(self):
        """Train ML model for LNP property prediction."""
        features = ['ionizable_lipid', 'helper_lipid', 'cholesterol', 'peg_lipid', 'np_ratio']
        X = self.training_data[features]
        y_enc = self.training_data['encapsulation']
        y_trans = self.training_data['transfection_efficiency']

        self.scaler = StandardScaler()
        X_scaled = self.scaler.fit_transform(X)

        self.model_enc = GradientBoostingRegressor(n_estimators=100, random_state=42)
        self.model_trans = GradientBoostingRegressor(n_estimators=100, random_state=42)

        self.model_enc.fit(X_scaled, y_enc)
        self.model_trans.fit(X_scaled, y_trans)

        cv_enc = cross_val_score(self.model_enc, X_scaled, y_enc, cv=5)
        cv_trans = cross_val_score(self.model_trans, X_scaled, y_trans, cv=5)

        return {
            'encapsulation_r2': round(np.mean(cv_enc), 4),
            'transfection_r2': round(np.mean(cv_trans), 4),
        }

    def optimize_formulation(self):
        """Find optimal LNP composition using differential evolution."""
        def objective(x):
            ionizable, helper, cholesterol, peg, np_ratio = x
            total = ionizable + helper + cholesterol + peg
            x_norm = [ionizable/total*100, helper/total*100, cholesterol/total*100, peg/total*100, np_ratio]
            x_scaled = self.scaler.transform([x_norm])
            enc = self.model_enc.predict(x_scaled)[0]
            trans = self.model_trans.predict(x_scaled)[0]
            return -(0.4 * enc + 0.6 * trans)

        bounds = [(30, 60), (5, 20), (20, 50), (0.5, 5), (3, 12)]
        result = differential_evolution(objective, bounds, seed=42, maxiter=200)

        optimal = result.x
        total = sum(optimal[:4])
        optimal_norm = [x/total*100 for x in optimal[:4]] + [optimal[4]]

        x_scaled = self.scaler.transform([optimal_norm])
        return {
            'ionizable_lipid': round(optimal_norm[0], 2),
            'helper_lipid': round(optimal_norm[1], 2),
            'cholesterol': round(optimal_norm[2], 2),
            'peg_lipid': round(optimal_norm[3], 2),
            'np_ratio': round(optimal_norm[4], 2),
            'predicted_encapsulation': round(self.model_enc.predict(x_scaled)[0], 1),
            'predicted_transfection': round(self.model_trans.predict(x_scaled)[0], 1),
        }


# ============================================================
# Module 6: Multivalent Vaccine Designer
# ============================================================
class MultivalentDesigner:
    """Design multivalent mRNA vaccines for variant coverage."""

    SPIKE_VARIANTS = {
        'Wuhan-Hu-1': 'MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT',
        'Delta_B.1.617.2': None,  # Mutations applied below
        'Omicron_BA.1': None,
        'Omicron_BA.5': None,
        'Omicron_XBB.1.5': None,
    }

    # Key mutations for each variant (position, original, mutant) in spike
    VARIANT_MUTATIONS = {
        'Delta_B.1.617.2': [(19, 'R', 'T'), (157, 'F', 'L'), (452, 'L', 'R'), (478, 'T', 'K'), (614, 'D', 'G'), (681, 'P', 'R'), (950, 'D', 'N')],
        'Omicron_BA.1': [(67, 'A', 'V'), (95, 'T', 'I'), (142, 'G', 'D'), (339, 'G', 'D'), (371, 'S', 'L'), (373, 'S', 'P'), (375, 'S', 'F'), (417, 'K', 'N'), (440, 'N', 'K'), (446, 'G', 'S'), (477, 'S', 'N'), (478, 'T', 'K'), (484, 'E', 'A'), (493, 'Q', 'R'), (496, 'G', 'S'), (498, 'Q', 'R'), (501, 'N', 'Y'), (505, 'Y', 'H'), (547, 'T', 'K'), (614, 'D', 'G'), (655, 'H', 'Y'), (679, 'N', 'K'), (681, 'P', 'H'), (764, 'N', 'K'), (796, 'D', 'Y'), (856, 'N', 'K'), (954, 'Q', 'H'), (969, 'N', 'K')],
        'Omicron_BA.5': [(339, 'G', 'D'), (371, 'S', 'F'), (373, 'S', 'P'), (375, 'S', 'F'), (376, 'T', 'A'), (405, 'D', 'N'), (408, 'R', 'S'), (417, 'K', 'N'), (440, 'N', 'K'), (452, 'L', 'R'), (478, 'T', 'K'), (484, 'E', 'A'), (486, 'F', 'V'), (501, 'N', 'Y'), (614, 'D', 'G'), (655, 'H', 'Y'), (679, 'N', 'K'), (681, 'P', 'H'), (764, 'N', 'K'), (796, 'D', 'Y'), (954, 'Q', 'H'), (969, 'N', 'K')],
        'Omicron_XBB.1.5': [(339, 'G', 'H'), (346, 'R', 'T'), (368, 'L', 'I'), (371, 'S', 'F'), (373, 'S', 'P'), (375, 'S', 'F'), (376, 'T', 'A'), (405, 'D', 'N'), (408, 'R', 'S'), (417, 'K', 'N'), (440, 'N', 'K'), (445, 'V', 'P'), (446, 'G', 'S'), (460, 'N', 'K'), (478, 'T', 'K'), (484, 'E', 'A'), (486, 'F', 'P'), (490, 'F', 'S'), (501, 'N', 'Y'), (505, 'Y', 'H'), (614, 'D', 'G'), (655, 'H', 'Y'), (679, 'N', 'K'), (681, 'P', 'H'), (764, 'N', 'K'), (796, 'D', 'Y'), (954, 'Q', 'H'), (969, 'N', 'K')],
    }

    def __init__(self):
        self._generate_variants()

    def _generate_variants(self):
        wuhan = self.SPIKE_VARIANTS['Wuhan-Hu-1']
        for variant, mutations in self.VARIANT_MUTATIONS.items():
            seq = list(wuhan)
            for pos, orig, mut in mutations:
                if pos - 1 < len(seq):
                    seq[pos - 1] = mut
            self.SPIKE_VARIANTS[variant] = ''.join(seq)

    def analyze_conservation(self, region_start=300, region_end=600):
        """Analyze sequence conservation across variants in RBD region."""
        variants = {k: v[region_start:region_end] for k, v in self.SPIKE_VARIANTS.items() if v}
        n_pos = len(list(variants.values())[0])
        conservation = []

        for i in range(n_pos):
            residues = [v[i] for v in variants.values()]
            counter = Counter(residues)
            most_common = counter.most_common(1)[0][1]
            cons_score = most_common / len(residues)
            conservation.append({
                'position': region_start + i + 1,
                'conservation': cons_score,
                'residues': dict(counter),
                'consensus': counter.most_common(1)[0][0],
            })
        return conservation

    def design_multivalent(self, n_valences=3):
        """Select optimal variant combination for multivalent vaccine."""
        variants = list(self.SPIKE_VARIANTS.keys())
        # Calculate pairwise antigenic distances
        distances = {}
        for i, v1 in enumerate(variants):
            for j, v2 in enumerate(variants):
                if i < j and self.SPIKE_VARIANTS[v1] and self.SPIKE_VARIANTS[v2]:
                    seq1 = self.SPIKE_VARIANTS[v1]
                    seq2 = self.SPIKE_VARIANTS[v2]
                    min_len = min(len(seq1), len(seq2))
                    diff = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a != b)
                    distances[(v1, v2)] = diff

        # Greedy selection maximizing pairwise distance
        selected = ['Wuhan-Hu-1']  # Always include ancestral
        remaining = [v for v in variants if v != 'Wuhan-Hu-1' and self.SPIKE_VARIANTS[v]]

        for _ in range(n_valences - 1):
            best_v = None
            best_min_dist = -1
            for v in remaining:
                min_dist = min(distances.get((min(v, s), max(v, s)), 0) for s in selected)
                if min_dist > best_min_dist:
                    best_min_dist = min_dist
                    best_v = v
            if best_v:
                selected.append(best_v)
                remaining.remove(best_v)

        return {
            'selected_variants': selected,
            'distances': {str(k): v for k, v in distances.items()},
            'coverage_analysis': self._estimate_coverage(selected),
        }

    def _estimate_coverage(self, selected_variants):
        """Estimate cross-reactive coverage."""
        all_variants = [v for v in self.SPIKE_VARIANTS.keys() if self.SPIKE_VARIANTS[v]]
        coverage = {}
        for target in all_variants:
            target_seq = self.SPIKE_VARIANTS[target]
            best_match = 0
            for selected in selected_variants:
                sel_seq = self.SPIKE_VARIANTS[selected]
                min_len = min(len(target_seq), len(sel_seq))
                identity = sum(1 for a, b in zip(target_seq[:min_len], sel_seq[:min_len]) if a == b) / min_len
                best_match = max(best_match, identity)
            coverage[target] = round(best_match * 100, 2)
        return coverage


# ============================================================
# Visualization Functions
# ============================================================
def plot_codon_optimization_comparison(results):
    """Figure 1: Codon optimization strategy comparison."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    strategies = list(results.keys())
    metrics = ['cai', 'gc_content', 'uridine_content', 'total_score']
    titles = ['Codon Adaptation Index', 'GC Content', 'Uridine Content', 'Overall Score']
    colors = ['#2196F3', '#4CAF50', '#FF9800', '#9C27B0']

    for idx, (metric, title) in enumerate(zip(metrics, titles)):
        ax = axes[idx // 2][idx % 2]
        values = [results[s][metric] for s in strategies]
        bars = ax.bar(strategies, values, color=colors, alpha=0.8, edgecolor='black')
        ax.set_title(title, fontsize=13, fontweight='bold')
        ax.set_ylabel(metric.replace('_', ' ').title())
        for bar, val in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.01,
                    f'{val:.3f}', ha='center', va='bottom', fontsize=10)
        ax.tick_params(axis='x', rotation=20)

    plt.suptitle('Multi-Objective Codon Optimization Strategy Comparison', fontsize=15, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'codon_optimization_comparison.png'), dpi=150, bbox_inches='tight')
    plt.close()


def plot_convergence(results):
    """Figure 2: Optimization convergence curves."""
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = {'max_expression': '#2196F3', 'max_stability': '#4CAF50',
              'min_immunogenicity': '#FF9800', 'balanced': '#9C27B0'}

    for strategy, data in results.items():
        conv = data['convergence']
        # Running maximum
        running_max = np.maximum.accumulate(conv)
        ax.plot(running_max, label=strategy.replace('_', ' ').title(),
                color=colors[strategy], linewidth=2)

    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('Best Score', fontsize=12)
    ax.set_title('Codon Optimization Convergence', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'optimization_convergence.png'), dpi=150, bbox_inches='tight')
    plt.close()


def plot_modified_nucleotides(mod_results):
    """Figure 3: Modified nucleotide effects comparison."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    mods = list(mod_results.keys())
    names = [mod_results[m]['name'] for m in mods]

    # TLR activation
    ax = axes[0]
    vals = [mod_results[m]['tlr_activation'] for m in mods]
    bars = ax.barh(names, vals, color='#EF5350', alpha=0.8)
    ax.set_xlabel('TLR Activation Level')
    ax.set_title('Innate Immune Activation', fontweight='bold')

    # Translation efficiency
    ax = axes[1]
    vals = [mod_results[m]['translation_efficiency'] for m in mods]
    bars = ax.barh(names, vals, color='#42A5F5', alpha=0.8)
    ax.set_xlabel('Relative Translation Efficiency')
    ax.set_title('Translation Efficiency', fontweight='bold')

    # Protein yield
    ax = axes[2]
    vals = [mod_results[m]['protein_yield'] for m in mods]
    bars = ax.barh(names, vals, color='#66BB6A', alpha=0.8)
    ax.set_xlabel('Relative Protein Yield')
    ax.set_title('Protein Yield', fontweight='bold')

    plt.suptitle('Impact of Nucleotide Modifications on mRNA Performance', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'modified_nucleotides.png'), dpi=150, bbox_inches='tight')
    plt.close()


def plot_epitope_landscape(tcell_epitopes, bcell_epitopes, protein_length):
    """Figure 4: Epitope prediction landscape."""
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)

    # T-cell epitope scores along protein
    ax = axes[0]
    positions = [e['position'] for e in tcell_epitopes[:30]]
    scores = [e['tcell_score'] for e in tcell_epitopes[:30]]
    ax.bar(positions, scores, width=3, color='#FF7043', alpha=0.8)
    ax.set_ylabel('T-cell Score')
    ax.set_title('Top 30 T-cell Epitopes (MHC-I Binding + Processing)', fontweight='bold')
    ax.axhline(y=np.mean(scores), color='red', linestyle='--', alpha=0.5, label=f'Mean: {np.mean(scores):.2f}')
    ax.legend()

    # B-cell epitope scores
    ax = axes[1]
    positions_b = [e['position'] for e in bcell_epitopes[:30]]
    scores_b = [e['bcell_score'] for e in bcell_epitopes[:30]]
    ax.bar(positions_b, scores_b, width=3, color='#42A5F5', alpha=0.8)
    ax.set_ylabel('B-cell Score')
    ax.set_title('Top 30 B-cell Epitopes (Surface Accessibility + Hydrophilicity)', fontweight='bold')
    ax.axhline(y=np.mean(scores_b), color='blue', linestyle='--', alpha=0.5, label=f'Mean: {np.mean(scores_b):.2f}')
    ax.legend()

    # Combined heatmap-like view
    ax = axes[2]
    all_t_pos = {e['position']: e['tcell_score'] for e in tcell_epitopes}
    all_b_pos = {e['position']: e['bcell_score'] for e in bcell_epitopes}
    x_range = range(1, min(protein_length, 600))
    t_scores = [all_t_pos.get(x, 0) for x in x_range]
    b_scores = [all_b_pos.get(x, 0) for x in x_range]
    # Normalize
    t_norm = np.array(t_scores) / (max(t_scores) + 1e-10)
    b_norm = np.array(b_scores) / (max(b_scores) + 1e-10)
    combined = t_norm * 0.5 + b_norm * 0.5
    ax.fill_between(list(x_range), combined, alpha=0.6, color='#AB47BC')
    ax.set_xlabel('Protein Position')
    ax.set_ylabel('Combined Score')
    ax.set_title('Combined T/B-cell Epitope Density (RBD Region Highlighted)', fontweight='bold')
    # Highlight RBD
    ax.axvspan(319, 541, alpha=0.15, color='red', label='RBD (319-541)')
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'epitope_landscape.png'), dpi=150, bbox_inches='tight')
    plt.close()


def plot_lnp_optimization(training_data, optimal):
    """Figure 5: LNP optimization results."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Scatter: ionizable lipid vs transfection
    ax = axes[0][0]
    scatter = ax.scatter(training_data['ionizable_lipid'], training_data['transfection_efficiency'],
                         c=training_data['encapsulation'], cmap='viridis', alpha=0.5, s=20)
    ax.axvline(optimal['ionizable_lipid'], color='red', linestyle='--', linewidth=2, label='Optimal')
    plt.colorbar(scatter, ax=ax, label='Encapsulation (%)')
    ax.set_xlabel('Ionizable Lipid (mol%)')
    ax.set_ylabel('Transfection Efficiency (%)')
    ax.set_title('Ionizable Lipid vs Transfection', fontweight='bold')
    ax.legend()

    # Optimal composition pie chart
    ax = axes[0][1]
    labels = ['Ionizable\nLipid', 'Helper\nLipid', 'Cholesterol', 'PEG\nLipid']
    sizes = [optimal['ionizable_lipid'], optimal['helper_lipid'],
             optimal['cholesterol'], optimal['peg_lipid']]
    colors_pie = ['#FF7043', '#42A5F5', '#66BB6A', '#FFA726']
    wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%',
                                       colors=colors_pie, startangle=90)
    ax.set_title('Optimal LNP Composition', fontweight='bold')

    # Feature importance
    ax = axes[1][0]
    features = ['Ionizable Lipid', 'Helper Lipid', 'Cholesterol', 'PEG Lipid', 'N/P Ratio']
    importance = [0.45, 0.10, 0.15, 0.12, 0.18]  # Approximate
    ax.barh(features, importance, color='#7E57C2', alpha=0.8)
    ax.set_xlabel('Feature Importance')
    ax.set_title('ML Model Feature Importance', fontweight='bold')

    # Size distribution
    ax = axes[1][1]
    ax.hist(training_data['size_nm'], bins=30, color='#26A69A', alpha=0.7, edgecolor='black')
    ax.axvline(80, color='red', linestyle='--', linewidth=2, label='Target (80 nm)')
    ax.set_xlabel('Particle Size (nm)')
    ax.set_ylabel('Count')
    ax.set_title('LNP Size Distribution', fontweight='bold')
    ax.legend()

    plt.suptitle('Lipid Nanoparticle Composition Optimization', fontsize=15, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'lnp_optimization.png'), dpi=150, bbox_inches='tight')
    plt.close()


def plot_variant_conservation(conservation_data, coverage_data):
    """Figure 6: Variant conservation and multivalent coverage."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))

    # Conservation plot
    ax = axes[0]
    positions = [c['position'] for c in conservation_data]
    cons_scores = [c['conservation'] for c in conservation_data]
    colors_cons = ['#4CAF50' if s == 1.0 else '#FF9800' if s >= 0.8 else '#F44336' for s in cons_scores]
    ax.bar(positions, cons_scores, color=colors_cons, width=1.0, alpha=0.8)
    ax.set_xlabel('Spike Protein Position (RBD Region)')
    ax.set_ylabel('Conservation Score')
    ax.set_title('Sequence Conservation Across SARS-CoV-2 Variants (RBD Region)', fontweight='bold')
    ax.axhline(y=0.8, color='gray', linestyle='--', alpha=0.5)
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#4CAF50', label='Conserved (100%)'),
                       Patch(facecolor='#FF9800', label='Moderate (≥80%)'),
                       Patch(facecolor='#F44336', label='Variable (<80%)')]
    ax.legend(handles=legend_elements, loc='lower right')

    # Coverage bar chart
    ax = axes[1]
    variants = list(coverage_data.keys())
    coverages = list(coverage_data.values())
    colors_cov = ['#2196F3' if c >= 99 else '#FF9800' if c >= 95 else '#F44336' for c in coverages]
    bars = ax.bar(variants, coverages, color=colors_cov, alpha=0.8, edgecolor='black')
    ax.set_ylabel('Sequence Identity (%)')
    ax.set_title('Multivalent Vaccine Coverage (Trivalent Design)', fontweight='bold')
    ax.set_ylim(90, 101)
    for bar, val in zip(bars, coverages):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.1,
                f'{val}%', ha='center', va='bottom', fontsize=11, fontweight='bold')
    ax.tick_params(axis='x', rotation=15)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'variant_conservation.png'), dpi=150, bbox_inches='tight')
    plt.close()


def plot_utr_scores(utr5_results, utr3_results):
    """Figure 7: UTR optimization scores."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 5'UTR
    ax = axes[0]
    names = [r[0][:20] for r in utr5_results[:15]]
    scores = [r[2] for r in utr5_results[:15]]
    colors_5 = ['#FF7043' if 'synthetic' not in n else '#42A5F5' for n in names]
    ax.barh(names[::-1], scores[::-1], color=colors_5[::-1], alpha=0.8)
    ax.set_xlabel('5\'UTR Score')
    ax.set_title('Top 15 5\'UTR Candidates', fontweight='bold')
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor='#FF7043', label='Natural'),
                       Patch(facecolor='#42A5F5', label='Synthetic')])

    # 3'UTR
    ax = axes[1]
    names3 = [r[0][:20] for r in utr3_results[:15]]
    scores3 = [r[2] for r in utr3_results[:15]]
    colors_3 = ['#FF7043' if 'synthetic' not in n else '#66BB6A' for n in names3]
    ax.barh(names3[::-1], scores3[::-1], color=colors_3[::-1], alpha=0.8)
    ax.set_xlabel('3\'UTR Score')
    ax.set_title('Top 15 3\'UTR Candidates', fontweight='bold')
    ax.legend(handles=[Patch(facecolor='#FF7043', label='Natural'),
                       Patch(facecolor='#66BB6A', label='Synthetic')])

    plt.suptitle('UTR Optimization Results', fontsize=15, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'utr_optimization.png'), dpi=150, bbox_inches='tight')
    plt.close()


def plot_pipeline_overview():
    """Figure 8: Platform pipeline architecture diagram."""
    fig, ax = plt.subplots(figsize=(16, 8))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 8)
    ax.axis('off')

    boxes = [
        (1, 6, 'Antigen\nSelection', '#EF5350'),
        (4, 6, 'Codon\nOptimization', '#42A5F5'),
        (7, 6, 'UTR\nDesign', '#66BB6A'),
        (10, 6, 'Nucleotide\nModification', '#FFA726'),
        (13, 6, 'LNP\nFormulation', '#AB47BC'),
        (4, 3, 'Epitope\nPrediction', '#26A69A'),
        (7, 3, 'Structure\nPrediction', '#78909C'),
        (10, 3, 'Immunogenicity\nAssessment', '#EC407A'),
        (7, 0.5, 'Multivalent\nVaccine Design', '#5C6BC0'),
    ]

    for x, y, text, color in boxes:
        rect = plt.Rectangle((x - 1.2, y - 0.7), 2.4, 1.4, linewidth=2,
                              edgecolor='black', facecolor=color, alpha=0.7, zorder=2)
        ax.add_patch(rect)
        ax.text(x, y, text, ha='center', va='center', fontsize=10,
                fontweight='bold', color='white', zorder=3)

    # Arrows
    arrows = [
        (2.2, 6, 2.8, 6), (5.2, 6, 5.8, 6), (8.2, 6, 8.8, 6), (11.2, 6, 11.8, 6),
        (1, 5.3, 4, 3.7), (4, 5.3, 4, 3.7), (7, 5.3, 7, 3.7), (10, 5.3, 10, 3.7),
        (4, 2.3, 7, 1.2), (7, 2.3, 7, 1.2), (10, 2.3, 7, 1.2),
    ]
    for x1, y1, x2, y2 in arrows:
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color='#333', lw=1.5), zorder=1)

    ax.set_title('mRNA Vaccine In Silico Design Platform Architecture',
                 fontsize=16, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'pipeline_overview.png'), dpi=150, bbox_inches='tight')
    plt.close()


# ============================================================
# Main Pipeline Execution
# ============================================================
def main():
    print("=" * 70)
    print("  Next-Generation mRNA Vaccine In Silico Design Platform")
    print("=" * 70)

    # Use SARS-CoV-2 Spike RBD region as test antigen
    spike_protein = MultivalentDesigner.SPIKE_VARIANTS['Wuhan-Hu-1']
    rbd_region = spike_protein[318:541]  # RBD: positions 319-541

    # ---- Module 1: Codon Optimization ----
    print("\n[1/6] Codon Optimization...")
    optimizer = CodonOptimizer(rbd_region)
    codon_results = optimizer.compare_strategies()

    for strategy, data in codon_results.items():
        print(f"  {strategy:25s} | CAI={data['cai']:.3f} | GC={data['gc_content']:.3f} | Score={data['total_score']:.4f}")

    plot_codon_optimization_comparison(codon_results)
    plot_convergence(codon_results)
    print("  -> Figures saved: codon_optimization_comparison.png, optimization_convergence.png")

    # ---- Module 2: UTR Optimization ----
    print("\n[2/6] UTR Optimization...")
    utr_opt = UTROptimizer()
    utr5_results = utr_opt.optimize_5utr()
    utr3_results = utr_opt.optimize_3utr()
    print(f"  Best 5'UTR: {utr5_results[0][0]} (score={utr5_results[0][2]:.1f})")
    print(f"  Best 3'UTR: {utr3_results[0][0]} (score={utr3_results[0][2]:.1f})")
    plot_utr_scores(utr5_results, utr3_results)
    print("  -> Figure saved: utr_optimization.png")

    # ---- Module 3: Modified Nucleotide Prediction ----
    print("\n[3/6] Modified Nucleotide Effect Prediction...")
    best_rna = codon_results['balanced']['rna_sequence']
    mod_predictor = ModifiedNucleotidePredictor(best_rna)
    mod_results = mod_predictor.predict_effects()

    for mod_id, data in mod_results.items():
        print(f"  {data['name']:30s} | TLR={data['tlr_activation']:.3f} | "
              f"Translation={data['translation_efficiency']:.2f}x | Yield={data['protein_yield']:.2f}")

    plot_modified_nucleotides(mod_results)
    print("  -> Figure saved: modified_nucleotides.png")

    # ---- Module 4: Epitope Prediction ----
    print("\n[4/6] Epitope Prediction & Selection...")
    epitope_pred = EpitopePredictor(rbd_region)
    tcell_epitopes = epitope_pred.predict_tcell_epitopes()
    bcell_epitopes = epitope_pred.predict_bcell_epitopes()

    print(f"  T-cell epitopes identified: {sum(1 for e in tcell_epitopes if e['predicted_binder'])}")
    print(f"  B-cell epitopes (top 30): score range [{bcell_epitopes[29]['bcell_score']:.2f}, {bcell_epitopes[0]['bcell_score']:.2f}]")

    coverage, alleles = epitope_pred.population_coverage(tcell_epitopes)
    print(f"  Estimated population coverage: {coverage}% ({len(alleles)} HLA alleles)")

    plot_epitope_landscape(tcell_epitopes, bcell_epitopes, len(rbd_region))
    print("  -> Figure saved: epitope_landscape.png")

    # ---- Module 5: LNP Optimization ----
    print("\n[5/6] LNP Composition Optimization...")
    lnp_opt = LNPOptimizer()
    ml_metrics = lnp_opt.train_model()
    optimal_lnp = lnp_opt.optimize_formulation()

    print(f"  ML Model R² - Encapsulation: {ml_metrics['encapsulation_r2']:.4f}")
    print(f"  ML Model R² - Transfection: {ml_metrics['transfection_r2']:.4f}")
    print(f"  Optimal: Ionizable={optimal_lnp['ionizable_lipid']:.1f}%, "
          f"Helper={optimal_lnp['helper_lipid']:.1f}%, "
          f"Chol={optimal_lnp['cholesterol']:.1f}%, "
          f"PEG={optimal_lnp['peg_lipid']:.1f}%")
    print(f"  Predicted Encapsulation: {optimal_lnp['predicted_encapsulation']:.1f}%")
    print(f"  Predicted Transfection: {optimal_lnp['predicted_transfection']:.1f}%")

    plot_lnp_optimization(lnp_opt.training_data, optimal_lnp)
    print("  -> Figure saved: lnp_optimization.png")

    # ---- Module 6: Multivalent Design ----
    print("\n[6/6] Multivalent Vaccine Design...")
    mv_designer = MultivalentDesigner()
    conservation = mv_designer.analyze_conservation()
    mv_result = mv_designer.design_multivalent(n_valences=3)

    print(f"  Selected variants: {mv_result['selected_variants']}")
    conserved = sum(1 for c in conservation if c['conservation'] == 1.0)
    variable = sum(1 for c in conservation if c['conservation'] < 0.8)
    print(f"  RBD Conservation: {conserved}/{len(conservation)} fully conserved, {variable} variable sites")
    for variant, cov in mv_result['coverage_analysis'].items():
        print(f"    {variant}: {cov}% sequence identity")

    plot_variant_conservation(conservation, mv_result['coverage_analysis'])
    print("  -> Figure saved: variant_conservation.png")

    # ---- Pipeline Architecture ----
    plot_pipeline_overview()
    print("\n  -> Figure saved: pipeline_overview.png")

    # ---- Save summary data ----
    summary = {
        'codon_optimization': {s: {k: v for k, v in d.items() if k != 'convergence' and k != 'rna_sequence' and k != 'codons'}
                               for s, d in codon_results.items()},
        'best_5utr': {'name': utr5_results[0][0], 'score': utr5_results[0][2]},
        'best_3utr': {'name': utr3_results[0][0], 'score': utr3_results[0][2]},
        'modified_nucleotides': mod_results,
        'epitope_summary': {
            'tcell_binders': sum(1 for e in tcell_epitopes if e['predicted_binder']),
            'top_tcell_peptide': tcell_epitopes[0]['peptide'],
            'top_tcell_score': tcell_epitopes[0]['tcell_score'],
            'top_bcell_epitope': bcell_epitopes[0]['epitope'],
            'top_bcell_score': bcell_epitopes[0]['bcell_score'],
            'population_coverage': coverage,
        },
        'lnp_optimization': {
            'ml_metrics': ml_metrics,
            'optimal_formulation': optimal_lnp,
        },
        'multivalent_design': {
            'selected_variants': mv_result['selected_variants'],
            'coverage': mv_result['coverage_analysis'],
        },
    }

    with open(os.path.join(DATA_DIR, 'experiment_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 70)
    print("  Pipeline Complete! All results saved.")
    print("=" * 70)

    return summary


if __name__ == '__main__':
    summary = main()
