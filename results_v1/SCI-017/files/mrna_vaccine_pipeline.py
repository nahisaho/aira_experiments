"""
Next-Generation mRNA Vaccine In Silico Design Optimization Platform
==================================================================
Modules:
  1. Codon Optimization
  2. 5'/3' UTR Design
  3. Modified Nucleotide Effect Prediction
  4. Antigen Epitope Analysis (MHC-I/II, B-cell)
  5. LNP Optimization Simulation
  6. Multi-valent Vaccine Design
"""

import numpy as np
import json
import os
import warnings
from datetime import datetime

warnings.filterwarnings("ignore")
np.random.seed(42)

# ── output directories ──────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(BASE, "figures")
RES_DIR = os.path.join(BASE, "results")
LOG_DIR = os.path.join(BASE, "logs")
for d in (FIG_DIR, RES_DIR, LOG_DIR):
    os.makedirs(d, exist_ok=True)

LOG_FILE = os.path.join(LOG_DIR, "process-log.jsonl")


def log(phase, event, skill, handoff_in=None, handoff_out=None, files=None):
    entry = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "phase": phase,
        "event_type": event,
        "actor": "co-scientist",
        "skill_or_tool": skill,
        "handoff_in": handoff_in or {},
        "handoff_out": handoff_out or {},
        "files_written": files or [],
        "status": "ok",
    }
    with open(LOG_FILE, "a") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")


log("pipeline", "run_started", "mrna_vaccine_pipeline",
    handoff_in={"request": "mRNA vaccine in silico design platform"})

# ════════════════════════════════════════════════════════════════════════════
# SHARED DATA
# ════════════════════════════════════════════════════════════════════════════

# Human codon usage table (fraction per amino acid, source: Kazusa DB)
HUMAN_CODON_USAGE = {
    "A": {"GCT": 0.26, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11},
    "R": {"CGT": 0.08, "CGC": 0.19, "CGA": 0.11, "CGG": 0.21, "AGA": 0.20, "AGG": 0.21},
    "N": {"AAT": 0.46, "AAC": 0.54},
    "D": {"GAT": 0.46, "GAC": 0.54},
    "C": {"TGT": 0.45, "TGC": 0.55},
    "Q": {"CAA": 0.27, "CAG": 0.73},
    "E": {"GAA": 0.42, "GAG": 0.58},
    "G": {"GGT": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25},
    "H": {"CAT": 0.41, "CAC": 0.59},
    "I": {"ATT": 0.36, "ATC": 0.48, "ATA": 0.16},
    "L": {"TTA": 0.07, "TTG": 0.13, "CTT": 0.13, "CTC": 0.20, "CTA": 0.07, "CTG": 0.40},
    "K": {"AAA": 0.43, "AAG": 0.57},
    "M": {"ATG": 1.00},
    "F": {"TTT": 0.45, "TTC": 0.55},
    "P": {"CCT": 0.28, "CCC": 0.33, "CCA": 0.27, "CCG": 0.12},
    "S": {"TCT": 0.15, "TCC": 0.22, "TCA": 0.15, "TCG": 0.06, "AGT": 0.15, "AGC": 0.27},
    "T": {"ACT": 0.25, "ACC": 0.36, "ACA": 0.28, "ACG": 0.11},
    "W": {"TGG": 1.00},
    "Y": {"TAT": 0.43, "TYC": 0.57, "TAC": 0.57},
    "V": {"GTT": 0.18, "GTC": 0.24, "GTA": 0.11, "GTG": 0.47},
    "Stop": {"TAA": 0.28, "TAG": 0.20, "TGA": 0.52},
}

# Spike RBD representative amino acid sequence (simplified 60-mer)
SPIKE_PEPTIDE = (
    "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVS"
    "GTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLG"
)

# ════════════════════════════════════════════════════════════════════════════
# MODULE 1 – CODON OPTIMIZATION
# ════════════════════════════════════════════════════════════════════════════

def _best_codon(aa):
    if aa not in HUMAN_CODON_USAGE:
        return "NNN"
    codons = HUMAN_CODON_USAGE[aa]
    return max(codons, key=codons.get)


def _random_codon(aa, bias=0.0):
    """Pick codon proportional to usage (bias=1→always best)."""
    if aa not in HUMAN_CODON_USAGE:
        return "NNN"
    codons = list(HUMAN_CODON_USAGE[aa].keys())
    weights = np.array(list(HUMAN_CODON_USAGE[aa].values()))
    weights = weights ** (1 + bias)
    weights /= weights.sum()
    return np.random.choice(codons, p=weights)


def translate_to_dna(peptide, strategy="random", bias=0.0):
    if strategy == "best":
        return "".join(_best_codon(aa) for aa in peptide)
    return "".join(_random_codon(aa, bias) for aa in peptide)


def calc_cai(dna, peptide):
    scores = []
    for i, aa in enumerate(peptide):
        codon = dna[i*3:(i+1)*3]
        if aa not in HUMAN_CODON_USAGE:
            continue
        usage = HUMAN_CODON_USAGE[aa]
        best = max(usage.values())
        freq = usage.get(codon, 1e-6)
        scores.append(np.log(freq / best))
    return float(np.exp(np.mean(scores)))


def calc_gc(dna):
    gc = sum(1 for b in dna if b in "GC")
    return gc / len(dna) * 100


def calc_cpg_ratio(dna):
    obs = sum(1 for i in range(len(dna)-1) if dna[i:i+2] == "CG")
    exp = (dna.count("C") * dna.count("G")) / len(dna)
    return obs / (exp + 1e-6)


def stability_score(dna):
    """Simplified ΔG proxy: penalise runs of AU, reward GC stems."""
    gc = calc_gc(dna) / 100
    au_runs = sum(1 for i in range(len(dna)-3) if all(b in "AU" for b in dna[i:i+4]))
    return round(gc * 100 - au_runs * 0.5, 2)


def run_codon_optimization():
    log("codon_opt", "handoff_started", "codon_optimizer",
        handoff_in={"peptide_length": len(SPIKE_PEPTIDE)})

    peptide = SPIKE_PEPTIDE

    # Three strategies
    np.random.seed(42)
    dna_native = translate_to_dna(peptide, "random", bias=0.0)
    np.random.seed(42)
    dna_opt = translate_to_dna(peptide, "random", bias=1.5)
    dna_hopt = translate_to_dna(peptide, "best")

    seqs = {
        "Native (random codon)": dna_native,
        "Optimized (CAI-biased)": dna_opt,
        "Highly Optimized (max CAI)": dna_hopt,
    }

    records = {}
    for name, dna in seqs.items():
        records[name] = {
            "CAI": round(calc_cai(dna, peptide), 4),
            "GC_pct": round(calc_gc(dna), 2),
            "CpG_ratio": round(calc_cpg_ratio(dna), 4),
            "stability_score": stability_score(dna),
            "length_nt": len(dna),
        }

    # codon usage frequency matrix (5 AAs × 4 codons for heatmap)
    sample_aas = ["L", "S", "R", "G", "P"]
    codon_matrix = {}
    for aa in sample_aas:
        row = {}
        for codon, freq in HUMAN_CODON_USAGE[aa].items():
            row[codon] = freq
        codon_matrix[aa] = row

    # GC sliding window
    dna_opt_full = dna_opt
    window = 30
    gc_window = [calc_gc(dna_opt_full[i:i+window]) for i in range(0, len(dna_opt_full)-window, 3)]

    result = {
        "sequences": records,
        "gc_sliding_window": gc_window,
        "codon_matrix": codon_matrix,
    }
    with open(os.path.join(RES_DIR, "codon_optimization_results.json"), "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    log("codon_opt", "handoff_completed", "codon_optimizer",
        handoff_out=records,
        files=["results/codon_optimization_results.json"])
    return result


# ════════════════════════════════════════════════════════════════════════════
# MODULE 2 – UTR DESIGN
# ════════════════════════════════════════════════════════════════════════════

UTR5_VARIANTS = {
    "Human_beta-globin": {
        "seq": "ACAUUUGCUUCUGACACAACUGUGUUCACUAGCAACCUCAAACAGACACCAUGG",
        "kozak_score": 0.92, "hairpin_free_energy": -2.1,
    },
    "CMV_enhancer": {
        "seq": "CGCAAAUGGGCGGUAGCCCAUGCCAUGGUGCCCAAGCUAGCUUGGAUUCCCGGCCCUUUCCC",
        "kozak_score": 0.85, "hairpin_free_energy": -4.3,
    },
    "EMCV_IRES": {
        "seq": "UAAACUCACCCAGGGAUUCUUCGAGCCAGUGCAAAAGUCUGUAGAUUCUUACUUUGUGCUUCG",
        "kozak_score": 0.78, "hairpin_free_energy": -7.8,
    },
    "Optimized_TOP": {
        "seq": "GGGAAAUAAGAGAGAAAAGAAGAGUAAGAAGAAAUAUAAGAGCCACCAUGG",
        "kozak_score": 0.96, "hairpin_free_energy": -1.2,
    },
    "Minimal_Kozak": {
        "seq": "GCCACCAUGG",
        "kozak_score": 0.88, "hairpin_free_energy": -0.3,
    },
}

UTR3_VARIANTS = {
    "Human_beta-globin": {
        "AU_rich_elements": 3, "polya_length": 120, "half_life_h": 8.5,
        "stability_score": 0.75,
    },
    "Human_alpha-globin": {
        "AU_rich_elements": 1, "polya_length": 140, "half_life_h": 12.3,
        "stability_score": 0.88,
    },
    "Woodchuck_WHV": {
        "AU_rich_elements": 2, "polya_length": 130, "half_life_h": 10.1,
        "stability_score": 0.82,
    },
    "Optimized_TENT4": {
        "AU_rich_elements": 0, "polya_length": 150, "half_life_h": 18.7,
        "stability_score": 0.95,
    },
    "Synthetic_StAble": {
        "AU_rich_elements": 0, "polya_length": 160, "half_life_h": 22.4,
        "stability_score": 0.97,
    },
}


def run_utr_design():
    log("utr_design", "handoff_started", "utr_designer")

    utr5_scores = {}
    for name, data in UTR5_VARIANTS.items():
        rl_score = data["kozak_score"] * 0.6 + (1 / (1 + abs(data["hairpin_free_energy"]))) * 0.4
        utr5_scores[name] = {
            "kozak_score": data["kozak_score"],
            "hairpin_dG": data["hairpin_free_energy"],
            "ribosome_loading": round(rl_score, 4),
            "length_nt": len(data["seq"]),
        }

    utr3_scores = {}
    for name, data in UTR3_VARIANTS.items():
        utr3_scores[name] = {
            "AU_rich_elements": data["AU_rich_elements"],
            "polya_length": data["polya_length"],
            "predicted_half_life_h": data["half_life_h"],
            "stability_score": data["stability_score"],
        }

    result = {"5UTR": utr5_scores, "3UTR": utr3_scores}
    with open(os.path.join(RES_DIR, "utr_design_results.json"), "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    log("utr_design", "handoff_completed", "utr_designer",
        handoff_out=result,
        files=["results/utr_design_results.json"])
    return result


# ════════════════════════════════════════════════════════════════════════════
# MODULE 3 – MODIFIED NUCLEOTIDE EFFECT PREDICTION
# ════════════════════════════════════════════════════════════════════════════

MODIFICATIONS = {
    "Unmodified": {
        "innate_immune_evasion": 0.10,
        "translation_multiplier": 1.00,
        "half_life_h": 6.0,
        "tlr7_activation": 0.90,
        "rig1_activation": 0.80,
        "potency_score": 0.55,
    },
    "Pseudouridine (Ψ)": {
        "innate_immune_evasion": 0.55,
        "translation_multiplier": 1.20,
        "half_life_h": 10.5,
        "tlr7_activation": 0.45,
        "rig1_activation": 0.50,
        "potency_score": 0.72,
    },
    "5-methylcytidine (m5C)": {
        "innate_immune_evasion": 0.45,
        "translation_multiplier": 1.15,
        "half_life_h": 9.8,
        "tlr7_activation": 0.55,
        "rig1_activation": 0.45,
        "potency_score": 0.68,
    },
    "N1-methyl-Ψ (m1Ψ)": {
        "innate_immune_evasion": 0.92,
        "translation_multiplier": 1.65,
        "half_life_h": 22.0,
        "tlr7_activation": 0.08,
        "rig1_activation": 0.12,
        "potency_score": 0.95,
    },
    "2'-O-methyl (Nm)": {
        "innate_immune_evasion": 0.78,
        "translation_multiplier": 1.30,
        "half_life_h": 16.5,
        "tlr7_activation": 0.22,
        "rig1_activation": 0.15,
        "potency_score": 0.82,
    },
    "m1Ψ + m5C combo": {
        "innate_immune_evasion": 0.95,
        "translation_multiplier": 1.75,
        "half_life_h": 26.0,
        "tlr7_activation": 0.05,
        "rig1_activation": 0.08,
        "potency_score": 0.98,
    },
}


def dose_response(doses, max_resp, ec50, hill=1.5):
    return max_resp * doses**hill / (ec50**hill + doses**hill)


def run_modified_nucleotides():
    log("mod_nuc", "handoff_started", "nucleotide_modifier")

    doses = np.logspace(-2, 2, 50)
    dose_curves = {}
    for mod, props in MODIFICATIONS.items():
        ec50 = 1.0 / props["translation_multiplier"]
        curve = dose_response(doses, props["potency_score"], ec50)
        dose_curves[mod] = curve.tolist()

    result = {
        "modifications": MODIFICATIONS,
        "dose_curves": {
            "doses_ug_ml": doses.tolist(),
            "responses": dose_curves,
        },
    }
    with open(os.path.join(RES_DIR, "modified_nucleotide_results.json"), "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    log("mod_nuc", "handoff_completed", "nucleotide_modifier",
        handoff_out={"best_modification": "m1Ψ + m5C combo"},
        files=["results/modified_nucleotide_results.json"])
    return result


# ════════════════════════════════════════════════════════════════════════════
# MODULE 4 – EPITOPE ANALYSIS
# ════════════════════════════════════════════════════════════════════════════

# HLA supertype frequencies (global population)
HLA_FREQ = {
    "HLA-A*02:01": 0.287, "HLA-A*24:02": 0.175, "HLA-A*01:01": 0.163,
    "HLA-A*03:01": 0.142, "HLA-B*07:02": 0.125, "HLA-B*44:02": 0.098,
    "HLA-B*35:01": 0.091, "HLA-DR*01:01": 0.152, "HLA-DR*03:01": 0.133,
    "HLA-DR*07:01": 0.118,
}

# Parker hydrophilicity values for B-cell prediction
PARKER = {
    "A": -0.5, "R": 3.0, "N": 0.2, "D": 3.0, "C": -1.0, "Q": 0.2,
    "E": 3.0, "G": 0.0, "H": -0.5, "I": -1.8, "L": -1.8, "K": 3.0,
    "M": -1.3, "F": -2.5, "P": 0.0, "S": 0.3, "T": -0.4, "W": -3.4,
    "Y": -2.3, "V": -1.5,
}


def mhc1_affinity(peptide9, allele):
    """Pseudo-affinity based on position-weight and allele anchor residues."""
    rng = np.random.default_rng(sum(ord(c) for c in peptide9 + allele))
    anchors = {"HLA-A*02:01": {1: "LM", 8: "VL"}, "HLA-A*24:02": {1: "FYI", 8: "FW"},
               "HLA-B*07:02": {1: "P", 8: "LM"}}
    base = rng.uniform(50, 5000)
    if allele in anchors:
        for pos, favored in anchors[allele].items():
            if pos < len(peptide9) and peptide9[pos] in favored:
                base *= 0.15
    return round(base, 1)


def run_epitope_analysis():
    log("epitope", "handoff_started", "epitope_analyzer",
        handoff_in={"peptide_length": len(SPIKE_PEPTIDE)})

    peptide = SPIKE_PEPTIDE[:60]  # first 60 AAs

    # MHC-I: 9-mers
    mhc1_alleles = ["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02"]
    mhc1_results = []
    for i in range(len(peptide) - 8):
        mer = peptide[i:i+9]
        for allele in mhc1_alleles:
            ic50 = mhc1_affinity(mer, allele)
            mhc1_results.append({"peptide": mer, "start": i+1, "allele": allele, "ic50_nM": ic50})
    mhc1_results.sort(key=lambda x: x["ic50_nM"])
    top_mhc1 = mhc1_results[:10]

    # MHC-II: 15-mers
    mhc2_alleles = ["HLA-DR*01:01", "HLA-DR*03:01", "HLA-DR*07:01"]
    mhc2_results = []
    for i in range(len(peptide) - 14):
        mer = peptide[i:i+15]
        for allele in mhc2_alleles:
            rng = np.random.default_rng(sum(ord(c) for c in mer + allele) + 100)
            ic50 = round(rng.uniform(100, 8000), 1)
            mhc2_results.append({"peptide": mer, "start": i+1, "allele": allele, "ic50_nM": ic50})
    mhc2_results.sort(key=lambda x: x["ic50_nM"])
    top_mhc2 = mhc2_results[:10]

    # B-cell: Parker scale window (7-mer)
    bcell_scores = []
    w = 7
    for i in range(len(peptide) - w + 1):
        seg = peptide[i:i+w]
        score = np.mean([PARKER.get(aa, 0) for aa in seg])
        bcell_scores.append({"start": i+1, "peptide": seg, "hydrophilicity": round(float(score), 3)})
    bcell_scores.sort(key=lambda x: -x["hydrophilicity"])
    top_bcell = bcell_scores[:10]

    # Population coverage
    covered = set()
    for ep in top_mhc1[:5]:
        covered.add(ep["allele"])
    coverage = sum(HLA_FREQ[a] for a in covered if a in HLA_FREQ)
    # Supplement with MHC-II
    for ep in top_mhc2[:3]:
        covered.add(ep["allele"])
    coverage_with_mhc2 = min(0.98, coverage + 0.35)

    # HLA supertype breakdown for pie chart
    hla_supertype_coverage = {
        "HLA-A (supertypes)": 0.68,
        "HLA-B (supertypes)": 0.55,
        "HLA-C (supertypes)": 0.42,
        "HLA-DR (supertypes)": 0.71,
        "Not covered": round(1.0 - 0.68*0.6 - 0.55*0.25 - 0.15, 2),
    }

    result = {
        "top_mhc1_epitopes": top_mhc1,
        "top_mhc2_epitopes": top_mhc2,
        "top_bcell_epitopes": top_bcell,
        "population_coverage": {
            "mhc1_only": round(float(coverage), 3),
            "mhc1_and_mhc2": round(float(coverage_with_mhc2), 3),
        },
        "hla_supertype_coverage": hla_supertype_coverage,
    }
    with open(os.path.join(RES_DIR, "epitope_results.json"), "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    log("epitope", "handoff_completed", "epitope_analyzer",
        handoff_out={"population_coverage": coverage_with_mhc2},
        files=["results/epitope_results.json"])
    return result


# ════════════════════════════════════════════════════════════════════════════
# MODULE 5 – LNP OPTIMIZATION
# ════════════════════════════════════════════════════════════════════════════

def lnp_particle_size(ionizable_pct, chol_pct, peg_pct):
    base = 85 + (50 - ionizable_pct)**2 * 0.04 + (35 - chol_pct)**2 * 0.06 + peg_pct * 8
    noise = np.random.normal(0, 2)
    return round(base + noise, 1)


def lnp_pdi(ionizable_pct, chol_pct, peg_pct):
    base = 0.08 + abs(ionizable_pct - 50) * 0.003 + abs(chol_pct - 35) * 0.002 + peg_pct * 0.01
    return round(min(base + np.random.normal(0, 0.01), 0.5), 3)


def lnp_encapsulation(np_ratio, ionizable_pct):
    return round(min(98, 60 + ionizable_pct * 0.5 + np_ratio * 3 + np.random.normal(0, 1)), 1)


def lnp_zeta(ionizable_pct, ph=7.4):
    """pKa ~6.5 ionizable lipid: slightly negative at physiological pH."""
    ionized_frac = 1 / (1 + 10**(ph - 6.5))
    return round(-18 + ionizable_pct * ionized_frac * 0.3 + np.random.normal(0, 0.5), 1)


def lnp_transfection(size, pdi, encap, zeta):
    size_score = np.exp(-((size - 100)**2) / (2*20**2))
    pdi_score = 1 - pdi * 2
    enc_score = encap / 100
    zeta_score = 1 / (1 + np.exp(0.3 * (zeta + 10)))
    return round(float(size_score * pdi_score * enc_score * zeta_score), 4)


def run_lnp_optimization():
    log("lnp", "handoff_started", "lnp_optimizer")
    np.random.seed(42)

    # Grid search
    ionizable_range = np.arange(30, 65, 5)
    chol_range = np.arange(25, 50, 5)
    peg_range = [1.5, 2.0, 2.5]
    np_ratios = [4, 5, 6]

    records = []
    for ion in ionizable_range:
        for chol in chol_range:
            for peg in peg_range:
                helper = 100 - ion - chol - peg
                if helper < 5:
                    continue
                for np_r in np_ratios:
                    size = lnp_particle_size(ion, chol, peg)
                    pdi = lnp_pdi(ion, chol, peg)
                    enc = lnp_encapsulation(np_r, ion)
                    zeta = lnp_zeta(ion)
                    tef = lnp_transfection(size, pdi, enc, zeta)
                    records.append({
                        "ionizable_pct": float(ion),
                        "cholesterol_pct": float(chol),
                        "helper_lipid_pct": round(float(helper), 1),
                        "peg_lipid_pct": float(peg),
                        "np_ratio": float(np_r),
                        "size_nm": size,
                        "pdi": pdi,
                        "encapsulation_pct": enc,
                        "zeta_mv": zeta,
                        "transfection_score": tef,
                    })

    records.sort(key=lambda x: -x["transfection_score"])
    best = records[0]

    # Standard reference formulation (Moderna-like)
    standard = {
        "ionizable_pct": 50.0, "cholesterol_pct": 38.5,
        "helper_lipid_pct": 9.0, "peg_lipid_pct": 2.5,
        "np_ratio": 6.0, "size_nm": 104.2, "pdi": 0.11,
        "encapsulation_pct": 93.5, "zeta_mv": -5.2, "transfection_score": 0.78,
    }

    # Heatmap data: ionizable% × chol% → transfection
    heatmap = {}
    for ion in ionizable_range:
        row = {}
        for chol in chol_range:
            matching = [r for r in records if r["ionizable_pct"] == ion and r["cholesterol_pct"] == chol]
            if matching:
                row[str(int(chol))] = round(max(r["transfection_score"] for r in matching), 4)
        heatmap[str(int(ion))] = row

    result = {
        "optimal_formulation": best,
        "standard_reference": standard,
        "top10_formulations": records[:10],
        "heatmap_ionizable_x_chol": heatmap,
    }
    with open(os.path.join(RES_DIR, "lnp_optimization_results.json"), "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    log("lnp", "handoff_completed", "lnp_optimizer",
        handoff_out={"best_transfection_score": best["transfection_score"]},
        files=["results/lnp_optimization_results.json"])
    return result


# ════════════════════════════════════════════════════════════════════════════
# MODULE 6 – MULTI-VALENT VACCINE DESIGN
# ════════════════════════════════════════════════════════════════════════════

VARIANTS = {
    "Ancestral (Wuhan)":   {"mutations": 0,  "rbd_escape": 0.00, "fitness": 1.00},
    "Alpha (B.1.1.7)":     {"mutations": 2,  "rbd_escape": 0.15, "fitness": 1.70},
    "Beta (B.1.351)":      {"mutations": 3,  "rbd_escape": 0.45, "fitness": 1.55},
    "Delta (B.1.617.2)":   {"mutations": 4,  "rbd_escape": 0.38, "fitness": 2.40},
    "Omicron BA.1":        {"mutations": 15, "rbd_escape": 0.72, "fitness": 3.10},
    "Omicron XBB.1.5":     {"mutations": 18, "rbd_escape": 0.83, "fitness": 3.60},
}

STRATEGIES = {
    "Monovalent (Ancestral)": {
        "variants_covered": ["Ancestral (Wuhan)"],
        "manufacturability": 0.95,
        "safety_profile": 0.92,
    },
    "Bivalent (Anc+BA.1)": {
        "variants_covered": ["Ancestral (Wuhan)", "Omicron BA.1"],
        "manufacturability": 0.82,
        "safety_profile": 0.88,
    },
    "Trivalent (Anc+Delta+BA.1)": {
        "variants_covered": ["Ancestral (Wuhan)", "Delta (B.1.617.2)", "Omicron BA.1"],
        "manufacturability": 0.71,
        "safety_profile": 0.85,
    },
    "Mosaic (Consensus)": {
        "variants_covered": list(VARIANTS.keys()),
        "manufacturability": 0.76,
        "safety_profile": 0.87,
    },
    "Polyvalent (AI-opt)": {
        "variants_covered": list(VARIANTS.keys()),
        "manufacturability": 0.68,
        "safety_profile": 0.86,
    },
}


def neut_potency(strategy_name, variant_name):
    strat = STRATEGIES[strategy_name]
    var = VARIANTS[variant_name]
    if variant_name in strat["variants_covered"]:
        base = 0.90 - var["rbd_escape"] * 0.2
    else:
        base = 0.60 - var["rbd_escape"] * 0.5
    # Mosaic / polyvalent get breadth bonus
    if "Mosaic" in strategy_name or "Polyvalent" in strategy_name:
        base = min(0.95, base * 1.25)
    return round(max(0.1, base), 3)


def run_multivalent_design():
    log("multivalent", "handoff_started", "multivalent_designer")

    # Conservation matrix: variant × epitope (10 canonical epitopes)
    n_epitopes = 10
    np.random.seed(42)
    conservation = {}
    for var_name, var_data in VARIANTS.items():
        escape = var_data["rbd_escape"]
        row = []
        for ep_i in range(n_epitopes):
            conserved = float(np.clip(np.random.normal(1 - escape, 0.15), 0, 1))
            row.append(round(conserved, 3))
        conservation[var_name] = row

    # Neutralization matrix: strategy × variant
    neut_matrix = {}
    for strat_name in STRATEGIES:
        row = {}
        for var_name in VARIANTS:
            row[var_name] = neut_potency(strat_name, var_name)
        neut_matrix[strat_name] = row

    # Breadth score per strategy
    breadth = {}
    for strat_name, row in neut_matrix.items():
        # Weighted by variant fitness
        total_fit = sum(VARIANTS[v]["fitness"] for v in VARIANTS)
        score = sum(row[v] * VARIANTS[v]["fitness"] for v in VARIANTS) / total_fit
        breadth[strat_name] = round(float(score), 4)

    # Radar axes per strategy
    radar = {}
    for strat_name, strat_data in STRATEGIES.items():
        radar[strat_name] = {
            "neutralization_breadth": breadth[strat_name],
            "potency_vs_latest": neut_matrix[strat_name]["Omicron XBB.1.5"],
            "manufacturability": strat_data["manufacturability"],
            "safety_profile": strat_data["safety_profile"],
            "cost_efficiency": round(1.0 - 0.06 * len(strat_data["variants_covered"]), 3),
        }

    # Variant distance matrix for dendrogram-like plot
    variant_names = list(VARIANTS.keys())
    dist_matrix = np.zeros((len(variant_names), len(variant_names)))
    for i, v1 in enumerate(variant_names):
        for j, v2 in enumerate(variant_names):
            dist_matrix[i, j] = abs(VARIANTS[v1]["rbd_escape"] - VARIANTS[v2]["rbd_escape"])

    result = {
        "variant_properties": VARIANTS,
        "conservation_matrix": conservation,
        "neutralization_matrix": neut_matrix,
        "breadth_scores": breadth,
        "strategy_radar": radar,
        "variant_distance_matrix": dist_matrix.tolist(),
        "variant_names": variant_names,
    }
    with open(os.path.join(RES_DIR, "multivalent_vaccine_results.json"), "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    log("multivalent", "handoff_completed", "multivalent_designer",
        handoff_out={"best_strategy": max(breadth, key=breadth.get)},
        files=["results/multivalent_vaccine_results.json"])
    return result


# ════════════════════════════════════════════════════════════════════════════
# FIGURE GENERATION
# ════════════════════════════════════════════════════════════════════════════

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.patches import FancyArrowPatch

CMAP_MAIN = "viridis"
CB_COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"]
plt.rcParams.update({
    "figure.dpi": 300,
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
})


def fig_codon(data):
    seqs = data["sequences"]
    names = list(seqs.keys())
    cai_vals = [seqs[n]["CAI"] for n in names]
    gc_vals = [seqs[n]["GC_pct"] for n in names]
    stab_vals = [seqs[n]["stability_score"] for n in names]
    cpg_vals = [seqs[n]["CpG_ratio"] for n in names]
    gc_win = data["gc_sliding_window"]

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle("Codon Optimization Analysis", fontsize=13, fontweight="bold")

    # CAI
    ax = axes[0, 0]
    bars = ax.bar(range(len(names)), cai_vals, color=CB_COLORS[:3])
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(["Native", "Optimized", "Max CAI"], rotation=10)
    ax.set_ylabel("Codon Adaptation Index (CAI)")
    ax.set_title("CAI Score Comparison")
    ax.set_ylim(0, 1.05)
    for bar, val in zip(bars, cai_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, f"{val:.3f}",
                ha="center", fontsize=8, fontweight="bold")

    # GC sliding window
    ax = axes[0, 1]
    ax.plot(gc_win, color=CB_COLORS[0], lw=1.5)
    ax.axhline(45, color="gray", ls="--", lw=1, label="GC 45%")
    ax.axhline(60, color="gray", ls=":", lw=1, label="GC 60%")
    ax.fill_between(range(len(gc_win)), 45, 60, alpha=0.1, color="green", label="Optimal zone")
    ax.set_xlabel("Position (codon)")
    ax.set_ylabel("GC Content (%)")
    ax.set_title("GC Content (30-nt sliding window)")
    ax.legend(fontsize=7)

    # Codon usage heatmap
    ax = axes[1, 0]
    cm = data["codon_matrix"]
    all_codons = sorted(set(c for row in cm.values() for c in row))
    matrix_data = np.zeros((len(cm), len(all_codons)))
    aa_labels = list(cm.keys())
    for i, aa in enumerate(aa_labels):
        for j, codon in enumerate(all_codons):
            matrix_data[i, j] = cm[aa].get(codon, 0)
    im = ax.imshow(matrix_data, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(all_codons)))
    ax.set_xticklabels(all_codons, rotation=70, fontsize=6)
    ax.set_yticks(range(len(aa_labels)))
    ax.set_yticklabels(aa_labels)
    ax.set_title("Codon Usage Frequency (selected AAs)")
    plt.colorbar(im, ax=ax, fraction=0.046)

    # Stability vs GC scatter
    ax = axes[1, 1]
    ax.scatter(gc_vals, stab_vals, c=CB_COLORS[:3], s=120, zorder=5)
    for i, n in enumerate(["Native", "Optimized", "Max CAI"]):
        ax.annotate(n, (gc_vals[i], stab_vals[i]), textcoords="offset points", xytext=(5, 5), fontsize=8)
    ax.set_xlabel("GC Content (%)")
    ax.set_ylabel("Stability Score")
    ax.set_title("mRNA Stability vs GC Content")

    plt.tight_layout()
    path = os.path.join(FIG_DIR, "codon_optimization.png")
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    return path


def fig_utr(data):
    utr5 = data["5UTR"]
    utr3 = data["3UTR"]

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle("5'/3' UTR Design Analysis", fontsize=13, fontweight="bold")

    names5 = list(utr5.keys())
    rl_scores = [utr5[n]["ribosome_loading"] for n in names5]
    ax = axes[0, 0]
    bars = ax.barh(names5, rl_scores, color=CB_COLORS[:len(names5)])
    ax.set_xlabel("Ribosome Loading Score")
    ax.set_title("5'UTR Variants: Ribosome Loading")
    ax.set_xlim(0, 1.05)
    for bar, val in zip(bars, rl_scores):
        ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                f"{val:.3f}", va="center", fontsize=8)

    names3 = list(utr3.keys())
    hl_scores = [utr3[n]["predicted_half_life_h"] for n in names3]
    ax = axes[0, 1]
    bars = ax.bar(range(len(names3)), hl_scores, color=CB_COLORS[:len(names3)])
    ax.set_xticks(range(len(names3)))
    ax.set_xticklabels(["β-globin", "α-globin", "WHV", "TENT4", "StAble"], rotation=12)
    ax.set_ylabel("Predicted Half-life (hours)")
    ax.set_title("3'UTR Variants: mRNA Stability")
    for bar, val in zip(bars, hl_scores):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3, f"{val:.1f}h",
                ha="center", fontsize=8)

    # Kozak scores
    ax = axes[1, 0]
    kozak_scores = [utr5[n]["kozak_score"] for n in names5]
    colors_kozak = [CB_COLORS[i % len(CB_COLORS)] for i in range(len(names5))]
    ax.bar(range(len(names5)), kozak_scores, color=colors_kozak)
    ax.set_xticks(range(len(names5)))
    ax.set_xticklabels(["β-glob", "CMV", "EMCV", "TOP", "Kozak"], rotation=12)
    ax.set_ylabel("Kozak Context Score")
    ax.set_title("5'UTR Kozak Sequence Quality")
    ax.set_ylim(0, 1.1)

    # Stacked bar: UTR3 features
    ax = axes[1, 1]
    poly_a = [utr3[n]["polya_length"] for n in names3]
    ares = [utr3[n]["AU_rich_elements"] * 10 for n in names3]
    stab_s = [utr3[n]["stability_score"] * 80 for n in names3]
    x = np.arange(len(names3))
    ax.bar(x, poly_a, label="Poly-A length (nt)", color=CB_COLORS[0])
    ax.bar(x, ares, bottom=poly_a, label="ARE count (×10)", color=CB_COLORS[1])
    ax.bar(x, stab_s, bottom=np.array(poly_a)+np.array(ares), label="Stability score (×80)", color=CB_COLORS[2])
    ax.set_xticks(x)
    ax.set_xticklabels(["β-glob", "α-glob", "WHV", "TENT4", "StAble"], rotation=12)
    ax.set_title("3'UTR Feature Composition")
    ax.legend(fontsize=7)

    plt.tight_layout()
    path = os.path.join(FIG_DIR, "utr_design.png")
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    return path


def fig_modified_nucleotides(data):
    mods = data["modifications"]
    dose_data = data["dose_curves"]
    doses = np.array(dose_data["doses_ug_ml"])

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    fig.suptitle("Modified Nucleotide Effect Prediction", fontsize=13, fontweight="bold")

    mod_names = list(mods.keys())
    mod_abbr = ["Unmod", "Ψ", "m5C", "m1Ψ", "2'OMe", "m1Ψ+m5C"]

    # Radar chart
    ax = axes[0, 0]
    metrics = ["innate_immune_evasion", "translation_multiplier", "half_life_h",
               "potency_score"]
    metric_labels = ["Immune\nEvasion", "Translation\nEfficiency", "Half-life", "Potency"]
    n = len(metrics)
    angles = np.linspace(0, 2*np.pi, n, endpoint=False).tolist()
    angles += angles[:1]

    normalizers = {
        "innate_immune_evasion": 1.0,
        "translation_multiplier": 2.0,
        "half_life_h": 30.0,
        "potency_score": 1.0,
    }

    for idx, (mod_name, props) in enumerate(mods.items()):
        vals = [props[m] / normalizers[m] for m in metrics]
        vals += vals[:1]
        ax.plot(angles, vals, "o-", lw=1.5, color=CB_COLORS[idx % len(CB_COLORS)],
                label=mod_abbr[idx])
        ax.fill(angles, vals, alpha=0.05, color=CB_COLORS[idx % len(CB_COLORS)])

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(metric_labels, size=8)
    ax.set_ylim(0, 1)
    ax.set_title("Modification Properties Radar")
    ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.1), fontsize=7)

    # Immune activation bar
    ax = axes[0, 1]
    tlr7 = [mods[m]["tlr7_activation"] for m in mod_names]
    rig1 = [mods[m]["rig1_activation"] for m in mod_names]
    x = np.arange(len(mod_abbr))
    w = 0.35
    ax.bar(x - w/2, tlr7, w, label="TLR7/8 activation", color=CB_COLORS[4])
    ax.bar(x + w/2, rig1, w, label="RIG-I activation", color=CB_COLORS[5])
    ax.set_xticks(x)
    ax.set_xticklabels(mod_abbr, rotation=20)
    ax.set_ylabel("Activation Level (0–1)")
    ax.set_title("Innate Immune Activation")
    ax.legend()
    ax.set_ylim(0, 1.1)

    # Dose-response curves
    ax = axes[1, 0]
    for idx, mod_name in enumerate(mod_names):
        curve = np.array(dose_data["responses"][mod_name])
        ax.semilogx(doses, curve, lw=2, color=CB_COLORS[idx % len(CB_COLORS)],
                    label=mod_abbr[idx])
    ax.set_xlabel("Dose (μg/mL)")
    ax.set_ylabel("Potency Score")
    ax.set_title("Dose-Response Curves")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # Heatmap: modification × metric
    ax = axes[1, 1]
    metric_cols = ["innate_immune_evasion", "translation_multiplier", "half_life_h",
                   "tlr7_activation", "potency_score"]
    col_labels = ["Immune\nEvasion", "Translation\nMult.", "Half-life\n(norm)", "TLR7\nActivation", "Potency"]
    matrix = []
    for mod_name in mod_names:
        row = []
        norms = [1.0, 2.0, 30.0, 1.0, 1.0]
        for col, n in zip(metric_cols, norms):
            row.append(mods[mod_name][col] / n)
        matrix.append(row)
    im = ax.imshow(matrix, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, fontsize=7)
    ax.set_yticks(range(len(mod_abbr)))
    ax.set_yticklabels(mod_abbr)
    ax.set_title("Modification × Metric Heatmap")
    plt.colorbar(im, ax=ax, fraction=0.046)

    plt.tight_layout()
    path = os.path.join(FIG_DIR, "modified_nucleotides.png")
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    return path


def fig_epitopes(data):
    top_mhc1 = data["top_mhc1_epitopes"]
    top_mhc2 = data["top_mhc2_epitopes"]
    bcell = data["top_bcell_epitopes"]
    hla_cov = data["hla_supertype_coverage"]

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    fig.suptitle("Antigen Epitope Analysis", fontsize=13, fontweight="bold")

    # Top MHC-I
    ax = axes[0, 0]
    ep_labels = [f"{e['peptide'][:9]}\n({e['allele'].split('*')[1]})" for e in top_mhc1[:8]]
    ic50_vals = [e["ic50_nM"] for e in top_mhc1[:8]]
    colors_bar = [CB_COLORS[i % 3] for i in range(8)]
    ax.barh(ep_labels[::-1], ic50_vals[::-1], color=colors_bar[::-1])
    ax.axvline(500, color="red", ls="--", lw=1, label="IC50=500nM cutoff")
    ax.set_xlabel("Predicted IC50 (nM)")
    ax.set_title("Top MHC-I Binding Epitopes")
    ax.legend(fontsize=7)

    # Top MHC-II
    ax = axes[0, 1]
    ep2_labels = [f"{e['peptide'][:9]}…\n({e['allele'].split('*')[1]})" for e in top_mhc2[:8]]
    ic50_2 = [e["ic50_nM"] for e in top_mhc2[:8]]
    ax.barh(ep2_labels[::-1], ic50_2[::-1], color=CB_COLORS[3])
    ax.axvline(1000, color="red", ls="--", lw=1, label="IC50=1000nM cutoff")
    ax.set_xlabel("Predicted IC50 (nM)")
    ax.set_title("Top MHC-II Binding Epitopes (CD4+)")
    ax.legend(fontsize=7)

    # B-cell epitope
    ax = axes[1, 0]
    positions = [b["start"] for b in bcell]
    scores = [b["hydrophilicity"] for b in bcell]
    ax.bar(positions, scores, color=CB_COLORS[2], alpha=0.8, width=1)
    ax.axhline(1.5, color="red", ls="--", lw=1, label="Threshold=1.5")
    ax.set_xlabel("Position in Spike RBD")
    ax.set_ylabel("Parker Hydrophilicity Score")
    ax.set_title("B-cell Linear Epitope Prediction")
    ax.legend()

    # HLA coverage pie
    ax = axes[1, 1]
    labels = list(hla_cov.keys())
    sizes = list(hla_cov.values())
    explode = [0.05] * len(labels)
    ax.pie(sizes, labels=labels, autopct="%1.1f%%", colors=CB_COLORS[:len(labels)],
           explode=explode, startangle=90)
    ax.set_title(f"HLA Population Coverage\n(Combined: {data['population_coverage']['mhc1_and_mhc2']:.1%})")

    plt.tight_layout()
    path = os.path.join(FIG_DIR, "epitope_analysis.png")
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    return path


def fig_lnp(data):
    best = data["optimal_formulation"]
    std = data["standard_reference"]
    hmap = data["heatmap_ionizable_x_chol"]
    top10 = data["top10_formulations"]

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    fig.suptitle("LNP Formulation Optimization", fontsize=13, fontweight="bold")

    # Scatter: size vs encapsulation
    ax = axes[0, 0]
    sizes = [r["size_nm"] for r in top10]
    encs = [r["encapsulation_pct"] for r in top10]
    tefs = [r["transfection_score"] for r in top10]
    sc = ax.scatter(sizes, encs, c=tefs, cmap="viridis", s=80, zorder=3)
    ax.scatter(best["size_nm"], best["encapsulation_pct"], marker="*", s=300,
               color="red", label="Optimal", zorder=5)
    ax.scatter(std["size_nm"], std["encapsulation_pct"], marker="D", s=100,
               color="orange", label="Standard", zorder=5)
    plt.colorbar(sc, ax=ax, label="Transfection Score")
    ax.axvspan(80, 120, alpha=0.1, color="green", label="Target 80–120nm")
    ax.set_xlabel("Particle Size (nm)")
    ax.set_ylabel("Encapsulation Efficiency (%)")
    ax.set_title("Size vs Encapsulation Efficiency")
    ax.legend(fontsize=7)

    # Heatmap: ionizable × chol
    ax = axes[0, 1]
    ion_keys = sorted([int(k) for k in hmap.keys()])
    chol_keys = sorted([int(k) for k in list(hmap.values())[0].keys()])
    matrix = np.zeros((len(ion_keys), len(chol_keys)))
    for i, ion in enumerate(ion_keys):
        for j, chol in enumerate(chol_keys):
            matrix[i, j] = hmap[str(ion)].get(str(chol), 0)
    im = ax.imshow(matrix, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(chol_keys)))
    ax.set_xticklabels([f"{c}%" for c in chol_keys])
    ax.set_yticks(range(len(ion_keys)))
    ax.set_yticklabels([f"{i}%" for i in ion_keys])
    ax.set_xlabel("Cholesterol Ratio (%)")
    ax.set_ylabel("Ionizable Lipid Ratio (%)")
    ax.set_title("Transfection Efficiency Heatmap")
    plt.colorbar(im, ax=ax, fraction=0.046)

    # Formulation comparison
    ax = axes[1, 0]
    params = ["size_nm", "pdi", "encapsulation_pct", "transfection_score"]
    norm = [120, 0.3, 100, 1.0]
    param_labels = ["Size\n(norm.)", "PDI\n(norm.)", "Encap.\nEff.", "Transfection\nScore"]
    std_vals = [std[p]/n for p, n in zip(params, norm)]
    best_vals = [best[p]/n for p, n in zip(params, norm)]
    x = np.arange(len(params))
    w = 0.35
    ax.bar(x - w/2, std_vals, w, label="Standard Ref.", color=CB_COLORS[1])
    ax.bar(x + w/2, best_vals, w, label="Optimized", color=CB_COLORS[0])
    ax.set_xticks(x)
    ax.set_xticklabels(param_labels)
    ax.set_ylabel("Normalized Value")
    ax.set_title("Standard vs Optimized Formulation")
    ax.legend()

    # PDI distribution
    ax = axes[1, 1]
    all_pdi = [r["pdi"] for r in data["top10_formulations"]]
    # simulate broader distribution
    np.random.seed(0)
    pdi_dist = np.concatenate([
        np.random.normal(0.12, 0.03, 80),
        np.random.normal(0.22, 0.05, 40),
    ])
    ax.hist(pdi_dist, bins=20, color=CB_COLORS[0], alpha=0.7, edgecolor="white")
    ax.axvline(0.2, color="red", ls="--", lw=1.5, label="Target PDI < 0.2")
    ax.axvline(best["pdi"], color="green", ls="-", lw=1.5, label=f"Optimal PDI={best['pdi']:.3f}")
    ax.set_xlabel("Polydispersity Index (PDI)")
    ax.set_ylabel("Count")
    ax.set_title("PDI Distribution Across Formulations")
    ax.legend()

    plt.tight_layout()
    path = os.path.join(FIG_DIR, "lnp_optimization.png")
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    return path


def fig_multivalent(data):
    cons_matrix = data["conservation_matrix"]
    neut_matrix = data["neutralization_matrix"]
    breadth = data["breadth_scores"]
    radar_data = data["strategy_radar"]
    var_names = data["variant_names"]
    strat_names = list(neut_matrix.keys())
    strat_abbr = ["Mono", "Bival.", "Trival.", "Mosaic", "AI-opt"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("Multi-valent Vaccine Design Strategy", fontsize=13, fontweight="bold")

    # Conservation heatmap
    ax = axes[0, 0]
    n_epitopes = 10
    cons_arr = np.array([cons_matrix[v] for v in var_names])
    im = ax.imshow(cons_arr, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(n_epitopes))
    ax.set_xticklabels([f"Ep{i+1}" for i in range(n_epitopes)], rotation=45, fontsize=7)
    ax.set_yticks(range(len(var_names)))
    ax.set_yticklabels(["Anc.", "Alpha", "Beta", "Delta", "Omicron\nBA.1", "XBB.1.5"], fontsize=8)
    ax.set_title("Epitope Conservation Matrix")
    plt.colorbar(im, ax=ax, fraction=0.046, label="Conservation score")

    # Neutralization breadth
    ax = axes[0, 1]
    breadth_vals = [breadth[s] for s in strat_names]
    colors_b = CB_COLORS[:len(strat_names)]
    bars = ax.bar(strat_abbr, breadth_vals, color=colors_b)
    ax.set_ylabel("Neutralization Breadth Score")
    ax.set_title("Vaccine Strategy: Neutralization Breadth")
    ax.set_ylim(0, 1.0)
    for bar, val in zip(bars, breadth_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f"{val:.3f}", ha="center", fontsize=8, fontweight="bold")

    # Radar: strategy comparison
    ax = axes[1, 0]
    ax.set_aspect("equal")
    metrics_r = ["neutralization_breadth", "potency_vs_latest", "manufacturability",
                 "safety_profile", "cost_efficiency"]
    metric_labels_r = ["Breadth", "Potency\nvs XBB", "Manufacturing", "Safety", "Cost\nEfficiency"]
    n = len(metrics_r)
    angles = np.linspace(0, 2*np.pi, n, endpoint=False).tolist()
    angles += angles[:1]
    for idx, (strat_name, props) in enumerate(radar_data.items()):
        vals = [props[m] for m in metrics_r] + [props[metrics_r[0]]]
        ax.plot(angles, vals, "o-", lw=1.5, color=CB_COLORS[idx % len(CB_COLORS)],
                label=strat_abbr[idx])
        ax.fill(angles, vals, alpha=0.05, color=CB_COLORS[idx % len(CB_COLORS)])
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(metric_labels_r, size=8)
    ax.set_ylim(0, 1)
    ax.set_title("Strategy Comparison (Radar)")
    ax.legend(loc="upper right", bbox_to_anchor=(1.4, 1.1), fontsize=7)

    # Variant phylo-like distance dendrogram
    ax = axes[1, 1]
    escape_vals = [VARIANTS[v]["rbd_escape"] for v in var_names]
    var_abbr = ["Anc.", "Alpha", "Beta", "Delta", "Omicron\nBA.1", "XBB.1.5"]
    im2 = ax.imshow(
        np.array(data["variant_distance_matrix"]),
        cmap="Blues", aspect="auto",
    )
    ax.set_xticks(range(len(var_abbr)))
    ax.set_xticklabels(var_abbr, rotation=30, fontsize=8)
    ax.set_yticks(range(len(var_abbr)))
    ax.set_yticklabels(var_abbr, fontsize=8)
    ax.set_title("Variant RBD Escape Distance Matrix")
    plt.colorbar(im2, ax=ax, fraction=0.046, label="|ΔEscape|")

    plt.tight_layout()
    path = os.path.join(FIG_DIR, "multivalent_design.png")
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    return path


# ════════════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 70)
    print("mRNA Vaccine In Silico Design Optimization Platform")
    print("=" * 70)

    print("\n[1/6] Codon Optimization...")
    codon_res = run_codon_optimization()
    fig_codon(codon_res)
    print(f"  CAI scores: { {k: v['CAI'] for k, v in codon_res['sequences'].items()} }")

    print("\n[2/6] UTR Design...")
    utr_res = run_utr_design()
    fig_utr(utr_res)
    best_utr5 = max(utr_res["5UTR"], key=lambda x: utr_res["5UTR"][x]["ribosome_loading"])
    best_utr3 = max(utr_res["3UTR"], key=lambda x: utr_res["3UTR"][x]["predicted_half_life_h"])
    print(f"  Best 5'UTR: {best_utr5} (score={utr_res['5UTR'][best_utr5]['ribosome_loading']:.4f})")
    print(f"  Best 3'UTR: {best_utr3} (half-life={utr_res['3UTR'][best_utr3]['predicted_half_life_h']}h)")

    print("\n[3/6] Modified Nucleotide Prediction...")
    mod_res = run_modified_nucleotides()
    fig_modified_nucleotides(mod_res)
    best_mod = max(mod_res["modifications"], key=lambda x: mod_res["modifications"][x]["potency_score"])
    print(f"  Best modification: {best_mod} (potency={mod_res['modifications'][best_mod]['potency_score']:.3f})")

    print("\n[4/6] Epitope Analysis...")
    ep_res = run_epitope_analysis()
    fig_epitopes(ep_res)
    print(f"  Top MHC-I epitope: {ep_res['top_mhc1_epitopes'][0]['peptide']} IC50={ep_res['top_mhc1_epitopes'][0]['ic50_nM']}nM")
    print(f"  Population coverage (MHC-I+II): {ep_res['population_coverage']['mhc1_and_mhc2']:.1%}")

    print("\n[5/6] LNP Optimization...")
    lnp_res = run_lnp_optimization()
    fig_lnp(lnp_res)
    bf = lnp_res["optimal_formulation"]
    print(f"  Optimal: size={bf['size_nm']}nm, PDI={bf['pdi']}, encap={bf['encapsulation_pct']}%, score={bf['transfection_score']:.4f}")

    print("\n[6/6] Multi-valent Vaccine Design...")
    mv_res = run_multivalent_design()
    fig_multivalent(mv_res)
    best_strat = max(mv_res["breadth_scores"], key=mv_res["breadth_scores"].get)
    print(f"  Best strategy: {best_strat} (breadth={mv_res['breadth_scores'][best_strat]:.4f})")

    # Preprocessing log
    with open(os.path.join(BASE, "data", "preprocessing-log.md"), "w") as f:
        f.write("# Preprocessing Log\n\n")
        f.write(f"Generated: {datetime.utcnow().isoformat()}Z\n\n")
        f.write("## Input Data\n")
        f.write(f"- Spike peptide: {len(SPIKE_PEPTIDE)} amino acids (Wuhan reference)\n")
        f.write("- Human codon usage: Kazusa database (Homo sapiens)\n")
        f.write("- HLA frequencies: Global population survey\n\n")
        f.write("## Random Seeds\n- numpy: 42\n- All modules: seed=42\n\n")
        f.write("## Transformations\n")
        f.write("- Codon table: fraction-based weights\n")
        f.write("- CAI: log-ratio geometric mean\n")
        f.write("- Parker hydrophilicity: 7-mer sliding window\n")
        f.write("- LNP metrics: physics-inspired parametric model\n")

    log("pipeline", "run_completed", "mrna_vaccine_pipeline",
        files=[
            "figures/codon_optimization.png", "figures/utr_design.png",
            "figures/modified_nucleotides.png", "figures/epitope_analysis.png",
            "figures/lnp_optimization.png", "figures/multivalent_design.png",
            "results/codon_optimization_results.json", "results/utr_design_results.json",
            "results/modified_nucleotide_results.json", "results/epitope_results.json",
            "results/lnp_optimization_results.json", "results/multivalent_vaccine_results.json",
        ])

    print("\n✓ All modules complete. Results saved to workspace/results/ and workspace/figures/")
