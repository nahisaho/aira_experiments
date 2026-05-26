"""
Biomedical Knowledge Graph Construction and Drug Repurposing via KGE Models.

This module builds a synthetic biomedical knowledge graph integrating
drug, gene, disease, pathway, and phenotype entities with relations
modeled after DrugBank, DisGeNET, STRING, and CTD data sources.
"""

import os
import json
import random
import numpy as np
import pandas as pd
import networkx as nx
from collections import defaultdict

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)))
DATA_DIR = os.path.join(OUTPUT_DIR, "data")
os.makedirs(DATA_DIR, exist_ok=True)


# ── Entity definitions ──────────────────────────────────────────────
DRUGS = [
    "Remdesivir", "Baricitinib", "Dexamethasone", "Hydroxychloroquine",
    "Ivermectin", "Tocilizumab", "Favipiravir", "Lopinavir",
    "Ritonavir", "Molnupiravir", "Nirmatrelvir", "Ruxolitinib",
    "Colchicine", "Fluvoxamine", "Metformin", "Aspirin",
    "Ibuprofen", "Azithromycin", "Chloroquine", "Sofosbuvir",
    "Ribavirin", "Interferon-beta", "Camostat", "Nafamostat",
    "Umifenovir", "Nitazoxanide", "Famotidine", "Losartan",
    "Atorvastatin", "Sirolimus"
]

GENES = [
    "ACE2", "TMPRSS2", "IL6", "TNF", "IFNG", "JAK1", "JAK2",
    "STAT3", "NFKB1", "TP53", "EGFR", "VEGFA", "MAPK1", "AKT1",
    "MTOR", "CASP3", "BCL2", "CXCL8", "CCL2", "IL1B",
    "IL10", "TLR4", "MYD88", "TGFB1", "CDK2", "PARP1",
    "HDAC1", "ESR1", "AR", "PPARG", "RdRp", "3CLpro", "PLpro",
    "Spike", "Nucleocapsid", "ORF3a", "ORF7a", "NSP1", "NSP3", "NSP12"
]

DISEASES = [
    "COVID-19", "SARS", "MERS", "Influenza", "Rheumatoid Arthritis",
    "Systemic Lupus Erythematosus", "Type 2 Diabetes", "Hypertension",
    "Alzheimer Disease", "Parkinson Disease", "Breast Cancer",
    "Lung Cancer", "Colorectal Cancer", "Asthma", "COPD",
    "Crohn Disease", "Ulcerative Colitis", "Multiple Sclerosis",
    "HIV/AIDS", "Hepatitis C"
]

PATHWAYS = [
    "JAK-STAT signaling", "NF-kB signaling", "PI3K-AKT signaling",
    "MAPK signaling", "TNF signaling", "Toll-like receptor signaling",
    "Cytokine-cytokine receptor interaction", "Apoptosis",
    "mTOR signaling", "HIF-1 signaling", "IL-17 signaling",
    "Chemokine signaling", "RIG-I-like receptor signaling",
    "NOD-like receptor signaling", "Complement and coagulation cascades"
]

PHENOTYPES = [
    "Cytokine storm", "Acute respiratory distress", "Pulmonary fibrosis",
    "Thrombocytopenia", "Lymphopenia", "Fever", "Hypoxia",
    "Endothelial dysfunction", "Cardiac injury", "Renal failure",
    "Coagulopathy", "Neuroinflammation", "Hepatotoxicity",
    "Immune evasion", "Viral replication"
]

RELATION_TYPES = {
    "drug_targets_gene": ("Drug", "Gene"),
    "drug_treats_disease": ("Drug", "Disease"),
    "gene_associated_disease": ("Gene", "Disease"),
    "gene_participates_pathway": ("Gene", "Pathway"),
    "disease_has_phenotype": ("Disease", "Phenotype"),
    "drug_interacts_drug": ("Drug", "Drug"),
    "gene_interacts_gene": ("Gene", "Gene"),
    "pathway_involves_phenotype": ("Pathway", "Phenotype"),
    "drug_inhibits_gene": ("Drug", "Gene"),
    "drug_upregulates_gene": ("Drug", "Gene"),
}


def build_triples():
    """Generate biologically-motivated triples for the KG."""
    triples = []

    # Drug-Gene targeting (DrugBank-inspired)
    covid_drug_gene = {
        "Remdesivir": ["RdRp", "NSP12"],
        "Baricitinib": ["JAK1", "JAK2"],
        "Dexamethasone": ["NFKB1", "IL6", "TNF"],
        "Hydroxychloroquine": ["TLR4", "ACE2"],
        "Tocilizumab": ["IL6"],
        "Favipiravir": ["RdRp"],
        "Lopinavir": ["3CLpro"],
        "Ritonavir": ["3CLpro", "PLpro"],
        "Molnupiravir": ["RdRp", "NSP12"],
        "Nirmatrelvir": ["3CLpro"],
        "Ruxolitinib": ["JAK1", "JAK2"],
        "Camostat": ["TMPRSS2"],
        "Nafamostat": ["TMPRSS2"],
        "Sofosbuvir": ["RdRp"],
        "Interferon-beta": ["IFNG", "STAT3"],
        "Colchicine": ["NFKB1", "CASP3"],
        "Fluvoxamine": ["IL6", "IL1B"],
        "Sirolimus": ["MTOR"],
        "Losartan": ["ACE2"],
        "Atorvastatin": ["AKT1", "MAPK1"],
        "Metformin": ["MTOR", "AKT1"],
        "Famotidine": ["3CLpro", "PLpro"],
    }
    for drug, genes in covid_drug_gene.items():
        for gene in genes:
            triples.append((drug, "drug_targets_gene", gene))

    # Additional random drug-gene targets
    for drug in DRUGS:
        n = random.randint(1, 3)
        for gene in random.sample(GENES, n):
            triples.append((drug, "drug_targets_gene", gene))

    # Drug-Disease treatments (known + speculative)
    known_treatments = {
        "Remdesivir": ["COVID-19"],
        "Baricitinib": ["COVID-19", "Rheumatoid Arthritis"],
        "Dexamethasone": ["COVID-19", "Asthma", "COPD"],
        "Tocilizumab": ["COVID-19", "Rheumatoid Arthritis"],
        "Hydroxychloroquine": ["Rheumatoid Arthritis", "Systemic Lupus Erythematosus"],
        "Ruxolitinib": ["Rheumatoid Arthritis"],
        "Metformin": ["Type 2 Diabetes"],
        "Losartan": ["Hypertension"],
        "Atorvastatin": ["Hypertension"],
        "Interferon-beta": ["Multiple Sclerosis", "Hepatitis C"],
        "Sofosbuvir": ["Hepatitis C"],
        "Ribavirin": ["Hepatitis C"],
    }
    for drug, diseases in known_treatments.items():
        for disease in diseases:
            triples.append((drug, "drug_treats_disease", disease))

    # Gene-Disease associations (DisGeNET-inspired)
    gene_disease = {
        "ACE2": ["COVID-19", "Hypertension"],
        "TMPRSS2": ["COVID-19"],
        "IL6": ["COVID-19", "Rheumatoid Arthritis", "Crohn Disease"],
        "TNF": ["Rheumatoid Arthritis", "Crohn Disease", "Ulcerative Colitis"],
        "JAK1": ["Rheumatoid Arthritis"],
        "JAK2": ["Rheumatoid Arthritis"],
        "TP53": ["Breast Cancer", "Lung Cancer", "Colorectal Cancer"],
        "EGFR": ["Lung Cancer", "Breast Cancer"],
        "VEGFA": ["Lung Cancer", "Colorectal Cancer"],
        "MTOR": ["Breast Cancer", "Type 2 Diabetes"],
        "NFKB1": ["COVID-19", "Rheumatoid Arthritis", "Asthma"],
        "STAT3": ["Lung Cancer", "Breast Cancer"],
        "MAPK1": ["Lung Cancer", "Colorectal Cancer"],
        "AKT1": ["Breast Cancer", "Type 2 Diabetes"],
        "BCL2": ["Breast Cancer", "Lung Cancer"],
        "IL1B": ["COVID-19", "Rheumatoid Arthritis"],
        "CXCL8": ["COVID-19", "COPD"],
        "CCL2": ["COVID-19", "Alzheimer Disease"],
        "TLR4": ["COVID-19", "Asthma"],
        "3CLpro": ["COVID-19"],
        "PLpro": ["COVID-19"],
        "RdRp": ["COVID-19"],
        "Spike": ["COVID-19", "SARS"],
        "Nucleocapsid": ["COVID-19", "SARS"],
    }
    for gene, diseases in gene_disease.items():
        for disease in diseases:
            triples.append((gene, "gene_associated_disease", disease))

    # Gene-Pathway (STRING/KEGG-inspired)
    gene_pathway = {
        "JAK1": ["JAK-STAT signaling"],
        "JAK2": ["JAK-STAT signaling"],
        "STAT3": ["JAK-STAT signaling", "PI3K-AKT signaling"],
        "NFKB1": ["NF-kB signaling", "TNF signaling"],
        "AKT1": ["PI3K-AKT signaling", "mTOR signaling"],
        "MAPK1": ["MAPK signaling"],
        "MTOR": ["mTOR signaling", "PI3K-AKT signaling"],
        "TNF": ["TNF signaling", "NF-kB signaling"],
        "TLR4": ["Toll-like receptor signaling", "NF-kB signaling"],
        "CASP3": ["Apoptosis"],
        "BCL2": ["Apoptosis"],
        "IL6": ["JAK-STAT signaling", "Cytokine-cytokine receptor interaction"],
        "IL1B": ["NF-kB signaling", "IL-17 signaling"],
        "CXCL8": ["Chemokine signaling", "IL-17 signaling"],
        "CCL2": ["Chemokine signaling"],
        "IFNG": ["JAK-STAT signaling", "Cytokine-cytokine receptor interaction"],
        "MYD88": ["Toll-like receptor signaling", "NF-kB signaling"],
    }
    for gene, pathways in gene_pathway.items():
        for pathway in pathways:
            triples.append((gene, "gene_participates_pathway", pathway))

    # Disease-Phenotype
    disease_phenotype = {
        "COVID-19": ["Cytokine storm", "Acute respiratory distress", "Lymphopenia",
                      "Hypoxia", "Coagulopathy", "Fever", "Viral replication"],
        "SARS": ["Acute respiratory distress", "Fever", "Lymphopenia"],
        "MERS": ["Acute respiratory distress", "Renal failure"],
        "Rheumatoid Arthritis": ["Cytokine storm", "Neuroinflammation"],
        "Asthma": ["Acute respiratory distress", "Hypoxia"],
        "COPD": ["Pulmonary fibrosis", "Hypoxia"],
        "Breast Cancer": ["Immune evasion"],
        "Lung Cancer": ["Immune evasion", "Pulmonary fibrosis"],
        "Type 2 Diabetes": ["Endothelial dysfunction"],
        "Hypertension": ["Endothelial dysfunction", "Cardiac injury"],
        "Alzheimer Disease": ["Neuroinflammation"],
    }
    for disease, phenotypes in disease_phenotype.items():
        for pheno in phenotypes:
            triples.append((disease, "disease_has_phenotype", pheno))

    # Gene-Gene interactions (STRING-inspired)
    ppi_pairs = [
        ("ACE2", "TMPRSS2"), ("JAK1", "JAK2"), ("JAK1", "STAT3"),
        ("JAK2", "STAT3"), ("NFKB1", "TNF"), ("NFKB1", "IL1B"),
        ("AKT1", "MTOR"), ("MAPK1", "AKT1"), ("TP53", "BCL2"),
        ("TP53", "CASP3"), ("IL6", "STAT3"), ("TLR4", "MYD88"),
        ("MYD88", "NFKB1"), ("EGFR", "MAPK1"), ("EGFR", "AKT1"),
        ("3CLpro", "NSP3"), ("RdRp", "NSP12"), ("Spike", "ACE2"),
    ]
    for g1, g2 in ppi_pairs:
        triples.append((g1, "gene_interacts_gene", g2))

    # Pathway-Phenotype
    pathway_phenotype = {
        "JAK-STAT signaling": ["Cytokine storm", "Immune evasion"],
        "NF-kB signaling": ["Cytokine storm", "Neuroinflammation"],
        "TNF signaling": ["Cytokine storm", "Fever"],
        "Apoptosis": ["Lymphopenia"],
        "mTOR signaling": ["Viral replication", "Immune evasion"],
        "Complement and coagulation cascades": ["Coagulopathy", "Thrombocytopenia"],
        "HIF-1 signaling": ["Hypoxia"],
        "Chemokine signaling": ["Neuroinflammation"],
    }
    for pathway, phenos in pathway_phenotype.items():
        for pheno in phenos:
            triples.append((pathway, "pathway_involves_phenotype", pheno))

    # Drug-Drug interactions
    ddi_pairs = [
        ("Lopinavir", "Ritonavir"), ("Remdesivir", "Baricitinib"),
        ("Hydroxychloroquine", "Azithromycin"), ("Nirmatrelvir", "Ritonavir"),
    ]
    for d1, d2 in ddi_pairs:
        triples.append((d1, "drug_interacts_drug", d2))

    # Drug inhibits/upregulates
    for drug in random.sample(DRUGS, 10):
        gene = random.choice(GENES[:20])
        triples.append((drug, "drug_inhibits_gene", gene))
    for drug in random.sample(DRUGS, 8):
        gene = random.choice(GENES[:20])
        triples.append((drug, "drug_upregulates_gene", gene))

    # Deduplicate
    triples = list(set(triples))
    return triples


def create_entity_type_map():
    """Create mapping of entity names to types."""
    etype = {}
    for d in DRUGS:
        etype[d] = "Drug"
    for g in GENES:
        etype[g] = "Gene"
    for dis in DISEASES:
        etype[dis] = "Disease"
    for p in PATHWAYS:
        etype[p] = "Pathway"
    for ph in PHENOTYPES:
        etype[ph] = "Phenotype"
    return etype


def save_triples(triples, filepath):
    """Save triples to TSV."""
    df = pd.DataFrame(triples, columns=["head", "relation", "tail"])
    df.to_csv(filepath, sep="\t", index=False)
    return df


def build_networkx_graph(triples):
    """Build NetworkX MultiDiGraph from triples."""
    G = nx.MultiDiGraph()
    etype = create_entity_type_map()
    for h, r, t in triples:
        G.add_node(h, entity_type=etype.get(h, "Unknown"))
        G.add_node(t, entity_type=etype.get(t, "Unknown"))
        G.add_edge(h, t, relation=r)
    return G


def compute_graph_statistics(G, triples):
    """Compute graph statistics."""
    stats = {
        "num_entities": G.number_of_nodes(),
        "num_triples": len(triples),
        "num_relations": len(set(r for _, r, _ in triples)),
        "num_drugs": sum(1 for _, d in G.nodes(data="entity_type") if d == "Drug"),
        "num_genes": sum(1 for _, d in G.nodes(data="entity_type") if d == "Gene"),
        "num_diseases": sum(1 for _, d in G.nodes(data="entity_type") if d == "Disease"),
        "num_pathways": sum(1 for _, d in G.nodes(data="entity_type") if d == "Pathway"),
        "num_phenotypes": sum(1 for _, d in G.nodes(data="entity_type") if d == "Phenotype"),
        "avg_degree": np.mean([d for _, d in G.degree()]),
        "density": nx.density(G),
    }

    # Data source contribution (simulated)
    stats["drugbank_triples"] = sum(1 for _, r, _ in triples if r in
                                     ["drug_targets_gene", "drug_treats_disease", "drug_interacts_drug"])
    stats["disgenet_triples"] = sum(1 for _, r, _ in triples if r == "gene_associated_disease")
    stats["string_triples"] = sum(1 for _, r, _ in triples if r in
                                   ["gene_interacts_gene", "gene_participates_pathway"])
    stats["ctd_triples"] = sum(1 for _, r, _ in triples if r in
                                ["drug_inhibits_gene", "drug_upregulates_gene",
                                 "disease_has_phenotype", "pathway_involves_phenotype"])

    return stats


if __name__ == "__main__":
    print("Building biomedical knowledge graph...")
    triples = build_triples()
    df = save_triples(triples, os.path.join(DATA_DIR, "kg_triples.tsv"))
    G = build_networkx_graph(triples)
    stats = compute_graph_statistics(G, triples)

    print(f"\n=== Knowledge Graph Statistics ===")
    for k, v in stats.items():
        print(f"  {k}: {v}")

    with open(os.path.join(DATA_DIR, "kg_stats.json"), "w") as f:
        json.dump(stats, f, indent=2, default=str)

    # Save entity type mapping
    etype = create_entity_type_map()
    with open(os.path.join(DATA_DIR, "entity_types.json"), "w") as f:
        json.dump(etype, f, indent=2)

    print(f"\nTriples saved to {DATA_DIR}/kg_triples.tsv")
    print(f"Statistics saved to {DATA_DIR}/kg_stats.json")
