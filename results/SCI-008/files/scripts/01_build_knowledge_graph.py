#!/usr/bin/env python3
"""
Step 1: Build a Biomedical Knowledge Graph for Drug Repurposing
Integrates data from DrugBank, DisGeNET, STRING, and CTD (simulated).
"""

import json
import os
import random
import csv
from datetime import datetime

import numpy as np
import pandas as pd
import networkx as nx

random.seed(42)
np.random.seed(42)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
LOG_DIR = os.path.join(BASE_DIR, "logs")
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

# ─── Entity definitions ───

DRUGS = [
    ("DB00945", "Aspirin"), ("DB00563", "Methotrexate"), ("DB01050", "Ibuprofen"),
    ("DB00503", "Ritonavir"), ("DB00608", "Chloroquine"), ("DB01611", "Hydroxychloroquine"),
    ("DB00641", "Simvastatin"), ("DB01076", "Atorvastatin"), ("DB00959", "Methylprednisolone"),
    ("DB00635", "Prednisone"), ("DB09065", "Baricitinib"), ("DB11817", "Remdesivir"),
    ("DB14761", "Molnupiravir"), ("DB16691", "Nirmatrelvir"), ("DB00631", "Clopidogrel"),
    ("DB01136", "Carvedilol"), ("DB00177", "Valsartan"), ("DB00966", "Telmisartan"),
    ("DB01197", "Captopril"), ("DB00691", "Finasteride"), ("DB00619", "Imatinib"),
    ("DB00530", "Erlotinib"), ("DB08901", "Ponatinib"), ("DB06290", "Simeprevir"),
    ("DB00696", "Ergotamine"), ("DB01048", "Abacavir"), ("DB00224", "Indinavir"),
    ("DB01601", "Lopinavir"), ("DB00932", "Tipranavir"), ("DB01264", "Darunavir"),
    ("DB00198", "Oseltamivir"), ("DB00558", "Zanamivir"), ("DB06817", "Ruxolitinib"),
    ("DB11979", "Eculizumab"), ("DB06273", "Tocilizumab"), ("DB00051", "Adalimumab"),
    ("DB01041", "Thalidomide"), ("DB00674", "Galantamine"), ("DB01065", "Melatonin"),
    ("DB00252", "Phenytoin"),
]

GENES = [
    ("ACE2", "Angiotensin Converting Enzyme 2"), ("TMPRSS2", "Transmembrane Serine Protease 2"),
    ("IL6", "Interleukin 6"), ("TNF", "Tumor Necrosis Factor"), ("IFNG", "Interferon Gamma"),
    ("JAK1", "Janus Kinase 1"), ("JAK2", "Janus Kinase 2"), ("STAT3", "Signal Transducer 3"),
    ("NFKB1", "Nuclear Factor Kappa B Subunit 1"), ("TP53", "Tumor Protein P53"),
    ("EGFR", "Epidermal Growth Factor Receptor"), ("VEGFA", "Vascular Endothelial GF A"),
    ("MAPK1", "Mitogen-Activated Protein Kinase 1"), ("AKT1", "AKT Serine/Threonine Kinase 1"),
    ("MTOR", "Mechanistic Target Of Rapamycin"), ("CASP3", "Caspase 3"),
    ("BCL2", "B-Cell Lymphoma 2"), ("IL1B", "Interleukin 1 Beta"),
    ("CXCL8", "C-X-C Motif Chemokine Ligand 8"), ("CCL2", "C-C Motif Chemokine Ligand 2"),
    ("TLR4", "Toll Like Receptor 4"), ("NLRP3", "NLR Family Pyrin Domain 3"),
    ("HMGCR", "HMG-CoA Reductase"), ("PTGS2", "Prostaglandin-Endoperoxide Synthase 2"),
    ("ADRB1", "Adrenoceptor Beta 1"), ("AGTR1", "Angiotensin II Receptor Type 1"),
    ("REN", "Renin"), ("AGT", "Angiotensinogen"), ("F2", "Coagulation Factor II"),
    ("SERPINE1", "Serpin Family E Member 1"), ("CTSL", "Cathepsin L"),
    ("BSG", "Basigin (CD147)"), ("DPP4", "Dipeptidyl Peptidase 4"),
    ("FURIN", "Furin"), ("CD4", "CD4 Molecule"),
]

DISEASES = [
    ("DOID:0080600", "COVID-19"), ("DOID:2841", "Asthma"), ("DOID:3393", "Coronary Artery Disease"),
    ("DOID:9352", "Type 2 Diabetes"), ("DOID:10763", "Hypertension"),
    ("DOID:2377", "Multiple Sclerosis"), ("DOID:10283", "Prostate Cancer"),
    ("DOID:3571", "Liver Cancer"), ("DOID:1612", "Breast Cancer"),
    ("DOID:219", "Colon Cancer"), ("DOID:332", "Parkinson Disease"),
    ("DOID:14330", "Alzheimer Disease"), ("DOID:9074", "Systemic Lupus Erythematosus"),
    ("DOID:7148", "Rheumatoid Arthritis"), ("DOID:0050117", "Acute Respiratory Distress"),
    ("DOID:0050741", "Pulmonary Fibrosis"), ("DOID:4", "Cardiovascular Disease"),
    ("DOID:0060903", "Thrombotic Disorder"), ("DOID:0080599", "Cytokine Storm"),
    ("DOID:552", "Pneumonia"),
]

PATHWAYS = [
    ("R-HSA-168256", "Immune System"), ("R-HSA-449147", "Signaling by Interleukins"),
    ("R-HSA-168164", "Toll-Like Receptor Cascades"), ("R-HSA-1280215", "Cytokine Signaling"),
    ("R-HSA-109582", "Hemostasis"), ("R-HSA-2262752", "Cellular Responses to Stress"),
    ("R-HSA-5653656", "Vesicle-mediated Transport"), ("R-HSA-1643685", "Disease"),
    ("R-HSA-162582", "Signal Transduction"), ("R-HSA-556833", "Metabolism of Lipids"),
    ("R-HSA-1474244", "Extracellular Matrix Organization"), ("R-HSA-74160", "Gene Expression"),
    ("R-HSA-1280218", "Adaptive Immune System"), ("R-HSA-168249", "Innate Immune System"),
    ("R-HSA-9006925", "Intracellular Signaling by Second Messengers"),
    ("hsa04668", "TNF Signaling Pathway"), ("hsa04630", "JAK-STAT Signaling Pathway"),
    ("hsa04151", "PI3K-Akt Signaling Pathway"), ("hsa04064", "NF-kappa B Signaling Pathway"),
    ("hsa04620", "Toll-like Receptor Signaling Pathway"),
]

PHENOTYPES = [
    ("HP:0001945", "Fever"), ("HP:0012735", "Cough"), ("HP:0002014", "Diarrhea"),
    ("HP:0002094", "Dyspnea"), ("HP:0001907", "Thromboembolism"),
    ("HP:0001919", "Acute Kidney Injury"), ("HP:0002321", "Vertigo"),
    ("HP:0100806", "Sepsis"), ("HP:0003326", "Myalgia"), ("HP:0002315", "Headache"),
    ("HP:0012115", "Hepatitis"), ("HP:0001250", "Seizures"),
    ("HP:0001649", "Tachycardia"), ("HP:0002090", "Pneumonia"),
    ("HP:0003493", "Antinuclear Antibody Positive"),
]

# ─── Build edges (relations) ───

def build_kg():
    triples = []

    # Drug-Gene (target) relations — DrugBank-like
    drug_gene_targets = {
        "DB00945": ["PTGS2", "NFKB1"],
        "DB00563": ["NFKB1", "TP53"],
        "DB01050": ["PTGS2"],
        "DB00503": ["CTSL", "FURIN"],
        "DB00608": ["ACE2", "TMPRSS2", "TLR4"],
        "DB01611": ["ACE2", "TMPRSS2", "TLR4", "NLRP3"],
        "DB00641": ["HMGCR", "NFKB1", "AKT1"],
        "DB01076": ["HMGCR", "NFKB1"],
        "DB00959": ["NFKB1", "IL6", "TNF"],
        "DB00635": ["NFKB1", "IL6", "TNF"],
        "DB09065": ["JAK1", "JAK2"],
        "DB11817": ["FURIN"],
        "DB14761": ["FURIN"],
        "DB16691": ["CTSL", "FURIN"],
        "DB00631": ["F2"],
        "DB01136": ["ADRB1"],
        "DB00177": ["AGTR1"],
        "DB00966": ["AGTR1", "NFKB1"],
        "DB01197": ["ACE2", "REN"],
        "DB00619": ["EGFR", "MAPK1", "AKT1"],
        "DB00530": ["EGFR"],
        "DB08901": ["EGFR", "VEGFA"],
        "DB06290": ["FURIN"],
        "DB01048": ["FURIN"],
        "DB00224": ["CTSL"],
        "DB01601": ["CTSL", "FURIN"],
        "DB00932": ["CTSL"],
        "DB01264": ["CTSL", "FURIN"],
        "DB00198": ["FURIN"],
        "DB00558": ["FURIN"],
        "DB06817": ["JAK1", "JAK2", "STAT3"],
        "DB11979": ["IL6"],
        "DB06273": ["IL6"],
        "DB00051": ["TNF"],
        "DB01041": ["TNF", "NFKB1", "VEGFA"],
        "DB00674": ["CASP3"],
        "DB01065": ["NLRP3", "NFKB1"],
        "DB00252": ["CASP3"],
        "DB00691": ["AGT"],
        "DB00696": ["EGFR"],
    }

    for drug_id, genes in drug_gene_targets.items():
        for gene in genes:
            triples.append((drug_id, "targets", gene))

    # Gene-Disease associations — DisGeNET-like
    gene_disease = {
        "ACE2": ["DOID:0080600", "DOID:10763", "DOID:4"],
        "TMPRSS2": ["DOID:0080600", "DOID:10283"],
        "IL6": ["DOID:0080600", "DOID:7148", "DOID:0080599", "DOID:552"],
        "TNF": ["DOID:7148", "DOID:9074", "DOID:0080599"],
        "IFNG": ["DOID:2377", "DOID:0080600"],
        "JAK1": ["DOID:7148", "DOID:0080600"],
        "JAK2": ["DOID:0080600", "DOID:1612"],
        "STAT3": ["DOID:3571", "DOID:1612", "DOID:0080600"],
        "NFKB1": ["DOID:0080600", "DOID:7148", "DOID:9074", "DOID:0080599"],
        "TP53": ["DOID:1612", "DOID:3571", "DOID:219", "DOID:10283"],
        "EGFR": ["DOID:1612", "DOID:3571"],
        "VEGFA": ["DOID:1612", "DOID:3393"],
        "MAPK1": ["DOID:1612", "DOID:0080600"],
        "AKT1": ["DOID:9352", "DOID:1612"],
        "MTOR": ["DOID:9352", "DOID:1612"],
        "CASP3": ["DOID:14330", "DOID:332"],
        "BCL2": ["DOID:1612", "DOID:219"],
        "IL1B": ["DOID:0080600", "DOID:7148", "DOID:0080599"],
        "CXCL8": ["DOID:0080600", "DOID:0050117"],
        "CCL2": ["DOID:0080600", "DOID:0050117"],
        "TLR4": ["DOID:0080600", "DOID:552"],
        "NLRP3": ["DOID:0080600", "DOID:0080599", "DOID:14330"],
        "HMGCR": ["DOID:3393", "DOID:4"],
        "PTGS2": ["DOID:219", "DOID:7148"],
        "ADRB1": ["DOID:4", "DOID:10763"],
        "AGTR1": ["DOID:10763", "DOID:4"],
        "REN": ["DOID:10763"],
        "AGT": ["DOID:10763", "DOID:0080600"],
        "F2": ["DOID:0060903", "DOID:0080600"],
        "SERPINE1": ["DOID:0060903", "DOID:0080600"],
        "CTSL": ["DOID:0080600"],
        "BSG": ["DOID:0080600"],
        "DPP4": ["DOID:9352", "DOID:0080600"],
        "FURIN": ["DOID:0080600"],
        "CD4": ["DOID:0080600"],
    }

    for gene, diseases in gene_disease.items():
        for disease in diseases:
            triples.append((gene, "associated_with", disease))

    # Gene-Gene interactions — STRING-like
    ppi_edges = [
        ("ACE2", "TMPRSS2"), ("ACE2", "AGT"), ("ACE2", "REN"), ("ACE2", "AGTR1"),
        ("IL6", "STAT3"), ("IL6", "JAK1"), ("IL6", "JAK2"), ("IL6", "NFKB1"),
        ("TNF", "NFKB1"), ("TNF", "CASP3"), ("TNF", "IL6"), ("TNF", "IL1B"),
        ("JAK1", "STAT3"), ("JAK2", "STAT3"), ("JAK1", "JAK2"),
        ("NFKB1", "IL1B"), ("NFKB1", "CXCL8"), ("NFKB1", "CCL2"),
        ("EGFR", "MAPK1"), ("EGFR", "AKT1"), ("MAPK1", "AKT1"),
        ("AKT1", "MTOR"), ("TP53", "BCL2"), ("TP53", "CASP3"),
        ("NLRP3", "IL1B"), ("NLRP3", "CASP3"), ("TLR4", "NFKB1"),
        ("HMGCR", "AKT1"), ("PTGS2", "NFKB1"), ("F2", "SERPINE1"),
        ("CTSL", "FURIN"), ("CTSL", "BSG"), ("FURIN", "ACE2"),
        ("IFNG", "JAK1"), ("IFNG", "STAT3"), ("DPP4", "AKT1"),
        ("ADRB1", "AGTR1"), ("AGT", "REN"), ("AGT", "AGTR1"),
        ("VEGFA", "MAPK1"), ("BCL2", "CASP3"), ("CCL2", "TLR4"),
    ]
    for g1, g2 in ppi_edges:
        triples.append((g1, "interacts_with", g2))
        triples.append((g2, "interacts_with", g1))

    # Gene-Pathway associations
    gene_pathway = {
        "IL6": ["R-HSA-449147", "R-HSA-1280215", "hsa04630"],
        "TNF": ["hsa04668", "R-HSA-1280215", "hsa04064"],
        "JAK1": ["hsa04630", "R-HSA-449147"],
        "JAK2": ["hsa04630"],
        "STAT3": ["hsa04630", "R-HSA-162582"],
        "NFKB1": ["hsa04064", "R-HSA-449147", "hsa04620"],
        "EGFR": ["hsa04151", "R-HSA-162582"],
        "MAPK1": ["hsa04151", "R-HSA-162582", "R-HSA-9006925"],
        "AKT1": ["hsa04151", "R-HSA-162582"],
        "MTOR": ["hsa04151"],
        "TLR4": ["hsa04620", "R-HSA-168164", "R-HSA-168249"],
        "NLRP3": ["R-HSA-168249", "R-HSA-1280215"],
        "CASP3": ["R-HSA-2262752"],
        "IL1B": ["R-HSA-449147", "R-HSA-1280215"],
        "ACE2": ["R-HSA-1643685"],
        "TMPRSS2": ["R-HSA-1643685"],
        "F2": ["R-HSA-109582"],
        "SERPINE1": ["R-HSA-109582"],
        "HMGCR": ["R-HSA-556833"],
        "PTGS2": ["R-HSA-556833", "hsa04668"],
        "FURIN": ["R-HSA-5653656"],
        "CTSL": ["R-HSA-5653656"],
        "BCL2": ["R-HSA-2262752"],
        "IFNG": ["R-HSA-1280215", "hsa04630"],
        "VEGFA": ["hsa04151"],
    }
    for gene, pathways in gene_pathway.items():
        for pw in pathways:
            triples.append((gene, "participates_in", pw))

    # Disease-Phenotype associations
    disease_phenotype = {
        "DOID:0080600": ["HP:0001945", "HP:0012735", "HP:0002094", "HP:0001907", "HP:0003326", "HP:0002090"],
        "DOID:552": ["HP:0001945", "HP:0012735", "HP:0002094", "HP:0002090"],
        "DOID:0080599": ["HP:0001945", "HP:0001649", "HP:0002094"],
        "DOID:0050117": ["HP:0002094", "HP:0001945"],
        "DOID:7148": ["HP:0003326", "HP:0003493"],
        "DOID:9074": ["HP:0003493", "HP:0001945"],
        "DOID:10763": ["HP:0002315", "HP:0001649"],
        "DOID:3393": ["HP:0002315", "HP:0001649"],
        "DOID:14330": ["HP:0002315", "HP:0001250"],
        "DOID:332": ["HP:0002321"],
        "DOID:9352": ["HP:0002094"],
        "DOID:1612": ["HP:0003326"],
        "DOID:0060903": ["HP:0001907"],
    }
    for disease, phenotypes in disease_phenotype.items():
        for ph in phenotypes:
            triples.append((disease, "has_phenotype", ph))

    # Known drug-disease (indication) — CTD-like
    drug_disease_known = {
        "DB00608": ["DOID:0080600"],  # Chloroquine (initial trials)
        "DB01611": ["DOID:0080600", "DOID:9074"],  # HCQ
        "DB11817": ["DOID:0080600"],  # Remdesivir
        "DB09065": ["DOID:7148", "DOID:0080600"],  # Baricitinib
        "DB06273": ["DOID:7148", "DOID:0080600"],  # Tocilizumab
        "DB00959": ["DOID:0050117", "DOID:7148"],  # Methylprednisolone
        "DB00635": ["DOID:7148", "DOID:9074"],  # Prednisone
        "DB00945": ["DOID:3393", "DOID:0060903"],  # Aspirin
        "DB00641": ["DOID:3393", "DOID:4"],  # Simvastatin
        "DB01076": ["DOID:3393"],  # Atorvastatin
        "DB00177": ["DOID:10763"],  # Valsartan
        "DB00966": ["DOID:10763"],  # Telmisartan
        "DB01197": ["DOID:10763"],  # Captopril
        "DB00619": ["DOID:1612"],  # Imatinib
        "DB00530": ["DOID:3571"],  # Erlotinib
        "DB01601": ["DOID:0080600"],  # Lopinavir
        "DB14761": ["DOID:0080600"],  # Molnupiravir
        "DB16691": ["DOID:0080600"],  # Nirmatrelvir
        "DB06817": ["DOID:0080600"],  # Ruxolitinib
        "DB00051": ["DOID:7148"],  # Adalimumab
        "DB00198": ["DOID:552"],  # Oseltamivir
    }
    for drug, diseases in drug_disease_known.items():
        for disease in diseases:
            triples.append((drug, "treats", disease))

    return triples


def main():
    print("=== Building Biomedical Knowledge Graph ===")
    triples = build_kg()

    # Save entity mappings
    entities = {}
    for d_id, d_name in DRUGS:
        entities[d_id] = {"type": "Drug", "name": d_name}
    for g_id, g_name in GENES:
        entities[g_id] = {"type": "Gene", "name": g_name}
    for d_id, d_name in DISEASES:
        entities[d_id] = {"type": "Disease", "name": d_name}
    for p_id, p_name in PATHWAYS:
        entities[p_id] = {"type": "Pathway", "name": p_name}
    for ph_id, ph_name in PHENOTYPES:
        entities[ph_id] = {"type": "Phenotype", "name": ph_name}

    with open(os.path.join(DATA_DIR, "entities.json"), "w") as f:
        json.dump(entities, f, indent=2)

    # Save triples as TSV
    df = pd.DataFrame(triples, columns=["head", "relation", "tail"])
    df.to_csv(os.path.join(DATA_DIR, "triples.tsv"), sep="\t", index=False)

    # Summary statistics
    G = nx.DiGraph()
    for h, r, t in triples:
        G.add_edge(h, t, relation=r)

    relation_counts = df["relation"].value_counts().to_dict()

    stats = {
        "total_entities": len(entities),
        "total_triples": len(triples),
        "entity_types": {
            "Drug": sum(1 for e in entities.values() if e["type"] == "Drug"),
            "Gene": sum(1 for e in entities.values() if e["type"] == "Gene"),
            "Disease": sum(1 for e in entities.values() if e["type"] == "Disease"),
            "Pathway": sum(1 for e in entities.values() if e["type"] == "Pathway"),
            "Phenotype": sum(1 for e in entities.values() if e["type"] == "Phenotype"),
        },
        "relation_types": relation_counts,
        "graph_density": nx.density(G),
        "avg_degree": np.mean([d for _, d in G.degree()]),
        "connected_components": nx.number_weakly_connected_components(G),
    }

    with open(os.path.join(DATA_DIR, "kg_stats.json"), "w") as f:
        json.dump(stats, f, indent=2, default=str)

    print(f"  Entities: {stats['total_entities']}")
    print(f"  Triples: {stats['total_triples']}")
    print(f"  Relations: {relation_counts}")
    print(f"  Graph density: {stats['graph_density']:.4f}")
    print(f"  Avg degree: {stats['avg_degree']:.2f}")
    print(f"  Weakly connected components: {stats['connected_components']}")

    # Log
    log_entry = {
        "timestamp": datetime.now().isoformat(),
        "phase": "data_preparation",
        "event_type": "kg_construction",
        "actor": "co-scientist",
        "skill_or_tool": "01_build_knowledge_graph.py",
        "handoff_out": stats,
        "files_written": ["data/entities.json", "data/triples.tsv", "data/kg_stats.json"],
        "status": "ok",
    }
    with open(os.path.join(LOG_DIR, "process-log.jsonl"), "a") as f:
        f.write(json.dumps(log_entry) + "\n")

    print("=== Knowledge Graph construction complete ===")


if __name__ == "__main__":
    main()
