"""
Biomedical Knowledge Graph Construction
Sources: DrugBank, DisGeNET, STRING, CTD (simulated with real biological structure)
Entities: Drugs, Diseases, Genes, Pathways, Phenotypes
"""

import json
import os
import random
import time
from pathlib import Path

import networkx as nx
import pandas as pd
import numpy as np

random.seed(42)
np.random.seed(42)

BASE = Path(__file__).parent.parent
DATA_DIR = BASE / "data"
LOG_FILE = BASE / "logs" / "process-log.jsonl"
DATA_DIR.mkdir(exist_ok=True)
(BASE / "logs").mkdir(exist_ok=True)


def log_event(event_type, details):
    entry = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "phase": "EXECUTE",
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": "co-scientist-drug-repurposing",
        **details,
    }
    with open(LOG_FILE, "a") as f:
        f.write(json.dumps(entry) + "\n")


log_event("run_started", {"script": "01_build_knowledge_graph.py"})

# ─────────────────────────────────────────────
# 1. Entity Definitions
# ─────────────────────────────────────────────

# Drugs (DrugBank-inspired)
DRUGS = {
    "DB00001": "Lepirudin",
    "DB00002": "Cetuximab",
    "DB00006": "Bivalirudin",
    "DB00007": "Leuprolide",
    "DB00014": "Goserelin",
    "DB00018": "Finasteride",
    "DB00026": "Basiliximab",
    "DB00028": "Human-Immunoglobulin",
    "DB00033": "Interferon-gamma",
    "DB00034": "Interferon-alfa",
    "DB00041": "Aldesleukin",
    "DB00058": "Alglucerase",
    "DB00059": "Pegaspargase",
    "DB00062": "Human-Serum-Albumin",
    "DB00067": "Vasopressin",
    "DB00072": "Trastuzumab",
    "DB00080": "Daptomycin",
    "DB00082": "Pegvisomant",
    "DB00083": "Botulinum-Toxin-A",
    "DB00091": "Cyclosporine",
    "DB00099": "Filgrastim",
    "DB00108": "Natalizumab",
    "DB00109": "Abciximab",
    "DB00111": "Human-Chorionic-Gonadotropin",
    "DB00114": "Pyridoxal-Phosphate",
    # COVID-related drugs
    "DB14443": "Remdesivir",
    "DB00207": "Azithromycin",
    "DB00795": "Sulfasalazine",
    "DB00993": "Azathioprine",
    "DB01076": "Atorvastatin",
    "DB01048": "Abacavir",
    "DB00783": "Estradiol",
    "DB00741": "Hydrocortisone",
    "DB01234": "Dexamethasone",
    "DB01029": "Irbesartan",
    "DB00177": "Valsartan",
    "DB00196": "Fluconazole",
    "DB01222": "Budesonide",
    "DB00991": "Oxaliplatin",
    "DB01232": "Saquinavir",
    "DB00932": "Lopinavir",
    "DB01601": "Lopinavir-Ritonavir",
    "DB01258": "Aliskiren",
    "DB00758": "Clopidogrel",
    "DB00001X": "Baricitinib",
    "DB00002X": "Tocilizumab",
    "DB00003X": "Sarilumab",
    "DB00004X": "Anakinra",
    "DB00005X": "Hydroxychloroquine",
    "DB00006X": "Ivermectin",
    "DB00007X": "Favipiravir",
    "DB00008X": "Molnupiravir",
    "DB00009X": "Paxlovid",
    "DB00010X": "Colchicine",
}

# Diseases (OMIM/MeSH-inspired)
DISEASES = {
    "MESH:D000755": "Anemia-Sickle-Cell",
    "MESH:D001249": "Asthma",
    "MESH:D001943": "Breast-Neoplasms",
    "MESH:D003550": "Cystic-Fibrosis",
    "MESH:D003920": "Diabetes-Mellitus",
    "MESH:D004827": "Epilepsy",
    "MESH:D006333": "Heart-Failure",
    "MESH:D006973": "Hypertension",
    "MESH:D007938": "Leukemia",
    "MESH:D008223": "Lymphoma",
    "MESH:D010300": "Parkinson-Disease",
    "MESH:D011289": "Prader-Willi-Syndrome",
    "MESH:D012174": "Retinitis-Pigmentosa",
    "MESH:D012559": "Schizophrenia",
    "MESH:D012600": "Scoliosis",
    "MESH:D013593": "Synovitis",
    "MESH:D014376": "Tuberculosis",
    "MESH:D015179": "Colorectal-Neoplasms",
    "MESH:D016180": "Lentiviruses",
    "MESH:D016360": "Molecular-Mimicry",
    # COVID-19
    "MESH:D000086382": "COVID-19",
    "MESH:D011014": "Pneumonia",
    "MESH:D012128": "Respiratory-Distress-Syndrome",
    "MESH:D018352": "Coronavirus-Infections",
    "MESH:D045169": "SARS",
    "MESH:D007249": "Inflammation",
    "MESH:D016638": "Critical-Illness",
    "MESH:D013927": "Thrombosis",
    "MESH:D012769": "Shock-Septic",
    "MESH:D011655": "Pulmonary-Embolism",
}

# Genes (HGNC-inspired)
GENES = {
    "HGNC:11892": "TNF",
    "HGNC:6018": "IL6",
    "HGNC:6116": "IL1B",
    "HGNC:7025": "MYC",
    "HGNC:9829": "TP53",
    "HGNC:3236": "EGFR",
    "HGNC:11998": "VEGFA",
    "HGNC:9588": "PTGS2",
    "HGNC:5465": "HMOX1",
    "HGNC:8975": "ACE2",  # COVID key gene
    "HGNC:15509": "TMPRSS2",  # COVID key gene
    "HGNC:7486": "NFKB1",
    "HGNC:10671": "STAT3",
    "HGNC:11876": "TLR4",
    "HGNC:4079": "FURIN",  # COVID related
    "HGNC:9948": "RELA",
    "HGNC:7127": "MMP9",
    "HGNC:3356": "F2",
    "HGNC:4882": "IFNG",
    "HGNC:5962": "IL10",
    "HGNC:6018": "IL6",
    "HGNC:6711": "CXCL8",
    "HGNC:1499": "CCL2",
    "HGNC:4128": "CSF2",
    "HGNC:11775": "TGFB1",
}

# Pathways (Reactome-inspired)
PATHWAYS = {
    "R-HSA-168928": "DDX58-IFIH1-mediated-induction-of-IFN",
    "R-HSA-1280215": "Cytokine-Signaling-in-Immune-System",
    "R-HSA-162582": "Signal-Transduction",
    "R-HSA-6798695": "Innate-Immune-System",
    "R-HSA-74160": "Gene-Expression",
    "R-HSA-5633007": "Regulation-of-TP53-Activity",
    "R-HSA-3232118": "SUMOylation-of-Transcription-Cofactors",
    "R-HSA-168256": "Immune-System",
    "R-HSA-109582": "Hemostasis",
    "R-HSA-392499": "Metabolism-of-proteins",
    "R-HSA-5357801": "Programmed-Cell-Death",
    "R-HSA-69306": "DNA-Repair",
    "R-HSA-400508": "Hormone-ligand-binding-receptors",
    "R-HSA-3700989": "Transcriptional-Regulation-by-Small-Molecules",
    "R-HSA-9612973": "Autophagy",
}

# Phenotypes (HPO-inspired)
PHENOTYPES = {
    "HP:0001945": "Fever",
    "HP:0002090": "Pneumonia",
    "HP:0002099": "Asthma",
    "HP:0001635": "Congestive-Heart-Failure",
    "HP:0000822": "Hypertension",
    "HP:0002018": "Nausea",
    "HP:0001744": "Splenomegaly",
    "HP:0002017": "Nausea-Vomiting",
    "HP:0011458": "Abdominal-symptom",
    "HP:0001988": "Recurrent-Hypoglycemia",
    "HP:0000716": "Depression",
    "HP:0012531": "Pain",
    "HP:0002098": "Respiratory-Distress",
    "HP:0001297": "Stroke",
    "HP:0001907": "Thromboembolism",
}

# ─────────────────────────────────────────────
# 2. Triple Construction
# ─────────────────────────────────────────────

triples = []


def add_triple(h, r, t):
    triples.append({"head": h, "relation": r, "tail": t})


# Drug–Disease (treats, contraindicated_in, investigated_for)
drug_disease_treats = [
    ("DB00091", "treats", "MESH:D003920"),
    ("DB00091", "treats", "MESH:D008223"),
    ("DB00033", "treats", "MESH:D007938"),
    ("DB00034", "treats", "MESH:D016180"),
    ("DB00072", "treats", "MESH:D001943"),
    ("DB00002", "treats", "MESH:D015179"),
    ("DB01076", "treats", "MESH:D006973"),
    ("DB01076", "treats", "MESH:D003920"),
    ("DB00741", "treats", "MESH:D001249"),
    ("DB01234", "treats", "MESH:D001249"),
    ("DB01234", "treats", "MESH:D007249"),
    ("DB01234", "treats", "COVID-19-treatment", ),  # known treatment
    ("DB00932", "investigated_for", "MESH:D018352"),
    ("DB14443", "treats", "MESH:D000086382"),  # Remdesivir → COVID-19
    ("DB14443", "investigated_for", "MESH:D045169"),
    ("DB00001X", "treats", "MESH:D000086382"),  # Baricitinib → COVID-19
    ("DB00002X", "treats", "MESH:D000086382"),  # Tocilizumab → COVID-19
    ("DB01234", "treats", "MESH:D000086382"),   # Dexamethasone → COVID-19
    ("DB00758", "investigated_for", "MESH:D013927"),
    ("DB00010X", "investigated_for", "MESH:D000086382"),  # Colchicine
    ("DB00009X", "treats", "MESH:D000086382"),   # Paxlovid
    ("DB00008X", "investigated_for", "MESH:D000086382"),
    ("DB00795", "investigated_for", "MESH:D007249"),
    ("DB00993", "treats", "MESH:D013593"),
    ("DB00067", "treats", "MESH:D006333"),
    ("DB00018", "treats", "MESH:D011289"),
]

for h, r, t in drug_disease_treats:
    if t in DISEASES or t.startswith("COVID"):
        if h in DRUGS:
            add_triple(h, r, t if t in DISEASES else "MESH:D000086382")

# Drug–Gene (targets, inhibits, activates, upregulates, downregulates)
drug_gene_rels = [
    ("DB14443", "inhibits", "HGNC:8975"),
    ("DB14443", "inhibits", "HGNC:7486"),
    ("DB00001X", "inhibits", "HGNC:10671"),  # Baricitinib → JAK/STAT3
    ("DB00002X", "inhibits", "HGNC:6018"),   # Tocilizumab → IL6
    ("DB00003X", "inhibits", "HGNC:6018"),   # Sarilumab → IL6
    ("DB00004X", "inhibits", "HGNC:6116"),   # Anakinra → IL1B
    ("DB01234", "downregulates", "HGNC:11892"),  # Dex → TNF
    ("DB01234", "downregulates", "HGNC:6018"),
    ("DB01076", "inhibits", "HGNC:9588"),
    ("DB01076", "downregulates", "HGNC:11892"),
    ("DB00072", "targets", "HGNC:3236"),
    ("DB00002", "targets", "HGNC:3236"),
    ("DB00091", "inhibits", "HGNC:7486"),
    ("DB00091", "inhibits", "HGNC:10671"),
    ("DB00795", "inhibits", "HGNC:9588"),
    ("DB00795", "inhibits", "HGNC:11892"),
    ("DB10029", "inhibits", "HGNC:8975"),   # ACE inhibitor
    ("DB01258", "inhibits", "HGNC:8975"),
    ("DB01029", "inhibits", "HGNC:8975"),
    ("DB00177", "inhibits", "HGNC:8975"),
    ("DB00933", "targets", "HGNC:4079"),   # FURIN inhibitor
    ("DB00010X", "inhibits", "HGNC:11892"),  # Colchicine → TNF
    ("DB00010X", "inhibits", "HGNC:6116"),
    ("DB00005X", "inhibits", "HGNC:11892"),
    ("DB00005X", "inhibits", "HGNC:7486"),
    ("DB00007X", "inhibits", "HGNC:15509"),  # Favipiravir → TMPRSS2
]

for h, r, t in drug_gene_rels:
    if h in DRUGS and t in GENES:
        add_triple(h, r, t)

# Gene–Disease (associated_with, causes, biomarker_of)
gene_disease_rels = [
    ("HGNC:11892", "associated_with", "MESH:D007249"),  # TNF → Inflammation
    ("HGNC:6018", "associated_with", "MESH:D000086382"),  # IL6 → COVID
    ("HGNC:6018", "associated_with", "MESH:D007249"),
    ("HGNC:6116", "associated_with", "MESH:D006333"),
    ("HGNC:6116", "associated_with", "MESH:D007249"),
    ("HGNC:10671", "associated_with", "MESH:D000086382"),  # STAT3 → COVID
    ("HGNC:8975", "associated_with", "MESH:D000086382"),   # ACE2 → COVID (receptor)
    ("HGNC:15509", "associated_with", "MESH:D000086382"),  # TMPRSS2 → COVID
    ("HGNC:4079", "associated_with", "MESH:D000086382"),   # FURIN → COVID
    ("HGNC:7486", "associated_with", "MESH:D000086382"),   # NFKB1 → COVID
    ("HGNC:7486", "associated_with", "MESH:D007938"),
    ("HGNC:9829", "associated_with", "MESH:D001943"),
    ("HGNC:9829", "associated_with", "MESH:D007938"),
    ("HGNC:3236", "associated_with", "MESH:D001943"),
    ("HGNC:11998", "associated_with", "MESH:D001943"),
    ("HGNC:7025", "associated_with", "MESH:D008223"),
    ("HGNC:9588", "associated_with", "MESH:D007249"),
    ("HGNC:7127", "associated_with", "MESH:D012128"),  # MMP9 → ARDS
    ("HGNC:3356", "associated_with", "MESH:D013927"),  # F2 → Thrombosis
    ("HGNC:4882", "associated_with", "MESH:D000086382"),  # IFNG
    ("HGNC:6711", "associated_with", "MESH:D000086382"),  # CXCL8
    ("HGNC:1499", "associated_with", "MESH:D000086382"),  # CCL2
    ("HGNC:11775", "associated_with", "MESH:D012128"),
]

for h, r, t in gene_disease_rels:
    if h in GENES and t in DISEASES:
        add_triple(h, r, t)

# Gene–Pathway (participates_in, regulates)
gene_pathway_rels = [
    ("HGNC:11892", "participates_in", "R-HSA-1280215"),
    ("HGNC:6018", "participates_in", "R-HSA-1280215"),
    ("HGNC:6116", "participates_in", "R-HSA-1280215"),
    ("HGNC:7486", "participates_in", "R-HSA-6798695"),
    ("HGNC:10671", "participates_in", "R-HSA-1280215"),
    ("HGNC:9829", "participates_in", "R-HSA-5633007"),
    ("HGNC:9829", "regulates", "R-HSA-5357801"),
    ("HGNC:3236", "participates_in", "R-HSA-162582"),
    ("HGNC:11998", "participates_in", "R-HSA-162582"),
    ("HGNC:8975", "participates_in", "R-HSA-109582"),
    ("HGNC:15509", "participates_in", "R-HSA-6798695"),
    ("HGNC:4079", "participates_in", "R-HSA-392499"),
    ("HGNC:4882", "participates_in", "R-HSA-168928"),
    ("HGNC:7127", "participates_in", "R-HSA-109582"),
    ("HGNC:3356", "participates_in", "R-HSA-109582"),
    ("HGNC:11775", "participates_in", "R-HSA-162582"),
    ("HGNC:5465", "participates_in", "R-HSA-9612973"),
    ("HGNC:9588", "participates_in", "R-HSA-74160"),
]

for h, r, t in gene_pathway_rels:
    if h in GENES and t in PATHWAYS:
        add_triple(h, r, t)

# Disease–Phenotype (has_phenotype)
disease_phenotype = [
    ("MESH:D000086382", "has_phenotype", "HP:0001945"),  # COVID → Fever
    ("MESH:D000086382", "has_phenotype", "HP:0002090"),  # COVID → Pneumonia
    ("MESH:D000086382", "has_phenotype", "HP:0002098"),  # COVID → Resp Distress
    ("MESH:D000086382", "has_phenotype", "HP:0001907"),  # COVID → Thromboembolism
    ("MESH:D007249", "has_phenotype", "HP:0001945"),
    ("MESH:D012128", "has_phenotype", "HP:0002098"),
    ("MESH:D001249", "has_phenotype", "HP:0002099"),
    ("MESH:D006333", "has_phenotype", "HP:0001635"),
    ("MESH:D006973", "has_phenotype", "HP:0000822"),
    ("MESH:D013927", "has_phenotype", "HP:0001907"),
    ("MESH:D011655", "has_phenotype", "HP:0001907"),
]

for h, r, t in disease_phenotype:
    if h in DISEASES and t in PHENOTYPES:
        add_triple(h, r, t)

# Drug–Pathway (modulates)
drug_pathway = [
    ("DB14443", "modulates", "R-HSA-6798695"),
    ("DB00001X", "modulates", "R-HSA-1280215"),
    ("DB00002X", "modulates", "R-HSA-1280215"),
    ("DB01234", "modulates", "R-HSA-6798695"),
    ("DB01076", "modulates", "R-HSA-162582"),
    ("DB00091", "modulates", "R-HSA-1280215"),
    ("DB00010X", "modulates", "R-HSA-6798695"),
    ("DB00005X", "modulates", "R-HSA-6798695"),
]

for h, r, t in drug_pathway:
    if h in DRUGS and t in PATHWAYS:
        add_triple(h, r, t)

# Gene–Gene (interacts_with via STRING)
gene_gene = [
    ("HGNC:11892", "interacts_with", "HGNC:6018"),
    ("HGNC:11892", "interacts_with", "HGNC:6116"),
    ("HGNC:6018", "interacts_with", "HGNC:10671"),
    ("HGNC:7486", "interacts_with", "HGNC:11892"),
    ("HGNC:7486", "interacts_with", "HGNC:6018"),
    ("HGNC:8975", "interacts_with", "HGNC:15509"),
    ("HGNC:9829", "interacts_with", "HGNC:7025"),
    ("HGNC:9829", "interacts_with", "HGNC:5633007"),
    ("HGNC:4882", "interacts_with", "HGNC:10671"),
    ("HGNC:4882", "interacts_with", "HGNC:11892"),
    ("HGNC:11775", "interacts_with", "HGNC:7127"),
    ("HGNC:3356", "interacts_with", "HGNC:7127"),
]

for h, r, t in gene_gene:
    if h in GENES and t in GENES:
        add_triple(h, r, t)

# Save triples
df_triples = pd.DataFrame(triples)
df_triples.to_csv(DATA_DIR / "kg_triples.tsv", sep="\t", index=False)
print(f"Total triples: {len(df_triples)}")
print(df_triples["relation"].value_counts())

# ─────────────────────────────────────────────
# 3. Build NetworkX graph for analysis
# ─────────────────────────────────────────────
G = nx.MultiDiGraph()
for _, row in df_triples.iterrows():
    G.add_edge(row["head"], row["tail"], relation=row["relation"])

# Node types
all_entities = {}
for k, v in DRUGS.items():
    all_entities[k] = {"name": v, "type": "drug"}
for k, v in DISEASES.items():
    all_entities[k] = {"name": v, "type": "disease"}
for k, v in GENES.items():
    all_entities[k] = {"name": v, "type": "gene"}
for k, v in PATHWAYS.items():
    all_entities[k] = {"name": v, "type": "pathway"}
for k, v in PHENOTYPES.items():
    all_entities[k] = {"name": v, "type": "phenotype"}

for node_id, attrs in all_entities.items():
    if node_id in G.nodes:
        G.nodes[node_id].update(attrs)

entity_df = pd.DataFrame([{"id": k, **v} for k, v in all_entities.items()])
entity_df.to_csv(DATA_DIR / "kg_entities.csv", index=False)

# Stats
node_types = {}
for n, d in G.nodes(data=True):
    t = d.get("type", "unknown")
    node_types[t] = node_types.get(t, 0) + 1

relation_counts = df_triples["relation"].value_counts().to_dict()

stats = {
    "total_triples": len(df_triples),
    "total_nodes": G.number_of_nodes(),
    "total_edges": G.number_of_edges(),
    "node_types": node_types,
    "relation_counts": relation_counts,
    "density": nx.density(G),
    "unique_relations": df_triples["relation"].nunique(),
}

with open(DATA_DIR / "kg_stats.json", "w") as f:
    json.dump(stats, f, indent=2)

print("\n=== Knowledge Graph Statistics ===")
for k, v in stats.items():
    print(f"  {k}: {v}")

log_event("file_written", {
    "files_written": ["data/kg_triples.tsv", "data/kg_entities.csv", "data/kg_stats.json"],
    "status": "ok",
    "details": stats
})

print("\n[✓] Knowledge graph built and saved.")
