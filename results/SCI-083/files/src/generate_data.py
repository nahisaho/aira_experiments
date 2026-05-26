"""
Synthetic data generation for metabolomics-microbiome integration study.
Simulates IBD cohort with paired metabolomics and 16S microbiome data.
"""
import numpy as np
import pandas as pd
import os

np.random.seed(42)

N_SAMPLES = 200  # 100 IBD, 100 healthy controls
N_METABOLITES = 150
N_TAXA = 80
N_PATHWAYS = 25

# Sample metadata
groups = np.array(['IBD'] * 100 + ['Control'] * 100)
sample_ids = [f'S{i:03d}' for i in range(N_SAMPLES)]

# Generate taxa names (genus-level)
genera = [
    'Bacteroides', 'Faecalibacterium', 'Roseburia', 'Eubacterium', 'Ruminococcus',
    'Blautia', 'Coprococcus', 'Prevotella', 'Bifidobacterium', 'Lactobacillus',
    'Akkermansia', 'Alistipes', 'Parabacteroides', 'Clostridium', 'Streptococcus',
    'Enterococcus', 'Escherichia', 'Veillonella', 'Dialister', 'Megamonas',
    'Sutterella', 'Bilophila', 'Desulfovibrio', 'Fusobacterium', 'Haemophilus',
    'Klebsiella', 'Proteus', 'Citrobacter', 'Lachnospira', 'Dorea',
    'Oscillospira', 'Butyricicoccus', 'Anaerostipes', 'Collinsella', 'Holdemanella',
    'Methanobrevibacter', 'Christensenella', 'Phascolarctobacterium', 'Acidaminococcus', 'Catenibacterium'
]
taxa_names = []
for g in genera:
    taxa_names.append(g)
    taxa_names.append(f'{g}_sp2')
taxa_names = taxa_names[:N_TAXA]

# Generate metabolite names
metabolite_classes = {
    'SCFA': ['Butyrate', 'Propionate', 'Acetate', 'Valerate', 'Isobutyrate', 'Isovalerate'],
    'BileAcid': ['Deoxycholic_acid', 'Lithocholic_acid', 'Ursodeoxycholic_acid', 'Chenodeoxycholic_acid',
                 'Cholic_acid', 'Taurocholic_acid', 'Glycocholic_acid', 'Taurodeoxycholic_acid'],
    'AminoAcid': ['Tryptophan', 'Phenylalanine', 'Tyrosine', 'Histidine', 'Leucine',
                  'Isoleucine', 'Valine', 'Glutamine', 'Glutamate', 'Alanine',
                  'Serine', 'Threonine', 'Proline', 'Glycine', 'Arginine'],
    'TryptophanMetab': ['Indole', 'Indole_3_acetic_acid', 'Kynurenine', 'Serotonin',
                        'Indole_3_propionic_acid', 'Quinolinic_acid', 'Xanthurenic_acid'],
    'Lipid': ['LysoPC_16_0', 'LysoPC_18_1', 'LysoPE_16_0', 'PC_34_1', 'PE_36_2',
              'Sphingomyelin_d18_1', 'Ceramide_d18_1_16_0', 'DG_34_1', 'TG_52_2'],
    'Vitamin': ['Riboflavin', 'Nicotinamide', 'Pantothenate', 'Biotin', 'Folate', 'Thiamine'],
    'Phenolic': ['Hippuric_acid', 'p_Cresol_sulfate', 'Phenylacetic_acid', 'Indoxyl_sulfate',
                 '4_Hydroxyphenylacetic_acid', '3_4_Dihydroxyphenylacetic_acid'],
}
metabolite_names = []
metabolite_class_labels = []
for cls, mets in metabolite_classes.items():
    for m in mets:
        metabolite_names.append(m)
        metabolite_class_labels.append(cls)
# Pad to N_METABOLITES
while len(metabolite_names) < N_METABOLITES:
    idx = len(metabolite_names)
    metabolite_names.append(f'Unknown_m{idx}')
    metabolite_class_labels.append('Unknown')
metabolite_names = metabolite_names[:N_METABOLITES]
metabolite_class_labels = metabolite_class_labels[:N_METABOLITES]

# --- Microbiome data (compositional, CLR-transformed for analysis) ---
# IBD: reduced Faecalibacterium, Roseburia; increased Escherichia, Fusobacterium
base_abundances = np.random.dirichlet(np.ones(N_TAXA) * 0.5, size=N_SAMPLES)

# Modify IBD samples
for i in range(100):  # IBD
    # Decrease beneficial
    for j, t in enumerate(taxa_names):
        if 'Faecalibacterium' in t or 'Roseburia' in t or 'Coprococcus' in t:
            base_abundances[i, j] *= 0.3
        if 'Escherichia' in t or 'Fusobacterium' in t or 'Klebsiella' in t:
            base_abundances[i, j] *= 3.0
    base_abundances[i] /= base_abundances[i].sum()

# CLR transform
from scipy.stats import gmean
def clr_transform(data):
    data = data + 1e-10
    gm = gmean(data, axis=1).reshape(-1, 1)
    return np.log(data / gm)

taxa_clr = clr_transform(base_abundances)

# --- Metabolomics data ---
# Correlated with microbiome + disease effect
metabolite_data = np.random.randn(N_SAMPLES, N_METABOLITES) * 0.5

# Create correlations between specific taxa and metabolites
# Faecalibacterium -> Butyrate (positive)
butyrate_idx = metabolite_names.index('Butyrate')
faecal_idx = taxa_names.index('Faecalibacterium')
metabolite_data[:, butyrate_idx] += taxa_clr[:, faecal_idx] * 0.7

# Roseburia -> Propionate
prop_idx = metabolite_names.index('Propionate')
ros_idx = taxa_names.index('Roseburia')
metabolite_data[:, prop_idx] += taxa_clr[:, ros_idx] * 0.6

# Escherichia -> Indoxyl_sulfate (positive in IBD)
indoxyl_idx = metabolite_names.index('Indoxyl_sulfate')
ecoli_idx = taxa_names.index('Escherichia')
metabolite_data[:, indoxyl_idx] += taxa_clr[:, ecoli_idx] * 0.5

# Bifidobacterium -> Acetate
acetate_idx = metabolite_names.index('Acetate')
bifido_idx = taxa_names.index('Bifidobacterium')
metabolite_data[:, acetate_idx] += taxa_clr[:, bifido_idx] * 0.5

# Disease effect on tryptophan pathway
trp_idx = metabolite_names.index('Tryptophan')
kyn_idx = metabolite_names.index('Kynurenine')
for i in range(100):
    metabolite_data[i, trp_idx] -= 1.0
    metabolite_data[i, kyn_idx] += 0.8

# Disease effect on bile acids
for i in range(100):
    da_idx = metabolite_names.index('Deoxycholic_acid')
    metabolite_data[i, da_idx] += 0.7

# Add noise
metabolite_data += np.random.randn(N_SAMPLES, N_METABOLITES) * 0.3

# --- Pathway definitions ---
pathway_names = [
    'Butyrate_biosynthesis', 'Propionate_biosynthesis', 'Acetate_biosynthesis',
    'Tryptophan_metabolism', 'Bile_acid_biotransformation', 'Folate_biosynthesis',
    'Indole_metabolism', 'Phenylpropanoid_degradation', 'Sphingolipid_metabolism',
    'Glycerophospholipid_metabolism', 'Branched_chain_amino_acid_metabolism',
    'Aromatic_amino_acid_metabolism', 'Glutamate_metabolism', 'Purine_metabolism',
    'Pyrimidine_metabolism', 'Sulfur_metabolism', 'Methane_metabolism',
    'Nitrogen_metabolism', 'TCA_cycle', 'Glycolysis',
    'Pentose_phosphate_pathway', 'Fatty_acid_biosynthesis', 'Biotin_metabolism',
    'Nicotinate_metabolism', 'Pantothenate_metabolism'
]

# Pathway-metabolite mapping
pathway_metabolite_map = {
    'Butyrate_biosynthesis': ['Butyrate', 'Acetate', 'Glutamate'],
    'Propionate_biosynthesis': ['Propionate', 'Valerate'],
    'Acetate_biosynthesis': ['Acetate', 'Isobutyrate'],
    'Tryptophan_metabolism': ['Tryptophan', 'Indole', 'Kynurenine', 'Serotonin', 'Indole_3_acetic_acid',
                              'Indole_3_propionic_acid', 'Quinolinic_acid', 'Xanthurenic_acid'],
    'Bile_acid_biotransformation': ['Deoxycholic_acid', 'Lithocholic_acid', 'Ursodeoxycholic_acid',
                                    'Chenodeoxycholic_acid', 'Cholic_acid'],
    'Indole_metabolism': ['Indole', 'Indoxyl_sulfate', 'Indole_3_acetic_acid', 'Indole_3_propionic_acid'],
    'Phenylpropanoid_degradation': ['Hippuric_acid', 'p_Cresol_sulfate', 'Phenylacetic_acid',
                                     '4_Hydroxyphenylacetic_acid'],
    'Sphingolipid_metabolism': ['Sphingomyelin_d18_1', 'Ceramide_d18_1_16_0'],
    'Glycerophospholipid_metabolism': ['LysoPC_16_0', 'LysoPC_18_1', 'LysoPE_16_0', 'PC_34_1', 'PE_36_2'],
    'Branched_chain_amino_acid_metabolism': ['Leucine', 'Isoleucine', 'Valine'],
    'Aromatic_amino_acid_metabolism': ['Tryptophan', 'Phenylalanine', 'Tyrosine'],
    'Glutamate_metabolism': ['Glutamine', 'Glutamate', 'Alanine'],
}

# Pathway-taxa mapping
pathway_taxa_map = {
    'Butyrate_biosynthesis': ['Faecalibacterium', 'Roseburia', 'Eubacterium', 'Butyricicoccus', 'Anaerostipes'],
    'Propionate_biosynthesis': ['Bacteroides', 'Prevotella', 'Alistipes', 'Phascolarctobacterium'],
    'Acetate_biosynthesis': ['Bifidobacterium', 'Blautia', 'Ruminococcus'],
    'Tryptophan_metabolism': ['Escherichia', 'Lactobacillus', 'Clostridium', 'Bacteroides'],
    'Bile_acid_biotransformation': ['Clostridium', 'Bacteroides', 'Eubacterium', 'Lactobacillus'],
    'Indole_metabolism': ['Escherichia', 'Clostridium', 'Bacteroides'],
    'Phenylpropanoid_degradation': ['Clostridium', 'Blautia', 'Dorea'],
}

# Save all data
os.makedirs('data', exist_ok=True)

taxa_df = pd.DataFrame(taxa_clr, index=sample_ids, columns=taxa_names)
taxa_df.to_csv('data/taxa_clr.csv')

met_df = pd.DataFrame(metabolite_data, index=sample_ids, columns=metabolite_names)
met_df.to_csv('data/metabolites.csv')

meta_df = pd.DataFrame({'SampleID': sample_ids, 'Group': groups})
meta_df.to_csv('data/metadata.csv', index=False)

rel_df = pd.DataFrame(base_abundances, index=sample_ids, columns=taxa_names)
rel_df.to_csv('data/taxa_relative.csv')

# Save pathway maps as JSON
import json
with open('data/pathway_metabolite_map.json', 'w') as f:
    json.dump(pathway_metabolite_map, f, indent=2)
with open('data/pathway_taxa_map.json', 'w') as f:
    json.dump(pathway_taxa_map, f, indent=2)

met_class_df = pd.DataFrame({'Metabolite': metabolite_names, 'Class': metabolite_class_labels})
met_class_df.to_csv('data/metabolite_classes.csv', index=False)

print(f"Generated data: {N_SAMPLES} samples, {N_METABOLITES} metabolites, {N_TAXA} taxa")
print(f"Groups: {np.unique(groups, return_counts=True)}")
