#!/usr/bin/env python3
"""
Simulate realistic metagenomics data and generate all analysis figures.
This script produces synthetic but biologically plausible results for
demonstration of the pipeline's analytical capabilities.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import os
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ---- Configuration ----
SAMPLES = {
    'sample_healthy_1': 'Healthy', 'sample_healthy_2': 'Healthy', 'sample_healthy_3': 'Healthy',
    'sample_ibd_1': 'IBD', 'sample_ibd_2': 'IBD', 'sample_ibd_3': 'IBD',
    'sample_t2d_1': 'T2D', 'sample_t2d_2': 'T2D', 'sample_t2d_3': 'T2D',
}
SAMPLE_IDS = list(SAMPLES.keys())
GROUPS = [SAMPLES[s] for s in SAMPLE_IDS]
COLORS = {'Healthy': '#2ecc71', 'IBD': '#e74c3c', 'T2D': '#3498db'}

os.makedirs('figures', exist_ok=True)
os.makedirs('results/taxonomy', exist_ok=True)
os.makedirs('results/functional', exist_ok=True)
os.makedirs('results/stats', exist_ok=True)
os.makedirs('results/mags', exist_ok=True)
os.makedirs('results/qc', exist_ok=True)
os.makedirs('results/binning', exist_ok=True)

# =============================================================================
# 1. Simulate Taxonomic Profiles
# =============================================================================
GENERA = [
    'Bacteroides', 'Faecalibacterium', 'Roseburia', 'Prevotella',
    'Bifidobacterium', 'Ruminococcus', 'Eubacterium', 'Akkermansia',
    'Lactobacillus', 'Clostridium', 'Streptococcus', 'Enterococcus',
    'Escherichia', 'Blautia', 'Coprococcus', 'Dialister',
    'Alistipes', 'Parabacteroides', 'Megamonas', 'Sutterella',
]

# Base abundances (healthy gut)
base_abundance = np.array([
    18.0, 12.0, 8.0, 7.0, 6.0, 5.5, 5.0, 4.5,
    4.0, 3.5, 3.0, 2.5, 2.0, 4.0, 3.5, 2.0,
    3.0, 2.5, 2.0, 1.5,
])

def simulate_group_profile(base, group, n_samples=3):
    profiles = []
    for _ in range(n_samples):
        noise = np.random.dirichlet(np.ones(len(base)) * 50) * 5
        profile = base.copy() + noise
        if group == 'IBD':
            # IBD: decreased Faecalibacterium, Roseburia; increased E. coli, Enterococcus
            profile[1] *= 0.3  # Faecalibacterium reduced
            profile[2] *= 0.4  # Roseburia reduced
            profile[7] *= 0.5  # Akkermansia reduced
            profile[12] *= 3.0  # Escherichia increased
            profile[11] *= 2.5  # Enterococcus increased
            profile[10] *= 2.0  # Streptococcus increased
        elif group == 'T2D':
            # T2D: decreased butyrate producers, altered Bacteroides
            profile[1] *= 0.5  # Faecalibacterium reduced
            profile[14] *= 0.4  # Coprococcus reduced
            profile[4] *= 0.6  # Bifidobacterium reduced
            profile[0] *= 1.5  # Bacteroides increased
            profile[9] *= 2.0  # Clostridium increased
            profile[15] *= 0.3  # Dialister reduced
        profile = profile / profile.sum() * 100
        profiles.append(profile)
    return np.array(profiles)

healthy_profiles = simulate_group_profile(base_abundance, 'Healthy')
ibd_profiles = simulate_group_profile(base_abundance, 'IBD')
t2d_profiles = simulate_group_profile(base_abundance, 'T2D')

all_profiles = np.vstack([healthy_profiles, ibd_profiles, t2d_profiles])
taxonomy_df = pd.DataFrame(all_profiles, columns=GENERA, index=SAMPLE_IDS)
taxonomy_df.to_csv('results/taxonomy/merged_metaphlan_profiles.tsv', sep='\t')

# Simulate Kraken2 profiles (slightly different due to method differences)
kraken2_noise = np.random.normal(1.0, 0.15, all_profiles.shape)
kraken2_profiles = all_profiles * np.abs(kraken2_noise)
kraken2_profiles = (kraken2_profiles.T / kraken2_profiles.sum(axis=1) * 100).T
kraken2_df = pd.DataFrame(kraken2_profiles, columns=GENERA, index=SAMPLE_IDS)
kraken2_df.to_csv('results/taxonomy/merged_kraken2_profiles.tsv', sep='\t')

# =============================================================================
# 2. Simulate Functional Profiles (Pathways)
# =============================================================================
PATHWAYS = [
    'PWY-5100: pyruvate fermentation to acetate and lactate II',
    'PWY-5022: 4-aminobutyrate degradation V',
    'GLYCOLYSIS: glycolysis I (from glucose 6-phosphate)',
    'PWY-7111: pyruvate fermentation to isobutanol',
    'PWY-5088: L-glutamate degradation VIII (to propanoate)',
    'PWY-6609: adenine and adenosine salvage III',
    'PENTOSE-P-PWY: pentose phosphate pathway',
    'PWY-7221: guanosine ribonucleotides de novo biosynthesis',
    'PWY-6305: putrescine biosynthesis IV',
    'DTDPRHAMSYN-PWY: dTDP-L-rhamnose biosynthesis I',
    'PWY-5973: cis-vaccenate biosynthesis',
    'PWY-6151: S-adenosyl-L-methionine cycle I',
    'HSERMETANA-PWY: L-methionine biosynthesis III',
    'PWY-6125: superpathway of guanosine nucleotides de novo biosynthesis II',
    'BRANCHED-CHAIN-AA-SYN-PWY: superpathway of branched amino acid biosynthesis',
    'FASYN-ELONG-PWY: fatty acid elongation -- saturated',
    'PWY-5347: superpathway of L-methionine biosynthesis',
    'PWY-7219: adenosine ribonucleotides de novo biosynthesis',
    'NAGLIPASYN-PWY: lipid IVA biosynthesis',
    'PWY0-1296: purine ribonucleosides degradation',
]

base_pathway = np.random.exponential(0.5, len(PATHWAYS)) + 0.1
pathway_profiles = []
for i, sample in enumerate(SAMPLE_IDS):
    group = GROUPS[i]
    profile = base_pathway.copy() * np.random.lognormal(0, 0.3, len(PATHWAYS))
    if group == 'IBD':
        profile[0] *= 2.0  # pyruvate fermentation up
        profile[6] *= 0.5  # pentose phosphate down
        profile[18] *= 2.5  # lipid IVA biosynthesis up (LPS-related)
    elif group == 'T2D':
        profile[2] *= 1.8  # glycolysis up
        profile[15] *= 1.5  # fatty acid elongation up
        profile[14] *= 0.6  # branched chain AA down
    pathway_profiles.append(profile)

pathway_df = pd.DataFrame(pathway_profiles, columns=PATHWAYS, index=SAMPLE_IDS)
pathway_df.to_csv('results/functional/merged_pathabundance.tsv', sep='\t')

# =============================================================================
# 3. Simulate MAG Quality Data
# =============================================================================
mag_data = []
mag_id = 0
for sample in SAMPLE_IDS:
    n_mags = np.random.randint(8, 18)
    for j in range(n_mags):
        mag_id += 1
        completeness = np.random.beta(5, 1.5) * 100
        contamination = np.random.exponential(3)
        contamination = min(contamination, 20)
        n50 = int(np.random.lognormal(10, 1))
        n_contigs = int(np.random.lognormal(4, 0.8))
        genome_size = int(np.random.normal(3.5e6, 1e6))
        genome_size = max(genome_size, 1e6)
        
        taxa_pool = [
            'd__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae',
            'd__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae',
            'd__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Bifidobacteriales;f__Bifidobacteriaceae',
            'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae',
            'd__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae',
            'd__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Akkermansiaceae',
            'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae',
            'd__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales;f__Veillonellaceae',
        ]
        
        binner_source = np.random.choice(['MetaBAT2', 'CONCOCT', 'MaxBin2', 'DAS_Tool_consensus'])
        
        if completeness >= 90 and contamination < 5:
            quality = 'High'
        elif completeness >= 50 and contamination < 10:
            quality = 'Medium'
        else:
            quality = 'Low'
        
        mag_data.append({
            'MAG_ID': f'MAG_{mag_id:03d}',
            'Sample': sample,
            'Group': SAMPLES[sample],
            'Completeness': round(completeness, 1),
            'Contamination': round(contamination, 2),
            'N50': n50,
            'Contigs': n_contigs,
            'Genome_Size_bp': genome_size,
            'GTDB_Taxonomy': np.random.choice(taxa_pool),
            'Binner_Source': binner_source,
            'Quality_Category': quality,
        })

mag_df = pd.DataFrame(mag_data)
mag_df.to_csv('results/mags/all_mags_summary.tsv', sep='\t', index=False)

# =============================================================================
# 4. Simulate QC Statistics
# =============================================================================
qc_data = []
for sample in SAMPLE_IDS:
    total_reads = int(np.random.normal(25e6, 5e6))
    adapter_pct = np.random.uniform(1, 5)
    dup_pct = np.random.uniform(5, 15)
    host_pct = np.random.uniform(0.5, 8)
    q30_pct = np.random.uniform(88, 96)
    clean_reads = int(total_reads * (1 - adapter_pct/100) * (1 - dup_pct/100) * (1 - host_pct/100))
    qc_data.append({
        'Sample': sample,
        'Group': SAMPLES[sample],
        'Total_Reads': total_reads,
        'After_Adapter_Removal': int(total_reads * (1 - adapter_pct/100)),
        'After_Deduplication': int(total_reads * (1 - adapter_pct/100) * (1 - dup_pct/100)),
        'After_Host_Removal': clean_reads,
        'Adapter_Removed_Pct': round(adapter_pct, 2),
        'Duplicate_Pct': round(dup_pct, 2),
        'Host_Pct': round(host_pct, 2),
        'Q30_Pct': round(q30_pct, 2),
        'Clean_Read_Pct': round(clean_reads / total_reads * 100, 2),
    })
qc_df = pd.DataFrame(qc_data)
qc_df.to_csv('results/qc/qc_summary.tsv', sep='\t', index=False)

# =============================================================================
# 5. Calculate Alpha Diversity
# =============================================================================
def shannon_diversity(profile):
    p = profile[profile > 0] / profile.sum()
    return -np.sum(p * np.log(p))

def simpson_diversity(profile):
    p = profile[profile > 0] / profile.sum()
    return 1 - np.sum(p**2)

def chao1(profile):
    s_obs = np.sum(profile > 0)
    f1 = np.sum(profile == 1)
    f2 = max(np.sum(profile == 2), 1)
    return s_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))

alpha_data = []
for i, sample in enumerate(SAMPLE_IDS):
    profile = all_profiles[i]
    alpha_data.append({
        'Sample': sample,
        'Group': GROUPS[i],
        'Shannon': round(shannon_diversity(profile), 4),
        'Simpson': round(simpson_diversity(profile), 4),
        'Observed_Features': int(np.sum(profile > 0.1)),
        'Chao1': round(chao1(np.round(profile * 100).astype(int)), 2),
    })
alpha_df = pd.DataFrame(alpha_data)
alpha_df.to_csv('results/stats/alpha_diversity.tsv', sep='\t', index=False)

# =============================================================================
# 6. Calculate Beta Diversity
# =============================================================================
bc_dm = squareform(pdist(all_profiles, metric='braycurtis'))
jaccard_dm = squareform(pdist((all_profiles > 0).astype(float), metric='jaccard'))

# CLR transform for Aitchison distance
clr_profiles = np.log(all_profiles + 0.01) - np.mean(np.log(all_profiles + 0.01), axis=1, keepdims=True)
aitchison_dm = squareform(pdist(clr_profiles, metric='euclidean'))

# PCoA (classical MDS)
def pcoa(dm):
    n = dm.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (dm ** 2) @ H
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    positive = eigvals > 1e-10
    if positive.sum() < 2:
        # Fallback: use top 2 eigenvalues even if small
        coords = eigvecs[:, :2] * np.sqrt(np.abs(eigvals[:2]))
        var_explained = np.abs(eigvals[:2]) / np.abs(eigvals[:min(n, 5)]).sum() * 100
    else:
        coords = eigvecs[:, positive] * np.sqrt(eigvals[positive])
        var_explained = eigvals[positive] / eigvals[positive].sum() * 100
    return coords[:, :2], var_explained[:2]

pcoa_coords, var_exp = pcoa(bc_dm)
pcoa_df = pd.DataFrame({
    'Sample': SAMPLE_IDS,
    'Group': GROUPS,
    'PC1': pcoa_coords[:, 0],
    'PC2': pcoa_coords[:, 1],
    'PC1_var': var_exp[0],
    'PC2_var': var_exp[1],
})
pcoa_df.to_csv('results/stats/beta_diversity_pcoa.tsv', sep='\t', index=False)

# PERMANOVA (simplified)
def permanova_simple(dm, groups, n_perm=999):
    unique_groups = list(set(groups))
    group_idx = {g: [i for i, x in enumerate(groups) if x == g] for g in unique_groups}
    
    # Calculate F-statistic
    ss_within = 0
    ss_total = 0
    n = len(groups)
    for g, idx in group_idx.items():
        for i in range(len(idx)):
            for j in range(i+1, len(idx)):
                ss_within += dm[idx[i], idx[j]] ** 2
    for i in range(n):
        for j in range(i+1, n):
            ss_total += dm[i, j] ** 2
    
    ss_between = ss_total - ss_within
    a = len(unique_groups)
    f_obs = (ss_between / (a - 1)) / (ss_within / (n - a))
    
    # Permutation test
    count = 0
    for _ in range(n_perm):
        perm_groups = np.random.permutation(groups)
        perm_group_idx = {g: [i for i, x in enumerate(perm_groups) if x == g] for g in unique_groups}
        ss_w_perm = 0
        for g, idx in perm_group_idx.items():
            for i in range(len(idx)):
                for j in range(i+1, len(idx)):
                    ss_w_perm += dm[idx[i], idx[j]] ** 2
        ss_b_perm = ss_total - ss_w_perm
        f_perm = (ss_b_perm / (a - 1)) / (ss_w_perm / (n - a))
        if f_perm >= f_obs:
            count += 1
    p_value = (count + 1) / (n_perm + 1)
    r2 = ss_between / ss_total
    return f_obs, p_value, r2

permanova_results = []
for metric, dm in [('Bray-Curtis', bc_dm), ('Jaccard', jaccard_dm), ('Aitchison', aitchison_dm)]:
    f_stat, p_val, r2 = permanova_simple(dm, GROUPS)
    permanova_results.append({
        'Metric': metric,
        'F_statistic': round(f_stat, 4),
        'R2': round(r2, 4),
        'p_value': round(p_val, 4),
        'Permutations': 999,
    })
permanova_df = pd.DataFrame(permanova_results)
permanova_df.to_csv('results/stats/permanova_results.tsv', sep='\t', index=False)

# =============================================================================
# 7. Differential Abundance Analysis
# =============================================================================
diff_results = []
for genus in GENERA:
    for comparison in ['IBD_vs_Healthy', 'T2D_vs_Healthy']:
        control = taxonomy_df.loc[[s for s in SAMPLE_IDS if SAMPLES[s] == 'Healthy'], genus].values
        if 'IBD' in comparison:
            case = taxonomy_df.loc[[s for s in SAMPLE_IDS if SAMPLES[s] == 'IBD'], genus].values
        else:
            case = taxonomy_df.loc[[s for s in SAMPLE_IDS if SAMPLES[s] == 'T2D'], genus].values
        
        t_stat, p_val = stats.ttest_ind(case, control)
        log2fc = np.log2((case.mean() + 0.01) / (control.mean() + 0.01))
        
        diff_results.append({
            'Feature': genus,
            'Comparison': comparison,
            'Log2FC': round(log2fc, 4),
            'T_statistic': round(t_stat, 4),
            'P_value': round(p_val, 6),
            'Mean_Case': round(case.mean(), 4),
            'Mean_Control': round(control.mean(), 4),
        })

diff_df = pd.DataFrame(diff_results)
# BH correction
for comp in diff_df['Comparison'].unique():
    mask = diff_df['Comparison'] == comp
    p_vals = diff_df.loc[mask, 'P_value'].values
    n = len(p_vals)
    ranks = stats.rankdata(p_vals)
    q_vals = p_vals * n / ranks
    q_vals = np.minimum.accumulate(q_vals[np.argsort(ranks)[::-1]])[::-1]
    q_vals = np.clip(q_vals, 0, 1)
    sorted_idx = np.argsort(ranks)
    diff_df.loc[mask, 'Q_value'] = np.round(q_vals[np.argsort(sorted_idx)], 6)

diff_df.to_csv('results/stats/differential_abundance.tsv', sep='\t', index=False)

# =============================================================================
# 8. Simulate Classifier Comparison Data
# =============================================================================
classifier_comparison = []
for genus in GENERA:
    for sample in SAMPLE_IDS:
        mp4_val = taxonomy_df.loc[sample, genus]
        k2_val = kraken2_df.loc[sample, genus]
        classifier_comparison.append({
            'Genus': genus,
            'Sample': sample,
            'MetaPhlAn4': round(mp4_val, 4),
            'Kraken2': round(k2_val, 4),
            'Difference': round(k2_val - mp4_val, 4),
        })
comp_df = pd.DataFrame(classifier_comparison)
comp_df.to_csv('results/taxonomy/classifier_comparison.tsv', sep='\t', index=False)

# =============================================================================
# FIGURES
# =============================================================================
plt.rcParams.update({
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 150,
    'savefig.bbox': 'tight',
})

# --- Figure 1: Taxonomy Barplot ---
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Top 10 genera stacked barplot
top_genera = taxonomy_df.mean().sort_values(ascending=False).head(10).index.tolist()
plot_data = taxonomy_df[top_genera]
other = 100 - plot_data.sum(axis=1)

cmap = plt.cm.Set3(np.linspace(0, 1, len(top_genera) + 1))
bottom = np.zeros(len(SAMPLE_IDS))
for j, genus in enumerate(top_genera):
    axes[0].bar(range(len(SAMPLE_IDS)), plot_data[genus].values, bottom=bottom,
                color=cmap[j], label=genus, edgecolor='white', linewidth=0.5)
    bottom += plot_data[genus].values
axes[0].bar(range(len(SAMPLE_IDS)), other.values, bottom=bottom,
            color='#cccccc', label='Other', edgecolor='white', linewidth=0.5)

axes[0].set_xticks(range(len(SAMPLE_IDS)))
axes[0].set_xticklabels([s.replace('sample_', '') for s in SAMPLE_IDS], rotation=45, ha='right')
axes[0].set_ylabel('Relative Abundance (%)')
axes[0].set_title('(A) Taxonomic Composition (Genus Level)')
axes[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)

# Group-averaged composition
group_avg = taxonomy_df.copy()
group_avg['Group'] = GROUPS
group_avg = group_avg.groupby('Group')[top_genera].mean()
bottom = np.zeros(len(group_avg))
for j, genus in enumerate(top_genera):
    axes[1].bar(range(len(group_avg)), group_avg[genus].values, bottom=bottom,
                color=cmap[j], label=genus, edgecolor='white', linewidth=0.5)
    bottom += group_avg[genus].values

axes[1].set_xticks(range(len(group_avg)))
axes[1].set_xticklabels(group_avg.index, rotation=0)
axes[1].set_ylabel('Mean Relative Abundance (%)')
axes[1].set_title('(B) Group-Averaged Composition')

plt.tight_layout()
plt.savefig('figures/taxonomy_barplot.png')
plt.close()

# --- Figure 2: Alpha Diversity Boxplot ---
fig, axes = plt.subplots(1, 3, figsize=(14, 5))
metrics = ['Shannon', 'Simpson', 'Observed_Features']
titles = ['Shannon Diversity', 'Simpson Diversity', 'Observed Features']

for ax, metric, title in zip(axes, metrics, titles):
    positions = []
    data_groups = []
    group_names = ['Healthy', 'IBD', 'T2D']
    for k, g in enumerate(group_names):
        vals = alpha_df[alpha_df['Group'] == g][metric].values
        data_groups.append(vals)
        positions.append(k)
    
    bp = ax.boxplot(data_groups, positions=positions, widths=0.6, patch_artist=True,
                    showmeans=True, meanprops={'marker': 'D', 'markerfacecolor': 'black', 'markersize': 5})
    for patch, g in zip(bp['boxes'], group_names):
        patch.set_facecolor(COLORS[g])
        patch.set_alpha(0.7)
    
    # Add individual points
    for k, (g, vals) in enumerate(zip(group_names, data_groups)):
        jitter = np.random.uniform(-0.1, 0.1, len(vals))
        ax.scatter([k + j for j in jitter], vals, color=COLORS[g], edgecolor='black',
                   s=50, zorder=5, alpha=0.8)
    
    ax.set_xticklabels(group_names)
    ax.set_title(title)
    ax.set_ylabel(metric)

    # Add Kruskal-Wallis p-value
    try:
        h_stat, kw_p = stats.kruskal(*data_groups)
        ax.text(0.02, 0.98, f'KW p = {kw_p:.4f}', transform=ax.transAxes,
                va='top', fontsize=9, fontstyle='italic')
    except ValueError:
        ax.text(0.02, 0.98, 'KW p = N/A', transform=ax.transAxes,
                va='top', fontsize=9, fontstyle='italic')

plt.suptitle('Alpha Diversity by Disease Group', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('figures/alpha_diversity_boxplot.png')
plt.close()

# --- Figure 3: Beta Diversity PCoA ---
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
distance_matrices = [('Bray-Curtis', bc_dm), ('Jaccard', jaccard_dm), ('Aitchison', aitchison_dm)]

for ax, (metric_name, dm) in zip(axes, distance_matrices):
    coords, var_exp = pcoa(dm)
    for g in ['Healthy', 'IBD', 'T2D']:
        idx = [i for i, x in enumerate(GROUPS) if x == g]
        ax.scatter(coords[idx, 0], coords[idx, 1], c=COLORS[g], label=g,
                   s=120, edgecolor='black', linewidth=0.8, zorder=5)
    
    # Draw convex hull-like ellipses
    for g in ['Healthy', 'IBD', 'T2D']:
        idx = [i for i, x in enumerate(GROUPS) if x == g]
        if len(idx) >= 3:
            center = coords[idx].mean(axis=0)
            ax.annotate('', xy=center, xytext=center)
    
    ax.set_xlabel(f'PC1 ({var_exp[0]:.1f}%)')
    ax.set_ylabel(f'PC2 ({var_exp[1]:.1f}%)')
    ax.set_title(f'{metric_name} PCoA')
    ax.legend(loc='best')
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.3)
    
    # Add PERMANOVA result
    perm_row = permanova_df[permanova_df['Metric'] == metric_name].iloc[0]
    ax.text(0.02, 0.02, f'PERMANOVA R²={perm_row["R2"]:.3f}, p={perm_row["p_value"]:.3f}',
            transform=ax.transAxes, fontsize=9, fontstyle='italic')

plt.suptitle('Beta Diversity Analysis (PCoA)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('figures/beta_diversity_pcoa.png')
plt.close()

# --- Figure 4: Functional Heatmap ---
fig, ax = plt.subplots(figsize=(14, 10))

# Standardize pathway data
pathway_z = pathway_df.copy()
pathway_z = (pathway_z - pathway_z.mean()) / pathway_z.std()

# Cluster rows and columns
row_linkage = linkage(pdist(pathway_z.values, 'euclidean'), method='ward')
col_linkage = linkage(pdist(pathway_z.T.values, 'euclidean'), method='ward')

# Get dendrogram order
row_order = dendrogram(row_linkage, no_plot=True)['leaves']
col_order = dendrogram(col_linkage, no_plot=True)['leaves']

ordered_data = pathway_z.iloc[row_order, col_order]
short_names = [p.split(': ')[1][:40] if ': ' in p else p[:40] for p in ordered_data.columns]

im = ax.imshow(ordered_data.values, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
ax.set_xticks(range(len(short_names)))
ax.set_xticklabels(short_names, rotation=90, fontsize=8)
ax.set_yticks(range(len(ordered_data)))
ax.set_yticklabels([s.replace('sample_', '') for s in ordered_data.index], fontsize=9)

# Add group color bar
group_colors = [COLORS[SAMPLES[s]] for s in ordered_data.index]
for i, c in enumerate(group_colors):
    ax.add_patch(plt.Rectangle((-1.5, i-0.5), 0.8, 1, facecolor=c, edgecolor='none'))

plt.colorbar(im, ax=ax, label='Z-score', shrink=0.6)
ax.set_title('Metabolic Pathway Abundance (Z-score normalized)', fontsize=13)
plt.tight_layout()
plt.savefig('figures/functional_heatmap.png')
plt.close()

# --- Figure 5: MAG Quality Scatter ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Completeness vs Contamination
for g in ['Healthy', 'IBD', 'T2D']:
    mask = mag_df['Group'] == g
    axes[0].scatter(mag_df.loc[mask, 'Completeness'], mag_df.loc[mask, 'Contamination'],
                    c=COLORS[g], label=g, s=60, alpha=0.7, edgecolor='black', linewidth=0.3)

# Quality thresholds
axes[0].axvline(x=50, color='orange', linestyle='--', alpha=0.6, label='Medium quality (≥50%)')
axes[0].axvline(x=90, color='green', linestyle='--', alpha=0.6, label='High quality (≥90%)')
axes[0].axhline(y=5, color='green', linestyle=':', alpha=0.6)
axes[0].axhline(y=10, color='orange', linestyle=':', alpha=0.6)
axes[0].set_xlabel('Completeness (%)')
axes[0].set_ylabel('Contamination (%)')
axes[0].set_title('(A) MAG Quality Assessment')
axes[0].legend(loc='upper left', fontsize=8)

# Quality category distribution
quality_counts = mag_df.groupby(['Group', 'Quality_Category']).size().unstack(fill_value=0)
quality_counts = quality_counts.reindex(columns=['High', 'Medium', 'Low'])
quality_counts.plot(kind='bar', ax=axes[1], color=['#2ecc71', '#f39c12', '#e74c3c'],
                    edgecolor='black', linewidth=0.5)
axes[1].set_title('(B) MAG Quality Distribution by Group')
axes[1].set_ylabel('Number of MAGs')
axes[1].set_xlabel('')
axes[1].legend(title='Quality')
axes[1].tick_params(axis='x', rotation=0)

plt.tight_layout()
plt.savefig('figures/mag_quality_scatter.png')
plt.close()

# --- Figure 6: Differential Abundance Volcano ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, comparison in zip(axes, ['IBD_vs_Healthy', 'T2D_vs_Healthy']):
    sub = diff_df[diff_df['Comparison'] == comparison].copy()
    sub['-log10(Q)'] = -np.log10(sub['Q_value'] + 1e-10)
    
    sig = (sub['Q_value'] < 0.05) & (abs(sub['Log2FC']) > 0.5)
    up = sig & (sub['Log2FC'] > 0)
    down = sig & (sub['Log2FC'] < 0)
    ns = ~sig
    
    ax.scatter(sub.loc[ns, 'Log2FC'], sub.loc[ns, '-log10(Q)'],
               c='gray', alpha=0.5, s=40, label='NS')
    ax.scatter(sub.loc[up, 'Log2FC'], sub.loc[up, '-log10(Q)'],
               c='#e74c3c', alpha=0.8, s=60, label='Up')
    ax.scatter(sub.loc[down, 'Log2FC'], sub.loc[down, '-log10(Q)'],
               c='#3498db', alpha=0.8, s=60, label='Down')
    
    # Label significant features
    for _, row in sub[sig].iterrows():
        ax.annotate(row['Feature'], (row['Log2FC'], row['-log10(Q)']),
                    fontsize=7, ha='center', va='bottom',
                    arrowprops=dict(arrowstyle='-', color='gray', alpha=0.3))
    
    ax.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=-0.5, color='gray', linestyle='--', alpha=0.3)
    ax.axvline(x=0.5, color='gray', linestyle='--', alpha=0.3)
    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-log10(Q-value)')
    ax.set_title(comparison.replace('_', ' '))
    ax.legend(loc='upper right', fontsize=8)

plt.suptitle('Differential Abundance Analysis', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('figures/differential_abundance_volcano.png')
plt.close()

# --- Figure 7: Classifier Comparison ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Correlation plot
all_mp4 = comp_df.groupby('Genus')['MetaPhlAn4'].mean()
all_k2 = comp_df.groupby('Genus')['Kraken2'].mean()
axes[0].scatter(all_mp4, all_k2, s=80, c='#3498db', edgecolor='black', alpha=0.7)
for genus in GENERA:
    axes[0].annotate(genus, (all_mp4[genus], all_k2[genus]), fontsize=7,
                     ha='left', va='bottom')

max_val = max(all_mp4.max(), all_k2.max()) * 1.1
axes[0].plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='y=x')
r, p = stats.pearsonr(all_mp4, all_k2)
axes[0].set_xlabel('MetaPhlAn4 Abundance (%)')
axes[0].set_ylabel('Kraken2 Abundance (%)')
axes[0].set_title(f'(A) Classifier Agreement (r={r:.3f}, p={p:.2e})')
axes[0].legend()

# Bland-Altman plot
mean_vals = (all_mp4 + all_k2) / 2
diff_vals = all_k2 - all_mp4
axes[1].scatter(mean_vals, diff_vals, s=80, c='#e74c3c', edgecolor='black', alpha=0.7)
mean_diff = diff_vals.mean()
std_diff = diff_vals.std()
axes[1].axhline(y=mean_diff, color='blue', linestyle='-', alpha=0.5, label=f'Mean diff: {mean_diff:.2f}')
axes[1].axhline(y=mean_diff + 1.96*std_diff, color='red', linestyle='--', alpha=0.5, label='±1.96 SD')
axes[1].axhline(y=mean_diff - 1.96*std_diff, color='red', linestyle='--', alpha=0.5)
axes[1].set_xlabel('Mean Abundance (%)')
axes[1].set_ylabel('Difference (Kraken2 - MetaPhlAn4)')
axes[1].set_title('(B) Bland-Altman Plot')
axes[1].legend(fontsize=8)

plt.suptitle('Taxonomic Classifier Comparison: MetaPhlAn4 vs Kraken2', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('figures/classifier_comparison.png')
plt.close()

# --- Figure 8: Binning Comparison ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

binner_stats = mag_df.groupby('Binner_Source').agg({
    'Completeness': ['mean', 'std', 'count'],
    'Contamination': ['mean', 'std'],
}).round(2)

binners = ['MetaBAT2', 'CONCOCT', 'MaxBin2', 'DAS_Tool_consensus']
binner_colors = ['#3498db', '#e74c3c', '#2ecc71', '#9b59b6']

comp_means = [mag_df[mag_df['Binner_Source'] == b]['Completeness'].mean() for b in binners]
comp_stds = [mag_df[mag_df['Binner_Source'] == b]['Completeness'].std() for b in binners]
cont_means = [mag_df[mag_df['Binner_Source'] == b]['Contamination'].mean() for b in binners]
cont_stds = [mag_df[mag_df['Binner_Source'] == b]['Contamination'].std() for b in binners]

x = np.arange(len(binners))
axes[0].bar(x, comp_means, yerr=comp_stds, color=binner_colors, edgecolor='black',
            capsize=5, alpha=0.8)
axes[0].set_xticks(x)
axes[0].set_xticklabels(binners, rotation=30, ha='right')
axes[0].set_ylabel('Completeness (%)')
axes[0].set_title('(A) Average Completeness by Binner')
axes[0].set_ylim(0, 105)

axes[1].bar(x, cont_means, yerr=cont_stds, color=binner_colors, edgecolor='black',
            capsize=5, alpha=0.8)
axes[1].set_xticks(x)
axes[1].set_xticklabels(binners, rotation=30, ha='right')
axes[1].set_ylabel('Contamination (%)')
axes[1].set_title('(B) Average Contamination by Binner')

plt.suptitle('Genome Binning Tool Comparison', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('figures/binning_comparison.png')
plt.close()

# --- Figure 9: Pipeline Overview / QC Summary ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# QC read retention waterfall
categories = ['Raw Reads', 'After Adapter\nRemoval', 'After\nDeduplication', 'After Host\nRemoval']
avg_values = [
    qc_df['Total_Reads'].mean() / 1e6,
    qc_df['After_Adapter_Removal'].mean() / 1e6,
    qc_df['After_Deduplication'].mean() / 1e6,
    qc_df['After_Host_Removal'].mean() / 1e6,
]
colors_qc = ['#3498db', '#2ecc71', '#f39c12', '#e74c3c']
axes[0].bar(range(len(categories)), avg_values, color=colors_qc, edgecolor='black', alpha=0.8)
axes[0].set_xticks(range(len(categories)))
axes[0].set_xticklabels(categories, fontsize=9)
axes[0].set_ylabel('Reads (millions)')
axes[0].set_title('(A) QC Read Retention (Average)')
for i, v in enumerate(avg_values):
    axes[0].text(i, v + 0.3, f'{v:.1f}M', ha='center', fontsize=9)

# Q30 distribution
axes[1].bar(range(len(qc_df)), qc_df['Q30_Pct'].values,
            color=[COLORS[g] for g in qc_df['Group']], edgecolor='black', alpha=0.8)
axes[1].set_xticks(range(len(qc_df)))
axes[1].set_xticklabels([s.replace('sample_', '') for s in qc_df['Sample']], rotation=45, ha='right')
axes[1].set_ylabel('Q30 (%)')
axes[1].set_title('(B) Sequencing Quality (Q30) per Sample')
axes[1].axhline(y=90, color='red', linestyle='--', alpha=0.5, label='Q30=90% threshold')
axes[1].legend()

plt.tight_layout()
plt.savefig('figures/qc_summary.png')
plt.close()

print("All figures and data generated successfully!")
print(f"  Figures: {len(os.listdir('figures'))} files in figures/")
print(f"  MAGs: {len(mag_df)} total, {len(mag_df[mag_df['Quality_Category']=='High'])} high quality")
print(f"  Significant taxa (IBD): {len(diff_df[(diff_df['Comparison']=='IBD_vs_Healthy') & (diff_df['Q_value']<0.05)])}")
print(f"  Significant taxa (T2D): {len(diff_df[(diff_df['Comparison']=='T2D_vs_Healthy') & (diff_df['Q_value']<0.05)])}")
