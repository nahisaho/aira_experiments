"""
Integrated metabolomics-microbiome analysis pipeline.
Implements: correlation network, sPLS integration, causal inference,
pathway enrichment, biomarker scoring, and IBD case study.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy import stats
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import ElasticNet, LogisticRegression
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import roc_auc_score, roc_curve, classification_report, confusion_matrix
from sklearn.feature_selection import mutual_info_classif
import json
import warnings
warnings.filterwarnings('ignore')

# --------------------------------------------------
# Load data
# --------------------------------------------------
taxa_df = pd.read_csv('data/taxa_clr.csv', index_col=0)
met_df = pd.read_csv('data/metabolites.csv', index_col=0)
meta_df = pd.read_csv('data/metadata.csv')
met_class_df = pd.read_csv('data/metabolite_classes.csv')
with open('data/pathway_metabolite_map.json') as f:
    pathway_met_map = json.load(f)
with open('data/pathway_taxa_map.json') as f:
    pathway_taxa_map = json.load(f)

groups = meta_df['Group'].values
le = LabelEncoder()
y = le.fit_transform(groups)  # 0=Control, 1=IBD

print("=" * 60)
print("INTEGRATED METABOLOMICS-MICROBIOME ANALYSIS PIPELINE")
print("=" * 60)

# ==================================================
# Module 1: Peak Annotation Simulation (Untargeted Metabolomics)
# ==================================================
print("\n--- Module 1: Peak Annotation Automation ---")

# Simulate annotation confidence levels
np.random.seed(42)
n_met = met_df.shape[1]
annotation_levels = np.random.choice(
    ['Level1_confirmed', 'Level2_putative', 'Level3_class', 'Level4_unknown'],
    size=n_met, p=[0.25, 0.35, 0.25, 0.15]
)
mz_values = np.random.uniform(100, 900, n_met)
rt_values = np.random.uniform(0.5, 25, n_met)
annotation_scores = np.where(
    annotation_levels == 'Level1_confirmed', np.random.uniform(0.85, 1.0, n_met),
    np.where(annotation_levels == 'Level2_putative', np.random.uniform(0.6, 0.85, n_met),
    np.where(annotation_levels == 'Level3_class', np.random.uniform(0.3, 0.6, n_met),
    np.random.uniform(0.05, 0.3, n_met)))
)

annot_df = pd.DataFrame({
    'Metabolite': met_df.columns,
    'mz': mz_values.round(4),
    'RT_min': rt_values.round(2),
    'Annotation_Level': annotation_levels,
    'Confidence_Score': annotation_scores.round(3),
    'Class': met_class_df['Class'].values
})
annot_df.to_csv('data/peak_annotations.csv', index=False)

level_counts = annot_df['Annotation_Level'].value_counts()
print(f"Annotation summary: {dict(level_counts)}")
print(f"Mean confidence: {annotation_scores.mean():.3f}")

# Figure 1: Annotation summary
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
level_counts.plot(kind='bar', ax=axes[0], color=['#2ecc71', '#3498db', '#e67e22', '#e74c3c'])
axes[0].set_title('Annotation Level Distribution')
axes[0].set_ylabel('Count')
axes[0].tick_params(axis='x', rotation=45)

axes[1].hist(annotation_scores, bins=20, color='#3498db', edgecolor='black', alpha=0.7)
axes[1].set_xlabel('Confidence Score')
axes[1].set_ylabel('Count')
axes[1].set_title('Annotation Confidence Distribution')

class_counts = annot_df['Class'].value_counts()
class_counts.plot(kind='barh', ax=axes[2], color='#9b59b6')
axes[2].set_xlabel('Count')
axes[2].set_title('Metabolite Class Distribution')

plt.tight_layout()
plt.savefig('figures/fig1_annotation_summary.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/fig1_annotation_summary.png")

# ==================================================
# Module 2: Correlation Network (Taxa-Metabolite)
# ==================================================
print("\n--- Module 2: Correlation Network ---")

# Compute Spearman correlations between top taxa and metabolites
top_taxa = taxa_df.columns[:30]
top_mets = met_df.columns[:50]

corr_matrix = np.zeros((len(top_taxa), len(top_mets)))
pval_matrix = np.zeros_like(corr_matrix)

for i, taxon in enumerate(top_taxa):
    for j, metab in enumerate(top_mets):
        r, p = spearmanr(taxa_df[taxon], met_df[metab])
        corr_matrix[i, j] = r
        pval_matrix[i, j] = p

# BH correction
from statsmodels.stats.multitest import multipletests
pvals_flat = pval_matrix.flatten()
_, pvals_adj, _, _ = multipletests(pvals_flat, method='fdr_bh')
pval_adj_matrix = pvals_adj.reshape(pval_matrix.shape)

# Build network: edges where |r| > 0.3 and FDR < 0.05
G = nx.Graph()
sig_edges = []
for i, taxon in enumerate(top_taxa):
    for j, metab in enumerate(top_mets):
        if abs(corr_matrix[i, j]) > 0.3 and pval_adj_matrix[i, j] < 0.05:
            G.add_node(taxon, node_type='taxon')
            G.add_node(metab, node_type='metabolite')
            G.add_edge(taxon, metab, weight=corr_matrix[i, j],
                       pval=pval_adj_matrix[i, j])
            sig_edges.append({
                'Taxon': taxon, 'Metabolite': metab,
                'Spearman_r': round(corr_matrix[i, j], 4),
                'FDR_q': round(pval_adj_matrix[i, j], 6)
            })

sig_edges_df = pd.DataFrame(sig_edges)
sig_edges_df.to_csv('data/significant_correlations.csv', index=False)
print(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# Figure 2: Correlation heatmap (top associations)
fig, axes = plt.subplots(1, 2, figsize=(18, 8))

# Heatmap of top correlations
top_corr = pd.DataFrame(corr_matrix[:15, :20], index=top_taxa[:15], columns=top_mets[:20])
sns.heatmap(top_corr, cmap='RdBu_r', center=0, ax=axes[0], vmin=-0.8, vmax=0.8,
            xticklabels=True, yticklabels=True)
axes[0].set_title('Taxa-Metabolite Spearman Correlations')
axes[0].tick_params(axis='x', rotation=90, labelsize=8)
axes[0].tick_params(axis='y', labelsize=8)

# Network visualization
if G.number_of_nodes() > 0:
    pos = nx.spring_layout(G, k=2, seed=42)
    node_colors = ['#e74c3c' if G.nodes[n].get('node_type') == 'taxon' else '#3498db'
                   for n in G.nodes()]
    edge_weights = [abs(G[u][v]['weight']) * 3 for u, v in G.edges()]
    edge_colors = ['red' if G[u][v]['weight'] < 0 else 'blue' for u, v in G.edges()]

    nx.draw(G, pos, ax=axes[1], node_color=node_colors, node_size=300,
            width=edge_weights, edge_color=edge_colors, alpha=0.7,
            with_labels=True, font_size=6)
    axes[1].set_title('Taxa-Metabolite Correlation Network')

plt.tight_layout()
plt.savefig('figures/fig2_correlation_network.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/fig2_correlation_network.png")

# ==================================================
# Module 3: sPLS Integration (mixOmics-style)
# ==================================================
print("\n--- Module 3: sPLS Integration ---")

scaler = StandardScaler()
X_taxa = scaler.fit_transform(taxa_df.values)
X_met = scaler.fit_transform(met_df.values)

# PLS regression (taxa -> metabolites)
n_components = 5
pls = PLSRegression(n_components=n_components, max_iter=500)
pls.fit(X_taxa, X_met)

X_scores = pls.x_scores_
Y_scores = pls.y_scores_

# Explained variance per component
x_var = []
y_var = []
for i in range(n_components):
    x_var.append(np.var(X_scores[:, i]) / np.var(X_taxa).sum() * 100)
    y_var.append(np.var(Y_scores[:, i]) / np.var(X_met).sum() * 100)

print(f"PLS components explained variance (X): {[f'{v:.1f}%' for v in x_var]}")
print(f"PLS components explained variance (Y): {[f'{v:.1f}%' for v in y_var]}")

# MelonnPan-style prediction: ElasticNet for metabolite prediction
print("\n--- MelonnPan-style Metabolite Prediction ---")
from sklearn.model_selection import KFold

melonnpan_results = []
kf = KFold(n_splits=5, shuffle=True, random_state=42)

# Predict top metabolites from taxa
target_metabolites = ['Butyrate', 'Propionate', 'Indoxyl_sulfate', 'Tryptophan', 'Kynurenine',
                      'Deoxycholic_acid', 'Acetate', 'Hippuric_acid', 'p_Cresol_sulfate', 'Serotonin']

for metab in target_metabolites:
    if metab in met_df.columns:
        y_met = met_df[metab].values
        scores = []
        for train_idx, test_idx in kf.split(X_taxa):
            en = ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=2000)
            en.fit(X_taxa[train_idx], y_met[train_idx])
            pred = en.predict(X_taxa[test_idx])
            r, _ = spearmanr(y_met[test_idx], pred)
            scores.append(r)
        mean_r = np.mean(scores)
        melonnpan_results.append({'Metabolite': metab, 'CV_Spearman_r': round(mean_r, 4)})

melonnpan_df = pd.DataFrame(melonnpan_results)
melonnpan_df.to_csv('data/melonnpan_prediction.csv', index=False)
print(f"MelonnPan prediction results:\n{melonnpan_df.to_string(index=False)}")

# Figure 3: PLS sample scores + MelonnPan
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

colors = ['#e74c3c' if g == 'IBD' else '#3498db' for g in groups]
axes[0].scatter(X_scores[:, 0], X_scores[:, 1], c=colors, alpha=0.6, s=40)
axes[0].set_xlabel(f'PLS Component 1 ({x_var[0]:.1f}%)')
axes[0].set_ylabel(f'PLS Component 2 ({x_var[1]:.1f}%)')
axes[0].set_title('sPLS Sample Scores (Taxa)')
from matplotlib.patches import Patch
axes[0].legend(handles=[Patch(color='#e74c3c', label='IBD'), Patch(color='#3498db', label='Control')])

axes[1].scatter(Y_scores[:, 0], Y_scores[:, 1], c=colors, alpha=0.6, s=40)
axes[1].set_xlabel(f'PLS Component 1 ({y_var[0]:.1f}%)')
axes[1].set_ylabel(f'PLS Component 2 ({y_var[1]:.1f}%)')
axes[1].set_title('sPLS Sample Scores (Metabolites)')
axes[1].legend(handles=[Patch(color='#e74c3c', label='IBD'), Patch(color='#3498db', label='Control')])

# MelonnPan bar plot
melonnpan_df_sorted = melonnpan_df.sort_values('CV_Spearman_r', ascending=True)
bars = axes[2].barh(melonnpan_df_sorted['Metabolite'], melonnpan_df_sorted['CV_Spearman_r'],
                     color=['#2ecc71' if v > 0.3 else '#e67e22' if v > 0 else '#e74c3c'
                            for v in melonnpan_df_sorted['CV_Spearman_r']])
axes[2].set_xlabel('Cross-validated Spearman r')
axes[2].set_title('MelonnPan Metabolite Prediction')
axes[2].axvline(x=0.3, color='gray', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('figures/fig3_spls_melonnpan.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/fig3_spls_melonnpan.png")

# ==================================================
# Module 4: Causal Inference (Granger Causality)
# ==================================================
print("\n--- Module 4: Causal Inference (Granger-style) ---")

# Simulate temporal data for Granger causality
np.random.seed(42)
n_timepoints = 50
# Faecalibacterium -> Butyrate (causal)
faecal_ts = np.cumsum(np.random.randn(n_timepoints) * 0.5) + 5
butyrate_ts = np.zeros(n_timepoints)
butyrate_ts[0] = 3
for t in range(1, n_timepoints):
    butyrate_ts[t] = 0.3 * butyrate_ts[t-1] + 0.5 * faecal_ts[t-1] + np.random.randn() * 0.3

# Escherichia -> Indoxyl_sulfate
ecoli_ts = np.cumsum(np.random.randn(n_timepoints) * 0.4) + 3
indoxyl_ts = np.zeros(n_timepoints)
indoxyl_ts[0] = 1
for t in range(1, n_timepoints):
    indoxyl_ts[t] = 0.2 * indoxyl_ts[t-1] + 0.4 * ecoli_ts[t-1] + np.random.randn() * 0.3

# Granger causality test
from statsmodels.tsa.stattools import grangercausalitytests

granger_results = []
pairs = [
    ('Faecalibacterium', faecal_ts, 'Butyrate', butyrate_ts),
    ('Escherichia', ecoli_ts, 'Indoxyl_sulfate', indoxyl_ts),
]
for taxon_name, taxon_ts, met_name, met_ts in pairs:
    data = np.column_stack([met_ts, taxon_ts])
    try:
        result = grangercausalitytests(data, maxlag=3, verbose=False)
        for lag in [1, 2, 3]:
            f_stat = result[lag][0]['ssr_ftest'][0]
            p_val = result[lag][0]['ssr_ftest'][1]
            granger_results.append({
                'Cause': taxon_name, 'Effect': met_name,
                'Lag': lag, 'F_stat': round(f_stat, 3), 'p_value': round(p_val, 6)
            })
    except Exception as e:
        print(f"Granger test failed for {taxon_name}->{met_name}: {e}")

granger_df = pd.DataFrame(granger_results)
granger_df.to_csv('data/granger_causality.csv', index=False)
print(f"Granger causality results:\n{granger_df.to_string(index=False)}")

# Mendelian Randomization simulation
print("\n--- Mendelian Randomization (Simulated) ---")
np.random.seed(42)
n_snps = 20
snp_effects_on_taxa = np.random.randn(n_snps) * 0.1
snp_effects_on_met = snp_effects_on_taxa * 0.6 + np.random.randn(n_snps) * 0.05  # causal path

# IVW estimator
beta_ivw = np.sum(snp_effects_on_taxa * snp_effects_on_met) / np.sum(snp_effects_on_taxa**2)
se_ivw = np.sqrt(np.sum((snp_effects_on_met - beta_ivw * snp_effects_on_taxa)**2) /
                 ((n_snps - 1) * np.sum(snp_effects_on_taxa**2)))
z_ivw = beta_ivw / se_ivw
p_ivw = 2 * (1 - stats.norm.cdf(abs(z_ivw)))

mr_results = {
    'Exposure': 'Faecalibacterium_abundance',
    'Outcome': 'Butyrate_level',
    'N_SNPs': n_snps,
    'Beta_IVW': round(beta_ivw, 4),
    'SE_IVW': round(se_ivw, 4),
    'Z_stat': round(z_ivw, 4),
    'P_value': round(p_ivw, 6),
    'Method': 'IVW'
}
print(f"MR IVW result: beta={beta_ivw:.4f}, SE={se_ivw:.4f}, p={p_ivw:.6f}")

# Figure 4: Causal inference
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Time series
axes[0, 0].plot(faecal_ts, label='Faecalibacterium', color='#e74c3c')
ax2 = axes[0, 0].twinx()
ax2.plot(butyrate_ts, label='Butyrate', color='#3498db')
axes[0, 0].set_xlabel('Time')
axes[0, 0].set_ylabel('Faecalibacterium', color='#e74c3c')
ax2.set_ylabel('Butyrate', color='#3498db')
axes[0, 0].set_title('Temporal Dynamics: Faecalibacterium → Butyrate')

axes[0, 1].plot(ecoli_ts, label='Escherichia', color='#e74c3c')
ax3 = axes[0, 1].twinx()
ax3.plot(indoxyl_ts, label='Indoxyl sulfate', color='#3498db')
axes[0, 1].set_xlabel('Time')
axes[0, 1].set_ylabel('Escherichia', color='#e74c3c')
ax3.set_ylabel('Indoxyl sulfate', color='#3498db')
axes[0, 1].set_title('Temporal Dynamics: Escherichia → Indoxyl sulfate')

# MR scatter
axes[1, 0].scatter(snp_effects_on_taxa, snp_effects_on_met, c='#9b59b6', alpha=0.7, s=60)
x_range = np.linspace(min(snp_effects_on_taxa), max(snp_effects_on_taxa), 100)
axes[1, 0].plot(x_range, beta_ivw * x_range, 'r-', label=f'IVW: β={beta_ivw:.3f}')
axes[1, 0].set_xlabel('SNP effect on Faecalibacterium')
axes[1, 0].set_ylabel('SNP effect on Butyrate')
axes[1, 0].set_title(f'MR Scatter Plot (p={p_ivw:.4f})')
axes[1, 0].legend()

# Granger F-stats
if len(granger_df) > 0:
    granger_pivot = granger_df.pivot_table(index=['Cause', 'Effect'], columns='Lag', values='F_stat')
    granger_pivot.plot(kind='bar', ax=axes[1, 1], colormap='viridis')
    axes[1, 1].set_title('Granger Causality F-statistics')
    axes[1, 1].set_ylabel('F-statistic')
    axes[1, 1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig('figures/fig4_causal_inference.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/fig4_causal_inference.png")

# ==================================================
# Module 5: Pathway Enrichment Analysis
# ==================================================
print("\n--- Module 5: Pathway Enrichment Analysis ---")

# Differential metabolites (IBD vs Control)
diff_results = []
for metab in met_df.columns:
    ibd_vals = met_df.loc[groups == 'IBD', metab]
    ctrl_vals = met_df.loc[groups == 'Control', metab]
    stat, pval = mannwhitneyu(ibd_vals, ctrl_vals, alternative='two-sided')
    fc = ibd_vals.mean() - ctrl_vals.mean()
    diff_results.append({'Metabolite': metab, 'MeanDiff': fc, 'U_stat': stat, 'p_value': pval})

diff_df = pd.DataFrame(diff_results)
_, diff_df['FDR_q'], _, _ = multipletests(diff_df['p_value'], method='fdr_bh')
diff_df.to_csv('data/differential_metabolites.csv', index=False)

sig_mets = set(diff_df.loc[diff_df['FDR_q'] < 0.05, 'Metabolite'].values)
print(f"Significant differential metabolites (FDR<0.05): {len(sig_mets)}")

# Pathway enrichment (Fisher's exact test)
all_mets = set(met_df.columns)
enrichment_results = []
for pathway, pathway_mets in pathway_met_map.items():
    pathway_set = set(pathway_mets) & all_mets
    if len(pathway_set) == 0:
        continue
    overlap = pathway_set & sig_mets
    a = len(overlap)
    b = len(pathway_set) - a
    c = len(sig_mets) - a
    d = len(all_mets) - a - b - c
    odds, pval = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
    enrichment_results.append({
        'Pathway': pathway,
        'Pathway_size': len(pathway_set),
        'Overlap': a,
        'Odds_ratio': round(odds, 3),
        'p_value': round(pval, 6),
        'Enrichment_ratio': round(a / max(len(pathway_set), 1), 3)
    })

enrich_df = pd.DataFrame(enrichment_results)
if len(enrich_df) > 0:
    _, enrich_df['FDR_q'], _, _ = multipletests(enrich_df['p_value'], method='fdr_bh')
enrich_df = enrich_df.sort_values('p_value')
enrich_df.to_csv('data/pathway_enrichment.csv', index=False)
print(f"Pathway enrichment results:\n{enrich_df[['Pathway', 'Overlap', 'Enrichment_ratio', 'p_value', 'FDR_q']].to_string(index=False)}")

# Differential taxa
diff_taxa_results = []
for taxon in taxa_df.columns:
    ibd_vals = taxa_df.loc[groups == 'IBD', taxon]
    ctrl_vals = taxa_df.loc[groups == 'Control', taxon]
    stat, pval = mannwhitneyu(ibd_vals, ctrl_vals, alternative='two-sided')
    fc = ibd_vals.mean() - ctrl_vals.mean()
    diff_taxa_results.append({'Taxon': taxon, 'MeanDiff': fc, 'p_value': pval})

diff_taxa_df = pd.DataFrame(diff_taxa_results)
_, diff_taxa_df['FDR_q'], _, _ = multipletests(diff_taxa_df['p_value'], method='fdr_bh')
diff_taxa_df.to_csv('data/differential_taxa.csv', index=False)
sig_taxa = set(diff_taxa_df.loc[diff_taxa_df['FDR_q'] < 0.05, 'Taxon'].values)
print(f"Significant differential taxa: {len(sig_taxa)}")

# Figure 5: Pathway enrichment + volcano
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Volcano plot (metabolites)
diff_df['neg_log10_q'] = -np.log10(diff_df['FDR_q'].clip(1e-20))
colors_vol = []
for _, row in diff_df.iterrows():
    if row['FDR_q'] < 0.05 and abs(row['MeanDiff']) > 0.3:
        colors_vol.append('#e74c3c' if row['MeanDiff'] > 0 else '#3498db')
    else:
        colors_vol.append('gray')
axes[0].scatter(diff_df['MeanDiff'], diff_df['neg_log10_q'], c=colors_vol, alpha=0.6, s=30)
axes[0].axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
axes[0].axvline(0.3, color='gray', linestyle='--', alpha=0.3)
axes[0].axvline(-0.3, color='gray', linestyle='--', alpha=0.3)
axes[0].set_xlabel('Mean Difference (IBD - Control)')
axes[0].set_ylabel('-log10(FDR q-value)')
axes[0].set_title('Metabolite Volcano Plot')

# Pathway enrichment bar
if len(enrich_df) > 0:
    top_enrich = enrich_df.head(10).sort_values('p_value', ascending=False)
    colors_e = ['#e74c3c' if q < 0.05 else '#95a5a6' for q in top_enrich['FDR_q']]
    axes[1].barh(top_enrich['Pathway'], -np.log10(top_enrich['p_value'].clip(1e-10)),
                 color=colors_e)
    axes[1].set_xlabel('-log10(p-value)')
    axes[1].set_title('Pathway Enrichment Analysis')
    axes[1].axvline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)

# Integrated pathway-taxa-metabolite network
axes[2].set_title('Integrated Pathway Map')
axes[2].text(0.5, 0.5, f'Sig. Metabolites: {len(sig_mets)}\nSig. Taxa: {len(sig_taxa)}\nEnriched Pathways: {len(enrich_df[enrich_df["FDR_q"]<0.05]) if len(enrich_df) > 0 else 0}',
             transform=axes[2].transAxes, ha='center', va='center', fontsize=14,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('figures/fig5_pathway_enrichment.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/fig5_pathway_enrichment.png")

# ==================================================
# Module 6: Disease Biomarker Integrated Scoring
# ==================================================
print("\n--- Module 6: Biomarker Integrated Scoring ---")

# Combine taxa + metabolite features
X_combined = np.hstack([X_taxa, X_met])
feature_names = list(taxa_df.columns) + list(met_df.columns)

# Feature importance via Random Forest
rf = RandomForestClassifier(n_estimators=200, random_state=42, max_depth=10)
rf.fit(X_combined, y)
importances = rf.feature_importances_

# Top features
top_k = 30
top_idx = np.argsort(importances)[-top_k:][::-1]
top_features = [(feature_names[i], importances[i], 'Taxa' if i < len(taxa_df.columns) else 'Metabolite')
                for i in top_idx]
top_feat_df = pd.DataFrame(top_features, columns=['Feature', 'Importance', 'Type'])
top_feat_df.to_csv('data/top_biomarkers.csv', index=False)
print(f"Top 10 biomarkers:\n{top_feat_df.head(10).to_string(index=False)}")

# Cross-validated AUC comparison
models = {
    'RF_Taxa_only': (RandomForestClassifier(n_estimators=100, random_state=42), X_taxa),
    'RF_Metabolites_only': (RandomForestClassifier(n_estimators=100, random_state=42), X_met),
    'RF_Integrated': (RandomForestClassifier(n_estimators=100, random_state=42), X_combined),
    'GB_Integrated': (GradientBoostingClassifier(n_estimators=100, random_state=42), X_combined),
    'LR_Integrated': (LogisticRegression(max_iter=1000, random_state=42), X_combined),
}

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
auc_results = {}
roc_data = {}

for name, (model, X) in models.items():
    aucs = []
    all_y_true = []
    all_y_prob = []
    for train_idx, test_idx in cv.split(X, y):
        model_clone = type(model)(**model.get_params())
        model_clone.fit(X[train_idx], y[train_idx])
        y_prob = model_clone.predict_proba(X[test_idx])[:, 1]
        aucs.append(roc_auc_score(y[test_idx], y_prob))
        all_y_true.extend(y[test_idx])
        all_y_prob.extend(y_prob)
    auc_results[name] = {'mean_AUC': np.mean(aucs), 'std_AUC': np.std(aucs)}
    fpr, tpr, _ = roc_curve(all_y_true, all_y_prob)
    roc_data[name] = (fpr, tpr)
    print(f"  {name}: AUC = {np.mean(aucs):.4f} ± {np.std(aucs):.4f}")

auc_df = pd.DataFrame(auc_results).T
auc_df.to_csv('data/model_performance.csv')

# Integrated biomarker score
# Top features selected
top_feature_idx = top_idx[:15]
X_top = X_combined[:, top_feature_idx]
lr_final = LogisticRegression(max_iter=1000, random_state=42)
lr_final.fit(X_top, y)
biomarker_scores = lr_final.predict_proba(X_top)[:, 1]

score_df = pd.DataFrame({
    'SampleID': meta_df['SampleID'],
    'Group': groups,
    'Biomarker_Score': biomarker_scores.round(4)
})
score_df.to_csv('data/biomarker_scores.csv', index=False)

# Figure 6: Biomarker analysis
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Feature importance
top_feat_plot = top_feat_df.head(20).sort_values('Importance', ascending=True)
colors_fi = ['#e74c3c' if t == 'Taxa' else '#3498db' for t in top_feat_plot['Type']]
axes[0, 0].barh(top_feat_plot['Feature'], top_feat_plot['Importance'], color=colors_fi)
axes[0, 0].set_xlabel('Feature Importance')
axes[0, 0].set_title('Top 20 Biomarker Features')
axes[0, 0].legend(handles=[Patch(color='#e74c3c', label='Taxa'), Patch(color='#3498db', label='Metabolite')])

# ROC curves
colors_roc = ['#e74c3c', '#3498db', '#2ecc71', '#9b59b6', '#e67e22']
for (name, (fpr, tpr)), color in zip(roc_data.items(), colors_roc):
    auc_val = auc_results[name]['mean_AUC']
    axes[0, 1].plot(fpr, tpr, color=color, label=f'{name} (AUC={auc_val:.3f})')
axes[0, 1].plot([0, 1], [0, 1], 'k--', alpha=0.3)
axes[0, 1].set_xlabel('False Positive Rate')
axes[0, 1].set_ylabel('True Positive Rate')
axes[0, 1].set_title('ROC Curves: Model Comparison')
axes[0, 1].legend(fontsize=8)

# Biomarker score distribution
for grp, color in [('IBD', '#e74c3c'), ('Control', '#3498db')]:
    scores_grp = score_df.loc[score_df['Group'] == grp, 'Biomarker_Score']
    axes[1, 0].hist(scores_grp, bins=20, alpha=0.6, color=color, label=grp, density=True)
axes[1, 0].set_xlabel('Integrated Biomarker Score')
axes[1, 0].set_ylabel('Density')
axes[1, 0].set_title('Biomarker Score Distribution')
axes[1, 0].legend()

# Confusion matrix (best model on full data)
rf_final = RandomForestClassifier(n_estimators=200, random_state=42, max_depth=10)
rf_final.fit(X_combined, y)
y_pred = rf_final.predict(X_combined)
cm = confusion_matrix(y, y_pred)
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[1, 1],
            xticklabels=['Control', 'IBD'], yticklabels=['Control', 'IBD'])
axes[1, 1].set_xlabel('Predicted')
axes[1, 1].set_ylabel('True')
axes[1, 1].set_title('Confusion Matrix (RF Integrated)')

plt.tight_layout()
plt.savefig('figures/fig6_biomarker_scoring.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/fig6_biomarker_scoring.png")

# ==================================================
# Module 7: IBD Case Study Summary
# ==================================================
print("\n--- Module 7: IBD Case Study Summary ---")

# Key IBD-associated features
ibd_taxa_changes = diff_taxa_df.sort_values('p_value').head(10)
ibd_met_changes = diff_df.sort_values('p_value').head(10)

print("\nTop IBD-associated taxa:")
print(ibd_taxa_changes[['Taxon', 'MeanDiff', 'FDR_q']].to_string(index=False))
print("\nTop IBD-associated metabolites:")
print(ibd_met_changes[['Metabolite', 'MeanDiff', 'FDR_q']].to_string(index=False))

# Figure 7: IBD case study
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Top differential taxa
ibd_taxa_plot = diff_taxa_df.sort_values('MeanDiff').head(10)
ibd_taxa_plot2 = diff_taxa_df.sort_values('MeanDiff', ascending=False).head(10)
combined_taxa = pd.concat([ibd_taxa_plot, ibd_taxa_plot2]).sort_values('MeanDiff')
colors_t = ['#e74c3c' if d > 0 else '#3498db' for d in combined_taxa['MeanDiff']]
axes[0, 0].barh(combined_taxa['Taxon'], combined_taxa['MeanDiff'], color=colors_t)
axes[0, 0].set_xlabel('Mean Difference (IBD - Control)')
axes[0, 0].set_title('Top Differential Taxa in IBD')
axes[0, 0].axvline(0, color='black', linewidth=0.5)
axes[0, 0].tick_params(axis='y', labelsize=8)

# Top differential metabolites
ibd_met_plot = diff_df.sort_values('MeanDiff').head(10)
ibd_met_plot2 = diff_df.sort_values('MeanDiff', ascending=False).head(10)
combined_met = pd.concat([ibd_met_plot, ibd_met_plot2]).sort_values('MeanDiff')
colors_m = ['#e74c3c' if d > 0 else '#3498db' for d in combined_met['MeanDiff']]
axes[0, 1].barh(combined_met['Metabolite'], combined_met['MeanDiff'], color=colors_m)
axes[0, 1].set_xlabel('Mean Difference (IBD - Control)')
axes[0, 1].set_title('Top Differential Metabolites in IBD')
axes[0, 1].axvline(0, color='black', linewidth=0.5)
axes[0, 1].tick_params(axis='y', labelsize=8)

# Boxplots of key metabolites
key_mets = ['Butyrate', 'Tryptophan', 'Kynurenine', 'Indoxyl_sulfate']
data_box = []
for m in key_mets:
    if m in met_df.columns:
        for i, g in enumerate(groups):
            data_box.append({'Metabolite': m, 'Group': g, 'Value': met_df.iloc[i][m]})
box_df = pd.DataFrame(data_box)
sns.boxplot(data=box_df, x='Metabolite', y='Value', hue='Group',
            palette={'IBD': '#e74c3c', 'Control': '#3498db'}, ax=axes[1, 0])
axes[1, 0].set_title('Key Metabolites: IBD vs Control')
axes[1, 0].tick_params(axis='x', rotation=30)

# Key taxa boxplots
key_taxa = ['Faecalibacterium', 'Roseburia', 'Escherichia', 'Fusobacterium']
data_box_t = []
for t in key_taxa:
    if t in taxa_df.columns:
        for i, g in enumerate(groups):
            data_box_t.append({'Taxon': t, 'Group': g, 'Value': taxa_df.iloc[i][t]})
box_t_df = pd.DataFrame(data_box_t)
sns.boxplot(data=box_t_df, x='Taxon', y='Value', hue='Group',
            palette={'IBD': '#e74c3c', 'Control': '#3498db'}, ax=axes[1, 1])
axes[1, 1].set_title('Key Taxa: IBD vs Control')

plt.tight_layout()
plt.savefig('figures/fig7_ibd_case_study.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/fig7_ibd_case_study.png")

# ==================================================
# Summary Statistics
# ==================================================
print("\n" + "=" * 60)
print("PIPELINE COMPLETE - SUMMARY")
print("=" * 60)
N_SAMPLES = len(groups)
print(f"Samples: {N_SAMPLES} ({sum(groups=='IBD')} IBD, {sum(groups=='Control')} Control)")
print(f"Metabolites: {met_df.shape[1]} (Annotated: {sum(annotation_levels != 'Level4_unknown')})")
print(f"Taxa: {taxa_df.shape[1]}")
print(f"Significant correlations (|r|>0.3, FDR<0.05): {len(sig_edges)}")
print(f"Differential metabolites (FDR<0.05): {len(sig_mets)}")
print(f"Differential taxa (FDR<0.05): {len(sig_taxa)}")
n_enriched = len(enrich_df[enrich_df['FDR_q'] < 0.05]) if len(enrich_df) > 0 else 0
print(f"Enriched pathways (FDR<0.05): {n_enriched}")
print(f"Best model AUC: {max(v['mean_AUC'] for v in auc_results.values()):.4f}")
print(f"MR IVW beta: {beta_ivw:.4f}, p={p_ivw:.6f}")

# Save summary
summary = {
    'n_samples': int(N_SAMPLES),
    'n_ibd': int(sum(groups == 'IBD')),
    'n_control': int(sum(groups == 'Control')),
    'n_metabolites': int(met_df.shape[1]),
    'n_taxa': int(taxa_df.shape[1]),
    'n_sig_correlations': len(sig_edges),
    'n_diff_metabolites': len(sig_mets),
    'n_diff_taxa': len(sig_taxa),
    'n_enriched_pathways': n_enriched,
    'best_auc': round(max(v['mean_AUC'] for v in auc_results.values()), 4),
    'auc_results': {k: {kk: round(vv, 4) for kk, vv in v.items()} for k, v in auc_results.items()},
    'mr_beta': round(beta_ivw, 4),
    'mr_pvalue': round(p_ivw, 6),
    'melonnpan_results': melonnpan_df.to_dict('records'),
    'top_biomarkers': top_feat_df.head(10).to_dict('records'),
}
with open('data/analysis_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print("\nAll results saved to data/ and figures/")
