#!/usr/bin/env python3
"""
Open Access / Open Data Impact Analysis Framework
===================================================
Quantitative analysis of the impact of open access and open data
on research communities using simulated bibliometric/altmetric data.

Modules:
1. OA Citation Advantage (OACA) - Causal estimation via PSM + DiD
2. Data Sharing & Reuse Patterns
3. Preprint Server Role Evaluation
4. FAIR Compliance Automated Assessment
5. Citizen Science & Outreach Impact
6. Life Sciences Open Data Case Study
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
FIGURES_DIR = Path(__file__).parent.parent / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

plt.rcParams.update({
    'figure.figsize': (10, 6),
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'figure.dpi': 150,
})

results = {}

# =============================================================================
# Module 1: OA Citation Advantage (OACA) — Causal Estimation
# =============================================================================
print("=" * 70)
print("Module 1: OA Citation Advantage (OACA) — PSM + DiD Estimation")
print("=" * 70)

N_PAPERS = 5000
years = np.random.choice(range(2015, 2025), N_PAPERS)
journal_prestige = np.random.uniform(0.5, 5.0, N_PAPERS)  # IF proxy
author_hindex = np.random.poisson(15, N_PAPERS)
num_authors = np.random.poisson(4, N_PAPERS) + 1
field = np.random.choice(['BioMed', 'Physics', 'CS', 'SocialSci', 'Engineering'], N_PAPERS)

# OA assignment with confounding
logit_oa = -1.0 + 0.15 * journal_prestige + 0.02 * author_hindex + 0.1 * (years - 2015)
prob_oa = 1 / (1 + np.exp(-logit_oa))
oa_status = np.random.binomial(1, prob_oa)

# Citation generation with OA causal effect
base_citations = (
    0.5 * journal_prestige ** 1.5 +
    0.3 * author_hindex +
    0.2 * num_authors +
    0.5 * (years - 2015)
)
oa_effect = oa_status * np.random.uniform(3.0, 8.0, N_PAPERS)  # True causal effect
noise = np.random.exponential(2.0, N_PAPERS)
citations = np.maximum(0, base_citations + oa_effect + noise).astype(int)

df_papers = pd.DataFrame({
    'year': years, 'journal_prestige': journal_prestige,
    'author_hindex': author_hindex, 'num_authors': num_authors,
    'field': field, 'oa_status': oa_status, 'citations': citations,
})

# Propensity Score Matching
X_psm = df_papers[['journal_prestige', 'author_hindex', 'num_authors', 'year']].values
y_psm = df_papers['oa_status'].values

lr = LogisticRegression(max_iter=1000)
lr.fit(X_psm, y_psm)
propensity_scores = lr.predict_proba(X_psm)[:, 1]
df_papers['propensity_score'] = propensity_scores

# Nearest-neighbor matching
treated_idx = np.where(y_psm == 1)[0]
control_idx = np.where(y_psm == 0)[0]

nn = NearestNeighbors(n_neighbors=1, metric='euclidean')
nn.fit(propensity_scores[control_idx].reshape(-1, 1))
distances, indices = nn.kneighbors(propensity_scores[treated_idx].reshape(-1, 1))
matched_control_idx = control_idx[indices.flatten()]

treated_citations = df_papers.iloc[treated_idx]['citations'].values
matched_control_citations = df_papers.iloc[matched_control_idx]['citations'].values

att = np.mean(treated_citations - matched_control_citations)
att_se = np.std(treated_citations - matched_control_citations) / np.sqrt(len(treated_idx))
att_ci = (att - 1.96 * att_se, att + 1.96 * att_se)
att_pct = att / np.mean(matched_control_citations) * 100

print(f"  ATT (Average Treatment Effect on Treated): {att:.2f} citations")
print(f"  95% CI: [{att_ci[0]:.2f}, {att_ci[1]:.2f}]")
print(f"  Percentage advantage: {att_pct:.1f}%")
print(f"  OA papers mean: {np.mean(treated_citations):.1f}, Matched control mean: {np.mean(matched_control_citations):.1f}")

results['oaca'] = {
    'ATT': round(att, 2), 'CI_lower': round(att_ci[0], 2), 'CI_upper': round(att_ci[1], 2),
    'pct_advantage': round(att_pct, 1),
    'n_treated': len(treated_idx), 'n_control': len(matched_control_idx),
}

# DiD analysis by year
did_results = []
for yr in range(2016, 2025):
    pre = df_papers[(df_papers['year'] < yr)]
    post = df_papers[(df_papers['year'] >= yr)]
    for period, subset in [('pre', pre), ('post', post)]:
        oa_mean = subset[subset['oa_status'] == 1]['citations'].mean()
        non_oa_mean = subset[subset['oa_status'] == 0]['citations'].mean()
        did_results.append({'cutoff_year': yr, 'period': period,
                           'oa_mean': oa_mean, 'non_oa_mean': non_oa_mean,
                           'diff': oa_mean - non_oa_mean})

df_did = pd.DataFrame(did_results)

# Figure 1: OACA by field
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

field_advantages = df_papers.groupby(['field', 'oa_status'])['citations'].mean().unstack()
field_advantages['advantage'] = field_advantages[1] - field_advantages[0]
field_advantages['pct'] = field_advantages['advantage'] / field_advantages[0] * 100

colors = sns.color_palette("Set2", len(field_advantages))
bars = axes[0].bar(field_advantages.index, field_advantages['pct'], color=colors, edgecolor='black', linewidth=0.5)
axes[0].set_ylabel('OA Citation Advantage (%)')
axes[0].set_xlabel('Research Field')
axes[0].set_title('(a) OA Citation Advantage by Field')
axes[0].axhline(y=att_pct, color='red', linestyle='--', alpha=0.7, label=f'Overall: {att_pct:.1f}%')
axes[0].legend()
for bar, val in zip(bars, field_advantages['pct']):
    axes[0].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.5,
                f'{val:.1f}%', ha='center', va='bottom', fontsize=9)

# Propensity score distribution
axes[1].hist(propensity_scores[y_psm == 1], bins=40, alpha=0.6, label='OA', density=True, color='steelblue')
axes[1].hist(propensity_scores[y_psm == 0], bins=40, alpha=0.6, label='Non-OA', density=True, color='salmon')
axes[1].set_xlabel('Propensity Score')
axes[1].set_ylabel('Density')
axes[1].set_title('(b) Propensity Score Distribution')
axes[1].legend()

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig1_oaca_analysis.png', bbox_inches='tight')
plt.close()
print("  → Saved fig1_oaca_analysis.png")

# Figure 2: DiD over time
fig, ax = plt.subplots(figsize=(10, 6))
yearly_stats = df_papers.groupby(['year', 'oa_status'])['citations'].agg(['mean', 'std', 'count']).reset_index()
for oa_val, label, color in [(1, 'Open Access', 'steelblue'), (0, 'Non-OA', 'salmon')]:
    subset = yearly_stats[yearly_stats['oa_status'] == oa_val]
    se = subset['std'] / np.sqrt(subset['count'])
    ax.plot(subset['year'], subset['mean'], 'o-', label=label, color=color, linewidth=2)
    ax.fill_between(subset['year'], subset['mean'] - 1.96*se, subset['mean'] + 1.96*se,
                    alpha=0.2, color=color)

ax.set_xlabel('Publication Year')
ax.set_ylabel('Mean Citations')
ax.set_title('Difference-in-Differences: OA vs Non-OA Citation Trajectories')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig2_did_trajectories.png', bbox_inches='tight')
plt.close()
print("  → Saved fig2_did_trajectories.png")


# =============================================================================
# Module 2: Data Sharing & Reuse Patterns
# =============================================================================
print("\n" + "=" * 70)
print("Module 2: Data Sharing & Reuse Patterns")
print("=" * 70)

N_DATASETS = 2000
ds_years = np.random.choice(range(2015, 2025), N_DATASETS)
ds_field = np.random.choice(['Genomics', 'Proteomics', 'Imaging', 'Clinical', 'Ecology'], N_DATASETS,
                            p=[0.35, 0.15, 0.15, 0.2, 0.15])
has_metadata = np.random.binomial(1, 0.7, N_DATASETS)
has_code = np.random.binomial(1, 0.4, N_DATASETS)
license_type = np.random.choice(['CC-BY', 'CC0', 'Restricted', 'None'], N_DATASETS,
                                 p=[0.35, 0.25, 0.25, 0.15])

# Reuse count depends on metadata, code, license, and field
base_reuse = np.random.poisson(3, N_DATASETS)
reuse_boost = (
    has_metadata * np.random.poisson(4, N_DATASETS) +
    has_code * np.random.poisson(3, N_DATASETS) +
    (np.isin(license_type, ['CC-BY', 'CC0'])).astype(int) * np.random.poisson(5, N_DATASETS)
)
reuse_count = base_reuse + reuse_boost
time_to_first_reuse = np.random.exponential(12, N_DATASETS)  # months

df_datasets = pd.DataFrame({
    'year': ds_years, 'field': ds_field, 'has_metadata': has_metadata,
    'has_code': has_code, 'license': license_type,
    'reuse_count': reuse_count, 'time_to_first_reuse_months': time_to_first_reuse,
})

# Normalized Reusability Index (NRI)
age_years = 2025 - df_datasets['year']
df_datasets['NRI'] = df_datasets['reuse_count'] / np.maximum(age_years, 1)

# Statistics
reuse_by_license = df_datasets.groupby('license')['reuse_count'].agg(['mean', 'median', 'std'])
reuse_by_meta = df_datasets.groupby('has_metadata')['reuse_count'].mean()
reuse_by_code = df_datasets.groupby('has_code')['reuse_count'].mean()

print(f"  Mean reuse count: {df_datasets['reuse_count'].mean():.1f}")
print(f"  With metadata: {reuse_by_meta[1]:.1f} vs without: {reuse_by_meta[0]:.1f}")
print(f"  With code: {reuse_by_code[1]:.1f} vs without: {reuse_by_code[0]:.1f}")
print(f"  CC-BY/CC0 mean reuse: {df_datasets[df_datasets['license'].isin(['CC-BY','CC0'])]['reuse_count'].mean():.1f}")
print(f"  Restricted/None mean reuse: {df_datasets[df_datasets['license'].isin(['Restricted','None'])]['reuse_count'].mean():.1f}")

results['data_sharing'] = {
    'mean_reuse': round(df_datasets['reuse_count'].mean(), 1),
    'metadata_effect': round(reuse_by_meta[1] - reuse_by_meta[0], 1),
    'code_effect': round(reuse_by_code[1] - reuse_by_code[0], 1),
    'open_license_mean': round(df_datasets[df_datasets['license'].isin(['CC-BY','CC0'])]['reuse_count'].mean(), 1),
}

# Figure 3: Data sharing patterns
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (a) Reuse by license
license_order = ['CC0', 'CC-BY', 'Restricted', 'None']
palette = {'CC0': '#2ecc71', 'CC-BY': '#3498db', 'Restricted': '#e74c3c', 'None': '#95a5a6'}
sns.boxplot(data=df_datasets, x='license', y='reuse_count', order=license_order,
            palette=palette, ax=axes[0, 0])
axes[0, 0].set_title('(a) Dataset Reuse by License Type')
axes[0, 0].set_ylabel('Reuse Count')

# (b) Reuse by field
field_order = df_datasets.groupby('field')['reuse_count'].mean().sort_values(ascending=False).index
sns.violinplot(data=df_datasets, x='field', y='NRI', order=field_order,
               palette='Set2', ax=axes[0, 1])
axes[0, 1].set_title('(b) Normalized Reusability Index by Field')
axes[0, 1].set_ylabel('NRI (reuses/year)')
axes[0, 1].tick_params(axis='x', rotation=15)

# (c) Metadata & code impact
categories = ['No Meta\nNo Code', 'Meta\nNo Code', 'No Meta\nWith Code', 'Meta\n& Code']
means = [
    df_datasets[(df_datasets['has_metadata']==0) & (df_datasets['has_code']==0)]['reuse_count'].mean(),
    df_datasets[(df_datasets['has_metadata']==1) & (df_datasets['has_code']==0)]['reuse_count'].mean(),
    df_datasets[(df_datasets['has_metadata']==0) & (df_datasets['has_code']==1)]['reuse_count'].mean(),
    df_datasets[(df_datasets['has_metadata']==1) & (df_datasets['has_code']==1)]['reuse_count'].mean(),
]
bar_colors = ['#e74c3c', '#f39c12', '#3498db', '#2ecc71']
bars = axes[1, 0].bar(categories, means, color=bar_colors, edgecolor='black', linewidth=0.5)
axes[1, 0].set_ylabel('Mean Reuse Count')
axes[1, 0].set_title('(c) Impact of Metadata & Code Availability')
for bar, val in zip(bars, means):
    axes[1, 0].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.3,
                   f'{val:.1f}', ha='center', va='bottom', fontsize=10)

# (d) Temporal trend
yearly_sharing = df_datasets.groupby('year').agg(
    mean_reuse=('reuse_count', 'mean'),
    pct_open_license=('license', lambda x: np.mean(np.isin(x, ['CC-BY', 'CC0'])) * 100),
    pct_metadata=('has_metadata', 'mean'),
).reset_index()
ax2 = axes[1, 1].twinx()
axes[1, 1].plot(yearly_sharing['year'], yearly_sharing['mean_reuse'], 'o-', color='steelblue', linewidth=2, label='Mean Reuse')
ax2.plot(yearly_sharing['year'], yearly_sharing['pct_open_license'], 's--', color='#2ecc71', linewidth=2, label='% Open License')
axes[1, 1].set_xlabel('Year')
axes[1, 1].set_ylabel('Mean Reuse Count', color='steelblue')
ax2.set_ylabel('% Open License', color='#2ecc71')
axes[1, 1].set_title('(d) Temporal Trends in Data Sharing')
lines1, labels1 = axes[1, 1].get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
axes[1, 1].legend(lines1 + lines2, labels1 + labels2, loc='upper left')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig3_data_sharing_patterns.png', bbox_inches='tight')
plt.close()
print("  → Saved fig3_data_sharing_patterns.png")


# =============================================================================
# Module 3: Preprint Server Role Evaluation
# =============================================================================
print("\n" + "=" * 70)
print("Module 3: Preprint Server Role Evaluation")
print("=" * 70)

N_PREPRINTS = 3000
preprint_server = np.random.choice(['bioRxiv', 'arXiv', 'medRxiv', 'SSRN', 'Other'], N_PREPRINTS,
                                    p=[0.35, 0.30, 0.15, 0.10, 0.10])
submit_year = np.random.choice(range(2018, 2025), N_PREPRINTS)
peer_review_days_base = np.random.lognormal(4.8, 0.5, N_PREPRINTS)  # ~120 days median
preprint_first = np.random.binomial(1, 0.6, N_PREPRINTS)

# Preprint-first reduces review time
review_reduction = preprint_first * np.random.uniform(15, 40, N_PREPRINTS)
peer_review_days = np.maximum(30, peer_review_days_base - review_reduction)
published = np.random.binomial(1, 0.75, N_PREPRINTS)
preprint_citations = np.random.poisson(5, N_PREPRINTS)
final_citations = np.where(published, preprint_citations + np.random.poisson(10, N_PREPRINTS),
                           preprint_citations)
community_comments = np.random.poisson(2, N_PREPRINTS)

df_preprints = pd.DataFrame({
    'server': preprint_server, 'year': submit_year,
    'preprint_first': preprint_first, 'review_days': peer_review_days,
    'published': published, 'preprint_citations': preprint_citations,
    'final_citations': final_citations, 'comments': community_comments,
})

# Analysis
mean_review_preprint_first = df_preprints[df_preprints['preprint_first'] == 1]['review_days'].mean()
mean_review_journal_first = df_preprints[df_preprints['preprint_first'] == 0]['review_days'].mean()
review_reduction_pct = (mean_review_journal_first - mean_review_preprint_first) / mean_review_journal_first * 100

t_stat, p_val = stats.ttest_ind(
    df_preprints[df_preprints['preprint_first'] == 1]['review_days'],
    df_preprints[df_preprints['preprint_first'] == 0]['review_days']
)

print(f"  Preprint-first review time: {mean_review_preprint_first:.0f} days")
print(f"  Journal-first review time: {mean_review_journal_first:.0f} days")
print(f"  Reduction: {review_reduction_pct:.1f}% (p={p_val:.2e})")
print(f"  Publication rate: {df_preprints['published'].mean()*100:.1f}%")

results['preprint'] = {
    'preprint_first_days': round(mean_review_preprint_first, 0),
    'journal_first_days': round(mean_review_journal_first, 0),
    'reduction_pct': round(review_reduction_pct, 1),
    'p_value': f"{p_val:.2e}",
    'publication_rate': round(df_preprints['published'].mean() * 100, 1),
}

# Figure 4: Preprint analysis
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (a) Review time comparison
sns.boxplot(data=df_preprints, x='preprint_first', y='review_days', ax=axes[0, 0],
            palette=['salmon', 'steelblue'])
axes[0, 0].set_xticklabels(['Journal-first', 'Preprint-first'])
axes[0, 0].set_ylabel('Peer Review Duration (days)')
axes[0, 0].set_title(f'(a) Review Duration (p={p_val:.2e})')
axes[0, 0].set_xlabel('')

# (b) By server
server_stats = df_preprints.groupby('server').agg(
    mean_review=('review_days', 'mean'),
    pub_rate=('published', 'mean'),
    mean_citations=('final_citations', 'mean')
).reset_index()
bars = axes[0, 1].bar(server_stats['server'], server_stats['mean_review'],
                       color=sns.color_palette("Set2", len(server_stats)), edgecolor='black', linewidth=0.5)
axes[0, 1].set_ylabel('Mean Review Time (days)')
axes[0, 1].set_title('(b) Review Duration by Preprint Server')
for bar, val in zip(bars, server_stats['mean_review']):
    axes[0, 1].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1,
                   f'{val:.0f}', ha='center', va='bottom', fontsize=9)

# (c) Temporal trend of preprint adoption
yearly_preprint = df_preprints.groupby('year').agg(
    n_preprints=('server', 'count'),
    pct_preprint_first=('preprint_first', 'mean'),
    mean_comments=('comments', 'mean')
).reset_index()
axes[1, 0].bar(yearly_preprint['year'], yearly_preprint['n_preprints'], alpha=0.6, color='steelblue', label='Count')
ax_twin = axes[1, 0].twinx()
ax_twin.plot(yearly_preprint['year'], yearly_preprint['pct_preprint_first']*100, 'ro-', linewidth=2, label='% Preprint-first')
axes[1, 0].set_xlabel('Year')
axes[1, 0].set_ylabel('Number of Preprints', color='steelblue')
ax_twin.set_ylabel('% Preprint-first', color='red')
axes[1, 0].set_title('(c) Preprint Adoption Trends')
lines1, labels1 = axes[1, 0].get_legend_handles_labels()
lines2, labels2 = ax_twin.get_legend_handles_labels()
axes[1, 0].legend(lines1 + lines2, labels1 + labels2, loc='upper left')

# (d) Citation advantage
cite_data = df_preprints.groupby(['preprint_first', 'published'])['final_citations'].mean().unstack()
cite_data.index = ['Journal-first', 'Preprint-first']
cite_data.columns = ['Unpublished', 'Published']
cite_data.plot(kind='bar', ax=axes[1, 1], color=['#e74c3c', '#2ecc71'], edgecolor='black', linewidth=0.5)
axes[1, 1].set_title('(d) Citation Impact by Preprint Strategy')
axes[1, 1].set_ylabel('Mean Citations')
axes[1, 1].tick_params(axis='x', rotation=0)
axes[1, 1].legend(title='Status')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig4_preprint_analysis.png', bbox_inches='tight')
plt.close()
print("  → Saved fig4_preprint_analysis.png")


# =============================================================================
# Module 4: FAIR Compliance Automated Assessment
# =============================================================================
print("\n" + "=" * 70)
print("Module 4: FAIR Compliance Automated Assessment")
print("=" * 70)

N_REPOS = 1500
repo_names = [f"repo_{i:04d}" for i in range(N_REPOS)]
repo_year = np.random.choice(range(2018, 2025), N_REPOS)
repo_field = np.random.choice(['Genomics', 'Climate', 'Physics', 'Chemistry', 'Social'], N_REPOS)

# FAIR sub-scores (0-1 scale)
f_score = np.clip(np.random.beta(3, 2, N_REPOS) + 0.05 * (repo_year - 2018), 0, 1)  # Findable
a_score = np.clip(np.random.beta(4, 2, N_REPOS) + 0.03 * (repo_year - 2018), 0, 1)  # Accessible
i_score = np.clip(np.random.beta(2, 3, N_REPOS) + 0.04 * (repo_year - 2018), 0, 1)  # Interoperable
r_score = np.clip(np.random.beta(2, 3, N_REPOS) + 0.03 * (repo_year - 2018), 0, 1)  # Reusable

fair_total = (f_score + a_score + i_score + r_score) / 4

# Sub-metrics
has_pid = (f_score > 0.5).astype(int)
has_standard_metadata = (f_score > 0.6).astype(int)
has_open_protocol = (a_score > 0.5).astype(int)
uses_vocab = (i_score > 0.4).astype(int)
has_license_file = (r_score > 0.4).astype(int)
has_provenance = (r_score > 0.5).astype(int)

df_fair = pd.DataFrame({
    'repo': repo_names, 'year': repo_year, 'field': repo_field,
    'F': f_score, 'A': a_score, 'I': i_score, 'R': r_score,
    'FAIR_total': fair_total,
    'has_pid': has_pid, 'has_standard_metadata': has_standard_metadata,
    'has_open_protocol': has_open_protocol, 'uses_vocab': uses_vocab,
    'has_license': has_license_file, 'has_provenance': has_provenance,
})

# Summary statistics
fair_means = df_fair[['F', 'A', 'I', 'R', 'FAIR_total']].mean()
print(f"  Mean FAIR score: {fair_means['FAIR_total']:.3f}")
print(f"  F={fair_means['F']:.3f}, A={fair_means['A']:.3f}, I={fair_means['I']:.3f}, R={fair_means['R']:.3f}")

fair_by_year = df_fair.groupby('year')[['F', 'A', 'I', 'R', 'FAIR_total']].mean()
fair_by_field = df_fair.groupby('field')[['FAIR_total']].mean().sort_values('FAIR_total', ascending=False)
print(f"  Top field: {fair_by_field.index[0]} ({fair_by_field.iloc[0, 0]:.3f})")

results['fair'] = {
    'mean_F': round(fair_means['F'], 3), 'mean_A': round(fair_means['A'], 3),
    'mean_I': round(fair_means['I'], 3), 'mean_R': round(fair_means['R'], 3),
    'mean_total': round(fair_means['FAIR_total'], 3),
}

# Figure 5: FAIR assessment
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (a) Radar chart - FAIR by field
from matplotlib.patches import FancyBboxPatch
fields_unique = df_fair['field'].unique()
fair_by_field_detail = df_fair.groupby('field')[['F', 'A', 'I', 'R']].mean()
x_pos = np.arange(len(['F', 'A', 'I', 'R']))
width = 0.15
for idx, fld in enumerate(fields_unique):
    vals = fair_by_field_detail.loc[fld].values
    axes[0, 0].bar(x_pos + idx * width, vals, width, label=fld, alpha=0.8)
axes[0, 0].set_xticks(x_pos + width * 2)
axes[0, 0].set_xticklabels(['Findable', 'Accessible', 'Interoperable', 'Reusable'])
axes[0, 0].set_ylabel('Score (0-1)')
axes[0, 0].set_title('(a) FAIR Scores by Domain')
axes[0, 0].legend(fontsize=8, ncol=2)
axes[0, 0].set_ylim(0, 1)

# (b) Temporal improvement
for col, label in [('F', 'Findable'), ('A', 'Accessible'), ('I', 'Interoperable'), ('R', 'Reusable')]:
    axes[0, 1].plot(fair_by_year.index, fair_by_year[col], 'o-', label=label, linewidth=2)
axes[0, 1].set_xlabel('Year')
axes[0, 1].set_ylabel('Mean Score')
axes[0, 1].set_title('(b) FAIR Score Trends Over Time')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

# (c) Distribution of total FAIR score
axes[1, 0].hist(df_fair['FAIR_total'], bins=50, color='steelblue', edgecolor='black', linewidth=0.5, alpha=0.7)
axes[1, 0].axvline(df_fair['FAIR_total'].mean(), color='red', linestyle='--', linewidth=2,
                   label=f"Mean={df_fair['FAIR_total'].mean():.3f}")
axes[1, 0].axvline(df_fair['FAIR_total'].median(), color='orange', linestyle='--', linewidth=2,
                   label=f"Median={df_fair['FAIR_total'].median():.3f}")
axes[1, 0].set_xlabel('FAIR Total Score')
axes[1, 0].set_ylabel('Count')
axes[1, 0].set_title('(c) Distribution of FAIR Compliance Scores')
axes[1, 0].legend()

# (d) Sub-metric compliance rates
metrics = ['has_pid', 'has_standard_metadata', 'has_open_protocol', 'uses_vocab', 'has_license', 'has_provenance']
metric_labels = ['PID', 'Std Metadata', 'Open Protocol', 'Std Vocab', 'License', 'Provenance']
rates = [df_fair[m].mean() * 100 for m in metrics]
bar_colors = plt.cm.RdYlGn(np.array(rates) / 100)
bars = axes[1, 1].barh(metric_labels, rates, color=bar_colors, edgecolor='black', linewidth=0.5)
axes[1, 1].set_xlabel('Compliance Rate (%)')
axes[1, 1].set_title('(d) FAIR Sub-metric Compliance Rates')
for bar, val in zip(bars, rates):
    axes[1, 1].text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2.,
                   f'{val:.1f}%', ha='left', va='center', fontsize=9)
axes[1, 1].set_xlim(0, 110)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig5_fair_assessment.png', bbox_inches='tight')
plt.close()
print("  → Saved fig5_fair_assessment.png")


# =============================================================================
# Module 5: Citizen Science & Outreach Impact
# =============================================================================
print("\n" + "=" * 70)
print("Module 5: Citizen Science & Outreach Impact")
print("=" * 70)

N_PROJECTS = 500
cs_field = np.random.choice(['Ecology', 'Astronomy', 'Health', 'Climate', 'Biodiversity'], N_PROJECTS)
n_participants = np.random.lognormal(6, 1.5, N_PROJECTS).astype(int)
n_publications = np.random.poisson(3, N_PROJECTS)
twitter_mentions = np.random.poisson(50, N_PROJECTS) * (1 + 0.3 * np.log1p(n_participants) / 10)
news_mentions = np.random.poisson(5, N_PROJECTS)
blog_mentions = np.random.poisson(8, N_PROJECTS)
policy_mentions = np.random.poisson(1, N_PROJECTS)
academic_citations = np.random.poisson(15, N_PROJECTS)
mendeley_readers = np.random.poisson(30, N_PROJECTS) * (1 + 0.2 * np.log1p(n_participants) / 10)

# Altmetric Attention Score (simplified)
altmetric_score = (twitter_mentions * 1 + news_mentions * 8 + blog_mentions * 5 +
                   policy_mentions * 10 + mendeley_readers * 0.5)

df_citizen = pd.DataFrame({
    'field': cs_field, 'participants': n_participants.astype(int),
    'publications': n_publications, 'twitter': twitter_mentions.astype(int),
    'news': news_mentions, 'blogs': blog_mentions, 'policy': policy_mentions,
    'citations': academic_citations, 'mendeley': mendeley_readers.astype(int),
    'altmetric_score': altmetric_score,
})

# Correlation analysis
corr_part_cite = stats.spearmanr(df_citizen['participants'], df_citizen['citations'])
corr_part_alt = stats.spearmanr(df_citizen['participants'], df_citizen['altmetric_score'])

print(f"  Mean participants: {df_citizen['participants'].mean():.0f}")
print(f"  Spearman(participants, citations): r={corr_part_cite.statistic:.3f}, p={corr_part_cite.pvalue:.2e}")
print(f"  Spearman(participants, altmetric): r={corr_part_alt.statistic:.3f}, p={corr_part_alt.pvalue:.2e}")
print(f"  Mean altmetric score: {df_citizen['altmetric_score'].mean():.0f}")

results['citizen_science'] = {
    'mean_participants': int(df_citizen['participants'].mean()),
    'corr_participants_citations': round(corr_part_cite.statistic, 3),
    'corr_participants_altmetric': round(corr_part_alt.statistic, 3),
    'mean_altmetric': round(df_citizen['altmetric_score'].mean(), 0),
}

# Figure 6: Citizen science impact
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (a) Participants vs Altmetric
scatter = axes[0, 0].scatter(np.log10(df_citizen['participants']), np.log10(df_citizen['altmetric_score']),
                             c=df_citizen['citations'], cmap='viridis', alpha=0.5, s=20)
plt.colorbar(scatter, ax=axes[0, 0], label='Academic Citations')
axes[0, 0].set_xlabel('log₁₀(Participants)')
axes[0, 0].set_ylabel('log₁₀(Altmetric Score)')
axes[0, 0].set_title(f'(a) Participation vs Outreach (ρ={corr_part_alt.statistic:.3f})')

# (b) Impact channels by field
channels = df_citizen.groupby('field')[['twitter', 'news', 'blogs', 'policy']].mean()
channels.plot(kind='bar', ax=axes[0, 1], color=['#1DA1F2', '#e74c3c', '#f39c12', '#2ecc71'],
              edgecolor='black', linewidth=0.5)
axes[0, 1].set_title('(b) Outreach Channels by Field')
axes[0, 1].set_ylabel('Mean Mentions')
axes[0, 1].tick_params(axis='x', rotation=30)
axes[0, 1].legend(fontsize=9)

# (c) Bibliometric vs Altmetric comparison
axes[1, 0].scatter(df_citizen['citations'], df_citizen['altmetric_score'],
                   alpha=0.4, s=20, color='steelblue')
z = np.polyfit(df_citizen['citations'], df_citizen['altmetric_score'], 1)
p = np.poly1d(z)
x_line = np.linspace(df_citizen['citations'].min(), df_citizen['citations'].max(), 100)
axes[1, 0].plot(x_line, p(x_line), 'r--', linewidth=2, label=f'y={z[0]:.1f}x+{z[1]:.0f}')
axes[1, 0].set_xlabel('Academic Citations')
axes[1, 0].set_ylabel('Altmetric Score')
axes[1, 0].set_title('(c) Bibliometrics vs Altmetrics')
axes[1, 0].legend()

# (d) Distribution of impact
impact_data = pd.DataFrame({
    'Metric': ['Citations', 'Twitter', 'News', 'Blogs', 'Policy', 'Mendeley'],
    'Mean': [df_citizen['citations'].mean(), df_citizen['twitter'].mean(),
             df_citizen['news'].mean(), df_citizen['blogs'].mean(),
             df_citizen['policy'].mean(), df_citizen['mendeley'].mean()],
})
# Normalize for comparison
impact_data['Normalized'] = impact_data['Mean'] / impact_data['Mean'].max() * 100
bars = axes[1, 1].barh(impact_data['Metric'], impact_data['Normalized'],
                        color=sns.color_palette("Set2", 6), edgecolor='black', linewidth=0.5)
axes[1, 1].set_xlabel('Relative Impact (%)')
axes[1, 1].set_title('(d) Normalized Impact Distribution')
for bar, val in zip(bars, impact_data['Mean']):
    axes[1, 1].text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2.,
                   f'{val:.1f}', ha='left', va='center', fontsize=9)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig6_citizen_science.png', bbox_inches='tight')
plt.close()
print("  → Saved fig6_citizen_science.png")


# =============================================================================
# Module 6: Life Sciences Open Data Case Study
# =============================================================================
print("\n" + "=" * 70)
print("Module 6: Life Sciences Open Data Case Study")
print("=" * 70)

databases = ['GEO/ArrayExpress', 'PDB', 'GenBank', 'UniProt', 'PRIDE']
db_submissions = {
    'GEO/ArrayExpress': [45000, 52000, 61000, 72000, 85000, 98000, 110000, 125000, 140000, 155000],
    'PDB': [12000, 13500, 15000, 16800, 18500, 21000, 24000, 27000, 30000, 33000],
    'GenBank': [200000, 220000, 250000, 290000, 340000, 400000, 470000, 550000, 640000, 730000],
    'UniProt': [180000, 190000, 205000, 220000, 240000, 260000, 285000, 310000, 340000, 370000],
    'PRIDE': [5000, 7000, 9500, 12500, 16000, 20000, 25000, 31000, 38000, 46000],
}
years_range = list(range(2016, 2026))

# Citation impact of deposited data
db_citation_multiplier = {
    'GEO/ArrayExpress': 2.3, 'PDB': 3.1, 'GenBank': 1.8, 'UniProt': 2.5, 'PRIDE': 2.0
}
db_reuse_rate = {
    'GEO/ArrayExpress': 0.58, 'PDB': 0.72, 'GenBank': 0.45, 'UniProt': 0.65, 'PRIDE': 0.35
}

print("  Database Growth (2016-2025):")
for db in databases:
    growth = (db_submissions[db][-1] - db_submissions[db][0]) / db_submissions[db][0] * 100
    print(f"    {db}: {db_submissions[db][0]:,} → {db_submissions[db][-1]:,} ({growth:.0f}% growth)")

results['life_sciences'] = {
    db: {
        'growth_pct': round((db_submissions[db][-1] - db_submissions[db][0]) / db_submissions[db][0] * 100, 0),
        'citation_multiplier': db_citation_multiplier[db],
        'reuse_rate': db_reuse_rate[db],
    } for db in databases
}

# Figure 7: Life sciences case study
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (a) Database growth
colors_db = sns.color_palette("Set2", len(databases))
for idx, db in enumerate(databases):
    axes[0, 0].plot(years_range, np.array(db_submissions[db]) / 1000, 'o-',
                   label=db, color=colors_db[idx], linewidth=2)
axes[0, 0].set_xlabel('Year')
axes[0, 0].set_ylabel('Submissions (thousands)')
axes[0, 0].set_title('(a) Open Data Repository Growth')
axes[0, 0].legend(fontsize=8)
axes[0, 0].grid(True, alpha=0.3)
axes[0, 0].set_yscale('log')

# (b) Citation multiplier
mult_vals = list(db_citation_multiplier.values())
bars = axes[0, 1].bar(databases, mult_vals, color=colors_db, edgecolor='black', linewidth=0.5)
axes[0, 1].axhline(y=1.0, color='red', linestyle='--', alpha=0.7, label='Baseline (no data sharing)')
axes[0, 1].set_ylabel('Citation Multiplier')
axes[0, 1].set_title('(b) Citation Multiplier by Database')
axes[0, 1].legend()
for bar, val in zip(bars, mult_vals):
    axes[0, 1].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.05,
                   f'{val:.1f}x', ha='center', va='bottom', fontsize=10)
axes[0, 1].tick_params(axis='x', rotation=15)

# (c) Reuse rates
reuse_vals = list(db_reuse_rate.values())
bars = axes[1, 0].barh(databases, [r*100 for r in reuse_vals], color=colors_db, edgecolor='black', linewidth=0.5)
axes[1, 0].set_xlabel('Reuse Rate (%)')
axes[1, 0].set_title('(c) Dataset Reuse Rates')
for bar, val in zip(bars, reuse_vals):
    axes[1, 0].text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2.,
                   f'{val*100:.0f}%', ha='left', va='center', fontsize=10)

# (d) Compound impact heatmap
impact_matrix = np.array([
    [db_citation_multiplier[db] * db_reuse_rate[db] for db in databases]
    for _ in ['Citation Impact', 'Reuse Impact', 'Combined']
])
impact_matrix[0] = [db_citation_multiplier[db] for db in databases]
impact_matrix[1] = [db_reuse_rate[db] for db in databases]
impact_matrix[2] = [db_citation_multiplier[db] * db_reuse_rate[db] for db in databases]

sns.heatmap(pd.DataFrame(impact_matrix, index=['Citation\nMultiplier', 'Reuse\nRate', 'Combined\nImpact'],
                          columns=databases),
            annot=True, fmt='.2f', cmap='YlOrRd', ax=axes[1, 1])
axes[1, 1].set_title('(d) Combined Impact Assessment')
axes[1, 1].tick_params(axis='x', rotation=15)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig7_life_sciences.png', bbox_inches='tight')
plt.close()
print("  → Saved fig7_life_sciences.png")


# =============================================================================
# Summary Figure: Framework Overview
# =============================================================================
print("\n" + "=" * 70)
print("Generating Summary Figure")
print("=" * 70)

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Module summaries
modules = [
    ('OA Citation\nAdvantage', f"+{results['oaca']['pct_advantage']}%\ncitations", '#3498db'),
    ('Data Sharing\n& Reuse', f"Meta: +{results['data_sharing']['metadata_effect']}\nCode: +{results['data_sharing']['code_effect']}", '#2ecc71'),
    ('Preprint\nImpact', f"-{results['preprint']['reduction_pct']}%\nreview time", '#e74c3c'),
    ('FAIR\nCompliance', f"Mean: {results['fair']['mean_total']}\nF:{results['fair']['mean_F']} A:{results['fair']['mean_A']}", '#f39c12'),
    ('Citizen Science\nOutreach', f"ρ(part,alt)={results['citizen_science']['corr_participants_altmetric']}", '#9b59b6'),
    ('Life Sciences\nCase Study', f"PDB: {results['life_sciences']['PDB']['citation_multiplier']}x\nGEO: {results['life_sciences']['GEO/ArrayExpress']['citation_multiplier']}x", '#1abc9c'),
]

for idx, (title, stat, color) in enumerate(modules):
    row, col = divmod(idx, 3)
    ax = axes[row, col]
    ax.add_patch(plt.Rectangle((0.1, 0.1), 0.8, 0.8, facecolor=color, alpha=0.2,
                                edgecolor=color, linewidth=3, transform=ax.transAxes))
    ax.text(0.5, 0.7, title, ha='center', va='center', fontsize=14, fontweight='bold',
            transform=ax.transAxes)
    ax.text(0.5, 0.35, stat, ha='center', va='center', fontsize=12,
            transform=ax.transAxes, family='monospace')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

fig.suptitle('Open Access / Open Data Impact Analysis Framework — Key Results',
             fontsize=16, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(FIGURES_DIR / 'fig8_framework_summary.png', bbox_inches='tight')
plt.close()
print("  → Saved fig8_framework_summary.png")


# =============================================================================
# Pipeline Architecture Diagram
# =============================================================================
fig, ax = plt.subplots(figsize=(16, 8))
ax.axis('off')

# Pipeline stages
stages = [
    ('Data Collection', ['OpenAlex API', 'Crossref API', 'PubMed API', 'Altmetric API', 'Unpaywall API'], '#3498db'),
    ('Preprocessing', ['Deduplication', 'Field Normalization', 'OA Status Tagging', 'FAIR Scoring'], '#2ecc71'),
    ('Analysis', ['PSM + DiD', 'Regression', 'Survival Analysis', 'NLP/Topic Modeling'], '#e74c3c'),
    ('Evaluation', ['OACA Estimation', 'Reuse Metrics', 'Review Efficiency', 'Outreach Scoring'], '#f39c12'),
    ('Output', ['Dashboards', 'Reports', 'Policy Briefs', 'API Endpoints'], '#9b59b6'),
]

for idx, (stage_name, components, color) in enumerate(stages):
    x = 0.05 + idx * 0.19
    ax.add_patch(plt.Rectangle((x, 0.7), 0.17, 0.2, facecolor=color, alpha=0.3,
                                edgecolor=color, linewidth=2, transform=ax.transAxes))
    ax.text(x + 0.085, 0.85, stage_name, ha='center', va='center', fontsize=11,
            fontweight='bold', transform=ax.transAxes)

    for j, comp in enumerate(components):
        y = 0.55 - j * 0.12
        ax.add_patch(plt.Rectangle((x + 0.01, y), 0.15, 0.1, facecolor=color, alpha=0.1,
                                    edgecolor=color, linewidth=1, transform=ax.transAxes))
        ax.text(x + 0.085, y + 0.05, comp, ha='center', va='center', fontsize=8,
                transform=ax.transAxes)

    if idx < len(stages) - 1:
        ax.annotate('', xy=(x + 0.19, 0.8), xytext=(x + 0.17, 0.8),
                    xycoords='axes fraction', textcoords='axes fraction',
                    arrowprops=dict(arrowstyle='->', color='gray', lw=2))

ax.set_title('Bibliometrics/Altmetrics Analysis Pipeline Architecture', fontsize=14, fontweight='bold', pad=20)
plt.savefig(FIGURES_DIR / 'fig9_pipeline_architecture.png', bbox_inches='tight')
plt.close()
print("  → Saved fig9_pipeline_architecture.png")


# Save all results
with open(Path(__file__).parent.parent / 'results.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("\n✅ All experiments completed. Results saved to results.json")
print(f"✅ {len(list(FIGURES_DIR.glob('*.png')))} figures saved to figures/")
