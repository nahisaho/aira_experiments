#!/usr/bin/env python3
"""
Cancer Proteogenomics Integrated Analysis Pipeline
===================================================
Simulates a CPTAC-style pancreatic ductal adenocarcinoma (PDAC) proteogenomics
analysis including:
  1. Variant peptide search
  2. RNA-protein expression discordance (translational regulation)
  3. Phosphoproteomics and kinase activity estimation (KSEA)
  4. Neoantigen candidate proteomics verification
  5. Multi-omics factor analysis (MOFA+) patient stratification
  6. CPTAC PDAC case study integration
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from sklearn.decomposition import NMF, PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import warnings
warnings.filterwarnings("ignore")

np.random.seed(42)
sns.set_theme(style="whitegrid", font_scale=1.1)

N_PATIENTS = 140
N_GENES = 5000
N_PHOSPHOSITES = 1200
N_KINASES = 45
SUBTYPES = ["Classical", "Basal-like", "Immunogenic"]

# ============================================================
# Helper: simulate CPTAC PDAC-like multi-omics data
# ============================================================
def simulate_cptac_pdac():
    """Generate synthetic CPTAC PDAC multi-omics data."""
    patients = [f"PDAC_{i:03d}" for i in range(N_PATIENTS)]
    genes = [f"Gene_{i}" for i in range(N_GENES)]

    # Assign subtypes with imbalanced proportions
    subtype_labels = np.random.choice(SUBTYPES, N_PATIENTS, p=[0.45, 0.35, 0.20])

    # RNA-seq (log2 TPM)
    rna = np.random.randn(N_PATIENTS, N_GENES) * 1.5 + 6
    # Add subtype-specific signals in first 300 genes
    for i, st in enumerate(SUBTYPES):
        mask = subtype_labels == st
        rna[mask, i*100:(i+1)*100] += np.random.uniform(1.5, 3.0, size=(mask.sum(), 100))

    # Proteomics (log2 intensity) – correlated with RNA but with noise
    protein = rna * 0.6 + np.random.randn(N_PATIENTS, N_GENES) * 1.2 + 2
    # Add translational regulation effects for specific gene sets
    translation_regulated = list(range(400, 600))
    protein[:, translation_regulated] += np.random.randn(N_PATIENTS, 200) * 2.0

    # Genomic variants (binary: 0/1 for each gene)
    mutation_rate = np.full(N_GENES, 0.03)
    # KRAS, TP53, SMAD4, CDKN2A high mutation rates for PDAC
    driver_genes = [0, 1, 2, 3]  # Gene_0=KRAS, Gene_1=TP53, etc.
    driver_names = {0: "KRAS", 1: "TP53", 2: "SMAD4", 3: "CDKN2A"}
    mutation_rate[0] = 0.92  # KRAS
    mutation_rate[1] = 0.72  # TP53
    mutation_rate[2] = 0.31  # SMAD4
    mutation_rate[3] = 0.25  # CDKN2A
    mutations = (np.random.rand(N_PATIENTS, N_GENES) < mutation_rate).astype(int)

    # Phosphoproteomics (log2 intensity)
    phospho_sites = [f"pSite_{i}" for i in range(N_PHOSPHOSITES)]
    phospho = np.random.randn(N_PATIENTS, N_PHOSPHOSITES) * 1.8 + 5
    # Subtype-specific phosphorylation patterns
    for i, st in enumerate(SUBTYPES):
        mask = subtype_labels == st
        phospho[mask, i*100:(i+1)*100] += np.random.uniform(1.0, 2.5, size=(mask.sum(), 100))

    rna_df = pd.DataFrame(rna, index=patients, columns=genes)
    protein_df = pd.DataFrame(protein, index=patients, columns=genes)
    mutation_df = pd.DataFrame(mutations, index=patients, columns=genes)
    phospho_df = pd.DataFrame(phospho, index=patients, columns=phospho_sites)

    clinical = pd.DataFrame({
        "Patient": patients,
        "Subtype": subtype_labels,
        "Stage": np.random.choice(["I", "II", "III", "IV"], N_PATIENTS, p=[0.1, 0.35, 0.35, 0.2]),
        "OS_months": np.random.exponential(18, N_PATIENTS).clip(1, 60),
        "Age": np.random.normal(65, 10, N_PATIENTS).astype(int).clip(35, 90),
    })
    clinical["OS_event"] = (clinical["OS_months"] < 24).astype(int)

    return rna_df, protein_df, mutation_df, phospho_df, clinical, driver_names


# ============================================================
# Module 1: Variant Peptide Search
# ============================================================
def variant_peptide_search(mutation_df, protein_df, driver_names):
    """Simulate variant peptide identification from genomic variants."""
    print("=" * 60)
    print("Module 1: Variant Peptide Search")
    print("=" * 60)

    # For each mutated gene, simulate variant peptide detection
    results = []
    all_genes_with_mutations = mutation_df.columns[mutation_df.sum(axis=0) > 0]

    for gene_idx in range(min(50, len(all_genes_with_mutations))):
        gene = all_genes_with_mutations[gene_idx]
        n_mutated = mutation_df[gene].sum()
        # Detection probability depends on protein abundance
        mean_abundance = protein_df[gene].mean()
        detect_prob = min(0.9, max(0.05, (mean_abundance - 2) / 10))
        detected = np.random.rand() < detect_prob
        peptide_score = np.random.uniform(20, 80) if detected else 0
        fdr = np.random.uniform(0.001, 0.05) if detected else 1.0

        gene_name = driver_names.get(int(gene.split("_")[1]), gene)
        results.append({
            "Gene": gene_name,
            "Mutated_Samples": n_mutated,
            "Variant_Peptide_Detected": detected,
            "Peptide_Score": round(peptide_score, 2),
            "FDR": round(fdr, 4),
            "Mean_Protein_Abundance": round(mean_abundance, 2),
        })

    vp_df = pd.DataFrame(results)
    detected = vp_df[vp_df["Variant_Peptide_Detected"]]
    print(f"  Total genes screened: {len(vp_df)}")
    print(f"  Variant peptides detected: {len(detected)} ({len(detected)/len(vp_df)*100:.1f}%)")
    print(f"  Mean peptide score (detected): {detected['Peptide_Score'].mean():.2f}")

    # Figure 1: Variant peptide detection overview
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Bar chart of top detected variant peptides
    top_detected = detected.nlargest(15, "Peptide_Score")
    axes[0].barh(top_detected["Gene"], top_detected["Peptide_Score"],
                 color=sns.color_palette("viridis", len(top_detected)))
    axes[0].set_xlabel("Peptide Identification Score")
    axes[0].set_title("Top 15 Variant Peptides by Score")
    axes[0].invert_yaxis()

    # Scatter: abundance vs detection
    colors = ["#e74c3c" if d else "#95a5a6" for d in vp_df["Variant_Peptide_Detected"]]
    axes[1].scatter(vp_df["Mean_Protein_Abundance"], vp_df["Mutated_Samples"],
                    c=colors, alpha=0.7, edgecolors="k", linewidth=0.5)
    axes[1].set_xlabel("Mean Protein Abundance (log2)")
    axes[1].set_ylabel("Number of Mutated Samples")
    axes[1].set_title("Variant Peptide Detection vs Protein Abundance")
    axes[1].legend(handles=[
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#e74c3c', markersize=8, label='Detected'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#95a5a6', markersize=8, label='Not Detected'),
    ])

    plt.tight_layout()
    plt.savefig("figures/fig1_variant_peptide.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  -> Saved figures/fig1_variant_peptide.png")

    return vp_df


# ============================================================
# Module 2: RNA-Protein Expression Discordance
# ============================================================
def rna_protein_discordance(rna_df, protein_df, clinical):
    """Analyze mRNA-protein expression discordance to infer translational regulation."""
    print("\n" + "=" * 60)
    print("Module 2: RNA-Protein Expression Discordance Analysis")
    print("=" * 60)

    correlations = []
    for gene in rna_df.columns[:1000]:
        r, p = stats.spearmanr(rna_df[gene], protein_df[gene])
        correlations.append({"Gene": gene, "Spearman_r": r, "p_value": p})

    corr_df = pd.DataFrame(correlations)
    corr_df["Category"] = pd.cut(corr_df["Spearman_r"],
                                  bins=[-1, 0.2, 0.5, 1.0],
                                  labels=["Low (<0.2)", "Medium (0.2-0.5)", "High (>0.5)"])

    median_corr = corr_df["Spearman_r"].median()
    low_corr = (corr_df["Spearman_r"] < 0.2).sum()
    print(f"  Median mRNA-protein correlation: {median_corr:.3f}")
    print(f"  Genes with low correlation (<0.2): {low_corr} ({low_corr/len(corr_df)*100:.1f}%)")
    print(f"  Genes with high correlation (>0.5): {(corr_df['Spearman_r'] > 0.5).sum()}")

    # Translational efficiency estimation
    te_scores = protein_df.values[:, :1000] - rna_df.values[:, :1000] * 0.6
    te_df = pd.DataFrame(te_scores, index=rna_df.index, columns=rna_df.columns[:1000])

    # Subtype-specific TE
    subtype_te = {}
    for st in SUBTYPES:
        mask = clinical["Subtype"].values == st
        subtype_te[st] = te_df.iloc[mask].mean()

    # Figure 2: RNA-protein discordance
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Distribution of correlations
    axes[0, 0].hist(corr_df["Spearman_r"], bins=50, color="#3498db", edgecolor="black", alpha=0.7)
    axes[0, 0].axvline(median_corr, color="red", linestyle="--", label=f"Median = {median_corr:.3f}")
    axes[0, 0].set_xlabel("Spearman Correlation (mRNA vs Protein)")
    axes[0, 0].set_ylabel("Number of Genes")
    axes[0, 0].set_title("Distribution of mRNA-Protein Correlation")
    axes[0, 0].legend()

    # Category pie chart
    cat_counts = corr_df["Category"].value_counts()
    axes[0, 1].pie(cat_counts, labels=cat_counts.index, autopct='%1.1f%%',
                    colors=["#e74c3c", "#f39c12", "#27ae60"])
    axes[0, 1].set_title("Correlation Category Distribution")

    # Scatter plot example (Gene_450 = translationally regulated)
    gene_ex = "Gene_450"
    axes[1, 0].scatter(rna_df[gene_ex], protein_df[gene_ex], alpha=0.5,
                        c=clinical["Subtype"].map({"Classical": "#3498db", "Basal-like": "#e74c3c", "Immunogenic": "#27ae60"}))
    r_val, _ = stats.spearmanr(rna_df[gene_ex], protein_df[gene_ex])
    axes[1, 0].set_xlabel(f"mRNA Expression (log2 TPM)")
    axes[1, 0].set_ylabel(f"Protein Expression (log2 intensity)")
    axes[1, 0].set_title(f"{gene_ex}: mRNA vs Protein (r={r_val:.3f})")

    # Heatmap: subtype-specific TE for top variable genes
    te_var = te_df.var(axis=0).nlargest(30)
    te_sub_mat = pd.DataFrame({st: subtype_te[st][te_var.index] for st in SUBTYPES})
    sns.heatmap(te_sub_mat.T, cmap="RdBu_r", center=0, ax=axes[1, 1],
                xticklabels=False, yticklabels=True)
    axes[1, 1].set_title("Translational Efficiency by Subtype\n(Top 30 Variable Genes)")
    axes[1, 1].set_xlabel("Genes")

    plt.tight_layout()
    plt.savefig("figures/fig2_rna_protein_discordance.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  -> Saved figures/fig2_rna_protein_discordance.png")

    return corr_df, te_df


# ============================================================
# Module 3: Phosphoproteomics & Kinase Activity (KSEA)
# ============================================================
def kinase_activity_ksea(phospho_df, clinical):
    """Perform KSEA-like kinase activity estimation from phosphoproteomics."""
    print("\n" + "=" * 60)
    print("Module 3: Phosphoproteomics & Kinase Activity Estimation (KSEA)")
    print("=" * 60)

    kinase_names = [
        "AKT1", "AKT2", "BRAF", "CDK1", "CDK2", "CDK4", "CDK6",
        "EGFR", "ERK1", "ERK2", "FGFR1", "GSK3B", "JAK2", "JNK1",
        "MAPK14", "MEK1", "MET", "MTOR", "PAK1", "PDGFRA", "PI3K",
        "PKA", "PKC_alpha", "PKC_delta", "PLK1", "RAF1", "RET",
        "ROCK1", "RSK2", "SRC", "TGFBR1", "VEGFR2", "ABL1",
        "AURORA_A", "AURORA_B", "CAMK2", "CK2", "DYRK1A", "FGFR2",
        "IGF1R", "IKK", "IRAK4", "LCK", "NEK2", "WEE1"
    ]

    # Simulate kinase-substrate relationships
    n_substrates_per_kinase = np.random.randint(10, 50, N_KINASES)
    kinase_substrate_map = {}
    for i, kinase in enumerate(kinase_names):
        substrates = np.random.choice(N_PHOSPHOSITES, n_substrates_per_kinase[i], replace=False)
        kinase_substrate_map[kinase] = substrates

    # KSEA scoring
    ksea_scores = np.zeros((N_PATIENTS, N_KINASES))
    for j, kinase in enumerate(kinase_names):
        substrate_idx = kinase_substrate_map[kinase]
        substrate_data = phospho_df.iloc[:, substrate_idx]
        # z-score enrichment
        mean_substrate = substrate_data.mean(axis=1)
        global_mean = phospho_df.mean(axis=1)
        global_std = phospho_df.std(axis=1)
        ksea_scores[:, j] = (mean_substrate - global_mean) / global_std

    # Add subtype-specific kinase patterns
    subtype_kinase_effects = {
        "Basal-like": ["EGFR", "SRC", "BRAF", "ERK1", "ERK2", "MEK1"],
        "Classical": ["AKT1", "MTOR", "PI3K", "CDK4", "CDK6"],
        "Immunogenic": ["JAK2", "JNK1", "IKK", "IRAK4"],
    }
    for st, kinases in subtype_kinase_effects.items():
        mask = clinical["Subtype"].values == st
        for k in kinases:
            if k in kinase_names:
                idx = kinase_names.index(k)
                ksea_scores[mask, idx] += np.random.uniform(1.0, 2.5, mask.sum())

    ksea_df = pd.DataFrame(ksea_scores, index=phospho_df.index, columns=kinase_names)

    # Statistics
    for st in SUBTYPES:
        mask = clinical["Subtype"].values == st
        top_kinases = ksea_df.iloc[mask].mean().nlargest(5)
        print(f"  {st} - Top kinases: {', '.join(top_kinases.index.tolist())}")

    # Figure 3: KSEA kinase activity
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Heatmap of kinase activities by subtype
    ksea_subtype = pd.DataFrame({
        st: ksea_df.iloc[clinical["Subtype"].values == st].mean()
        for st in SUBTYPES
    })
    top_variable_kinases = ksea_df.var().nlargest(25).index
    sns.heatmap(ksea_subtype.loc[top_variable_kinases], cmap="RdBu_r", center=0,
                ax=axes[0], annot=True, fmt=".2f")
    axes[0].set_title("Mean Kinase Activity Score by Subtype\n(Top 25 Variable Kinases)")
    axes[0].set_xlabel("Subtype")

    # Box plot of key kinases
    key_kinases = ["EGFR", "MTOR", "JAK2", "SRC", "AKT1", "ERK1"]
    plot_data = []
    for k in key_kinases:
        for i, st in enumerate(SUBTYPES):
            mask = clinical["Subtype"].values == st
            for val in ksea_df.loc[mask, k]:
                plot_data.append({"Kinase": k, "Subtype": st, "Activity": val})
    plot_df = pd.DataFrame(plot_data)
    sns.boxplot(data=plot_df, x="Kinase", y="Activity", hue="Subtype", ax=axes[1],
                palette={"Classical": "#3498db", "Basal-like": "#e74c3c", "Immunogenic": "#27ae60"})
    axes[1].set_title("Kinase Activity Scores by Subtype")
    axes[1].set_ylabel("KSEA z-score")
    axes[1].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig("figures/fig3_kinase_activity.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  -> Saved figures/fig3_kinase_activity.png")

    return ksea_df


# ============================================================
# Module 4: Neoantigen Candidate Verification
# ============================================================
def neoantigen_verification(mutation_df, protein_df, rna_df, vp_df):
    """Verify neoantigen candidates using proteomics data."""
    print("\n" + "=" * 60)
    print("Module 4: Neoantigen Candidate Proteomics Verification")
    print("=" * 60)

    # Simulate HLA binding predictions
    hla_types = ["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02", "HLA-B*35:01",
                 "HLA-C*07:01", "HLA-C*04:01"]
    neoantigen_candidates = []

    detected_vps = vp_df[vp_df["Variant_Peptide_Detected"]].head(20)
    for _, row in detected_vps.iterrows():
        for hla in np.random.choice(hla_types, np.random.randint(1, 4), replace=False):
            binding_affinity = np.random.exponential(200)
            immunogenicity = np.random.uniform(0, 1)
            ms_validated = np.random.rand() < 0.35
            neoantigen_candidates.append({
                "Gene": row["Gene"],
                "HLA_Allele": hla,
                "Binding_Affinity_nM": round(binding_affinity, 1),
                "Immunogenicity_Score": round(immunogenicity, 3),
                "Strong_Binder": binding_affinity < 500,
                "MS_Validated": ms_validated,
                "Peptide_Score": row["Peptide_Score"],
            })

    neo_df = pd.DataFrame(neoantigen_candidates)
    strong_binders = neo_df[neo_df["Strong_Binder"]]
    validated = neo_df[neo_df["MS_Validated"]]

    print(f"  Total neoantigen candidates: {len(neo_df)}")
    print(f"  Strong binders (<500 nM): {len(strong_binders)} ({len(strong_binders)/len(neo_df)*100:.1f}%)")
    print(f"  MS-validated neoantigens: {len(validated)} ({len(validated)/len(neo_df)*100:.1f}%)")
    print(f"  Strong binders & MS-validated: {len(neo_df[(neo_df['Strong_Binder']) & (neo_df['MS_Validated'])])}")

    # Figure 4: Neoantigen verification
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Binding affinity distribution
    axes[0].hist(neo_df["Binding_Affinity_nM"], bins=30, color="#9b59b6",
                 edgecolor="black", alpha=0.7)
    axes[0].axvline(500, color="red", linestyle="--", label="Strong binding threshold")
    axes[0].set_xlabel("Binding Affinity (nM)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("HLA Binding Affinity Distribution")
    axes[0].legend()

    # Immunogenicity vs binding affinity
    colors = ["#27ae60" if v else "#e74c3c" for v in neo_df["MS_Validated"]]
    axes[1].scatter(neo_df["Binding_Affinity_nM"], neo_df["Immunogenicity_Score"],
                    c=colors, alpha=0.6, edgecolors="k", linewidth=0.5)
    axes[1].set_xlabel("Binding Affinity (nM)")
    axes[1].set_ylabel("Immunogenicity Score")
    axes[1].set_title("Neoantigen Characterization")
    axes[1].legend(handles=[
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#27ae60', markersize=8, label='MS Validated'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#e74c3c', markersize=8, label='Not Validated'),
    ])

    # HLA allele distribution
    hla_counts = neo_df.groupby("HLA_Allele")["MS_Validated"].agg(["sum", "count"])
    hla_counts.columns = ["Validated", "Total"]
    hla_counts["Not_Validated"] = hla_counts["Total"] - hla_counts["Validated"]
    hla_counts[["Validated", "Not_Validated"]].plot(kind="bar", stacked=True, ax=axes[2],
                                                      color=["#27ae60", "#bdc3c7"])
    axes[2].set_title("Neoantigen Validation by HLA Allele")
    axes[2].set_ylabel("Count")
    axes[2].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig("figures/fig4_neoantigen.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  -> Saved figures/fig4_neoantigen.png")

    return neo_df


# ============================================================
# Module 5: Multi-Omics Factor Analysis (MOFA+)
# ============================================================
def mofa_analysis(rna_df, protein_df, phospho_df, clinical):
    """Perform MOFA+-like multi-omics factor decomposition."""
    print("\n" + "=" * 60)
    print("Module 5: Multi-Omics Factor Analysis (MOFA+)")
    print("=" * 60)

    # Use top variable features from each omics
    n_features = 500
    rna_var = rna_df.var(axis=0).nlargest(n_features).index
    prot_var = protein_df.var(axis=0).nlargest(n_features).index
    phospho_var = phospho_df.var(axis=0).nlargest(min(n_features, N_PHOSPHOSITES)).index

    # Concatenate and standardize
    combined = pd.concat([
        rna_df[rna_var].add_prefix("RNA_"),
        protein_df[prot_var].add_prefix("Prot_"),
        phospho_df[phospho_var].add_prefix("Phospho_"),
    ], axis=1)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(combined)

    # Factor decomposition via NMF-like approach (shifted data for NMF)
    X_shifted = X_scaled - X_scaled.min() + 0.01
    n_factors = 10
    nmf = NMF(n_components=n_factors, random_state=42, max_iter=500)
    W = nmf.fit_transform(X_shifted)  # Patient x Factor
    H = nmf.components_                # Factor x Feature

    # Also do PCA for visualization
    pca = PCA(n_components=n_factors, random_state=42)
    W_pca = pca.fit_transform(X_scaled)
    explained_var = pca.explained_variance_ratio_

    # Cluster patients using factor scores
    kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(W_pca[:, :5])
    sil_score = silhouette_score(W_pca[:, :5], clusters)

    factor_df = pd.DataFrame(W_pca[:, :n_factors],
                              index=rna_df.index,
                              columns=[f"Factor_{i+1}" for i in range(n_factors)])
    factor_df["Cluster"] = clusters
    factor_df["Subtype"] = clinical["Subtype"].values

    print(f"  Variance explained (first 5 factors): {explained_var[:5].sum()*100:.1f}%")
    print(f"  Silhouette score (3 clusters): {sil_score:.3f}")

    # Cross-tabulate clusters vs subtypes
    ct = pd.crosstab(factor_df["Cluster"], factor_df["Subtype"])
    print(f"  Cluster-Subtype cross-tabulation:\n{ct}")

    # Variance explained per omics per factor
    n_rna_feat = len(rna_var)
    n_prot_feat = len(prot_var)
    n_phospho_feat = len(phospho_var)
    omics_var = np.zeros((n_factors, 3))
    for f in range(n_factors):
        loadings = pca.components_[f]
        omics_var[f, 0] = np.sum(loadings[:n_rna_feat] ** 2)
        omics_var[f, 1] = np.sum(loadings[n_rna_feat:n_rna_feat+n_prot_feat] ** 2)
        omics_var[f, 2] = np.sum(loadings[n_rna_feat+n_prot_feat:] ** 2)
    omics_var = omics_var / omics_var.sum(axis=1, keepdims=True)

    # Figure 5: MOFA+ analysis
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # PCA scatter colored by subtype
    subtype_colors = {"Classical": "#3498db", "Basal-like": "#e74c3c", "Immunogenic": "#27ae60"}
    for st in SUBTYPES:
        mask = factor_df["Subtype"] == st
        axes[0, 0].scatter(W_pca[mask, 0], W_pca[mask, 1],
                           c=subtype_colors[st], label=st, alpha=0.6, edgecolors="k", linewidth=0.3)
    axes[0, 0].set_xlabel(f"Factor 1 ({explained_var[0]*100:.1f}%)")
    axes[0, 0].set_ylabel(f"Factor 2 ({explained_var[1]*100:.1f}%)")
    axes[0, 0].set_title("Multi-Omics Factor Analysis: Patient Space")
    axes[0, 0].legend()

    # Variance explained bar chart
    axes[0, 1].bar(range(1, n_factors+1), explained_var * 100, color="#3498db", edgecolor="black")
    axes[0, 1].set_xlabel("Factor")
    axes[0, 1].set_ylabel("Variance Explained (%)")
    axes[0, 1].set_title("Variance Explained per Factor")

    # Omics contribution per factor
    x = np.arange(n_factors)
    w = 0.6
    bottom_rna = np.zeros(n_factors)
    bottom_prot = omics_var[:, 0]
    bottom_phospho = omics_var[:, 0] + omics_var[:, 1]
    axes[1, 0].bar(x, omics_var[:, 0], w, label="Transcriptomics", color="#3498db")
    axes[1, 0].bar(x, omics_var[:, 1], w, bottom=bottom_prot, label="Proteomics", color="#e74c3c")
    axes[1, 0].bar(x, omics_var[:, 2], w, bottom=bottom_phospho, label="Phosphoproteomics", color="#27ae60")
    axes[1, 0].set_xlabel("Factor")
    axes[1, 0].set_ylabel("Relative Contribution")
    axes[1, 0].set_title("Omics Contribution per Factor")
    axes[1, 0].legend()
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels([f"F{i+1}" for i in range(n_factors)])

    # Survival by cluster
    for c in sorted(factor_df["Cluster"].unique()):
        mask = (factor_df["Cluster"] == c).values
        os_vals = clinical.loc[mask, "OS_months"].sort_values()
        survival = np.arange(1, len(os_vals) + 1)[::-1] / len(os_vals)
        axes[1, 1].step(os_vals, survival, label=f"Cluster {c} (n={mask.sum()})", linewidth=2)
    axes[1, 1].set_xlabel("Overall Survival (months)")
    axes[1, 1].set_ylabel("Survival Probability")
    axes[1, 1].set_title("Kaplan-Meier by MOFA+ Cluster")
    axes[1, 1].legend()
    axes[1, 1].set_xlim(0, 60)

    plt.tight_layout()
    plt.savefig("figures/fig5_mofa_analysis.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  -> Saved figures/fig5_mofa_analysis.png")

    return factor_df, explained_var, sil_score


# ============================================================
# Module 6: CPTAC PDAC Case Study Integration
# ============================================================
def cptac_case_study(rna_df, protein_df, mutation_df, phospho_df, clinical,
                     ksea_df, corr_df, factor_df):
    """Integrated CPTAC PDAC case study analysis."""
    print("\n" + "=" * 60)
    print("Module 6: CPTAC PDAC Integrated Case Study")
    print("=" * 60)

    # Multi-omics integration summary
    # Subtype-specific differential analysis
    results = {}
    for st in SUBTYPES:
        mask = clinical["Subtype"].values == st
        rest = ~mask

        # Differential protein expression
        diff_prot = []
        for gene in protein_df.columns[:500]:
            t_stat, p_val = stats.ttest_ind(protein_df.loc[mask, gene], protein_df.loc[rest, gene])
            fc = protein_df.loc[mask, gene].mean() - protein_df.loc[rest, gene].mean()
            diff_prot.append({"Gene": gene, "log2FC": fc, "p_value": p_val, "t_stat": t_stat})

        diff_df = pd.DataFrame(diff_prot)
        diff_df["padj"] = diff_df["p_value"] * len(diff_df)  # Bonferroni
        diff_df["padj"] = diff_df["padj"].clip(upper=1.0)
        sig_up = ((diff_df["log2FC"] > 0.5) & (diff_df["padj"] < 0.05)).sum()
        sig_down = ((diff_df["log2FC"] < -0.5) & (diff_df["padj"] < 0.05)).sum()
        results[st] = {"sig_up": sig_up, "sig_down": sig_down, "diff_df": diff_df}
        print(f"  {st}: {sig_up} up-regulated, {sig_down} down-regulated proteins (|log2FC|>0.5, padj<0.05)")

    # Figure 6: Integrated case study
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # Volcano plot for Basal-like subtype
    diff = results["Basal-like"]["diff_df"]
    diff["-log10p"] = -np.log10(diff["padj"].clip(lower=1e-50))
    colors = np.where((diff["log2FC"].abs() > 0.5) & (diff["padj"] < 0.05),
                       np.where(diff["log2FC"] > 0, "#e74c3c", "#3498db"), "#cccccc")
    axes[0, 0].scatter(diff["log2FC"], diff["-log10p"], c=colors, alpha=0.5, s=15)
    axes[0, 0].axhline(-np.log10(0.05), color="gray", linestyle="--", alpha=0.5)
    axes[0, 0].axvline(0.5, color="gray", linestyle="--", alpha=0.5)
    axes[0, 0].axvline(-0.5, color="gray", linestyle="--", alpha=0.5)
    axes[0, 0].set_xlabel("log2 Fold Change")
    axes[0, 0].set_ylabel("-log10(adjusted p-value)")
    axes[0, 0].set_title("Volcano Plot: Basal-like vs Others")

    # Mutation landscape
    driver_genes = ["Gene_0", "Gene_1", "Gene_2", "Gene_3"]
    driver_labels = ["KRAS", "TP53", "SMAD4", "CDKN2A"]
    mut_freq = mutation_df[driver_genes].mean() * 100
    bar_colors = ["#e74c3c", "#3498db", "#27ae60", "#f39c12"]
    axes[0, 1].bar(driver_labels, mut_freq, color=bar_colors, edgecolor="black")
    axes[0, 1].set_ylabel("Mutation Frequency (%)")
    axes[0, 1].set_title("Driver Gene Mutation Frequency (PDAC)")
    for i, v in enumerate(mut_freq):
        axes[0, 1].text(i, v + 1, f"{v:.1f}%", ha="center", fontweight="bold")

    # Multi-omics integration heatmap (top features per subtype)
    # Select top variable genes from protein data
    top_genes = protein_df.var().nlargest(50).index
    sorted_idx = np.argsort(clinical["Subtype"].values)
    mat = protein_df.iloc[sorted_idx][top_genes].values
    mat_z = (mat - mat.mean(axis=0)) / mat.std(axis=0)
    im = axes[1, 0].imshow(mat_z.T, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    axes[1, 0].set_xlabel("Patients (sorted by subtype)")
    axes[1, 0].set_ylabel("Top 50 Variable Proteins")
    axes[1, 0].set_title("Protein Expression Heatmap")
    plt.colorbar(im, ax=axes[1, 0], label="z-score")

    # Subtype summary
    summary_data = {
        "Subtype": SUBTYPES,
        "N_patients": [sum(clinical["Subtype"] == st) for st in SUBTYPES],
        "Median_OS": [clinical.loc[clinical["Subtype"] == st, "OS_months"].median() for st in SUBTYPES],
        "Mean_Age": [clinical.loc[clinical["Subtype"] == st, "Age"].mean() for st in SUBTYPES],
    }
    summary_df = pd.DataFrame(summary_data)
    axes[1, 1].axis("off")
    table = axes[1, 1].table(
        cellText=summary_df.values,
        colLabels=summary_df.columns,
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.2, 1.8)
    axes[1, 1].set_title("Patient Subtype Summary", pad=20)

    plt.tight_layout()
    plt.savefig("figures/fig6_cptac_case_study.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  -> Saved figures/fig6_cptac_case_study.png")

    return results


# ============================================================
# Main Pipeline
# ============================================================
def main():
    print("=" * 60)
    print("  Cancer Proteogenomics Integrated Analysis Pipeline")
    print("  CPTAC Pancreatic Ductal Adenocarcinoma Case Study")
    print("=" * 60)
    print()

    # Simulate data
    print("Generating synthetic CPTAC PDAC data...")
    rna_df, protein_df, mutation_df, phospho_df, clinical, driver_names = simulate_cptac_pdac()
    print(f"  Patients: {N_PATIENTS}, Genes: {N_GENES}, Phosphosites: {N_PHOSPHOSITES}")
    print(f"  Subtypes: {dict(zip(*np.unique(clinical['Subtype'], return_counts=True)))}")
    print()

    # Module 1
    vp_df = variant_peptide_search(mutation_df, protein_df, driver_names)

    # Module 2
    corr_df, te_df = rna_protein_discordance(rna_df, protein_df, clinical)

    # Module 3
    ksea_df = kinase_activity_ksea(phospho_df, clinical)

    # Module 4
    neo_df = neoantigen_verification(mutation_df, protein_df, rna_df, vp_df)

    # Module 5
    factor_df, explained_var, sil_score = mofa_analysis(rna_df, protein_df, phospho_df, clinical)

    # Module 6
    results = cptac_case_study(rna_df, protein_df, mutation_df, phospho_df, clinical,
                                ksea_df, corr_df, factor_df)

    print("\n" + "=" * 60)
    print("  Pipeline Complete - All figures saved to figures/")
    print("=" * 60)

    # Return key metrics for report
    return {
        "n_patients": N_PATIENTS,
        "n_genes": N_GENES,
        "n_phosphosites": N_PHOSPHOSITES,
        "vp_detected": len(vp_df[vp_df["Variant_Peptide_Detected"]]),
        "vp_total": len(vp_df),
        "median_corr": corr_df["Spearman_r"].median(),
        "low_corr_pct": (corr_df["Spearman_r"] < 0.2).mean() * 100,
        "sil_score": sil_score,
        "explained_var_5": explained_var[:5].sum() * 100,
        "neo_total": len(neo_df),
        "neo_validated": len(neo_df[neo_df["MS_Validated"]]),
    }


if __name__ == "__main__":
    main()
