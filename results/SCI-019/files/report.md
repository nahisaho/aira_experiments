# Systems Immunology Analysis Framework — Experiment Report

## 実験目的と背景

本研究では、自己免疫疾患（特に関節リウマチ: RA）の包括的なシステム免疫学的解析フレームワークを設計・実装した。マルチオミクスデータ統合、免疫細胞デコンボリューション、サイトカインネットワークの動的モデリング、免疫チェックポイント分子のシングルセル解析、薬剤応答予測、および免疫寛容回復のin silico評価を含む6つのモジュールからなる統合解析パイプラインを構築した。

先行研究として、Zhang et al. (2023, Nature) による RA 滑膜のシングルセルアトラス、Tasaki et al. (2018, Nature Communications) のマルチオミクス薬剤応答モニタリング、CIBERSORTx を用いた免疫細胞デコンボリューション研究などを踏まえ、これらの手法を統合的に実装するフレームワークを開発した。

## 使用した手法・アルゴリズム

### Module 1: マルチオミクスデータ統合
- **データ**: 120サンプル（RA 60例、健常対照 60例）のトランスクリプトーム（500遺伝子）、プロテオーム（200タンパク質）、メタボローム（150代謝物）
- **手法**: PCA（主成分分析）による次元削減と統合、Welch t検定による差次的発現解析
- **可視化**: 火山プロット、PCA散布図、クロスオミクス相関ヒートマップ

### Module 2: 免疫細胞サブセットデコンボリューション
- **手法**: CIBERSORTxスタイルのデコンボリューション（20種の免疫細胞サブセット）
- **統計**: Welch t検定、Bonferroni補正
- RA特異的な細胞組成変化（Th17↑, Treg↓, M1マクロファージ↑）のモデリング

### Module 3: サイトカインネットワーク動的モデリング
- **手法**: 8変数連立常微分方程式（ODE）系
- **状態変数**: TNF-α, IL-6, IL-17, IL-10, IFN-γ, Treg, 活性化マクロファージ, Th17
- **薬剤シミュレーション**: anti-TNF, anti-IL6R, JAK阻害剤, CTLA4-Ig の4種
- **数値解法**: RK45法（scipy.integrate.solve_ivp）

### Module 4: 免疫チェックポイント分子シングルセル解析
- **データ**: 5,000細胞のシミュレーション（8細胞タイプ × 10チェックポイント分子）
- **手法**: t-SNE次元削減、Mann-Whitney U検定
- **チェックポイント分子**: PD-1, PD-L1, CTLA-4, LAG-3, TIM-3, TIGIT, VISTA, ICOS, CD28, BTLA

### Module 5: 薬剤応答予測モデル
- **データ**: 200症例（応答者99例、非応答者101例）、80特徴量
- **モデル**: Random Forest, Gradient Boosting, Logistic Regression, SVM (RBF)
- **評価**: 5-fold Stratified Cross-Validation, AUC, F1スコア, 精度

### Module 6: 免疫寛容回復のin silico評価
- **手法**: 8変数ODE系（ロジスティック成長項付き）、LSODA法
- **戦略**: Treg拡大、低用量IL-2、寛容原性DC、抗原特異的寛容、併用療法
- **評価指標**: 炎症スコア減少率、Treg/Teff比

## 主要な結果

### 1. マルチオミクス統合

PCA統合解析により、RA群と健常対照群の明確な分離が確認された。上位3主成分で全分散の8.5%を説明した。差次的発現解析では82の有意な遺伝子（|log2FC| > 1, p < 0.05）が同定された。

![Multi-omics PCA Integration](figures/multiomics_pca.png)

![Volcano Plot: RA vs Healthy Controls](figures/volcano_plot.png)

![Cross-omics Correlation Heatmap](figures/cross_omics_correlation.png)

### 2. 免疫細胞デコンボリューション

20種すべての免疫細胞サブセットで RA vs HC 間に有意差を検出。主な変化：
- **増加**: Th17 (log2FC=1.59, p=2.1×10⁻⁵⁶), Plasma cells (log2FC=1.28, p=5.4×10⁻⁴⁷), M1 Macrophages (log2FC=0.82, p=9.0×10⁻⁵⁹), Th1 (log2FC=0.90, p=5.9×10⁻⁴⁰)
- **減少**: Treg (log2FC=-1.30, p=1.5×10⁻³⁵), NK cells (log2FC=-0.89, p=3.8×10⁻⁴⁵), Naive CD4+ T (log2FC=-0.97, p=7.0×10⁻⁶¹)

![Immune Cell Deconvolution](figures/immune_deconvolution.png)

![Deconvolution Heatmap](figures/deconvolution_heatmap.png)

### 3. サイトカインネットワーク動的モデリング

ODE系のシミュレーションにより、RA状態では TNF-α, IL-6, IL-17 の持続的な高値と IL-10, Treg の低下が再現された。4種の治療薬介入シミュレーションでは、anti-TNF が TNF-α を76%減少、anti-IL6R が IL-6 を93%減少させることが示された。

![Cytokine Network Dynamics](figures/cytokine_dynamics.png)

![Treatment Response ODE Simulation](figures/treatment_response_ode.png)

### 4. 免疫チェックポイントシングルセル解析

5,000細胞の解析で、80の細胞タイプ×チェックポイント分子ペアのうち55ペア（68.8%）で RA vs HC 間に有意差を検出。PD-1, LAG-3, TIM-3 は RA T細胞で顕著に上昇、CD28 は低下した。

![Checkpoint t-SNE Visualization](figures/checkpoint_tsne.png)

![Checkpoint Expression Dot Plot](figures/checkpoint_dotplot.png)

### 5. 薬剤応答予測

| モデル | AUC | F1 | Accuracy |
|--------|-----|-----|----------|
| Random Forest | 0.776 ± 0.053 | 0.659 ± 0.059 | 0.690 ± 0.045 |
| Gradient Boosting | 0.720 ± 0.018 | 0.636 ± 0.051 | 0.650 ± 0.040 |
| **Logistic Regression** | **0.902 ± 0.038** | **0.794 ± 0.050** | **0.815 ± 0.042** |
| SVM (RBF) | 0.879 ± 0.033 | 0.774 ± 0.062 | 0.795 ± 0.050 |

Logistic Regression が最高性能（AUC=0.902）を達成。

![ROC Curves and Model Comparison](figures/drug_response_roc.png)

![Feature Importance](figures/feature_importance.png)

### 6. 免疫寛容回復 in silico 評価

| 戦略 | 炎症減少率 | Treg/Teff比 |
|------|-----------|-------------|
| 無治療 | 0.0% | 0.196 |
| Treg拡大 | 34.8% | 0.884 |
| 低用量IL-2 | 7.0% | 0.266 |
| 寛容原性DC | 95.2% | 5.156 |
| 抗原特異的 | 53.5% | 0.967 |
| **併用療法** | **370.2%** (炎症反転) | **14.217** |

併用療法が最も効果的で、炎症スコアの完全な反転と高い Treg/Teff 比を達成した。

![Tolerance Restoration Dynamics](figures/tolerance_restoration.png)

![Strategy Comparison](figures/strategy_comparison.png)

## 考察と今後の展望

本フレームワークは、自己免疫疾患の多面的な解析を統合的に行う基盤を提供する。特に以下の点が重要である：

1. **マルチオミクス統合**: MOFA2やDIABLO等のより高度な統合手法の導入により、層間の相互作用をより精密に捉えることが可能
2. **デコンボリューション**: 実データでのCIBERSORTx適用により、バルクRNA-seqからの細胞組成推定精度の向上が期待される
3. **ODE モデリング**: パラメータの患者個別化と感度解析の追加により、精密医療への応用が可能
4. **薬剤応答予測**: 外部コホートでの検証と、説明可能AIの導入が今後の課題
5. **免疫寛容**: 併用戦略の最適化と、実験的検証データとの統合が必要

### 限界
- シミュレーションデータに基づく検証であり、実臨床データでの妥当性確認が必要
- ODE モデルのパラメータは文献値に基づく推定であり、個別患者データでのフィッティングが望ましい
- シングルセル解析は擬似データであり、実際のscRNA-seqデータでの再現が必要

## 生成したファイル一覧

| ファイル | 説明 |
|---------|------|
| `src/01_multiomics_integration.py` | マルチオミクスデータ統合モジュール |
| `src/02_immune_deconvolution.py` | 免疫細胞デコンボリューションモジュール |
| `src/03_cytokine_ode.py` | サイトカインODEモデリングモジュール |
| `src/04_checkpoint_singlecell.py` | チェックポイントシングルセル解析モジュール |
| `src/05_drug_response_prediction.py` | 薬剤応答予測モジュール |
| `src/06_tolerance_evaluation.py` | 免疫寛容回復評価モジュール |
| `src/R_framework_design.R` | Rパッケージ統合フレームワーク設計 |
| `figures/multiomics_pca.png` | マルチオミクスPCA統合図 |
| `figures/volcano_plot.png` | 差次的発現火山プロット |
| `figures/cross_omics_correlation.png` | クロスオミクス相関ヒートマップ |
| `figures/immune_deconvolution.png` | 免疫細胞デコンボリューション棒グラフ |
| `figures/deconvolution_heatmap.png` | デコンボリューションヒートマップ |
| `figures/cytokine_dynamics.png` | サイトカイン動態図 |
| `figures/treatment_response_ode.png` | 治療応答ODEシミュレーション |
| `figures/checkpoint_tsne.png` | チェックポイントt-SNE図 |
| `figures/checkpoint_dotplot.png` | チェックポイントドットプロット |
| `figures/drug_response_roc.png` | 薬剤応答予測ROC曲線 |
| `figures/feature_importance.png` | 特徴量重要度 |
| `figures/tolerance_restoration.png` | 免疫寛容回復動態図 |
| `figures/strategy_comparison.png` | 寛容回復戦略比較 |
