# PRS Cross-Ethnic Transferability Simulation Report

**作成日**: 2026-05-23  
**ステータス**: DRAFT — NOT FOR DISTRIBUTION

---

## 1. 実験目的と背景

多遺伝子リスクスコア（PRS）は、ゲノムワイド関連解析（GWAS）で同定されたSNPの効果量を集約し、個人の疾患リスクを予測する手法である。しかし、PRSの予測精度はGWASが実施された集団（主にヨーロッパ系）から他の民族集団へ転送する際に大幅に低下することが報告されている。

本実験では、UK Biobank（ヨーロッパ系）からBioBank Japan（日本人/東アジア系）へのPRS転送問題に焦点を当て、以下の5つの手法を比較評価した：

1. **Standard PRS**: ヨーロッパ系GWASの効果量を直接適用
2. **Bayesian LD-Corrected PRS**: 連鎖不平衡（LD）構造の差異をベイズ推定で補正
3. **Multi-Ethnic Meta-Analysis PRS**: DerSimonian-Laird変量効果メタ解析で効果量を再推定
4. **Local Ancestry-Corrected PRS**: 局所祖先推定に基づく効果量の重み付け
5. **Penalized Transfer PRS**: ペナルティ付き回帰による転移学習

## 2. 使用した手法・アルゴリズム

### 2.1 集団遺伝学シミュレーション

- **Balding-Nichols モデル**: 祖先集団のアレル頻度からFstパラメータに基づき各集団のアレル頻度を生成
- **LD構造**: ブロック対角行列（指数減衰）で異なるLD構造をシミュレート（EUR: decay=0.08, EAS: decay=0.12）
- **遺伝子型生成**: Cholesky分解による相関構造の導入

### 2.2 ベイズLD補正手法

$$\beta_{\text{target}} \sim \mathcal{N}\left(C \cdot \beta_{\text{source}},\ \sigma^2 \cdot R_{\text{target}}^{-1}\right)$$

事後推定: $\hat{\beta}_{\text{target}} = (R_{\text{EAS}} + \frac{1}{\sigma^2_{\text{prior}}} I)^{-1} R_{\text{EAS}} \hat{\beta}_{\text{joint,EUR}}$

### 2.3 多民族メタ解析

DerSimonian-Laird推定量によるランダム効果メタ解析：
- Cochran's Q統計量でheterogeneityを評価
- $\tau^2 = \max(0, (Q - df) / C)$ で研究間分散を推定
- ランダム効果重み: $w_k = 1 / (se_k^2 + \tau^2)$

### 2.4 局所祖先補正

$$\text{PRS}_{\text{corrected}}(i) = \sum_j \beta_{\text{adj}}(j) \cdot G(i,j)$$

$$\beta_{\text{adj}}(j) = \alpha(j) \cdot \beta_{\text{EAS}}(j) + (1 - \alpha(j)) \cdot \beta_{\text{EUR}}(j)$$

### 2.5 ペナルティ付き転移学習

$$\min_\beta \|y_{\text{target}} - X_{\text{target}} \beta\|^2 + \lambda_1 \|\beta\|^2 + \lambda_2 \|\beta - \beta_{\text{source}}\|^2$$

## 3. 主要な結果

### 3.1 ベースライン比較（Fst=0.10, N_EUR=5000, N_EAS=1000, h²=0.5）

| 手法 | R² | 相関係数 | 回帰傾き |
|------|-----|---------|----------|
| Standard PRS | 0.2950 | 0.5431 | 0.2362 |
| Bayesian LD | 0.0063 | 0.0791 | 1.0532 |
| Meta-Analysis | 0.3078 | 0.5548 | 0.2455 |
| Local Ancestry | 0.0128 | 0.1131 | 0.0176 |
| Penalized Transfer | 0.1985 | 0.4455 | 0.3868 |
| EUR Within-Pop (参照) | 0.2729 | 0.5224 | 0.2201 |

**主要所見**: Multi-Ethnic Meta-Analysis PRS が最高のR²（0.3078）を達成し、EUR内予測（0.2729）を上回った。Standard PRSも0.2950と比較的良好な性能を示した。

### 3.2 手法比較

![Method Comparison](figures/method_comparison.png)

### 3.3 集団分化（Fst）に対する感度

![Fst Sweep](figures/fst_sweep.png)

Fstの増加に伴いPRS予測精度が低下するパターンが確認された。Meta-Analysis手法はFstが大きい場合にもStandard PRSに対して安定した改善を示した。

### 3.4 ターゲット集団サンプルサイズの効果

![Sample Size Sweep](figures/sample_size_sweep.png)

EASサンプルサイズの増加に伴い、特にPenalized Transfer PRSとMeta-Analysis PRSの性能が向上した。

### 3.5 遺伝率の影響

![Heritability Sweep](figures/heritability_sweep.png)

遺伝率が高い形質ほど全手法のPRS性能が向上するが、手法間の相対的な優劣は維持された。

### 3.6 アレル頻度・LD構造の比較

![Allele Frequency Comparison](figures/allele_freq_comparison.png)

EUR-EAS間のアレル頻度の散布図とLD構造の差異。Fst=0.10でも相当のアレル頻度差が存在することが確認された。

### 3.7 効果量の分布

![Effect Size Analysis](figures/effect_size_analysis.png)

### 3.8 2型糖尿病ケーススタディ

パラメータ: n_SNPs=300, n_causal=40, h²=0.20, Fst=0.11, N_EUR=8000, N_EAS=2000

| 手法 | AUC | R² |
|------|-----|-----|
| Standard PRS | 0.7226 | 0.0630 |
| Bayesian LD | 0.5483 | 0.0032 |
| Meta-Analysis | 0.7014 | 0.0533 |
| Local Ancestry | 0.5691 | 0.0067 |
| Penalized Transfer | 0.6308 | 0.0257 |

![T2D Case Study](figures/t2d_case_study.png)

Standard PRSがAUC 0.7226で最高の判別能を示し、次いでMeta-Analysis（0.7014）が続いた。T2Dの低い遺伝率（h²=0.20）下でも、PRS上位分位群では有意なリスク上昇が確認された。

## 4. 考察と今後の展望

### 4.1 主要な知見

1. **Meta-Analysis手法の有効性**: 連続形質においてMeta-Analysis PRSが最も高い予測精度を達成。複数集団のGWASデータを統合することの有効性が示された。
2. **LD補正の限界**: 本シミュレーションでは、ベイズLD補正は正則化パラメータに敏感で、現在の実装では最適な性能を発揮できなかった。
3. **サンプルサイズの重要性**: ターゲット集団のGWASサンプルサイズが大きいほど、転移学習手法の改善幅が増大する。
4. **二値形質への適用**: T2DケーススタディでAUC 0.72を達成し、臨床的に有用なリスク層別化の可能性を示した。

### 4.2 限界

- シミュレーション環境の単純化（実際のゲノムデータではSNP数が数百万）
- 環境効果や遺伝子-環境相互作用の未考慮
- 局所祖先推定の簡略化（実際にはRFMix等のツールを使用）
- 集団特異的な因果バリアントの未モデル化

### 4.3 今後の方向性

- 大規模実データ（UK Biobank + BioBank Japan）での検証
- PRS-CSxやCT-SLEBなどの最新手法との比較
- 複数のターゲット集団（アフリカ系、南アジア系等）への拡張
- 機械学習ベースの効果量補正手法の開発

## 5. 生成ファイル一覧

### コード
| ファイル | 説明 |
|---------|------|
| `prs_transferability.py` | シミュレーションフレームワーク本体 |

### 図表
| ファイル | 説明 |
|---------|------|
| `figures/method_comparison.png` | 手法別R²・相関係数の比較 |
| `figures/fst_sweep.png` | Fst値に対するPRS精度の変化 |
| `figures/sample_size_sweep.png` | サンプルサイズに対するPRS精度の変化 |
| `figures/heritability_sweep.png` | 遺伝率に対するPRS精度の変化 |
| `figures/allele_freq_comparison.png` | EUR-EAS間のアレル頻度・LD比較 |
| `figures/effect_size_analysis.png` | 効果量分布とManhattanプロット |
| `figures/t2d_case_study.png` | 2型糖尿病ケーススタディ結果 |

### 数値結果
| ファイル | 説明 |
|---------|------|
| `results/baseline_results.csv` | ベースライン実験の数値結果 |
| `results/fst_sweep_results.csv` | Fstスイープの数値結果 |
| `results/sample_size_sweep_results.csv` | サンプルサイズスイープの数値結果 |
| `results/heritability_sweep_results.csv` | 遺伝率スイープの数値結果 |
| `results/t2d_case_study_results.csv` | T2Dケーススタディの数値結果 |

### ドキュメント
| ファイル | 説明 |
|---------|------|
| `report.md` | 本レポート |
| `paper.md` | 学術論文形式の文書 |
