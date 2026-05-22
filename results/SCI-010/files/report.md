# DRAFT — NOT FOR DISTRIBUTION

# ADC ペイロード・リンカー最適化プラットフォーム報告書

- 作成対象: `workspace/adc_platform.py`
- 乱数シード: 42
- 対象系: HER2 高発現腫瘍を想定した ADC 設計評価

## 実験目的と背景

本実装の目的は、ADC（Antibody-Drug Conjugate）の payload-linker 設計を、**DAR 分布、リンカー切断、バイスタンダー効果、血漿安定性と腫瘍放出の最適化、PK/PD、Monte Carlo 感度解析、HER2 指向ケーススタディ**の 7 モジュールで統合評価できる計算基盤として構築することです。

特に T-DXd 類似条件を意識しつつ、DAR・cleavable linker・腫瘍内放出・血漿曝露のトレードオフを数値的に比較し、どの設計変数が有効性と安全性を強く左右するかを可視化しました。

## 使用した手法・アルゴリズムの概要

1. **DAR 分布モデリング**  
   DAR 0–8 を二項分布と Poisson 近似で記述し、10,000 分子の Monte Carlo サンプリングを実施しました。各 DAR 種について clearance rate、therapeutic index、hydrophobicity penalty を算出しました。

2. **リンカー切断 ODE**  
   酸感受性加水分解、cathepsin B による Michaelis-Menten 切断、還元型 disulfide 切断を Plasma と Tumor/Lysosome 条件で 24 時間追跡しました。

3. **バイスタンダー効果 PDE**  
   1D reaction-diffusion 方程式を有限差分法で解き、腫瘍半径 1 mm、100 分割で自由薬物の空間分布を評価しました。

4. **安定性最適化**  
   `J = efficacy_score - toxicity_score` を目的関数とし、`k_tumor / k_plasma >= 50` の制約付きグリッド探索と Pareto front 抽出を行いました。

5. **PK/PD モデル**  
   Central / Peripheral / Tumor の 3 compartment に、腫瘍内 target binding と cell kill を組み込んだ ODE 系を 21 日サイクルで解きました。

6. **Monte Carlo 感度解析**  
   Latin Hypercube Sampling により 1000 条件を生成し、CL, Vc, k_release_tumor, EC50, DAR を ±30% 変動させました。

7. **HER2-targeted ケーススタディ**  
   `DAR=4 cleavable`、`DAR=8 cleavable`、`DAR=8 non-cleavable` を比較しました。

## 各モジュールの主要な結果と数値

### 1. DAR 分布
- 平均 DAR: **6.247**
- DAR 標準偏差: **1.172**
- DAR 3–4 の治療域に入る割合: **7.35%**
- 高 DAR 条件では平均積載量は大きい一方、治療域に入る粒子の割合は限定的でした。
- Figure: `figures/01_dar_distribution.png`

### 2. リンカー切断
- 酸感受性リンカーの腫瘍側 24 h 放出率: **0.998**
- Cathepsin B モデルの腫瘍/血漿選択性比: **7.12**
- 還元型 disulfide の血漿 24 h 放出率: **0.152**
- 腫瘍内放出は酸性・酵素条件で強く促進され、血漿側では比較的抑制されました。
- Figure: `figures/02_linker_cleavage_kinetics.png`

### 3. バイスタンダー拡散
- 24 h 推定バイスタンダー半径: **0.515 mm**
- 24 h 中心濃度: **0.005123 a.u.**
- 24 h 辺縁濃度: **0.000112 a.u.**
- 腫瘍中心から辺縁へ濃度勾配が形成され、膜透過性 payload の空間的伝播を再現しました。
- Figure: `figures/03_bystander_diffusion.png`

### 4. 血漿安定性 vs 腫瘍放出最適化
- 最適グリッド点: `k_plasma = 0.0010 /day`, `k_tumor = 10.0000 /day`
- 最良目的関数 `J`: **0.9990**
- 低 plasma cleavage / 高 tumor cleavage の領域が最適設計に対応しました。
- Figure: `figures/04_optimization_landscape.png`

### 5. PK/PD
- ADC plasma AUC: **490.447**
- Free drug plasma AUC: **4.159**
- Tumor drug AUC: **322.759**
- Day 21 viable cell fraction: **2.983e-07**
- 名目条件では強い腫瘍制御が予測され、腫瘍内放出が PD 効果を支配しました。
- Figure: `figures/05_pkpd_simulation.png`

### 6. Monte Carlo 感度解析
- 成功シミュレーション数: **1000 / 1000**
- 平均 tumor AUC: **325.261**
- Day 21 viable cell fraction 平均: **2.797e-07**
- 95% bootstrap CI: **[2.792e-07, 2.802e-07]**
- 感度係数では `k_release_tumor` と `DAR` が最も強く、`EC50` と `Vc` も cell kill に影響しました。
- Figure: `figures/06_monte_carlo_sensitivity.png`

### 7. HER2-targeted ADC ケーススタディ

| Scenario | Therapeutic index | Tumor AUC | Plasma exposure | Bystander radius (mm) | Day21 viable fraction |
|---|---:|---:|---:|---:|---:|
| DAR=4 cleavable | 174.725 | 110.275 | 2.079 | 0.354 | 8.266e-07 |
| DAR=8 cleavable | 58.994 | 322.759 | 4.159 | 0.515 | 2.983e-07 |
| DAR=8 non-cleavable | 27.767 | 37.979 | 1.040 | 0.283 | 5.091e-07 |

- 治療指数では **DAR=4 cleavable** が最良でした。
- 腫瘍 AUC と bystander radius は **DAR=8 cleavable** が最大でした。
- **DAR=8 non-cleavable** は血漿曝露は低いものの、腫瘍内放出と空間効果が小さく、総合性能は限定的でした。
- Figure: `figures/07_case_study_comparison.png`

## 考察と今後の展望

- 本解析では、**高 DAR は腫瘍内薬物負荷を増やす一方、治療指数は中程度 DAR で改善しうる**ことが示されました。
- 一方、**cleavable linker を伴う DAR=8 設計は腫瘍 AUC と bystander spread を最大化**し、不均一 HER2 発現腫瘍に有利となる可能性があります。
- 感度解析より、**k_release_tumor と DAR が優先的に最適化すべき設計変数**であることが示唆されました。
- 限界として、本モデルは 1D 拡散、単回 21 日サイクル、簡略化された target turnover、規格化濃度単位を採用しており、実験データによる較正が必要です。
- 今後は、(1) 実測 cathepsin 活性や pH 分布の導入、(2) 複数投与サイクル化、(3) 腫瘍不均一性と細胞集団分布の導入、(4) 実測 PK データに基づく Bayesian calibration が有用です。

## 生成したファイル一覧

### Figures
- `figures/01_dar_distribution.png`
- `figures/02_linker_cleavage_kinetics.png`
- `figures/03_bystander_diffusion.png`
- `figures/04_optimization_landscape.png`
- `figures/05_pkpd_simulation.png`
- `figures/06_monte_carlo_sensitivity.png`
- `figures/07_case_study_comparison.png`

### Results
- `results/dar_analysis.csv`
- `results/linker_kinetics.csv`
- `results/optimization_results.csv`
- `results/pkpd_timecourse.csv`
- `results/monte_carlo_results.csv`
- `results/case_study_summary.csv`
- `results/statistical-summary.md`
- `results/summary_metrics.json`

### Data
- `data/dar_monte_carlo_samples.csv`
- `data/linker_terminal_release_summary.csv`
- `data/bystander_diffusion_profiles.csv`
- `data/optimization_landscape.csv`
- `data/monte_carlo_sensitivity_coefficients.csv`
- `data/preprocessing-log.md`

### Logs
- `logs/process-log.jsonl`
