# 実験レポート: 無細胞タンパク質合成（CFPS）生産性最適化フレームワーク

---

## 1. 実験目的と背景

### 1.1 目的

本研究では、無細胞タンパク質合成（Cell-Free Protein Synthesis; CFPS）システムの生産性を最大化するための統合計算フレームワークを設計・実装する。具体的には以下の6課題を扱った：

1. **転写-翻訳連成モデル**（リソース競合を考慮したODEシステム）
2. **エネルギー再生系の比較**（クレアチンリン酸、PEP、マルトース）
3. **Mg²⁺/K⁺/ポリアミン濃度の最適化マップ**
4. **mRNA安定性とリボソーム負荷の予測モデル**
5. **バッチ→半連続→連続系のスケールアップ設計**
6. **膜タンパク質発現（ナノディスク統合）のケーススタディ**

### 1.2 背景

CFPS技術は、細胞の生存制約から解放されたタンパク質生産プラットフォームである。近年、ワクチン製造、膜タンパク質構造解析、非天然アミノ酸組み込みなど多様な応用が急速に進んでいる（Zhu et al. 2025; Warfel et al. 2023）。しかし、エネルギー枯渇、イオン濃度最適化、リボソーム飽和などの複数の制約が同時に存在するため、システム全体の最適化は依然として困難である。本研究では、ODEモデルとベイズ最適化を統合することでこの問題に取り組む。

---

## 2. 先行研究サーベイの結果

### 2.1 特定された主要論文（2020年以降）

| # | タイトル | 著者 | 年 | DOI |
|---|---------|------|----|----|
| 1 | AI-driven high-throughput droplet screening of cell-free gene expression | Zhu et al. | 2025 | 10.1038/s41467-025-58139-0 |
| 2 | Cell-Free Protein Synthesis as a Method to Rapidly Screen ML-Generated Protease Variants | Thornton et al. | 2025 | 10.1021/acssynbio.5c00062 |
| 3 | A Low-Cost, Thermostable CFPS Platform (maltodextrin energy substrate) | Warfel et al. | 2023 | 10.1021/acssynbio.2c00392 |
| 4 | A highly efficient human cell-free translation system | Aleksashin et al. | 2023 | 10.1261/rna.079825.123 |
| 5 | Breakthrough in K. phaffii CFPS: AOX1 promoter drives T7-independent expression | Zhang et al. | 2025 | 10.3724/abbs.2025115 |
| 6 | CFPS platform for pyrrolysine-based ncAA incorporation | Ranji Charna et al. | 2022 | 10.1002/biot.202200096 |
| 7 | High-yield CHO CECF system for membrane proteins | Thoring et al. | 2017 | 10.1038/s41598-017-12188-8 |
| 8 | Combining mechanistic and ML models for metabolic optimization | Zhang et al. | 2020 | 10.1038/s41467-020-17910-1 |

### 2.2 先行研究の主要知見

- **Zhu et al. 2025（DropAI）**: マイクロ流体液滴とAIを組み合わせ、4倍のコスト削減を達成。E. coli CFEシステムの組成を転移学習でB. subtilis系に適用可能。
- **Warfel et al. 2023**: マルトデキストリンをエネルギー基質・凍結乾燥保護剤として使用。室温で4週間安定、ワクチン1回分約0.50ドル。
- **Zhang et al. 2025（K. phaffii）**: K⁺グルタミン酸とMg²⁺グルタミン酸の相乗効果によりGFP収量596 mg/L達成（CFPS最高記録）。
- **Thoring et al. 2017（CHO CECF）**: CHO細胞溶解物ベースのCECF系で膜タンパク質最大980 µg/mL達成。

### 2.3 先行研究の課題・限界

1. 個別コンポーネントの最適化に終始し、系全体の統合最適化が欠如
2. ODEモデルとベイズ最適化の組み合わせが未検討
3. エネルギー系、イオン条件、スケールアップ、膜タンパク質発現の比較が不十分
4. 高次元パラメータ空間の探索に必要なデータ量（>500サンプル）が実験的に困難

---

## 3. 使用した手法・アルゴリズムの概要

### 3.1 ODEモデル（7状態変数）

7つの状態変数を定義した連成ODE系：

| 変数 | 意味 | 単位 |
|------|------|------|
| D | DNAテンプレート | nM |
| M | mRNA | nM |
| R | 遊離リボソーム | nM |
| RL | リボソーム-mRNA複合体 | nM |
| P | タンパク質 | nM |
| E | エネルギー（ATP当量） | mM |
| AA | アミノアシルtRNA池 | 比率 |

主要方程式：

$$\frac{dM}{dt} = \frac{k_{tx,max} \cdot D \cdot E/(E+E_{thresh})}{1 + D/K_{m}} - k_{dm} \cdot M$$

$$\frac{dP}{dt} = k_p \cdot RL \cdot AA \cdot \frac{E}{E+E_{thresh}} - k_{dp} \cdot P$$

- 数値積分: Radau法（scipy.integrate.solve_ivp）
- 相対誤差 $10^{-6}$, 絶対誤差 $10^{-9}$

### 3.2 ベイズ最適化

- サロゲートモデル: ガウス過程（Matérn 5/2カーネル）
- 獲得関数: Expected Improvement（EI）
- 探索空間: 5次元（Mg²⁺, K⁺, スペルミン, エネルギー再生速度, リボソーム量）
- 初期ランダム評価: 12回 → BO反復: 30回

### 3.3 交差検証

- データセット: ODEモデル評価 + 18% CVのガウスノイズ（150サンプル）
- 手法: 5分割交差検証（KFold, shuffle=True）
- 評価指標: R², RMSE

---

## 4. 主要な結果と数値

### 4.1 エネルギー再生系比較

![Figure 1: エネルギー系比較](figures/fig1_energy_comparison.png)

**表1: エネルギー系別タンパク質収量（4時間シミュレーション）**

| エネルギー系 | 初期ATP (mM) | 最大再生速度 (mM/s) | 最終収量 (nM) | 換算収量 (µg/mL)* |
|------------|------------|-----------------|------------|-----------------|
| クレアチンリン酸 | 25 | 0.025 | **3,866** | 104 |
| PEP | 30 | 0.018 | 2,843 | 77 |
| マルトース | 20 | 0.010 | 1,594 | 43 |

*GFP（27 kDa）換算

**解釈:** クレアチンリン酸が最高収量を達成するが、高速な酵素失活と無機リン酸蓄積が制限因子となる。マルトースは緩やかな減衰プロファイルを示し、長時間連続系（CECF）に適する。

---

### 4.2 イオン濃度最適化マップ

![Figure 2: イオン最適化マップ](figures/fig2_ion_optimization_map.png)

**表2: スペルミン濃度別最適イオン条件**

| スペルミン (mM) | 最適Mg²⁺ (mM) | 最適K⁺ (mM) | ピーク相対収量 (%) |
|--------------|-------------|------------|-----------------|
| 0.0 | 10.1 | 101 | 100 |
| 0.5 | 10.1 | 101 | 127 |
| 1.0 | 10.1 | 101 | 132 |
| 2.0 | 10.1 | 101 | 118 |

**解釈:** 全スペルミン濃度でMg²⁺≈10 mM、K⁺≈100 mMが最適値。スペルミン1 mMで最大32%の収量向上。2 mM超では阻害効果が現れる（文献値と一致）。

---

### 4.3 mRNA安定性とリボソーム負荷

![Figure 3: mRNA安定性とリボソーム負荷](figures/fig3_mrna_ribosome.png)

- **ポリソーム占有率90%**: mRNA半減期が占有率10%の場合と比較して約3.4倍延長
- **リボソーム飽和点**: 総リボソーム40 nMに対して、mRNA濃度~20–40 nMでリボソーム負荷が飽和

---

### 4.4 スケールアップ設計比較

![Figure 4: スケールアップ比較](figures/fig4_scaleup.png)

**表3: 運転モード別性能比較**

| モード | 運転時間 | 最終収量 (nM) | 収量 (µg/mL) | バッチ比 |
|------|---------|------------|-----------|--------|
| バッチ | 2時間 | 1,605 | 43.3 | 1.0× |
| 半連続 | 8時間 | **6,389** | 172.5 | **3.98×** |
| 連続（定常状態） | 24時間 | 1,269 SS | 34.3 | 0.79× SS |

*SS = 定常状態瞬時濃度（容積生産性は連続系が優位）

---

### 4.5 膜タンパク質発現ケーススタディ

![Figure 5: 膜タンパク質発現](figures/fig5_membrane_protein.png)

- **最適ナノディスク濃度**: 5.0 µM
- **最適界面活性剤濃度**: 0.31 mM  
- **最大機能的収量**: 1,234 nM（~33 µg/mL）
- ナノディスク非存在下では凝集による分解が主要制限因子

---

### 4.6 ベイズ最適化結果

![Figure 6: ベイズ最適化](figures/fig6_bayesian_optimization.png)

**表4: ベイズ最適化で特定された最適条件**

| パラメータ | 最適値 |
|---------|-------|
| Mg²⁺ | 8.84 mM |
| K⁺ | 107 mM |
| スペルミン | 1.14 mM |
| エネルギー再生速度 | 最大値付近 |
| リボソーム量 | 最大値付近 |
| **予測最大収量** | **3,919 nM** |

- 収束: 約25–30反復で安定化
- 初期ランダム探索(12点) → EIガイド探索(30点) = 計42評価

---

### 4.7 GPサロゲートモデルの交差検証

![Figure 7: 交差検証](figures/fig7_cross_validation.png)

**表5: 5分割交差検証結果（n=150）**

| 評価指標 | 平均 ± 標準偏差 |
|--------|--------------|
| R² | 0.255 ± 0.266 |
| RMSE (nM) | 547 ± 46 |

---

## 5. 考察と今後の展望

### 5.1 結果の解釈

**エネルギー系**: 文献と一致して、クレアチンリン酸が短時間反応で最高収量を示すが、マルトースのほうが長時間系に適する。

**イオン最適化**: Mg²⁺≈10 mM, K⁺≈100 mMという結果は、Zhang et al. [2025]やJewett研究室の報告と一致しており、モデルの妥当性を示す。

**スケールアップ**: 半連続系の4倍収量向上は文献値（72倍の改善例あり：Jackson et al. 2014）より保守的だが、本シミュレーションでは単純な2時間サイクルの4回反復のみを考慮している。最適化されたCECF設計ではさらなる改善が期待できる。

### 5.2 ⚠️ 自己批判的評価

**1. 合成データへの依存性**  
交差検証は「ODEモデル + ノイズ」というデータのみで実施した。これはGPがODEの近似を学習する循環的検証であり、実験データに対する汎化能力を示すものではない。真の検証には、実際のCFPS実験データ（例：Jewett研究室の公開データセット、K. phaffiiの596 mg/L達成実験）に対するフィッティングが必要である。

**2. 低いR²（0.255±0.266）の意味**  
この低いR²と大きな標準偏差は意図的に現実的なノイズ（CV≈18%）を加えた結果であり、5次元空間での150サンプルによるGP回帰の本質的困難を正直に反映している。高精度予測には500サンプル以上が必要と推定される（Zhu et al.がO(10⁴)液滴を使用していることと整合）。

**3. パラメータ設定の不確実性**  
ODEのレート定数は文献値と内部整合性から設定したが、実験フィッティングを行っていない。特に`k_dm`（mRNA分解速度）、`k_rl`（リボソーム結合速度）は系に大きく依存し、E. coliとCHO系で1桁以上異なる可能性がある。

**4. 過度に楽観的な収量値の可能性**  
クレアチンリン酸系での3,866 nM（~104 µg/mL）は文献の典型的なバッチ値（20–100 µg/mL）の上限に相当し、非現実的ではないが、エネルギー枯渇・副産物蓄積の詳細な動態が簡略化されているため、実測値では10–30%低下する可能性がある。

**5. 実世界への一般化可能性**  
- 本モデルはE. coliベースのCFPSを想定しているが、CHOやコムギ胚芽系では翻訳機構が根本的に異なる
- ナノディスク積分モデルは定性的パラメータに基づいており、標的膜タンパク質ごとにキャリブレーションが必要
- pH変動、転写因子の活性低下、プロテアーゼ活性などの実験的変数が未考慮

### 5.3 今後の展望

1. **実験データによるパラメータ推定**: ベイズ推定（MCMC）によるODEパラメータの事後分布推定
2. **高次元最適化**: DNA濃度、RNAポリメラーゼ量、シャペロン濃度、pHも加えた10次元以上の探索
3. **ハイブリッドモデル**: ODEから計算したリボソーム占有率などの特徴量をGPの入力に加えた物理ガイドML
4. **実験自動化との統合**: DropAI [Zhu et al. 2025]のような液滴マイクロ流体と本フレームワークを統合し、自動実験-最適化ループを構築
5. **多目的最適化**: 収量・コスト・スケーラビリティのパレート最適化

---

## 6. 生成したファイル一覧

| ファイル名 | 種別 | 説明 |
|----------|------|------|
| `cfps_model.py` | Python | CFPS ODEモデル・全解析スクリプト |
| `results_summary.json` | JSON | 主要数値結果のまとめ |
| `figures/fig1_energy_comparison.png` | 図 | エネルギー再生系比較 |
| `figures/fig2_ion_optimization_map.png` | 図 | イオン濃度最適化マップ |
| `figures/fig3_mrna_ribosome.png` | 図 | mRNA安定性・リボソーム負荷 |
| `figures/fig4_scaleup.png` | 図 | スケールアップ設計比較 |
| `figures/fig5_membrane_protein.png` | 図 | 膜タンパク質発現ケーススタディ |
| `figures/fig6_bayesian_optimization.png` | 図 | ベイズ最適化収束・結果 |
| `figures/fig7_cross_validation.png` | 図 | 交差検証パリティプロット |
| `paper.md` | 論文 | 学術論文形式のレポート（英語） |
| `report.md` | レポート | 実験全結果レポート（本ファイル） |

---

## 参考文献

1. Zhu J et al. AI-driven high-throughput droplet screening of cell-free gene expression. *Nat Commun.* 2025. DOI: 10.1038/s41467-025-58139-0
2. Thornton EL et al. Cell-Free Protein Synthesis as a Method to Rapidly Screen Machine Learning-Generated Protease Variants. *ACS Synth Biol.* 2025. DOI: 10.1021/acssynbio.5c00062
3. Warfel KF et al. A Low-Cost, Thermostable, Cell-Free Protein Synthesis Platform. *ACS Synth Biol.* 2023. DOI: 10.1021/acssynbio.2c00392
4. Aleksashin NA, Chang ST, Cate JHD. A highly efficient human cell-free translation system. *RNA.* 2023. DOI: 10.1261/rna.079825.123
5. Zhang Y et al. Breakthrough in K. phaffii CFPS. *Acta Biochim Biophys Sin.* 2025. DOI: 10.3724/abbs.2025115
6. Ranji Charna A et al. CFPS platform for pyrrolysine-based ncAA incorporation. *Biotechnol J.* 2022. DOI: 10.1002/biot.202200096
7. Thoring L et al. High-yield production of "difficult-to-express" proteins in CHO CECF system. *Sci Rep.* 2017. DOI: 10.1038/s41598-017-12188-8
8. Zhang J et al. Combining mechanistic and ML models for metabolic optimization. *Nat Commun.* 2020. DOI: 10.1038/s41467-020-17910-1
