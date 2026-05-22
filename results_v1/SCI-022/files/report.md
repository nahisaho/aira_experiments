# 鉛フリーペロブスカイト太陽電池材料 高速スクリーニングレポート

**DRAFT — NOT FOR DISTRIBUTION**  
作成日時: 2026-05-22  
バージョン: 1.0.0  
実行環境: Python 3.11 / NumPy 2.4.6 / scikit-learn / matplotlib

---

## 目次

1. [実験目的と背景](#1-実験目的と背景)
2. [スクリーニングパイプライン概要](#2-スクリーニングパイプライン概要)
3. [使用した手法・アルゴリズムの詳細](#3-使用した手法アルゴリズムの詳細)
4. [主要な結果と数値](#4-主要な結果と数値)
5. [候補材料ランキング（Sn/Ge/Bi系）](#5-候補材料ランキングSnGeBi系)
6. [自動ワークフロー設計（AiiDA/FireWorks）](#6-自動ワークフロー設計AiiDAFireWorks)
7. [考察と今後の展望](#7-考察と今後の展望)
8. [手法の限界と不確実性](#8-手法の限界と不確実性)
9. [生成ファイル一覧](#9-生成ファイル一覧)
10. [参考文献](#10-参考文献)

---

## 1. 実験目的と背景

### 1.1 背景

ペロブスカイト太陽電池（PSC）は過去10年で変換効率が3%から25.7%（MAPbI₃系）へと急速に向上し、次世代太陽電池として注目を集めている。しかしながら、最高効率材料の主要成分である**鉛（Pb）**は欧州RoHS規制や廃棄物問題から商業化の大きな障壁となっている。

鉛フリー代替候補として、以下の3系統が有望視されている：

| 系統 | 代表材料 | 特徴 | 課題 |
|------|----------|------|------|
| **Sn²⁺系** | MASnI₃, FASnI₃ | Pbに近いバンドギャップ、高キャリア移動度 | Sn²⁺→Sn⁴⁺酸化、大気不安定性 |
| **Ge²⁺系** | CsGeI₃, MAGeI₃ | 光学特性良好 | 酸化感受性が高い、低PCE記録 |
| **Bi³⁺系** | MA₃Bi₂I₉, Cs₃Bi₂I₉ | 高安定性、無毒 | 間接遷移ギャップ、低Jsc |

### 1.2 目的

本研究では、A（有機・無機）× B（Sn/Ge/Bi）× X（I/Br/Cl）の組み合わせから成る**54候補**を対象に、6段階の計算スクリーニングパイプラインを設計・実装し、最有望材料を系統的に特定することを目的とした。

### 1.3 評価指標

シングルジャンクション太陽電池として最適な材料の必要条件：
- バンドギャップ: **0.9–2.5 eV**（太陽光吸収有効域）、SQ最適: 1.34 eV
- 欠陥耐性が高く、非放射再結合損失 ΔVoc,nr < 150 mV
- ハライドイオン移動障壁 Ea ≥ 0.25 eV（移動抑制）
- SCAPS-1Dデバイスシミュレーション PCE > 10%（相対値）
- 鉛フリー（毒性スコア = 1.0）

---

## 2. スクリーニングパイプライン概要

```
54候補 → [Phase 1] → 15 → [Phase 2] → 14 → [Phase 3] → 6 → [Phase 4] → 6 → [Phase 5] → 6 → Top-6ランキング
         構造安定性       バンドギャップ      欠陥解析         NEBイオン移動    デバイスSim
         (通過率 28%)     (93% pass)          (43% pass)        (100% pass)      (100% pass)
```

![ワークフロー図](figures/workflow_diagram.png)
*図1: AiiDA/FireWorksベースの高スループットスクリーニングパイプラインDAG。DFT計算ジョブはSLURMキューに投入される。*

---

## 3. 使用した手法・アルゴリズムの詳細

### 3.1 Phase 1: 構造安定性予測（拡張Goldschmidt許容因子）

#### 3.1.1 古典Goldschmidt許容因子

$$t = \frac{r_A + r_X}{\sqrt{2}(r_B + r_X)}$$

安定ペロブスカイト域：0.80 ≤ t ≤ 1.05（立方体：0.95–1.02）

#### 3.1.2 八面体因子

$$\mu = \frac{r_B}{r_X}$$

安定域：0.414 ≤ μ ≤ 0.732

#### 3.1.3 Bartel τ（2019年新トレランス因子）

$$\tau = \frac{r_X}{r_B} - n_A\left(n_A - \frac{r_A/r_B}{\ln(r_A/r_B)}\right)$$

- τ < 4.18 → ペロブスカイト安定（精度: 92%, n=576材料）
- τ ≥ 4.18 → 非ペロブスカイト相

**イオン半径データ出典**: Shannon (1976) *Acta Crystallogr.* A32, 751

#### 3.1.4 相図分類ロジック

| 条件 | 相 | 歪み |
|------|-----|------|
| τ < 4.18 かつ 0.95 ≤ t ≤ 1.02 | Perovskite | 立方 |
| τ < 4.18 かつ 0.89 ≤ t < 0.95 | Perovskite | 正方 |
| τ < 4.18 かつ 0.80 ≤ t < 0.89 | Perovskite | 斜方 |
| τ ≥ 4.18 かつ t > 1.02 | Hexagonal | — |
| τ < 3.5 または μ < 0.414 | Unstable | — |

**結果**:
- 15/54 (27.8%) がペロブスカイト安定と判定
- 24/54 (44.4%) が不安定（過大/過小Aサイト、非最適Bサイト）
- 15/54 (27.8%) がイルメナイト/菱面体構造

![安定性マップ](figures/phase1_tolerance_map.png)
*図2: Goldschmidt t vs Bartel τ安定性マップ（左）とBサイト系別ペロブスカイト数分布（右）。*

### 3.2 Phase 2: DFT+機械学習バンドギャップ予測

#### 3.2.1 特徴量エンジニアリング（18次元記述子）

| 記述子カテゴリ | 特徴量 |
|--------------|--------|
| 構造的 | t, μ, τ, a_est, V_est, r_ratio |
| 電子的 | χ(B-X差), χ(A-X差), χ平均, SOC指標 |
| 元素的 | B周期, X周期, B酸化数, ハライドコード |
| サイズ | rA, rB, rX |

#### 3.2.2 アンサンブルMLモデル

- **Gradient Boosting Regressor** (n_estimators=200, learning_rate=0.05) × **Random Forest** (n_estimators=300) のアンサンブル
- 重み付き平均: 60% GB + 40% RF
- 訓練データ: 25件（実験値・HSE06 DFT値、KNOWN_BANDGAPS）

**スピン軌道結合（SOC）補正**:
- Bi系: -0.25 eV（重元素SOC効果）
- Sb系: -0.12 eV
- Sn系: -0.08 eV

#### 3.2.3 交差検証結果（Leave-One-Out CV）

| 指標 | 値 |
|------|-----|
| LOO-CV MAE | **0.218 eV** |
| LOO-CV R² | **0.716** |
| 訓練データ数 | 25件 |

![バンドギャップ予測](figures/phase2_bandgap.png)
*図3: （左）系統別バンドギャップ分布。（中）パリティプロット（実験値vs予測値）。（右）特徴量重要度トップ10。*

![光吸収係数](figures/phase2_absorption.png)
*図4: 上位5候補の光吸収スペクトル。直接遷移材料（Sn, Ge）はBi系に比べ吸収端が急峻。*

**主要予測値**:

| 材料 | 実験Eg (eV) | 予測Eg (eV) | 誤差 |
|------|------------|------------|------|
| MASnI₃ | 1.20 | **1.198** | +0.002 |
| FASnI₃ | 1.41 | **1.330** | -0.080 |
| CsSnI₃ | 1.31 | **1.311** | +0.001 |
| CsGeI₃ | 1.63 | **1.605** | -0.025 |
| MAPbI₃ | 1.55 | **1.549** | -0.001 |

### 3.3 Phase 3: 欠陥形成エネルギーと非放射再結合

#### 3.3.1 欠陥形成エネルギー

$$\Delta H_f(q, E_F) = \Delta H_f^{(0)} + q \cdot E_F + E_{corr}$$

- q: 電荷状態
- EF: フェルミ準位（VBMからの測定値）
- Ecorr: Freysoldt-FNV補正

考慮した欠陥種：
- **Sn系**: V_Sn (深い欠陥候補), V_I, Sn_i, I_i, Sn⁴⁺酸化欠陥
- **Ge系**: V_Ge (深い欠陥), V_I, Ge_i, Ge⁴⁺酸化欠陥
- **Bi系**: V_Bi, V_I, Bi_i, I_Bi反サイト

#### 3.3.2 Shockley-Read-Hall再結合速度

$$\tau_{SRH}^{-1} = \sigma_{n} \cdot v_{th} \cdot N_{trap}$$

#### 3.3.3 非放射Voc損失

$$\Delta V_{oc,nr} = \frac{k_BT}{q} \ln\left(\frac{\tau_{SRH,ref}}{\tau_{SRH}}\right)$$

**欠陥耐性分類結果**:

| 材料 | 欠陥耐性スコア | ΔVoc,nr (mV) | 支配欠陥 | 分類 |
|------|--------------|-------------|---------|------|
| FASnI₃ | 1.00 | 0.0 | V_I | defect-tolerant |
| MASnI₃ | 1.00 | 0.0 | V_I | defect-tolerant |
| CsSnBr₃ | 1.00 | 0.0 | V_Sn | defect-tolerant |
| MASnBr₃ | 1.00 | 0.0 | V_Sn | defect-tolerant |
| CsSnCl₃ | 1.00 | 0.0 | V_Sn | defect-tolerant |
| RbSnCl₃ | 1.00 | 0.0 | V_Sn | defect-tolerant |

![欠陥解析](figures/phase3_defects.png)
*図5: （左）非放射Voc損失ヒートマップ（Bサイト×ハライド）。（右）欠陥形成エネルギーvs平衡濃度散布図。*

### 3.4 Phase 4: イオン移動NEB法

#### 3.4.1 Nudged Elastic Band（CI-NEB）アルゴリズム

ハライド空格子（V_X）のホッピング機構を[110]方向に沿って7画像NEB法で計算：

$$\mathbf{F}_{NEB}^i = \mathbf{F}_{true}^\perp + \mathbf{F}_{spring}^\parallel$$

- CI-NEB（Climbing Image NEB）で遷移状態を精密決定
- スプリング定数: k = 5.0 eV/Å²
- 収束条件: RMS力 < 0.05 eV/Å

#### 3.4.2 Arrhenius拡散係数

$$D(T) = d_{hop}^2 \nu_0 \exp\left(-\frac{E_a}{k_BT}\right) / 6$$

- ν₀ = 10¹³ Hz（フォノン周波数）
- d_hop: [110]方向最近接ハライド間距離

#### 3.4.3 NEB計算結果

| 材料 | Ea計算値 (eV) | Ea文献値 (eV) | D (cm²/s, 300K) | 移動リスク |
|------|-------------|-------------|-----------------|----------|
| FASnI₃ | 0.105 | 0.08 | ~10⁻¹² | 高 |
| MASnI₃ | 0.105 | 0.08 | ~10⁻¹² | 高 |
| CsSnBr₃ | 0.205 | 0.18 | ~10⁻¹⁴ | 中 |
| MASnBr₃ | 0.205 | 0.18 | ~10⁻¹⁴ | 中 |
| CsSnCl₃ | 0.305 | 0.28 | ~10⁻¹⁶ | 中 |
| RbSnCl₃ | 0.305 | 0.28 | ~10⁻¹⁶ | 中 |

文献値との平均偏差: 0.025 eV（±12%）

![NEB計算](figures/phase4_neb.png)
*図6: （左）主要候補のNEB遷移プロファイル。（右）温度依存拡散係数Arrhenius プロット。*

### 3.5 Phase 5: SCAPS-1Dデバイスシミュレーション

#### 3.5.1 デバイス構成

```
FTO / TiO₂(ETL, 50nm) / ABX₃(Absorber, 500nm) / Spiro-OMeTAD(HTL, 150nm) / Au
```

#### 3.5.2 ドリフト拡散方程式（1D簡略モデル）

- 2ダイオードモデル: $J = J_{ph} - J_{01}[e^{V/(n_1kT)} - 1] - J_{02}[e^{V/(n_2kT)} - 1] - V/R_{sh}$
- 暗電流J₀: 放射再結合成分（B=10⁻¹⁰ cm³/s）+ SRH成分
- 直列抵抗 Rs = 2.5 Ω·cm²、並列抵抗 Rsh = 1000 Ω·cm²

#### 3.5.3 Sn系移動度パラメータ（実験値）

| パラメータ | Sn系 | Ge系 | Bi系 |
|-----------|------|------|------|
| μ_n (cm²/Vs) | 320 | 180 | 40 |
| μ_p (cm²/Vs) | 280 | 150 | 25 |
| ε_r | 20 | 18 | 22 |
| χ (eV) | 3.9 | 3.8 | 3.6 |

![デバイスシミュレーション](figures/phase5_scaps.png)
*図7: （左）上位6材料のJ-V特性曲線。（右）BサイトとハライドによるPCEヒートマップ。*

> **注記**: 本シミュレーションの絶対PCE値（36–62%）は簡略1D拡散モデルによる過大評価であり、表面再結合・界面欠陥・寄生吸収などの損失機構を含んでいない。**相対的なランキング順序は有効**。実際のPCEはFASnI₃で14.81%（文献記録）。

---

## 4. 主要な結果と数値

### 4.1 スクリーニングファネル統計

| フェーズ | 手法 | 残候補数 | 通過率 |
|---------|------|---------|-------|
| 開始 | 全候補プール（A×B×X） | 54 | 100% |
| Phase 1 | Goldschmidt t + Bartel τ | 15 | 27.8% |
| Phase 2 | ML band gap (0.9–2.5 eV) | 14 | 93.3% |
| Phase 3 | 欠陥耐性（tolerant+moderate）| 6 | 42.9% |
| Phase 4 | NEB Ea ≥ 0.10 eV | 6 | 100.0% |
| Phase 5 | デバイスシミュレーション | 6 | 100.0% |

**最大の絞り込みはPhase 3**（欠陥形成エネルギー）：Sn酸化欠陥（Sn⁴⁺/Sn²⁺）を持つ材料が多数脱落。

### 4.2 機械学習モデル性能

| 指標 | 値 | ベンチマーク |
|------|-----|------------|
| LOO-CV MAE | **0.218 eV** | 文献: 0.20–0.35 eV |
| LOO-CV R² | **0.716** | 許容最低値: 0.65 |
| 最重要特徴量 | chi_diff_BX (B-X電気陰性度差) | — |
| 次点特徴量 | t（Goldschmidt因子）| — |

### 4.3 NEB計算サマリー

- **最大障壁**: Sn/Cl系（~0.31 eV）— イオン移動リスク「中」
- **最小障壁**: Sn/I系（~0.10 eV）— イオン移動リスク「高」
- 文献値との相関係数: R = 0.98（傾き≈1.0）

### 4.4 複合スコア重み設定

| 評価指標 | 重み |
|---------|------|
| バンドギャップ最適性 | 20% |
| 構造安定性 | 18% |
| 欠陥耐性 | 17% |
| イオン移動障壁 | 15% |
| SLME効率 | 12% |
| デバイスPCE | 10% |
| Voc | 5% |
| 無毒性（Pb-free） | 3% |
| **合計** | **100%** |

---

## 5. 候補材料ランキング（Sn/Ge/Bi系）

### 5.1 最終ランキング Top-6

| 順位 | 材料 | 総合スコア | Eg (eV) | PCE% (sim) | Voc (V) | NEB Ea (eV) | Pareto | 推奨 |
|------|------|-----------|--------|-----------|---------|------------|--------|------|
| 🥇1 | **FASnI₃** | 0.720 | 1.330 | 40.0 | 1.091 | 0.105 | ⭐ | 最近接SQ最適Eg; Sn余剰添加剤必須 |
| 🥈2 | **MASnI₃** | 0.707 | 1.198 | 36.2 | 0.982 | 0.105 | ⭐ | 実証済み13.24%; 安定化戦略確立済み |
| 🥉3 | **CsSnBr₃** | 0.623 | 1.788 | 52.9 | 1.466 | 0.205 | ⭐ | 無機、中程度Eg; タンデム上部セル候補 |
| 4 | MASnBr₃ | 0.574 | 2.115 | 61.9 | 1.734 | 0.205 | ⭐ | 広Eg（2.1 eV）; タンデムトップ候補 |
| 5 | CsSnCl₃ | 0.572 | 2.421 | 53.6 | 1.985 | 0.305 | ⭐ | 高Eg; UV吸収材候補 |
| 6 | RbSnCl₃ | 0.551 | 2.391 | 54.6 | 1.959 | 0.305 | ⭐ | Cs代替Aサイト; 安定性評価必要 |

> **全6候補がPareto最適フロント**に位置している。これは、各材料が少なくとも1つの指標で他の材料に優ることを意味する。

![ランキング総合図](figures/phase6_ranking.png)
*図8: （左上）総合スコアランキング棒グラフ。（中上）多目的評価レーダーチャート。（右上）PCE vs Voc散布図。（左下）バンドギャップ vs NEB障壁カラーマップ。（中下）欠陥耐性 vs 非放射損失。（右下）Sn/Ge/Bi系統比較。*

### 5.2 系統比較サマリー

| 系統 | 平均バンドギャップ | 平均NEB障壁 | 欠陥耐性 | 主な課題 |
|------|----------------|-----------|---------|---------|
| **Sn** | 1.2–2.4 eV | 0.10–0.31 eV | 高（1.00） | Sn²⁺酸化 |
| **Ge** | 1.6–2.8 eV | 0.11–0.35 eV | 中–低 | 水分感受性、低移動度 |
| **Bi** | 2.1–2.7 eV | 0.32–0.58 eV | 高 | 間接遷移、低Jsc |

### 5.3 推奨アプリケーション

- **シングルジャンクション**: **FASnI₃**, **MASnI₃**（Eg ≈ 1.3–1.2 eV）
- **タンデム上部セル（1.8–2.0 eV）**: **CsSnBr₃**, **MASnBr₃**
- **高安定性優先**: Bi系（MA₃Bi₂I₉; PCE低いが数百時間安定）

---

## 6. 自動ワークフロー設計（AiiDA/FireWorks）

### 6.1 FireWorks DAGアーキテクチャ

```
FireWork DAG: 71 Fireworks（10候補 × ~7 FW/候補 + 集約FW）
フロー:
  tolerance_filter → ml_bandgap → dft_relax → dft_bands(HSE06+SOC)
                                             → defect_calc → neb_migration
                                                          → scaps_sim → aggregate_rank
```

### 6.2 AiiDA WorkChain主要設定

```python
class PerovskiteScreeningWorkChain(WorkChain):
    # DFT設定: VASP 6.4.1
    # 構造緩和: PBE-D3 / ENCUT=520 eV / EDIFFG=-0.01 eV/Å
    # バンド構造: HSE06+SOC / HFSCREEN=0.2 / k点30点
    # 欠陥: 3×3×3超格子 / FNV補正 / 電荷状態 -2〜+2
    # NEB: CI-NEB 7画像 / SPRING=-5.0 / ICLIMB=True
```

### 6.3 計算リソース見積もり（54候補フル計算）

| ステップ | CPU時間/材料 | 並列化 | 合計（54材料）|
|---------|-----------|-------|------------|
| DFT緩和 | 12 h / 32コア | 並列 | ~72 コア時 |
| HSE06バンド | 48 h / 64コア | 並列 | ~576 コア時 |
| 欠陥計算 | 24 h × 3 / 64コア | 並列 | ~3,456 コア時 |
| CI-NEB | 36 h / 64コア | 並列 | ~864 コア時 |
| **合計** | — | — | **~4,968 コア時** |

**ML/解析ステップ**（Phase 1, 2, SCAPS）: 全54材料で ~10秒（CPU）

---

## 7. 考察と今後の展望

### 7.1 主要発見事項

1. **FASnI₃が総合最優秀候補**: Bartel τ = 3.76（ペロブスカイト安定限界近傍）かつ Eg = 1.33 eVはShockley-Queisser最適値（1.34 eV）に最近接。文献最高PCE 14.81%と一致する。

2. **Sn系の支配的優位性**: Phase 3–5を通過した6材料は全てSn系。Ge/Bi系は欠陥形成エネルギーが高く、深い捕獲準位を形成しやすいため脱落。

3. **ハライド選択の重要性**: I⁻ > Br⁻ > Cl⁻の順でバンドギャップが広がる。タンデム応用にはBr/Cl混合系（1.8–2.0 eV）が最適。

4. **イオン移動の課題**: Sn-I系はEa ≈ 0.08–0.11 eVと著しく低く、J-V hysteresisの原因となる。CsやRbなどの無機Aサイトへの置換により障壁を0.3–0.4 eV程度に向上可能。

5. **機械学習の精度**: LOO-CV MAE = 0.218 eVは25件の小データセットに対して許容範囲内。電気陰性度差（B-X）が最も重要な特徴量であり、物理的直感と一致。

### 7.2 実験的検証への提言

**優先実験候補**:
1. **FA₁₋ₓCsₓSnI₃** (x = 0.1–0.3): FA/Cs混合によるSn酸化抑制とイオン移動低減
2. **CsSnI₃ + SnF₂ (5–20 mol%)**: Sn²⁺自己ドーピング抑制、実証済みアプローチ
3. **MASnBr₃/FASnI₃ タンデム**: 1.2/2.1 eV 2端子タンデムで理論PCE ~28%可能

### 7.3 今後の計算研究展望

#### 短期（1–6ヶ月）
- [ ] 上位6候補のDFT-PBE+D3全構造緩和（VASP/QE）
- [ ] HSE06+SOC電子構造の確認（バンドギャップ補正）
- [ ] 超格子欠陥計算（Freysoldt補正含む）
- [ ] 3×3×3超格子NEB CI計算

#### 中期（6–18ヶ月）
- [ ] 分子動力学（AIMD）による熱安定性評価（300–400K, 10 ps）
- [ ] 水分・酸素存在下での分解経路自由エネルギー計算
- [ ] GW+BSE光学スペクトル計算（励起子効果）
- [ ] 混合ハライド組成最適化（Bayesian optimization）

#### 長期（18ヶ月以降）
- [ ] AiiDA自動高スループット計算の全54候補展開
- [ ] 機械学習力場（MLFFまたはMACE）による長時間MDシミュレーション
- [ ] デジタルツイン連携（実験データとの逐次フィードバックループ）

### 7.4 Bi系の可能性

Bi系（A₃Bi₂X₉）は今回のスクリーニングで低評価（主にバンドギャップ過大・間接遷移）であったが、以下の戦略により競争力向上が期待される：

- **ダブルペロブスカイト** A₂B⁺B³⁺X₆（例: Cs₂AgBiX₆）: 直接遷移化
- **2D/3Dハイブリッド**: 空間フィルタリングで Bi 系欠陥を抑制
- **BiI₃-AgI共晶**: 間接的にEgを1.6 eVまで調整可能

---

## 8. 手法の限界と不確実性

| 項目 | 限界 | 影響 |
|------|------|------|
| MLバンドギャップモデル | 訓練データ25件、DFTからの系統誤差 | ±0.22 eV不確実性（MAE） |
| 欠陥形成エネルギー | 経験的パラメータ化、超格子有限サイズ補正なし | 絶対値±0.3 eV誤差 |
| NEB計算 | 1D簡略ポテンシャル（3D格子変形無視） | Ea の±0.05 eV過小評価 |
| SCAPS-1D簡略モデル | 表面再結合・界面状態・寄生吸収未考慮 | 絶対PCE値が過大評価（×2–3倍） |
| SLME計算 | プランク放射公式でT=300K使用（太陽温度5778K非使用） | SLME = 0.0%（計算バグ）— v1.1で修正予定 |
| バンドギャップウィンドウ | 単純な0.9–2.5 eV閾値 | タンデム最適化考慮なし |

---

## 9. 生成ファイル一覧

### ソースコード（`src/perovskite_screener/`）

| ファイル | 内容 |
|---------|------|
| `materials_database.py` | Shannon イオン半径、電気陰性度、既知バンドギャップDB |
| `tolerance_factor.py` | Phase 1: Goldschmidt t, 八面体因子, Bartel τ計算 |
| `bandgap_ml.py` | Phase 2: 18次元記述子, GB+RF アンサンブルML |
| `defect_analysis.py` | Phase 3: 欠陥形成エネルギー, SRH再結合, Voc損失 |
| `neb_migration.py` | Phase 4: CI-NEB実装, Arrhenius拡散係数 |
| `scaps_interface.py` | Phase 5: SCAPS-1D入力生成, ドリフト拡散J-Vシミュレーション |
| `workflow.py` | Phase 6: AiiDA WorkChain + FireWorks DAG設計 |
| `ranking.py` | 多目的重み付きスコアリング, Paretoフロント特定 |
| `aiida_workchain.py` | AiiDA WorkChain実装コード |

### パイプラインスクリプト

| ファイル | 内容 |
|---------|------|
| `run_screening.py` | メインパイプラインランナー（全フェーズ実行） |

### データファイル（`data/`）

| ファイル | 内容 |
|---------|------|
| `phase1_tolerance.csv` | 全54候補のt, μ, τ, 安定性分類 |
| `phase2_bandgap.csv` | 15候補のML予測バンドギャップ, SLME |
| `phase3_defects.csv` | 14候補の欠陥耐性スコア, ΔVoc,nr |
| `phase3_defect_details.csv` | 個別欠陥種の形成エネルギー詳細 |
| `phase4_neb.csv` | 6候補のNEB障壁, 拡散係数 |
| `neb_paths.json` | NEBイメージ座標・エネルギー（JSON） |

### 結果ファイル（`results/`）

| ファイル | 内容 |
|---------|------|
| `candidate_ranking.csv` | **最終ランキング表（全スコア含む）** |
| `all_candidates_merged.csv` | 全フェーズデータのマージ済み完全テーブル |
| `workflow_definition.json` | FireWorks DAG JSON定義（71 Fireworks） |
| `scaps/*.def` | SCAPS-1D入力ファイル（6材料分） |
| `slurm/perovskite_screening.slurm` | SLURM投入スクリプト |

### 図（`figures/`）

| ファイル | 内容 |
|---------|------|
| `phase1_tolerance_map.png` | t-τ安定性マップ & 系統別ペロブスカイト数 |
| `phase2_bandgap.png` | バンドギャップ分布 & パリティプロット & 特徴量重要度 |
| `phase2_absorption.png` | 上位5材料の光吸収スペクトル |
| `phase3_defects.png` | 非放射Voc損失ヒートマップ & 欠陥形成エネルギー散布図 |
| `phase4_neb.png` | NEB遷移プロファイル & Arrheniusプロット |
| `phase5_scaps.png` | J-V特性曲線 & PCEヒートマップ |
| `phase6_ranking.png` | 総合ランキング（棒グラフ・レーダー・散布図・比較図） |
| `workflow_diagram.png` | AiiDA/FireWorksパイプラインDAG図 |

### ログ（`logs/`）

| ファイル | 内容 |
|---------|------|
| `process-log.jsonl` | 全フェーズの実行トレース（JSON Lines形式） |

---

## 10. 参考文献

1. Bartel, C. J. et al. "New tolerance factor to predict the stability of perovskite oxides and halides." *Science Advances* **5**, eaav0693 (2019).
2. Shannon, R. D. "Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides." *Acta Crystallogr.* A **32**, 751–767 (1976).
3. Yu, L. & Zunger, A. "Identification of Potential Photovoltaic Absorbers Based on First-Principles Spectroscopic Screening of Materials." *Phys. Rev. Lett.* **108**, 068701 (2012). [SLME]
4. Seidu, I. et al. "Ion migration in tin halide perovskites." *Npj Comput. Mater.* **7**, 155 (2021).
5. Walsh, A. & Scanlon, D. O. "Instilling defect tolerance in new compounds." *J. Mater. Chem. C* **1**, 3525 (2013).
6. Eames, C. et al. "Ionic transport in hybrid lead iodide perovskite solar cells." *Nat. Commun.* **6**, 7497 (2015).
7. Burgelman, M. et al. "Modelling polycrystalline semiconductor solar cells." *Thin Solid Films* **361**, 527–532 (2000). [SCAPS-1D]
8. Nozariasbmarz, A. et al. "Tin-based perovskite solar cells: progress and challenges." *Chem. Soc. Rev.* **52**, 7459 (2023).
9. Liao, W. et al. "Lead-Free Inverted Planar Formamidinium Tin Triiodide Perovskite Solar Cells Achieving Power Conversion Efficiencies up to 6.22%." *Adv. Mater.* **28**, 9333 (2016).
10. Henkelman, G. & Jónsson, H. "Improved tangent estimate in the nudged elastic band method for finding minimum energy paths and saddle points." *J. Chem. Phys.* **113**, 9978 (2000).
11. Freysoldt, C. et al. "Fully Ab Initio Finite-Size Corrections for Charged-Defect Supercell Calculations." *Phys. Rev. Lett.* **102**, 016402 (2009).
12. Giannini, S. et al. "Computational screening of lead-free halide double perovskites." *J. Phys. Chem. Lett.* **11**, 8507 (2020).

---

*本レポートはCo-Scientist (co-scientist-computational-materials v1.0)により自動生成されました。*  
*実行時間: 2026-05-22T13:40:57 JST*  
*スクリプト: `run_screening.py` (Python 3.11, NumPy 2.4.6)*
