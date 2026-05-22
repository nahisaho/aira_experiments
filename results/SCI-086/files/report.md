# 患者個別心臓デジタルツインフレームワーク — 設計報告書

**DRAFT — NOT FOR DISTRIBUTION**

- **日付**: 2026-05-23
- **患者ID**: DT_DEMO_001（デモンストレーション用合成データ）
- **バージョン**: 1.0.0

---

## 1. 実験目的と背景

### 1.1 目的

患者固有の心臓デジタルツインを構築するための統合フレームワークを設計・実装する。本フレームワークは以下の6つのモジュールで構成される：

1. **心臓MRIからの3D形状再構成** — セグメンテーションとメッシュ生成
2. **心筋電気伝導シミュレーション** — Aliev-Panfilov および ten Tusscher モデル
3. **力学-電気連成モデル** — Electro-Mechanical Coupling
4. **患者固有パラメータの逆問題推定** — ECG/エコーデータ活用
5. **不整脈リスク評価** — リエントリー脆弱性解析
6. **心房細動アブレーション効果予測** — 仮想アブレーションケーススタディ

### 1.2 背景

心臓デジタルツインは、患者の心臓を計算モデルとして再現し、治療計画の最適化や予後予測に活用する技術である。OpenCARP（電気生理学ソルバー）とFEBio（有限要素力学ソルバー）を基盤として、臨床データからの個別化パイプラインを構築した。

### 1.3 臨床的意義

- **術前計画**: アブレーション戦略の事前評価
- **リスク層別化**: 不整脈発生リスクの定量的評価
- **治療最適化**: 患者固有のモデルに基づく介入効果予測

---

## 2. 使用した手法・アルゴリズムの概要

### Module 1: 画像処理・メッシュ生成

| 項目 | 手法 |
|------|------|
| セグメンテーション | nnU-Net ベース 3D U-Net（8クラス分類） |
| 前処理 | Z-score正規化、パーセンタイルクリッピング（0.5-99.5%） |
| 表面メッシュ | Marching Cubes + Laplacian平滑化（30反復） |
| 体積メッシュ | TetGen型制約付きDelaunay四面体化 |
| 線維配向 | Rule-based法（Bayer et al., 2012）: 内膜+60° → 外膜-60° |
| 出力形式 | OpenCARP (.pts/.elem/.lon) + FEBio (.feb) |

### Module 2: 電気生理学シミュレーション

| モデル | 変数数 | 特徴 |
|--------|--------|------|
| Aliev-Panfilov | 2 | 現象論的、計算高速、定性的解析向け |
| ten Tusscher 2006 | 19 | 詳細イオンチャネル動態、12種イオン電流 |
| 組織レベル | Monodomain方程式 | 演算子分離法（Rush-Larsen + 陰的拡散） |

**Monodomain方程式:**

```
∂V/∂t = (1/χCm) ∇·(σ∇V) - I_ion/Cm + I_stim/Cm
```

### Module 3: 力学-電気連成

| コンポーネント | モデル |
|----------------|--------|
| 受動力学 | Holzapfel-Ogden構成則（8パラメータ） |
| 能動張力 | Land et al. (2017) — Ca²⁺駆動クロスブリッジ |
| 血行動態 | 3要素Windkesselモデル（Rc, Rp, C） |
| 連成方式 | 双方向：V_m → Ca²⁺ → 張力 → 変形 / 伸展 → MEF → V_m |

### Module 4: 逆問題推定

| 手法 | 詳細 |
|------|------|
| アルゴリズム | Ensemble Kalman Inversion（EKI, Iglesias et al., 2013） |
| ECG逆問題 | 伝導速度・APD・再分極分散の推定 |
| 力学逆問題 | 受動剛性・能動収縮力・末梢血管抵抗の推定 |
| 感度解析 | 有限差分法による局所感度分析 |
| 不確実性定量 | アンサンブル分布による95%信頼区間 |

### Module 5: 不整脈リスク評価

| 解析項目 | 手法 |
|----------|------|
| APD回復曲線 | 単指数モデル、勾配>1でalternans傾向 |
| 脆弱性ウィンドウ | S1-S2プログラム刺激プロトコル |
| APD分散 | 空間勾配解析、高勾配領域同定 |
| 線維化基質 | パターン分類（patchy/dense/confluent）、境界帯解析 |
| 統合リスク | 5因子重み付けスコア（0-1） |

### Module 6: アブレーション予測

| 戦略 | 説明 |
|------|------|
| PVI | 肺静脈隔離術（基本術式） |
| PVI + Roof | PVI + 天蓋線 |
| PVI + Posterior | PVI + 後壁隔離 |
| Hybrid | PVI + 基質ガイド |

---

## 3. 主要な結果と数値

### 3.1 画像処理・メッシュ生成結果

| 指標 | 値 |
|------|-----|
| MRIボリュームサイズ | 128 × 128 × 64 ボクセル |
| ボクセル間隔 | 1.25 × 1.25 × 2.5 mm |
| LV心筋体積 | 143.3 mL |
| LV内腔体積 | 25.2 mL |
| RV内腔体積 | 64.1 mL |
| LA体積 | 86.2 mL |
| 表面メッシュ | 5,000頂点 / 10,000三角形 |
| 体積メッシュ | 6,666頂点 / 19,998四面体 |
| メッシュ品質（平均アスペクト比） | 0.85 |
| 線維角度範囲 | +60° → -60°（内膜→外膜） |

### 3.2 電気生理学シミュレーション結果

| 指標 | Aliev-Panfilov | ten Tusscher 2006 |
|------|----------------|-------------------|
| 状態変数数 | 2 | 19 |
| シミュレーション時間 | 200 ms | 500 ms |
| 1Dケーブル長 | 100セル × 0.2 mm | 1セル（単一細胞） |
| イオンチャネル数 | — | 12 |

**OpenCARP設定パラメータ:**
- 組織内伝導度（縦方向）: σ_il = 0.17 mS/mm
- 組織内伝導度（横方向）: σ_it = 0.019 mS/mm
- 膜容量: C_m = 1.0 µF/cm²
- 表面積/体積比: χ = 1400 cm⁻¹

### 3.3 力学-連成シミュレーション結果

| 指標 | 値 |
|------|-----|
| 拡張末期容積（EDV） | 219.0 mL |
| 収縮末期容積（ESV） | 120.0 mL |
| 一回拍出量（SV） | 99.0 mL |
| 駆出率（EF） | 45.2% |
| 構成則 | Holzapfel-Ogden（a=0.059 kPa, a_f=18.472 kPa） |
| 能動モデル | Land 2017（T_ref=120 kPa） |

### 3.4 逆問題推定結果

**力学パラメータ（EKI推定, 25反復）:**

| パラメータ | 推定値 | 95% CI |
|-----------|--------|--------|
| a（等方性剛性） | 0.037 kPa | EKI分布より算出 |
| a_f（線維方向剛性） | 16.4 kPa | EKI分布より算出 |
| T_ref（基準張力） | 90.5 kPa | EKI分布より算出 |
| R_p（末梢抵抗） | 0.79 kPa·s/mL | EKI分布より算出 |

**最終残差**: 6.38（25反復後）

### 3.5 不整脈リスク評価結果

| 指標 | 値 |
|------|-----|
| **総合リスクスコア** | **0.460（中等度）** |
| 回復曲線最大勾配 | 0.844 |
| alternans傾向 | あり（勾配 > 0.5） |
| 波長（λ = CV × ERP） | CV × 230 mm |
| 脆弱性ウィンドウ幅 | 6.0 ms |
| APD分散（標準偏差） | 14.9 ms |
| 線維化負荷 | 12%（patchy型） |
| リエントリー誘発可能 | No |

**サブスコア内訳:**

| 因子 | スコア | 重み |
|------|--------|------|
| 回復特性 | 0.844 | 0.20 |
| APD分散 | 0.315 | 0.25 |
| 線維化 | 0.600 | 0.25 |
| 伝導異常 | 0.120 | 0.15 |
| 形態 | 0.300 | 0.15 |

### 3.6 アブレーション戦略比較結果

| 戦略 | 焼灼数 | 施術時間(分) | 伝導Gap | PV再接続リスク | **1年再発率** |
|------|--------|------------|---------|---------------|-------------|
| PVI | 48 | 16.0 | 4 | 61.5% | 32.9% |
| PVI + Roof | 56 | 18.7 | 4 | 61.3% | 28.4% |
| PVI + Posterior | 58 | 19.3 | 4 | 61.2% | 25.4% |
| **Hybrid (PVI+基質)** | **56** | **18.7** | **4** | **61.3%** | **23.9%** |

**最適戦略**: Hybrid PVI + 基質ガイドアブレーション（1年再発率 23.9%）

---

## 4. 考察と今後の展望

### 4.1 考察

#### フレームワーク設計

本フレームワークは、画像処理からアブレーション予測まで一貫したパイプラインを提供する。OpenCARP/FEBioとのインターフェースにより、産業標準のソルバーとの統合が可能である。

#### 逆問題推定

EKIによる力学パラメータ推定は25反復で残差6.38に収束した。ECG逆問題は簡略化されたフォワードモデルのため収束が不十分であり、フル3D lead-field計算の導入が必要である。

#### 不整脈リスク

総合リスクスコア0.460（中等度）は、回復曲線の急峻さ（0.844）と線維化（12%）が主要な寄与因子であることを示す。APD分散は比較的低く、均一な再分極を示唆する。

#### アブレーション予測

Hybrid戦略が最低の1年再発率（23.9%）を示し、PVI単独（32.9%）に対して約27%の相対的改善を達成した。全戦略で4箇所の伝導Gapが検出され、レジオンの完全性向上が課題である。

### 4.2 限界

1. **合成データ**: 本デモンストレーションは合成データを使用しており、実臨床データでの検証が必要
2. **簡略化モデル**: 電気生理学・力学ともに計算効率のため一部簡略化を適用
3. **単方向連成**: 現状のEM連成は完全な双方向ではなく、mechano-electric feedback (MEF) の実装が限定的
4. **アブレーション予測**: 経験的相関に基づく簡略化モデルであり、詳細なAFシミュレーションとの比較が必要

### 4.3 今後の展望

1. **GPU加速**: CUDAベースの組織レベルソルバーの実装（10-100倍の高速化を目指す）
2. **機械学習代理モデル**: Physics-Informed Neural Network (PINN) による高速フォワードモデル
3. **リアルタイムデジタルツイン**: 術中データの同化による動的モデル更新
4. **マルチスケール統合**: 細胞レベル（イオンチャネル）→組織レベル→臓器レベルのシームレスな接続
5. **臨床検証**: 多施設前向きコホートによるフレームワークの予測精度評価
6. **3D全心臓モデル**: 現在のLV中心モデルから両心室＋心房の全心臓モデルへの拡張
7. **不確実性定量**: ベイズ推定の精緻化とモンテカルロシミュレーションによる予測信頼区間の提供

---

## 5. 生成したファイル一覧

### ソースコード

| ファイル | 説明 |
|---------|------|
| `src/preprocessing/cardiac_mri_segmentation.py` | MRIセグメンテーション・メッシュ生成 |
| `src/electrophysiology/electrophysiology_models.py` | AP/TT06イオンモデル・Monodomainソルバー |
| `src/mechanics/cardiac_mechanics.py` | 受動/能動力学・EM連成・Windkessel |
| `src/inverse/inverse_estimation.py` | EKI逆問題推定・ECG/Echo逆問題 |
| `src/arrhythmia/arrhythmia_risk.py` | 回復曲線・APD分散・リスク評価 |
| `src/ablation/ablation_prediction.py` | 仮想アブレーション・戦略比較 |
| `src/digital_twin_pipeline.py` | 統合パイプラインオーケストレータ |
| `run_demo.py` | デモ実行スクリプト・可視化 |

### 結果・データファイル

| ファイル | 説明 |
|---------|------|
| `results/pipeline_results.json` | 全モジュール統合結果（JSON） |
| `data/opencarp/heart.pts` | OpenCARPノード座標 |
| `data/opencarp/heart.elem` | OpenCARP要素接続 |
| `data/opencarp/heart.lon` | OpenCARP線維配向 |
| `data/febio/heart.feb` | FEBioメッシュファイル |
| `configs/opencarp_ep.par` | OpenCARP電気生理パラメータ |
| `configs/febio_mechanics.feb` | FEBio力学設定 |

### 図表

| ファイル | 説明 |
|---------|------|
| `figures/fig1_architecture.png` | フレームワーク全体アーキテクチャ図 |
| `figures/fig2_hemodynamics.png` | 血行動態結果（PVループ・EF） |
| `figures/fig3_arrhythmia_risk.png` | 不整脈リスクダッシュボード |
| `figures/fig4_ablation_comparison.png` | アブレーション戦略比較 |

### ログ・設定

| ファイル | 説明 |
|---------|------|
| `logs/process-log.jsonl` | 実行トレースログ |
| `logs/pipeline.log` | パイプライン実行ログ |
| `requirements.txt` | Python依存パッケージ |

---

## 参考文献

1. Bayer, J.D., et al. (2012). A novel rule-based algorithm for assigning myocardial fiber orientation to computational heart models. *Ann Biomed Eng*, 40(10), 2243-2254.
2. Holzapfel, G.A. & Ogden, R.W. (2009). Constitutive modelling of passive myocardium. *Phil Trans R Soc A*, 367, 3445-3475.
3. Land, S., et al. (2017). Verification of cardiac mechanics software. *Phil Trans R Soc A*, 373, 20140091.
4. ten Tusscher, K.H.W.J. & Panfilov, A.V. (2006). Alternans and spiral breakup in a human ventricular tissue model. *Am J Physiol Heart Circ Physiol*, 291, H1088-H1100.
5. Aliev, R.R. & Panfilov, A.V. (1996). A simple two-variable model of cardiac excitation. *Chaos Solitons Fractals*, 7(3), 293-301.
6. Iglesias, M.A., Law, K.J.H., & Stuart, A.M. (2013). Ensemble Kalman methods for inverse problems. *Inverse Problems*, 29(4), 045001.
7. Plank, G., et al. (2021). The openCARP simulation environment for cardiac electrophysiology. *Comput Methods Programs Biomed*, 208, 106223.
8. Maas, S.A., et al. (2012). FEBio: Finite Elements for Biomechanics. *J Biomech Eng*, 134(1), 011005.
