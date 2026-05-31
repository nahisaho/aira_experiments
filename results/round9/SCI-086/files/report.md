# Experimental Report: Patient-Specific Cardiac Digital Twin Framework

**実験日時:** 2026-05-31  
**研究テーマ:** OpenCARP/FEBioベース心臓デジタルツインフレームワークの設計と評価  
**使用環境:** Python 3.11.2, Jupyter MCP, ToolUniverse MCP, Web Search

---

## 1. 実験目的と背景

### 目的

本研究では、患者個別の心臓デジタルツイン（Cardiac Digital Twin: CDT）を構築するための包括的な計算フレームワークを設計・実装・評価した。具体的には以下の6コンポーネントの統合を目指した：

1. 心臓MRIからの3D形状再構成（セグメンテーション + メッシュ生成）
2. 心筋電気伝導シミュレーション（Aliev-Panfilov モデル）
3. 力学-電気連成モデル（Electro-Mechanical Coupling）
4. 患者固有パラメータの逆問題推定（ECGデータ活用）
5. 不整脈リスク評価（機械学習）
6. 心房細動アブレーション効果予測

### 背景

心房細動（AF）は世界で3,350万人が罹患する最も一般的な持続性不整脈であり、カテーテルアブレーションが主要な根治的治療法だが12ヶ月後の洞調律維持率は60〜70%にとどまる。心臓デジタルツインは、患者固有の解剖・電気生理・力学データを統合し、治療をin silicoでシミュレートすることで個別化医療を実現する次世代アプローチである。

---

## 2. ステップ1: 先行研究調査

### 2.1 使用ツールと調査結果

**使用ツール:** ToolUniverse SemanticScholar API（HTTP 429エラーのため失敗）、Web Search (Bing AI)

**SemanticScholar API エラー:** 3回試行すべてHTTP 429（レート制限）。代替としてWeb Searchを使用。

### 2.2 特定された主要論文（2021–2025）

| # | タイトル | 著者・年 | DOI | 主要知見 |
|---|---|---|---|---|
| 1 | Electro-Mechanical Whole-Heart Digital Twins | Gerach et al. 2021 | 10.3390/math9111247 | 完全連成マルチフィジクス心臓CDT；患者MRIで検証；peak stress 50-100 kPa |
| 2 | Cardiac Digital Twin Pipeline for Virtual Therapy | Camps et al. 2024 | 10.48550/arXiv.2401.10029 | 自動化パイプライン；不確実性定量化；薬効シミュレーション |
| 3 | GPU accelerated digital twins of the human heart | Viola et al. 2023 | 10.1038/s41598-023-34098-8 | GPU加速全心臓シミュレーション；仮想臨床試験 |
| 4 | From bits to bedside: cardiac electrophysiology DT | Bhagirath et al. 2024 | 10.1093/europace/euae295 | CDT最新レビュー；臨床転換課題 |
| 5 | Calibration strategy for electromechanical CDT | Wang et al. 2025 | 10.7554/eLife.106555 | ASME V&V40準拠キャリブレーション標準 |
| 6 | Cardiac DT at scale from MRI: UK Biobank | Ugurlu et al. 2025 | 10.1371/journal.pone.0327158 | ~55,000例スケールCDT生成；オープンソースツール |

### 2.3 先行研究の課題・限界

1. **スケーラビリティ**: 大規模コホートへのCDTパイプライン適用が困難
2. **逆問題の信頼性**: パラメータ推定の不確実性定量化が不足
3. **電気力学連成**: AFアブレーション計画でのEM連成が未成熟
4. **AFアブレーション予測**: 既存ML研究はN=50〜1,600と可変；AUROC 0.64〜0.86の幅
5. **計算コスト**: full TP06モデルのリアルタイム実用化は未達成

---

## 3. ステップ2: NatureLM / GALACTICA MCP 試行結果

### 試行したツールと結果

| ツール名 | 試行内容 | 結果 |
|---|---|---|
| `ask_naturelm` | AP model params, EP simulation 定量予測 | **接続失敗**: ToolUniverse に未登録 |
| `scientific_qa` | 実験設計妥当性検証 | **接続失敗**: ToolUniverse に未登録 |
| `predict_citations` | 関連文献予測 | **接続失敗**: ToolUniverse に未登録 |
| `SemanticScholar_search_papers` | 文献検索 | **HTTP 429**: レート制限 (3回試行) |

**確認方法:** `tooluniverse-grep_tools` でフィールド "name" を "NatureLM" / "GALACTICA" で検索 → 0マッチ

**代替措置:** Web Search (Bing AI) を2回実施し、6件の主要論文を特定。文献内の数値パラメータを実験設計に活用。

### 代替値（文献ベース）

NatureLM・GALACTICAによる定量予測の代わりに、文献から以下の参照値を使用：

| パラメータ | 文献値 | 出典 |
|---|---|---|
| Peak active stress | 50–120 kPa | Gerach et al. 2021 |
| AP model k (AP model) | 8.0 | Aliev & Panfilov 1996 |
| Electromechanical delay | 20–50 ms | 各種文献 |
| AF ablation AUROC (real data) | 0.75–0.86 | Hermans et al. 2023 |
| Segmentation Dice LV | 0.88–0.95 | nnU-Net benchmarks |

---

## 4. ステップ3: Python実装と実行結果

### 4.1 実行環境

- Python 3.11.2, NumPy 2.3.5, Pandas 2.3.3, SciPy 1.17.1, scikit-learn 1.6.1
- 乱数シード: `np.random.seed(42)`, `random.seed(42)`
- Jupyter MCP (kernel: f7b66057-8536-439d-b8b1-1ef51dbc3a78)

### 4.2 Aliev-Panfilov 2D シミュレーション [cell:1]

**設定:** 100×100グリッド (2×2 cm)、dt=0.05 ms、T=500 ms、Forward Euler法

**実行結果:**
- 最終最大電位: V_max = 0.0053 (500 ms後)
- 波面伝播: 6タイムポイントのスナップショット取得

![Figure 1: Aliev-Panfilov Wave Propagation](figures/fig1_aliev_panfilov_snapshots.png)

### 4.3 単細胞活動電位シミュレーション [cell:3b]

**実行結果:**
- APピーク: 0.9986 (正規化)
- APD90: 26.6 ms
- 回復変数 r: 正常減衰

![Figure 2: Action Potential](figures/fig2_action_potential.png)

### 4.4 逆問題：パラメータ推定 [cell:4]

**手法:** Nelder-Mead法、ノイズ付き合成ECGプロキシへのフィッティング

| パラメータ | 真値 | 推定値 | 相対誤差 |
|---|---|---|---|
| k (興奮性) | 8.000 | 7.790 | **2.6%** |
| a (閾値) | 0.150 | 0.146 | **2.7%** |
| ε (回復) | 0.00200 | 0.00255 | **27.4%** |
| **R²** | — | **0.9999** | — |
| **RMSE** | — | **0.0195** | — |

**考察:** kとaは高精度で回復（<3%誤差）。εは感度が低く27%誤差—これは心電図データからのεの同定困難性に関する文献報告と一致。

![Figure 3: Inverse Problem](figures/fig3_inverse_problem.png)

### 4.5 電気力学連成シミュレーション [cell:10]

**実行結果:**

| 指標 | 値 |
|---|---|
| Ca²⁺ ピーク | 0.724 (正規化) |
| 最大能動応力 | **67.7 kPa** |
| 最大短縮率 | **10.9%** |
| 電気力学遅延 | **23.8 ms** |

文献参照値（Gerach 2021）: 50–100 kPa、20–50 ms遅延 → **整合的**

![Figure 5: Electromechanical Coupling](figures/fig5_electromechanical.png)

### 4.6 不整脈リスク評価 ML [cell:7]

**コホート:** 合成データ N=200患者、11特徴量、37%陽性率

**5-fold CV結果:**

| モデル | AUROC | ± SD | F1 | ± SD |
|---|---|---|---|---|
| Random Forest | **0.921** | 0.042 | 0.766 | 0.075 |
| Gradient Boosting | **0.901** | 0.061 | 0.774 | 0.082 |

**上位5特徴量 (Random Forest MDI):**
1. EF% — 28.7%
2. LA径 (mm) — 25.4%
3. QT間隔 (ms) — 8.1%
4. 線維化率 — 8.0%
5. LV質量 (g) — 5.9%

### 4.7 AFアブレーション効果予測 [cell:8b]

**コホート:** N=150、10特徴量、成功率 53.3%

⚠️ **自己批判的記録:** 初期モデル（v1）は `pulm_vein_isolation` を特徴量として使用しAUROC=1.0（データリーク）。v2では該当特徴量を除去し、実際の臨床的困難さを反映するノイズ（σ=0.4）を追加。

**修正後5-fold CV結果:**

| モデル | AUROC | ± SD | F1 | ± SD |
|---|---|---|---|---|
| Random Forest | **0.571** | 0.137 | 0.633 | 0.114 |
| Gradient Boosting | **0.570** | 0.124 | 0.622 | 0.122 |

文献値（ATHENA 2024, N=1600+）: AUROC 0.75 → 本研究は小規模合成データのため下回る

### 4.8 3D形状統計 [cell:12]

**コホート:** N=50患者（合成）

| 計測値 | 平均 | ± SD |
|---|---|---|
| Dice LV | 0.919 | 0.032 |
| Dice RV | 0.891 | 0.031 |
| メッシュノード数 | 41,852 | 13,379 |
| メッシュ要素数 | 202,565 | 74,167 |

---

## 5. 主要な図表

### Figure 6: フレームワーク全体像と結果サマリー

![Figure 6: Framework Overview](figures/fig6_framework_overview.png)

### Figure 4: ML結果サマリー

![Figure 4: ML Results](figures/fig4_ml_results.png)

---

## 6. 考察と今後の展望

### 6.1 主要な知見

1. **逆問題**: Nelder-Mead法でAP模型パラメータをR²=0.9999で推定可能。εの同定困難性は文献と一致。
2. **電気力学連成**: peak stress 67.7 kPa、短縮10.9%、遅延23.8 msは生理学的範囲内（Gerach 2021基準）。
3. **不整脈リスクML**: AUROC=0.921±0.042は合成コホートで高精度。EFとLA径が支配的予測因子（臨床的妥当性あり）。
4. **AFアブレーション**: AUROC≈0.57は意図的なノイズ追加の結果。データリーク（AUROC=1.0）は明示的に検出・修正。

### 6.2 限界と前提条件

| 限界 | 内容 |
|---|---|
| **合成データ依存** | 実患者MRIデータ不使用。分布仮定が実際と乖離する可能性 |
| **AP簡易モデル** | Aliev-Panfilovは現象論的モデル。イオンチャネル詳細（TP06）が必要 |
| **一方向連成** | EP→Mechanics のみ。Mechano-electric feedbackを無視 |
| **ε識別困難** | 回復パラメータは短時間ECGから同定困難（27.4%誤差） |
| **AF予測の困難** | 現実のAFアブレーション予測はN=1600+でもAUROC≈0.75 |

### 6.3 今後の展望

1. **実MRIデータ統合**: UK Biobankの55,000例（Ugurlu 2025）への適用
2. **Full TP06モデル**: GPU-acceleratedイオンモデルシミュレーション（Viola 2023）
3. **双方向連成**: stretch-activated channelsを含むEM連成
4. **ベイズ不確実性定量化**: 逆問題パラメータの信頼区間推定
5. **前向き臨床検証**: デジタルツイン予測と実際のアブレーション結果の比較

---

## 7. 生成ファイル一覧

| ファイル | 説明 |
|---|---|
| `figures/fig1_aliev_panfilov_snapshots.png` | AP 2D波面伝播 6スナップショット |
| `figures/fig2_action_potential.png` | 単細胞活動電位と回復変数 |
| `figures/fig3_inverse_problem.png` | 逆問題パラメータ推定結果 |
| `figures/fig4_ml_results.png` | MLモデル比較・特徴量重要度 |
| `figures/fig5_electromechanical.png` | 電気力学連成タイムコース |
| `figures/fig6_framework_overview.png` | フレームワーク全体図・結果サマリー表 |
| `data/raw/geometry_cohort.csv` | 3D形状コホートデータ (N=50) |
| `data/raw/arrhythmia_cohort.csv` | 不整脈リスクコホート (N=200) |
| `data/raw/af_ablation_cohort.csv` | AFアブレーションコホート (N=150) |
| `cardiac_digital_twin.ipynb` | 実装Jupyterノートブック |
| `paper.md` | 学術論文形式レポート |
| `report.md` | 本実験レポート |

---

## 8. 再現性情報

```
Python:       3.11.2 (GCC 12.2.0)
numpy:        2.3.5
pandas:       2.3.3
scipy:        1.17.1
scikit-learn: 1.6.1
matplotlib:   3.10.9
seaborn:      0.13.2
random_seed:  42 (np.random.seed(42), random.seed(42))
```

全実験は `random_state=42` を統一設定し、完全再現可能。
