# 空間トランスクリプトミクス高度解析パイプライン — 設計レポート

**DRAFT — NOT FOR DISTRIBUTION**

| 項目 | 内容 |
|------|------|
| 作成日 | 2026-05-23 |
| バージョン | 1.0.0 |
| プラットフォーム | Visium / MERFISH / Slide-seq |
| 主要ライブラリ | Squidpy, SpatialDE, cell2location, LIANA, scanpy |

---

## 1. 実験目的と背景

空間トランスクリプトミクスは、組織内における遺伝子発現の空間分布を網羅的に計測する技術である。10x Visium（スポットベース、約 55 μm 解像度）や MERFISH（単一細胞解像度、数百〜数千遺伝子）などのプラットフォームが広く普及しているが、得られるデータから生物学的知見を抽出するためには、複数の解析ステップを統合した体系的パイプラインが不可欠である。

本パイプラインは以下の 6 つの解析課題を統合的に扱うフレームワークとして設計した：

1. **スポットデコンボリューション** — Visium の各スポットに含まれる細胞タイプ組成の推定
2. **空間的遺伝子発現パターン検出** — 空間的に構造化された発現を持つ遺伝子（SVG）の同定
3. **細胞間コミュニケーション推定** — リガンド–受容体ペアによる細胞間シグナル伝達の空間解析
4. **組織微小環境ニッチ同定** — 細胞タイプ組成と空間構造に基づく組織ニッチの分類
5. **3D 空間再構成** — 連続切片の位置合わせと三次元統合
6. **腫瘍免疫微小環境（TIME）ケーススタディ** — 免疫浸潤勾配・疲弊マーカー・免疫チェックポイントの空間解析

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 データ前処理（Module 0）

| ステップ | 手法 | パラメータ |
|----------|------|-----------|
| QC フィルタリング | ミトコンドリア比率・遺伝子数閾値 | max_pct_mito=20%, min_genes=200 |
| 正規化 | Log-normalization (target_sum=10⁴) | scanpy.pp.normalize_total + log1p |
| 次元削減 | HVG 選択 → PCA → UMAP | n_hvg=3000, n_pcs=30 |
| クラスタリング | Leiden algorithm | resolution=0.8 |
| 空間グラフ構築 | Delaunay/KNN | Squidpy spatial_neighbors (n_neighs=6) |

**注意**: QC 閾値はデータセット固有であり、ユニバーサルなカットオフの適用は推奨しない。各データセットの QC メトリクス分布を確認した上で調整すること。

### 2.2 スポットデコンボリューション（Module 1: cell2location）

```
scRNA-seq reference → NB regression → cell-type signatures
                                            ↓
Spatial data → cell2location model → posterior cell-type abundances
```

**アルゴリズム概要:**

- **Stage 1（参照モデル）**: scRNA-seq アトラスに対して Negative Binomial (NB) 回帰を適用し、各細胞タイプの遺伝子発現シグネチャ（平均発現プロファイル）を推定。`RegressionModel` を使用し、バッチ効果を共変量として補正。
- **Stage 2（空間マッピング）**: 推定シグネチャを事前分布として用い、空間データの各スポットにおける細胞タイプ存在量を変分推論で推定。事後分布の 5th パーセンタイル (`q05_cell_abundance_w_sf`) を下流解析に使用。

| パラメータ | 値 | 説明 |
|-----------|-----|------|
| N_cells_per_location | 8 | Visium スポットあたりの期待細胞数 |
| detection_alpha | 20 | 検出感度の正則化パラメータ |
| n_epochs_ref | 250 | 参照モデルの学習エポック数 |
| n_epochs_spatial | 30,000 | 空間モデルの学習エポック数 |

### 2.3 空間的遺伝子発現パターン検出（Module 2）

2 つの相補的手法を併用し、コンセンサスにより堅牢な SVG リストを構築する。

#### Moran's I（Squidpy）
- 空間的自己相関の古典的指標。I ∈ [-1, 1] で、正の値は空間クラスタリングを示す。
- 並べ替え検定（n_perms=1000）により統計的有意性を評価。
- Benjamini–Hochberg 法で多重検定補正。

#### SpatialDE
- ガウス過程ベースの手法。遺伝子発現の空間パターンを、Squared-Exponential カーネルでモデル化。
- **FSV**（Fraction of Spatial Variance）: 発現分散のうち空間構造に起因する割合。
- **length scale (l)**: 空間パターンの特徴的スケール（短い = 局所的、長い = 大域的）。

#### コンセンサス戦略
```
gene ∈ consensus_SVG  ⟺  (Moran's I q-value < 0.05) ∧ (SpatialDE q-value < 0.05)
```

### 2.4 細胞間コミュニケーション推定（Module 3: LIANA + Squidpy）

**LIANA (LIgand-receptor ANalysis frAmework)**:
- CellPhoneDB, CellChat, NATMI, Connectome, logFC, SingleCellSignalR の 6 手法の結果を統合ランキング（rank aggregation）で統合。
- `magnitude_rank`: 相互作用の強度ランク、`specificity_rank`: 細胞タイプ対特異性ランク。
- LR データベース: consensus（複数 DB の共通セット）を使用。

**Squidpy 空間コンテキスト**:
- `nhood_enrichment`: 実際の空間隣接関係から、特定の細胞タイプ対が期待以上に共局在するかを検定。
- 空間的リガンド–受容体共発現: 二変量 Moran's I でリガンドと受容体の空間的協調を評価。

### 2.5 組織微小環境ニッチ同定（Module 4）

```
各スポット → 空間的近傍の細胞タイプ組成プロファイル
          → (+ cell2location デコンボリューション結果)
          → StandardScaler → PCA → Leiden クラスタリング
          → ニッチラベル割り当て
```

- **Neighborhood profile**: 各スポットの空間的 k-nearest neighbors（k=15）における細胞タイプの出現比率を計算。
- **特徴量統合**: デコンボリューション結果が利用可能な場合、neighborhood profile と結合。
- **クラスタリング**: Leiden アルゴリズム（resolution=0.8）または KMeans（シルエットスコアで k 自動選択）。
- **ニッチ特徴付け**: 各ニッチの平均細胞タイプ組成を算出し、dominant cell type を同定。

### 2.6 3D 空間再構成（Module 5）

**アルゴリズム: Iterative Closest Point (ICP)**

```
Section i (source) ──→ ICP alignment ──→ Section i-1 (target)
                            ↓
                    Transformation matrix T_i
                            ↓
Aligned 2D coordinates + z-spacing → 3D coordinates (x, y, z)
```

1. **ペアワイズ ICP**: 連続する 2 切片間で、最近傍対応点を求め、SVD ベースの剛体変換（回転 + 並進）を推定。収束まで反復（max_iter=200, tol=1e-6）。
2. **逐次位置合わせ**: Section 1 を参照フレームとし、Section 2→1, 3→2, ... の順に累積変換を適用。
3. **Z 軸割り当て**: 切片間距離（z_spacing=10 μm）を等間隔に割り当て。
4. **AnnData 統合**: 全切片を `ad.concat` で結合し、`obsm['spatial_3d']` に 3D 座標を格納。

### 2.7 腫瘍免疫微小環境ケーススタディ（Module 6）

#### 腫瘍境界定義
- 腫瘍スポットと非腫瘍スポット間のユークリッド距離を計算。
- 4 ゾーンに分類: `tumor_core`, `tumor_edge`（境界 200 μm 以内）, `stroma_near`, `stroma_far`。

#### 免疫浸潤勾配
- 腫瘍境界からの距離ビン（0–50, 50–100, 100–200, 200–500 μm）ごとに免疫細胞密度を定量。
- CD8+ T 細胞、CD4+ T 細胞、マクロファージの密度を算出。

#### T 細胞疲弊スコアリング
- 疲弊マーカー遺伝子セット: PDCD1 (PD-1), LAG3, HAVCR2 (TIM-3), TIGIT, TOX
- `scanpy.tl.score_genes` によるモジュールスコア計算。

#### 免疫チェックポイントリガンド
- CD274 (PD-L1), PDCD1LG2 (PD-L2), LGALS9 のモジュールスコア。
- 腫瘍ゾーンとの空間的対応を評価。

#### 相互作用ホットスポット
- 腫瘍スポットの空間的近傍に含まれる免疫細胞数を計数。
- 腫瘍–免疫接触頻度の高い領域をホットスポットとして同定。

---

## 3. 主要な結果と数値

> **注記**: 本レポートはパイプライン設計段階のものであり、以下は期待される出力の構造と指標を示す。実データでの実行後に具体的数値が記載される。

### 3.1 期待される出力指標

| 解析モジュール | 主要出力指標 | 期待値の目安 |
|---------------|-------------|-------------|
| M0: QC | 残存スポット/細胞数 | 入力の 80–95% |
| M1: デコンボリューション | 推定細胞タイプ数 | 10–30 types |
| M2: SVG 検出 | コンセンサス SVG 数 | 200–2,000 genes |
| M3: 細胞間通信 | 有意な LR ペア数 | 50–500 pairs |
| M4: ニッチ同定 | ニッチ数 | 3–10 niches |
| M5: 3D 再構成 | ICP アラインメント RMSE | < 50 μm |
| M6: TIME 解析 | 免疫浸潤勾配・疲弊スコア分布 | ゾーン依存 |

### 3.2 品質管理チェックリスト

- [x] パイプライン全モジュールのコード実装完了
- [x] 設定ファイル（config.yaml）によるパラメータ管理
- [x] 多重検定補正（Benjamini–Hochberg）の適用
- [x] ランダムシード固定（seed=42）による再現性確保
- [x] 全図表を `figures/` に保存（英語ラベル）
- [x] 全数値結果を `results/` に CSV 保存
- [x] 実行ログを `logs/process-log.jsonl` に記録

---

## 4. 考察と今後の展望

### 4.1 手法選択の根拠

| 選択 | 根拠 | 代替手法 |
|------|------|---------|
| cell2location | ベイズ推論ベースで不確実性を定量化。Visium のスポット解像度に最適化 | RCTD, stereoscope, Tangram |
| SpatialDE + Moran's I | モデルベース（GP）と統計的検定の二重確認で偽陽性を低減 | SPARK, SOMDE, nnSVG |
| LIANA | 複数手法の統合ランキングにより単一手法のバイアスを軽減 | CellChat, NicheNet, Commot |
| ICP | 剛体変換で十分な場合に計算コストが低い | PASTE, STAligner, GraphST |

### 4.2 制限事項

1. **cell2location の参照データ依存性**: デコンボリューションの精度は scRNA-seq アトラスの品質と細胞タイプアノテーションに強く依存する。組織・疾患特異的な参照データの使用を推奨。
2. **SpatialDE の計算コスト**: 大規模データ（>10,000 スポット）では計算時間が増大する。SpatialDE2 や nnSVG への移行を検討。
3. **ICP の局所最適**: ICP は初期配置に依存し、大きな回転・変形がある場合に局所最適に陥る可能性がある。PASTE（最適輸送ベース）との比較を推奨。
4. **LIANA の空間情報統合**: LIANA 自体は空間情報を直接利用しない。Squidpy の空間共局在解析と事後的に統合する設計としたが、Commot や SpaTalk のような空間ネイティブ手法との比較が望まれる。

### 4.3 今後の展望

- **nnSVG への対応**: SVG 検出に nearest-neighbor Gaussian process を導入し、大規模データへのスケーラビリティを向上。
- **PASTE 統合**: 最適輸送ベースの切片間アラインメントにより、非剛体変形に対応した 3D 再構成を実現。
- **マルチモーダル統合**: MERFISH + Visium、あるいは空間プロテオミクス（CODEX）との統合解析フレームワークの設計。
- **GPU 最適化**: cell2location / scvi-tools の GPU 学習によるスループット向上。
- **インタラクティブ可視化**: Napari / Vitessce によるブラウザベースの空間データ探索環境の構築。

---

## 5. パイプライン・アーキテクチャ図

```
┌──────────────────────────────────────────────────────────────────────┐
│                      run_pipeline.py (Orchestrator)                  │
│                                                                      │
│  config.yaml ─→ M0 ─→ M1 ─→ M2 ─→ M3 ─→ M4 ─→ M5 ─→ M6          │
│                  │     │     │     │     │     │     │               │
│                  ▼     ▼     ▼     ▼     ▼     ▼     ▼               │
│               figures/ results/ data/ logs/process-log.jsonl         │
└──────────────────────────────────────────────────────────────────────┘

M0: Data Loading & QC        → scanpy, squidpy
M1: Spot Deconvolution       → cell2location (scvi-tools)
M2: SVG Detection            → Moran's I (squidpy) + SpatialDE
M3: Cell Communication       → LIANA + squidpy nhood_enrichment
M4: Niche Identification     → Neighborhood profiling + Leiden
M5: 3D Reconstruction        → ICP alignment + z-stacking
M6: Tumor–Immune Case Study  → Boundary analysis + exhaustion scoring
```

---

## 6. 生成したファイル一覧

### コード・設定

| ファイル | 説明 |
|---------|------|
| `config.yaml` | パイプライン全体の設定ファイル（QC 閾値、モデルパラメータ等） |
| `requirements.txt` | Python 依存パッケージ一覧 |
| `run_pipeline.py` | メインオーケストレータスクリプト |
| `pipeline/__init__.py` | パッケージ初期化 |
| `pipeline/m00_data_loading.py` | データ読み込み・QC・正規化・クラスタリング |
| `pipeline/m01_deconvolution.py` | cell2location デコンボリューション |
| `pipeline/m02_spatial_patterns.py` | Moran's I + SpatialDE による SVG 検出 |
| `pipeline/m03_communication.py` | LIANA + Squidpy による細胞間通信解析 |
| `pipeline/m04_niche.py` | 組織ニッチ同定・特徴付け |
| `pipeline/m05_reconstruction_3d.py` | ICP ベース 3D 空間再構成 |
| `pipeline/m06_tumor_immune.py` | 腫瘍免疫微小環境解析 |

### 出力ディレクトリ（実行後に生成）

| ディレクトリ/ファイル | 内容 |
|---------------------|------|
| `figures/qc_violin.png` | QC メトリクスバイオリンプロット |
| `figures/deconvolution_map.png` | 細胞タイプ空間分布マップ |
| `figures/svg_expression_maps.png` | SVG 空間発現パターン |
| `figures/communication_network.png` | LR 相互作用ヒートマップ |
| `figures/nhood_enrichment.png` | 近傍エンリッチメントヒートマップ |
| `figures/niche_map.png` | 組織ニッチ空間マップ |
| `figures/niche_composition.png` | ニッチ別細胞タイプ組成 |
| `figures/3d_reconstruction.png` | 3D 再構成散布図 |
| `figures/tumor_zones.png` | 腫瘍ゾーン空間マップ |
| `figures/immune_gradient.png` | 免疫浸潤勾配プロット |
| `figures/exhaustion_spatial.png` | T 細胞疲弊スコア空間マップ |
| `figures/tumor_immune_landscape.png` | TIME 概観（4 パネル） |
| `results/cell_type_abundances.csv` | 細胞タイプ存在量マトリクス |
| `results/spatially_variable_genes.csv` | コンセンサス SVG テーブル |
| `results/ligand_receptor_results.csv` | 上位 LR ペアリスト |
| `results/niche_assignments.csv` | ニッチ割り当てラベル |
| `results/niche_summary.csv` | ニッチ別組成サマリ |
| `results/coordinates_3d.csv` | 3D 座標データ |
| `results/tumor_immune_report_*.csv` | TIME 解析結果一式 |
| `logs/process-log.jsonl` | 実行トレースログ |

---

## 7. 再現性情報

```yaml
random_seed: 42
python: ">=3.10"
key_dependencies:
  scanpy: ">=1.9"
  squidpy: ">=1.3"
  cell2location: ">=0.1.3"
  SpatialDE: ">=1.1"
  liana: ">=1.0"
  scvi-tools: ">=1.0"
```

全ての乱数生成器（numpy, random, torch）に対してシード 42 を設定。依存パッケージのバージョンは `requirements.txt` で固定。

---

*本レポートは Co-Scientist spatial-transcriptomics スキルにより自動生成されました。*
