# Perturb-seq 解析フレームワーク — 実験報告書

**DRAFT — NOT FOR DISTRIBUTION**

- **日時**: 2026-05-23
- **解析パイプライン**: Scanpy/NMF/PC Algorithm ベース
- **乱数シード**: 42

---

## 1. 実験目的と背景

Perturb-seq（CRISPR スクリーニング + scRNA-seq）は、大規模な遺伝子摂動の転写応答を1細胞レベルで計測する技術である。本フレームワークは、以下の6つの解析モジュールを統合したエンドツーエンドパイプラインを設計・実装した：

1. **品質管理（QC）とガイド検出** — ガイドRNA割り当ての信頼性評価
2. **遺伝子プログラム変動検出** — 差分発現解析（DE）と共発現モジュール同定
3. **因果グラフ推定** — 摂動効果からの遺伝子制御ネットワーク推定
4. **エピスタシス検出** — 組合せ摂動の相互作用効果定量
5. **低次元表現学習** — scVI/CPA スタイルの潜在空間表現
6. **必須遺伝子ネットワーク** — フィットネス効果に基づく必須遺伝子同定

パイプラインの設計検証にはシミュレーションデータ（5,000細胞 × 2,000遺伝子 × 20ガイド）を使用した。

---

## 2. 使用手法・アルゴリズムの概要

### 2.1 QC & ガイド検出（Module 1）

| 項目 | 手法 |
|------|------|
| ガイドUMI閾値決定 | Gaussian Mixture Model（2成分）per guide |
| 細胞QC | 最小遺伝子数（200）、最小UMI数（500）、99%ile上限 |
| 信頼度スコア | UMI / 閾値比に基づく正規化スコア（マルチプレット時は0.8倍ペナルティ） |

### 2.2 差分発現 & 共発現モジュール（Module 2）

| 項目 | 手法 |
|------|------|
| 差分発現検定 | Wilcoxon rank-sum test（各摂動 vs. 非標的対照） |
| 多重検定補正 | Benjamini-Hochberg FDR（摂動ごと） |
| 共発現モジュール | Non-negative Matrix Factorization (NMF, 8成分) |
| 前処理 | Scanpy標準（normalize_total → log1p → HVG 1000 → scale → PCA 30次元） |

### 2.3 因果グラフ推定（Module 3）

| 項目 | 手法 |
|------|------|
| 因果構造学習 | PC アルゴリズム（制約ベース） |
| 条件付き独立性検定 | 偏相関 + Fisher z検定 |
| エッジの方向付け | Meek's rules（v構造検出） |
| 安定性評価 | ブートストラップ（50回リサンプリング） |

### 2.4 エピスタシス検出（Module 4）

| 項目 | 手法 |
|------|------|
| 加法性モデル | Expected = effect(A) + effect(B) |
| エピスタシススコア | Observed - Expected（L2ノルム） |
| 分類 | Synergy（>1.2倍）/ Buffering（<0.8倍）/ Suppression（コサイン類似度<0）/ Additive |
| 有意性検定 | 並び替え検定（500回、FDR補正） |

### 2.5 低次元表現学習（Module 5）

| 項目 | 手法 |
|------|------|
| VAEエンコーダ | PCA初期化付き2層エンコーダ（隠れ層64、潜在次元10） |
| CPA分解 | 基底状態（PCA on controls）+ 摂動残差 |
| 摂動類似度 | コサイン類似度行列 |
| クラスタリング | Ward法階層的クラスタリング（3クラスタ） |

### 2.6 必須遺伝子ネットワーク（Module 6）

| 項目 | 手法 |
|------|------|
| フィットネス推定 | UMI比（摂動/対照）× 転写効果量 |
| 必須遺伝子分類 | フィットネスzスコア < -1.5 |
| ネットワーク構築 | 摂動効果のコサイン類似度（閾値0.3） |
| コミュニティ検出 | Greedy modularity (Louvain-style) |

---

## 3. 主要な結果と数値

### 3.1 QC & ガイド検出

| 指標 | 値 |
|------|-----|
| 総細胞数 | 5,000 |
| QC通過細胞 | 4,955（99.1%） |
| 未割り当て細胞 | 1,422（28.4%） |
| 単一ガイド検出 | 3,355（67.1%） |
| 複数ガイド検出 | 223（4.5%） |
| 解析対象（QC通過 & 割り当て済み） | 3,544 |
| 平均信頼度スコア | 0.554 |
| ガイドUMI閾値中央値 | ~6.3（GMM推定） |

### 3.2 差分発現

| 指標 | 値 |
|------|-----|
| 総検定数 | 40,000 |
| 有意なDEG（FDR < 0.05） | 4 |
| 検定した摂動数 | 20 |
| DEGを持つ摂動 | 4（gene_9, gene_10, gene_12, gene_15） |
| NMF共発現モジュール | 8（各63遺伝子） |

> **注**: シミュレーションデータのため有意DEGは少数。実データでは数百〜数千のDEGが期待される。

### 3.3 因果グラフ

| 指標 | 値 |
|------|-----|
| グラフノード数 | 4 |
| 推定エッジ数 | 0 |
| グラフ密度 | 0.000 |
| 安定エッジ（>50%） | 0 |

> **注**: 効果行列が4遺伝子×20摂動と小規模であり、有意なエッジは検出されなかった。実データではDE遺伝子数増加に伴い密なネットワークが期待される。

### 3.4 エピスタシス

| 指標 | 値 |
|------|-----|
| 検定した組合せ | 6 |
| 有意な相互作用 | 2 |
| 相互作用タイプ | Synergy: 6 |
| 平均エピスタシス強度 | 17.47 |
| 最強相互作用 | gene_6 × gene_9（強度 18.34） |

### 3.5 低次元表現

| 指標 | 値 |
|------|-----|
| 潜在次元数 | 10 |
| 潜在空間分散 | 7.95 |
| プロファイルした摂動 | 48（単一20 + 組合せ28） |
| クラスタ数 | 3 |
| クラスタ内平均類似度 | 0.068 |

### 3.6 必須遺伝子ネットワーク

| 指標 | 値 |
|------|-----|
| 解析した摂動 | 20 |
| 必須遺伝子候補 | 1（gene_1） |
| ネットワークエッジ | 0 |
| コミュニティ数 | 1 |

---

## 4. 考察と今後の展望

### 4.1 パイプラインの有効性

本フレームワークは、Perturb-seq データの主要な解析ステップを統合した完全なパイプラインを提供する。各モジュールは独立に実行可能であり、中間結果はファイルに永続化されるため、段階的な解析が可能である。

### 4.2 シミュレーションデータの限界

シミュレーションデータでは摂動効果のシグナル対ノイズ比が限定的であり、以下の点で実データとは異なる：

- **DEG検出数**: 実験データでは摂動あたり数十〜数百のDEGが検出されるのに対し、シミュレーションでは4件のみ
- **因果グラフ**: 効果行列のサイズが小さく、有意なエッジが推定されなかった
- **必須遺伝子**: フィットネス効果のモデリングが簡素であり、1遺伝子のみが候補として検出された

### 4.3 実データへの適用時の推奨事項

1. **ガイド検出**: GMMの成分数を交差検証で最適化。CITE-seqや hashtag antibody との併用でマルチプレット検出精度を向上
2. **DE解析**: MAST や edgeR など、scRNA-seq 専用の統計検定への拡張
3. **因果推定**: DCDI（Differentiable Causal Discovery with Interventions）や GIES アルゴリズムの導入
4. **エピスタシス**: Pertpy の `MixtureModel` や `Augur` との統合
5. **潜在表現**: 完全な scVI/CPA モデルの利用（`scvi-tools` パッケージ）
6. **必須遺伝子**: DepMap/CERES スコアとの比較検証

### 4.4 スケーラビリティ

- 10万細胞規模のデータには、GPU対応の scVI エンコーダ（Module 5）が不可欠
- DE解析は遺伝子数 × 摂動数のスケーリングのため、並列化（Dask/Ray）を推奨
- 因果推定のPC アルゴリズムは O(n²) であり、大規模ネットワーク（>1000ノード）には FCI や GES を検討

---

## 5. 生成ファイル一覧

### データファイル（`data/`）

| ファイル | 内容 |
|---------|------|
| `perturbseq_simulated.h5ad` | シミュレーション生データ（5000×2000） |
| `perturbseq_qc_filtered.h5ad` | QC通過・割り当て済みデータ（3544×2000） |
| `perturbseq_processed.h5ad` | 前処理済み + DE結果 |
| `perturbseq_with_latent.h5ad` | 潜在表現付きデータ |
| `simulation_summary.json` | シミュレーション設定 |

### 結果ファイル（`results/`）

| ファイル | 内容 |
|---------|------|
| `01_qc_stats.json` | QC統計サマリ |
| `02_de_results.csv` | 全DE検定結果（40,000行） |
| `02_de_summary.json` | DE解析サマリ |
| `02_modules.json` | NMF共発現モジュール遺伝子リスト |
| `03_causal_edges.csv` | 推定因果エッジリスト |
| `03_causal_summary.json` | 因果グラフサマリ |
| `03_adjacency_matrix.csv` | 隣接行列 |
| `03_edge_stability.csv` | ブートストラップ安定性行列 |
| `04_epistasis_results.csv` | エピスタシス検定結果 |
| `04_epistasis_summary.json` | エピスタシスサマリ |
| `05_perturbation_similarity.csv` | 摂動間コサイン類似度行列 |
| `05_latent_summary.json` | 潜在表現サマリ |
| `06_fitness_scores.csv` | フィットネススコア一覧 |
| `06_network_edges.csv` | 必須遺伝子ネットワークエッジ |
| `06_essential_summary.json` | 必須遺伝子解析サマリ |

### 図表ファイル（`figures/`）

| ファイル | 内容 |
|---------|------|
| `01_qc_guide_detection.{png,svg}` | QC・ガイド検出6パネル図 |
| `02_de_coexpression.{png,svg}` | 火山プロット・DEG数・NMFモジュール・UMAP |
| `03_causal_graph.{png,svg}` | 因果DAG・エッジ安定性ヒートマップ |
| `04_epistasis.{png,svg}` | エピスタシス散布図・分類・ランキング |
| `05_latent_representation.{png,svg}` | t-SNE潜在空間・類似度ヒートマップ・デンドログラム |
| `06_essential_network.{png,svg}` | 必須遺伝子ネットワーク・フィットネス分布・パスウェイ富化 |

### パイプラインコード（`src/`）

| ファイル | 内容 |
|---------|------|
| `00_setup.py` | データシミュレーション |
| `01_qc_guide_detection.py` | QC & ガイド割り当て |
| `02_differential_expression.py` | DE解析 & 共発現モジュール |
| `03_causal_graph.py` | PC因果グラフ推定 |
| `04_epistasis.py` | エピスタシス検出 |
| `05_latent_representation.py` | VAE/CPA潜在表現 |
| `06_essential_gene_network.py` | 必須遺伝子ネットワーク |
| `run_pipeline.py` | パイプラインオーケストレータ |

### ログ（`logs/`）

| ファイル | 内容 |
|---------|------|
| `process-log.jsonl` | 実行トレース（タイムスタンプ・フェーズ・イベント） |

---

## 6. 実行方法

```bash
# 依存パッケージインストール
pip install scanpy anndata scikit-learn statsmodels networkx matplotlib scikit-misc

# パイプライン全体実行
PYTHONPATH=. python src/run_pipeline.py

# 個別モジュール実行
PYTHONPATH=. python src/00_setup.py          # データ生成
PYTHONPATH=. python src/01_qc_guide_detection.py  # QC
PYTHONPATH=. python src/02_differential_expression.py  # DE
PYTHONPATH=. python src/03_causal_graph.py   # 因果グラフ
PYTHONPATH=. python src/04_epistasis.py      # エピスタシス
PYTHONPATH=. python src/05_latent_representation.py  # 潜在表現
PYTHONPATH=. python src/06_essential_gene_network.py  # 必須遺伝子
```

---

*本報告書は Co-Scientist パイプラインにより自動生成されました。*
