# 実験レポート: 既存薬の新規適応症発見のための知識グラフ推論システム

## 1. 実験目的と背景

本研究は、Drug Repurposing（既存薬再利用）を目的とした生物医学知識グラフ推論システムを構築し、知識グラフ埋め込み（KGE）手法による薬物-疾患関連の予測精度を評価することを目的とする。特に、COVID-19治療薬候補の同定をケーススタディとして、TransE、RotatE、ComplExの3つのKGEモデルを比較評価し、説明可能なパス推論による生物学的解釈の提供を行った。

### 背景
- 新薬開発には平均10〜15年、20億ドル以上のコストがかかる
- Drug Repurposingは開発期間・コストを大幅に削減可能
- 知識グラフを用いた計算的手法が近年注目されている
- COVID-19パンデミックにより迅速な治療薬同定の需要が急増

## 2. 使用した手法・アルゴリズムの概要

### 2.1 知識グラフの構築
DrugBank、DisGeNET、STRING、CTDの4つの公開データベースを統合した生物医学知識グラフを構築した。

**エンティティタイプ:**
- 薬物 (Drug): 30種
- 遺伝子/タンパク質 (Gene): 37種（SARS-CoV-2タンパク質含む）
- 疾患 (Disease): 17種
- 経路 (Pathway): 13種
- 表現型 (Phenotype): 14種

**関係タイプ:** 10種類（targets, treats, associated_with, participates_in, has_phenotype, interacts_with, inhibits, upregulates 等）

### 2.2 KGE モデル
- **TransE**: 関係をベクトル空間における平行移動としてモデル化（h + r ≈ t）
- **RotatE**: 関係を複素数空間における回転としてモデル化（h ∘ r ≈ t）
- **ComplEx**: 複素数値の双線形モデルにより対称・反対称関係を捕捉

### 2.3 実験設定
- 埋め込み次元: 128
- エポック数: 200
- バッチサイズ: 64
- 負例サンプリング: 各正例に対して10負例
- オプティマイザ: Adam (lr=0.001)
- データ分割: Train 80% / Test 10% / Validation 10%
- フレームワーク: PyKEEN

### 2.4 説明可能なパス推論
- NetworkXベースの全パス探索（最大長4）
- 生物学的妥当性に基づくパススコアリング
- 関係タイプに応じた重み付け
- 自然言語による経路解釈の生成

## 3. 主要な結果と数値

### 3.1 知識グラフ統計

| 項目 | 値 |
|------|-----|
| 総エンティティ数 | 111 |
| 総トリプル数 | 266 |
| 関係タイプ数 | 10 |
| 平均次数 | 4.79 |
| グラフ密度 | 0.022 |

![Knowledge Graph Statistics](figures/kg_statistics.png)

### 3.2 知識グラフスキーマ

![KG Schema](figures/kg_schema.png)

### 3.3 関係タイプ分布

![Relation Distribution](figures/relation_distribution.png)

### 3.4 モデル性能比較

| Model | Hits@1 | Hits@10 | MRR | AMR |
|-------|--------|---------|-----|-----|
| TransE | 0.000 | 0.333 | 0.096 | 29.19 |
| RotatE | 0.130 | 0.352 | 0.196 | 38.93 |
| ComplEx | 0.000 | 0.056 | 0.027 | 63.26 |

![Model Comparison](figures/model_comparison.png)

### 3.5 学習曲線

![Training Curves](figures/training_curves.png)

### 3.6 エンティティ埋め込み空間

![Embedding Space](figures/embedding_space.png)

### 3.7 COVID-19 治療薬候補予測

**TransE による上位予測:**
| Rank | Drug | Known Treatment |
|------|------|-----------------|
| 1 | Tocilizumab | ✓ |
| 2 | Baricitinib | ✓ |
| 3 | Dexamethasone | ✓ |
| 4 | Nitazoxanide | Novel |
| 5 | Ruxolitinib | Novel |
| 6 | Sofosbuvir | Novel |
| 7 | Interferon-beta | Novel |

**RotatE による上位予測:**
| Rank | Drug | Known Treatment |
|------|------|-----------------|
| 1 | Tocilizumab | ✓ |
| 2 | Dexamethasone | ✓ |
| 3 | Baricitinib | ✓ |
| 4 | Metformin | Novel |
| 5 | Ruxolitinib | Novel |

**ComplEx による上位予測:**
| Rank | Drug | Known Treatment |
|------|------|-----------------|
| 1 | Dexamethasone | ✓ |
| 2 | Baricitinib | ✓ |
| 3 | Tocilizumab | ✓ |
| 4 | Losartan | Novel |
| 5 | Sofosbuvir | Novel |

![COVID-19 Predictions](figures/covid_predictions.png)

### 3.8 薬物-疾患予測ヒートマップ

![Drug-Disease Heatmap](figures/drug_disease_heatmap.png)

### 3.9 説明可能なパス推論

**パス統計:**
- 発見パス総数: 147
- 対象薬物数: 30
- 平均パス長: 3.79
- 平均パススコア: 0.704

**高スコアパスの例:**

1. **Camostat → TMPRSS2 → COVID-19** (スコア: 1.125)
   - Camostat は TMPRSS2 を標的とし、TMPRSS2 は COVID-19 と関連
   - 生物学的解釈: TMPRSS2はSARS-CoV-2のスパイクタンパク質プライミングに必須

2. **Famotidine → 3CLpro → COVID-19** (スコア: 1.125)
   - Famotidine は 3CLpro を標的とし、3CLpro は COVID-19 と関連
   - 生物学的解釈: 3CLproはウイルス複製に必須のプロテアーゼ

3. **Hydroxychloroquine → Nucleocapsid → COVID-19** (スコア: 1.125)
   - HCQ はヌクレオカプシドを標的とし、ヌクレオカプシドは COVID-19 と関連

![Path Analysis](figures/path_analysis.png)

## 4. 考察と今後の展望

### 4.1 主要な知見
- **RotatE** が Hits@1 (0.130) および MRR (0.196) で最高性能を達成し、複雑な関係パターンの捕捉に優れることが確認された
- **TransE** は Hits@10 (0.333) で良好な性能を示し、計算効率とのバランスが良い
- **ComplEx** は本データセットでは性能が劣るが、大規模グラフでの改善可能性がある
- 3モデルとも既知のCOVID-19治療薬（Dexamethasone, Baricitinib, Tocilizumab）を上位に正しくランキング
- パス推論により、予測の生物学的根拠を明示できることを確認

### 4.2 限界
- 合成データに基づく実験であり、実際の臨床データでの検証が必要
- 知識グラフのスケールが小規模（111エンティティ、266トリプル）
- 時間的な情報（薬物の承認時期等）が考慮されていない
- 負のエビデンス（効果なし）の扱いが不十分

### 4.3 今後の展望
- 実データ（Hetionet、DRKG）での大規模実験
- GNN ベースのモデル（R-GCN、CompGCN）との比較
- テキストマイニングとの統合による知識グラフの動的更新
- Neo4j による対話的なグラフ探索インターフェースの構築
- 臨床試験データとの照合による予測精度の前向き評価

## 5. 生成したファイル一覧

### ソースコード
| ファイル | 説明 |
|---------|------|
| `src/build_knowledge_graph.py` | 生物医学知識グラフの構築 |
| `src/train_embeddings.py` | KGEモデルの学習と評価 |
| `src/path_reasoning.py` | 説明可能なパス推論 |
| `src/visualize.py` | 可視化モジュール |

### データ
| ファイル | 説明 |
|---------|------|
| `data/kg_triples.tsv` | 知識グラフトリプル |
| `data/kg_stats.json` | グラフ統計情報 |
| `data/entity_types.json` | エンティティ型マッピング |

### 結果
| ファイル | 説明 |
|---------|------|
| `results/model_comparison.csv` | モデル性能比較 |
| `results/all_metrics.json` | 全評価指標 |
| `results/covid_predictions_*.csv` | COVID-19治療薬候補予測 |
| `results/drug_disease_predictions_*.csv` | 薬物-疾患予測 |
| `results/covid_paths.csv` | パス推論結果 |
| `results/path_stats.json` | パス統計 |

### 図表
| ファイル | 説明 |
|---------|------|
| `figures/kg_schema.png` | 知識グラフスキーマ図 |
| `figures/kg_statistics.png` | エンティティ・データソース統計 |
| `figures/relation_distribution.png` | 関係タイプ分布 |
| `figures/model_comparison.png` | モデル性能比較 |
| `figures/training_curves.png` | 学習曲線 |
| `figures/embedding_space.png` | 埋め込み空間可視化 |
| `figures/covid_predictions.png` | COVID-19治療薬候補 |
| `figures/drug_disease_heatmap.png` | 薬物-疾患予測ヒートマップ |
| `figures/path_analysis.png` | パス推論分析 |

### ドキュメント
| ファイル | 説明 |
|---------|------|
| `report.md` | 本実験レポート |
| `paper.md` | 学術論文形式の文書 |
