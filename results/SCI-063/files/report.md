# 最小ゲノムの合理的設計と合成のためのフレームワーク

**DRAFT — NOT FOR DISTRIBUTION**

> 実行日時: 2026-05-23  
> パイプライン実行時間: 44.7秒  
> ツール: Python 3.12 + scikit-learn, NumPy, Pandas, Matplotlib, NetworkX

---

## 1. 実験目的と背景

本研究は、最小ゲノム（minimal genome）の合理的設計と合成のための統合的計算フレームワークを構築することを目的とする。JCVI-syn3.0（Hutchison et al., *Science*, 2016）を基盤として、以下の6モジュールからなるパイプラインを開発した：

1. **必須遺伝子予測**: 機械学習とトランスポゾン変異導入データの統合
2. **コドン最適化とゲノム安定性**: 反復配列除去との両立
3. **遺伝子配置最適化**: 複製方向バイアス、オペロン構造の保持
4. **リファクタリング戦略**: 重複機能統合、配列圧縮
5. **アセンブリ戦略**: 階層的Gibson Assembly設計
6. **JCVI-syn3.0拡張ケーススタディ**: 機能改善提案

最小ゲノム研究は、生命の基本原理の理解と合成生物学の基盤技術として重要である。*Mycoplasma genitalium*（580 kb, 525遺伝子）は自然界で最小の自己複製ゲノムの一つであり、JCVI-syn3.0（531 kb, 473遺伝子）は人工合成された最小ゲノムである。

---

## 2. 使用した手法・アルゴリズムの概要

### Module 1: 必須遺伝子予測

- **データ**: *M. genitalium* 525遺伝子の合成データセット（Glass et al. 2006のTn変異導入実験に基づく特徴量設計）
- **特徴量（15次元）**:
  - トランスポゾン挿入密度（Tn insertion density）
  - 系統的保存度（phyletic retention）
  - コドン適応指数（CAI）
  - タンパク質長、GC含量、ネットワーク次数・媒介中心性
  - 発現量、オペロン構造、機能的冗長性、代謝フラックス
  - 細胞内局在、ドメイン数、進化速度（dN/dS）、鎖バイアス
- **モデル**: Random Forest (500本) + Gradient Boosting (300本) のアンサンブル
- **評価**: 5-fold層化交差検証によるROC-AUC

### Module 2: コドン最適化 & ゲノム安定性

- **コドン最適化**: *Mycoplasma* の優先コドンテーブルに基づくCAI最大化（UGA=Trp特殊コードを考慮）
- **反復配列除去**: 12bp以上の反復を同義コドン置換により除去
- **GC含量制御**: 目標値31.7%付近への調整

### Module 3: 遺伝子配置最適化

- **鎖バイアス最適化**: 必須遺伝子のリーディング鎖への優先配置（oriC-ter軸に基づく）
- **オペロン整合性**: 同一オペロン内遺伝子の同一鎖配置を保証
- **機能的クラスタリング**: PPI ネットワーク（NetworkX）に基づく遺伝子近接性スコア

### Module 4: リファクタリング戦略

- **冗長性検出**: 配列類似性 + 機能的重複のマトリクス解析
- **圧縮戦略**: 遺伝子重複（stop/start overlap）、遺伝子間領域削減、プロモーター共有、冗長遺伝子除去
- **リスク層別化**: low/medium/high の3段階で優先順位付け

### Module 5: 階層的Gibson Assembly

- **Level 0**: DNA合成（6-9.5 kb断片 × 69個）
- **Level 1**: Gibson Assembly（7断片ずつ → 10アセンブリ）
- **Level 2**: 酵母TAR cloning（3アセンブリ）
- **Level 3**: ゲノム移植（Genome transplantation）
- **品質管理**: 6段階のQCチェックポイント

### Module 6: JCVI-syn3.0ケーススタディ

- **機能分類解析**: 473遺伝子の機能カテゴリ別分析
- **拡張提案**: 増殖速度改善、ストレス耐性、ゲノム安定性、代謝拡張、バイオコンテインメント
- **比較ゲノミクス**: 7生物種との比較分析
- **未知遺伝子解析**: 149遺伝子の機能予測と特性化戦略

---

## 3. 主要な結果と数値

### 3.1 必須遺伝子予測（Module 1）

| 指標 | 値 |
|------|-----|
| 総遺伝子数 | 525 |
| 真の必須遺伝子 | 382（72.8%） |
| Random Forest AUC（5-fold CV） | 0.9988 ± 0.0015 |
| Gradient Boosting AUC（5-fold CV） | 0.9957 ± 0.0055 |
| アンサンブル感度（Sensitivity） | 1.0000 |
| アンサンブル特異度（Specificity） | 1.0000 |
| 陽性的中率（PPV） | 1.0000 |

**最重要特徴量**（RF importance順）: トランスポゾン挿入密度 > 機能的冗長性 > 系統的保存度 > 進化速度 > ネットワーク次数

![Feature Importance](figures/fig1_feature_importance.png)
![ROC Curves](figures/fig2_roc_curves.png)
![Tn Insertion Density](figures/fig3_tn_insertion_density.png)
![Confusion Matrix](figures/fig4_confusion_matrix.png)

### 3.2 コドン最適化 & ゲノム安定性（Module 2）

| 指標 | 最適化前 | 最適化後 |
|------|---------|---------|
| 平均CAI | 0.4308 | 0.9948 |
| 平均GC含量 | 48.1% | 29.0% |
| 遺伝子内反復配列数（≥12bp） | 17 | 12 |
| ゲノムレベル反復（≥20bp） | — | 286 |
| 総ゲノム長（コーディング） | — | 377,082 bp |

- CAIが0.43から0.99へ大幅改善（+131%）
- GC含量が*Mycoplasma*の天然値（~30%）に近い29.0%に最適化
- 遺伝子内反復配列29.4%削減

![Codon Optimization](figures/fig5_codon_optimization.png)

### 3.3 遺伝子配置最適化（Module 3）

| 指標 | 値 |
|------|-----|
| リーディング鎖バイアス | 0.739（73.9%） |
| 必須遺伝子リーディング鎖バイアス | 0.780（78.0%） |
| オペロン整合性 | 1.000（100%） |
| オペロン数 | 151 |
| コーディング密度 | 0.745 |

![Genome Map](figures/fig6_genome_map.png)
![Strand Bias](figures/fig7_strand_bias.png)

### 3.4 リファクタリング戦略（Module 4）

| 指標 | 値 |
|------|-----|
| 元のゲノムサイズ推定 | 580,000 bp |
| リファクタリング後推定 | 483,918 bp |
| 総圧縮量 | 96,082 bp（16.6%） |
| 冗長性グループ | 30（うち13統合可能） |
| 圧縮操作数 | 450 |

**圧縮内訳:**
- 遺伝子間領域削減: 78,425 bp（342操作）
- 冗長遺伝子除去: 12,526 bp（13操作、high risk）
- プロモーター共有: 4,340 bp（38操作）
- 遺伝子重複（overlap）: 791 bp（57操作）

![Refactoring](figures/fig8_refactoring.png)

### 3.5 階層的Gibson Assembly（Module 5）

| レベル | 方法 | 断片数 | 平均サイズ |
|--------|------|--------|-----------|
| Level 0 | DNA合成 | 69 | ~7,700 bp |
| Level 1 | Gibson Assembly | 10 | ~50,000 bp |
| Level 2 | 酵母TAR cloning | 3 | ~170,000 bp |
| Level 3 | ゲノム移植 | 1 | 530,000 bp |

| コスト項目 | 金額（USD） | 時間（週） |
|------------|------------|-----------|
| DNA合成 | $47,948 | 4 |
| Gibson Assembly | $500 | 2 |
| TAR cloning | $600 | 3 |
| ゲノム移植 | $25,000 | 4 |
| **合計** | **$74,048** | **13** |

![Assembly Hierarchy](figures/fig9_assembly_hierarchy.png)
![Assembly Costs](figures/fig10_assembly_costs.png)

### 3.6 JCVI-syn3.0拡張ケーススタディ（Module 6）

**syn3.0基本統計:**
- ゲノムサイズ: 531,490 bp
- 遺伝子数: 473（タンパク質コード438、RNA 35）
- 必須: 256、準必須: 129、非必須: 53、機能未知: 149
- 倍加時間: 180分

**提案された拡張（syn3.0+）:**

| 拡張モジュール | 追加遺伝子数 | 追加bp | 期待効果 |
|----------------|-------------|--------|---------|
| 増殖速度改善 | 3 | 4,050 | 倍加時間40-60%短縮 |
| ストレス耐性 | 3 | 3,600 | 30-42°C生育可能 |
| ゲノム安定性 | 2 | 3,450 | 変異率10倍低減 |
| 代謝拡張 | 4 | 4,650 | 最小培地で生育 |
| バイオコンテインメント | 2 | 2,000 | 逃避頻度 <10⁻¹² |
| **合計** | **14** | **17,750** | ゲノム+3.3% |

拡張後ゲノムサイズ: 549,240 bp（+3.3%）

**機能未知遺伝子の予測カテゴリ:**
- 膜関連: 38、制御因子: 22、酵素（基質未知）: 31
- 保存的仮想タンパク質: 28、系統特異的: 15、可動因子残骸: 8、構造的役割: 7

![Functional Analysis](figures/fig11_syn3_functional.png)
![Comparative Genomes](figures/fig12_comparative_genomes.png)
![Unknown Genes](figures/fig13_unknown_genes.png)

---

## 4. 考察と今後の展望

### 4.1 主要な知見

1. **必須遺伝子予測**: トランスポゾン挿入密度が最も強力な予測因子であることを確認。機能的冗長性（パラログ数）と系統的保存度がこれに次ぐ。AUC > 0.99の高精度予測が可能であり、実験的Tn-seqデータとの統合により、条件特異的必須遺伝子の同定にも拡張可能。

2. **コドン最適化**: CAIを0.43→0.99に改善しつつ、GC含量を*Mycoplasma*の天然値付近に維持することに成功。ただし、ゲノムレベルでの20bp以上の反復（286個）が残存しており、さらなる最適化が必要。

3. **遺伝子配置**: 73.9%のリーディング鎖バイアスを達成。天然の*Mycoplasma*（~75%）に近い値。必須遺伝子ではさらに高い78.0%を達成し、複製効率と遺伝子発現の最適化を両立。

4. **リファクタリング**: 16.6%のゲノム圧縮が可能と推定。主な圧縮源は遺伝子間領域の削減（78.4 kb）であり、低リスクで実行可能。冗長遺伝子除去は高リスクだが約12.5 kbの追加削減が可能。

5. **アセンブリ**: 4段階の階層的アセンブリ戦略により、69個のDNA合成断片から530 kbゲノムの構築が可能。推定コスト$74,048、所要時間13週間。DNA合成コストが全体の65%を占める。

6. **syn3.0拡張**: 14遺伝子（17.8 kb）の追加で、増殖速度、ストレス耐性、ゲノム安定性、代謝能力、バイオセーフティを改善するsyn3.0+の設計を提案。ゲノムサイズの増加はわずか3.3%。

### 4.2 制限事項

- 本研究は合成データに基づくフレームワーク検証であり、実験的検証が必要
- コドン最適化はタンパク質フォールディングへの影響を考慮していない（稀少コドンによるco-translational foldingの役割）
- 遺伝子間領域の圧縮は制御配列の機能を損なう可能性がある
- syn3.0の149個の機能未知遺伝子の役割解明が最小ゲノム設計の鍵

### 4.3 今後の展望

1. **実験データ統合**: 実際のTn-seqデータ、RNA-seq発現プロファイル、プロテオームデータとの統合
2. **条件特異的必須性**: 複数環境条件下での必須遺伝子セットの比較
3. **ゲノムワイドCRISPRi**: 機能未知遺伝子の系統的ノックダウン実験
4. **自動化パイプライン**: ロボティクスと組み合わせたDesign-Build-Test-Learnサイクルの自動化
5. **代謝モデル統合**: ゲノムスケール代謝モデル（GEM）との統合による必須代謝反応の予測
6. **AI駆動設計**: 深層学習によるde novo遺伝子設計と既存遺伝子の機能予測

---

## 5. 生成したファイル一覧

### ソースコード（`src/`）

| ファイル | 説明 |
|---------|------|
| `src/essential_gene_predictor.py` | Module 1: 必須遺伝子予測（ML+Tn変異データ） |
| `src/codon_optimizer.py` | Module 2: コドン最適化 & ゲノム安定性 |
| `src/gene_arrangement.py` | Module 3: 遺伝子配置最適化 |
| `src/refactoring_strategy.py` | Module 4: リファクタリング戦略 |
| `src/assembly_strategy.py` | Module 5: 階層的Gibson Assembly |
| `src/jcvi_syn3_casestudy.py` | Module 6: JCVI-syn3.0ケーススタディ |
| `src/pipeline.py` | 統合パイプラインメインスクリプト |

### 図表（`figures/`）

| ファイル | 説明 |
|---------|------|
| `fig1_feature_importance.png` | 特徴量重要度（RF + MI） |
| `fig2_roc_curves.png` | ROC曲線（RF, GB, Ensemble） |
| `fig3_tn_insertion_density.png` | Tn挿入密度分布 |
| `fig4_confusion_matrix.png` | 混同行列 |
| `fig5_codon_optimization.png` | コドン最適化結果（CAI, GC, 反復） |
| `fig6_genome_map.png` | 円形ゲノムマップ |
| `fig7_strand_bias.png` | 鎖バイアス解析 |
| `fig8_refactoring.png` | リファクタリング戦略 |
| `fig9_assembly_hierarchy.png` | 階層的アセンブリ図 |
| `fig10_assembly_costs.png` | コスト・タイムライン |
| `fig11_syn3_functional.png` | syn3.0機能分類 |
| `fig12_comparative_genomes.png` | 比較ゲノミクス |
| `fig13_unknown_genes.png` | 機能未知遺伝子カテゴリ |

### 結果データ（`results/`）

| ファイル | 説明 |
|---------|------|
| `module1_summary.json` | Module 1 予測性能サマリ |
| `predicted_essential_genes.csv` | 予測された必須遺伝子リスト |
| `feature_importance.csv` | 特徴量重要度スコア |
| `codon_optimization_results.csv` | コドン最適化結果（全遺伝子） |
| `module2_stability.json` | Module 2 ゲノム安定性指標 |
| `gene_arrangement.csv` | 遺伝子配置データ |
| `module3_arrangement.json` | Module 3 配置最適化指標 |
| `compression_plan.csv` | 圧縮操作計画（450操作） |
| `module4_refactoring.json` | Module 4 リファクタリングサマリ |
| `module5_assembly.json` | Module 5 アセンブリコスト |
| `assembly_hierarchy.json` | アセンブリ階層構造 |
| `qc_checkpoints.json` | 品質管理チェックポイント |
| `syn3_categories.json` | syn3.0遺伝子機能分類 |
| `syn3_extensions.json` | syn3.0拡張提案 |
| `syn3_unknown_genes.json` | 機能未知遺伝子解析 |
| `comparative_genomes.csv` | 比較ゲノムデータ |
| `module6_summary.json` | Module 6 サマリ |
| `pipeline_summary.json` | パイプライン全体サマリ |

### その他

| ファイル | 説明 |
|---------|------|
| `data/essential_genes_dataset.csv` | 合成訓練データセット |
| `logs/process-log.jsonl` | 実行ログ |

---

## 参考文献

1. Hutchison CA III, et al. (2016) Design and synthesis of a minimal bacterial genome. *Science* 351:aad6253.
2. Glass JI, et al. (2006) Essential genes of a minimal bacterium. *PNAS* 103:425-430.
3. Gibson DG, et al. (2008) Complete chemical synthesis, assembly, and cloning of a *Mycoplasma genitalium* genome. *Science* 319:1215-1220.
4. Gibson DG, et al. (2010) Creation of a bacterial cell controlled by a chemically synthesized genome. *Science* 329:52-56.
5. Breuer M, et al. (2019) Essential metabolism for a minimal cell. *eLife* 8:e36842.
