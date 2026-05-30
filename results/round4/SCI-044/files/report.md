# RNA二次構造予測実験レポート

## 実験目的と背景

RNA二次構造の正確な予測は、遺伝子発現調節・ウイルス複製・医薬品設計において基盤的重要性を持つ。Turner最近接隣接エネルギーモデルによる動的計画法（DP）は半世紀にわたって標準的手法であるが、（1）疑似結び目（pseudoknot）を含む構造に対する計算量の爆発、（2）化学プローブ（SHAPE/DMS）データの統合、（3）多種アライメント（MSA）由来の共変情報の活用、という三つの課題が依然として残る。本研究では、これら全要素を統合した **HybridFold** アルゴリズムを設計・実装し、合成ベンチマークと SARS-CoV-2 5'UTR ケーススタディで評価した。

---

## 先行研究調査（ToolUniverse MCP使用）

### 使用ツールと試行状況

| ツール名 | 試行結果 |
|---------|---------|
| `SemanticScholar_search_papers` | 最初の3リクエストで API 400/429 エラー（レート制限）。その後のリクエストも断続的にエラー |
| `PubMed_search_articles` | 正常動作。5クエリ成功 |
| `Crossref_search_works` | 正常動作。2クエリ成功 |

Semantic Scholar API はレート制限（1 req/sec without API key）により断続的に失敗した。PubMed・Crossref を主要検索ソースとして使用した。これは科学的透明性の観点から記録する。

### 発見した主要論文（2020年以降、5件以上）

| # | タイトル | 著者 | 年 | DOI | 主要知見 |
|---|---------|------|----|----|---------|
| 1 | UFold: Fast and Accurate RNA Secondary Structure Prediction with Deep Learning | Fu et al. | 2020 | 10.1101/2020.08.17.254896 | U-Netアーキテクチャでドットブラケット行列を直接予測。従来MFE法より高速 |
| 2 | CParty: hierarchically constrained partition function of RNA pseudoknots | Gray et al. | 2024 | 10.1093/bioinformatics/btae748 | HFoldに基づく疑似結び目分配関数。O(n³)時間・O(n²)空間。SARS-CoV-2治療標的解析 |
| 3 | SparseRNAFolD: optimized sparse RNA pseudoknot-free folding with dangle consideration | Gray et al. | 2024 | 10.1186/s13015-024-00256-4 | スパース化MFE予測。ダングル貢献を考慮しつつ時間・空間効率を改善 |
| 4 | DinoKnot: Duplex Interaction of Nucleic Acids With PseudoKnots | Newman et al. | 2024 | 10.1109/TCBB.2024.3362308 | 疑似結び目を含む核酸二重鎖相互作用予測。SARS-CoV-2ゲノムで検証 |
| 5 | AliNA: deep learning program for RNA secondary structure prediction | Nasaev et al. | 2023 | 10.1002/minf.202300113 | アライメントベースのデータ拡張で非相同配列に対する汎化性を改善 |
| 6 | Machine learning modeling of RNA structures | Wu et al. | 2023 | 10.1093/bib/bbad210 | RNA構造ML手法の包括的レビュー。熱力学原理とDLの統合の重要性を指摘 |
| 7 | DivideFold+: AI-based tool for RNA secondary structure prediction | Omnes et al. | 2026 | 10.1016/j.jmb.2026.169865 | 長鎖RNAを独立サブドメインに分割してDL予測。データ拡張戦略を提案 |
| 8 | Diverse database and machine learning model (eFold) | de Lajarte et al. | 2026 | 10.1126/sciadv.adz4967 | EvoformerインスパイアのRNA二次構造予測。化学プローブデータで訓練 |
| 9 | Unveiling hidden structural patterns in the SARS-CoV-2 genome | Ziesel & Jabbari | 2024 | 10.1371/journal.pone.0298164 | SARS-CoV-2の40の構造的ゲノム領域を同定。複数予測手法の比較 |
| 10 | The trRosettaRNA server for RNA structure prediction | Wang et al. | 2026 | 10.1038/s41596-026-01356-8 | end-to-end DLによるRNA 3D構造予測サーバー |

### 先行研究の課題・限界

1. **疑似結び目の計算複雑性**: 一般的な疑似結び目を含む構造予測はNP完全。CPartyでもO(n³)に抑えるためhierarchical仮定が必要
2. **実験データの不足**: eFoldが指摘する通り、実験的に検証されたRNA構造データが不足
3. **化学プローブ統合の最適化**: SHAPE/DMS疑似エネルギー変換パラメータの最適化は配列依存性が高い
4. **MSA共変情報の定量的活用**: DCA/MI情報をMFE計算に統合する効率的な枠組みが未確立
5. **配列長スケーラビリティ**: 1000 nt以上の長鎖RNAに対してO(n³)でも実用的でない場合がある

---

## 使用した手法・アルゴリズムの概要

### HybridFold アーキテクチャ

```
入力: RNA配列 S (n nt) + SHAPE/DMS反応性 r + MSA M
  │
  ├── [1] SHAPE正規化 (2-8%法)
  │         r_norm → 疑似エネルギー ΔG_SHAPE = 1.8·ln(r+1) − 0.6
  │
  ├── [2] MSA共変情報 (相互情報量)
  │         MI(i,j) = H(i) + H(j) − H(i,j)
  │
  ├── [3] 統合重みベクトル
  │         w[i] = ΔG_SHAPE[i] − 0.5·⟨MI(i,·)⟩
  │
  ├── [4] Turner MFE-DP (O(n³))
  │         V[i][j] = min{ヘアピン, スタッキング, 内部ループ} + w[i] + w[j]
  │         W[i][j] = min{W[i][j-1], min_k{W[i][k-1] + V[k][j]}}
  │
  ├── [5] 階層的疑似結び目検出
  │         コア構造G(pseudoknot-free)を取得 → 未対合領域を再折り畳み
  │
  └── 出力: ドットブラケット記法, MFE (kcal/mol), 疑似結び目座標
```

### ベースライン比較対象

| 手法 | アルゴリズム | 計算量 | 特徴 |
|------|------------|-------|------|
| **Nussinov** | 最大塩基対DP | O(n³) | エネルギーモデルなし |
| **Turner-MFE** | 最小自由エネルギーDP | O(n³) | Turner最近接エネルギー |
| **Turner+SHAPE** | MFE + 化学プローブ拘束 | O(n³) | 本研究の中間手法 |
| **HybridFold** | MFE + SHAPE + MI共変 | O(n³) | 本研究の提案手法 |

---

## 主要な結果と数値

### 1. 5-fold 交差検証ベンチマーク（合成RNA、n=40/60/80、各fold 20サンプル）

| 手法 | F1 (mean±SD) | Precision (mean±SD) | Recall (mean±SD) | 平均時間 (ms) |
|------|-------------|---------------------|-----------------|-------------|
| Nussinov | 0.344 ± 0.317 | 0.285 ± 0.264 | 0.437 ± 0.399 | 6.77 ± 4.82 |
| Turner-MFE | 0.354 ± 0.332 | 0.308 ± 0.293 | 0.416 ± 0.387 | 43.87 ± 37.53 |
| Turner+SHAPE | **0.385 ± 0.349** | **0.343 ± 0.317** | 0.440 ± 0.391 | 44.05 ± 37.63 |
| HybridFold | **0.385 ± 0.349** | 0.342 ± 0.316 | **0.441 ± 0.392** | 56.01 ± 43.60 |

![Figure 1: ベンチマーク結果](figures/benchmark_results.png)

**解釈**: SHAPE統合により F1 が +8.8% 向上（0.354→0.385）。MSA共変情報の追加は限定的（標準偏差の範囲内）。

### 2. 配列長スケーリング

![Figure 2: 計算量スケーリング](figures/length_scaling.png)

すべての手法が O(n³) の理論計算量に従うことが確認された。線形-線形スケールでは n=300 nt で Turner-MFE が約 2000 ms を要する。

### 3. SARS-CoV-2 5'UTR 構造予測（最初の100 nt）

| 手法 | 予測塩基対数 | MFE (kcal/mol) |
|------|-----------|---------------|
| Nussinov | 35 | — |
| Turner-MFE | 19 | −18.3 |
| Turner+SHAPE+DMS | 20 | −16.9 |
| Turner+Pseudoknot | 28 | −27.9 |

ドットブラケット記法（Turner+SHAPE+DMS）:
```
(((.((((.(...)))))...............))........
         ↑SL1                    ↑SL2領域
```

![Figure 3: SARS-CoV-2 5'UTR アーク図](figures/sars_structure.png)

### 4. 疑似結び目検出ベンチマーク（H型疑似結び目合成配列）

| 配列長 | HybridFold F1 | Turner-MFE F1 |
|-------|--------------|--------------|
| n=50 | 0.042 ± 0.042 | 0.048 ± 0.048 |
| n=70 | 0.016 ± 0.026 | 0.018 ± 0.029 |
| n=90 | 0.019 ± 0.022 | 0.022 ± 0.025 |

![Figure 4: 疑似結び目ベンチマーク](figures/pseudoknot_benchmark.png)

**注**: 両手法とも正確な疑似結び目ペア回収率は低い（F1 <0.05）。これは疑似結び目予測の本質的難しさを反映する。

### 5. SHAPE データ品質感度分析

| ノイズレベル (σ) | F1 (mean±SD) |
|--------------|-------------|
| 0.0 (完全) | 0.899 ± 0.035 |
| 0.1 | 0.869 ± 0.031 |
| 0.2 | 0.852 ± 0.036 |
| 0.3 | 0.848 ± 0.037 |
| 0.5 | 0.816 ± 0.049 |
| 0.8 (高ノイズ) | 0.800 ± 0.045 |

![Figure 5: SHAPE感度分析](figures/shape_sensitivity.png)

SHAPE反応性ノイズが増加しても F1 の低下は緩やか（σ=0 vs σ=0.8: Δ=0.099）。

### 6. MSA共変情報ヒートマップ（SARS-CoV-2 5'UTR）

![Figure 6: 相互情報量ヒートマップ](figures/covariation_heatmap.png)

---

## 考察と今後の展望

### 主要知見
1. **SHAPE統合の有効性**: Turner+SHAPE は Turner-MFE に対して F1 で平均 +8.8% の改善。ノイズ耐性も高い（σ=0.8 でも F1=0.80 維持）
2. **疑似結び目予測の困難さ**: 階層的アプローチでも正確なペア回収率は低い（<5%）。これは CParty（2024）や DinoKnot（2024）の知見と整合する
3. **MSA共変情報の貢献**: 今回の実装では SHAPE に対して限定的な追加効果。より多くのMSA配列（50以上）と DCA（Direct Coupling Analysis）の使用が推奨される
4. **スケーラビリティ**: すべての手法が O(n³) に従い、n>300 では実用的な時間内に収まらない

### 限界
- 合成データを使用（実験的に検証された構造データベースによる評価が必要）
- 疑似エネルギーパラメータ（m=1.8, b=-0.6）は Mathews 2009 の値を固定使用
- 簡略化した Turner エネルギー行列（全パラメータのサブセット）

### 今後の課題
1. LinearFold/LinearPartition 的な線形近似による O(n) アルゴリズムへの拡張
2. RNA-FM・EternaFold 等の大規模言語モデルとの統合
3. RNApdbee/Rfam データセットによる本格的評価
4. SARS-CoV-2全ゲノムへのスライディングウィンドウ適用

---

## 生成したファイル一覧

| ファイル | 内容 |
|---------|-----|
| `src/rna_structure.py` | HybridFoldコアアルゴリズム（Turner DP, 疑似結び目, SHAPE統合, MI計算） |
| `src/experiment.py` | ベンチマーク実験スクリプト |
| `figures/benchmark_results.png` | 5-fold CVベンチマーク結果（Figure 1） |
| `figures/length_scaling.png` | 計算量スケーリング（Figure 2） |
| `figures/sars_structure.png` | SARS-CoV-2 5'UTRアーク構造図（Figure 3） |
| `figures/pseudoknot_benchmark.png` | 疑似結び目ベンチマーク（Figure 4） |
| `figures/shape_sensitivity.png` | SHAPE感度分析（Figure 5） |
| `figures/covariation_heatmap.png` | MI共変情報ヒートマップ（Figure 6） |
| `report.md` | 本レポート |
| `paper.md` | 学術論文形式文書 |
