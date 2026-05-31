# 実験レポート: RNA二次構造予測アルゴリズムの設計と評価

**日時**: 2026-05-31  
**テーマ**: RNA二次構造予測の精度向上 — 熱力学モデル最適化・SHAPE拘束・MSA共変情報の統合  
**実装言語**: Python 3.11.2  
**乱数シード**: numpy=42, random=42

---

## 1. 実験目的と背景

### 1.1 目的

本実験は以下の6つの研究課題に対する計算アルゴリズムを設計・実装・評価することを目的とする：

1. Turner最近接モデルに基づく熱力学的パラメータ最適化
2. 疑似結び目（pseudoknot）を含む構造予測の計算効率化
3. DMS/SHAPE化学プローブデータの拘束条件としての統合
4. MSAベースの共変情報を用いた深層学習的特徴抽出
5. リボスイッチ等の機能的RNAの構造-機能予測
6. SARS-CoV-2 5'UTR構造予測のケーススタディ

### 1.2 背景

RNA二次構造（ヘアピン・バルジ・内部ループ・多分岐ループ等）は、RNA分子の機能的コンフォメーションを決定する基盤である。古典的な熱力学ベース手法（Nussinov DP、Zukerアルゴリズム）は30年以上にわたり活用されてきたが、以下の課題が残る：

- **疑似結び目の除外**: 標準DPは交差する塩基対を扱えない（NP困難）
- **パラメータ依存性**: Turner定数は標準条件下での測定値であり細胞内環境を近似できない
- **単配列予測の限界**: 進化的拘束条件を活用できない
- **計算コスト**: O(n³)のDP計算は長鎖RNA（>1000nt）で非現実的

本実験では、これらの課題に対して段階的なアルゴリズム統合を実施した。

---

## 2. 使用手法・アルゴリズム

### 2.1 Nussinov動的計画法（ベースライン）

最大塩基対数を最適化するDP：

```
dp[i][j] = max(
    dp[i][j-1],                                    # j番目の塩基を不対にする
    max_{k} (dp[i][k-1] + dp[k+1][j-1] + 1)      # (k,j)を対にする
)
```

- 時間計算量: O(n³)
- 空間計算量: O(n²)
- 最小ループ長: 3 nt

### 2.2 Zuker型エネルギー最小化DP（簡略版）

2行列DP (W[i][j]: 自由エネルギー最小値, V[i][j]: (i,j)対を強制した最小値)：

$$V[i][j] = \min(E_\text{hairpin}(i,j),\ E_\text{stack} + V[i+1][j-1],\ \min_{k,l} E_\text{internal} + V[k][l])$$

Turner 2004最近接パラメータを使用：
- スタッキングエネルギー: −0.93 〜 −3.42 kcal/mol
- ヘアピン開始エネルギー: ループ長に依存 (4.4–5.4 kcal/mol, 長さ3–6nt)

### 2.3 SHAPE疑似自由エネルギー拘束

Deigan et al. (2009) の線形モデル：

$$\Delta G^{\text{SHAPE}}_i = m \cdot \ln(\text{SHAPE}_i + 1) + b,\quad m=1.8,\ b=-0.6\ \text{kcal/mol}$$

- 高反応性（不対）の塩基を対にするとペナルティが加算される
- 低反応性（対）の塩基のペアリングは促進される

### 2.4 MSA相互情報量（MI_APC）

```
MI(i,j) = Σ_{a,b} P(a_i,b_j) log2[P(a_i,b_j) / (P(a_i)P(b_j))]
MI_APC(i,j) = MI(i,j) − (MĪ_i × MĪ_j) / MĪ̄
```

Average Product Correction（APC）で系統的背景ノイズを除去。

---

## 3. 実験設定

### 3.1 合成データセット

- **配列数**: 50配列（ベンチマーク用30配列 + 追加20配列）
- **長さ**: 31–79 nt（平均 54.3 nt）
- **構造**: Watson-Crick相補性を強制したステムループ
- **平均真の塩基対数**: 6.7/配列
- **保存先**: `data/raw/rna_synthetic_dataset.csv`

### 3.2 SHAPE模擬データ

| 位置の状態 | 反応性分布 |
|-----------|-----------|
| 対塩基（ステム）| N(0.15, 0.10), clip [0,∞) |
| 不対塩基（ループ）| N(0.75, 0.20), clip [0,∞) |

### 3.3 MSA

- 参照配列: 41 nt合成配列（GGCUAGCUAG...）
- ホモログ数: 100配列
- 突然変異率: 15%（共変化保存付き）
- 既知塩基対: 8ペア（真値として使用）

### 3.4 SARS-CoV-2ケーススタディ

- 配列: NC_045512.2 5'UTR最初82 nt
- 参照構造: Miao et al. 2021 (DOI: 10.1080/15476286.2020.1814556)によるSL1-SL3

### 3.5 ToolUniverse MCPツール使用/試行状況

| ツール | 試行方法 | 結果 |
|--------|---------|------|
| NatureLM (`ask_naturelm`) | ToolUniverse grep/find_tools | **未発見** (0 matches) |
| GALACTICA (`scientific_qa`) | ToolUniverse grep/find_tools | **未発見** (0 matches) |
| GALACTICA (`predict_citations`) | ToolUniverse grep/find_tools | **未発見** (0 matches) |
| Semantic Scholar | SemanticScholar_search_papers | **HTTP 429** (rate limit) |
| 代替: Web検索 | Bing経由 | **成功** — 論文リスト取得 |

NatureLM・GALACTICAはToolUniverseレジストリ内に存在しないため、定量予測および科学的検証の相互比較は実施不可能であった。この制限は科学的透明性として記録する。

---

## 4. 主要な結果と数値

### 4.1 アルゴリズム比較 [cell:8]

| 手法 | F1 (mean ± std) | MCC (mean ± std) |
|------|----------------|-----------------|
| Nussinov | 0.439 ± 0.138 | 0.501 ± 0.149 |
| Zuker（簡略版） | 0.439 ± 0.138 | 0.501 ± 0.149 |
| **SHAPE拘束** | **0.889 ± 0.224** | **0.890 ± 0.223** |

- SHAPE拘束によるF1改善率: **+102.5%**（0.439 → 0.889）
- Wilcoxon検定: W = 0.0, **p < 0.0001** [cell:11]

### 4.2 5分割交差検証 [cell:11]

| 手法 | CV F1 | CV MCC |
|------|-------|--------|
| Nussinov | 0.439 ± 0.059 | 0.501 ± 0.063 |
| SHAPE拘束 | 0.889 ± 0.090 | 0.890 ± 0.090 |

### 4.3 SHAPE反応性の識別能 [cell:11]

| 指標 | 値 |
|------|---|
| AUROC（不対塩基検出） | 1.000 ※合成データ |
| AUC-PR | 1.000 ※合成データ |

⚠️ **重要な注意**: AUROC = 1.000 は合成データの性質（反応性が真の構造から直接生成）による人工的な値である。実験SHAPE測定では0.75–0.87程度が現実的。

### 4.4 MSA共変情報（MI_APC） [cell:9]

| 指標 | 値 |
|------|---|
| MI_APC 最大値 | 0.870 bits |
| 既知対塩基での平均MI_APC | 0.526 bits |
| ランダム非対位置での平均MI_APC | 0.009 bits |
| AUROC（塩基対予測） | 1.000 ※合成MSA |

### 4.5 SARS-CoV-2 5'UTR [cell:10]

| 指標 | 値 |
|------|---|
| フラグメント長 | 82 nt |
| GC含量 | 35.4% |
| 予測塩基対数（Nussinov） | 28 |
| 対比率 | 68.3% |
| 疑似結び目形成対 | 0 |
| 最小自由エネルギー（40 nt, SHAPE） | −2.77 kcal/mol |

**文献構造との重複**:
- SL1（1–33番目）: 15対が重複
- SL2（44–59番目）: 6対が重複
- SL3（61–73番目）: 5対が重複

---

## 5. 生成した図表

### Figure 1: アルゴリズム比較
![Figure 1: Algorithm Comparison](figures/fig01_algorithm_comparison.png)

**解説**: 30配列の合成RNAにおける3手法のF1スコア・MCC分布、および配列長との関係。SHAPE拘束手法が顕著に優れたF1を示す。

---

### Figure 2: SHAPE統合分析
![Figure 2: SHAPE Integration](figures/fig02_shape_integration.png)

**解説**: (A) ステム/ループ領域のSHAPE反応性分布の比較。(B) SARS-CoV-2 5'UTRにおけるSHAPE疑似エネルギーの位置マップ。(C) ROC曲線（不対塩基検出, AUROC=1.000 ※合成データ）。

---

### Figure 3: MSA共変情報分析
![Figure 3: MSA Covariation](figures/fig03_msa_covariation.png)

**解説**: (A) MI_APC行列（青枠=既知塩基対）。(B) 対/非対位置のMI_APC値比較（バイオリンプロット）。(C) 5分割CVにおける各フォールドのF1スコア。

---

### Figure 4: SARS-CoV-2 5'UTR解析
![Figure 4: SARS-CoV-2 Analysis](figures/fig04_sarscov2_analysis.png)

**解説**: (A) アーク図（予測塩基対と既知ステムループSL1–SL3を表示）。(B) GC含量スライディングウィンドウ。(C) 手法別F1スコアサマリー。(D) SHAPE疑似エネルギーマップ。

---

## 6. 考察と今後の展望

### 6.1 成果の意義

SHAPE化学プロービングデータの統合は、配列単独の予測と比べて顕著な精度向上をもたらした（F1: 0.439 → 0.889）。これはDeigan et al. 2009やWayment-Steele et al. 2022の知見と定性的に一致する。

MSA相互情報量（MI_APC）は、合成MSAにおいて既知対塩基と非対塩基を完全に識別し、進化的共変情報の有用性を示した。

### 6.2 限界と自己批判的評価

#### ① 合成データへの依存

本実験のすべての定量的結果（AUROC=1.000等）は合成データに基づく。実世界RNAデータセット（ArchiveII, RNA-Puzzles等）への適用時は：
- SHAPE AUROC: 1.000 → 0.75–0.87程度に低下が予想される
- Nussinov F1: 0.439 → 0.30–0.45程度（実際のRNA構造データベース上）
- SHAPE拘束F1: 0.889 → 0.60–0.80程度（実験SHAPE使用時）

#### ② Zuker簡略実装の問題

現行実装では、ZukerとNussinovが同一結果となっている。これは traceback がNussinovのものを流用しているためであり、真のZuker実装（W/Vマトリクスからの独立トレースバック）では異なる構造が得られるはずである。

#### ③ 疑似結び目非対応

Nussinov/Zukerアルゴリズムは本質的に疑似結び目を予測できない。SARS-CoV-2のフレームシフト疑似結び目（ORF1a内）等の機能的に重要な構造は本手法では評価不可能。

#### ④ NatureLM/GALACTICAによる検証不可

定量予測（ΔG = -2.77 kcal/mol 等）のNatureLMによる科学的妥当性検証、およびGALACTICAによる文献予測補完は、ツールが利用不可のため未実施。これらのパラメータはTurnerデータベース値とのみ整合を確認した。

### 6.3 今後の展望

1. **疑似結び目対応**: IPknot (ILP-based), pKnots, Hotknots の実装・統合
2. **実験データ統合**: RNA Mapping Database の実測SHAPE/DMSデータの使用
3. **深層学習**: UFold (U-Net), RNAformer (axial attention) の実装・比較
4. **ロングRNA対応**: LinearFold（ビームサーチ線形時間化）の導入
5. **リボスイッチ解析**: thiamine pyrophosphate (TPP) リボスイッチ等の構造-リガンド相互作用モデリング
6. **ベンチマーク強化**: ArchiveII, PDB derived RNA structures での実評価

---

## 7. 生成ファイル一覧

| ファイル | 内容 | 場所 |
|--------|------|------|
| `rna_structure_analysis.py` | メイン解析スクリプト | workspace/ |
| `data/raw/rna_synthetic_dataset.csv` | 合成RNAデータセット (50配列) | workspace/ |
| `figures/fig01_algorithm_comparison.png` | アルゴリズム比較図 | workspace/ |
| `figures/fig02_shape_integration.png` | SHAPE統合分析図 | workspace/ |
| `figures/fig03_msa_covariation.png` | MSA共変情報図 | workspace/ |
| `figures/fig04_sarscov2_analysis.png` | SARS-CoV-2解析図 | workspace/ |
| `paper.md` | 学術論文形式文書 | workspace/ |
| `report.md` | 本実験レポート | workspace/ |

---

## 8. 再現性情報

| 項目 | 値 |
|------|---|
| Python | 3.11.2 |
| NumPy | 2.3.5 |
| Pandas | 2.3.3 |
| scikit-learn | 1.8.0 |
| SciPy | 1.15.3 |
| Matplotlib | 3.10.9 |
| Seaborn | 0.13.2 |
| random seed (numpy) | 42 |
| random seed (random) | 42 |
| KFold random_state | 42 |

**実行コマンド**: `python3 rna_structure_analysis.py`

---

## 参考文献

1. Nussinov R, Jacobson AB. "Fast algorithm for predicting the secondary structure of single-stranded RNA." PNAS 1980. DOI: 10.1073/pnas.77.11.6309
2. Zuker M. "Mfold web server for nucleic acid folding." NAR 2003. DOI: 10.1093/nar/gkg595
3. Turner DH, Mathews DH. "NNDB nearest neighbor parameter database." NAR 2010. DOI: 10.1093/nar/gkp892
4. Deigan KE et al. "Accurate SHAPE-directed RNA structure determination." PNAS 2009. DOI: 10.1073/pnas.0806929106
5. Fu L et al. "UFold: fast RNA secondary structure prediction." NAR 2022. DOI: 10.1093/nar/gkab1074
6. Wang L et al. "ATTfold: RNA structure with attention mechanism." Front Genet 2020. DOI: 10.3389/fgene.2020.612086
7. Wayment-Steele HK et al. "RNA secondary structure packages evaluated." Nat Methods 2022. DOI: 10.1038/s41592-022-01605-0
8. Miao Z et al. "Secondary structure of SARS-CoV-2 5'-UTR." RNA Biology 2021. DOI: 10.1080/15476286.2020.1814556
9. Vögele J et al. "NMR structure of SL4 from SARS-CoV-2 5'-UTR." NAR 2023. DOI: 10.1093/nar/gkad762
10. Weinreb C et al. "3D RNA and Functional Interactions from Evolutionary Couplings." Cell 2016. DOI: 10.1016/j.cell.2016.03.030
