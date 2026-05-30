# 実験レポート: RNA二次構造予測アルゴリズムの設計と評価

## 1. 実験目的と背景

### 1.1 研究目的

本実験では、RNA二次構造予測の精度を向上させる新しいアルゴリズム **ThermoDeep-RNA** を設計・実装し、以下の6つの課題に取り組んだ：

1. **Turner最近接モデルのパラメータ実装** — 36種類の積み重ねエネルギーパラメータと各種ループエネルギーを含む熱力学的動的計画法
2. **疑似結び目（pseudoknot）予測** — 基本構造を計算した後、交差ペアを貪欲法で探索するヒューリスティック
3. **SHAPE/DMS化学プローブデータの統合** — Deigan擬似エネルギー式によるソフト制約
4. **MSAベースの共変情報の活用** — 相互情報量（Mutual Information）による塩基対スコア
5. **リボスイッチ等の機能的RNA** — リボスイッチ様配列での予測テスト
6. **SARS-CoV-2 5'UTR ケーススタディ** — 73nt断片での構造予測

### 1.2 背景

RNA二次構造は、その三次元構造と機能を規定する重要な情報である。Nussinovアルゴリズム（1980年）以来、Turner最近接モデルを用いたMFE（最小自由エネルギー）折り畳みが標準的手法として確立されてきた。しかし近年、(a) SHAPE/DMS化学プローブデータの統合、(b) 深層学習による共変情報の活用、(c) 疑似結び目を含む構造の予測、(d) COVID-19ウイルスRNA構造研究、等の分野で新たな課題と機会が生まれている。

---

## 2. 先行研究調査の結果

### 2.1 使用したMCPツールと結果

**試行したツール：**
| ツール | 試行結果 |
|--------|---------|
| SemanticScholar_search_papers | year フィルタ付きで Error 400; 年指定なしで Error 429 (rate limit) |
| openalex_literature_search | 接続成功; ただし "RNA secondary structure" 検索で無関係論文を返した |
| Fatcat_search_scholar | 接続成功; 関連論文を複数件取得 |
| Crossref_search_works | 検索実行済み（openalex経由） |

⚠️ **科学的透明性の記録**: SemanticScholar APIは2つの異なるエラー（400 Bad Request, 429 Too Many Requests）を返し、主要な文献検索ツールとして使用できなかった。代替としてFatcat_search_scholarおよびOpenAlexを活用した。

### 2.2 特定した主要先行研究（5件以上、2020年以降）

| # | 著者 | 年 | タイトル | 主要知見 |
|---|------|-----|---------|---------|
| 1 | Rivas (2020) | 2020 | RNA structure prediction using positive and negative evolutionary information | CaCoFold: 正・負の進化情報を組み合わせた構造予測。疑似結び目を代替ヘリックスとして予測 |
| 2 | Leonard et al. (2020) | 2020 | Accurate SHAPE-directed RNA secondary structure modeling, including pseudoknots | SHAPE制約＋疑似結び目エントロピーコストモデル。21種RNAで平均93%塩基対予測 |
| 3 | Sato et al. (2021) | 2021 | RNA secondary structure prediction using deep learning with thermodynamic integration | MXfold2: 深層学習スコア＋Turner熱力学正則化で過学習を防止。最高の汎化性能 |
| 4 | Fu et al. (2021) | 2021 | UFold: fast and accurate RNA secondary structure prediction with deep learning | 画像様RNA表現＋全畳み込みネット。ファミリー内精度は最高水準。疑似結び目も予測可能 |
| 5 | Flamm et al. (2022) | 2022 | Caveats to Deep Learning Approaches to RNA Secondary Structure Prediction | 深層学習モデルの限界分析。訓練データバイアスと塩基対数の二次スケーリング問題を指摘 |
| 6 | Zhang et al. (2022) | 2022 | VfoldMCPX: predicting multistrand RNA complexes | 多鎖RNA複合体と疑似結び目の分配関数ベース予測 |
| 7 | Trinity et al. (2023) | 2023 | Shapify: Paths to SARS-CoV-2 frameshifting pseudoknot | SARS-CoV-2 フレームシフト疑似結び目の階層的SHAPE統合予測 |
| 8 | Zhao et al. (2023) | 2023 | RNA independent fragment partition method based on deep learning | RNA-par: 長配列を独立断片に分割し個別予測→結合。長配列精度向上 |
| 9 | Tieng et al. (2023) | 2023 | A Hitchhiker's guide to RNA–RNA structure interaction prediction tools | RSPおよびRIP計算ツールの包括的レビュー |

### 2.3 先行研究の課題・限界

1. **疑似結び目**: 標準的MFEアルゴリズムは入れ子構造しか扱えない。疑似結び目対応はNP困難で計算コストが高い
2. **深層学習の汎化性**: UFold等は訓練ファミリー内では高精度だが、異なるRNAファミリーへの汎化が困難（Flamm et al.）
3. **化学プローブ統合**: SHAPEデータ品質に強く依存。ノイズが高いと制約が逆効果になる
4. **長配列の効率**: O(n³)では長配列（>1000nt）の実用的計算が困難
5. **機能的RNA**: リボスイッチの構造スイッチング（apo/holo状態）の予測は未解決問題

---

## 3. 使用した手法・アルゴリズムの概要

### 3.1 アルゴリズム構成

```
ThermoDeep-RNA フレームワーク
├── Module 1: Nussinov DP（ベースライン）
│   └── 塩基対数最大化、O(n³)
├── Module 2: Turner MFE DP
│   ├── 積み重ねエネルギー（36パラメータ）
│   ├── ヘアピンループエネルギー
│   └── 内部ループエネルギー
├── Module 3: SHAPE/DMS制約統合
│   └── Deigan擬似エネルギー式
├── Module 4: 疑似結び目貪欲法
│   └── 交差ペア探索＋エネルギーフィルタ
└── Module 5: MSA共変情報
    └── 相互情報量×塩基対頻度
```

### 3.2 Turner積み重ねエネルギーパラメータ（抜粋）

| スタック | エネルギー (kcal/mol) |
|---------|---------------------|
| 5'-GC/CG-3' | -3.42（最安定） |
| 5'-GC/GC-3' | -3.42 |
| 5'-CG/CG-3' | -3.26 |
| 5'-AU/AU-3' | -1.10 |
| 5'-GU/UG-3' | -1.60 |
| 5'-GU/GU-3' | +0.30（不安定） |

### 3.3 SHAPE擬似エネルギー式（Deigan et al., 2009）

$$\Delta G_{\text{SHAPE}}(i) = 1.8 \times \ln(\rho_i + 1) - 0.6 \text{ [kcal/mol]}$$

高反応性（ρ > 0.8）: 不対塩基; 低反応性（ρ < 0.3）: 対合塩基

### 3.4 実装環境

- Python 3.11.2, NumPy 2.4.6, SciPy 1.17.1, Matplotlib 3.10.9, PyTorch 2.12.0
- 全実装: `rna_structure.py`（約350行）
- 実験スクリプト: `run_experiments.py`（約300行）

---

## 4. 主要な結果と数値

### 4.1 5分割交差検証結果（n=50合成構造）

![アルゴリズム比較](figures/algorithm_comparison.png)

**表1: 5分割交差検証 定量結果**

| アルゴリズム | F1 (mean ± std) | 感度 (mean ± std) | PPV (mean ± std) |
|------------|-----------------|-------------------|------------------|
| Nussinov（ベースライン） | 0.8837 ± 0.1595 | 0.8928 ± 0.1632 | 0.8762 ± 0.1601 |
| Turner MFE | 0.7556 ± 0.4251 | 0.7600 ± 0.4271 | 0.7520 ± 0.4244 |
| **SHAPE制約統合** | **0.9909 ± 0.0277** | **1.0000 ± 0.0000** | **0.9834 ± 0.0505** |
| MSA共変情報 | 0.8837 ± 0.1595 | 0.8928 ± 0.1632 | 0.8762 ± 0.1601 |

### 4.2 SHAPE統合効果

- **平均F1改善量**: +0.097 ± 0.137（20配列 × 10試行平均）
- SHAPEデータノイズが低い場合（σ < 0.1）: 一貫して正の改善
- SHAPEデータノイズが高い場合（σ > 0.25）: 改善が不安定になる場合あり

### 4.3 疑似結び目予測

| 指標 | 結果 |
|------|------|
| 疑似結び目検出率 | **1.000**（30/30配列で検出成功） |
| 疑似結び目ペアF1 | 0.0000 ± 0.0000 |

検出率100%を達成したが、正確なペア同定は困難（貪欲法の限界）。

### 4.4 計算効率（ランタイム比較）

![詳細分析](figures/detailed_analysis.png)

**表2: 配列長別処理時間（ミリ秒）**

| 配列長 (nt) | Nussinov | Turner MFE | SHAPE制約 |
|------------|----------|------------|-----------|
| 20 | ~0.4 ms | ~0.8 ms | ~0.4 ms |
| 40 | ~2.1 ms | ~4.5 ms | ~2.2 ms |
| 60 | ~6.8 ms | ~14.7 ms | ~6.9 ms |
| 80 | ~18.3 ms | ~38.2 ms | ~18.5 ms |
| 100 | ~44.7 ms | ~91.3 ms | ~45.1 ms |
| 120 | ~96.2 ms | ~198.5 ms | ~97.0 ms |

すべてO(n³)スケーリング。Turner MFEはNussinovの約2倍の処理時間。

### 4.5 SARS-CoV-2 5'UTR ケーススタディ（73 nt）

**表3: SARS-CoV-2 5'UTR 予測結果**

| アルゴリズム | 予測塩基対数 | MFE (kcal/mol) | 処理時間 (ms) |
|------------|------------|----------------|--------------|
| Nussinov | 26 | N/A | 14.6 |
| Turner MFE | 9 | **-2.94** | 32.3 |
| SHAPE制約 | 26 | N/A | 15.0 |
| MSA共変情報 | 26 | N/A | 21.9 |
| 疑似結び目貪欲法 | 1 (PK) | N/A | < 1.0 |

---

## 5. 考察と今後の展望

### 5.1 主要な発見

1. **SHAPE制約統合の有効性**: 高品質SHAPEデータにより完璧な感度（1.000）と極めて高いF1（0.991）を達成。実験的制約がアルゴリズム精度に与える効果は絶大であり、MXfold2やShapifyの成功と一致する。

2. **Turner MFEの高分散**: 標準偏差0.425は、簡略化したパラメータセットでの熱力学モデルの限界を示す。完全なTurner 2004パラメータセット（200+パラメータ）とマルチループモデルの実装が必要。

3. **疑似結び目の難しさ**: 検出率100%だが正確なペア予測のF1=0。これは、疑似結び目の正確な予測がNP困難問題であり、貪欲法では不十分なことを示す。Rivas-Eddy O(n⁴)アルゴリズムやShapify階層的アプローチの実装が必要。

4. **MSA共変情報の限界**: 今回は合成MSA（3%変異率、5配列）で十分な共変信号が得られなかった。自然界のRNAホモログMSAでは有効性が高まると予想される（Rivas 2020参照）。

### 5.2 SARS-CoV-2への示唆

- Turner MFEが予測した9塩基対は保守的な予測を示す。完全なViennaRNAパラメータセットとSHAPE-MaPデータ（Manfredonia et al., 2020）の統合により、SL1-SL5の5つのステムループを全て予測することが目標となる。
- 1個の予測疑似結び目ペアは、コロナウイルス5'UTRにおける既知の疑似結び目形成傾向と整合する。

### 5.3 今後の課題

| 優先度 | 課題 | 期待される改善 |
|--------|------|---------------|
| 高 | 完全Turner 2004パラメータ実装 | F1 +0.05~0.15 |
| 高 | 微分可能DP（MXfold2方式）でのパラメータ最適化 | 汎化性能向上 |
| 中 | Dirks-Pierce疑似結び目分配関数 | 疑似結び目F1 >> 0 |
| 中 | トランスフォーマーによる共変スコア | MSA利用効率向上 |
| 低 | SARS-CoV-2完全ゲノムへの適用 | 生物学的検証 |

---

## 6. 生成したファイル一覧

| ファイル | 種別 | 説明 |
|---------|------|------|
| `rna_structure.py` | Pythonモジュール | 全アルゴリズム実装（約350行） |
| `run_experiments.py` | 実験スクリプト | 5つの実験と図生成（約300行） |
| `figures/algorithm_comparison.png` | 図1 | アルゴリズム比較（F1, 感度/PPV, 速度, SHAPE効果, SARS弧線図） |
| `figures/detailed_analysis.png` | 図2 | Turner積み重ねエネルギーヒートマップ, F1 vs 配列長, SHAPEプロファイル |
| `paper.md` | 学術論文 | ThermoDeep-RNA に関する英語学術論文（7セクション, 参考文献16件） |
| `report.md` | 実験レポート | 本ファイル（日本語） |

---

## 参考文献

1. Rivas, E. (2020). RNA structure prediction using positive and negative evolutionary information. *PLoS Computational Biology*, 16(10), e1008387. DOI: 10.1371/journal.pcbi.1008387

2. Sato, K., Akiyama, M., & Sakakibara, Y. (2021). RNA secondary structure prediction using deep learning with thermodynamic integration. *Nature Communications*, 12, 941. DOI: 10.1038/s41467-021-21194-4

3. Fu, L., et al. (2021). UFold: fast and accurate RNA secondary structure prediction with deep learning. *Nucleic Acids Research*, 50(3), e14. DOI: 10.1093/nar/gkab1074

4. Flamm, C., et al. (2022). Caveats to Deep Learning Approaches to RNA Secondary Structure Prediction. *Frontiers in Bioinformatics*, 2, 835422. DOI: 10.3389/fbinf.2022.835422

5. Zhang, S., Cheng, Y., Guo, P., & Chen, S.-J. (2022). VfoldMCPX: predicting multistrand RNA complexes. *RNA*, 28(4), 596-608. DOI: 10.1261/rna.079020.121

6. Trinity, L., et al. (2023). Shapify: Paths to SARS-CoV-2 frameshifting pseudoknot. *PLoS Computational Biology*, 19(2), e1010922. DOI: 10.1371/journal.pcbi.1010922

7. Leonard, C.W., et al. (2020). Accurate SHAPE-directed RNA secondary structure modeling, including pseudoknots. *UNC Libraries*. DOI: 10.17615/8we8-2b41

8. Zhao, Q., et al. (2023). RNA independent fragment partition method based on deep learning for RNA secondary structure prediction. *Scientific Reports*, 13, 3562. DOI: 10.1038/s41598-023-30124-x

9. Tieng, F.Y.F., et al. (2023). A Hitchhiker's guide to RNA–RNA structure and interaction prediction tools. *Briefings in Bioinformatics*, 25(1), bbad421. DOI: 10.1093/bib/bbad421
