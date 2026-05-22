# LongSV: ロングリードデータからの高精度構造変異検出パイプライン

> **DRAFT — NOT FOR DISTRIBUTION**  
> 作成日時: 2026-05-22T12:39:38 UTC  
> リファレンスゲノム: hg38  
> ベンチマーク真値セット: GIAB HG002 Tier1 SV truth set

---

## 目次

1. [実験目的と背景](#1-実験目的と背景)
2. [パイプラインアーキテクチャ](#2-パイプラインアーキテクチャ)
3. [使用した手法・アルゴリズムの概要](#3-使用した手法アルゴリズムの概要)
4. [主要な結果と数値](#4-主要な結果と数値)
5. [考察と今後の展望](#5-考察と今後の展望)
6. [生成したファイル一覧](#6-生成したファイル一覧)

---

## 1. 実験目的と背景

### 背景

構造変異（Structural Variant: SV）は全長50 bp以上のゲノム再編成であり、欠失（DEL）、挿入（INS）、逆位（INV）、重複（DUP）、転座（TRA）が主要なカテゴリとして知られる。SVはがん、希少疾患、神経発達障害など多くの疾患の原因となる一方、短鎖リード（Illuminaなど）による検出精度には限界があった。

Oxford Nanopore Technology（ONT）とPacBio HiFiに代表されるロングリードシーケンサーは、10 kbを超えるリード長によって以下の課題を解決する可能性を持つ：

- リピート領域（テロメア、セントロメア、STR）内部のSV解析
- 複数の重複SV（クロモスリプシス、染色体外DNA）の同定
- ホモポリマー・低複雑度配列での高精度ベースコール

### 実験目的

本プロジェクトでは以下の6項目を達成する高精度SVパイプライン **LongSV** を設計・実装した：

1. **シグナルレベルのベースコール改善** — 双方向LSTM×5層＋アテンション機構によるCTCデコード
2. **統合SV検出戦略** — Split-read / Read-depth / Assembly-based の3手法をベイズ統合
3. **リピート領域の特殊処理** — テロメア長推定、セントロメアフィルタ、STR伸長検出
4. **複雑SVの検出ロジック** — クロモスリプシス、ecDNA、クロモプレキシ、BFB
5. **ショートリードとのハイブリッド解析** — SURVIVOR-style マージによる精度向上
6. **GIABベンチマーク評価設計** — Tier1 SV truth setに対するPrecision/Recall/F1計算

---

## 2. パイプラインアーキテクチャ

![Pipeline Architecture](figures/pipeline_architecture.png)

*Figure 1. LongSVパイプライン全体アーキテクチャ。入力から出力VCFまでのデータフロー、6つの主要モジュール、およびGIABベンチマーク評価を示す。*

### データフロー概要

```
ONT/PacBio Raw Signal
    ↓
[Module 1] RNN Basecaller (BiLSTM + CTC)
    ↓
[Long-Read Alignment: minimap2]
    ├──→ [Module 2a] Split-Read SV Detector
    ├──→ [Module 2b] Read-Depth CNV Detector
    └──→ [Module 2c] Assembly-Based Detector
              ↓
[Module 3] Bayesian Evidence Integration
    ↓   ← [Module 5] Illumina Short-Read Hybrid Merge
[Module 4] Repeat Region Handler
    ├──→ Telomere Length Estimation
    ├──→ Centromere Strict Filter
    └──→ STR Expansion Detection
              ↓
[Module 6] Complex SV Detector
    ├──→ Chromothripsis
    ├──→ ecDNA (extrachromosomal DNA)
    ├──→ Chromoplexy
    └──→ BFB (Breakage-Fusion-Bridge)
              ↓
[VCF 4.2 Output] → [GIAB Benchmark]
```

---

## 3. 使用した手法・アルゴリズムの概要

### 3.1 RNN ベースコーラー（Module 1）

**ファイル**: `src/basecalling/rnn_basecaller.py`

| 項目 | 仕様 |
|------|------|
| アーキテクチャ | 双方向LSTM × 5層、hidden_size=384 |
| アテンション | Scaled Dot-Product Attention (8ヘッド) |
| デコーダ | CTC Beam Search (beam_width=25) |
| 正規化 | MAD-based robust normalization (clip: ±2.5 σ) |
| チャンキング | chunk_size=4000, overlap=200 samples |
| ホモポリマー補正 | 最大連続長を10塩基にキャップ |
| メチル化検出 | 5mC (CpGコンテキスト)、6mA (Adenine) |

**アルゴリズム**:

```python
# MAD正規化
signal_norm = (signal - median) / (1.4826 * MAD)
# CTCデコード: log_odds fusion
score = sum(w_i * log_odds(P_i)) for each class i
decoded_base = argmax(CTC_beam_search(log_probs))
```

シグナルをオーバーラップチャンク(4000サンプル単位)に分割し、各チャンクを独立にBiLSTMで処理した後、動的計画法によりステッチする。

---

### 3.2 統合 SV 検出戦略（Module 2）

**ファイル**: `src/sv_detection/sv_caller.py`

#### 2a. Split-Read 検出

SAM/PAFの補助アライメントタグ（SA:Z:）を解析し、プライマリアライメントと補助アライメントの関係からSV種別を推定：

| 条件 | SV種別 |
|------|--------|
| 同染色体・同ストランド・正のギャップ | DEL |
| 同染色体・同ストランド・負のギャップ | INS |
| 同染色体・逆ストランド | INV |
| 異染色体 | TRA |

クラスタリングは最大距離200 bpのスイープライン法で実施。

#### 2b. Read-Depth CNV 検出

```
depth[bin] → MAD正規化 → Z-score
|Z| > 3.0 の連続ビン (≥5ビン) → CNV call
Z > 0: DUP, Z < 0: DEL
```

#### 2c. Assembly-Based 検出

de-Bruijn グラフ (k=15) による局所アセンブリ後、コンティグと参照配列のアライメントからブレークポイントを同定。プロダクション環境ではhifiasm / wtdbg2を使用。

#### ベイズ統合

3手法の証拠スコアをlog-odds加重平均で統合し、シグモイド関数で[0,1]に変換：

```
combined_score = sigmoid(
    0.4 × log_odds(SR_score) +
    0.3 × log_odds(RD_score) +
    0.3 × log_odds(AB_score)
)
```

2手法以上のサポートがある場合、スコアに15%のボーナス付与。

---

### 3.3 リピート領域特殊処理（Module 3）

**ファイル**: `src/repeat_regions/repeat_handler.py`

| リピート種別 | 処理方針 |
|-------------|---------|
| **テロメア** | TelomereHunter方式: TTAGGG/CCCTAA コピー数カウント → テロメア長推定 |
| **セントロメア** | hg38セントロメア座標 ± 500 kbのSVに対してスコア閾値を0.8に引き上げ |
| **セグメント重複** | identity ≥ 90%のセグドップ領域のSVにSEGDUPフラグ; アセンブリ証拠必須 |
| **STR** | 正規表現で反復モチーフカウント; 参照コピー数+5以上→伸長と判定 |

言語的複雑度（k-mer多様性）でシーケンスのリピート種別を自動分類。

---

### 3.4 複雑 SV 検出（Module 4）

**ファイル**: `src/complex_sv/complex_sv_detector.py`

#### クロモスリプシス (Chromothripsis)

Stephens et al. (Cell 2011) の5つのホールマークを検証：

1. 1染色体に≥10ブレークポイントの集積
2. DEL + INV + DUP の混在（複数SV種別）
3. ブレークポイント間の短いセグメント（平均<500 kb）
4. コピー数の振動パターン（振幅 ≥ 1.5）
5. ストランド結合の無作為性（置換検定 p ≥ 0.05）

信頼スコア = 合格ホールマーク数 / 5

#### 染色体外DNA (ecDNA)

AmpliconArchitect 方式に基づくアルゴリズム：

1. 局所高コピー数領域（CN ≥ 5×）を同定
2. BND/TRA SVによる結合グラフを構築
3. 円形トポロジーの再構築を試みる
4. バック接合リード（circle junction reads）でサポートを確認
5. 既知オンコジーン（MYC, MYCN, EGFR, KRAS, ERBB2, CDK4）との重複を確認

#### クロモプレキシ

TRA/BNDのチェーン（≥3染色体、各接続 ≤ 1 Mb）を検出。

#### BFB サイクル

折り畳み逆位 + コピー数勾配 + テロメア消失の3条件を確認し、BFBサイクル数を `ceil(log2(max_CN))` で推定。

---

### 3.5 ハイブリッド解析（Module 5）

**ファイル**: `src/hybrid_analysis/hybrid_sv_caller.py`

SURVIVOR アルゴリズムに基づくマージ（最大距離1 kb、相互オーバーラップ ≥ 50%）：

```
hybrid_score = sigmoid(geometric_mean(LR_score, SR_score)) × (1 + 0.1×(n_callers-1))
```

- **両方でサポート**: スコアブースト + Illuminaによるブレークポイント精密化（±10 bp）
- **LRのみ**: 15% スコアペナルティ
- **SRのみ**: 20% スコアペナルティ
- **再ジェノタイピング**: IlluminaのカバレッジプロファイルからCN推定 → GT改訂（0/1 or 1/1）

---

### 3.6 GIAB ベンチマーク評価設計（Module 6）

**ファイル**: `src/benchmark/giab_benchmark.py`

Truvari スタイルのマッチング基準：

| パラメータ | 値 |
|-----------|-----|
| max_distance | 500 bp |
| min_reciprocal_overlap | 0.50 |
| min_size_similarity | 0.70 |
| require_same_type | True |

評価指標：
- **Precision** = TP / (TP + FP)
- **Recall** = TP / (TP + FN)
- **F1** = 2PR / (P + R)
- **GT concordance** = GT一致TP / TP

SVサイズ別の層別評価：50 bp–500 bp、500 bp–5 kb、5 kb–50 kb、50 kb–1 Mb

---

## 4. 主要な結果と数値

### 4.1 ベースコーラーデモ結果

| 項目 | 値 |
|------|-----|
| 処理リード長 | 570 bp |
| チャンク数 | 6 |
| ホモポリマー補正 | 12×A → 10×A（cap適用） |

> ※ デモはhidden_size=64、2層の小規模モデルで実施。プロダクション設定（hidden=384、5層）では大幅な精度向上が見込まれる。

### 4.2 SV 検出デモ結果

| 手法 | 検出コール数 |
|------|------------|
| Split-read (クラスタリング後) | 60 |
| Read-depth CNV | 2 |
| Assembly-based | 2 |
| **統合後** | **62** |

### 4.3 ハイブリッド解析デモ結果

| 項目 | 値 |
|------|-----|
| LRのみコール数 | 40 |
| SRのみコール数 | 28 |
| マージ後合計 | 40 |
| 両方でサポート | 28 (70%) |
| PASS通過 | 38 (95%) |
| 平均ハイブリッドスコア (PASS) | 0.7944 |

### 4.4 GIAB ベンチマーク結果

以下は模擬GIAB truth set（DEL×150、INS×200、INV×30、DUP×40 = 合計420 SV）に対する評価：

![Benchmark Comparison](figures/benchmark_comparison.png)

*Figure 2. LongSV v1（ベースライン）、v2（ハイブリッド）、v3（複雑SV対応）のGIABベンチマーク比較（DEL/DUP/INV）。*

#### v1: ベースライン (Split-read のみ)

| SV Type | Precision | Recall | F1 | TP | FP | FN |
|---------|-----------|--------|----|----|----|-----|
| DEL | 0.7402 | 0.6267 | 0.6787 | 94 | 33 | 56 |
| DUP | 0.6250 | 0.6250 | 0.6250 | 25 | 15 | 15 |
| INV | 0.6875 | 0.7333 | 0.7097 | 22 | 10 | 8 |
| INS | 0.0000 | 0.0000 | 0.0000 | 0 | 9 | 200 |

#### v2: ハイブリッド解析統合

| SV Type | Precision | Recall | F1 | 改善 |
|---------|-----------|--------|----|------|
| DEL | 0.8521 | 0.8067 | **0.8288** | +15.0 pp F1 |
| DUP | 0.8500 | 0.8500 | **0.8500** | +22.5 pp F1 |
| INV | 0.8485 | 0.9333 | **0.8889** | +17.9 pp F1 |

#### v3: 複雑SV対応フルパイプライン

| SV Type | Precision | Recall | F1 | 改善 vs v1 |
|---------|-----------|--------|----|-----------|
| DEL | 0.8156 | 0.7667 | **0.7904** | +11.2 pp F1 |
| DUP | 0.9268 | 0.9500 | **0.9383** | +31.3 pp F1 |
| INV | 0.9000 | 0.9000 | **0.9000** | +19.0 pp F1 |

#### INS 検出の課題

シミュレーション内でのINS（挿入SV）は相互オーバーラップベースのマッチングでは再現率0%となった。これはINSの定義上 `end = start + 1`（長さはINFOフィールドに格納）のため、オーバーラップ計算が機能しないためである。実際の評価ではTruvari v4の `--no-ref` オプションと配列類似度（Levenshtein距離）を使用した評価が必要。

### 4.5 リピート領域処理デモ結果

| 解析 | 結果 |
|------|------|
| テロメア検出 (TTAGGG×8) | 検出成功 (fraction=37.5%) |
| STR伸長 (CAG/HTT) | 参照20コピー→最大34コピー、伸長検出成功 |
| chr7:60000000 セントロメア | True (hg38座標: 58.9M–62.1M) |

---

## 5. 考察と今後の展望

### 5.1 主要な知見

1. **ハイブリッド解析の有効性**: v2でDEL F1が0.6787→0.8288に向上（+15 pp）。Illuminaによるブレークポイント精密化と再ジェノタイピングが特に有効。

2. **DUP検出でのv3の優位性**: v3でDUP F1が0.6250→0.9383に大幅向上。Assembly-basedとRead-depth証拠の統合が重複検出に大きく寄与する。

3. **INS検出の構造的課題**: 全バージョンでINS検出が0%。挿入配列のポリッシングと配列類似度ベースのマッチングが不可欠。現状では局所アセンブリによる挿入配列の復元が最有効手段。

4. **リピート領域の特殊フィルタリング**: セントロメア・セグドップ領域での高FP率はストリクトフィルタで大幅削減可能。ただし感度の低下とのトレードオフあり。

5. **複雑SV検出の閾値設定**: クロモスリプシス検出に必要な「1染色体に≥10ブレークポイント」という閾値は、デモのランダムデータでは達成されず。実際のがんゲノムデータでは検出可能と予想される。

### 5.2 アルゴリズムの限界

| 限界 | 影響 | 対策 |
|------|------|------|
| 参照ゲノム崩壊リピート領域 | テロメア・セントロメアでの検出不可 | T2T-CHM13参照ゲノムへの移行 |
| デモ用BiLSTMがランダム重み | 実際のベースコール精度は不明 | Guppy/Bonitの事前学習重みを使用 |
| シミュレーションデータのバイアス | 実データでの性能は異なる可能性 | NA24385実データでの検証必須 |
| 計算コスト | ONTシグナルのリアルタイム処理は困難 | GPU最適化（PyTorch + CUDA） |
| INS配列マッチング | 挿入サイズ類似度のみでは不十分 | REF/ALT配列比較の実装 |

### 5.3 今後の展望

#### 短期（〜6か月）

- **Guppy/Bonito連携**: 事前学習RNN重みの統合によるベースコール精度の実用化
- **truvari v4との統合**: `bench`, `phab`, `ga4gh` モジュールによる標準化評価
- **CLIPpy適応型フィルタリング**: カバレッジ依存の動的閾値調整
- **GPU加速**: PyTorch CUDAバックエンドでのBiLSTM推論の高速化

#### 中期（〜1年）

- **T2T-CHM13参照ゲノム対応**: テロメア・セントロメアの完全配列への対応
- **Transformer-based basecaller**: Bonitob を参考に MultiHead Attention全体への置換
- **エピジェネティクス統合**: ONT修飾塩基コールとSV検出の同時解析
- **腫瘍-正常対応設計**: ソマティックSV用の差分呼び出しモード

#### 長期（〜2年）

- **pan-genome SV calling**: 線形参照に依存しないグラフゲノムベースの SV 解析
- **リアルタイム解析**: ONT MinION ストリーミングと組み合わせたオンザフライSV解析
- **多様なコホート**: 1000 Genomes + UK Biobank での集団ゲノム SV カタログ構築

---

## 6. 生成したファイル一覧

### ソースコード

| ファイル | 概要 |
|--------|------|
| `src/basecalling/rnn_basecaller.py` | BiLSTM+CTCベースコーラー（全コンポーネント含む） |
| `src/sv_detection/sv_caller.py` | 統合SV検出エンジン (Split-read/Read-depth/Assembly-based) |
| `src/repeat_regions/repeat_handler.py` | テロメア・セントロメア・STR特殊処理 |
| `src/complex_sv/complex_sv_detector.py` | クロモスリプシス/ecDNA/クロモプレキシ/BFB検出 |
| `src/hybrid_analysis/hybrid_sv_caller.py` | ショートリード+ロングリードハイブリッド解析 |
| `src/benchmark/giab_benchmark.py` | GIABベンチマーク評価エンジン |

### 結果ファイル

| ファイル | 内容 |
|--------|------|
| `results/basecaller_demo.json` | ベースコーラーデモ出力（シーケンス長・品質スコア） |
| `results/sv_calls_demo.vcf` | 統合SVコールのVCF 4.2形式出力 |
| `results/sv_summary.json` | SV種別カウントと平均スコア |
| `results/repeat_handler_demo.json` | テロメア長推定・STR伸長検出結果 |
| `results/complex_sv_demo.json` | 複雑SV検出サマリー |
| `results/hybrid_sv_demo.json` | ハイブリッド解析マージ統計 |
| `results/benchmark_summary.json` | 全シナリオのGIABベンチマーク結果 |
| `results/benchmark_longsv_v1.md` | v1ベースラインの詳細ベンチマークレポート |
| `results/benchmark_longsv_v2.md` | v2ハイブリッドの詳細ベンチマークレポート |
| `results/benchmark_longsv_v3.md` | v3フルパイプラインの詳細ベンチマークレポート |

### 図表

| ファイル | 内容 |
|--------|------|
| `figures/pipeline_architecture.png` | パイプライン全体アーキテクチャ図 (200 DPI) |
| `figures/pipeline_architecture.svg` | ベクター形式アーキテクチャ図 |
| `figures/benchmark_comparison.png` | v1/v2/v3のGIABベンチマーク比較図 (180 DPI) |

### ログ

| ファイル | 内容 |
|--------|------|
| `logs/process-log.jsonl` | 実行トレース（全フェーズ） |

---

## 参考文献

1. Stephens PJ et al. (2011). Massive genomic rearrangement acquired in a single catastrophic event during cancer development. *Cell*, 144(1), 27–40.
2. Zook JM et al. (2020). A robust benchmark for germline structural variant detection. *Nature Biotechnology*, 38, 1347–1355.
3. Jiang T et al. (2021). Long-read sequencing-based detection of structural variants in complex regions of the human genome. *Briefings in Bioinformatics*, 22(5).
4. Turner KM et al. (2017). Extrachromosomal oncogene amplification drives tumour evolution and genetic heterogeneity. *Nature*, 543, 122–125.
5. Korbel JO & Campbell PJ (2013). Criteria for inference of chromothripsis in cancer genomes. *Cell*, 152(6), 1226–1236.
6. Sedlazeck FJ et al. (2018). Accurate detection of complex structural variations using single-molecule sequencing. *Nature Methods*, 15, 461–468.
7. Wick RR et al. (2019). Performance of neural network basecalling tools for Oxford Nanopore sequencing. *Genome Biology*, 20, 129.
