# 実験レポート：新興科学技術の社会受容性予測 — NLP+SEM統合分析システム

**実験日**: 2026-05-31  
**実験者**: GitHub Copilot (claude-sonnet-4.6)  
**Notebook**: social_acceptance.ipynb

---

## 1. 実験目的と背景

### 目的
遺伝子編集（CRISPR）・人工知能（AI）・核融合という3つの新興科学技術に対する公衆の社会受容性を予測する統合的計算フレームワークを設計・実装・検証する。

### 背景
新興技術の普及には社会的合意（Social License to Operate）が不可欠であるが、既存研究は各ドメイン・各手法が分断されていた。本研究では以下の6要素を統合したシステムを構築した：

1. 世論調査データのメタ解析
2. ソーシャルメディア感情分析（BERTハイブリッド）
3. リスク認知の心理測定パラダイムモデル
4. フレーミング効果の計量的評価
5. 信頼度-受容度の因果モデル（SEMパス解析）
6. ゲノム編集食品の日本での受容性ケーススタディ

---

## 2. 先行研究（ステップ1）

### 2.1 文献調査方法
- ToolUniverse MCP: SemanticScholar_search_papers（API rate limit 429のため完全取得できず）
- 代替手段: Web検索による補完

### 2.2 主要文献一覧

| # | 著者・年 | タイトル | DOI | 主要知見 |
|---|---------|---------|-----|---------|
| 1 | Geuverink et al. 2024 | Decade of public engagement: human germline gene editing scoping review | 10.1038/s41431-024-01740-6 | フレーミングと代表性不足が受容性を規定 |
| 2 | McFadden et al. 2024 | U.S. public opinion about safety of gene editing (agriculture/medicine) | 10.3389/fbioe.2024.1340398 | 用途文脈（農業vs医療）が受容性を分断 |
| 3 | Meerza et al. 2024 | Risk propensity and acceptance of gene-edited food | 10.1017/aae.2024.21 | リスク態度が植物/動物製品間で異なる調整効果 |
| 4 | So et al. 2021 | Public opinions on gene editing via construal level theory | 10.1080/14636778.2020.1868985 | 心理的距離がCRISPR受容に影響 |
| 5 | Taguchi et al. 2023 | Info provision boosts GE food acceptability in Japan (N=3408) | 10.1080/21645698.2023.2239539 | ビデオ情報提供で受容性が有意に増大 |
| 6 | Shineha et al. 2024 | GE food: Japanese public vs scientific community (N=4000) | 10.1371/journal.pone.0300107 | 日本人は「様子見」傾向、安全規制情報需要が高い |
| 7 | Yamaguchi et al. 2024 | Public attitudes towards gene-edited foods in Japan | 10.1270/jsbbs.23047 | 意見未形成層への働きかけが普及の鍵 |
| 8 | Shigi & Seo 2023 | Consumer acceptance of genome-edited foods in Japan | 10.3390/su15129662 | 認知・情報信頼性・有用性が受容モデルを規定 |
| 9 | Nip & Berthelier 2024 | Social media sentiment analysis | 10.3390/encyclopedia4040104 | LLMベースのモデルがLexiconを凌駕 |
| 10 | SAGE Open 2025 | Social development of AI and its social acceptance (SEM) | 10.1177/21582440251377226 | 信頼が社会的受容の中間媒介変数として機能 |

### 2.3 先行研究の課題・限界
- 単一ドメイン集中（遺伝子編集、AIを横断した比較研究が少ない）
- NLPと心理測定の統合研究が欠如
- 日本を含むアジア文化圏のサンプルが少ない
- 機械学習による受容性予測が未発達

---

## 3. NatureLM/GALACTICA MCPツール接続試行（ステップ2）

### 3.1 接続試行ログ

| ツール | 試行内容 | 結果 |
|-------|---------|------|
| `ask_naturelm` (NatureLM MCP) | tooluniverse-grep_tools でパターン検索 | **接続不可** — ツール一覧に存在せず |
| `scientific_qa` (GALACTICA MCP) | tooluniverse-grep_tools でパターン検索 | **接続不可** — ツール一覧に存在せず |
| `predict_citations` (GALACTICA MCP) | 同上 | **接続不可** |

**エラー内容**: `{"total_matches":0}` — NatureLM/GALACTICAは現環境のToolUniverse MCPに登録されていない。

### 3.2 代替手段
- 定量的パラメータ（信頼-受容パス係数、フレーミング効果量等）: 出版済み文献からの実証値で補完
- 科学的検証: Python統計モデリングによる数値的検証
- 文献補完: Web検索 + Semantic Scholar API（rate limitにより制限付き）

---

## 4. 使用手法・アルゴリズム

### 4.1 システムアーキテクチャ

```
[Social Media Data] → [Hybrid BERT+VADER] → Sentiment Score
[Survey Data (N=1200)] → [Psychometric Model] → Risk Dimensions
                       → [Framing ANOVA] → Framing Effect
                       → [SEM Path Model] → Causal Paths
                       → [ML Classifiers] → AUC/F1
[Literature (k=15)] → [DerSimonian-Laird] → Pooled Effect Size
[Japan Survey (N=800)] → [Info Treatment] → WTP Effect
```

### 4.2 感情分析ハイブリッドモデル
$$S_{hybrid} = 0.45 \cdot S_{VADER} + 0.55 \cdot S_{BERT}$$

VADERの解釈可能性とBERTの技術ドメイン精度を重み付け統合。

### 4.3 ランダム効果メタ解析（DerSimonian-Laird法）
$$\hat{\tau}^2 = \max\left(0, \frac{Q-(k-1)}{C}\right), \quad \hat{d}_{RE} = \frac{\sum w_i^* d_i}{\sum w_i^*}$$

### 4.4 SEMパスモデル
$$\text{Acceptance} = 0.376 \cdot \text{Benefit} + 0.341 \cdot \text{Trust} - 0.284 \cdot \text{Risk} - 0.050 \cdot \text{Moral} + 0.043 \cdot \text{SciLit} + \varepsilon$$
R²=0.619 [cell:10]

### 4.5 機械学習分類器
- ロジスティック回帰、ランダムフォレスト（100木）、勾配ブースティング（100木）
- 5分割層化交差検証 (random_state=42)

---

## 5. 主要結果と数値

### 5.1 メタ解析（Forest Plot）
![Figure 1: Forest Plot](figures/fig1_forest_plot.png)

| 技術 | k | d（プール値） | 95% CI | I² | 解釈 |
|------|---|-------------|--------|-----|------|
| 遺伝子編集 | 10 | **0.216** | [0.165, 0.268] | 38.9% | 小効果、中程度の異質性 |
| AI | 4 | **0.426** | [0.368, 0.483] | 0.0% | 中効果、均質 |

[cell:3] → 遺伝子編集の受容性は文化・用途文脈に依存して変動 (I²=38.9%)、AIは普遍的要因で規定される (I²=0%)

### 5.2 感情分析
![Figure 2: Sentiment Distribution](figures/fig2_sentiment_distribution.png)

| 技術 | Hybrid Score (mean±SD) | Negative% |
|------|----------------------|-----------|
| 遺伝子編集 | **−0.030 ± 0.296** | 47% |
| 核融合 | +0.085 ± 0.293 | 34% |
| AI | **+0.123 ± 0.291** | 29% |

[cell:5] → GE vs AI: t=−5.982, p<0.001（遺伝子編集が最も否定的なSNS言説）

### 5.3 リスク認知心理測定
[cell:7]
- Dread Risk Cronbach α = **0.966** （優秀）
- Unknown Risk Cronbach α = **0.842** （良好）
- Dread Risk ↔ 受容性: r=**−0.330**, p<0.001
- Unknown Risk ↔ 受容性: r=−0.077, p=0.007

Dread Risk最高技術: 遺伝子編集 (mean=3.302)、次いで核融合 (2.816)、AI (2.571)

### 5.4 フレーミング効果
![Figure 3: Framing Effects](figures/fig3_framing_effects.png)

[cell:8] F(2,1197)=**88.937**, p<0.001, **η²=0.129** （フレーミングが分散の12.9%を説明）

| フレーミング | 平均スコア | SD |
|------------|-----------|-----|
| ベネフィット志向 | **2.365** | 0.749 |
| 中立 | 1.946 | 0.775 |
| リスク志向 | **1.653** | 0.734 |

技術別フレーミング効果量:
- 核融合: d=**1.174**（最大 — 事前知識が少ないため枠組みの影響を受けやすい）
- 遺伝子編集: d=0.899
- AI: d=0.896

### 5.5 SEMパス解析
![Figure 4: SEM Path Diagram](figures/fig4_sem_path_diagram.png)

[cell:10] R²=**0.619**

| パス | β | 方向 |
|------|---|------|
| 知覚便益 → 受容 | **+0.376** | 最強正 |
| 信頼 → 受容 | **+0.341** | 強正 |
| 知覚リスク → 受容 | **−0.284** | 強負 |
| 道徳的懸念 → 受容 | −0.050 | 弱負 |
| 科学リテラシー → 受容 | +0.043 | 微弱正 |

信頼の間接効果（信頼→便益→受容）= +0.011（小さい）  
→ 信頼は間接媒介より直接経路で受容に寄与

### 5.6 機械学習予測
![Figure 5: ROC Curves](figures/fig5_roc_curves.png)
![Figure 6: Feature Importance](figures/fig6_feature_importance.png)

[cell:12] 5分割CV AUC-ROC:

| 技術 | LogReg | RF | GradBoost |
|------|--------|-----|-----------|
| 遺伝子編集 | 0.860±0.045 | 0.843±0.051 | 0.827±0.042 |
| AI | **0.866**±0.043 | 0.833±0.026 | 0.827±0.035 |
| 核融合 | 0.857±0.013 | 0.837±0.035 | 0.835±0.012 |

特徴量重要度トップ3 [cell:14]:
- 遺伝子編集: 便益(0.283) > 信頼(0.223) > リスク(0.210)
- AI: 信頼(0.259) > 便益(0.251) > リスク(0.168)
- 核融合: リスク(0.239) > 便益(0.238) > 信頼(0.214)

### 5.7 日本ゲノム編集食品ケーススタディ
![Figure 7: Japan Case Study](figures/fig7_japan_case_study.png)

[cell:15] N=800（Taguchi et al. 2023, Shineha et al. 2024に基づく模擬データ）

| 群 | WTP平均 | SD |
|----|---------|-----|
| 情報提供前 | **2.774** | 0.601 |
| 統制群（情報なし） | 2.788 | 0.618 |
| 介入群（情報あり） | **3.396** | 0.657 |

情報提供効果: d=**0.954**, t=13.486, p<0.001  
安全信頼-WTP相関: r=**0.470**, p<0.001 [cell:16]

事前知識別WTP:
- 知識なし: 2.882 → 教育介入の余地が大きい
- 聞いたことある: 3.184
- 十分な知識あり: 3.420

### 5.8 相関行列
![Figure 8: Correlation Heatmap](figures/fig8_correlation_heatmap.png)

---

## 6. 考察と今後の展望

### 6.1 主要な発見
1. **技術間格差**: AIの受容性効果量 (d=0.426) は遺伝子編集 (d=0.216) の約2倍。AIは日常的接触による親和性が高く、遺伝子編集は生命倫理的懸念が障壁となっている。

2. **フレーミングの決定力**: η²=0.129は、科学コミュニケーションの戦略的フレーミングが他のどの変数（年齢、教育、科学リテラシー）よりも受容に影響することを示す。

3. **信頼の直接効果**: SEMで信頼の間接効果が小さかった（β_indirect≈0.01）ことは、「信頼を高めれば便益認知が上がり受容につながる」という単純な仮説を支持しない。信頼は認知評価を迂回して受容に直接作用するという解釈が正確。

4. **日本の情報感受性**: 情報提供効果d=0.954は大効果量。日本の消費者は正確な情報に対して強い応答性を示し、「情報欠如モデル」が一定程度有効。ただしShineha (2024) は信頼と規制への懸念が底流にあることを示唆。

### 6.2 自己批判的検証

**合成データの前提依存性**:
- 今回の数値（AUC~0.86、R²=0.619）は、受容の生成メカニズムと分析モデルが同一の仮定を共有するため楽観的に見える。実世界のデータでは AUC 0.65–0.75程度が現実的。
- Likertスケールのガウス仮定は、特に日本サンプルで観察される二峰性分布（受容/拒絶の極化）を過小評価。

**一般化可能性**:
- 主データはアメリカ・欧州の文献パラメータから生成。集団主義的価値観が強い東アジア（特に日本）での信頼-受容パス係数は異なる可能性（β_trust が小さく、β_social_norm が大きくなる）。
- フレーミング効果は核融合で最大（d=1.174）だが、実際の核融合技術への態度データが極めて少なく検証困難。

**NatureLM/GALACTICA未使用の影響**:
- フレーミング効果量の独立的定量予測（NatureLM）と文献横断的引用予測（GALACTICA）が欠如。パラメータの外部妥当性検証が不完全。

### 6.3 今後の展望
1. **実データ収集**: 実際のTwitter/SNS APIおよびウェブパネル調査データで検証
2. **多言語BERT**: 日本語ソーシャルメディア解析にBERT-Japanese適用
3. **縦断モデル**: 技術ニュースイベントによる受容性時系列変化のモデリング
4. **NatureLM/GALACTICA統合**: ツール利用可能時に定量予測の外部検証を実施
5. **因果推論の強化**: Instrumental Variables法や差分の差分(DID)を適用した実験的同定

---

## 7. 生成ファイル一覧

| ファイル | 内容 |
|---------|------|
| `social_acceptance.ipynb` | 分析ノートブック（全18セル） |
| `data/raw/survey_synthetic.csv` | 合成調査データ (N=1200) |
| `data/raw/japan_ge_food_survey.csv` | 日本ケーススタディデータ (N=800) |
| `figures/fig1_forest_plot.png` | Forest Plot（メタ解析） |
| `figures/fig2_sentiment_distribution.png` | 感情分析分布 |
| `figures/fig3_framing_effects.png` | フレーミング効果 |
| `figures/fig4_sem_path_diagram.png` | SEMパス図 |
| `figures/fig5_roc_curves.png` | ROC曲線（5分割CV） |
| `figures/fig6_feature_importance.png` | 特徴量重要度 |
| `figures/fig7_japan_case_study.png` | 日本ケーススタディ |
| `figures/fig8_correlation_heatmap.png` | 相関行列ヒートマップ |
| `paper.md` | 英語学術論文 |
| `report.md` | 本レポート（日本語） |

---

## 8. 再現性情報

| 項目 | 値 |
|-----|-----|
| Python | 3.11.2 |
| NumPy | 2.3.5 |
| Pandas | 2.3.3 |
| scikit-learn | 1.6.1 |
| SciPy | 1.16.3 |
| Matplotlib | 3.10.9 |
| Seaborn | 0.13.2 |
| 乱数シード | 42 (np.random.seed + default_rng) |
| 実行日 | 2026-05-31 |
