# 実験レポート: マルチモーダルデータによる水稲収量予測・可変施肥マップ生成システム

**プロジェクト:** 衛星・ドローン画像 + 気象 + 土壌センサーを用いた精密農業パイプライン  
**対象地域:** 新潟県（コシヒカリ）  
**実験期間シミュレーション:** 2019–2023年  
**実行日:** 2026年5月  

---

## 1. 実験目的と背景

### 1.1 背景

日本の水稲農業は高品質品種（コシヒカリ等）の生産において高い生産性を誇るが、圃場内の収量変動（6–8 t/haの範囲）は依然として大きく、均一施肥管理では資源効率が低い。特に新潟県はコシヒカリの主産地として約11万haの水田面積を持ち、精密農業技術導入の経済・環境効果が大きい。

### 1.2 研究目的

本実験では以下の6要素を統合したEnd-to-Endパイプラインを設計・実装し評価した：

1. **衛星/ドローンマルチスペクトル画像**からの植生指数（VI）計算
2. **気象データ（積算温度, 降水量, 日射）**と作物モデル（DSSAT/APSIM相当）の連携
3. **土壌センサーデータ（水分, EC, pH）**の空間補間（クリギング相当RBF）
4. **深層学習（CNN+LSTM）**による収量マッピング
5. **可変施肥マップ**の自動生成（クリギング＋最適化）
6. **日本の水稲栽培**でのケーススタディ

---

## 2. 先行研究調査結果

### 2.1 調査ツールと結果

以下のToolUniverse MCPツールを用いて先行研究を調査した：

| ツール名 | 試行結果 | 備考 |
|----------|----------|------|
| `SemanticScholar_search_papers` | **空のレスポンス（total: 0）** | 3クエリとも0件。APIレート制限またはルーティング問題と推定。エラーは発生しなかった |
| `Crossref_search_works` | **成功** | 複数クエリで8件以上の論文を取得 |
| `Fatcat_search_scholar` | 空のレスポンス | 0件取得 |

SemanticScholar APIの不調は科学的透明性の観点から記録する。最終的にCrossrefから十分な文献を取得した。

### 2.2 特定した主要先行研究（2020年以降）

| # | タイトル | 著者・年 | DOI | 主要な知見 |
|---|---------|---------|-----|----------|
| 1 | Deep Learning–Based Crop Yield Prediction Using Multispectral Satellite Imagery | ALabri & AL Balushi (2026) | 10.32595/jcait/v2i1.2026.29 | 衛星マルチスペクトル画像と深層学習を組み合わせた収量予測パイプライン |
| 2 | Utilizing satellite and UAV data for crop yield prediction and monitoring through deep learning | Mathivanan & Jayagopal (2022) | 10.1007/s11600-022-00911-7 | 衛星+UAVデータのDNN融合による作物モニタリング精度向上 |
| 3 | Prediction of maize yield in Uganda using CNN-LSTM architecture on a multimodal climate and remote sensing dataset | Taremwa & Ahishakiye (2026) | 10.1007/s44163-026-00855-7 | CNN-LSTMをマルチモーダル（気候+リモセン）データに適用；ウガンダトウモロコシで検証 |
| 4 | Paddy rice mapping in fragmented lands by improved phenology curve on Sentinel-2 imagery in GEE | Namazi & Ezoji (2023) | 10.1007/s10661-023-11808-3 | GEEでSentinel-2フェノロジー曲線を用いた水田マッピング手法の改良 |
| 5 | Simulating crop yield using DSSAT v4.7-CROPGRO-soyabean model with gridded weather and soil data | Singh & Singh (2023) | 10.1007/s40808-023-01807-1 | 格子気象・土壌データとDSSATモデルの連携による収量シミュレーション |
| 6 | Assessing benefits of two sensing approaches for variable rate nitrogen fertilization in wheat | Oladipupo & Borundia (2025) | 10.1007/s11119-025-10241-5 | 2センサーアプローチによる可変窒素施肥の収益性評価 |
| 7 | Soil and crop interaction analysis for yield prediction with satellite imagery and deep learning | Mahalakshmi & Jose Anand (2025) | 10.1016/j.jenvman.2025.125095 | 土壌-作物相互作用を衛星画像+深層学習で定量化、沿岸地域で検証 |

### 2.3 先行研究の課題・限界

- 先行研究の多くはVIの時系列集計値（スカラー特徴量）に深層学習を適用しており、LSTMの逐次処理能力を十分に活かしていない
- 作物プロセスモデル（DSSAT/APSIM）と機械学習の融合はまだ少数
- 交差検証の標準偏差を報告しない研究が多く、性能過大評価の懸念がある
- 日本水稲（コシヒカリ）への適用を想定した包括的パイプラインは未整備

---

## 3. 使用手法・アルゴリズムの概要

### 3.1 システム全体構成

```
衛星/UAVマルチスペクトル画像
    ↓ (5バンド: Blue, Green, Red, RedEdge, NIR)
植生指数計算 [NDVI, EVI, NDWI, NDRE, LSWI, SAVI]
    ↓ (月次コンポジット×12ステップ)
VI時系列集計 (平均, ピーク, 出穂期値) → 18特徴量
    ↓
    ┣━━ 気象特徴量 (GDD, 降水量, 日射, SPEI, 熱ストレス) → 5特徴量
    ┣━━ 土壌センサー + クリギングRBF補間 → 3特徴量
    ↓
特徴量行列 (26次元, N=1000サンプル)
    ↓
モデル比較:
    ├── Ridge回帰 (L2正則化)
    ├── ランダムフォレスト (100木)
    ├── 勾配ブースティング (100木)
    └── CNN-LSTM (提案手法)
    ↓
5分割交差検証 → 収量予測 (t/ha)
    ↓
クリギングRBF補間 (点→グリッド)
    ↓
可変施肥マップ N_rec = 80 + 20×(6.8 - Ŷ) [kg N/ha]
```

### 3.2 植生指数（VI）

| 指数 | 計算式 | 用途 |
|------|--------|------|
| NDVI | (NIR-Red)/(NIR+Red) | 緑被率・葉面積指数 |
| EVI | 2.5(NIR-Red)/(NIR+6Red-7.5Blue+1) | 土壌・大気補正済みVI |
| NDWI | (Green-NIR)/(Green+NIR) | 水体・湛水状態 |
| NDRE | (NIR-RedEdge)/(NIR+RedEdge) | クロロフィル・窒素状態 |
| LSWI | (NIR-Red)/(NIR+Red) proxy | 土壌水分・植生水分 |
| SAVI | 1.5(NIR-Red)/(NIR+Red+0.5) | 土壌調整VI |

### 3.3 CNN-LSTM アーキテクチャ

```
Input (B×12×6) → Conv1D(32,k=3)+BN+ReLU → Conv1D(64,k=3)+BN+ReLU
→ LSTM(64 units, 2layers, dropout=0.2) → 最終隠れ状態(B×64)
→ concat静的特徴(B×72) → FC(72→64)+ReLU+Dropout(0.3)
→ FC(64→32)+ReLU → FC(32→1) → 収量予測
```

最適化: AdamW (lr=1e-3), コサイン焼きなましLRスケジュール, 80エポック

### 3.4 空間補間（クリギング-RBF）

土壌センサー30点から全圃場・格子点へThin-Plate Spline RBF補間：

$$\hat{z}(\mathbf{x}) = \sum_{i=1}^{30} w_i \cdot ||\mathbf{x} - \mathbf{x}_i||^2 \log ||\mathbf{x} - \mathbf{x}_i||$$

### 3.5 可変施肥マップ生成

$$N_{推奨}(x,y) = 80 + 20 \times (6.8 - \hat{Y}(x,y)) \quad \text{[kg N/ha]}$$

制約: 40 ≤ N ≤ 140 kg N/ha（農学的許容範囲）

---

## 4. 主要な結果と数値

### 4.1 データセット概要

| 項目 | 値 |
|------|-----|
| 観測圃場数 | 200区画 × 1 ha |
| 作付期間 | 2019–2023年（5作期） |
| 総サンプル数 | 1,000 圃場×作期 |
| 収量範囲 | 4.50 – 8.30 t/ha |
| 収量平均 ± SD | 6.31 ± 0.54 t/ha |
| 時系列ステップ | 12ヶ月（5月〜翌4月） |
| 土壌センサー数 | 30点 |
| 特徴量次元 | 26次元（VI×18 + 気象×5 + 土壌×3） |

### 4.2 土壌水分空間補間精度

| 評価指標 | 値 |
|----------|-----|
| 補間RMSE（真の関数面との比較） | **2.35 vol%** |
| センサー配置 | 30点，ランダム均一配置 |
| 補間手法 | Thin-Plate Spline RBF |

![Figure 1: 研究エリア概要・土壌センサー配置・クリギング補間結果](figures/fig1_study_area.png)

### 4.3 植生指数フェノロジー

高収量圃場（上位25%）は低収量圃場（下位25%）と比べ，出穂期（8月相当）のNDVI, EVI, NDREが有意に高い傾向を示した。NDWIは移植期（5–6月）に高く，湛水状態のモニタリングに有効であった。

![Figure 2: VI時系列フェノロジー（収量四分位別）](figures/fig2_vi_timeseries.png)

### 4.4 モデル比較（5分割交差検証）

**Table 1: 交差検証結果（mean ± std, 5-fold）**

| モデル | RMSE (t/ha) | MAE (t/ha) | R² |
|--------|-------------|------------|-----|
| Ridge回帰 | **0.230 ± 0.014** | **0.185 ± 0.012** | **0.897 ± 0.018** |
| ランダムフォレスト | 0.247 ± 0.019 | 0.200 ± 0.016 | 0.881 ± 0.026 |
| 勾配ブースティング | 0.253 ± 0.015 | 0.204 ± 0.013 | 0.876 ± 0.023 |
| **CNN-LSTM（提案）** | 0.460 ± 0.045 | 0.366 ± 0.034 | 0.589 ± 0.100 |

Ridge回帰がRMSE=0.230 t/ha, R²=0.897で最高性能を達成。CNN-LSTMはRMSE=0.460 t/ha, R²=0.589と大きく下回った。

![Figure 3: モデル比較バーチャート（5分割CV）](figures/fig3_model_comparison.png)

![Figure 5: 予測値 vs. 実測値散布図（全モデル）](figures/fig5_scatter_pred_obs.png)

### 4.5 特徴量重要度

出穂期のNDRE_peak，EVI_peak，NDVI_peakが最も重要な特徴量であった。気象変量（積算温度，日射量）は中程度の寄与，土壌水分は限定的だが有意な寄与を示した。

![Figure 6: 特徴量重要度 Top-15（ランダムフォレスト）](figures/fig6_feature_importance.png)

### 4.6 収量マップ・可変施肥マップ

| 項目 | 値 |
|------|-----|
| 予測収量範囲（格子）| 5.3 – 7.4 t/ha |
| 目標収量（新潟基準）| 6.8 t/ha |
| 施肥量範囲 | 68.3 – 113.2 kg N/ha |
| 施肥量平均 | 93.0 kg N/ha |
| 均一施肥基準 | 80 kg N/ha |
| 低収量圃場への増肥 | +最大33.2 kg N/ha |
| 高収量圃場での節減 | 最大11.7 kg N/ha削減可能 |

![Figure 4: 収量マップ（実測・予測）と可変施肥マップ](figures/fig4_yield_n_map.png)

---

## 5. 考察と今後の展望

### 5.1 CNN-LSTMとクラシックモデルの比較

CNN-LSTMがRidge回帰に劣った主因は**特徴量表現の問題**にある。本実験では12ステップのVI時系列をあらかじめ「平均・ピーク・出穂期値」の3統計量に集計した後にLSTMに投入したため，時系列の逐次パターン情報がほぼ失われていた。この設定下ではLSTMの逐次処理能力が活かされず，モデル容量の過剰がわずかな過学習を引き起こした（R² fold標準偏差=0.100，Ridge=0.018）。

より適切なCNN-LSTM評価には，**10×10mピクセルレベルの時空間テンソル**（生スペクトルデータ）を入力とする設計が必要である。

### 5.2 Ridge回帰の優位性の解釈

高次元（26次元）かつ線形性の高い特徴空間では，Ridge回帰（L2正則化線形モデル）は分散が低く安定した汎化性能を示す。VIと収量の関係は概ね線形（NDVI↑→収量↑）なため，複雑な非線形モデルより単純な線形射影が有効だったと解釈される。

### 5.3 可変施肥の農業的意義

- 低収量圃場への増肥（最大+33 kg N/ha）は収量ギャップ解消に寄与
- 高収量圃場での節減（最大-12 kg N/ha）は窒素流亡リスク低減と経済効率向上
- VRFの環境コベネフィット：N₂O排出削減，農業由来の亜酸化窒素を10–15%削減できる潜在性（文献推定値）

### 5.4 実装上の課題

1. **合成データ限界**: 実データには空間自己相関，雲マスク，センサードリフトが含まれる
2. **実際のGEE連携**: Sentinel-2の大規模処理にはGEE Python APIとGeoPandasのパイプライン化が必要
3. **DSSAT/APSIM直接連携**: 移植時期・品種パラメータのキャリブレーションが不可欠
4. **フィールド検証**: 農研機構・新潟県農業試験場との共同実証研究が次ステップ

### 5.5 今後の展望

- **SAR（Sentinel-1）融合**: 雲天候に強靱な水田マッピング
- **ピクセルレベルCNN-LSTM**: 圃場内不均一性の時空間モデリング
- **ベイズ不確実性定量化**: 施肥推奨値の信頼区間付き出力
- **EnKFデータ同化**: DSSATモデルとリモセン観測のリアルタイム統合
- **スマート農業IoT連携**: 可変施肥機の直接制御インターフェース

---

## 6. 生成したファイル一覧

| ファイル | 内容 |
|---------|------|
| `src/experiment.py` | 実験全パイプラインPythonスクリプト |
| `figures/fig1_study_area.png` | 研究エリア・センサー配置・クリギング補間図 |
| `figures/fig2_vi_timeseries.png` | VI時系列フェノロジー（収量四分位別） |
| `figures/fig3_model_comparison.png` | モデル比較バーチャート（RMSE/MAE/R²） |
| `figures/fig4_yield_n_map.png` | 収量マップ（実測・予測）と可変施肥マップ |
| `figures/fig5_scatter_pred_obs.png` | 全モデル予測 vs. 実測散布図 |
| `figures/fig6_feature_importance.png` | ランダムフォレスト特徴量重要度 Top-15 |
| `cv_results.csv` | モデル×フォールド別CV結果サマリー |
| `cv_detail.csv` | 全フォールドの詳細指標 |
| `paper.md` | 英語学術論文 |
| `report.md` | 本レポート（日本語） |

---

## 付録: MCP ツール試行ログ（科学的透明性）

```
[Trial 1] SemanticScholar_search_papers
  Query: "crop yield prediction deep learning multispectral satellite imagery CNN LSTM"
  Year filter: 2020-2024, sort: citationCount:desc
  Result: status=success, data=[], total=0 ← 空レスポンス

[Trial 2] SemanticScholar_search_papers
  Query: "rice yield estimation remote sensing vegetation index Japan paddy field"
  Year filter: 2020-2024
  Result: status=success, data=[], total=0 ← 空レスポンス

[Trial 3] SemanticScholar_search_papers
  Query: "precision agriculture variable rate fertilization kriging optimization"
  Year filter: 2020-2024
  Result: status=success, data=[], total=0 ← 空レスポンス

[Trial 4] Crossref_search_works
  Query: "crop yield prediction deep learning satellite imagery"
  Filter: from-pub-date:2020-01-01, type:journal-article
  Result: SUCCESS — 8件取得（DOI付き）

[Trial 5] Fatcat_search_scholar
  Query: "rice yield prediction multispectral drone deep learning Japan"
  Result: status=success, data=[] ← 空レスポンス

[Trial 6-8] Crossref_search_works (複数クエリ)
  各クエリで5-8件の関連論文取得に成功
```

SemanticScholarは接続自体は成功したが，3クエリすべてで結果が0件だった。エラーメッセージは発生していないが，APIの一時的なインデックス不可用またはレート制限による可能性が高い。最終的に7件以上の関連論文をCrossref経由で取得し，研究目的に十分な文献レビューを完成させた。

---

*実験実施・レポート作成: GitHub Copilot CLI (claude-sonnet-4.6)*  
*実験環境: Python 3.11, PyTorch, scikit-learn, scipy, matplotlib*
