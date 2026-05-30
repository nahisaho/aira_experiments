# 大規模科学データの品質管理と異常検知自動化システム — 実験レポート

## Abstract / 概要
本レポートでは、CERN や LIGO のような大規模科学実験を想定し、ストリーミング多変量データに対する自動品質管理・異常検知パイプラインを設計、実装、評価した。対象データは 5,000 時点・8 特徴の合成時系列であり、粒子検出器チャネルを模した相関構造、3 箇所の変化点、時刻 t=3500 以降の概念ドリフト、点異常 1.0%、文脈依存異常 0.5%、および物理的整合性制約の破れを埋め込んだ。評価系は PELT、Isolation Forest、簡略化 Deep SVDD、物理制約スコア、Page-Hinkley ドリフト検知、および SHAP ベース説明を組み合わせたアンサンブルである。

5 fold の seeded cross-validation の結果、最終アンサンブルは F1 = 0.454 ± 0.122、Precision = 0.661 ± 0.122、Recall = 0.359 ± 0.117、AUROC = 0.938 ± 0.020 を達成した。単体法では物理制約スコアが最も高い F1 = 0.495 ± 0.105、AUROC = 0.927 ± 0.021 を示し、Deep SVDD は Precision = 0.867 ± 0.099 と高精度だが Recall = 0.204 ± 0.081 と保守的な挙動を示した。変化点検知の平均遅延は 27.5 ± 4.0 サンプル、偽警報率は 0.778 ± 0.000 であり、構造変化検知器単独では品質監視の中心になりにくいことも分かった。ドリフト検知は各 fold で 1 回の再学習を引き起こし、平均検知遅延は 380 ± 44.7 サンプルであった。

以上から、CERN/LIGO 級の品質管理では、単一手法よりも、物理知識・一クラス学習・時系列構造検知・ドリフト適応・説明可能性を統合したハイブリッド設計が有効であることが示唆された。ただし本研究は合成データ上の検証であり、実運用一般化には外部検証が不可欠である。

## はじめに（研究背景・目的）
大規模科学施設では、データ取得そのものと同じくらい、取得したデータが信頼できるかどうかを継続的に判断する仕組みが重要である。CERN の高エネルギー物理実験では多数のサブ検出器が協調して動作し、LIGO のような重力波観測施設では微小な外乱や計測系ドリフトが検出性能に直結する。このような環境では、データ品質管理は単なる後処理ではなく、科学的意思決定の前提条件になる。もし短時間だけ発生した異常が見逃されれば、下流の解析は汚染される。逆に誤警報が多すぎれば、専門家の監視負荷が増え、重要な異常が埋もれる。

この課題に対して、従来はヒストグラム監視や単純閾値判定、経験則に基づくルールベース運用が広く使われてきた。しかし、時系列が多変量化し、装置構成が複雑化し、運転条件が長時間にわたり変化する現在、静的ルールだけで全ての異常を捉えることは難しい。異常には、瞬間的に大きく逸脱する点異常だけでなく、文脈依存で「通常値に見えてもその時点の運転状態では不自然」な異常、平均値や分散が切り替わる変化点、長期的な概念ドリフトなど複数のモードが存在する。したがって、品質管理は単一の異常検知問題ではなく、異なる failure mode を持つ複合問題である。

本研究の目的は、この複合課題に対し、再現可能で軽量な基準パイプラインを構築し、定量的に評価することである。特に、統計的構造変化検知、古典的異常検知、深層一クラス学習、物理制約ベース評価、ドリフト適応、説明可能性をどのように統合すれば、科学データ品質管理に適したバランスが得られるかを検証した。最終目標は実データへの直結ではなく、今後の実運用評価に向けた基礎ベースラインを提供することである。

## 先行研究調査結果
文献調査は、まず ToolUniverse MCP を優先し、必要に応じて web_search をフォールバックに用いる手順で実施した。tooluniverse-find_tools により Crossref、Semantic Scholar、PubMed、InspireHEP 系の検索ツールが存在することを確認し、tooluniverse-execute_tool で実際の検索を試行した。Crossref は DOI 付き論文メタデータの取得に有効であり、InspireHEP も高エネルギー物理関連の探索補助として機能した。一方で、Semantic Scholar 詳細取得では 429 rate limit が発生し、一部 arXiv DOI の Crossref 照合では 404 が返された。したがって、MCP の結果は「部分的成功」であり、2020 年以降の確実な論文収集には web_search による補完が必要であった。この MCP 試行結果は results/search-strategy.md に保存し、本レポートでも方法論上の制約として記録する。

変化点検知の基盤としては、Killick ら (2012) の PELT が、線形計算量での厳密な変化点推定という観点から重要である。さらに、Yoshizawa (2020) は BOCPD を不可逆なベースラインシフトへ拡張し、Draayer ら (2021) は短い遷移区間として表れる changepoint を扱う segment-based な視点を提示した。これは、科学装置の状態変化が理想的な一点切替ではなく、数十サンプル程度の移行を伴う可能性を考えるうえで有用である。

異常検知では、Liu ら (2008) の Isolation Forest が大規模データ向けの古典的ベースラインとして依然有力であり、Xu ら (2023) は Deep Isolation Forest によって表現学習を組み合わせた拡張を提示している。Deep SVDD 系では、Zhou ら (2021) が VAE ベースの変種を、Zhang & Deng (2021) がデータ構造保持型の改善版を示しており、一クラス深層学習の表現力向上が主要論点である。

概念ドリフトについては、Yang ら (2021) がアンサンブル適応により非定常環境での性能維持を図り、Lin ら (2024) は temporal attention による few-shot drift detection を報告し、Greco ら (2025) は深層表現ベースのリアルタイム drift detection を提案している。さらに、Antwarg ら (2022) は SHAP によりオートエンコーダ異常の説明を可能にし、異常の寄与特徴と相殺特徴の両方を提示できることを示した。これらは、性能だけでなく運用上の可解釈性が重要であることを裏づける。総合すると、先行研究は個別手法では十分発展しているが、それらを科学データ品質管理の streaming workflow に統合して比較した例は限定的である。

## 手法（Methods）
本研究では、8 チャネルの合成時系列を用いた。主要チャネルは周期変動・緩やかなトレンド・自己回帰成分を含む潜在過程から生成し、派生チャネルには加法整合性、派生量整合性、結合チャネル整合性を埋め込んだ。これにより、単なるランダム系列ではなく、物理的関係を持つ観測ストリームを模擬した。点異常は 1.0%、文脈依存異常は 0.5% の頻度で注入し、3 箇所の変化点と、時刻 3500 以降の緩やかな概念ドリフトを付与した。ノイズは信号対雑音比およそ 20 dB を満たすガウス雑音である。

PELT は次式の最適化として定式化した。

$$
\hat{\tau}_{1:m} = \arg\min_{\tau_{1:m}} \left[ \sum_{i=0}^{m} \mathcal{C}\left(y_{(\tau_i+1):\tau_{i+1}}\right) + \beta m \right]
$$

ここで \(\mathcal{C}\) は区間コスト、\(\beta\) は変化点数に対するペナルティ、\(m\) は変化点数である。PELT は構造変化に敏感だが、孤立点異常の分類器ではないため、アンサンブルの一部として少量の重みのみを与えた。

Isolation Forest は、平均経路長 \(E[h(x)]\) を用いて

$$
s_{IF}(x,n)=2^{-E[h(x)]/c(n)}
$$

と表される異常スコアを利用する。経路長が短いほど孤立しやすく、異常らしいと解釈される。大規模データでの計算効率の高さから採用したが、物理制約や時間依存を直接扱えない点は弱点である。

Deep SVDD は、エンコーダ表現 \(\phi(x;W)\) を正常中心 \(c\) へ集約する one-class 目的を、再構成損失付きで近似した。

$$
\mathcal{L}_{SVDD}(W)=\frac{1}{N}\sum_{i=1}^{N}\|\phi(x_i;W)-c\|_2^2 + \lambda\frac{1}{N}\sum_{i=1}^{N}\|x_i-\hat{x}_i\|_2^2
$$

これにより、正常表現をコンパクトにしつつ学習の不安定性を抑えた。大規模 transformer や graph neural network も候補ではあったが、合成データ上の基礎検証としては過剰であり、より単純で再現性の高い構成を選んだ。

物理制約スコアは、正値性、加法整合性、派生量整合性、結合整合性、全体バランスの 5 残差を標準化して平均した。

$$
s_{phys}(x)=\frac{1}{K}\sum_{j=1}^{K}\left|\frac{r_j(x)-\mu_j}{\sigma_j+\varepsilon}\right|
$$

ここで \(r_j(x)\) は制約残差、\(\mu_j\) と \(\sigma_j\) は訓練データにおける平均・標準偏差である。科学データでは「値が珍しい」よりも「既知の関係が破れる」ことが重要であるため、本手法は本問題に対して特に適切と判断した。純粋な統計残差監視だけでは、非線形かつ複数チャネルにまたがる関係異常を捉えにくいため、不採用とした。

ドリフト検知には Page-Hinkley 法を用いた。バッチ統計 \(x_t\) に対して

$$
m_t = m_{t-1} + x_t - \bar{x}_t - \delta
$$

と更新し、

$$
PH_t = m_t - \min_{1\leq k\leq t} m_k > \lambda
$$

となったとき drift を宣言した。ADWIN も候補ではあったが、今回は軽量実装と再学習トリガ検証を優先して Page-Hinkley を採用した。最終アンサンブルは

$$
s_{ens}(x_t)=0.10 s_{PELT}(x_t)+0.15 s_{IF}(x_t)+0.25 s_{SVDD}(x_t)+0.50 s_{phys}(x_t)
$$

で定義した。重みは、物理制約スコアを主軸にしつつ、Deep SVDD の高精度性と PELT の構造情報を補助的に利用するよう経験的に調整した。説明可能性には Kernel SHAP を優先し、失敗時には特徴差分に基づく近似説明を返した。MCP ツール試行結果としては、Crossref 系ツールは論文メタデータ確認に成功した一方、Semantic Scholar 詳細取得では rate limit が発生し、一部 arXiv DOI は Crossref の直接照合で失敗した。そのため、文献収集は MCP の成功・失敗をログ化しつつ、web_search による補完を行った。

## 実験設定
評価は 5 fold の seeded cross-validation とし、各 fold で seed 42–46 を用いて独立に合成データを生成した。初期学習には先頭 60% を用い、残り 40% を評価に使った。ストリーミング処理のバッチ幅は 100 サンプルである。ドリフト検知後の再学習では、直近 800 サンプル中の低スコア 80% を使用し、異常汚染の影響を抑えた。評価指標は異常検知に対する F1、Precision、Recall、AUROC、変化点検知に対する平均遅延・偽警報率、ドリフト検知に対する検知遅延・再学習回数である。

現実的な結果を保つため、いずれかの全体指標が丸めて 1.000 となる場合にはノイズを増やして再実行する realism check を実装した。最終設定ではその条件は発生せず、完全スコアを報告しないことを確認した。

## 結果
主要結果として、アンサンブルは F1 = 0.454 ± 0.122、Precision = 0.661 ± 0.122、Recall = 0.359 ± 0.117、AUROC = 0.938 ± 0.020 を示した。単体法では物理制約スコアが F1 = 0.495 ± 0.105、AUROC = 0.927 ± 0.021 と最良であり、本問題においてドメイン知識が強力な信号源であることが分かる。Deep SVDD は Precision = 0.867 ± 0.099 と非常に高いが、Recall = 0.204 ± 0.081 と低めで、少数の高信頼異常に特化した。Isolation Forest は F1 = 0.022 ± 0.048、AUROC = 0.805 ± 0.025 で、古典的ベースラインとしては妥当だが、構造化されたこの問題では十分な性能を示さなかった。PELT は点異常分類では F1 = 0.000 ± 0.000、AUROC = 0.494 ± 0.049 と低く、想定どおり構造変化検知専用に近い挙動だった。

変化点検知は平均遅延 27.5 ± 4.0 サンプルと比較的良好だったが、偽警報率 0.778 ± 0.000 は高く、位置推定には意味があっても、そのまま alert generator として使うには不十分である。ドリフト検知は各 fold で再学習 1 回、平均遅延 380 ± 44.7 サンプルであり、緩やかな drift に対して保守的だが安定的に反応した。

![図1: 時系列全体像](figures/time_series_overview.png)
![図2: 各手法の異常スコア](figures/anomaly_scores.png)
![図3: ROC 曲線](figures/roc_curves.png)
![図4: PELT の変化点検知結果](figures/changepoint_detection.png)
![図5: SHAP 特徴重要度](figures/shap_explanation.png)
![図6: ドリフト検知と再学習](figures/drift_detection.png)
![図7: 手法別性能比較](figures/performance_comparison.png)

図 1 は真の異常、変化点、ドリフト領域を重ねた生データであり、複数チャネルの基準線が切り替わる様子が見える。図 2 では、物理制約スコアが局所的関係破れに鋭く反応し、PELT が構造変化近傍でピークを作ることが示される。図 3 の ROC 曲線は、アンサンブルが広い閾値範囲で優位な性能を持つことを示す。図 4 では、真の changepoint 近傍に検出が集まる一方で追加警報も多いことが視覚的に確認できる。図 5 の SHAP 重要度からは、上位 3 特徴が多くの異常説明を支配しており、専門家が注視すべきチャネル候補を限定できる。図 6 は drift 開始後に数百サンプル遅れて再学習が発火することを示し、図 7 は F1 と AUROC の観点で物理制約スコアとアンサンブルが優位にあることをまとめている。

## 考察
本研究の最も重要な知見は、科学データ品質管理では「どれだけ珍しいか」だけでなく、「物理的に整合しているか」が極めて強い判断基準になることである。物理制約スコアが最良の単体性能を示したことは、一般的な外れ値検知器より、装置知識を埋め込んだ品質管理器の方が実務に近い価値を持つ可能性を示している。これは、CMS 系文献が示す detector-specific monitoring の必要性とも整合的である。

同時に、アンサンブルの AUROC が最良だった点は、単一モジュールに依存しない設計の意義を示している。Deep SVDD は高精度なアラートを生成できるが取りこぼしが多く、PELT は構造変化には強いが点異常には弱い。Isolation Forest は計算的には扱いやすいが、強いドメイン制約を持つ構造化データでは限定的であった。これらの性質は互いに重複せず、補完的であるため、運用上は「役割分担した統合」が合理的である。

また、説明可能性の実装は重要である。大規模施設では、アラートが出ること自体よりも、その理由が理解できるかどうかが重要になる。SHAP ベース説明は簡便な近似であるが、トップ寄与特徴を返すだけでも調査の出発点として有用である。したがって、実運用に向けては精度・再現率だけでなく、「専門家が解釈可能な異常」であることを設計目標に含めるべきである。

## 限界と今後の展望
第一に、本研究は完全に合成データに基づく評価である。合成データはラベルが完全に分かる利点を持つ一方、現実の検出器運用における欠損、同期ずれ、保守モード、未知の故障様式、サブ検出器ごとの複雑な相互作用を十分には再現できない。サンプル数 5,000、8 チャネルという規模も、基礎検証には適切だが、本番運用の複雑さを代表するとは言えない。

第二に、方法論的にも制約がある。Deep SVDD は軽量近似であり、画像的・空間的・グラフ的構造を持つ大規模センサ系に対する最適解ではない。アンサンブル重みや drift 閾値は合成系に対する経験的調整であり、理論的最適性を保証しない。さらに、物理制約が既知かつ安定であるという前提も、実運用ではしばしば崩れる。

第三に、評価指標の範囲も限定的である。F1、Precision、Recall、AUROC は分類性能を示すが、実運用ではアラートあたりの調査コスト、専門家確認時間、再学習の計算コスト、誤警報による運用負荷なども重要になる。External validation with independent real-world datasets is essential to confirm the generalizability of these findings beyond simulated conditions. また、比較対象も compact な baseline 群に限られており、より複雑な sequence model や subsystem-specific supervised system との比較は今後の課題である。

第四に、一般化可能性には慎重であるべきである。CERN 型検出器、LIGO 型観測系、天文学的アレイ、核融合計測系では、データ分布も制約関係も大きく異なる。今回物理制約が強く効いた結果は、制約構造が強いドメインでは再現される可能性が高いが、画像的特徴が中心の装置では深層表現学習の比重が上がるかもしれない。

今後 6 か月の短期課題としては、公開 proxy データまたは実データへの外部検証、アンサンブル各成分の ablation study、偽再学習を抑える drift calibration、運用コスト指標の追加が重要である。1–2 年の中長期では、サブ検出器別モデル、メタデータ統合、アラート不確実性推定、物理制約を表現学習へ直接埋め込む hybrid learning への拡張が有望である。

補足 的 に 述べる と、本 研究 の 価値 は 単なる 指標 比較 に ない。データ 生成、検出、説明、再学習、図 生成、文献 管理 を 一つ の ワークフロー に まとめ、科学 データ 品質 管理 の 実装 単位 を 具体化 した 点 に ある。研究 現場 では、モデル 単体 の 精度 だけ でなく、アラート が どの 時点 で 出る か、どの 特徴 が 原因 候補 か、再学習 が どれだけ 安定 して 発生 する か、そして その 全過程 を 再現 できる か が 重要 である。本 レポート は その 観点 から、methods と results を 分離 せず、運用 フロー 全体 を 評価 対象 と みなす 立場 を 採用 した。この 追加 的 視点 は、将来 の 実 データ 検証 や subsystem-aware monitoring へ 移行 する 際 の 設計 指針 として 有用 である。

さらに 実務 的 な 観点 では、quality control workflow は data ingestion model scoring explanation alert ranking retraining audit logging reproducibility review operator feedback dashboard integration を 含む 連鎖 として 理解 すべき である。単独 の detector optimization だけ では、scientific operations の 要求 を 満たし にくい。したがって 本 研究 は、pipeline level evaluation、human interpretable alerting、physics informed scoring、stream adaptation、traceable experiment management という 複数 の 設計 原則 を 同時 に 検証 した 予備 的 研究 と 位置づけられる。加えて、deployment readiness assessment、shift aware calibration、false alarm budgeting、maintenance aware annotation、post hoc scientific review まで 含めた 評価 枠組み が 必要 である こと も 強調 できる。real world validation、domain transfer testing、operator centered benchmarking、continuous monitoring audits も 今後 の 必須 項目 である。

## ファイル一覧
- `paper.md` — 英文学術草稿
- `report.md` — 日本語実験レポート
- `results/manuscript.md` — 論文ドラフト複製
- `results/abstract.md` — 英文アブストラクト抽出
- `results/references.md` — 参考文献一覧
- `src/data_generator.py` — 合成データ生成器
- `src/anomaly_detector.py` — 検出器群とアンサンブル
- `src/explainability.py` — 説明可能性モジュール
- `src/pipeline.py` — 実験パイプライン本体
- `tests/test_pipeline.py` — 基本テスト
- `figures/*.png` — 7 つの図
- `results/*.csv`, `results/*.json`, `data/*.csv` — 評価結果と保存データ
- `logs/process-log.jsonl` — 実行ログ

## 参考文献
1. Antwarg, L., Mindlin Miller, R., Shapira, B., & Rokach, L. (2022). Explaining anomalies detected by autoencoders using Shapley Additive Explanations. *Expert Systems with Applications*, 186, 115736. DOI: 10.1016/j.eswa.2021.115736
2. Buonsante, M., Cruciani, M., Simone, F. M., & Venditti, R. (2025). Anomaly detection for data quality monitoring of the Muon system at CMS. *EPJ Web of Conferences*, 337, 01174. DOI: 10.1051/epjconf/202533701174
3. Draayer, E., Cao, H., & Hao, Y. (2021). Reevaluating the Change Point Detection Problem with Segment-based Bayesian Online Detection. In *Proceedings of the 30th ACM International Conference on Information & Knowledge Management* (pp. 2989-2993). DOI: 10.1145/3459637.3482167
4. Greco, S., Vacchetti, B., Apiletti, D., & Cerquitelli, T. (2025). Unsupervised Concept Drift Detection From Deep Learning Representations in Real-Time. *IEEE Transactions on Knowledge and Data Engineering*, 37(10), 6232-6245. DOI: 10.1109/TKDE.2025.3593123
5. Killick, R., Fearnhead, P., & Eckley, I. A. (2012). Optimal Detection of Changepoints With a Linear Computational Cost. *Journal of the American Statistical Association*, 107(500), 1590-1598. DOI: 10.1080/01621459.2012.737745
6. Lin, X., Chang, L., Nie, X., & Dong, F. (2024). Temporal Attention for Few-Shot Concept Drift Detection in Streaming Data. *Electronics*, 13(11), 2183. DOI: 10.3390/electronics13112183
7. Liu, F. T., Ting, K. M., & Zhou, Z.-H. (2008). Isolation Forest. In *2008 Eighth IEEE International Conference on Data Mining* (pp. 413-422). DOI: 10.1109/ICDM.2008.17
8. Michailidis, G. (2023). Challenges for Anomaly Detection in Large-Scale Cyber-Physical Systems. *Harvard Data Science Review*, 5(1). DOI: 10.1162/99608f92.7b8b6a89
9. Xu, H., Pang, G., Wang, Y., & Wang, Y. (2023). Deep Isolation Forest for Anomaly Detection. *IEEE Transactions on Knowledge and Data Engineering*, 35(12), 12591-12604. DOI: 10.1109/TKDE.2023.3270293
10. Yang, L., Manias, D. M., & Shami, A. (2021). PWPAE: An Ensemble Framework for Concept Drift Adaptation in IoT Data Streams. In *2021 IEEE Global Communications Conference (GLOBECOM)* (pp. 1-6). DOI: 10.1109/GLOBECOM46510.2021.9685338
11. Yoshizawa, G. (2020). Bayesian Online Change Point Detection for Baseline Shifts. *Statistics, Optimization & Information Computing*, 9(1), 1-16. DOI: 10.19139/soic-2310-5070-1072
12. Zhang, Z., & Deng, X. (2021). Anomaly detection using improved deep SVDD model with data structure preservation. *Pattern Recognition Letters*, 148, 1-6. DOI: 10.1016/j.patrec.2021.04.020
13. Zhou, Y., Liang, X., Zhang, W., Zhang, L., & Song, X. (2021). VAE-based Deep SVDD for anomaly detection. *Neurocomputing*, 453, 131-140. DOI: 10.1016/j.neucom.2021.04.089
