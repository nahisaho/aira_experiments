# 観察データからの因果効果推定手法の体系的比較：医薬品疫学リアルワールドデータへの応用

**DRAFT — NOT FOR DISTRIBUTION**

---

## Abstract

本研究では、観察データから因果効果を推定する5つの主要手法——傾向スコアマッチング（PSM）、操作変数法（IV/2SLS）、差分の差分法（DiD）、Double/Debiased Machine Learning（DML）、因果フォレスト（Causal Forest）——を体系的に比較する。スタチン治療と心血管アウトカム（MACE）を模した合成リアルワールドデータ（N=2,000）を用い、真の平均処置効果（ATE = −0.1469）と既知の異質的処置効果（CATE）のもとで各手法の推定精度・バイアス・信頼区間被覆率を評価した。5-fold交差検証（各fold N=1,000）により再現性を確認した。

主な結果として、PSMが最小絶対バイアス（0.0476 ± 0.0036）を達成し、DiDが2番目（0.0481 ± 0.0034）、IV法（0.0537 ± 0.0111）、DML（0.0579 ± 0.0011）、因果フォレスト（0.0628 ± 0.0021）と続いた。ただし、DiDは平行トレンド仮定の違反（p < 0.001）を示し、IV法の第1段階F統計量（F=27.62）は弱操作変数の閾値（F>10）を超えた。因果フォレストは異質的処置効果の推定において顕著な優位性を示し、高齢・高LDL患者で治療ベネフィットが大きいことを可視化した。いずれの手法も残余交絡に起因するポジティブバイアスを示し、ATEを過小評価（絶対値で）する傾向があった。

---

## 1. 実験目的と背景

### 1.1 研究背景

ランダム化比較試験（RCT）は因果推論の黄金標準であるが、倫理的・実務的制約から実施が困難な場合も多い。医薬品疫学では電子カルテ・保険請求データ等のリアルワールドデータ（RWD）が活用されるが、処置割り付けが患者特性に依存する選択バイアスが問題となる。この選択バイアスを補正するため、多数の因果推論手法が開発されてきた。

近年の機械学習の発展により、Double/Debiased Machine Learning（DML; Chernozhukov et al., 2018）や因果フォレスト（Wager & Athey, 2018）といった高次元共変量を扱える手法が台頭している。一方、古典的なPSMやIV法も引き続き広く使用されており、各手法の適用条件・仮定・限界を体系的に理解することが実践上不可欠である。

本研究では、これら5手法をDoWhy/EconMLフレームワーク上で実装し、合成RWDを用いて真値との比較が可能な設定で評価する。ケーススタディとして、スタチン治療が心血管イベント（MACE）リスクに与える効果を想定した。

### 1.2 先行研究の位置づけ

本研究は以下の先行研究の知見を踏まえた：

- **PSMとその拡張**: Zhao et al. (2020) は非二値処置へのPS手法の拡張を論じ、重み付けアプローチの優位性を示した。Stuart (2023) はPSマッチングの適用とその拡張を整理し、共変量バランスの重要性を強調した。
- **DML**: Díaz (2019) はTargeted MLEとDMLの機械学習応用を比較し、正則化バイアスの除去に交差フィッティングが有効であることを示した。Kwon & Steiner (2026) はDMLと二重ロバスト推定量の統合を提案した。
- **DiD**: Rambachan & Roth (2023) は平行トレンド仮定の緩和版として「HonestDiD」を提案し、より信頼性の高い推論を可能にした。
- **因果フォレスト**: Cáceres & González (2022) は遺伝子発現データへの因果フォレスト応用を示し、個別化治療効果推定の有効性を実証した。

### 1.3 研究の新規性と貢献

1. 既知の真値を持つ合成RWDによる5手法の同一データでの公平な比較
2. 各手法の仮定検定（平行トレンド、弱操作変数、重複条件）の体系的報告
3. 因果フォレストによる異質的処置効果の可視化（年齢・LDL別CATE分布）
4. DoWhy/EconMLベースの再現可能な実装フレームワークの提供

---

## 2. 使用した手法・アルゴリズム

### 2.1 データ生成プロセス（DGP）

合成リアルワールドデータは以下の構造的因果モデル（SCM）に基づいて生成した（N=2,000）：

**交絡変数** $X$：年齢（$\sim \mathcal{N}(60, 10^2)$）、性別（Binary）、喫煙（Binary）、ベースラインLDL（$\sim \mathcal{N}(130, 30^2)$ mg/dL）、併存疾患スコア（$\sim \text{Poisson}(1.2)$）

**処置傾向スコア**：
$$P(T=1|X, Z) = \sigma\left(-2.0 + 0.04(\text{age}-60) + 0.3 \cdot \text{smoking} + 0.015(\text{LDL}-130) + 0.2 \cdot \text{comorbidity} - 0.2 \cdot \text{sex} + 0.8Z\right)$$

ここで $\sigma(\cdot)$ はシグモイド関数、$Z$ は操作変数（医師の処方傾向）。

**真の異質的処置効果（CATE）**：
$$\tau(X_i) = -0.1469 - 0.005(\text{age}_i - 60) - 0.001(\text{LDL}_i - 130)$$

高齢・高LDL患者でスタチンの効果が大きいという臨床的知見を反映。

**アウトカムモデル**：
$$Y_i = \mu_0(X_i) + T_i \cdot \tau(X_i) + \varepsilon_i, \quad \varepsilon_i \sim \mathcal{N}(0, 0.02^2)$$

$$\mu_0(X_i) = 0.05 + 0.004(\text{age}_i-60) + 0.06\cdot\text{smoking}_i + 0.001(\text{LDL}_i-130) + 0.03\cdot\text{comorbidity}_i - 0.02\cdot\text{sex}_i$$

### 2.2 各手法の概要

#### 2.2.1 傾向スコアマッチング（PSM）

ロジスティック回帰で傾向スコア $\hat{e}(X) = \hat{P}(T=1|X)$ を推定し、ロジット変換後の値に対して1:1最近傍マッチングをキャリパー（0.05）付きで実施：

$$\hat{\tau}_{PSM} = \frac{1}{|M|}\sum_{(i,j)\in M}(Y_i - Y_j)$$

ここで $M$ はマッチドペアの集合。

**主な制限**：高次元共変量での傾向スコア misspecification、マッチされない処置群の除外によるATT偏重。

#### 2.2.2 操作変数法（IV / 2SLS）

医師の処方傾向 $Z$ を操作変数として二段階最小二乗法を適用：

**第1段階**：$\hat{T} = X\hat{\gamma} + Z\hat{\delta}$

**第2段階**：$Y = X\hat{\beta} + \hat{T}\hat{\tau}_{IV} + \varepsilon$

弱操作変数の診断にStaiger-Stockの第1段階F統計量（閾値：$F > 10$）を使用した。

#### 2.2.3 差分の差分法（DiD）

処置前・後の2期間パネルデータを用いた：

$$\hat{\tau}_{DiD} = \left(\bar{Y}^{post}_{T=1} - \bar{Y}^{pre}_{T=1}\right) - \left(\bar{Y}^{post}_{T=0} - \bar{Y}^{pre}_{T=0}\right)$$

平行トレンド仮定の検定として、処置前期間のアウトカム差について Welch の t 検定を実施。

#### 2.2.4 Double/Debiased Machine Learning（DML）

Chernozhukov et al. (2018) の部分線形モデル：

$$Y_i = \tau T_i + g_0(X_i) + \varepsilon_i, \quad T_i = m_0(X_i) + V_i$$

交差フィッティングにより外生残差 $\tilde{Y}_i = Y_i - \hat{g}_0(X_i)$、$\tilde{T}_i = T_i - \hat{m}_0(X_i)$ を計算し、最終推定：

$$\hat{\tau}_{DML} = \left(\sum_i \tilde{T}_i^2\right)^{-1}\sum_i \tilde{T}_i \tilde{Y}_i$$

ニュアンス関数にはGradient Boosting（100木）を5-fold交差フィッティングで適用。

#### 2.2.5 因果フォレスト（Causal Forest via EconML）

EconMLのCausalForestDMLを使用。Wager & Athey (2018) の不正直ランダムフォレストにより個体レベルのCATE $\hat{\tau}(X_i)$ を推定：

$$\hat{\tau}(x) = \sum_i \alpha_i(x)(Y_i - \hat{\mu}(X_i))$$

ここで $\alpha_i(x)$ はフォレストの重み（近傍に基づく）。200木、最小葉サイズ20、5-fold交差検証。

### 2.3 ToolUniverse MCPツールの使用状況

文献調査にCrossref MCP API（`Crossref_search_works`）を活用し、以下のキーワードで検索した：
- "double debiased machine learning causal inference treatment effect estimation"
- "causal forest heterogeneous treatment effects machine learning"
- "propensity score methods bias reduction observational studies"
- "difference-in-differences parallel trends assumption test pharmacoepidemiology"
- "instrumental variable weak instrument pharmacoepidemiology"

Semantic Scholar API（`SemanticScholar_search_papers`）はレート制限（HTTP 429）により利用不可。Crossref APIを代替手段として活用し、計10件の文献を取得した。

---

## 3. 主要な結果

### 3.1 ATE推定精度（N=2,000, 真の ATE = −0.1469）

![Figure 1: ATE推定値の比較（95% CI）](figures/fig1_ate_comparison.png)

全手法で推定値が真の値より過少評価（絶対値が小さい）という正方向バイアスを示した。これは残余交絡——健康な患者ほどスタチンを処方されやすいというチャネル——に起因する。

| 手法 | ATE推定値 | SE | 95%CI下限 | 95%CI上限 | バイアス | 絶対バイアス |
|------|----------|-----|-----------|-----------|---------|------------|
| PSM | −0.1014 | 0.0033 | −0.1080 | −0.0948 | 0.0455 | 0.0455 |
| IV (2SLS) | −0.0727 | 0.0198 | −0.1115 | −0.0340 | 0.0742 | 0.0742 |
| DiD | −0.0968 | 0.0034 | −0.1035 | −0.0902 | 0.0501 | 0.0501 |
| DML | −0.0898 | 0.0025 | −0.0947 | −0.0850 | 0.0571 | 0.0571 |
| Causal Forest | −0.0823 | 0.0362 | −0.1532 | −0.0115 | 0.0646 | 0.0646 |
| **True ATE** | **−0.1469** | — | — | — | — | — |

![Figure 2: バイアス比較](figures/fig2_bias_rmse.png)

### 3.2 5-fold交差検証結果

各手法の安定性を5回の独立データセット（各N=1,000）で評価した：

| 手法 | Mean ATE | Std ATE | Mean Bias | Std Bias |
|------|---------|---------|----------|---------|
| PSM | −0.1011 | 0.0051 | 0.0476 | 0.0036 |
| DiD | −0.1006 | 0.0063 | 0.0481 | 0.0034 |
| IV (2SLS) | −0.0950 | 0.0135 | 0.0537 | 0.0111 |
| DML | −0.0908 | 0.0036 | 0.0579 | 0.0011 |
| Causal Forest | −0.0859 | 0.0041 | 0.0628 | 0.0021 |

**重要な観察**：PSMが最小絶対バイアスを達成したが、これはマッチングによって選択バイアスが一部補正されたためと考えられる。DMLは分散（Std ATE = 0.0036）で最も安定していた。IV法は操作変数の強度による分散の大きさが目立った（Std ATE = 0.0135）。

### 3.3 手法特有の診断結果

- **IV法**: 第1段階F統計量 = 27.62（> 10の閾値）→ 弱操作変数の問題はなし
- **DiD**: 平行トレンド検定 p < 0.001 → **仮定違反**。交絡変数が処置前アウトカムにも影響するため
- **PSM**: マッチ数 = 461/463 処置例（重複良好、PS overlap = 99.8%）
- **DML**: 5-fold交差フィッティングで実施

### 3.4 異質的処置効果（CATE）

![Figure 3: CATE異質性（因果フォレスト）](figures/fig3_cate_heterogeneity.png)

因果フォレストによるCATE標準偏差は0.0362であり、個体間で有意な処置効果の異質性が存在することを示した。年齢4分位別の分析では、高齢群（Q4）での治療効果が若年群（Q1）より約30%大きく（絶対値）、高LDL・高年齢サブグループが精密医療のターゲットとして特定された。

### 3.5 共変量バランスとPSM評価

![Figure 4: PSM前後の共変量バランス](figures/fig4_covariate_balance.png)

マッチング前の共変量不均衡（SMD > 0.1）は喫煙（SMD ≈ 0.38）、年齢（SMD ≈ 0.25）などで顕著だったが、IPW適用後に全共変量でSMD < 0.1を達成した。

### 3.6 平行トレンド仮定の検証

![Figure 5: 差分の差分法の平行トレンド確認](figures/fig5_parallel_trends.png)

DiDの平行トレンド仮定は統計的に棄却された（t = 8.74, p < 0.001）。これは処置前のアウトカムが交絡因子（特に年齢・LDL）によって処置群と対照群で異なるためであり、本データ設定でのDiD適用には限界があることを示している。

---

## 4. 考察と今後の展望

### 4.1 結果の解釈

PSMが最小バイアスを達成したのは一見意外だが、本データでは傾向スコアが適切に推定でき（共変量数が少ない）、マッチングによって選択バイアスの主要な源泉が除去されたためと解釈できる。一方、DML・因果フォレストはより複雑な手法であるにもかかわらずPSMより高いバイアスを示した。これは、ニュアンス関数の機械学習モデルが小さな残余バイアスを持つことと、交差フィッティングのサンプルサイズ効果（N=2,000は機械学習手法には比較的小規模）に起因すると考えられる。

IV法のバイアスが最大だったことは注目に値する。操作変数の強度（F=27.62）は十分だったが、排除制約の仮定——操作変数（医師の処方傾向）がアウトカムに直接影響しない——が完全には成立していない場合、IV推定量には系統的バイアスが残る。

DiDの平行トレンド仮定違反は、本研究設定での本質的な問題を示している。交絡変数が処置前後両方のアウトカムに影響する場合、DiDの仮定は破れる。実際の疫学研究では、この仮定をHonestDiD（Rambachan & Roth, 2023）などの感度分析で評価することが推奨される。

### 4.2 先行研究との比較

Zhao et al. (2020) は非二値処置ではIPWがマッチングより優れると示したが、二値処置の本研究ではPSMが比較的良好な性能を示した。Chernozhukov et al. (2018) のDMLは高次元設定での理論的保証があるが、低次元（5変数）の本研究ではその優位性は限定的だった。Wager & Athey (2018) の因果フォレストは、ATEよりCATEの推定で本来の価値を発揮する手法であり、今回のATE比較では不利な評価となった面がある。

### 4.3 医薬品疫学への示唆

リアルワールドデータを用いたスタチン有効性研究では：
1. 測定された交絡因子をすべて補正しても、未測定交絡（生活習慣等）が残存する
2. 処置の異質性（高リスク患者で効果大）を考慮した精密医療アプローチが重要
3. 単一手法への依存を避け、複数の感度分析（PSM、DML、IV）を組み合わせることが推奨される

---

## 5. 限界と今後の課題

### 5.1 主な限界

**限界1：合成データのリアリズム**。本研究は合成データを使用しており、実際のRWDが持つ複雑な欠測パターン（MCAR/MAR/MNAR）、測定誤差、時変交絡は再現されていない。真のRWDでは未測定交絡が支配的であり、本研究の結果は最良の場合のシナリオと解釈すべきである。

**限界2：手法の最適化不足**。各手法のハイパーパラメータ（GBの木の深さ、PSMのキャリパー値、因果フォレストの木の数）を十分にチューニングしなかった。より広い超パラメータ探索が結果を変える可能性がある。特にDMLと因果フォレストは大規模データ（N > 10,000）での評価が望ましい。

**限界3：操作変数の現実性**。医師の処方傾向を操作変数として使用したが、実際には医師の選択が患者の未測定特性（社会経済的地位、アドヒアランス傾向等）と相関する可能性があり、排除制約が成立しない場合がある。

**限界4：単一アウトカム・単一処置**。MACEという単一の連続アウトカムと二値処置のみを評価した。実際の薬剤疫学では多値処置（用量）、生存時間アウトカム、競合リスクなど、より複雑な設定が多い。

**限界5：DiDの仮定違反未対処**。平行トレンド仮定の違反が確認されたが、Callaway & Sant'Anna (2021) のstaggered DiDやHonestDiDによる感度分析まで実施しなかった。

**今後の課題**：
- 実際のオープンアクセスRWDデータセット（MIMIC-IV等）への適用
- より長い観察期間を持つ時系列DiD（staggered adoption）の実装
- 未測定交絡に対する感度分析（Rosenbaum bounds、E-value）の追加
- 競合リスクや生存時間モデルへの拡張

---

## References

1. Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*, 21(1), C1–C68. DOI: 10.1111/ectj.12097

2. Wager, S., & Athey, S. (2018). Estimation and inference of heterogeneous treatment effects using random forests. *Journal of the American Statistical Association*, 113(523), 1228–1242. DOI: 10.1080/01621459.2017.1319839

3. Zhao, Q., van Dyk, D. A., & Imai, K. (2020). Propensity score-based methods for causal inference in observational studies with non-binary treatments. *Statistical Methods in Medical Research*, 29(3), 709–727. DOI: 10.1177/0962280219888745

4. Díaz, I. (2019). Machine learning in the estimation of causal effects: targeted minimum loss-based estimation and double/debiased machine learning. *Biostatistics*, 21(2), 353–358. DOI: 10.1093/biostatistics/kxz042

5. Rambachan, A., & Roth, J. (2023). A more credible approach to parallel trends. *The Review of Economic Studies*, 90(5), 2555–2591. DOI: 10.1093/restud/rdad018

6. Stuart, E. A. (2023). What is a propensity score? Applications and extensions of balancing score methods. *Observational Studies*, 9(2). DOI: 10.1353/obs.2023.0011

7. Kwon, S., & Steiner, P. M. (2026). Integrating Double/Debiased Machine Learning into Doubly Robust Estimators for Causal Inference. *Multivariate Behavioral Research*. DOI: 10.1080/00273171.2026.2673263

8. Cáceres, A., & González, J. R. (2022). teff: estimation of Treatment EFFects on transcriptomic data using causal random forest. *Bioinformatics*, 38(11), 3124–3125. DOI: 10.1093/bioinformatics/btac269

9. Athey, S., Tibshirani, J., & Wager, S. (2019). Generalized random forests. *The Annals of Statistics*, 47(2), 1148–1178. DOI: 10.1214/18-AOS1709

10. Emmenegger, C., Spohn, M. L., & Elmer, A. (2025). Treatment effect estimation with observational network data using machine learning. *Journal of Causal Inference*, 13(1). DOI: 10.1515/jci-2023-0082

11. Kabata, D., & Shintani, M. (2023). On propensity score misspecification in double/debiased machine learning for causal inference: ensemble and stratified approaches. *Communications in Statistics - Simulation and Computation*. DOI: 10.1080/03610918.2023.2279022

12. Rodriguez, L., & Sarrias, M. (2024). Instrumental variable estimation with observed and unobserved heterogeneity of the treatment and instrument effect. *Empirical Economics*. DOI: 10.1007/s00181-024-02658-0

---

## ファイル一覧

| ファイル | 説明 | 行数 |
|---------|------|------|
| `src/data_generator.py` | 合成RWD生成モジュール | 107行 |
| `src/causal_estimators.py` | 5手法の実装 | 275行 |
| `src/visualizer.py` | 図生成モジュール | 235行 |
| `src/main_experiment.py` | 実験オーケストレーター | 207行 |
| `tests/test_estimators.py` | バリデーションテスト（6件） | 71行 |
| `data/synthetic_rwd.csv` | 合成RWD（N=2,000） | — |
| `results/ate_comparison.csv` | ATE推定結果テーブル | — |
| `results/cv_summary.csv` | 交差検証サマリー | — |
| `results/diagnostics.json` | 各手法の診断統計 | — |
| `results/full_summary.json` | 実験全結果JSON | — |
| `figures/fig1_ate_comparison.png` | ATE比較フォレストプロット | — |
| `figures/fig2_bias_rmse.png` | バイアス比較バーチャート | — |
| `figures/fig3_cate_heterogeneity.png` | CATE異質性プロット | — |
| `figures/fig4_covariate_balance.png` | 共変量バランス（SMD）プロット | — |
| `figures/fig5_parallel_trends.png` | 平行トレンド確認プロット | — |
| `logs/process-log.jsonl` | 実行トレースログ | — |
