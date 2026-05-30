# 実験レポート：食事成分と腸内細菌叢の相互作用予測システムバイオロジーフレームワーク

---

## 1. 実験目的と背景

### 1.1 研究目的

本実験では、食事成分から腸内細菌叢の動態と代謝産物（短鎖脂肪酸：SCFA）生成を予測するマルチスケール計算フレームワークを設計・実装した。具体的な目標は以下の通り：

1. SHIME（模擬腸管生態系）ベースの消化吸収動態モデルから、食事種別ごとの大腸基質供給量を推定
2. 一般化Lotka-Volterra（gLV）方程式による5種腸内細菌群集の長期動態シミュレーション
3. 化学量論的収率モデルによるSCFA（酢酸・プロピオン酸・酪酸）フラックス予測
4. 3種食事パターン（西洋食・地中海食・高繊維食）における菌叢組成・SCFA産生の比較
5. プロバイオティクス・プレバイオティクス介入後の菌叢応答シミュレーション
6. 発酵食品摂取による抗生物質後菌叢回復のケーススタディ

### 1.2 背景

腸内微生物叢は宿主の代謝、免疫、神経系機能において重要な役割を担い、SCFAは特に腸管上皮バリア維持・免疫調節・全身エネルギー代謝に不可欠なシグナル分子として機能する。食事組成は腸内細菌叢の組成と機能を決定する最も強力な因子の一つであるが、マクロ栄養素の摂取量から菌叢動態・代謝産物産生を定量的に予測するメカニスティックなフレームワークは未だ確立されていない。本研究ではSHIME消化モデル、gLV生態学モデル、SCFA収率モデルの3層統合フレームワークを実装し、その予測性能と限界を評価した。

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 モジュール1：SHIME消化吸収モデル

**構造：** 3コンパートメント（胃・小腸・大腸）の逐次移行モデル

**数理式：** 各コンパートメント $c$、基質 $s$ について：
$$\frac{dS_{c,s}}{dt} = -\left(k_{c,s} + \frac{1}{\tau_c}\right) S_{c,s} + \frac{(1-a_{c-1,s})}{\tau_{c-1}} S_{c-1,s}$$

**パラメータ：**
- 滞留時間 $\tau$: 胃2.0h、小腸3.5h、大腸36.0h
- 消化速度定数 $k$: 食品成分・コンパートメント別に設定
- 吸収率 $a$: 小腸での炭水化物吸収率70%、タンパク65%、脂質75%

**出力：** 1日3食摂取を仮定した大腸基質到達量 [g/h]

### 2.2 モジュール2：gLV細菌群集モデル

**モデル種：** *Bacteroides*, *Lactobacillus*, *Bifidobacterium*, *Faecalibacterium*, *Ruminococcus*

**数理式：**
$$\frac{dN_i}{dt} = N_i \left[ r_i + \sum_{j} A_{ij} N_j + \nu_i^{\text{carb}} \bar{S}_{\text{carb}} + \nu_i^{\text{prot}} \bar{S}_{\text{prot}} \right]$$

**パラメータ設定根拠：** Clark et al. (2021)、Liu et al. (2022) の発表値範囲から設定。繊維分解菌（*Bifidobacterium*, *Faecalibacterium*, *Ruminococcus*）の栄養素感受性を高く設定。

**繊維による自己抑制調節：**
$$A_{ii}^{\text{fibre}} = A_{ii} \cdot \max(0.3,\ 1 - 0.05 \cdot \bar{S}_{\text{carb}})$$

### 2.3 モジュール3：SCFAフラックスモデル

**化学量論的収率行列：** 5種 × 3種SCFA

| 菌種 | 酢酸 | プロピオン酸 | 酪酸 |
|------|------|------------|------|
| *Bacteroides* | 2.50 | 1.20 | 0.30 |
| *Lactobacillus* | 1.80 | 0.40 | 0.60 |
| *Bifidobacterium* | 2.00 | 0.30 | 0.80 |
| *Faecalibacterium* | 0.80 | 0.20 | **3.50** |
| *Ruminococcus* | 1.20 | 0.80 | **2.80** |

**フラックス計算：**
$$F_k(t) = \sum_{i} Y_{ik} \cdot N_i(t) \cdot b_{\text{fibre}}$$

### 2.4 交差検証

- **gLV軌跡予測：** 5分割時系列交差検証（訓練80% → 後半20%予測）
- **SCFA収率モデル：** 15%ガウスノイズを添加した合成データでの5分割交差検証

---

## 3. 主要な結果と数値

### 3.1 SHIME消化モデル：食事種別大腸基質供給量

![Fig. 1 — SHIME消化モデル](figures/fig1_shime_digestion.png)

**大腸への日平均炭水化物到達量：**
- 西洋食: 2.85 g/h
- 地中海食: 3.56 g/h  
- 高繊維食: 4.63 g/h（西洋食比 +63%）

小腸での吸収（炭水化物70%）後に大腸へ到達する未消化繊維量が食事パターンによって大きく異なり、菌叢の基質利用可能量に直接影響する。

---

### 3.2 長期菌叢動態シミュレーション（90日間）

![Fig. 2 — gLV長期菌叢動態](figures/fig2_glv_longterm.png)

すべての食事条件において約10日以内に安定した定常状態に収束した。

**定常状態バイオマス（g/L）：**

| 菌種 | 西洋食 | 地中海食 | 高繊維食 |
|------|--------|----------|---------|
| *Bacteroides* | 1.776 | 1.963 | 2.259 |
| *Lactobacillus* | 1.994 | 2.246 | 2.636 |
| *Bifidobacterium* | 2.571 | 2.972 | 3.597 |
| *Faecalibacterium* | 2.496 | 2.868 | 3.456 |
| *Ruminococcus* | 2.311 | 2.612 | 3.094 |

繊維分解菌3種（*Bifidobacterium*, *Faecalibacterium*, *Ruminococcus*）は高繊維食で最大39%バイオマスが増加。これはLiu et al. (2022)が報告した実験結果と方向性が一致する。

---

### 3.3 SCFAフラックス予測

![Fig. 3 — SCFAフラックス予測](figures/fig3_scfa_flux.png)

**定常状態SCFAフラックス（mmol/h）：**

| 食事 | 酢酸 | プロピオン酸 | 酪酸 |
|------|------|------------|------|
| 西洋食 | 17.94 | 6.05 | 18.99 |
| 地中海食 | 20.32 | 6.81 | 30.33 |
| 高繊維食 | 24.06 | 8.01 | 51.80 |

**主要知見：** 酪酸フラックスが食事依存性を最も強く示し、高繊維食では西洋食比2.73倍に増加（18.99 → 51.80 mmol/h）。酢酸（+34%）・プロピオン酸（+32%）は比較的変動が小さい。

---

### 3.4 定常状態比較ヒートマップ

![Fig. 7 — 定常状態比較](figures/fig7_diet_comparison.png)

バイオマスとSCFAフラックスの全食事・全種比較を可視化。*Bifidobacterium*と*Faecalibacterium*が繊維量に対して最も感受性が高く、SCFA産生への寄与も大きい。

---

### 3.5 プロバイオティクス・プレバイオティクス摂動

![Fig. 4 — 摂動実験](figures/fig4_perturbation.png)

**観察された回復動態：**
- プロバイオティクス（*Lactobacillus* 1.5 g/L追加）：対象種は即時増加後3〜5日以内に定常状態へ収束
- プロバイオティクス（*Bifidobacterium* 1.5 g/L追加）：同様のパターンで4〜6日以内に回復
- プレバイオティクス（FOS/GOS模擬、*Bifidobacterium*・*Faecalibacterium* 80%増加）：最も大きな一時的変動（最大1.7倍）を示し、競合排除によって約5〜7日で平衡に戻る

**重要な知見：** すべての摂動条件において食事依存的な定常状態への収束が観察された。これは腸内微生物叢の生態学的レジリエンスを示しており、プロバイオティクス/プレバイオティクスが一時的な変動を引き起こすが長期的な群集構成は食事が決定する、という先行研究の見解と一致する。

---

### 3.6 発酵食品摂取ケーススタディ（抗生物質後回復）

![Fig. 5 — 発酵食品ケーススタディ](figures/fig5_fermented_food.png)

**実験設計：** 抗生物質摂取をシミュレート（全菌種を重度に減少）後、地中海食のみ（コントロール）vs. 地中海食+発酵食品（介入群）の28日間回復を比較。

**主要結果（28日時点）：**

| 指標 | コントロール | 発酵食品群 |
|------|------------|----------|
| Shannon多様度 H' | 1.598 | 1.579 |
| *Bifidobacterium* バイオマス (g/L) | 2.972 | 3.635 |
| *Faecalibacterium* バイオマス (g/L) | 2.868 | 3.543 |
| 初期（抗生物質後）H' | 0.937 | 0.937 |

**解釈：** コントロール群は約3日以内に事前定常状態に完全回復（H'=1.598）。発酵食品群は*Bifidobacterium*（+22.3%）・*Faecalibacterium*（+23.5%）が有意に高いが、特定菌種への偏りからShannon多様度はわずかに低下（ΔH'=−0.019）。これは発酵食品が多様性よりも特定の健康関連菌種の定着に寄与することを示唆する。

---

### 3.7 モデル検証（交差検証）

![Fig. 6 — 交差検証結果](figures/fig6_cv_results.png)

**交差検証性能（5分割）：**

| モデル | 平均R² | ±標準偏差 | 平均RMSE | ±標準偏差 |
|-------|--------|----------|---------|---------|
| gLV軌跡予測 | 0.673 | 0.043 | 0.256 g/L | 0.033 |
| SCFA収率回帰 | 0.915 | 0.005 | — | — |

**⚠️ 自己批判的評価：**
- **gLV予測 R²=0.673** は合成データ上の値であり、実データへの適用性は保証されない。Joseph et al. (2024)が示すように、FBAベースの菌種間相互作用予測は実験データと相関しないことが多い。
- **SCFA R²=0.915** は同一モデル生成データからのパラメータ回収であり、楽観的バイアスが含まれる。実データでは大幅に低下する可能性がある。
- **合成データの前提条件依存性：** すべての定量値は設定パラメータに依存しており、実世界への一般化には個人差・宿主免疫・環境要因等の追加モデリングが必要。

---

## 4. 考察と今後の展望

### 4.1 フレームワークの意義と限界

本フレームワークは食事→消化→菌叢動態→代謝産物産生の連鎖を単一の計算パイプラインで再現することに成功した。食事パターンによる定性的な予測（繊維多→酪酸産生増→健康関連菌種増加）は先行研究と一致し、実験仮説生成ツールとして有用性がある。

**主な限界：**
1. **5種モデルの単純化：** 実際の腸内細菌叢は数百〜数千のOTUで構成される。高次の種間相互作用・キーストーン種効果が欠如。
2. **静的基質近似：** 1日平均値を使用しており、食事タイミング・食後変動・個人差（腸管通過時間等）が無視されている。
3. **gLV線形性仮定：** 実際の相互作用は非線形・文脈依存。特にクロスフィーディング（乳酸→酪酸産生）の動的機構が欠如。
4. **合成データ検証：** 実験時系列データでの検証が行われておらず、予測値の臨床的意味は不明。
5. **宿主側効果の不在：** 免疫調節、腸管蠕動、胆汁酸分泌等の宿主因子が含まれていない。

### 4.2 実世界への一般化可能性

⚠️ **重要な留意事項：** 本フレームワークの定量的予測値（例：酪酸フラックス50mmol/h）を実際のヒト生理値と直接比較することは適切でない。合成シミュレーションにおける数値はパラメータ選択に強く依存し、実際の個体差（遺伝、年齢、健康状態、抗生物質使用歴等）を反映していない。

本研究の主な貢献は絶対値の予測ではなく、**食事パターンの変化が菌叢動態・SCFA産生に与える相対的影響の方向性と大まかなスケール**を示すことにある。

### 4.3 MICOMおよびgapseqとの統合展望

MICOM（Diener et al. 2020）やgapseq（Zimmermann et al. 2021）などのゲノムスケール代謝モデリングツールとの統合が最も重要な次のステップである。具体的には：

- **AGORA2データベース（7,302菌株）** を利用したSCFAフラックスのFBAベース予測への置き換え
- **SHIME出力をMICOMの培地組成入力** として使用することで、より精密な基質-フラックス関係のモデリングが可能
- **個人ゲノムデータ（16S/メタゲノム）** との統合による個別化予測
- **gapseq** による腸内細菌ゲノムからの代謝ネットワーク自動再構築と統合

### 4.4 今後の実験的検証

1. Liu et al. (2022)のマウスイヌリン/レジスタントスターチ時系列データを使用したgLVパラメータフィッティング
2. 臨床介入試験データ（例：DIETFITS試験）を使用したSCFA予測の検証
3. フローレアクター実験（SHIME実機）との比較によるモデル較正

---

## 5. 生成したファイル一覧

| ファイル名 | 種別 | 内容 |
|-----------|------|------|
| `simulate_gut.py` | Pythonスクリプト | メインシミュレーションコード |
| `figures/fig1_shime_digestion.png` | 図 | SHIME消化モデル：3食事条件の大腸基質時系列 |
| `figures/fig2_glv_longterm.png` | 図 | gLV長期動態（90日）：3食事条件の菌叢軌跡 |
| `figures/fig3_scfa_flux.png` | 図 | SCFAフラックス時系列（3食事条件） |
| `figures/fig4_perturbation.png` | 図 | プロバイオティクス/プレバイオティクス摂動・回復 |
| `figures/fig5_fermented_food.png` | 図 | 抗生物質後回復ケーススタディ（対照vs発酵食品） |
| `figures/fig6_cv_results.png` | 図 | 5分割交差検証結果 |
| `figures/fig7_diet_comparison.png` | 図 | 定常状態比較ヒートマップ |
| `paper.md` | 論文 | 英語学術論文形式のまとめ |
| `report.md` | レポート | 本実験レポート |

---

## 6. 参考文献

1. Fan Y, Pedersen O. (2020). Gut microbiota in human metabolic health and disease. *Nature Reviews Microbiology*, 19, 55–71. https://doi.org/10.1038/s41579-020-0433-9

2. Silva YP, Bernardi A, Frozza RL. (2020). The Role of Short-Chain Fatty Acids From Gut Microbiota in Gut-Brain Communication. *Frontiers in Endocrinology*, 11, 25. https://doi.org/10.3389/fendo.2020.00025

3. Liu H, et al. (2022). Ecological dynamics of the gut microbiome in response to dietary fiber. *The ISME Journal*, 16, 2458–2470. https://doi.org/10.1038/s41396-022-01253-4

4. Heinken A, et al. (2023). Genome-scale metabolic reconstruction of 7,302 human microorganisms for personalized medicine. *Nature Biotechnology*, 41, 1320–1331. https://doi.org/10.1038/s41587-022-01628-0

5. Quinn-Bohmann N, et al. (2024). Microbial community-scale metabolic modelling predicts personalized short-chain fatty acid production profiles in the human gut. *Nature Microbiology*, 9, 1700–1712. https://doi.org/10.1038/s41564-024-01728-4

6. Clark RL, et al. (2021). Design of synthetic human gut microbiome assembly and butyrate production. *Nature Communications*, 12, 3254. https://doi.org/10.1038/s41467-021-22938-y

7. Marco ML, et al. (2021). ISAPP consensus statement on fermented foods. *Nature Reviews Gastroenterology & Hepatology*, 18, 196–208. https://doi.org/10.1038/s41575-020-00390-5

8. Joseph C, et al. (2024). Predicting microbial interactions with FBA: an evaluation. *BMC Bioinformatics*, 25, 87. https://doi.org/10.1186/s12859-024-05651-7

9. Jansma J, El Aidy S. (2021). Understanding the host-microbe interactions using metabolic modeling. *Microbiome*, 9, 16. https://doi.org/10.1186/s40168-020-00955-1

10. de Vos WM, et al. (2022). Gut microbiome and health: mechanistic insights. *Gut*, 71, 1020–1032. https://doi.org/10.1136/gutjnl-2021-326789
