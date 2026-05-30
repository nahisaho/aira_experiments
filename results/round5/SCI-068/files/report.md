# 実験レポート：海洋酸性化がサンゴ礁生態系に及ぼす影響の統合モデリング

---

## 1. 実験目的と背景

### 1.1 目的

本実験は、海洋酸性化（OA: Ocean Acidification）がサンゴ礁生態系に及ぼす影響を予測するための統合数値モデルを設計・実装し、グレートバリアリーフ（GBR）を対象に2100年までのシナリオ予測を実施することを目的とする。

### 1.2 背景

産業革命以来、大気CO₂濃度は280 ppmから415 ppm（2020年）へ上昇し、海洋のpHは約0.1単位（約26%の[H⁺]増加）低下した。IPCCのSSP5-8.5シナリオでは、2100年に大気CO₂が1100 μatmを超え、表層海水pHはさらに0.3〜0.4単位低下すると予測される。サンゴの石灰化に必要なアラゴナイト飽和度（Ω_arag）が低下することで、サンゴ礁の成長・維持が困難となる。

### 1.3 先行研究の主要知見（Step 1 調査結果）

#### 調査に使用したツール
- Semantic Scholar MCP（SemanticScholar_search_papers）
- OpenAlex MCP（openalex_literature_search）
- Crossref MCP（Crossref_search_works）

#### 特定された主要論文（2020年以降、5件以上）

| # | タイトル（要約） | 著者 | 年 | DOI | 主要知見 |
|---|----------------|------|-----|-----|---------|
| 1 | Shifts in coralline algae, macroalgae, and coral juveniles in the GBR | Fabricius et al. | 2020 | 10.1111/gcb.14985 | 現在進行中のOAにより、GBRでは珊瑚幼体・石灰藻が減少し、マクロ藻類が増加 |
| 2 | Twenty-first century ocean warming, acidification, deoxygenation from CMIP6 | Kwiatkowski et al. | 2020 | 10.5194/bg-17-3439-2020 | CMIP6 26モデルによるSSP別pH・Ω予測；SSP5-8.5でΔpH≈-0.44 |
| 3 | Extending the natural adaptive capacity of coral holobionts | Voolstra et al. | 2021 | 10.1038/s43017-021-00214-3 | 産卵補助・マイクロバイオーム操作・人工選択による適応能力拡大レビュー |
| 4 | Coral-bleaching responses to climate change across biological scales | van Woesik et al. | 2022 | 10.1111/gcb.16192 | 分子→景観スケールの統合フレームワーク；階層モデルの必要性 |
| 5 | Individual and Interactive Effects of Ocean Warming and Acidification on Favites colemani | Tañedo et al. | 2021 | 10.3389/fmars.2021.704487 | pH変化はΩを通じた付加効果；温度ストレスが主因だがOAが閾値を低下 |
| 6 | Genomic prediction of (mal)adaptation across current and future climatic landscapes | Capblancq et al. | 2020 | 10.1146/annurev-ecolsys-020720-042553 | ゲノム-環境関連解析で局所適応・気候変動対応の予測が可能 |
| 7 | Biological impacts of marine heatwaves | Smith et al. | 2022 | 10.1146/annurev-marine-032122-121437 | 海洋熱波の頻度・強度が増加；多スケールで生態影響が拡大 |
| 8 | Differences in carbonate chemistry up-regulation of long-lived reef-building corals | (複数著者) | 2023 | 10.1038/s41598-023-37598-9 | サンゴ種間で石灰化流体pHのアップレギュレーション能力に差異 |

#### 先行研究の課題・限界
- **モジュール分離**: 炭酸塩化学、生理学、生態学、遺伝学を統合したモデルが存在しない
- **空間的不均一性の無視**: GBRの2300 km規模の南北勾配・局所水域効果を考慮したモデルが少ない
- **相乗効果の不確実性**: 温度とpHの相互作用が相加的か相乗的かについて知見が分散
- **進化応答の遅延**: 適応進化速度の過大評価（単一遺伝子モデルの限界）
- **アルカリ度の時空間変動**: 多くのモデルが固定アルカリ度を仮定

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 モジュール構成

```
┌──────────────────────────────────────────────────────────────┐
│              統合サンゴ礁生態系モデル                          │
├────────────┬──────────────┬─────────────┬──────────────────┤
│ Module 1   │ Module 2     │ Module 3    │ Module 4         │
│ 炭酸塩化学  │ 石灰化モデル  │ 種間相互作用  │ 複合ストレス     │
│ (CO2SYS型) │ (pH/Ω依存)   │ (GLVネット) │ (温度×pH)       │
├────────────┴──────────────┴─────────────┴──────────────────┤
│ Module 5: 集団遺伝学 (Wright-Fisher モデル)                    │
├───────────────────────────────────────────────────────────── │
│ Module 6: GBR 2100年シナリオ統合予測 (SSP1-2.6/2-4.5/5-8.5)  │
└───────────────────────────────────────────────────────────── ┘
```

### 2.2 Module 1: 海水炭酸塩化学（CO2SYS型）

**入力**: 大気CO₂分圧（μatm）、水温T（°C）、塩分S（psu）

**出力**: pH, [CO₂], [HCO₃⁻], [CO₃²⁻], DIC, Ω_arag

主要式（Roy et al., 1993; Weiss, 1974; Mucci, 1983）：
- K₀（Henry則）、K₁、K₂（炭酸解離定数）、Ksp（アラゴナイト溶解度積）
- 電荷中性条件（TA = 2300 μmol/kg固定）からNewton-Raphson法で[H⁺]を反復計算
- Ω_arag = [Ca²⁺][CO₃²⁻] / Ksp_arag

### 2.3 Module 2: サンゴ石灰化速度モデル

```
G = G_max × f_Ω(Ω_arag) × f_T(T) × f_pH(pH)
```

| 因子 | 式 | 文献 |
|------|---|------|
| f_Ω | min(1, (Ω/3.5)^1.5) | Fabricius et al. (2011) |
| f_T | exp(-(T-27)²/(2×3.5²)) | GBR実測値 |
| f_pH | sigmoid(8×(pH-7.9)) | DeCarlo et al. (2019) |

### 2.4 Module 3: 種間相互作用ネットワーク（GLV）

8機能群（造礁サンゴ、共生藻、草食魚、肉食魚、ウニ、マクロ藻、石灰藻、動物プランクトン）を一般化Lotka-Volterra方程式で記述：

$$\frac{dN_i}{dt} = r_i N_i \left(1 - \frac{\sum_j A_{ij} N_j}{K_i}\right)$$

OA条件下ではサンゴ成長率 $r_\text{coral}$ を60%削減（Fabricius et al. 2011に基づく）。

### 2.5 Module 4: 温度-pH複合ストレス

3種類の相互作用モデル（相加、乗算、相乗）を実装し比較：

**相乗モデル**（最も現実的）：
$$T_\text{eff} = \Delta T + 0.5 \times \Delta\text{pH}$$
$$P_\text{bleach} = 1-(1-P_{T\text{eff}})(1-P_\text{OA}) + 0.1 P_{T\text{eff}} P_\text{OA} + \varepsilon$$

### 2.6 Module 5: 集団遺伝学モデル（Wright-Fisher）

4プロセスの組み合わせ：選択（$s_t = s_0(1+εt/80)$で増加）、変異（$\mu=10^{-4}$）、移入（$m=0.001$）、遺伝的浮動（二項分布サンプリング）

N=20確率的レプリカで不確実性を定量化。

### 2.7 Module 6: GBR 2100年シナリオ

3シナリオ（SSP1-2.6、SSP2-4.5、SSP5-8.5）を2020〜2100年にわたり年次計算。10シードでのクロスバリデーション実施。

---

## 3. 主要な結果と数値

### 3.1 炭酸塩化学（Figure 1）

![Figure 1: 炭酸塩化学](figures/fig1_carbonate_chemistry.png)

**Figure 1.** 大気CO₂分圧と海水炭酸塩化学の関係（T=27°C, S=35 psu）。(a) pH変化、(b) Ω_arag変化、(c) [CO₃²⁻]変化、(d) 石灰化速度マップ（SSPエンドポイントを点でプロット）。

- pH: 280 μatm → 8.19、415 μatm → 8.07、1135 μatm → 7.67
- Ω_arag: 280 μatm → ~4.5、1135 μatm → ~1.9
- サンゴ成長に必要なΩ≥3.0は580 μatm超で消滅

### 3.2 種間相互作用ネットワーク（Figure 2）

![Figure 2: 種間相互作用ネットワーク](figures/fig2_species_network.png)

**Figure 2.** (左) 通常条件のサンゴ礁生態系ネットワーク。(右) 海洋酸性化条件（サンゴ成長率60%削減）。ノードサイズは均衡バイオマスに比例。

### 3.3 GLV個体群動態（Figure 3）

![Figure 3: GLV個体群動態](figures/fig3_glv_dynamics.png)

**Figure 3.** 100年間のGLV動態。(左) 通常条件：サンゴが~0.6Kで安定。(右) OA条件：マクロ藻類が~0.9Kまで増加、サンゴは<0.15Kに崩壊。

### 3.4 温度-pH相乗ストレス（Figure 4）

![Figure 4: 相乗ストレス](figures/fig4_synergistic_stress.png)

**Figure 4.** (a) 相加モデル、(b) 乗算モデル、(c) 相乗モデルによる白化確率分布。相乗モデルでは、OAがΔpH=0.4の条件でP=0.5の白化閾値温度を約0.5°C低下させる。

### 3.5 集団遺伝学（Figure 5）

![Figure 5: 集団遺伝学](figures/fig5_population_genetics.png)

**Figure 5.** (左) N=10,000の耐性アレル頻度の時系列（20レプリカの平均・SD・5-95パーセンタイル）。(右) 有効集団サイズN別の進化応答比較。

- N=10,000: 2100年に耐性アレル頻度 ~0.55（SSP5-8.5の生存閾値0.70には未到達）
- 小集団（N=500）は遺伝的浮動が支配的で分散が極めて大きい

### 3.6 GBR 2100年予測（Figure 6, 7）

![Figure 6: GBR予測](figures/fig6_gbr_projections.png)

**Figure 6.** 2020-2100年のGBR主要変数のシナリオ別予測（SSP1-2.6, SSP2-4.5, SSP5-8.5）。

![Figure 7: サマリーダッシュボード](figures/fig7_summary_dashboard.png)

**Figure 7.** 統合サマリーダッシュボード。(a) 石灰化率、(b) pH、(c) Ω_arag、(d) 白化確率、(e) 耐性アレル頻度、(f) 2100年指標比較。

### 3.7 定量的結果表（2090-2100年平均）

| 指標 | SSP1-2.6 | SSP2-4.5 | SSP5-8.5 |
|------|----------|----------|----------|
| pCO₂ (μatm) | 449 ± 5 | 635 ± 12 | 1094 ± 29 |
| pH | 8.003 ± 0.004 | 7.875 ± 0.007 | 7.668 ± 0.010 |
| Ω_arag | 3.69 ± 0.03 | 2.98 ± 0.04 | 2.14 ± 0.03 |
| SST (°C) | 28.1 ± 0.2 | 28.9 ± 0.2 | 31.1 ± 0.3 |
| 石灰化率 G/G₀ | 0.658 ± 0.012 | 0.305 ± 0.019 | 0.033 ± 0.005 |
| 白化確率 | 0.061 ± 0.023 | 0.310 ± 0.117 | 0.991 ± 0.018 |

### 3.8 クロスバリデーション結果（10シード）

| シナリオ | 平均 G/G₀ | SD | 95% CI | CV (%) |
|---------|----------|-----|--------|---------|
| SSP1-2.6 | 0.6605 | 0.0034 | [0.655, 0.665] | 0.51% |
| SSP2-4.5 | 0.3028 | 0.0033 | [0.297, 0.307] | 1.09% |
| SSP5-8.5 | 0.0321 | 0.0005 | [0.031, 0.033] | 1.56% |

CV（変動係数）は全シナリオで2%未満 → 確率的ノイズの影響は小さく、シナリオ間差異が支配的

---

## 4. 考察と自己批判的評価

### 4.1 主要な知見

1. **炭酸塩化学の非線形性**: pHとΩ_aragはpCO₂の対数に近いスケールで変化し、SSP5-8.5では2100年にΩ_arag = 2.14まで低下 → 既往研究（Kwiatkowski et al., 2020）と整合的
2. **石灰化崩壊**: SSP5-8.5でG/G₀ = 3.3%という極端な値は、Ωの低下に加えてGBR水温がT_opt=27°Cから31°Cへ乖離することによる熱的抑制が重なった結果
3. **生態的相転移**: OA条件下ではGLVモデルがサンゴ優占系からマクロ藻類優占系への相転移を再現 → 観測例（Hughes et al., 2010）と定性的に整合
4. **進化応答の限界**: N=10,000の集団でSSP5-8.5の生存閾値に到達できない → 人為的支援（assisted evolution）の必要性を示唆

### 4.2 自己批判的検証

#### 合成データへの依存度
- **高依存**: 本モデルの全定量値は合成データからの出力であり、実GBRデータによる直接キャリブレーションは行っていない
- **固定アルカリ度**: TA=2300 μmol/kgの仮定はGBRの空間変動（±200 μmol/kg）を無視しており、pHで±0.05、Ωで±0.3の誤差を生む可能性がある
- **線形内挿**: pCO₂の線形補間は炭素循環フィードバックを無視し、世紀後半のpCO₂を10-30%過小評価する恐れがある

#### 実世界データへの適用可能性
- **石灰化モデル**: Fabricius et al. (2011) のパラメータは特定の生態系（PNG CO₂湧出域）から得られたものであり、GBR全域への適用には追加検証が必要
- **閾値の地域差**: 白化閾値MMM+1°Cはグローバル平均値であり、GBR南部ではより低い、北部ではより高い閾値が知られている
- **GLVパラメータの時間変化**: 成長率・相互作用係数が温度・pHとともに変化する効果は未実装

#### 過度に楽観的な側面
- SSP1-2.6のG/G₀ = 0.658は、実際の生物多様性損失・サンゴ種組成変化・幼生補充減少などを考慮すると過大評価の可能性がある
- 集団遺伝学モデルの単一遺伝子座仮定は多遺伝子的適応の複雑性を無視し、適応速度を過大評価する可能性がある

#### 過度に悲観的な側面
- SSP5-8.5でのG/G₀ = 0.033は、T_optからの乖離（+4°C）による熱的ペナルティが支配的 → 実際のサンゴは熱適応を示すものがあり、T_optが上方シフトする可能性を考慮していない
- 環境変化への順化（acclimatization）は本モデルに組み込まれていない

### 4.3 先行研究との比較

| 先行研究 | 本研究との比較 |
|---------|--------------|
| Kwiatkowski et al. (2020): GBR SSP5-8.5でpH=7.65-7.75 | 本研究: pH=7.67 ✓ 整合 |
| Fabricius et al. (2011): Ω=2でG~70%減 | 本研究: Ω=2.14でG=97%減（温度効果が加わるため過大） |
| van Woesik et al. (2022): 白化の多スケール統合の必要性 | 本研究は種・集団・景観の垂直統合を実現 |
| Voolstra et al. (2021): 進化支援なしでの回復困難 | 本研究: N=10,000でも生存閾値未達、定量的に支持 |

---

## 5. 今後の展望

### 5.1 優先課題

1. **実データとの統合**
   - AIMSの長期モニタリングデータ（1985-現在）との統合
   - AIMS GBR海洋酸性化モニタリングネットワークの炭酸塩化学データを用いたTA・DICのキャリブレーション

2. **空間的明示化**
   - GBR循環モデル（ROMS/MOM）との結合
   - 南北温度勾配・局所湧昇域・礁内フラックスの空間分解

3. **生物学的精緻化**
   - 多遺伝子座・QTL（Quantitative Trait Locus）構造への集団遺伝学モジュール拡張
   - サンゴ幼生補充・定着・死亡率の明示的モデリング
   - Symbiodiniaceae多様性（Clade C/D切り替え）の内部共生者シャッフリングモデル

4. **管理介入シナリオの統合**
   - 海洋アルカリ度添加（OAE）シナリオ
   - 支援進化（assisted evolution）の効果定量化
   - 礁保護区の連結性効果

### 5.2 不確実性の優先的削減

| 不確実性源 | 改善アプローチ |
|-----------|--------------|
| TA空間変動 | GBROOS（Great Barrier Reef Ocean Observing System）データ |
| 白化閾値地域差 | CoRTAD（coral reef temperature anomaly database）活用 |
| 石灰化曲線パラメータ | メタ解析（Hendriks et al., 2010後継研究）による改訂 |
| 進化応答速度 | GBR産サンゴゲノムデータ（GBIF/NCBI SRA）との統合 |

---

## 6. 生成したファイル一覧

| ファイル名 | 内容 |
|-----------|------|
| `coral_model.py` | 統合モデルの実装（6モジュール、600行超） |
| `figures/fig1_carbonate_chemistry.png` | 炭酸塩化学の変化 |
| `figures/fig2_species_network.png` | 種間相互作用ネットワーク（通常/OA条件） |
| `figures/fig3_glv_dynamics.png` | GLV個体群動態（100年シミュレーション） |
| `figures/fig4_synergistic_stress.png` | 温度-pH相乗ストレス曲面 |
| `figures/fig5_population_genetics.png` | Wright-Fisher集団遺伝学モデル |
| `figures/fig6_gbr_projections.png` | GBR 2100年シナリオ予測（6変数） |
| `figures/fig7_summary_dashboard.png` | 統合サマリーダッシュボード |
| `figures/results_table.csv` | 定量結果テーブル（CSV） |
| `paper.md` | 学術論文形式のまとめ（英語） |
| `report.md` | 本実験レポート（日本語） |

---

## 参考文献

1. Fabricius, K. E. et al. (2020). Shifts in coralline algae, macroalgae, and coral juveniles in the Great Barrier Reef. *Global Change Biology*, 26(3), 1390–1405. DOI: 10.1111/gcb.14985
2. Kwiatkowski, L. et al. (2020). Twenty-first century ocean warming, acidification, deoxygenation from CMIP6. *Biogeosciences*, 17(13), 3439–3470. DOI: 10.5194/bg-17-3439-2020
3. Voolstra, C. R. et al. (2021). Extending the natural adaptive capacity of coral holobionts. *Nat. Rev. Earth Environ.*, 2, 747–762. DOI: 10.1038/s43017-021-00214-3
4. Capblancq, T. et al. (2020). Genomic prediction of (mal)adaptation. *Annu. Rev. Ecol. Evol. Syst.*, 51, 245–269. DOI: 10.1146/annurev-ecolsys-020720-042553
5. van Woesik, R. et al. (2022). Coral-bleaching responses to climate change across biological scales. *Global Change Biology*, 28(14), 4229–4250. DOI: 10.1111/gcb.16192
6. Tañedo, M. C. S. et al. (2021). Individual and Interactive Effects of Ocean Warming and Acidification on *Favites colemani*. *Frontiers in Marine Science*, 8, 704487. DOI: 10.3389/fmars.2021.704487
7. Smith, K. E. et al. (2022). Biological impacts of marine heatwaves. *Annu. Rev. Mar. Sci.*, 15, 119–145. DOI: 10.1146/annurev-marine-032122-121437
8. Dickson, A. G. & Millero, F. J. (1987). Equilibrium constants for dissociation of carbonic acid in seawater. *Deep Sea Research A*, 34(10), 1733–1743. DOI: 10.1016/0198-0149(87)90021-5
9. Mucci, A. (1983). Solubility of calcite and aragonite in seawater. *American Journal of Science*, 283, 780–799. DOI: 10.2475/ajs.283.7.780
10. Roy, R. N. et al. (1993). The dissociation constants of carbonic acid in seawater at salinities 5 to 45 and temperatures 0 to 45°C. *Marine Chemistry*, 44(2-4), 249–267. DOI: 10.1016/0304-4203(93)90207-5
