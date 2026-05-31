# 実験レポート: 海洋酸性化がサンゴ礁生態系に及ぼす影響の統合予測モデル

**プロジェクト**: CO2ReefSys — サンゴ礁生態系統合シミュレーションフレームワーク  
**作成日**: 2026-05-31  
**作成者**: GitHub Copilot (Claude Sonnet 4.6)  
**実行環境**: Python 3.11.2, NumPy 2.3.5, SciPy 1.17.1

---

## 1. 実験目的と背景

### 目的

大気中CO₂濃度の上昇による海洋酸性化がグレートバリアリーフ（GBR）のサンゴ礁生態系に及ぼす影響を、2100年までのIPCC AR6 SSPシナリオに基づいて統合的に予測するモデル（CO2ReefSys）を設計・実装する。

### 研究背景

海洋酸性化は産業革命前からのpHが約0.11単位下降（~8.18 → ~8.07）した現象で、現在の進行速度は過去5500万年間で前例がない（IPCC AR6, 2021）。特にサンゴ礁にとって危機的なのはアラゴナイト飽和度（Ω_arag）の低下であり、Ω_arag < 3では健全なサンゴ礁形成が困難になる。加えて、海水温上昇（白化）と酸性化の複合ストレスは相乗効果（Synergy）を示すことが複数の研究で確認されている（Allison et al., 2021）。

---

## 2. 使用した手法・アルゴリズムの概要

CO2ReefSysは6つのモジュールを統合した計算フレームワークである：

| モジュール | 手法 | 検証状況 |
|-----------|------|---------|
| **1. 炭酸塩化学** | 熱力学的平衡定数（K0: Weiss 1974, K1/K2: Lueker 2000, Ksp: Mucci 1983） | 既知値と一致（pH=8.025, Ω=3.37 @420μatm） |
| **2. 石灰化速度モデル** | 非対称ガウス温度応答 × Ω依存べき乗則 | 種別パラメータを文献から設定 |
| **3. 生態系ネットワーク** | 一般化Lotka-Volterra ODE（7種） | 基準状態での安定性を確認 |
| **4. 複合ストレスモデル** | 相乗効果指数（Synergy Index） | Allison et al. (2021)と定性的一致 |
| **5. 集団遺伝学** | Wright-Fisher模型（有限集団 + 方向性選択） | 50レプリケートで収束確認 |
| **6. シナリオ予測** | SSP1-2.6/SSP2-4.5/SSP5-8.5の炭酸塩・生態応答計算 | IPCC AR6準拠の強制力設定 |

### 先行研究調査（ToolUniverse Semantic Scholar使用）

Semantic Scholar MCPツールにより以下の論文を特定：

1. **Spreter et al. (2022)** — Arabian Sea上昇流によるOA影響、栄養塩の部分的緩和効果
2. **Noonan et al. (2025)** — 段階的OA下での群集変化（Acropora > Porites > Pavona 感受性順）
3. **Allison et al. (2021)** — 温度とOAの非加算的（相乗的）相互作用
4. **Fuller et al. (2020)** — A. milleporaゲノムワイド関連解析、白化耐性の多遺伝子スコア
5. **González-Espinosa & Donner (2021)** — Random Forestモデル（精度0.834）、雲被覆が白化を緩和
6. **Boonnam et al. (2022)** — SVM（88.85%精度）、pH・SSTが白化の主要予測因子
7. **Jagadeesh & Pradhan (2025)** — XGBoost（R²=0.85）、機械学習による白化予測

### NatureLM / GALACTICA MCPの試行結果

- **NatureLM MCP**: `tooluniverse-grep_tools`でパターン"NatureLM"を検索 → **0件ヒット（利用不可）**
- **GALACTICA MCP**: `tooluniverse-grep_tools`でパターン"GALACTICA"を検索 → **0件ヒット（利用不可）**
- **代替手段**: 文献値による直接パラメータ検証、Semantic Scholarによる追加文献調査

---

## 3. 主要な結果と数値

### 3.1 炭酸塩化学計算結果

現在条件（T=26°C, S=35, pCO₂=420 μatm, TA=2300 μmol/kg）：

| パラメータ | 値 |
|-----------|-----|
| pH (total scale) | **8.025** |
| Ω_aragonite | **3.37** |
| DIC | 2004 μmol kg⁻¹ |
| [CO₃²⁻] | 211 μmol kg⁻¹ |

SSP5-8.5条件（850 μatm）: pH = 7.765、Ω_arag = 2.06（ΔpH = −0.260 units、ΔΩ = −1.31）

![Figure 1: 炭酸塩化学と石灰化](figures/fig1_carbonate_chemistry.png)

**図1.** 大気pCO₂に対する海水炭酸塩化学とサンゴ石灰化速度の変化。(A) pH低下曲線、(B) アラゴナイト飽和度、(C) 炭酸イオン濃度、(D) 種別石灰化速度対Ω曲線、(E) 石灰化速度対pCO₂（種比較）、(F) Acroporaの温度応答（4シナリオ）。

### 3.2 種別石灰化速度

産業革命前から高排出シナリオ（SSP5-8.5）への変化：

| 種 | 前産業時代（280 μatm） | 現在（420 μatm） | SSP5-8.5（850 μatm） | 前産業比変化 |
|----|----------------------|-----------------|---------------------|------------|
| Acropora millepora | 486.7 | 332.9 | 129.4 | **−73.5%** |
| Porites lobata | 230.2 | 169.5 | 79.2 | **−65.6%** |

（単位: mmol CaCO₃ m⁻² d⁻¹、T=26°C）

### 3.3 生態系動態モデル

7種一般化Lotka-Volterra ODEの結果（2020〜2100年）：

| シナリオ | 最終総サンゴ被覆率 | 基準線からの変化 |
|---------|------------------|----------------|
| 基準線（無気候変動） | **35.1%** | — |
| OAのみ（SSP5-8.5） | **3.8%** | **−89.1%** |
| 複合（OA + +3°C昇温） | **1.3%** | **−96.3%** |

初期被覆率45.0%から基準線では生態系は競争平衡に向かって安定（35.1%）。
OA+昇温複合では藻類被覆率が増加（OAのみ: 9.1%）、ウニ個体数が37.0%（基準線19.0%）に急増。

![Figure 2: 生態系動態](figures/fig2_ecosystem_dynamics.png)

**図2.** 7種サンゴ礁生態系の動態シミュレーション。(A) 基準線、(B) OAのみ（SSP5-8.5）、(C) OA+昇温複合、(D) 2100年時点の種別サンゴ被覆率、(E) 総サンゴ被覆率の推移、(F) 種間相互作用行列。

### 3.4 温度-pH複合ストレスの相乗効果

T=28°C、Ω=2.50（SSP2-4.5シナリオ2060年頃の予測値）における石灰化：

- Acropora: コントロール比 **37.6%**（相乗指数 +0.46：超加算的ストレス）
- Porites: **67.3%**（相対的に耐性）
- Pavona: **56.5%**
- Montipora: **47.9%**

T=29°C（SSP5-8.5 2080年頃）ではAcroporaが **13.8%** に低下。

![Figure 3: 複合ストレス](figures/fig3_combined_stress.png)

**図3.** 4種のサンゴに対する温度-pH複合ストレス応答曲面。各パネルは石灰化速度の等値線図（% of max）。

### 3.5 GBR 2100年予測シナリオ

| シナリオ | pCO₂ 2100 (μatm) | 水温 (°C) | pH 2100 | Ω_arag 2100 | 石灰化率（2020比） | 礁健全度指数 |
|---------|------------------|---------|---------|------------|-----------------|------------|
| **SSP1-2.6** | 403 | 27.5 | 8.037 | 3.62 | **90.4%** | **88.1%** |
| **SSP2-4.5** | 728 | 28.5 | 7.823 | 2.53 | **32.2%** | **61.4%** |
| **SSP5-8.5** | 1262 | 29.5 | 7.614 | 1.73 | **7.5%** | **33.0%** |

礁健全度指数（RHI）= pH寄与（30%）× Ω_arag寄与（30%）× 温度寄与（40%）の加重複合指数。

![Figure 5: GBR 2100年予測](figures/fig5_gbr_projections.png)

**図5.** 3つのSSPシナリオ下での2020〜2100年GBR予測。(A) pCO₂軌跡、(B) 水温軌跡（白化閾値表示）、(C) pH予測、(D) アラゴナイト飽和度、(E) Acropora石灰化速度、(F) 複合礁健全度指数。

### 3.6 集団遺伝学・局所適応

Wright-Fisher模型（150世代、50レプリケート）の結果：

| 有効集団サイズ (Ne) | 熱耐性対立遺伝子最終頻度 | OA耐性対立遺伝子最終頻度 | 対立遺伝子消失率（熱耐性） |
|--------------------|------------------------|------------------------|--------------------------|
| 100 | 0.203 ± 0.358 | 0.253 ± 0.341 | **72%** |
| 1000 | 0.445 ± 0.184 | 0.307 ± 0.192 | 4% |
| 5000 | 0.461 ± 0.083 | 0.321 ± 0.080 | **0%** |

小さな孤立礁（Ne=100）では熱耐性対立遺伝子が72%の確率で消失（遺伝的浮動による）。大規模礁（Ne=5000）では消失確率0%、熱耐性頻度は平均0.461（±0.083）。

![Figure 4: 集団遺伝学](figures/fig4_population_genetics.png)

**図4.** サンゴ局所適応の集団遺伝学的シミュレーション。(A) 熱耐性対立遺伝子動態、(B) OA耐性対立遺伝子動態、(C) 固定確率対Ne、(D) 適応確率ヒートマップ、(E) 最終頻度分布、(F) 環境変化速度対適応確率。

### 3.7 機械学習予測モデル

800サンプルの合成データ（10%ガウスノイズ付加）による5分割交差検証：

| モデル | R²（テスト）± SD | R²（訓練）± SD | RMSE（テスト） |
|-------|----------------|--------------|-------------|
| Random Forest | **0.902 ± 0.009** | 0.979 ± 0.001 | 3.25 ± 0.36 |
| Gradient Boosting | **0.928 ± 0.019** | 0.997 ± 0.000 | 2.75 ± 0.38 |
| Ridge Regression | 0.620 ± 0.045 | 0.635 ± 0.010 | 6.35 ± 0.25 |

特徴量重要度（Random Forest）：Ω_arag（~35%）> 気温（~28%）> 白化歴（~15%）> 光量（~10%）

Gradient Boosting の訓練R²=0.997に対するテストR²=0.928のギャップは軽度の過学習を示唆する。

![Figure 6: 機械学習モデル](figures/fig6_ml_model.png)

**図6.** 機械学習モデルの性能。(A) 特徴量重要度（Random Forest）、(B) モデル比較（5分割CV）、(C) 予測値対シミュレーション値、(D) 残差分布。

### 3.8 総括ダッシュボード

![Figure 7: 統合サマリーダッシュボード](figures/fig7_summary_dashboard.png)

**図7.** 統合モデル予測サマリーダッシュボード。SSPシナリオ別の主要指標と礁生存可能性フェーズ空間を表示。

---

## 4. 考察と今後の展望

### 4.1 主要な考察

**OAと昇温の相乗効果**: T=28°C以上でOAの影響が指数関数的に増幅される。Allison et al. (2021)の実験的知見と定性的に一致。SSP2-4.5シナリオでは2050〜2060年頃に50%石灰化閾値（T≈27.5°C、Ω≈2.7）を超える可能性が高い。

**生態系崩壊の非線形性**: 生態系モデルでは−89.1%〜−96.3%という急激な崩壊が生じた。これはLotka-Volterra系の臨界閾値を超えた際の相転移的挙動を反映しており、実際のサンゴ礁でも白化頻度の増加による累積ストレスが非線形崩壊を引き起こす可能性がある（Hughes et al., 2019）。

**進化的レスキュー**: 大規模礁（Ne≥5000）では適応対立遺伝子頻度が0.461に達するが、これは固定（頻度>0.95）には至らず、多数の個体が適応対立遺伝子を持ちながら礁全体が適応する「軟選択（soft selective sweep）」シナリオに対応する。Fuller et al. (2020)の多遺伝子スコア（polygenic score）アプローチと整合的。

### 4.2 自己批判的評価

1. **合成データへの依存**: 全学習データが同一の力学モデルから生成されているため、ML R²（0.928）は実世界予測精度を過大評価している可能性が高い。AIMS（オーストラリア海洋科学研究所）の長期モニタリングデータによる外部検証が必須。

2. **Lotka-Volterra系の過度の単純化**: 海草、CCA（石灰藻）、サンゴ幼生着底基盤のダイナミクスが含まれていない。算出された−96.3%崩壊は上限推定値として解釈すべきである。

3. **単遺伝子座仮定**: 白化耐性は6.8百万以上のSNPに分散した多遺伝子形質（Fuller et al., 2020）。単一二対立遺伝子座のWright-Fisher模型は進化速度を過大評価する可能性がある。

4. **NatureLM/GALACTICA MCPの不使用**: 定量予測の外部検証として本来使用すべきNatureLMとGALACTICAが利用不可であった。文献値による代替検証を行ったが、完全な代替にはなり得ない。

### 4.3 今後の展望

- **空間明示的モデル**: GBR北部・中部・南部の緯度勾配と海流による局所OAの不均一性を組み込んだ空間モデルへの拡張
- **CO2SYS/Atlantis統合**: フルバージョンのCO2SYS（pressure補正、Mg²⁺/Ca²⁺効果）とAtlantis生態系モデルとのカップリング
- **実データ検証**: AIMS 1995〜2023年の長期モニタリングデータを用いた後方検証（backcast validation）
- **共生微生物**: Symbiodiniaceaeのクレード変化（shuffling/switching）による部分適応の組み込み
- **マイクロ電解質**: サンゴの石灰化流体pHアップレギュレーション（proton pump）の明示的モデリング

---

## 5. 生成したファイル一覧

### Python実行スクリプト（Jupyter相当セル）
| セル番号 | 内容 | 主要出力 |
|---------|------|---------|
| Cell 1 | 環境セットアップ、乱数シード固定 | — |
| Cell 2/3 | 炭酸塩化学モジュール + Fig.1生成 | `figures/fig1_carbonate_chemistry.png`, `data/raw/calcification_rates.csv` |
| Cell 4b | 生態系ネットワークODE + Fig.2生成 | `figures/fig2_ecosystem_dynamics.png`, `data/raw/ecosystem_dynamics.csv` |
| Cell 5 | 複合ストレス解析 + Fig.3生成 | `figures/fig3_combined_stress.png`, `data/raw/synergy_analysis.csv` |
| Cell 6 | 集団遺伝学 + Fig.4生成 | `figures/fig4_population_genetics.png`, `data/raw/population_genetics.csv` |
| Cell 7 | GBR 2100年予測 + Fig.5生成 | `figures/fig5_gbr_projections.png`, `data/raw/gbr_2100_scenarios.csv` |
| Cell 8 | 機械学習モデル + Fig.6生成 | `figures/fig6_ml_model.png`, `data/raw/ml_training_data.csv` |
| Cell 9 | 総括ダッシュボード + Fig.7生成 | `figures/fig7_summary_dashboard.png` |

### 生成された図
```
figures/
├── fig1_carbonate_chemistry.png    — 炭酸塩化学と石灰化速度
├── fig2_ecosystem_dynamics.png     — 生態系動態（7種ODE）
├── fig3_combined_stress.png        — 温度-pH複合ストレス応答曲面
├── fig4_population_genetics.png    — 集団遺伝学・局所適応
├── fig5_gbr_projections.png        — GBR 2100年シナリオ予測
├── fig6_ml_model.png               — 機械学習モデル性能
└── fig7_summary_dashboard.png      — 統合サマリーダッシュボード
```

### 生成されたデータ
```
data/raw/
├── calcification_rates.csv         (1200行: 100 pCO₂ × 3温度 × 4種)
├── ecosystem_dynamics.csv          (2400行: 800タイムステップ × 3シナリオ)
├── synergy_analysis.csv            (24行: 6 T-Ω組合せ × 4種)
├── population_genetics.csv         (16行: 4 Ne × 4選択係数)
├── gbr_2100_scenarios.csv          (3行: SSPシナリオサマリー)
└── ml_training_data.csv            (800行: ML訓練データ)
```

### 論文・レポート
```
paper.md    — 学術論文形式（英語）
report.md   — 本実験レポート（日本語）
```

---

## 6. 再現性情報

```
Python: 3.11.2
numpy: 2.3.5
pandas: 2.3.3
scipy: 1.17.1
scikit-learn: 1.6.1
matplotlib: 3.10.9
seaborn: 0.13.2
xgboost: 3.2.0
lightgbm: 4.6.0
rdkit: 2026.3.2

乱数シード: np.random.seed(42), random.seed(42) — 全セルで統一
```

---

## 参考文献

1. Allison, N., et al. (2021). Resolving the interactions of ocean acidification and temperature on coral calcification media pH. *Coral Reefs*, 40, 1719–1729. https://doi.org/10.1007/s00338-021-02170-2
2. Boonnam, N., et al. (2022). Coral Reef Bleaching under Climate Change: Prediction Modeling and Machine Learning. *Sustainability*, 14(10), 6161. https://doi.org/10.3390/su14106161
3. Fuller, Z., et al. (2020). Population genetics of the coral Acropora millepora. *Science*, 369(6501), eaba4674. https://doi.org/10.1126/science.aba4674
4. González-Espinosa, P.C., & Donner, S. (2021). Cloudiness reduces the bleaching response of coral reefs. *Global Change Biology*, 27(16), 3837–3849. https://doi.org/10.1111/gcb.15676
5. Jagadeesh, M., & Pradhan, U. (2025). Coral Reef Bleaching Prediction: A Machine Learning Approach. https://doi.org/10.1109/ICSCDS65426.2025.11166756
6. Noonan, S., et al. (2025). Progressive changes in coral reef communities with increasing ocean acidification. *Commun. Biol.*, 8, 414. https://doi.org/10.1038/s42003-025-08889-w
7. Spreter, P.M., et al. (2022). Calcification response of reef corals to seasonal upwelling. *Biogeosciences*, 19, 3559–3573. https://doi.org/10.5194/bg-19-3559-2022
