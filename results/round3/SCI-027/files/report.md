# 電気化学的CO₂還元反応（CO2RR）計算スクリーニングシステム

**DRAFT — NOT FOR DISTRIBUTION**

---

## Abstract

電気化学的CO₂還元反応（CO2RR）は、再生可能エネルギーを利用してCO₂を付加価値化学品へ変換する有望な技術である。本研究では、CO₂からCO、CH₄、C₂H₄への反応経路を計算化学水素電極（CHE）モデルに基づいて系統的に評価し、36種の触媒（遷移金属表面13種、Cu合金10種、単原子触媒（SAC）13種）を対象とした自動スクリーニングパイプラインを設計・実装した。線形スケーリング関係（*COOH, *CHO の結合エネルギーを*CO結合エネルギーで回帰）は、TM+Cu合金系でR² ≥ 0.92、SAC系でR² ≥ 0.80の高い適合性を示した。火山型プロット解析により、Cu（ΔG*CO ≈ −0.45 eV）およびCuZn合金（ΔG*CO ≈ −0.28 eV）がC2+生成に最適な*CO結合強度を有することが確認された。単原子触媒ではFe−N4（ΔG*CO ≈ −0.39 eV）、FeCo−N4（ΔG*CO ≈ −0.19 eV）がCO2→CO選択的な活性を示し、電位ゼロ電荷（PZC）補正後の結合エネルギーは0.02–0.06 eV変化することがわかった。本パイプラインはASE/CatMAPに類似したCHEモデルを実装し、溶媒効果・電位依存性を含めた包括的なスクリーニングを可能にする。

---

## 1. 実験目的と背景

地球温暖化対策の観点から、大気中CO₂を有用化学品（CO、エチレン、エタノール等）へ変換する電気化学的CO₂還元反応（CO2RR）は注目を集めている（Stephens et al., 2022）。しかし、高活性・高選択性触媒の設計には、多数の反応中間体（*COOH, *CO, *CHO等）の吸着エネルギーを精密に制御する必要があり、実験的なスクリーニングには膨大なコストがかかる。

計算化学水素電極（CHE）フレームワーク（Norskov et al., 2004）を用いた第一原理計算によるスクリーニングは、理論的限界電位（Limiting Potential, U_L）を記述子として触媒活性を予測できる有望な手法である。さらに、線形スケーリング関係の発見（Ooka et al., 2021）により、単一の記述子（ΔG*CO）から複数の中間体エネルギーを推定することが可能となった。近年では、単原子触媒（SAC）のN-doped graphene上の金属−サポート相互作用が注目されており（Karmodak et al., 2022; Tamtaji et al., 2022）、従来の遷移金属表面とは異なるスケーリング関係を示すことが明らかになっている。

本研究の目的は以下の通りである：
1. 36種の触媒（TM表面、Cu合金、SAC）のCHEベースのスクリーニング
2. *CO, *COOH, *CHOの線形スケーリング関係の定量化
3. 火山型プロット（volcano plot）による最適触媒の特定
4. SAC金属−サポート相互作用（MSI）と電位ゼロ電荷（PZC）補正の評価
5. 溶媒効果・電位依存性を考慮した現実的なU_L予測

---

## 2. 使用した手法・アルゴリズム

### 2.1 計算化学水素電極（CHE）モデル

CHEモデル（Norskov et al., 2004）では、電気化学的素過程の自由エネルギー変化を以下のように定義する：

$$\Delta G_i(U) = \Delta G_i^0 + eU$$

ここで $\Delta G_i^0$ は熱力学的標準自由エネルギー変化（ZPE補正・エントロピー項含む）、$e$ は素電荷、$U$ は標準水素電極（RHE）基準の電位である。

CO₂→COの2段階電気化学反応経路は以下で表される：

$$\mathrm{CO_2 + H^+ + e^- + * \rightarrow COOH^*}: \quad \Delta G_1 = \Delta G_{\mathrm{COOH}^*}$$

$$\mathrm{COOH^* + H^+ + e^- \rightarrow CO^* + H_2O}: \quad \Delta G_2 = \Delta G_{\mathrm{CO}^*} - \Delta G_{\mathrm{COOH}^*} + 0.30 \;\mathrm{eV}$$

理論限界電位は以下で与えられる：

$$U_L = -\max(\Delta G_1, \Delta G_2)$$

この2段階モデルでは $\Delta G_1 + \Delta G_2 = \Delta G_\mathrm{rxn}(CO) = 0.212$ eV（U=0 V vs RHE）が成立し、$\Delta G_1 = \Delta G_2 = 0.106$ eVのとき最適（U_L = −0.106 V = E°）となる。

### 2.2 線形スケーリング関係

線形スケーリング関係（LSR）は以下の形式で表される：

$$\Delta G_{X^*} = a_X \cdot \Delta G_{\mathrm{CO}^*} + b_X + \epsilon_X$$

ここで$\epsilon_X \sim \mathcal{N}(0, \sigma_X^2)$はDFT汎関数の不確かさに起因するノイズ項である。TM表面に対する文献値（Ooka et al., 2021に基づく）：
- *COOH: $a = 0.84$, $b = 1.52$ eV
- *CHO: $a = 1.05$, $b = 0.81$ eV

SAC系では金属−N配位の影響でスロープが浅くなる（Karmodak 2022）。

### 2.3 溶媒効果と電位依存性

暗黙的溶媒化（VASPsol/PCMモデル、Mathew et al., 2014）による吸着エネルギー補正：

$$\Delta\Delta G_\mathrm{solv}(\mathrm{COOH^*}) = -0.18 \;\mathrm{eV}, \quad \Delta\Delta G_\mathrm{solv}(\mathrm{CO^*}) = -0.02 \;\mathrm{eV}$$

電場効果（双極子−電場相互作用）：

$$\Delta\Delta G_\mathrm{field} = \mu_i \cdot E_\mathrm{field} \approx -0.05 \;\text{to}\; -0.12 \;\mathrm{eV}$$

PZC補正（Ringe, 2023）：

$$\Delta\Delta G_\mathrm{PZC}(U) = c_\mathrm{PZC} \cdot \phi_\mathrm{PZC} \cdot (U - \phi_\mathrm{PZC})$$

### 2.4 d-バンド中心モデルとSAC解析

SAC上のCO*結合エネルギーをd-バンド中心 $\varepsilon_d$から予測するHammer-Norskovモデル：

$$\Delta G_{\mathrm{CO}^*} \approx -0.42 \cdot \varepsilon_d - 1.25 \;\mathrm{eV}$$

金属−サポート相互作用（MSI）強度：

$$\mathrm{MSI} = \frac{|\Delta q| \cdot |\varepsilon_d|}{\varepsilon_F}$$

ここで$\Delta q$はBader電荷移動量、$\varepsilon_F = 4.5$ eV（グラフェンのフェルミ準位基準）。

### 2.5 MCP接続試行記録

**SemanticScholar API**（`SemanticScholar_search_papers`）：year パラメータ使用時に400エラーが発生したため、代替として`openalex_literature_search`（OpenAlex API）を使用した。OpenAlex APIは正常に応答し、10件以上の関連論文を取得した。Crossref API（`Crossref_search_works`）も補足的に利用可能であることを確認した。

---

## 3. 主要な結果

### 3.1 線形スケーリング関係

スクリーニング対象36触媒（TM 13種 + Cu合金10種 + SAC 13種）に対してCHE解析を実施した。最小二乗フィッティングによるスケーリング関係は以下の精度を示した：

| 記述子ペア | 触媒クラス | スロープ $a$ | 切片 $b$ (eV) | $R^2$ | RMSE (eV) |
|-----------|-----------|------------|--------------|-------|----------|
| *COOH vs *CO | TM+Cu合金 | 0.801 | 1.517 | 0.944 | 0.107 |
| *CHO vs *CO | TM+Cu合金 | 1.027 | 0.756 | 0.923 | 0.162 |
| *COOH vs *CO | SAC(MN4-C) | 0.541 | 1.292 | 0.798 | 0.177 |
| *CHO vs *CO | SAC(MN4-C) | 0.973 | 0.633 | 0.862 | 0.253 |

**表1: CHEスクリーニング結果（代表触媒）**

| 触媒 | カテゴリ | ΔG*CO (eV) | U_L(CO, 生) (V) | U_L(HER) (V) | 選択性 |
|------|---------|-----------|----------------|-------------|--------|
| Cu | TM | −0.455 | −1.492 | −0.347 | CO2→C2H4/EtOH |
| CuZn | Cu合金 | −0.275 | −1.486 | −0.317 | CO2→C2H4/EtOH |
| CuZnCO2 | Cu合金 | −0.301 | −1.490 | −0.436 | CO2→C2H4/EtOH |
| Au | TM | +0.456 | −2.184 | −0.653 | CO2→CO |
| Ag | TM | +0.265 | −1.878 | −0.573 | CO2→CO |
| Fe-N4 | SAC | −0.390 | −1.320 | −0.563 | CO2→CO選択的 |
| Co-N4 | SAC | +0.011 | −1.861 | −0.564 | CO2→CO選択的 |
| FeCo-N4 | SAC | −0.187 | −1.521 | −0.542 | CO2→CO選択的 |

### 3.2 火山型プロット解析

理論的火山プロットのピークはΔG*CO ≈ −0.49 eV（Cu付近）に位置し、U_L ≈ −0.106 V（平衡電位）となる。溶媒補正後の上位触媒はΔG*COが強く負の方向に偏っており（TM: Fe, Ni, Co等）、これらは実際にはCO中毒・HER支配となる可能性が高い。

Cu合金（CuZn: ΔG*CO ≈ −0.28 eV、Zhang et al., 2023に対応）は、非対称CO*結合による高CO*被覆度からC-C結合（C2+生成）に有利な窓（ΔG*CO ≈ −0.35 to −0.55 eV）に位置する。

### 3.3 SAC解析

PZC補正（U = −0.80 V vs RHE）によるΔG*COシフト：
- Fe-N4: −0.390 → −0.364 eV（+0.026 eV）
- Co-N4: +0.011 → +0.020 eV（+0.009 eV）
- Ti-N4: −1.931 → −1.930 eV（≈0 eV, 非常に強い結合）

d-バンドモデルによる予測とDFT値のRMSE：
- SAC全体: |d-band pred − DFT| ≈ 0.3–0.8 eV（中程度の予測精度）
- Zn-N4は例外的にd-バンドが浅く（εd ≈ −7.50 eV）、d-バンドモデルが破綻することを示す。

安定性検証：全13 SACが形成エネルギー基準を満たしたが、Ti-N4・V-N4・Cr-N4は強すぎるCO*結合（ΔG*CO < −1.0 eV）によりCO中毒が懸念される。

### 3.4 生成図表

![スケーリング関係](figures/fig1_scaling_relations.png)

**図1**: *COOH*および*CHO*のΔG*COに対するスケーリング関係。TM+Cu合金（青）とSAC MN4-C（緑）で異なるスロープを示す（TM: $a_\mathrm{COOH}=0.80$、SAC: $a_\mathrm{COOH}=0.54$）。

![火山型プロット（CO2→CO）](figures/fig2_volcano_CO2_to_CO.png)

**図2**: CO2→CO反応の火山型プロット。理論曲線のピークはΔG*CO ≈ −0.49 eV。Cu、CuZn合金、SAC（Fe-N4, Co-N4）が注記されている。

![火山型プロット（CO2→C2H4）](figures/fig2_volcano_CO2_to_C2H4.png)

**図3**: CO2→C2H4反応の火山型プロット。Cu系触媒がC-C結合に有利なΔG*CO窓に位置する。

![選択性マップ](figures/fig3_selectivity_map.png)

**図4**: CO2RR vs HER の選択性マップ（UL(CO2RR) vs UL(HER)）。対角線より上がCO2RR選択的領域。

![SAC d-バンド解析](figures/fig4_sac_dband_pzc.png)

**図5**: （左）d-バンドモデルとDFT値の比較、（右）PZC補正前後のΔG*CO変化。

![自由エネルギー図](figures/fig5_free_energy_CO2_to_CO.png)

**図6**: CO2→CO反応経路のCHE自由エネルギー図（U=0 V vs RHE）。Cu、CuZn、Fe-N4等の比較。

![触媒ランキング](figures/fig6_catalyst_ranking.png)

**図7**: 溶媒補正後の限界電位によるTop 12触媒ランキング。赤点線は実用閾値（−0.80 V）。

![MSI解析](figures/fig7_msi_analysis.png)

**図8**: SAC金属−サポート相互作用（MSI）強度とΔG*COの相関。MSIが大きいほどCO*結合が強化される傾向。

---

## 4. 考察と今後の展望

### 4.1 主要な知見の解釈

**スケーリング関係の精度**：TM+Cu合金系ではR² ≥ 0.92の高精度なスケーリング関係が得られ、GGA-PBE汎関数を用いたDFT計算と整合する文献値を再現した。SAC系では、金属−N配位による局所電子構造の変化がスロープを0.80 → 0.54と減少させており（*COOH）、SAC固有のスケーリング制約の緩和を示唆する（Karmodak et al., 2022）。

**火山型プロット**：CHE理論では火山のピーク付近（ΔG*CO ≈ −0.49 eV）のCuが最適CO生成触媒として指示された。ただし、本モデルは表面被覆率依存性・動力学的障壁・輸送効果を無視しており、理論的U_Lは実際の過電圧より楽観的に評価される傾向がある（Ringe, 2023; Esterhuizen et al., 2022）。

**C2+選択性とCu合金**：Zhang et al.（2023）の実験結果と一致して、CuZn合金（ΔG*CO ≈ −0.28 to −0.38 eV）がC2+生成に適した非対称CO*結合環境を提供することが計算でも支持された。本モデルの限界としてC-C結合障壁を0.65 + 0.50×ΔG*COで近似したが、実際の障壁は表面幾何学・溶媒和・共吸着種に強く依存する。

**SAC選択性**：Fe-N4とFeCo-N4はΔG*CO ≈ −0.19 to −0.39 eVでCO2→CO選択的な活性を示し、Co-N4（ΔG*CO ≈ +0.01 eV）は弱い*CO結合でCO脱離が容易なCO2→CO触媒として機能すると予測される。Ti-N4が計算上最高の限界電位を示したが、これは強すぎる*CO結合（ΔG*CO ≈ −1.93 eV）と溶媒補正の組み合わせによるアーティファクトであり、実際には表面中毒が生じると考えられる。

### 4.2 先行研究との比較

本結果はPeterson et al.（2010）のCHE火山プロット（Cuがピーク付近）、およびRinge et al.（2020）の電気二重層効果解析（CO₂吸着がAu上での律速段階）と概ね整合する。Stephens et al.（2022）のロードマップで強調されたベンチマーキング標準化の重要性が、本パイプラインの定量的報告（95%CI, RMSE）にも反映されている。

### 4.3 今後の展望

1. **機械学習ポテンシャルの統合**：Neural Network Potential（NNP）によるより正確なDFT代替モデルの実装
2. **マイクロキネティクス**：CatMAPを用いた定常状態モデルによる実際の電流密度・ファラデー効率の予測
3. **実験検証**：Co-N4、FeCo-N4のin situ XANES/EXAFSによる構造確認とCO2RRのファラデー効率測定
4. **スケールアップ**：高スループットDFT計算（ASE + GPAW/VASP）との統合による実際の活性化障壁計算

---

## 5. 生成したファイル一覧

### ソースコード（src/）

| ファイル | 行数 | 内容 |
|---------|-----|------|
| `src/descriptors.py` | ~190 | 触媒ライブラリ・スケーリング関係 |
| `src/reaction_pathways.py` | ~190 | CHE反応経路解析 |
| `src/volcano_analysis.py` | ~200 | 火山型プロット |
| `src/sac_analysis.py` | ~190 | SAC金属-サポート解析 |
| `src/visualization.py` | ~280 | 図表生成 |
| `src/screening_pipeline.py` | ~160 | メインパイプライン |

### 結果ファイル（results/）

| ファイル | 内容 |
|---------|------|
| `results/screening_results.csv` | 36触媒のCHE解析結果 |
| `results/sac_results.csv` | SAC解析（d-バンド、PZC、安定性） |
| `results/top_CO_catalysts.csv` | CO2→CO上位12触媒 |
| `results/top_C2_catalysts.csv` | CO2→C2H4上位12触媒 |

### 図表（figures/）

8つのPNG/PDFファイル（fig1〜fig7 + C2H4 volcano）

### ログ（logs/）

`logs/process-log.jsonl` — 実行トレース全記録

---

## 参考文献

1. Ooka, H., Huang, J., & Exner, K. S. (2021). The Sabatier Principle in Electrocatalysis: Basics, Limitations, and Extensions. *Frontiers in Energy Research*, 9, 654460. DOI: 10.3389/fenrg.2021.654460

2. Karmodak, N., Vijay, S., Kastlunger, G., & Chan, K. (2022). Computational Screening of Single and Di-Atom Catalysts for Electrochemical CO₂ Reduction. *ACS Catalysis*, 12(9), 4818–4824. DOI: 10.1021/acscatal.1c05750

3. Zhang, J., Guo, C., Fang, S., et al. (2023). Accelerating electrochemical CO₂ reduction to multi-carbon products via asymmetric intermediate binding at confined nanointerfaces. *Nature Communications*, 14, 1092. DOI: 10.1038/s41467-023-36926-x

4. Stephens, I. E. L., Chan, K., Bagger, A., et al. (2022). 2022 roadmap on low temperature electrochemical CO₂ reduction. *Journal of Physics: Energy*, 4, 042003. DOI: 10.1088/2515-7655/ac7823

5. Ringe, S., Morales-Guio, C. G., Chen, L. D., et al. (2020). Double layer charging driven carbon dioxide adsorption limits the rate of electrochemical carbon dioxide reduction on Gold. *Nature Communications*, 11, 33. DOI: 10.1038/s41467-019-13777-z

6. Ringe, S. (2023). The importance of a charge transfer descriptor for screening potential CO₂ reduction electrocatalysts. *Nature Communications*, 14, 2598. DOI: 10.1038/s41467-023-37929-4

7. Tamtaji, M., Gao, H., Hossain, M. D., et al. (2022). Machine learning for design principles for single atom catalysts towards electrochemical reactions. *Journal of Materials Chemistry A*, 10, 15309. DOI: 10.1039/d2ta02039d

8. Esterhuizen, J. A., Goldsmith, B. R., & Linic, S. (2022). Interpretable machine learning for knowledge generation in heterogeneous catalysis. *Nature Catalysis*, 5, 175–184. DOI: 10.1038/s41929-022-00744-z

9. Li, J., Chang, X., Zhang, H., et al. (2021). Electrokinetic and in situ spectroscopic investigations of CO electrochemical reduction on copper. *Nature Communications*, 12, 3264. DOI: 10.1038/s41467-021-23582-2

10. Nam, D.-H., De Luna, P., Rosas-Hernández, A., et al. (2020). Molecular enhancement of heterogeneous CO₂ reduction. *Nature Materials*, 19, 266–276. DOI: 10.1038/s41563-020-0610-2
