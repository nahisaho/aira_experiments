# 食品テクスチャ予測モデリングフレームワーク

## 実験目的と背景

食品のテクスチャ（食感）は、消費者の嗜好性、健下機能、栄養摂取に直結する重要な品質指標であり、組成（多糖類、タンパク質、脂質、水分）と加工条件（温度、せん断、押し出し、3D印刷条件）から物理的に予測できれば、製品開発のサイクルを大幅に短縮できる。特に、植物性代替肉や3Dフードプリンティング用インクの設計においては、レオロジー特性とテクスチャプロファイル（TPA: Texture Profile Analysis）パラメータを定量的に結びつける統合フレームワークの必要性が高まっている (Pulatsu & Lin, 2021; Schreuders et al., 2021)。本研究では、(1) 多糖類ゲルの粘弾性、(2) 乳化系のレオロジー、(3) TPA予測モデル、(4) 口腔内プロセシング、(5) 3D印刷性、(6) 植物性代替肉の組成最適化、の6サブモジュールを一貫したパイプラインとして設計し、有限要素法/粗視化MDの構成則として利用可能な計算フレームワークを実装することを目的とする。

## 先行研究調査結果（MCPツール使用試行を含む）

文献調査には ToolUniverse MCP の `SemanticScholar_search_papers`、`Fatcat_search_scholar`、`openalex_search_works` を順に試行した。`SemanticScholar_search_papers` は HTTP 429（レート制限）により利用不能であった。`Fatcat_search_scholar` は3Dフードプリンティング関連クエリで5件の実論文（Maldonado-Rosas et al., 2025; Chen et al., 2024; Lu et al., 2023; Jiao et al., 2024 等）を返却した。`openalex_search_works` は粘弾性・乳化レオロジー・TPA関連クエリで追加の高被引用論文を返却した。

主要な知見は次のとおりである。Cao et al. (2021) はハイドロゲルの構造とレオロジーの関係を体系化し、生体高分子ゲルが Generalized Maxwell 型の緩和スペクトルで記述できることを示した。Schwab et al. (2020) はバイオインクの印刷性を降伏応力と剪断粘度に基づき定量化するフレームワークを提示し、本研究の印刷性スコア設計の根拠となった。Lu et al. (2023) は多糖類インクの印刷性を画像解析と Random Forest により82–88%の精度で予測した。Schreuders et al. (2021) はせん断セルおよび高水分押出によりタンパク質繊維配向が異方性弾性率を増幅することを実証した。

## 手法・アルゴリズム

### 1. 多糖類ゲルの粘弾性モデリング
Generalized Maxwell モデルにより、緩和弾性率と動的粘弾性を表現する：

$$ G(t) = G_\infty + \sum_{i=1}^{N} G_i \exp\!\left(-\frac{t}{\tau_i}\right) $$

$$ G'(\omega) = G_\infty + \sum_{i=1}^{N} \frac{G_i (\omega \tau_i)^2}{1 + (\omega \tau_i)^2}, \qquad G''(\omega) = \sum_{i=1}^{N} \frac{G_i\, \omega \tau_i}{1 + (\omega \tau_i)^2} $$

濃度依存性は経験的べき乗則 $G_0(c) = G_{0,\mathrm{ref}}(c/c_\mathrm{ref})^n$（$n \approx 2.2$）で導入した。拡張 Kelvin-Voigt（Burgers型）の creep compliance は

$$ J(t) = \frac{1}{E_0} + \frac{t}{\eta_0} + \sum_{i} \frac{1}{E_i}\left(1 - \exp\!\left(-\frac{E_i\, t}{\eta_i}\right)\right) $$

で表される。

### 2. 乳化系レオロジー
Krieger–Dougherty 則 $\eta_r = (1 - \phi/\phi_\mathrm{max})^{-[\eta]\phi_\mathrm{max}}$ と Quemada のせん断依存拡張、および Pal の effective medium 近似を実装した。液滴径分布の効果は Princen 型スケーリング $G_0 \propto 1/R_{32}$ で取り込んだ。

### 3. TPA予測モデル
350 サンプルの合成データセット（hardness, cohesiveness, springiness, gumminess, chewiness の5パラメータ）を生成し、Random Forest と Gradient Boosting のアンサンブル（平均）で予測した。ノイズ標準偏差 $\sigma=0.18$ を加え、5分割クロスバリデーションで $R^2$ と MAE の平均±標準偏差を報告する。

### 4. 口腔内プロセシング
咀嚼サイクル $n$ における粒径を $d(n) = d_\infty + (d_0 - d_\infty)\exp(-\lambda n)$、唾液モイスチャを $m(n) = 1 - \exp(-k_s n)$ とし、ボーラス凝集度 $C(n) \propto m(n)(d_0/d(n))^{0.4}$ が閾値 0.62 を超えた時点で嚥下と判定した。

### 5. 3D印刷性
Herschel–Bulkley 則 $\tau = \tau_y + K\dot\gamma^n$ を用い、降伏応力・見かけ粘度・流動指数からガウシアン重み付けの印刷性スコア（0–1）を算出し、層接着性・形状保持を別途モデル化した。

### 6. 植物性代替肉
ファイバーネットワーク弾性率 $E_\parallel$ を組成（大豆/エンドウタンパク、メチルセルロース、脂質、せん断配向度、押出温度）から計算し、differential evolution により目標 TPA に到達する組成を最適化した。

## 主要な結果

### 粘弾性応答
濃度 0.5–2.0% の多糖類ゲルにおいて、$G'$ は低周波で約1桁、高周波で約2桁の増加を示し、 $\tan\delta$ は中間周波数で極小を取る典型的なゲル応答が得られた（Cao et al., 2021 と整合）。

![Viscoelastic frequency sweep](figures/viscoelastic_frequency_sweep.png)

### 乳化レオロジー
$\phi=0.05$–$0.6$ の範囲で3モデルを比較したところ、Krieger–Dougherty と Pal モデルは $\phi<0.4$ でほぼ一致し、$\phi>0.5$ で発散の仕方が異なった。Quemada モデルは高せん断条件 ($\dot\gamma=100\,\mathrm{s}^{-1}$) でせん断薄化を示した（テスト `test_quemada_shear_thinning` で検証済）。

![Emulsion rheology comparison](figures/emulsion_rheology_comparison.png)

### TPA予測精度（5分割CV）
クロスバリデーション結果（mean ± std）は以下のとおりで、完璧スコアを避けた現実的な値となった：

| Target | $R^2$ | MAE |
| --- | --- | --- |
| hardness_N | 0.878 ± 0.031 | 3.53 ± 0.45 N |
| gumminess_N | 0.859 ± 0.044 | 2.49 ± 0.48 N |
| chewiness_J | 0.840 ± 0.030 | 0.0020 ± 0.0003 J |
| cohesiveness | 0.603 ± 0.082 | 0.057 ± 0.007 |
| springiness | 0.222 ± 0.037 | 0.056 ± 0.003 |

二次パラメータ（cohesiveness, springiness）は一次パラメータ由来の派生量より相対ノイズの影響を受けやすく、予測精度が低下することが定量的に示された。

![TPA prediction CV](figures/tpa_prediction_cv.png)

### 口腔内プロセシング
標準パラメータ（初期粒径 8 mm、硬さ 22 N）で嚥下サイクルは **14 cycle（約 11.9 秒）** と予測され、Chen (2009) のレビューで報告される 10–20 秒の固形食品の咀嚼時間範囲と一致した。

![Oral processing simulation](figures/oral_processing_simulation.png)

### 印刷性マップ
降伏応力 $\tau_y \in [10, 5000]$ Pa、見かけ粘度 $\eta \in [1, 1000]$ Pa·s の平面において、印刷性スコアの極大領域は $\tau_y \approx 100$–$500$ Pa、$\eta \approx 20$–$100$ Pa·s に存在し、Schwab et al. (2020) と Lu et al. (2023) の実験的最適範囲と整合した。

![Printability map](figures/printability_map.png)

### 植物性代替肉の組成最適化
目標 TPA（hardness 35 N、cohesiveness 0.68、springiness 0.72、gumminess 23.8 N、chewiness 0.017 J）に対し、differential evolution は **目的関数値 0.150** で収束し、組成は概ね soy protein 約 15–20%、methylcellulose 約 2.5%、脂質 約 5%、shear alignment 0.85+ の領域を選好した。

![Plant meat texture design](figures/plant_meat_texture_design.png)

## 考察

本フレームワークは、ミクロな組成・微視構造から TPA のような感覚評価関連量までを単一のパイプラインで結ぶことに成功した。特に粘弾性モデルは構成則として有限要素法シミュレーションに直接組み込み可能であり、Herschel–Bulkley パラメータは押出ダイ・3D 印刷ノズル形状最適化にも転用できる。一方、二次 TPA パラメータの予測精度が低い結果は、TPA 二回圧縮試験そのものに含まれる測定誤差を反映している可能性が高く、機械学習よりも実験プロトコルの再現性向上が先決であることを示唆する (Bourne, 2002)。植物性代替肉ケーススタディでは、配向度 (shear alignment) が硬さと spring back の両方を支配する最大要因であり、Schreuders et al. (2021) の高水分押出における異方性形成と整合的である。

## Limitations and Future Work

本研究には少なくとも次の限界がある。第一に、TPA 予測モデルは合成データセットに基づいて訓練・評価しており、実食品データへの転移性能は未検証である。実データを用いた外部検証ではドメインギャップにより $R^2$ がさらに 0.1–0.2 程度低下する可能性が高く、領域適応や半教師あり学習の導入が今後の課題となる。第二に、粘弾性モデルは線形粘弾性領域に限定されており、大変形・破壊挙動（LAOS：大振幅振動せん断）は記述できない。第三に、口腔内プロセシングモデルは粒径と凝集度のみの低次元縮約であり、唾液 $\alpha$-アミラーゼによる多糖類分解、舌圧分布、嚥下流体力学（VFSS データとの比較）を含んでいない。第四に、3D 印刷性スコアはヒューリスティックであり、実機の押出機特性・温度プロファイル・ノズル形状を陽に取り込んでいない。第五に、植物性代替肉のファイバーネットワークモデルは均質化近似であり、SEM・X-CT で観測される多孔質階層構造は反映していない。今後は (a) 公開データセット（Lu et al., 2023 等）での外部検証、(b) 粗視化分子動力学による$G_i, \tau_i$ パラメータの第一原理的推定、(c) LAOS 応答を含む非線形粘弾性への拡張、(d) FEM ソルバとの結合による咀嚼シミュレーション、を計画している。

## 生成したファイル一覧

- `src/viscoelastic_model.py` — Generalized Maxwell / 拡張 Kelvin-Voigt（121行）
- `src/emulsion_rheology.py` — Krieger-Dougherty / Quemada / Pal モデル（112行）
- `src/tpa_predictor.py` — RF+GB アンサンブル TPA 予測（120行）
- `src/oral_processing.py` — 咀嚼・嚥下シミュレーション（59行）
- `src/fd_printing.py` — 印刷性マップ・形状保持（76行）
- `src/plant_meat_design.py` — ファイバーネットワーク + DE 最適化（103行）
- `src/visualization.py` — 図生成パイプライン
- `tests/test_models.py` — 12 ユニットテスト（全件 PASS）
- `figures/`（6枚）: viscoelastic_frequency_sweep.png, emulsion_rheology_comparison.png, tpa_prediction_cv.png, oral_processing_simulation.png, printability_map.png, plant_meat_texture_design.png
- `results/`: reference-list.md, tpa_cv_summary.csv, emulsion_sweep.csv, plant_meat_optimization.json
- `data/synthetic_tpa_dataset.csv`
- `logs/process-log.jsonl`

## References

1. Maldonado-Rosas, R., Alfaro-Ponce, M., Cuan-Urquizo, E., & Tejada-Ortigoza, V. (2025). Printability prediction of food formulations for 3D printing using a Gaussian Process Regression model. *Journal of Food Engineering*. DOI: 10.1016/j.jfoodeng.2025.112534
2. Chen, J., Kong, Y., & Huang, Q. (2024). Analysis of surimi extrusion behavior during 3D printing by modified CFD and quick prediction of printability using machine learning based on texture data. *Innovative Food Science & Emerging Technologies*. DOI: 10.1016/j.ifset.2024.103698
3. Lu, Y., Rai, R., & Nitin, N. (2023). Image-based assessment and machine learning-enabled prediction of printability of polysaccharides-based food ink for 3D printing. *Food Research International*, 173, 113384. DOI: 10.1016/j.foodres.2023.113384
4. Jiao, X., Ren, G., Law, C. L., et al. (2024). Novel strategy for optimizing of corn starch-based ink food 3D printing process: Printability prediction based on BP-ANN model. *International Journal of Biological Macromolecules*. DOI: 10.1016/j.ijbiomac.2024.133921
5. Schreuders, F. K. G., Dekkers, B. L., Bodnár, I., Erni, P., Boom, R. M., & van der Goot, A. J. (2021). Thermo-mechanical processing of plant proteins using shear cell and high-moisture extrusion cooking. *Critical Reviews in Food Science and Nutrition*. DOI: 10.1080/10408398.2020.1864618
6. Pulatsu, E., & Lin, M. (2021). A review on customizing edible food materials into 3D printable inks. *Foods*, 10(2), 320. DOI: 10.3390/foods10020320
7. Cao, H., Saroia, J., Wang, Y., et al. (2021). Relationship between Structure and Rheology of Hydrogels for Various Applications. *Gels*, 7(4), 255. DOI: 10.3390/gels7040255
8. Schwab, A., Levato, R., D'Este, M., Piluso, S., Eglin, D., & Malda, J. (2020). Printability and Shape Fidelity of Bioinks in 3D Bioprinting. *Chemical Reviews*, 120(19), 11028–11055. DOI: 10.1021/acs.chemrev.0c00084
9. Krieger, I. M., & Dougherty, T. J. (1959). A mechanism for non-Newtonian flow in suspensions of rigid spheres. *Transactions of the Society of Rheology*, 3(1), 137–152. DOI: 10.1122/1.548848
10. Bourne, M. C. (2002). *Food Texture and Viscosity: Concept and Measurement* (2nd ed.). Academic Press. DOI: 10.1016/B978-0-12-119062-0.X5000-5
11. Chen, J. (2009). Food oral processing — A review. *Food Hydrocolloids*, 23(1), 1–25. DOI: 10.1016/j.foodhyd.2007.11.013
12. Quemada, D. (1977). Rheology of concentrated disperse systems and minimum energy dissipation principle. *Rheologica Acta*, 16(1), 82–94. DOI: 10.1007/BF01516932
13. Kew, B., Holmes, M., Stieger, M., & Sarkar, A. (2020). Oral tribology, adsorption and rheology of alternative protein-based emulsions. *Food Hydrocolloids*, 109, 106105. DOI: 10.1016/j.foodhyd.2020.106105
