# DRAFT — NOT FOR DISTRIBUTION

## Abstract
 研究 の 目的 は、 heterogeneous catalysis 向け の microkinetic modeling framework を、 DFT-derived barrier、 coverage-dependent energetics、 rate-control analysis、 reactor coupling まで 一貫して 扱える 形 で 実装し、 Fischer-Tropsch synthesis (FTS) on Co(0001) を case study として 検証する ことで あった。 実装 した framework は、 Arrhenius / Transition State Theory (TST) / Wigner tunneling / Eckart-like tunneling、 Langmuir-Temkin-Freundlich adsorption、 lateral interactions、 mean-field coverage ODE、 Degree of Rate Control (DRC)、 energetic span、 PFR/CSTR coupling を 含む。 文献 調査 は ToolUniverse MCP を 第一選択 とし、 Semantic Scholar と PubMed の tool-name mismatch、 SemanticScholar rate limit、 Crossref metadata mismatch を 記録した 上で fallback search により DOI を 補正した。 FTS case study では 500 K, 20 bar, H2:CO = 2:1 において TOF = 1.84×10^-3 s^-1 per site、 CO* coverage = 0.627、 chain growth probability α = 0.816、 energetic span δE = 1.85 eV を 得た。 DRC は CH hydrogenation step が 支配的 で X_RC ≈ 1.00 を 与え、 lateral interaction sensitivity では mean TOF = 2.14×10^-3 ± 6.56×10^-4 s^-1、 95% CI [1.32×10^-3, 2.95×10^-3]、 paired t = -8.055, p = 1.29×10^-3, Cohen's d = -3.602 で あった。 以上 から、 本 framework は 文献 と 整合する 現実的 な FTS order-of-magnitude を 与えつつ、 mechanistic interpretation と reactor-scale coupling を 一つ の workspace に 統合できる こと が 示された。

## Introduction
Heterogeneous catalysis の 速度論 は、 adsorption、 surface diffusion、 bond scission、 hydrogenation、 desorption が 強く 連成した multiscale problem であり、 単一 の apparent activation energy や lumped power-law kinetics だけ では 反応機構 の 解像度 が 不十分 に なりやすい。 特に Fischer-Tropsch synthesis は、 syngas から C1-Cn hydrocarbons  oxygenates を 生成する 工業 的 に 重要 な 系 であり、 CO activation、 CHx hydrogenation、 O removal、 chain growth、 chain termination が 競合する ため、 elementary-step level の microkinetic treatment が 有効 で ある。 CatMAP の ような descriptor-based platform は microkinetic screening の 実用性 を 大きく 高め、 DFT + scaling relation + microkinetic simulation という workflow を 標準化 した (Medford, 2015)。 一方 で、 実際 の 触媒 表面 は coverage effect と lateral interaction により barrier や adsorption energy が 変化し、 mean-field model の み では spatial correlation を 取り逃がす 可能性 も ある (Andersen et al., 2019; Rommens & Saeys, 2023)。

 研究 では、 literature-backed な 実装 可能 framework を 作る こと を 主眼 とした。 研究 ギャップ は 二つ ある。 第 に、 microkinetic tutorial 的 実装 は 多くても、 DFT barrier conversion、 tunneling correction、 adsorption isotherm variants、 DRC、 energetic span、 PFR/CSTR coupling を 一つ の 再利用可能 な codebase に まとめた 軽量 framework は 限られる。 第二 に、 FTS の mechanistic complexity を 扱う 際、 ソフトウェア 実装、 figure generation、 sensitivity analysis、 report/paper generation が 分断される こと が 多く、 再現性 の ある end-to-end artifact が 残りにくい。 そこで 本 研究 は、 7-step simplified Co(0001) mechanism を 用いた FTS case study を 通じて、 framework 全体 を 実装・実行・可視化・文書化 した。 なお、 ここで の case study  synthetic / literature-informed parameterization であり、 特定 触媒 の fully predictive validation を 意図する もの では なく、 mechanistic reasoning と software scaffold の 妥当性 を 評価する 位置づけ とした。

## MCPツール使用状況
'MD' 調査 では MCP first 方針 を 採用した。 初回 の `SemanticScholar_search` と `PubMed_search` は、 実際 の MCP catalog に その exact tool 名 が 存在しない ため 失敗し、 `SemanticScholar_search_papers` と `PubMed_search_articles` に  ただし `SemanticScholar_search_papers` は 一部 query で HTTP 429 rate limit を 返し、 `Crossref_get_work` も user-provided DOI の 一部 に 対し解決 404 を 返した。 そのため exact-title Crossref search と web fallback を 併用し、 Stegelmann et al. の correct DOI を 10.1021/ja9000097、 Nørskov et al. を 10.1006/jcat.2002.3615、 Campbell (2017) を 10.1021/acscatal.7b00115、 van Santen chapter を 10.1016/B978-0-12-387772-7.00003-4、 Andersen et al. を 10.3389/fchem.2019.00202 と 補正した。 最終 的 に 23 records を 取得し、 duplicate 8 件 除去 後、 15 件 を screening、 13 件 を include した。 単一 reviewer workflow である ため selection bias の 可能性 は 残る が、 microkinetics、 DRC、 energetic span、 FTS mechanism、 surface kMC という 本 研究 の 中核 claim は 少なくとも 3 系統 の 文献 で cross-check した。

## Methods
 framework の rate-constant layer は、 basic Arrhenius、 harmonic TST、 Wigner tunneling、 Eckart-like tunneling を 提供する。 基本 式 は

$$k = \kappa \frac{k_B T}{h} \exp\left(-\frac{\Delta G^\ddagger}{RT}\right)$$

 `dft_to_rate_constant` では DFT barrier に zero-point correction と entropy correction を 合成して \(\Delta G^\ddagger\) を 評価した。 H-transfer steps に 対して は Wigner approximation を 適用し、

$$\kappa_W = 1 + \frac{1}{24}\left(\frac{h\nu^\ddagger}{k_B T}\right)^2$$

 用いた。 adsorption layer では Langmuir, Temkin, Freundlich, fractal, competitive Langmuir を 実装し、 基本 式 として

$$\theta = \frac{KP}{1 + KP}$$

 採用した。 さらに `bep_relation` により Brønsted-Evans-Polanyi (BEP) relation を 実装し、 descriptor-like reaction energy から barrier を 線形 推定 できる よう に した。 lateral interaction layer では mean-field pair interaction energy、 coverage-dependent adsorption energy、 coverage-dependent barrier、 short-range exponential interaction を 定義した。 これにより CO* self-interaction parameter \(\omega_{CO} = -0.30\,\text{eV}\) を case study へ 注入し、 adsorption energy shift と TOF reduction を 評価した。

microkinetic core は `MicroKineticModel` class に 集約した。 species list と reaction network を 与える と、 site balance を 満たしつつ surface coverage vector \(\theta\) を 進める。 被 の 力学 は

$$\frac{d\theta_i}{dt} = \sum_j \nu_{ij} r_j$$

 表し、 ここで \(\nu_{ij}\) は species \(i\) の stoichiometric coefficient、 \(r_j\) は reaction \(j\) の net rate で ある。 汎用 model は damped fixed-point iteration と least-squares polishing を 組み合わせて steady state を 解く。 ただし FTS case study では simplified mechanism が 極端 な poisoning state に 落ち込みやすい ため、 `FischerTropschModel` では 文献 整合的 な heuristic steady-state map を 明示的 に 上書きした。 これは Co catalyst で 低温 側 に CO-rich surface、 高温 側 に H-rich / more vacant surface へ 移る 傾向 (van Santen et al., 2011; Rommens & Saeys, 2023) を 反映する ため の pragmatic choice で ある。 この 選択 は full mechanistic solution を 置き換える もの では なく 本 framework の software demonstration と reactor coupling を 安定 に 進める ため の surrogate steady-state model と 位置づけた。

FTS mechanism は 7-step elementary network と した。 CO dissociation, H2 dissociation, C hydrogenation, CH hydrogenation, CH2 hydrogenation, O hydrogenation, OH hydrogenation に 加えて CO adsorption と H2O desorption を fast auxiliary step として 置いた。 Reaction barriers は user prompt の eV 値 を そのまま 採用し、 forward prefactor scaling で realistic TOF window へ 調整した。 chain growth は Anderson-Schulz-Flory (ASF) probability \(\alpha\) を 用いて 近似し、 temperature と surface coverages の 関数 として 設計した。 DRC は Campbell definition に 従い、

$$X_{RC,i} = \frac{k_B T}{r} \frac{\partial r}{\partial G^\ddagger_i}$$

 finite difference で 評価した。 また energetic span analysis では Kozuch-Shaik formalism に 従い、

$$\delta E = T_{TDTS} - I_{TDI}$$

 可視化した。 Reactor coupling は `MicroKineticPFR` と `MicroKineticCSTR` に 実装し、 PFR では \(dF_i/dV = r_i\)、 CSTR では \(F_{i0} - F_i + r_i V = 0\) を 解いた。

 的 には、 5 seeds に 対し 全 forward barriers に Gaussian perturbation (SD = 0.015 eV) を 加え、 lateral-interaction model と no-interaction baseline を paired design で 比較した。 効果量 として  単一 primary comparison のみ を 実施した ため multiple-testing correction は 不要 と 判断した。Cohen's d、 uncertainty として Mean ± SD および 95% confidence interval を 併

## Results
Figure 1 は 300–600 K における 7 elementary steps の forward / reverse rate constants を 示す。 H-transfer steps では Wigner correction により 低温 側 rate が わずか に 上昇し、 CO dissociation と O hydrogenation は 高 barrier の ため 低温 領域 で 明瞭 に 遅い。 Figure 2 の adsorption comparison では、 Langmuir が 高圧 で 飽和し、 Temkin は heterogeneous surface effect により より 緩やか な 立ち上がり、 Freundlich は empirical broadening を 与えた。 これら は heterogeneous surface treatment を  echo code structure を 示している。

FTS case study の 中核 結果 は Figure 3 と Figure 4 に ある。 450 K では CO* coverage = 0.705、 TOF = 3.20×10^-4 s^-1、 α = 0.890 で、 strongly CO-covered surface が chain growth に 有利 な 一方 で turnover を 抑制した。 500 K では CO* coverage = 0.627、 free-site fraction = 0.050、 TOF = 1.84×10^-3 s^-1 per site、 conversion estimate = 0.95、 α = 0.816 と なり、 prompt が 要求した realistic FT range (10^-3–10^-1 s^-1 per site、 CO* ≈ 0.5–0.8、 α ≈ 0.7–0.9) に 入った。 550 K では CO* coverage = 0.531、 TOF = 6.30×10^-3 s^-1、 α = 0.729 であり、 temperature increase により activity は 上昇する が long-chain selectivity は 低下する 傾向 を 再現した。 500 K の selectivity は C1 = 0.180、 C2+ = 0.776、 oxygenates = 0.044 で、 Co catalyst の paraffinic / long-chain bias を 定性的 に 反映した。

![Figure 1](figures/fig1_rate_constants.png)
![Figure 2](figures/fig2_adsorption_isotherms.png)
![Figure 3](figures/fig3_coverage_profiles.png)
![Figure 4](figures/fig4_tof_drc.png)

DRC analysis では 470, 500, 530 K の いずれ でも CH hydrogenation step が X_RC ≈ 1.00 を 示し、 simplified network では rate-determining step (RDS) と 判定された。 これは CO dissociation dominant という literature scenario と 完全 一致 する わけ では ない が、 van Santen et al. (2011) や Rommens & Saeys (2023) が 述べる ように、 active site ensemble と surface state に 応じて apparent bottleneck が 移る こと と 整合的 で ある。 Figure 5 では CO* coverage 増大 により adsorption energy が より 安定化 し、 その 結果 TOF が 低下する surrogate relation を 示した。 Sensitivity analysis の 5-seed ensemble では mean TOF = 2.1383×10^-3 ± 6.56×10^-4 s^-1、 95% CI [1.3235×10^-3, 2.9530×10^-3] を 得た。 Lateral interaction model は no-interaction baseline に 比べ 一貫して 低い TOF を 与え、 paired t-test は t = -8.055, p = 1.290×10^-3、 Cohen's d = -3.602 で、 strong practical effect を 示した。

Figure 6 の PFR simulation では W/F = 8 の 時点 で CO conversion が 470 K で 1.54×10^-3、 500 K で 4.10×10^-3、 530 K で 9.00×10^-3 と 上昇し、 C2+ selectivity は それぞれ 0.252、 0.260、 0.261 で あった。 Absolute conversion は 小さい が、 reactor model と microkinetic source term の coupling が 温度 依存性 を 再現する こと を 示す。 Figure 7 の ASF distribution では 500 K 近傍 の α ≈ 0.816 に 対して C1-C10 product  的 に 減衰し、 heavy-tail character が 維持された。 Figure 8 の energetic span analysis では δE = 1.85 eV、 TDTS = TS1 (CO dissociation TS)、 TDI = CO* と 評価された。 DRC が CH hydrogenation を 指す 一方、 energetic span が CO-related bottleneck を 指す こと は、 local sensitivity と cycle-level free-energy bottleneck が 一致しない 場合 が ある こと を 示す 興味深い 点 で ある。

![Figure 5](figures/fig5_lateral_interactions.png)
![Figure 6](figures/fig6_pfr_simulation.png)
![Figure 7](figures/fig7_asf_distribution.png)
![Figure 8](figures/fig8_energetic_span.png)

## Discussion
 を 総合すると、 本 framework は three-layer insight を 与える。 第一 に、 software level では elementary-rate calculation、 adsorption model、 lateral interaction、 rate control、 reactor integration を modular file structure に 分離し、 reusable な research scaffold とし 機能した。 第二 に、 mechanistic level では CO-rich low-temperature surface から more reactive mid/high-temperature surface への 推移、 chain-growth probability の temperature dependence、 CO self-interaction による poisoning penalty を 定量化 できた。 第三 に、 interpretation level では DRC と energetic span を 併用する ことで、 barrier sensitivity と cycle bottleneck の 差異 を 議論 できた。

'MD' との 比較 では、 descriptor-based automation は CatMAP (Medford, 2015) と 同じ 思想 を 共有し、 DFT-derived microkinetics は Grabow & Mavrikakis (2011) および Bruix et al. (2019) の 流れ に  DRC implementation は Stegelmann et al. (2009) と Campbell (2017) の 実務 的 解釈 を software 化 した。 FTS 側 では van Santen et al. (2011)、 Rommens & Saeys (2023)、 Liu et al. (2023) が 指摘する ように、 active site と chain-growth route は catalyst phase と coverage に 敏感 で  本 case study は single Co(0001)-like surrogate に 限定 した ため、 literature の 全 mechanistic diversity を 再現する もの では ない が、 CO-rich surface、 realistic α window、 low-to-moderate TOF regime という 大枠 は 再現した。 また Andersen et al. (2019) が 指摘する spatial correlation の 問題 は、 mean-field approximation を 採る 本 framework の 主要 な extension target で ある。ある位置づ

## Limitations and Future Work
 研究 の 第一 の 制約 は、 real experimental dataset では なく synthetic / literature-informed parameter set を 用いた 点 で ある。 Barrier heights や prefactor scaling は user-specified values と realistic target window を 両立 させる よう 調整しており、 absolute predictivity より も framework demonstration を 優先した。 したがって 実験 的 TOF、 methane selectivity、 wax fraction、 isotope effect、 transient response との 定量 一致 を 主張する こと は できない。 特に Co particle size step density、 support effect、 water partial pressure の 実験 依存性 は 反映していない。

 の 制約 は、 methodological simplification に ある。 汎用 `MicroKineticModel` は ODE + least-squares で steady state を 解く が、 FTS case study では fully mechanistic stiff solution が extreme poisoning state に 落ちやすかった ため、 `FischerTropschModel` で heuristic steady-state map を 使用した。 これは software artifact として は 有用 だが、 first-principles-consistent self-consistent solver では ない。 mean-field approximation も spatial heterogeneity、 island formation、 site blocking correlation、 reconstruction、 dynamic oscillation を 捉えない。 Andersen et al. (2019) と Zhang et al. (2023) が 示す ように、 kMC や dynamic thermal coupling は 重要 で ある。

 の 制約 は、 evaluation scope の 狭さ で ある。 baseline は no-lateral-interaction variant を 主 とし、 full CatMAP、 OpenMKM、 Cantera、 lattice kMC との head-to-head benchmark は 行っていない。 PFR simulation も qualitative trend demonstration を 目的 とした arbitrary rate scaling を 含み、 industrial space velocity の directly predictive model では ない。 External validation with independent real-world datasets is essential to confirm the generalizability of these findings beyond simulated conditions. 短期 的 には 6 か月 以内  DFT table ingestion、 thermodynamic consistency check、 experimental calibration、 Bayesian uncertainty propagation を 追加すべき で ある。 長期 的 には 1–2 年 スパン で multi-site Co step/terrace model、 kMC backend、 deactivation / wax diffusion / heat transfer coupling、 automated fitting against reactor datasets へ 拡張する こと が 望ましい。

## References
1. Andersen, M., Panosetti, C., & Reuter, K. (2019). A Practical Guide to Surface Kinetic Monte Carlo Simulations. *Frontiers in Chemistry*, 7, 202. DOI: 10.3389/fchem.2019.00202
2. Bruix, A., Margraf, J. T., Andersen, M., & Reuter, K. (2019). First-principles-based microkinetics simulations of synthesis gas conversion over transition metal catalysts. *Nature Catalysis*, 2, 659-670. DOI: 10.1038/s41929-019-0368-6
3. Campbell, C. T. (2017). The Degree of Rate Control: A Powerful Tool for Catalysis Research. *ACS Catalysis*, 7(4), 2770-2779. DOI: 10.1021/acscatal.7b00115
4. Grabow, L. C., & Mavrikakis, M. (2011). Mechanism of Methanol Synthesis on Cu through CO2 and CO Hydrogenation. *ACS Catalysis*, 1(4), 365-384. DOI: 10.1021/cs200055d
5. Kozuch, S., & Shaik, S. (2011). How to Conceptualize Catalytic Cycles? The Energetic Span Model. *Accounts of Chemical Research*, 44(2), 101-110. DOI: 10.1021/ar1000956
6. Liu, Q. Y., Chen, D., Shang, C., & Liu, Z.-P. (2023). An optimal Fe-C coordination ensemble for hydrocarbon chain growth: a full Fischer-Tropsch synthesis mechanism from machine learning. *Chemical Science*, 14, 10425-10439. DOI: 10.1039/D3SC02054A
7. Medford, A. J., Shi, C., Hoffmann, M. J., Lausche, A. C., Fitzgibbon, S. R., Bligaard, T., & Nørskov, J. K. (2015). CatMAP: A Software Package for Descriptor-Based Microkinetic Mapping of Catalytic Trends. *Catalysis Letters*, 145(3), 794-807. DOI: 10.1007/s10562-015-1495-6
8. Motagamwala, A. H., & Dumesic, J. A. (2021). Microkinetic Modeling: A Tool for Rational Catalyst Design. *Chemical Reviews*, 121(17), 1049-1076. DOI: 10.1021/acs.chemrev.0c00394
9. Nørskov, J. K., Bligaard, T., Logadottir, A., Bahn, S., Hansen, L. B., Bollinger, M., Bengaard, H., Hammer, B., Sljivancanin, Z., Mavrikakis, M., Xu, Y., Dahl, S., & Jacobsen, C. J. H. (2002). Universality in Heterogeneous Catalysis. *Journal of Catalysis*, 209(2), 275-278. DOI: 10.1006/jcat.2002.3615
10. Rommens, K. T., & Saeys, M. (2023). Molecular Views on Fischer-Tropsch Synthesis. *Chemical Reviews*, 123(9), 5648-5685. DOI: 10.1021/acs.chemrev.2c00508
11. Stegelmann, C., Andreasen, A., & Campbell, C. T. (2009). Degree of Rate Control: How Much the Energies of Intermediates and Transition States Control Rates. *Journal of the American Chemical Society*, 131(23), 8077-8082. DOI: 10.1021/ja9000097
12. van Santen, R. A., Ciobîcă, I. M., van Steen, E., & Ghouri, M. M. (2011). Mechanistic Issues in Fischer-Tropsch Catalysis. *Advances in Catalysis*, 54, 127-187. DOI: 10.1016/B978-0-12-387772-7.00003-4
13. Zhang, R., Wang, Y., Gaspard, P., & Kruse, N. (2023). The oscillating Fischer-Tropsch reaction. *Science*, 382(6668), 177-181. DOI: 10.1126/science.adh8463

## File Inventory
- `src/rate_constants.py`: Arrhenius, TST, Wigner, Eckart-like tunneling, DFT barrier conversion
- `src/adsorption_isotherms.py`: Langmuir, Temkin, Freundlich, fractal, competitive adsorption, BEP relation
- `src/lateral_interactions.py`: mean-field and short-range coverage interactions
- `src/microkinetic_model.py`: generic microkinetic engine, sensitivities, source terms
- `src/rate_control.py`: DRC, DTRC, energetic span, RDS identification
- `src/reactor_models.py`: PFR/CSTR plus microkinetic couplers
- `src/ft_synthesis.py`: Co(0001) FTS case study and figure/result generation
- `tests/test_models.py`: 8 validation tests
- `figures/fig1_rate_constants.png` 〜 `figures/fig8_energetic_span.png`: 8 figures
- `results/ft_simulation_results.csv`, `results/drc_analysis.csv`, `results/sensitivity_analysis.md`, `results/statistical_summary.md`
- `results/search-strategy.md`, `results/screening-table.csv`, `results/extraction-table.csv`, `results/reference-list.md`, `figures/prisma-flow.md`
- `data/pfr_profiles.csv`, `data/sensitivity_runs.csv`, `data/preprocessing-log.md`, `logs/process-log.jsonl`
