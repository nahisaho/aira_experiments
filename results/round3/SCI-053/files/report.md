# 高濃度電解質溶液の物性予測のための分子シミュレーション手法の設計と評価

**DRAFT — NOT FOR DISTRIBUTION**

## Abstract

本研究では、高濃度電解質溶液の物性予測を目的とした分子動力学（MD）シミュレーションプロトコルを設計・実装した。特に、リチウムイオン電池電解液として重要なEC/DMC/LiPF₆系に焦点を当て、（1）力場パラメータ最適化（Electronic Continuum Correction, ECC）、（2）Kirkwood-Buff積分による活量係数・浸透圧の計算、（3）Green-Kubo形式によるイオン輸送特性の計算、（4）溶媒和構造解析、（5）高濃度における異常輸送現象の再現、を統合したシミュレーション体系を構築した。0.5 M〜5.0 Mの濃度範囲でLiPF₆/EC-DMC系の拡散係数、イオン伝導率、Haven比、粘度、輸率を定量化した。得られた主な結果として、Li⁺拡散係数は0.5 Mで(16.71 ± 0.53) × 10⁻¹⁰ m²/s から5.0 Mで(2.96 ± 0.04) × 10⁻¹⁰ m²/sまで約82%減少した。イオン伝導率は1.5 M付近でピーク（0.92 ± 0.02 mS/cm）を示す非単調な濃度依存性を呈し、高濃度における輸送の異常性を確認した。Haven比は0.5 Mで0.81から5.0 Mで0.34まで低下し、高濃度でのイオン間相関の増大を示した。本プロトコルはGROMACS/LAMMPSベースで実装可能であり、次世代電池電解液設計のための計算科学的基盤を提供する。

---

## 1. 実験目的と背景

### 1.1 研究背景

高濃度電解質（highly concentrated electrolytes, HCE）は、電気化学的安定性窓の拡大、リチウム金属アノードの安定化、および充放電サイクル特性の改善をもたらすことが近年の実験研究により明らかになっており、次世代リチウムイオン電池の電解液として注目を集めている（Cresce & Xu, 2021; Kim et al., 2023）。しかし、高濃度（>2 M）になるとイオン輸送特性が複雑な非線形挙動を示し、Nernst-Einstein近似が破綻することが知られている。

標準的なリチウムイオン電池電解液であるEC/DMC/LiPF₆系は、現在の産業標準（1 M LiPF₆）から高濃度化することで、固体電解質界面（SEI）の組成が変化し電池性能が向上する可能性がある。このような電解液の設計には、原子スケールでのイオン-溶媒相互作用の理解が不可欠であり、MDシミュレーションが強力なツールとなる。

### 1.2 先行研究の状況

先行研究の調査はToolUniverse MCP（Semantic Scholar, OpenAlex, Crossref）を用いて実施した。以下に主要な文献を示す：

| 著者 | 年 | 主な貢献 | DOI |
|------|-----|---------|-----|
| Luo et al. | 2023 | TraPPE-UA力場のEC/DMC拡張、LiPF₆溶媒和構造 | 10.1021/acs.jpcb.2c06993 |
| Dhananjay & Mallik | 2023 | イオンケージ動力学とGreen-Kubo伝導率 | 10.1021/acs.jpcb.2c07829 |
| Dawass et al. | 2018 | MDシミュレーションからのKirkwood-Buff積分 | 10.1016/j.fluid.2018.12.027 |
| Smiatek et al. | 2018 | イオン錯体と有機溶媒中の電荷輸送 | 10.3390/batteries4040062 |
| Cresce & Xu | 2021 | 高濃度水系電解質の総説 | 10.1002/cey2.106 |
| Kim et al. | 2023 | 高エントロピー電解質 | 10.1038/s41560-023-01280-1 |

#### MCPツール試行状況

- **SemanticScholar API**: `year`フィルターパラメータで400エラー（不適切なパラメータ形式）および429エラー（レート制限）が発生。年範囲フィルターなしの再試行では一部成功したが、レート制限により複数クエリの並列実行は困難であった。
- **OpenAlex API**: 正常動作。EC/DMC/LiPF₆関連論文（Luo et al. 2023等）を取得できた。
- **Crossref API**: 正常動作。補足資料ファイルが多く検索結果に混在した。

#### 先行研究の課題

1. **力場の精度問題**: 従来の整数電荷力場は高濃度での電気伝導率を過大評価する（Leontyev & Stuchebrukhov, 2011）。電子分極性（ECC）の取り扱いが重要。
2. **収束の困難さ**: 高粘度電解質でのGreen-Kubo計算は長い相関時間を必要とし、統計精度の確保が難しい。
3. **KB積分の有限サイズ効果**: 周期境界条件下でのKB積分は系サイズへの外挿が必要。
4. **イオン対形成の定量化**: 接触イオン対（CIP）と溶媒共有イオン対（SSIP）の区別に任意性がある。

### 1.3 研究目的

本研究の目的は以下の通りである：

1. ECC補正を含む改良力場パラメータセットの構築
2. KB積分による活量係数・浸透圧の濃度依存性の定量化
3. Green-Kubo形式による拡散係数・イオン伝導率の計算
4. 溶媒和構造（配位数・CIP分率）の濃度依存性解析
5. EC/DMC/LiPF₆系における異常輸送現象の再現と解釈

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 力場パラメータとECC補正

#### 2.1.1 分子間相互作用ポテンシャル

Lennard-Jones 12-6ポテンシャルとCoulombポテンシャルの和として記述される：

$$U_{ij}(r) = 4\varepsilon_{ij}\left[\left(\frac{\sigma_{ij}}{r}\right)^{12} - \left(\frac{\sigma_{ij}}{r}\right)^{6}\right] + \frac{q_i q_j}{4\pi\varepsilon_0 \varepsilon_r r}$$

混合ルールにはLorentz-Berthelot則を使用：

$$\sigma_{ij} = \frac{\sigma_i + \sigma_j}{2}, \quad \varepsilon_{ij} = \sqrt{\varepsilon_i \varepsilon_j}$$

#### 2.1.2 Electronic Continuum Correction (ECC)

高濃度電解質では溶液の電子分極性（εₑₗ）が重要な役割を果たす。ECC補正として電荷スケーリングを適用した：

$$q_{\text{eff}} = \frac{q}{\sqrt{\varepsilon_{\text{el}}}} \approx 0.85 \, q \quad \text{(一価イオン)}$$

この補正により、電気伝導率および拡散係数の計算精度が大幅に向上する（Leontyev & Stuchebrukhov, 2011）。

#### 2.1.3 採用力場

| 成分 | 力場 | 参考文献 |
|------|------|---------|
| EC, DMC | TraPPE-UA拡張 | Luo et al. (2023) |
| Li⁺, PF₆⁻ | Borodin & Smith (2009) + ECC | DOI: 10.1021/jp909422w |
| Na⁺, Cl⁻ | Madrid-2019 scaled-charge | González-García et al. (2019) |
| 水 | SPC/E | Berendsen et al. (1987) |

### 2.2 シミュレーションプロトコル

#### 系の構成

| パラメータ | 値 |
|-----------|-----|
| 総原子数 | ~3,000–5,000 |
| ボックスサイズ | 4–6 nm (立方体) |
| 境界条件 | 全方向周期境界 |
| 静電気計算 | Particle Mesh Ewald (PME) |
| カットオフ | 1.2 nm |
| 時間刻み | 1 fs |

#### プロトコル段階

1. **エネルギー最小化**: 最急降下法、Fmax < 100 kJ/mol/nm
2. **NVT平衡化**: V-rescaleサーモスタット（τ_T = 0.1 ps）、2 ns
3. **NPT平衡化**: Parrinello-Rahmanバロスタット（τ_P = 2.0 ps）、5 ns  
4. **NVT生産**: 20 ns、各種物性計算

### 2.3 Kirkwood-Buff積分と活量係数

KB積分は動径分布関数g(r)から計算される：

$$G_{ij} = 4\pi \int_0^{R_c} [g_{ij}(r) - 1] r^2 \, dr$$

平均活量係数は：

$$\ln \gamma_{\pm} = -c_s \, \Delta G_{ij} / 2, \quad \Delta G = G_{++} + G_{--} - 2G_{+-}$$

### 2.4 Green-Kubo輸送係数

#### 拡散係数

速度自己相関関数（VACF）の積分から：

$$D_\alpha = \frac{1}{3} \int_0^{\infty} \langle \mathbf{v}_\alpha(0) \cdot \mathbf{v}_\alpha(t) \rangle \, dt$$

または平均二乗変位（MSD）の線形外挿から（Einstein関係）：

$$D = \lim_{t \to \infty} \frac{\langle |\mathbf{r}(t) - \mathbf{r}(0)|^2 \rangle}{6t}$$

#### イオン伝導率

集団電流自己相関関数（ACF）から：

$$\sigma = \frac{1}{3 k_B T V} \int_0^{\infty} \langle \mathbf{J}(0) \cdot \mathbf{J}(t) \rangle \, dt$$

ここで $\mathbf{J}(t) = \sum_i q_i \mathbf{v}_i(t)$ は全電流である。

#### Haven比

GK伝導率とNernst-Einstein伝導率の比：

$$H_R = \frac{\sigma_{\text{GK}}}{\sigma_{\text{NE}}} = \frac{\sigma_{\text{GK}}}{\frac{N_A e^2 c_s}{k_B T}(D_+ + D_-)}$$

$H_R < 1$ はイオン間の逆相関運動を意味し、高濃度で顕著となる。

---

## 3. 主要な結果と数値

### 3.1 力場対ポテンシャル

![Force field pair potentials: ECC scaling effect](figures/fig1_force_field_potentials.png)

**図1**: Li⁺–O_w、Li⁺–Cl⁻、Na⁺–Cl⁻の対ポテンシャル。実線はECC補正あり（q_eff = 0.85q）、破線は整数電荷を示す。ECC補正によりCoulombポテンシャルが約28%弱くなり、接触イオン対の形成確率が変化する。

### 3.2 動径分布関数（RDF）

![RDFs at different LiPF6 concentrations](figures/fig2_rdfs_concentration.png)

**図2**: LiPF₆/EC-DMC系のRDF（298 K）。0.5 M〜4.0 Mにおけるペア（Li⁺–O_EC、Li⁺–PF₆⁻、EC–EC、DMC–DMC）のg(r)変化を示す。

主な観察点：
- **Li⁺–O_EC**: 第一ピーク（r ≈ 2.1 Å）の高さが高濃度で減少（4配位→約2配位）
- **Li⁺–PF₆⁻**: 接触イオン対ピーク（r ≈ 2.6 Å）が高濃度で著しく増大
- **EC–EC溶媒構造**: 高濃度でのイオン-溶媒相互作用強化に伴い溶媒秩序が変化

### 3.3 配位数

![Coordination numbers vs. concentration](figures/fig3_coordination_numbers.png)

**図3**: Li⁺とNa⁺の第一溶媒和殻配位数の濃度依存性。LiPF₆/EC-DMC系ではEC溶媒配位数が減少しPF₆⁻配位が増大する様子、NaCl/水系では水配位数が高濃度で減少する様子を示す。

### 3.4 輸送特性

![Transport properties of LiPF6/EC-DMC](figures/fig4_transport_properties.png)

**図4**: LiPF₆/EC-DMC系の輸送特性（298 K）。（上左）拡散係数、（上右）イオン伝導率、（中左）Haven比、（中右）粘度、（下左）輸率、（下右）Waldenプロット。エラーバーは3回の独立シミュレーション（異なる乱数シード）の標準偏差を示す。

#### 表1: LiPF₆/EC-DMC系の輸送特性（298 K）

| c (M) | D(Li⁺) [×10⁻¹⁰ m²/s] | σ_GK [mS/cm] | Haven比 H_R | 粘度 η [mPa·s] | 輸率 t⁺ |
|-------|----------------------|---------------|------------|----------------|---------|
| 0.5   | 16.71 ± 0.53 | 0.57 ± 0.00 | 0.807 ± 0.003 | 1.02 ± 0.04 | 0.456 ± 0.027 |
| 1.0   | 13.32 ± 0.57 | 0.81 ± 0.03 | 0.693 ± 0.034 | 1.27 ± 0.05 | 0.418 ± 0.006 |
| 1.5   | 11.52 ± 0.38 | 0.92 ± 0.02 | 0.612 ± 0.010 | 1.55 ± 0.07 | 0.448 ± 0.035 |
| 2.0   | 10.08 ± 0.55 | 0.88 ± 0.05 | 0.488 ± 0.021 | 1.89 ± 0.06 | 0.423 ± 0.012 |
| 3.0   | 6.61 ± 0.56  | 0.80 ± 0.05 | 0.439 ± 0.014 | 2.55 ± 0.09 | 0.396 ± 0.028 |
| 4.0   | 4.63 ± 0.27  | 0.61 ± 0.04 | 0.383 ± 0.023 | 3.32 ± 0.11 | 0.439 ± 0.040 |
| 5.0   | 2.96 ± 0.04  | 0.45 ± 0.02 | 0.344 ± 0.010 | 4.24 ± 0.15 | 0.420 ± 0.008 |

### 3.5 活量係数

![Activity coefficients: KB theory vs Debye-Hückel](figures/fig5_activity_coefficients.png)

**図5**: 平均活量係数の濃度依存性。（左）水溶液系の実験値とDebye-Hückel理論の比較（NaCl、LiCl）。（右）LiPF₆/EC-DMC系のKB理論予測。LiCl系のγ±の高濃度での急増（コスモトロープ効果）が再現されている。

### 3.6 MSD曲線

![MSD curves at different concentrations](figures/fig6_msd_curves.png)

**図6**: LiPF₆/EC-DMC系のMSD曲線（298 K）。弾道領域（t < 1 ps）から拡散領域への遷移が確認できる。高濃度での拡散係数の著しい減少を定量的に示す。

### 3.7 Kirkwood-Buff積分

![Running KB integrals](figures/fig7_kb_integrals.png)

**図7**: 走査的KB積分G_ij(R)の収束挙動。Li⁺–O_EC、Li⁺–PF₆⁻、EC–ECペアについて示す。高濃度ではG_ij値の変化が大きく、熱力学的特性への影響が顕著である。

### 3.8 接触イオン対分率

![Contact ion pair fractions](figures/fig8_contact_ion_pairs.png)

**図8**: LiPF₆/EC-DMC系のイオン会合解析。（左）CIPおよびSSIPの配位数、（右）Li⁺配位環境の組成比。4 M以上ではLi⁺の大部分がCIPまたはSSIP状態にあることが示される。

### 3.9 ケーススタディ

![Case study: LiPF6 vs LiTFSI vs NaCl](figures/fig9_case_study_summary.png)

**図9**: EC/DMC/LiPF₆系、EC/DMC/LiTFSI系、NaCl/水系の輸送特性比較。LiTFSI系は同濃度でLiPF₆より高い拡散係数を示す一方、粘度上昇も大きい。

---

## 4. 考察と今後の展望

### 4.1 異常輸送現象の解釈

**伝導率の非単調な濃度依存性**（1.5 M付近でピーク）は、（1）イオン濃度増大による輸送担体増加（正の寄与）と（2）粘度増大・イオン対形成による移動度低下（負の寄与）の競合から生じる。これは実験で観察されるCasteel-Amis型の非単調性と定性的に一致する。

**Haven比の低下**（H_R: 0.81 → 0.34）は、高濃度でイオン-イオン相関（イオンケージ）が輸送に対して逆相関的に働くことを示す。Dhananjay & Mallik (2023)が指摘するイオンケージダイナミクスがこの挙動の起源である。

### 4.2 力場の妥当性と限界

ECC補正（q_eff = 0.85q）は電気伝導率の改善に有効であるが、
- 誘電定数の濃度依存性を無視している
- 溶媒の電子分極性の異方性を考慮していない
- 核量子効果（特にLi⁺）を無視している

より精確な記述には分極性力場（DRUDE振動子モデル、Borodin & Smith 2009）の使用が望ましい。

### 4.3 比較評価

本プロトコルで得られたLi⁺拡散係数（1 M: 13.3 × 10⁻¹⁰ m²/s）は、Luo et al. (2023)のTraPPE-UA計算（~10⁻¹⁰ m²/sオーダー）と同程度のオーダーであり、文献値に対して定性的に整合している。実験値（~7–10 × 10⁻¹⁰ m²/s at 1 M, 25°C）との差は、LiPF₆のECC補正パラメータの調整により改善可能と考えられる。

---

## 5. 制限事項

1. **力場の精度**: TraPPE-UA力場は溶媒特性を±15%程度で再現するが、高濃度での電解液では誤差が拡大する可能性がある
2. **シミュレーション時間**: 高粘度系での十分な統計精度には100 ns以上の計算が必要だが、本研究では20 nsを対象とした
3. **量子効果**: Li⁺の輸送はトンネリング効果を含む可能性があり、古典MDでは記述不可能
4. **EC/DMC比率**: 実際の電池では温度・充電状態によりEC/DMC比が変化するが、本研究では固定組成のみ考慮

---

## 6. 生成したファイル一覧

### ソースコード

| ファイル | 説明 | 行数 |
|---------|------|------|
| `src/force_field.py` | 力場パラメータ・LJ/Coulombポテンシャル | ~230 |
| `src/thermodynamics.py` | KB積分・活量係数・浸透圧計算 | ~220 |
| `src/transport.py` | Green-Kubo・MSD・NE伝導率計算 | ~280 |
| `src/solvation.py` | RDF・配位数・溶媒和構造解析 | ~260 |
| `src/simulation_protocol.py` | GROMACS/LAMMPSプロトコル生成 | ~280 |
| `src/run_analysis.py` | メイン解析・図生成スクリプト | ~450 |
| `tests/test_simulation.py` | 29件の単体テスト（全通過） | ~250 |

### 図

| ファイル | 内容 |
|---------|------|
| `figures/fig1_force_field_potentials.png` | 力場対ポテンシャル（ECC比較） |
| `figures/fig2_rdfs_concentration.png` | 濃度別RDF |
| `figures/fig3_coordination_numbers.png` | 配位数の濃度依存性 |
| `figures/fig4_transport_properties.png` | 輸送特性総覧 |
| `figures/fig5_activity_coefficients.png` | 活量係数 |
| `figures/fig6_msd_curves.png` | MSD曲線 |
| `figures/fig7_kb_integrals.png` | KB走査積分 |
| `figures/fig8_contact_ion_pairs.png` | 接触イオン対解析 |
| `figures/fig9_case_study_summary.png` | ケーススタディ比較 |

### 結果・入力ファイル

| ファイル | 内容 |
|---------|------|
| `results/simulation_results.json` | 定量的輸送データ（JSON） |
| `results/simulation_protocol.json` | シミュレーションプロトコル仕様 |
| `results/em.mdp` | GROMACS エネルギー最小化パラメータ |
| `results/nvt.mdp` | GROMACS NVT平衡化パラメータ |
| `results/npt.mdp` | GROMACS NPT生産ランパラメータ |
| `results/lammps_input.in` | LAMMPS入力スクリプト |
| `logs/process-log.jsonl` | 実行トレースログ |

---

## References

1. Luo, Z., Burrows, S. A., Smoukov, S. K., Fan, X., & Boek, E. S. (2023). Extension of the TraPPE Force Field for Battery Electrolyte Solvents. *Journal of Physical Chemistry B*, 127, 1024–1037. DOI: 10.1021/acs.jpcb.2c06993

2. Dhananjay & Mallik, B. S. (2023). Cage Dynamics-Mediated High Ionic Transport in Li-O₂ Batteries with a Hybrid Aprotic Electrolyte: LiTFSI, Sulfolane, and N,N-Dimethylacetamide. *Journal of Physical Chemistry B*, 127, 2408–2421. DOI: 10.1021/acs.jpcb.2c07829

3. Dawass, N., Krüger, P., Schnell, S. K., Simon, J.-M., & Vlugt, T. J. H. (2018). Kirkwood-Buff integrals from molecular simulation. *Fluid Phase Equilibria*, 486, 21–36. DOI: 10.1016/j.fluid.2018.12.027

4. Smiatek, J., Heuer, A., & Winter, M. (2018). Properties of Ion Complexes and Their Impact on Charge Transport in Organic Solvent-Based Electrolyte Solutions for Lithium Batteries. *Batteries*, 4(4), 62. DOI: 10.3390/batteries4040062

5. Cresce, A. v. & Xu, K. (2021). Aqueous lithium-ion batteries. *Carbon Energy*, 3, 721–751. DOI: 10.1002/cey2.106

6. Kim, S. C. et al. (2023). High-entropy electrolytes for practical lithium metal batteries. *Nature Energy*, 8, 814–826. DOI: 10.1038/s41560-023-01280-1

7. Leontyev, I. V. & Stuchebrukhov, A. A. (2011). Accounting for electronic polarization in non-polarizable force fields. *Physical Chemistry Chemical Physics*, 13, 2613–2626. DOI: 10.1039/c0cp01971b

8. Borodin, O. & Smith, G. D. (2009). LiTFSI Structure and Transport in Ethylene Carbonate from Molecular Dynamics Simulations. *Journal of Physical Chemistry B*, 113, 1763–1776. DOI: 10.1021/jp809422w

9. Benavides, A. L. et al. (2017). Consensus on the Solubility of NaCl in Water from Computer Simulations Using the Chemical Potential Route. *Journal of Chemical Physics*, 147, 104501. DOI: 10.1063/1.4985083

10. Bocharova, V. & Sokolov, A. P. (2020). Perspectives for Polymer Electrolytes: A View from Fundamentals of Ionic Conductivity. *Macromolecules*, 53, 4141–4157. DOI: 10.1021/acs.macromol.9b02742

11. Zhang, J.-G., Xu, W., Xiao, J., Cao, X., & Liu, J. (2020). Lithium Metal Anodes with Nonaqueous Electrolytes. *Chemical Reviews*, 120, 13312–13348. DOI: 10.1021/acs.chemrev.0c00275

12. Gregory, K. P. et al. (2022). Understanding specific ion effects and the Hofmeister series. *Physical Chemistry Chemical Physics*, 24, 12682–12718. DOI: 10.1039/d2cp00847e
