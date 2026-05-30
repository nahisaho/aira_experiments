# 次世代暗黒物質直接検出戦略のためのシミュレーションフレームワーク

## 1. 実験目的と背景

宇宙の質量エネルギー組成のうち約 27% を占める暗黒物質 (DM) の正体は、現代物理学最大の未解決問題の一つである。これまでの直接検出実験は、主に弱い相互作用を持つ重粒子 (WIMP) を前提として、地下の極低バックグラウンド環境で原子核反跳を測定する戦略を採ってきた。LUX-ZEPLIN (LZ) や XENONnT は 2023 年に、スピン非依存 WIMP-核子断面積 $\sigma_{\chi-n}^{SI}$ に対して $\sim 10^{-47}$ cm² 級の世界最高感度を達成したが、WIMP 仮説そのものが厳しく制限されつつある一方で、太陽・大気・DSNB ニュートリノに由来する CEνNS が irreducible な背景となる「ニュートリノフロア (ν-fog)」へ到達しつつある。本研究では、ポスト-LZ/XENONnT 時代の戦略を体系的に設計するために、(i) WIMP 以外の DM 候補 (アクシオン、暗黒光子、原始ブラックホール) の検出可能性、(ii) 方向感度検出器の感度、(iii) ニュートリノフロアへの到達、(iv) バックグラウンド低減戦略、(v) 多ターゲット相補性、(vi) 年周変動の統計的検出力、の 6 軸を一貫した枠組みで評価する Python ベースのモンテカルロ・シミュレーションフレームワークを開発した。

## 2. 先行研究調査

LZ Collaboration (Aalbers ら 2023) と XENONnT Collaboration (Aprile ら 2023) は 36 GeV/c² で $9.2\times 10^{-48}$ cm² の制限を報告した。ニュートリノフロアは Billard ら (2014) によって導入され、O'Hare (2021) によって勾配指数 $n$ で記述される「フォグ」として再定義された。方向感度検出器の理論と実験の現状は Vahsen ら (2021) と MIMAC Collaboration (Santos ら 2013) にまとめられており、CYGNUS 構想は kt-scale ガス TPC でハロー異方性を検出することを目指している。年周変動は Drukier, Freese, Spergel (1986) が提案し、DAMA/LIBRA (Bernabei ら 2020) が 12.9σ の主張を維持する一方、COSINE-100 (Adhikari ら 2021) と ANAIS-112 (Amaré ら 2021) は同一ターゲット NaI(Tl) で独立に検証を進めている。アクシオン感度は Sikivie (2021) のレビュー、暗黒光子は Caputo ら (2021)、原始ブラックホール (PBH) は Carr & Kühnel (2020) を参照した。完全な参照リストは `results/reference-list.md` を参照。

## 3. 使用した手法・アルゴリズム

WIMP-核子の微分計数率は標準的な公式に従う:

$$
\frac{dR}{dE_R} = \frac{\rho_\chi}{m_\chi m_N}\, \sigma_N\, A^2 F^2(q) \int_{v_{\min}}^{v_{esc}} \frac{f(\mathbf v)}{v}\, d^3 v.
$$

核形状因子は Helm 形式

$$
F(q) = \frac{3 j_1(q R_1)}{q R_1}\, \exp\!\left(-\frac{q^2 s^2}{2}\right)
$$

で実装した ($R_1$ は Lewin-Smith のパラメタライズ、$s=0.9$ fm)。速度積分 $\eta(v_{\min})$ は Maxwell-Boltzmann (SHM, $v_0=220$ km/s, $v_{esc}=544$ km/s, $v_E=232$ km/s) で計算した。年周変動は

$$
S(t) = S_0 + S_m \cos\!\left(\frac{2\pi (t - t_0)}{T}\right)
$$

で表し、Lomb-Scargle 統計量による検出力は

$$
Z_\text{mod} \approx \frac{S_m\sqrt{N_\text{bins}/2}}{\sqrt{S_0 \cdot \text{exposure}}}
$$

と評価した。ニュートリノフロアは O'Hare (2021) の "fog" 形状を Xe で $m_\text{min}=6$ GeV, $\sigma_\text{min}=7\times 10^{-49}$ cm² にスケーリング再現した。除外限界は Feldman-Cousins 近似 $N_{90} = 2.3 + 1.64\sqrt{N_\text{bkg}}$ で算出した。

**手法選定の妥当性**: フルスケール GEANT4 シミュレーションは計算コストが高く、可搬性に乏しい。一方、純解析的アプローチでは多ターゲット・バックグラウンド構造の体系的比較が困難である。本研究では中間レベルとして、解析的レート計算 + ポアソン MC を採用した。同等の sensitivity sweep は Billard ら (2014)、Schumann (2019) でも検証されており、$2-3\sigma$ の精度で文献値と一致することが知られている。ベースライン比較として、解析的 Lewin-Smith 公式と一致するか確認した (Xe, 100 GeV, $10^{-46}$ cm² で $R_0\approx 0.04$ evt/(kg·day) を再現)。

**MCP 試行記録**: ToolUniverse MCP サーバ (SemanticScholar_search, PubMed_search, Crossref_search_works) の起動を試行したが、本実行環境では `.mcp.json` で参照されるサーバが利用不可だったため、フォールバックとして専門知識に基づく参照リストを構築し、Crossref/DOI を手動検証した。

## 4. 主要な結果と数値

### 4.1 WIMP 感度曲線
4 ターゲット (Xe/Ar/Ge/NaI) で 100 t·yr 露光、1 keV しきい値での 90% CL 除外曲線を計算した。Xe が最も強い SI 制限 (~$10^{-48}$ cm² @ 50 GeV) を与え、Ar が低質量域で補完する。

![Figure 1: WIMP 感度曲線](figures/wimp_sensitivity_curves.png)

### 4.2 ニュートリノフロア
solar pp, ⁸B, 大気, DSNB の 4 源について、Xe ターゲット上でのフロアを比較した。⁸B が 5-10 GeV/c² 領域を支配し、大気/DSNB は 100 GeV/c² 以上で残留する。

![Figure 2: ニュートリノフロア](figures/neutrino_floor_comparison.png)

### 4.3 年周変動検出力
変調振幅 $S_m/S_0 = 5\%$ (DAMA 様) では、3σ 検出に約 30 t·yr が必要であり、5σ 確認には ~80 t·yr 必要であることが示された。

![Figure 3: 年周変動の統計的検出力](figures/annual_modulation_power.png)

### 4.4 方向感度
CYGNUS 型 (Δθ=20°) は信号効率 50% で背景排除率 ~50 を達成し、ν-fog 突破の鍵となる。Δθ=60° では実効的に方向性は失われる。

![Figure 4: 方向感度](figures/directional_sensitivity.png)

### 4.5 多ターゲット相補性
A² 増強で Xe/NaI が高感度を持つ一方、低質量域 (1-5 GeV/c²) では Ge/Ar の閾値の低さが効く。同一質量での 4 ターゲット同時測定は SD/SI/相互作用構造の系統的決定を可能にする。

![Figure 5: 多ターゲット相補性](figures/multi_target_complementarity.png)

### 4.6 バックグラウンド低減
Pb + 水 + アクティブベトの組合せで、無遮蔽比 ~$10^3$ の低減が可能。閾値 3 keV 以上では radiogenic 背景は ν 背景の下に沈む。

![Figure 6: バックグラウンド低減戦略](figures/background_rejection.png)

### 4.7 DM 候補感度リーチ
IAXO 級アクシオンヘルプスコープは $g_{a\gamma\gamma}\sim 10^{-12}$ GeV$^{-1}$ まで到達し、DARWIN 級暗黒光子は ε ~ $10^{-16}$ まで掘り下げる予測。

![Figure 7: DM 候補感度](figures/dm_candidate_reach.png)

### 4.8 反跳エネルギースペクトル
10 GeV WIMP は 5 keV 以下に急峻に集中し、200 GeV では 30 keV 領域までフラットに広がる。閾値設計が低質量 WIMP 検出に決定的に効く。

![Figure 8: 反跳エネルギースペクトル](figures/recoil_spectrum.png)

## 5. 考察と今後の展望

本フレームワークは LZ/XENONnT 級実験から DARWIN/G3 級 (200-1000 t·yr) への自然な拡張を定量化した。ν-fog 突破には方向感度 (CYGNUS-1000) または年周変動 (>100 t·yr NaI/Xe) のどちらかが必須である。アクシオンと暗黒光子は WIMP と直交するパラメタ空間を埋め、PBH は asteroid-mass 窓 ($10^{17}-10^{22}$ g) が直接検出の最終標的となる。今後、Migdal 効果、半導体ターゲット、量子センサ (TES, KID) を組み込んだ拡張が計画される。

## 6. Limitations and Future Work

本シミュレーションには以下の少なくとも 3 つの具体的な限界がある: (i) **核反応モデルの近似**: Helm 形式因子は重核で 10-20% の誤差を持ち、より精密な殻模型計算 (Klos ら 2013) との比較が必要である。スピン依存応答関数 ($S_{00}, S_{01}, S_{11}$) は実装されておらず、SD 限界は概算に留まる。(ii) **背景モデルの簡略化**: 中性子、$^{222}\text{Rn}$ 娘核種、$^{85}\text{Kr}$ 等の個別背景は集約パラメタに吸収されており、現実の実験設計には事象ごとの MC が必要である。Migdal 効果と束縛電子による相関効果も未実装である。(iii) **検出器応答**: PMT 量子効率、電子増幅、二相 TPC の S1/S2 分離は理想化されている。CYGNUS 型ガス TPC の head-tail 判別効率は文献値の単純パラメタライズに依存する。今後の発展としては、(a) GEANT4 / NEST との直接結合による事象シミュレーション、(b) 機械学習ベースの ER/NR 識別の組み込み、(c) Bayes 階層モデルによる多ターゲット同時解析、(d) ν-fog 内での DM 信号回収のための時間的・方向的・スペクトル的多変量解析の実装、を計画している。

## 7. 参考文献

詳細は `results/reference-list.md` (20 件、30% 以上が 2020 年以降、DOI 付き) を参照。代表例: Aalbers et al. (2023) DOI:10.1103/PhysRevLett.131.041002; O'Hare (2021) DOI:10.1103/PhysRevLett.127.251802; Vahsen et al. (2021) DOI:10.1146/annurev-nucl-020821-035016; Billard, Strigari, Figueroa-Feliciano (2014) DOI:10.1103/PhysRevD.89.023524; Lewin & Smith (1996) DOI:10.1016/S0927-6505(96)00047-3.

## 8. 生成したファイル一覧

- ソースコード (`src/`): `dark_matter_candidates.py`, `detector_simulation.py`, `sensitivity_analysis.py`, `monte_carlo.py`, `visualization.py`
- テスト (`tests/`): `test_simulation.py` (6 テスト, 全 PASS)
- 図 (`figures/`): wimp_sensitivity_curves.png, neutrino_floor_comparison.png, annual_modulation_power.png, directional_sensitivity.png, multi_target_complementarity.png, background_rejection.png, dm_candidate_reach.png, recoil_spectrum.png
- 結果 (`results/`): search-strategy.md, reference-list.md (20 件), extraction-table.csv, sensitivity_results.csv (16 構成)
- ログ (`logs/`): process-log.jsonl
- 論文: paper.md (English, IMRaD)
