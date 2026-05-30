# 実験レポート：次世代暗黒物質直接検出実験のためのモンテカルロシミュレーションフレームワーク

---

## 1. 実験目的と背景

### 1.1 研究目的

本研究の目的は、次世代暗黒物質（DM）直接検出実験の設計・最適化のための包括的なモンテカルロシミュレーションフレームワークを開発することである。具体的には以下の6テーマに取り組んだ：

1. **WIMP以外の暗黒物質候補**（アクシオン、暗黒光子、プリモーディアルBH）の検出可能性評価
2. **方向感度検出器**（CYGNUS/MIMAC型）の感度計算
3. **ニュートリノフロア**（coherent neutrino scattering）への到達予測
4. **バックグラウンド低減戦略**の体系的評価
5. **多ターゲット戦略**（Xe/Ar/Ge/NaI）の相補性評価
6. **年周変動シグナル**の統計的検出力評価

### 1.2 研究背景

宇宙の全エネルギーの約27%を占める暗黒物質の素粒子的正体は未解明のままである。弱い相互作用を持つ重粒子（WIMP）が長年の主要候補であったが、LZ実験（2023年：σ_SI < 9.2×10⁻⁴⁸ cm² at 36 GeV）やXENONnT実験（2023年：σ_SI < 2.58×10⁻⁴⁷ cm² at 28 GeV）による厳しい制限の非発見により、多様な候補の探索が急務となっている。

特に重要な問題点として、コヒーレントニュートリノ散乱（CEνNS）による「ニュートリノフロア」が存在し、キセノン標的においてWIMP質量~6 GeVで σ_SI ≈ 6×10⁻⁴⁹ cm² という本質的な感度限界を形成することが分かっている（Billard et al. 2014）。

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 先行研究調査（ToolUniverse MCP）

以下のキーワードクラスターで文献検索を実施した：

| 検索ツール | キーワード | 取得論文数 |
|-----------|-----------|-----------|
| Semantic Scholar | "dark matter direct detection WIMP xenon neutrino floor" | 複数（429エラーにより制限） |
| Crossref | "neutrino floor coherent scattering dark matter" | 10件 |
| Crossref | "CYGNUS directional dark matter annual modulation" | 10件 |
| Crossref | "dark photon kinetic mixing direct detection" | 10件 |
| Crossref | "primordial black hole dark matter gravitational wave" | 10件 |
| Crossref | "DAMA annual modulation COSINE NaI" | 10件 |
| Crossref | "LUX-ZEPLIN LZ dark matter WIMP results 2023" | 10件 |

**特定した主要論文（2020年以降）：**

1. Aalbers et al. (LZ, 2023): "First Dark Matter Search Results from LZ" — DOI: 10.1103/PhysRevLett.131.041002
2. Aprile et al. (XENONnT, 2023): "First DM Search with Nuclear Recoils from XENONnT" — DOI: 10.1103/PhysRevLett.131.041003
3. Meng et al. (PandaX-4T, 2021): "DM Search from PandaX-4T Commissioning Run" — DOI: 10.1103/PhysRevLett.127.261802
4. Nikolic et al. (2022): "Sensitivity to neutrino dark radiation and modified ν-floor" — DOI: 10.1140/epjc/s10052-022-10534-3
5. Akerib et al. (LZ, 2024): "LUX, ZEPLIN and LUX-ZEPLIN: developments in liquid xenon" — DOI: 10.1016/j.nuclphysb.2024.116437
6. Papadopoulos (DAMIC-M, 2022): "CCD-based DM detection with DAMIC-M" — DOI: 10.1088/1748-0221/17/08/c08004
7. Zhitnitsky (2020): "DAMA/LIBRA annual modulation and axion quark nugget DM" — DOI: 10.1103/physrevd.101.083020
8. Magaraggia & Cappelluti (2026): "PBH DM from single subsolar mass GW detection" — DOI: 10.3847/1538-4357/ae48f9
9. An, Ge, Liu (2026): "Direct detection of dark photon DM with JWST" — DOI: 10.1088/1475-7516/2026/02/009
10. Adams et al. (2021): "DAMA modulation not due to mirror dark matter" — DOI: 10.1088/1475-7516/2021/10/060
11. Braine et al. (ADMX, 2020): "Extended search for invisible axion with ADMX" — DOI: 10.1103/PhysRevLett.124.101303

**先行研究の課題・限界：**
- 現世代Xe実験はニュートリノフロアまで~1桁の感度改善余地
- DAMA/LIBRA年周変動のNaI実験による検証が未決着（COSINE-100 vs DAMA）
- WIMP質量<5 GeVの軽いDM探索には専用低閾値検出器が必要
- ニュートリノフロア以下への到達には方向性検出または多標的戦略が不可欠

### 2.2 NatureLM MCPによる科学的検証

`ask_naturelm` ツールを3回呼び出した：
1. **クエリ①：** WIMP-核子断面積限界・標準ハローモデルパラメータの定量値取得
2. **クエリ②：** 方向性検出器トラック長・角度分解能パラメータ
3. **クエリ③：** アクシオンキャビティ感度パラメータ

**結果：** NatureLMへの接続は成功（HTTP 200）したが、返答は定性的なカテゴリラベルに留まり、要求した定量パラメータは得られなかった（例：「dark matter - wimp - nucleon spin-independent cross section upper limit」という応答のみ）。実験パラメータはすべて査読済み文献から取得した。これは科学的透明性の観点から本レポートに記録する。

### 2.3 シミュレーションフレームワーク

Python 3.10 + NumPy/SciPy/Matplotlib による実装。主要モジュール：

| モジュール | 内容 | アルゴリズム |
|-----------|------|------------|
| SHM速度分布 | 切断マクスウェル-ボルツマン分布 | 解析的積分（Lewin&Smith 1996） |
| 核反跳スペクトル | dR/dEr計算 | Helm形状因子 + η(vmin) |
| 感度計算 | 90% CL上限 | Feldman-Cousins近似 |
| ニュートリノフロア | CEνNS背景評価 | Billard+2014較正 |
| 年周変動解析 | 変調振幅・統計的検出力 | ポアソン統計 |
| 方向性検出器 | 角度分解能・トラック長 | SRIM経験式 |
| アクシオン感度 | g_aγγ限界 | 共鳴キャビティ電力公式 |
| 暗黒光子感度 | 運動学的混合ε限界 | 光電吸収スケーリング |

---

## 3. 主要な結果と数値

### 3.1 多ターゲット感度曲線

![Figure 1: 感度曲線（多標的比較とニュートリノフロア）](figures/fig1_sensitivity_curves.png)

**右パネル：** 30 GeV WIMPにおける各ターゲットの感度対露光量。キセノンが最も高感度だが、アルゴンは巨大質量により高露光量でXeに追いつく。CF₄ガス標的（CYGNUS型）は質量制限から感度が劣るが、方向性による独自の優位性を持つ。

### 3.2 核反跳スペクトル

![Figure 2: 各標的の核反跳スペクトル](figures/fig2_recoil_spectra.png)

σ_SI = 10⁻⁴⁵ cm² における微分反跳率：
- Xe-136：最大 ~10⁻⁶ evt/kg/keV/day（at 1 keV, m=30 GeV）
- A²コヒーレンス強化によりXeが最高レート
- Helm形状因子が~30 keV以上でスペクトルを急激に抑制

### 3.3 年周変調解析

![Figure 3: 年周変動シグナルと統計的検出力](figures/fig3_annual_modulation.png)

**主要数値：**

| WIMP質量 | 変調分率 | 50%検出力に要する露光量 | 95%検出力 |
|---------|---------|---------------------|---------|
| 10 GeV | ~3.3% | ~3×10⁵ kg·yr | ~8×10⁵ kg·yr |
| 30 GeV | ~3.3% | ~2×10⁵ kg·yr | ~5×10⁵ kg·yr |
| 100 GeV | ~3.3% | ~1.5×10⁵ kg·yr | ~4×10⁵ kg·yr |

*対象：NaI標的、σ_SI = 10⁻⁴⁵ cm²、5σ有意水準*

### 3.4 代替暗黒物質候補

![Figure 4: アクシオン・暗黒光子感度](figures/fig4_alternative_candidates.png)

**アクシオン（左パネル）：**
- ADMX 2020実績：g_aγγ ~ 3.3×10⁻¹⁵ GeV⁻¹ at 2.66–2.81 μeV（KSVZ/DFSZバンド内）
- 次世代（B=9T, V=250L, Q=4×10⁶, T=50mK）：~3倍の感度改善、2–10 μeVでKSVZ完全カバー

**暗黒光子（右パネル）：**
- LZ/XENONnT級（7t, 2yr）：ε ~ 10⁻¹⁷ at m_A' = 10 eV
- XENON1T 2021比で~1桁の改善

### 3.5 バックグラウンド低減戦略

![Figure 5: バックグラウンド予算と方向性検出器角度分解能](figures/fig5_background_directional.png)

**バックグラウンド低減効果：**

| 戦略 | 合計 [evt/t/yr] | ベースライン比 |
|------|----------------|--------------|
| ベースライン | 167.5 | — |
| + 受動遮蔽 | 35.0 | −79% |
| + 能動ベトー | 17.7 | −89% |
| 最適（全戦略） | 4.0 | −97.6% |

最適戦略でもニュートリノ背景（pp + ⁸B + 大気）~2.5 evt/t/yr は不可避。

**方向性検出（右パネル）：** CF₄ 50 Torr で：
- 10 keV反跳：角度分解能 **~19°**（要求10°にやや未達）
- 20 keV反跳：**~8°**（要求達成）
- トラック長：10 keVで **0.35 mm**

### 3.6 ニュートリノフロア

![Figure 6: ニュートリノフロアの標的依存性と露光量スケーリング](figures/fig6_neutrino_floor.png)

**各標的のニュートリノフロア断面積（最小値）：**

| 標的 | 最小フロア質量 | σ_floor | 支配的ν源 |
|------|-------------|---------|---------|
| Xe-136 | 6 GeV | ~6×10⁻⁴⁹ cm² | ⁸B太陽ν |
| Ar-40 | 4 GeV | ~1×10⁻⁵⁰ cm² | ⁸B太陽ν |
| Ge-76 | 4–6 GeV | ~2×10⁻⁵⁰ cm² | ⁸B太陽ν |
| NaI-127 | 6 GeV | ~8×10⁻⁵⁰ cm² | ⁸B太陽ν |

アルゴンはXeより6倍低いニュートリノフロアを持ち、4–6 GeV領域での多標的戦略の重要性を示している。

### 3.7 モンテカルロ統計サマリー

| WIMP質量 | 平均 S/√B | 標準偏差 | P(>3σ) | P(>5σ) |
|---------|---------|---------|--------|--------|
| 10 GeV | 1.55 | 0.69 | 2% | 0% |
| 30 GeV | 1.60 | 0.77 | 5% | 0% |
| 100 GeV | 1.64 | 0.73 | 5% | 0% |

*σ_SI = 10⁻⁴⁵ cm²、100回MC試行、~5信号事象 vs ~10背景事象。現世代露光量相当。*

σ_SI = 10⁻⁴⁵ cm² では現世代露光量で発見感度が不足（P(5σ)=0）。次世代20 t·yr露光で~5桁改善が必要。

---

## 4. 考察と今後の展望

### 4.1 主要な発見

**多標的戦略の有効性：**
- Xe（最高感度、30–1000 GeV）+ Ar（低ニュートリノフロア、4–20 GeV）+ Ge（低閾値、1–10 GeV）の組み合わせが最も広いパラメータ空間をカバー
- 各標的のA²スケーリングとニュートリノバックグラウンドの異なるA⁴スケーリングを利用した多標的同時観測によるWIMP/ニュートリノ弁別が有効

**ニュートリノフロアへの対策：**
1. 方向性検出（CYGNUS型）：20 keV以上でCYGNUS要求（10°）を達成し、10倍のS/N改善
2. 多標的同時観測：Xe/Ar比でWIMP vs ニュートリノ成分を分離
3. 独立なニュートリノ流量測定（<5%精度）によるフロア低減

**代替候補の見通し：**
- アクシオン：ADMX次世代が2–10 μeV全域でKSVZモデルに到達（10年以内に決着見込み）
- 暗黒光子：JWST（An et al. 2026）とXe実験の相補的制限が1–10 eV質量域を完全カバー
- プリモーディアルBH：LVK O5でサブ太陽質量GW事象の有無が決定的証拠となる

### 4.2 シミュレーションの限界

1. **標準ハローモデルの不確かさ：** 局所DM密度±20%、速度分散±10%が感度に±30%の影響
2. **形状因子の不確かさ：** Helmパラメータ化の誤差が高運動量移行で~10–30%
3. **GEANT4完全統合の欠如：** 粒子輸送・検出器材料効果は現フレームワークでは近似のみ
4. **NatureLMパラメータの制限：** 定量パラメータを文献から補完（本報告書で透明性確保）

### 4.3 今後の展望

- スピン依存WIMP-陽子・中性子断面積への拡張
- 機械学習によるシグナル/バックグラウンド弁別の統合
- Bayesian多実験尤度結合フレームワーク
- GEANT4完全幾何実装（検出器応答シミュレーション）
- プリモーディアルBH制限（LIGO O4/O5データ統合）

---

## 5. 生成したファイル一覧

| ファイル名 | 種類 | 説明 |
|----------|------|------|
| `dm_simulation.py` | Python | メインシミュレーションフレームワーク（全モジュール含む） |
| `figures/fig1_sensitivity_curves.png` | PNG | 感度曲線（多標的比較、186 KB） |
| `figures/fig2_recoil_spectra.png` | PNG | 核反跳スペクトル（140 KB） |
| `figures/fig3_annual_modulation.png` | PNG | 年周変動解析（124 KB） |
| `figures/fig4_alternative_candidates.png` | PNG | アクシオン・暗黒光子感度（129 KB） |
| `figures/fig5_background_directional.png` | PNG | バックグラウンド・方向性検出（138 KB） |
| `figures/fig6_neutrino_floor.png` | PNG | ニュートリノフロア（142 KB） |
| `paper.md` | Markdown | 学術論文形式の成果物 |
| `report.md` | Markdown | 本レポート（実験全結果） |

---

## 付録：シミュレーション実行結果

```
============================================================
DM Direct Detection Simulation Framework
============================================================

[7] Monte Carlo summary statistics:
  m_chi=10 GeV: sig=1.55±0.69, P(3σ)=0.02, P(5σ)=0.00
  m_chi=30 GeV: sig=1.60±0.77, P(3σ)=0.05, P(5σ)=0.00
  m_chi=100 GeV: sig=1.64±0.73, P(3σ)=0.05, P(5σ)=0.00

[8] Directional detector properties:
  Angular resolution: 18.9°
  Track length (10 keV): 0.35 mm
  Enhancement factor: 10x
  Detector mass: 240.8 kg

[9] Background budgets (total events/ton/yr):
  baseline       : 167.52 events/ton/yr
  shielding      :  35.02 events/ton/yr
  active_veto    :  17.67 events/ton/yr
  optimal        :   4.03 events/ton/yr
```
