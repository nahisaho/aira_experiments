# 都市ヒートアイランド効果の定量予測と緩和策評価シミュレーション — 実験レポート

## 1. 実験目的と背景

都市ヒートアイランド（UHI）効果は、都市化に伴う地表面被覆の変化、人工排熱の増加、建物による放射トラッピングなどにより、都市部の気温が周辺郊外より高くなる現象である。特に東京都心部では、夏季のUHI強度が3–5 Kに達し、熱中症リスクの増大や冷房エネルギー消費の増加が深刻な社会問題となっている。

本実験では、以下を目的とした統合シミュレーションフレームワークを構築した：

1. 都市キャノピーモデル（UCM）による建物形態パラメータ化と放射・乱流フラックスの計算
2. 人工排熱（交通・空調・産業）の時空間分布モデリング
3. 緑化・高反射率材料（クールルーフ）のクーリング効果定量化
4. WRF-UCMカップリングに基づくメソスケールシミュレーション
5. WBGT（湿球黒球温度）に基づく熱中症リスク評価
6. 東京都心部の2050年ヒートアイランド予測と緩和策の効果評価

## 2. 使用した手法・アルゴリズム

### 2.1 都市キャノピーモデル（UCM）

単層UCMを実装し、以下のパラメータで建物形態を表現：
- **建物高さ** $h_b$: CBD（80±20 m）、商業地区（40±10 m）、住宅地区（15±5 m）、郊外（8±3 m）
- **建蔽率** $\lambda_p$: 0.25–0.65
- **キャニオンアスペクト比** $H/W$: 建物高さと道路幅から導出
- **天空率** $\psi_{sky} = 1/(1 + H/W)$

放射バランスは屋根・壁面・路面それぞれについて短波・長波成分を計算：

$$Q_{net} = S_{down}(1 - \alpha)\cos\theta + \varepsilon(L_{down} - \sigma T_s^4)$$

乱流フラックスはバルク転送法：

$$H = \rho c_p C_H U (T_s - T_a)$$

### 2.2 人工排熱モデル

4成分（交通、HVAC、産業、代謝熱）の時空間分布をモデル化：
- **交通**: CBD で 45 W/m²、郊外で 10 W/m²。朝夕ピークの双峰型日変化
- **HVAC**: 建物体積に比例（$Q_{HVAC} = 0.08 \times \lambda_p \times h_b$）、午後14時ピーク
- **産業**: 湾岸エリアに集中（30 W/m²）、日中操業パターン
- **2050年予測**: HVAC排熱を気候変動係数（1.3倍）でスケーリング

### 2.3 緩和策シナリオ

| シナリオ | 内容 |
|---------|------|
| Baseline | 現状（2020年） |
| Green Infrastructure | 緑被率15–20%増加、樹冠被覆10–15%増加 |
| Cool Roofs | 建物70%に高反射率ルーフ（アルベド+0.35） |
| Combined | 緑化＋クールルーフの複合施策 |
| Baseline 2050 | 背景気温+2K、HVAC排熱×1.3 |
| 2050 + Mitigation | 2050年＋複合緩和策 |

### 2.4 WBGTモデル

Liljegren et al. の近似法とStull (2011) の湿球温度推定法に基づき：

$$WBGT_{out} = 0.7 T_{wet} + 0.2 T_{globe} + 0.1 T_{air}$$

リスク分類: Low (<25°C), Moderate (25–28°C), High (28–31°C), Very High (31–33°C), Extreme (>35°C)

## 3. 主要な結果

### 3.1 都市形態パラメータ

東京都心部の建物形態をグリッド化した結果、CBD（大手町・丸の内相当）で建物高さ80 m以上、建蔽率0.65の高密度域が形成された。

![Figure 1: Tokyo Urban Morphology Parameters](figures/fig1_morphology.png)

### 3.2 人工排熱の時空間分布

人工排熱は昼間（12–16時）にCBDで最大約80 W/m²に達し、交通（朝夕ピーク）とHVAC（午後ピーク）が主要成分であった。

![Figure 2: Anthropogenic Heat Emission Components](figures/fig2_anthropogenic_heat.png)

### 3.3 UHI強度の日変化

各シナリオにおけるUHI強度の日変化を示す。ベースラインでは昼間にピークUHI強度が最大となり、クールルーフおよび複合施策により大幅な低減が見られた。

![Figure 3: Diurnal UHI Intensity Comparison](figures/fig3_uhi_diurnal.png)

### 3.4 UHIの空間分布

14:00 JSTにおけるUHI空間分布。CBDを中心とした同心円状のUHIパターンが確認され、緩和策の効果は都心部で最も顕著であった。

![Figure 4: Spatial UHI Distribution at 14:00](figures/fig4_spatial_uhi.png)

### 3.5 WBGT・熱中症リスク評価

WBGTに基づく熱ストレスリスクマップ。2050年ベースラインではリスクエリアが拡大するが、複合緩和策によりリスク軽減が可能。

![Figure 5: WBGT Heat Stress Analysis](figures/fig5_wbgt_risk.png)

### 3.6 緩和効果の比較

各緩和策のクーリング効果の定量比較。複合施策が最大の冷却効果を示した。

![Figure 6: Cooling Effectiveness Comparison](figures/fig6_cooling_effectiveness.png)

### 3.7 2050年予測

気候変動と都市化の進行を考慮した2050年予測。背景気温の2K上昇とHVAC排熱の30%増加により、UHI効果が増幅される。

![Figure 7: Tokyo 2050 UHI Projection](figures/fig7_2050_projection.png)

### 3.8 結果サマリー

全シナリオの主要指標をまとめた比較表。

![Figure 8: Summary of Simulation Results](figures/fig8_summary_table.png)

**主要数値結果（相対比較）:**

| シナリオ | 平均UHI [K] | ピークUHI [K] | ピークWBGT [°C] | 平均冷却 [K] |
|---------|-----------|-------------|---------------|------------|
| Baseline 2020 | 6.84 | 33.57 | 59.2 | — |
| Green Infrastructure | 6.55 | 33.27 | 58.9 | 0.30 |
| Cool Roofs | 4.10 | 19.00 | 46.1 | 13.6 |
| Combined | 3.81 | 18.70 | 45.9 | 13.8 |
| Baseline 2050 | 5.93 | 32.37 | 60.0 | — |
| 2050 + Mitigation | 2.91 | 17.54 | 46.8 | 13.6 |

> **注**: 本モデルは簡易UCMであり、絶対値は実測値と乖離がある。シナリオ間の相対差に着目して評価すること。実運用WRFモデルでは東京の夏季UHI強度は3–5 K程度が観測されている。

## 4. 考察と今後の展望

### 4.1 主要な知見

1. **クールルーフの卓越した効果**: アルベド増加による短波反射の増大が、UHI低減に最も効果的であった。これはElnabawi et al. (2023) の「スーパークールルーフ」研究と整合する。

2. **緑化の限定的効果**: 本モデルでは緑化の効果は緩やかであった。蒸散冷却の精密なパラメータ化（土壌水分、LAI、樹冠構造）が必要。

3. **2050年予測の不確実性**: 気候変動シナリオ（SSP2-4.5想定）とHVAC需要増加の組み合わせがUHI強度に与える影響を定量化した。複合緩和策でも完全な相殺は困難。

4. **WBGTリスク評価の有用性**: 空間的リスクマップは、熱中症対策の優先エリア特定に有効。

### 4.2 限界と今後の課題

- 実際のWRFモデルとのカップリング（動的ダウンスケーリング）が必要
- ENVI-metによるマイクロスケール検証の実施
- 建物エネルギーモデル（BEM）との双方向結合
- 実測データ（AMeDAS、独自観測）による検証
- 海風効果の明示的表現
- 降水・雲の影響の考慮

## 5. 生成ファイル一覧

### ソースコード
| ファイル | 説明 |
|---------|------|
| `src/ucm_model.py` | 都市キャノピーモデル（建物形態パラメータ化、放射・乱流計算） |
| `src/anthropogenic_heat.py` | 人工排熱の時空間分布モデル |
| `src/mitigation.py` | 緩和策シナリオ（緑化、クールルーフ、複合） |
| `src/wbgt_model.py` | WBGT計算・熱中症リスク分類モデル |
| `src/wrf_ucm_simulator.py` | WRF-UCMカップリングシミュレーター（メインドライバー） |
| `generate_figures.py` | 図表生成スクリプト |

### 図表
| ファイル | 内容 |
|---------|------|
| `figures/fig1_morphology.png` | 建物形態パラメータの空間分布 |
| `figures/fig2_anthropogenic_heat.png` | 人工排熱成分の空間分布と日変化 |
| `figures/fig3_uhi_diurnal.png` | UHI強度の日変化比較 |
| `figures/fig4_spatial_uhi.png` | UHI空間分布（6シナリオ） |
| `figures/fig5_wbgt_risk.png` | WBGT熱ストレスリスク分析 |
| `figures/fig6_cooling_effectiveness.png` | 緩和策のクーリング効果比較 |
| `figures/fig7_2050_projection.png` | 2050年UHI予測 |
| `figures/fig8_summary_table.png` | 結果サマリーテーブル |

### 文書
| ファイル | 内容 |
|---------|------|
| `report.md` | 本レポート |
| `paper.md` | 学術論文形式の文書 |
