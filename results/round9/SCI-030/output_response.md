Now let me try NatureLM and GALACTICA tools, and set up the Jupyter notebook:Now let me create the paper.md and report.md files:すべてのステップが完了しました。以下に結果をまとめます。

---

## 完了サマリー

### ステップ1: 先行研究調査
Webサーチ経由で関連論文7件を特定（Semantic Scholar APIはレート制限429で一時不可）：
| 論文 | DOI |
|-----|-----|
| Lei et al. 2023 (THM-EGS) | 10.1038/s41598-023-41745-7 |
| Gudala et al. 2023 (SC-CO₂ vs Water) | 10.1115/1.4055538 |
| Xiao & Li 2024 (DFN-EGS) | 10.1002/ese3.1800 |
| Zhang et al. 2024 (THM coupling) | 10.1093/jge/gxae041 |
| Seismic Kakkonda 2024 | 10.3124/segj.77.24 |
| Mori et al. 2023 (Kakkonda monitoring) | 10.1038/s41598-023-35159-8 |

### ステップ2: NatureLM/GALACTICA
両ツールともToolUniverseに未登録（接続失敗）→ IAPWS-IF97（標準EOS）で代替。Methodsに記録済み。

### ステップ3: Jupyter実行結果（主要数値）
| 指標 | 値 | セル |
|-----|---|------|
| 平均熱出力 | **71.6 MW** | [cell:7c] |
| 平均電気出力（η=20%） | **14.3 MW** | [cell:7c] |
| 30年累積熱エネルギー | **18,816 GWh** | [cell:7c] |
| 熱的减衰 | **59.7 °C** | [cell:7c] |
| 最大ΔCS (r=100m) | **4.19 MPa** | [cell:9] |
| 最適坑井間隔 | **800 m** | [cell:6] |
| DFN P₂₁密度 | **0.0044 m/m²** | [cell:3] |

### 生成ファイル
- ������ **`paper.md`** — 学術論文（英語、Abstract 250語超、全セクション）
- 📄 **`report.md`** — 実験レポート（日本語）
- 🖼 `figures/fig01_dfn_map.png` — DFN亀裂ネットワーク
- 🖼 `figures/fig02_thm_results.png` — THM 30年結果（6パネル）
- 🖼 `figures/fig03_seismicity_risk.png` — 誘発地震リスク解析
- 🖼 `figures/fig04_overview.png` — 統合フレームワーク概要