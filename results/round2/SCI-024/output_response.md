I have sufficient literature data. Now let me implement the computational framework and generate figures.All figures generated. Now let me create the academic paper and report.すべてのファイルが正常に生成されました。以下に結果をまとめます。

---

## 完了報告

### ✅ 作成ファイル

| ファイル | 内容 |
|---------|------|
| `paper.md` | 英語学術論文（9セクション、参考文献9件、表4件） |
| `report.md` | 実験全結果レポート（日本語、計算コード付き） |
| `figures/fig1_band_structure.png` | バンド構造（トポロジカル vs. 自明相） |
| `figures/fig2_z2_wilson_loop.png` | Wilson Loop / Z₂不変量（WCC追跡） |
| `figures/fig3_surface_states.png` | 表面ディラック錐 + スピンテクスチャ |
| `figures/fig4_phase_diagram.png` | SOC–M₀位相図 |
| `figures/fig5_screening.png` | Bi₂Se₃類似体スクリーニング（12物質） |
| `figures/fig6_workflow.png` | QE/Wannier90/Z2Pack統合ワークフロー |

---

### ������ 主要結果

| 項目 | 値 |
|------|-----|
| バンドギャップ（Bi₂Se₃モデル） | **305 meV**（実験300 meV±1.7%） |
| 表面ディラック速度 vF | **3.33 eV·Å**（実験3.0–3.6 eV·Å） |
| Z₂不変量（バンド反転基準） | **ν=1**（M₀<0）/ **ν=0**（M₀>0） |
| スクリーニング結果 | 11/12物質がZ₂=1（トポロジカル） |
| **最有望新規候補** | **ZrHgSeBr₂**（NatureLM予測、ギャップ420 meV） |

### 🔬 NatureLM MCP 結果
- `predict_material_composition`: Bi-Se型とZrHgSeBr₂型を予測 ✅
- `ask_naturelm`: バンド反転エネルギー0.16 eV、臨界SOC 0.25 eV取得 ✅  
- `predict_property (band_gap/SOC)`: **非対応エラー** → `ask_naturelm`で代替

### ������ 先行研究（ToolUniverse MCP、8件特定）
ElcoroらのMagnetic TQC (Nature Commun. 2021)、Pengらの磁性空間群分類 (PRB 2022)、Zhangらの歪み位相転移 Li₂AuBi (Nano Lett. 2024) など計9件のDOI付き文献を収録。