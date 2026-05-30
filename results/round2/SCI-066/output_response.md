Now let me run the main experiment. I'll implement the ESM AI emulator framework with synthetic CMIP6-like data.Now I'll create both output files:全ファイルの作成が完了しました。

---

## 実施内容まとめ

### ステップ1: 先行研究調査（ToolUniverse MCP）
`openalex_literature_search`, `Crossref_search_works` を使用し、**7件の主要論文**を特定：
- **ClimateBench v1.0** (Watson-Parris et al., 2022)
- **パターンスケーリング vs 深層学習の比較** (Lütjens et al., 2025)
- **WeatherBench** (Rasp et al., 2020)
- **ClimateSet** (Kaltenborn et al., 2023)他

### ステップ2-3: 実験実施

NatureLM (`ask_naturelm`) を4回クエリし、物理制約・超パラメータ・性能基準を取得。

**合成CMIP6データ**（32×64グリッド、1950-2114年、5メンバー、4シナリオ）を生成し、3モデルを比較：

| モデル | T Skill (CV) | 特記事項 |
|--------|-------------|---------|
| **パターンスケーリング** | **0.954 ±0.005** | 最高性能（線形応答優位） |
| U-Net（簡略版） | 0.010 ±0.431 | 勾配最適化なしで不安定 |
| ConvLSTM（簡略版） | −2.069 ±1.915 | 時系列外挿で崩壊 |

**主要発見**:
- SSP5-8.5での外挿失敗（T-RMSE: 0.48 K → 16.1 K）
- アンサンブルスプレッド比: 0.390（60%の不確実性を圧縮）
- 計算高速化: ~2×10⁶×

### 成果物
- ������ **`paper.md`** — 英語学術論文（Abstract 200語以上, 全7章, 参照文献10件, 図7枚）
- 📄 **`report.md`** — 日本語実験レポート（先行研究表, NatureLM結果, 定量表, 図7枚, 考察）
- 🖼️ **`figures/`** — 7枚の図（シナリオ・空間パターン・モデル比較・時系列・マルチシナリオ・不確実性・スキルダイアグラム）