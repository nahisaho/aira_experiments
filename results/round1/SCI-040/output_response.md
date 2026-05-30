Now I have sufficient literature. Let me proceed with the implementation.Good. Now let me create the main inversion framework and experiment scripts.The PyMC graph has issues with the forward model. Let me switch to a more robust approach using `scipy` and `emcee`-style MCMC.Good. Now let me rewrite using emcee for MCMC:All experiments completed successfully. Let me verify the figures and results:All experiments ran successfully. Now let me create the report.md and paper.md files.Now let me clean up the old file and verify everything:All files are in place. Here's a summary of what was completed:

## 完了サマリー

**ステップ1（先行研究調査）**: Web検索で8件以上の関連論文を特定（Liao et al. 2023, Nikkhoo & Rivalta 2022/2023, Anderson & Gu 2024, Sun & Tang 2021, Holt et al. 2022, Fukui et al. 2013 等）

**ステップ2（実験計画）**: 6つの実験を設計 — モデル比較、MCMC逆解析、統合インバージョン、カルマンフィルタ、粘弾性補正、ケーススタディ

**ステップ3（実験実施）**: 全実験を `volcano_inversion_v2.py` で実行完了

### 主要結果
- **MCMC**: 全パラメータを相対誤差<1%で回復（深さ: 4999±4m, 真値5000m）
- **統合インバージョン**: 残差がノイズレベルと同等（GNSS ~3.4mm, InSAR ~4.9mm）
- **カルマンフィルタ**: 体積変化RMSE = 2.02×10⁵ m³
- **粘弾性**: η=10¹⁸で67%、η=10¹⁹で41%の変位増幅（10年）
- **ケーススタディ**: 桜島・阿蘇とも高精度回復

### 生成ファイル
- `report.md`, `paper.md` — レポートと学術論文
- `volcano_inversion_v2.py` — フレームワーク本体
- `figures/` — 9枚の図（全てreport.md/paper.mdに埋め込み済み）
- `results.json` — 数値結果