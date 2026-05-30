Now let me find specific papers with DOIs for proper citation.Good, I now have solid references. Let me proceed to implement the full pipeline and generate figures.Now let me create the main analysis pipeline and generate all figures in one comprehensive script.Pipeline executed successfully. Now let me create the report.md and paper.md files.Now let me verify all files are in place and commit.全ステップ完了しました。以下が成果物のサマリーです：

## 先行研究調査（8件、DOI付き）
- Dominissini+ (2012, Nature) — MeRIP-seq原著
- Meyer & Jaffrey (2014, Nat Rev MCB) — m6Aレビュー
- Meyer (2019, Nat Methods) — DART-seq
- Huang+ (2020, Signal Transduct Target Ther) — m6Aと腫瘍
- Pratanwanich+ (2021, Nat Biotechnol) — xPore/nanopore
- Leger+ (2021, Nat Commun) — Nanocompore
- Jiang+ (2021, Mol Cancer) — がんm6Aレギュレーター
- Zhou+ (2026, GPB) — exomePeak2

## 実験結果
| 項目 | 結果 |
|------|------|
| MeRIP-seq感度/精度/F1 | 0.994 / 0.510 / 0.674 |
| Nanopore ML AUC | 1.000 |
| 差分修飾遺伝子 | 22 (hyper:10, hypo:12) |
| がん遺伝子Δm6A | +0.535 (p=1.2×10⁻⁹⁴) |

## 生成ファイル
- `src/epitranscriptome_pipeline.py` — 統合パイプライン
- `figures/fig1-8_*.png` — 8枚の図表
- `report.md` — 実験レポート
- `paper.md` — 学術論文（References含む）