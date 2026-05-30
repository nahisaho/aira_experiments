ToolUniverseのSemanticScholar等を使って先行研究を並行検索します。先行研究が揃いました。次にPythonシミュレーション実験を実行します。全7図が生成されました。次にpaper.mdとreport.mdを作成します。全ファイルが正常に生成されました。以下に実験の要約を示します。

---

## 実験完了サマリー

### ステップ1: 先行研究調査
- **SemanticScholar**: API 400エラー（全クエリ失敗）
- **Crossref / OpenAlex**: ✅ 成功 → **10件の主要論文**を同定（2020–2024年）
  - Yuan et al. 2021 (*Nature*, 引用1174件), Ngou et al. 2022 (*Plant Cell*, 引用932件)等

### ステップ2-3: 実験（6モジュール）

| モデル | 手法 | 主要結果 |
|---|---|---|
| 受容体結合 | ODE (RK45) | ETI信号=PTI信号の**7.1倍**（0.64 vs 4.53 a.u.） |
| MAPKカスケード | ODE (Michaelis-Menten) | Hill係数≈2.3の超感度スイッチ；ETI→~8 min半最大活性化 |
| SA/JAクロストーク | ODE (8変数) | 活物寄生性: PR1 max 47.03；壊死栄養性: PDF1.2 max 3.52 |
| WRKY/TGA GRN | Boolean Network (20ノード) | PTI/ETI両方が固定点アトラクター（PR1=ON, HR=ON）に収束 |
| 共進化ゲーム理論 | 複製子動力学 + Red Queen ODE | NLRアレイが最安定戦略；Red Queen周期~65世代 |
| イネいもち病 | ODE (12変数) + 交差検証n=5 | 抵抗性菌体負荷: **0.0059±0.0008** vs 感受性: **0.0676±0.0028**（11.5倍低減） |

### ステップ4: 成果物
- **`paper.md`** — 英語学術論文（Abstract 260語以上、References 12件、図7枚埋め込み）
- **`report.md`** — 日本語実験レポート（全セクション・図7枚埋め込み）
- **`figures/`** — 図1〜7（PNG、計7ファイル）
- **`src/plant_immunity_simulation.py`** — 全シミュレーションコード