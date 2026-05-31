Now let me set up the Jupyter notebook and implement the QKD simulations in parallel:Now I have all the numerical results. Let me create the paper and report files:タスクが完了しました。以下が成果物のサマリーです：

---

## 完了した実験・成果物

### ステップ1: 先行研究調査（Crossref API 使用）
6件の論文を特定（2020年以降）：Su (2020) BB84解析、Mizutani et al. (2025) 有限鍵解析、Krawec (2023) Twin-field QKD、Yu (2025) 量子リピータ網、Ghosal et al. (2025) リピータプロトコル、Popp et al. (2025) 蒸留プロトコル

### ステップ2: NatureLM/GALACTICA
- **両ツールとも ToolUniverse レジストリに存在せず**（grep 結果 0 件）
- 代替として文献値・解析モデルで対応。paper.md の Methods セクションに記録済み

### ステップ3: Jupyter 実行結果（主要数値）

| 分析 | 主要結果 |
|-----|---------|
| BB84 有限鍵 [cell:1] | QBER=5%, N=10⁶: **0.2035** bits/signal（漸近値比95.3%） |
| E91 CHSH [cell:2] | QBER=1%: S=**2.7719**（Bell不等式違反確認） |
| 量子リピータ [cell:3] | 200km, n=4: F=**0.8833**, T_mem=**6.40 ms** |
| 蒸留 [cell:4] | F₀=0.90 → F=0.99: **7ラウンド, 195.6ペア** |
| 東京ネットワーク [cell:7] | 平均鍵レート**1316.9 kbps**、ボトルネック514.7 kbps |
| MC統計 [cell:9] | KS検定 p=**1.56×10⁻⁸³**（QBER 3% vs 5%） |

### 成果物ファイル
- **`paper.md`** — 学術論文形式（Abstract 200語以上、全セクション、参考文献9件）
- **`report.md`** — 日本語実験レポート（全数値・図表・考察）
- **7枚の図** (`figures/fig1_bb84_finite_key.png` ～ `fig7_comprehensive_summary.png`)
- **`data/raw/tokyo_qkd_simulation.csv`** — 東京QKDネットワーク生データ