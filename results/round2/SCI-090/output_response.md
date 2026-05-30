Now I have enough literature data. Let me create the simulation framework and output files.All 6 figures generated. Now creating the output files.完了しました。以下が成果物の概要です。

---

## 完了報告

### 📊 シミュレーション結果サマリー

| 指標 | ベースライン | ZEB最適化 | 改善率 |
|---|---|---|---|
| 年間EUI | 294.5 kWh/m²/年 | 147.5 kWh/m²/年 | **−49.9%** |
| CO₂排出量 | 138.4 kg/m²/年 | 69.0 kg/m²/年 | **−50.1%** |
| 昼光自律率(DA) | 62.3% | 74.2% | +19.1% |
| 夏季換気(ACH) | 4.5 | 12.4 | +175.6% |
| IFC変換率 | — | **95.6%** | — |

### ������ NatureLM使用記録（透明性）
- 4クエリ実施。ZEB熱パラメータ・東京エネルギー負荷でCOP<1.0・消費量0桁違いなど**非現実的な数値**が返答されたため定量値として不採用。定性的傾向確認と昼光指標説明のみ参照。全試行をMethodsに記録。

### 📁 生成ファイル

| ファイル | 内容 |
|---|---|
| `paper.md` | 学術論文（英語、358行、参考文献10件）|
| `report.md` | 実験レポート（日本語、292行）|
| `figures/fig1〜fig6.png` | 6図（IFC変換・熱負荷・CFD・昼光・ダッシュボード・感度解析）|
| `bim_simulation.py` | Pythonシミュレーション本体 |

### 感度解析の主要発見
**窓U値が最大影響因子**（Pearson r=0.805）→ ZEB投資優先度の明確化