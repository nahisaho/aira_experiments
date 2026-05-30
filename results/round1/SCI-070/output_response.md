Now I have sufficient literature. Let me proceed with the experiment implementation.The optimization is slow. Let me optimize the code and re-run.All experiments completed and 7 figures generated. Now let me create the report and paper documents.すべての成果物が完成しました。以下に要約します：

## 完了した作業

### ステップ1: 先行研究調査
12件の文献を特定（DOI付き10件）。Costanza et al. (2014), Chaplin-Kramer et al. (2019), Hein et al. (2020), Bagstad et al. (2021), Villa et al. (2014) 等、InVEST/ARIES/SEEA-EA の主要論文をカバー。

### ステップ2: 実験計画
先行研究の手法を参考に、6モジュール統合パイプラインを設計。改良点として、選択実験によるWTP統合、6種割引率の体系的比較、SEEA-EA準拠会計を追加。

### ステップ3: 実験実施 — 主要結果
- **InVEST空間分析**: 50×50グリッド、5土地利用、5サービスモデル
- **WTP推定**: 生物多様性 ¥1,552、水質 ¥1,149、景観 ¥986/世帯/年
- **割引率**: 市場利子率(5%)とStern型(1.4%)でNPVに2.6倍の差
- **里山**: 基準ES価値 ¥726,200/ha/年、農地転換で71.3%損失

### 生成ファイル
| ファイル | 内容 |
|---------|------|
| `report.md` | 全結果レポート（図7枚埋め込み） |
| `paper.md` | 学術論文（References 12件、DOI付き） |
| `experiment.py` | 実験コード |
| `figures/` | 7枚のPNG図表 |