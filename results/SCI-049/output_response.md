Now let me create the complete anomaly detection system. I'll build all components in parallel:Now let me create the main experiment runner and visualization:The key is `top_rules` not `decision_rules` in the serialized results. Let me fix:All experiments ran and all 6 figures generated. Now let me create the report:全実験が完了し、すべてのファイルを生成しました。

## 実行サマリー

**6モジュール構成の異常検知システム**を設計・実装・検証しました：

| モジュール | 主要結果 |
|---|---|
| **変化点検出** | PELTが4/4の変化点を100%検出 |
| **多変量外れ値** | Isolation Forest F1=0.672（Deep SVDD F1=0.506を上回る） |
| **物理制約** | η範囲制約で1.26%の非物理的イベントを捕捉 |
| **ドリフト検出** | ADWIN 136件、Page-Hinkley 3件検出、再訓練8回発火 |
| **説明可能性** | 根本原因: pT(23.2%), missing_ET(21.2%)が主因 |
| **ストリーミング** | 37,449 events/sec、4層階層アーキテクチャ設計 |

**生成ファイル**: ソースコード7本、図表6枚、結果JSON、データCSV 2本、report.md — 詳細は `report.md` をご参照ください。