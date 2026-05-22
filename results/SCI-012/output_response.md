Now let me implement all the components:Now let me create the main experiment runner and visualization script:The simulation ran but Potjans and WM need parameter tuning. Let me fix the external drive calculations:The Potjans model needs current calibration. Let me fix it with proper drive calculation and vectorize the inner loop:All experiments ran successfully. Now let me generate the final `report.md`:全実験が完了しました。以下が成果のサマリーです。

---

## 実験結果サマリー

### ������ ニューロンモデル比較
| モデル | 発火率 | ISI CV | 計算速度 |
|--------|--------|--------|---------|
| Hodgkin-Huxley | 70 Hz | 0.005 | ~100ms/s-sim |
| Izhikevich-RS | 25 Hz | 0.233 | ~2ms/s-sim (**50×高速**) |
| AdEx | 15 Hz | 0.344 | ~3ms/s-sim |

### ������ シナプス可塑性
- STDP + ホメオスタティックスケーリングが協調動作
- 最終平均荷重: 0.3004（双峰分布に収束）

### ⚡ 並列アーキテクチャ
- **100,000 ニューロン × 2,000 万シナプス** を 7 秒未満で実行（Numba JIT）
- CUDA バックエンド（自動選択）により GPU 環境では 10–50× 加速

### 🧬 Potjans-Diesmann モデル
- E/I 発火率比 ≈ 1:3.4（文献の定性的パターンを再現）
- MI(L2/3E→L4E) = **0.281 bits**, TE = **0.036 bits**

### 🎯 作業記憶タスク
- Assembly 0: 符号化 28.2 Hz → 遅延期 **6.9 Hz（持続発火）**
- 実験値(Funahashi 1989) との比較 z スコア: **−0.37**（統計的整合）

生成ファイル: `report.md`、9 図、6 結果 JSON、`logs/process-log.jsonl` — 全て `workspace/` に保存済み。