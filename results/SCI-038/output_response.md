Now I have sufficient literature. Let me build the full simulation system.Now let me create the main simulation code:Need to fix two issues: the net simulation casting error and the RPO ΔV values. Let me fix both:Now fix the RPO initial conditions - the issue is that v0 should be reasonable:All simulations completed successfully. Now let me verify the figures and create the report and paper:Now let me verify the workspace and commit:すべて完了しました。以下が成果物の概要です：

## 成果物

**シミュレーション** (`adr_simulation.py`) — 6モジュール統合ADRミッション設計システム：
1. **ターゲット選定**: 30デブリから衝突リスク×除去効果で上位5件選定
2. **軌道最適化**: GA(200世代)で巡回順序最適化 → 総ΔV 4,331 m/s
3. **RPO**: Hill方程式で510mからのランデブー → ΔV 0.833 m/s
4. **タンブリング推定**: 平均誤差 0.50 deg/s (7.5%相対誤差)
5. **捕獲機構**: アーム(30s)/ネット(5.0s)/ハープーン(3.2s)全成功
6. **コスト最適化**: 全120順列探索 → $66.0M / 230日

**生成図**: 8枚（figures/配下に全て格納）

**文書**:
- `report.md` — 全結果・手法・考察の日本語レポート（図埋め込み済み）
- `paper.md` — 英語学術論文（12件のDOI付き参考文献、Abstract〜References完備）