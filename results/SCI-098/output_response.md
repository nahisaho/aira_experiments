大規模なシミュレーションフレームワークを構築します。複数のモジュールを並行して作成していきます。次にメインのシミュレーション実行スクリプトと可視化を作成します。シミュレーションを実行します。シミュレーションは完了しましたが、物理計算の単位変換にいくつかの問題があります。修正してから図表を生成します。Now fix the annual modulation rate calculation and other issues:Now fix the sensitivity floor, neutrino floor calculation, directional detector, and annual modulation scaling:Now fix the neutrino floor calculation and the directional detector:Now fix the directional detector to handle no-H/T case with quadrupole test:Fix the simulation runner to handle inf values and cap sensitivity limits:Also fix the multi-target sensitivity capping and the annual modulation exposure:Results are improving. Let me debug the multi-target response (showing 0 events) and then generate figures:Found the root cause — the momentum transfer formula is wrong. `q_mev` calculation has incorrect unit conversion:Also fix the same bug in the background model's CEνNS cross section:Rates are now physically correct (~10⁻⁵ events/keV/kg/day at σ=10⁻⁴⁶). Re-running the full simulation:Excellent! Results are now physically correct. Let me generate the figures:Now create the process log and report:最終的なファイル構成を確認します：すべて正常に生成されました。成果物をまとめます：

---

## 完了サマリー

暗黒物質直接検出の次世代戦略シミュレーションフレームワーク **DMDDSF v1.0** を構築しました。

### 主要な結果

| 項目 | 結果 |
|------|------|
| **WIMP感度** | DARWIN: σ_min = 2.7×10⁻⁴⁸ cm² (m_χ=33 GeV) |
| **アクシオン** | g_ae = 10⁻¹³で~104事象（発見可能境界） |
| **方向性検出** | H/T有り15°で3σ検出に202事象 |
| **ニュートリノフロア** | Xe: ~3.5×10⁻⁴⁷ cm²（究極exposure） |
| **BG低減** | 最適組合せでGe 73.5%、Xe 69.1%削減 |
| **多ターゲット** | Xe(DARWIN)が12,749事象で圧倒的（m_χ=50 GeV） |
| **年周変動** | DAMA-like信号は1年で28σ確認可能 |

### 生成ファイル
- **ソースコード**: 7モジュール（`src/`配下）
- **結果データ**: 8 JSON（`results/`）
- **図表**: 7×2（PNG+SVG、`figures/`）
- **レポート**: `report.md`（全手法・結果・考察）