Now I have sufficient literature. Let me build the full pipeline.Two small fixes needed:All 10 tests pass. Now run the full pipeline:Pipeline complete. Now let me read results and write the reports:Now I have all data. Let me write the reports:完了しました。以下が実験の全成果です。

---

## Experiment Complete: CO₂RR Computational Screening Pipeline

### Key Scientific Findings

1. **線形スケーリング関係**: TM+Cu合金系で *COOH scaling: $a = 0.80$, $R^2 = 0.944$（RMSE = 0.107 eV）。SAC MN4-C系では$a = 0.54$と浅くなり、N4配位が*COOHと*COの結合を部分的に分離することを示す。

2. **火山型プロット**: CO₂→CO 理論ピークはΔG*CO ≈ −0.49 eV（Cu付近、U_L = E° = −0.106 V）。Cu合金（CuZn, ΔG*CO ≈ −0.28 eV）はC2+生成の最適窓（ΔG*CO ≈ −0.3〜−0.5 eV）に位置し、Zhang et al. 2023の実験と整合。

3. **SAC選択性**: Fe-N4（ΔG*CO = −0.39 eV）とFeCo-N4（ΔG*CO = −0.19 eV）がCO2→CO選択的活性を示す。PZC補正は+0.009〜+0.026 eVのシフトをもたらす（U = −0.8 V）。

4. **溶媒効果**: *COOH*の暗黙的溶媒化補正（−0.18 eV）と電場補正（−0.12 eV）がU_Lを0.3〜0.4 V改善する。

5. **d-バンドモデル**: SAC全体のRMSEは〜0.3〜0.8 eVと中程度。Zn-N4は例外的に浅いd-バンド（εd ≈ −7.5 eV）でモデルが破綻し、非d-バンド機構の重要性を示す。

### 主要図表

![Volcano Plot](figures/fig2_volcano_CO2_to_CO.png)

### 成果物一覧

| 種別 | ファイル | 内容 |
|------|---------|------|
| ソースコード | `src/*.py` | 6モジュール（~1,200行） |
| レポート | `report.md` | 1,013語（日本語） |
| 論文 | `paper.md` | 3,630語（英語、IMRaD） |
| 図表 | `figures/` | 8図×PNG+PDF = 16ファイル |
| 結果 | `results/` | 4 CSV（36触媒の全スクリーニング結果） |
| ログ | `logs/process-log.jsonl` | 実行トレース |
| テスト | `tests/test_pipeline.py` | 10テスト（全パス） |

### Limitations

- CHEモデルは熱力学的限界電位のみを計算し、動力学的活性化障壁（±0.3〜0.8 eV）を無視
- 被覆率依存効果・横方向相互作用は未考慮（C2+予測に±0.3 eV誤差の可能性）
- GGA-PBEはCO吸着エネルギーを過大評価する傾向（Au/Ag系で顕著）