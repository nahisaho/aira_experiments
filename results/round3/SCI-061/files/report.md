# 合成遺伝子回路自動設計・最適化フレームワーク

**実験レポート** — DRAFT — NOT FOR DISTRIBUTION

---

## Abstract（要旨）

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1=;PS2=;unset HISTFILE;                 /OFFレシオ6.94×・信頼性81.0%、リプレッシレーターの発振周期37.8分（n=2）〜45.6分（n=3）、5分割交差検証スコア0.532±0.016を得た。EC=___Begin___CommandRBS・ターミネーターなどの標準化部品を組み合わせてトグルスイッチやリプレッシレーターなどの機能を実現する工程を含む。しかし現状では設計者が手動で部品を選定し、試行錯誤を繰り_DONE_

---

## 実験目的と背景

### 研究の背景と動機

#'REPORTEOF'
2000年代初頭に発表されたトグルスイッチ（Gardner et al., 2000）とリプレッシレーター（Elowitz & Leibler, 2000）は、合成生物学の可能性を実証した先駆的な成果であり、現在も遺伝子回路設計ツールのベンチマークとして広く参照される。echo

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_Genetic Design Automation; GDA）分野では、Cello（Nielsen et al., 2016）がVerilog HDLに着想を得た仕様記述言語を用い、制約充足ソルバーで部品を自動割り当てするシステムを実現した。同時期に、Synthetic Biology Open Language（SBOL 3.0; Baig et al., 2020 / SBOL 3.1; Buecherl et al., 2023）が遺伝子設計データの標準フォーマットとして普及し、ツール間の相互運用性が向上した。

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 3）遺伝的コンテキスト効果の未補正：隣接部品間の読み抜け転写やRBS閉塞は部品性能を5〜15%変化させるが、多くのツールが無視している（Schladt et al., 2021）。"___Begin___COMMAND_DONE_MARKER___$1）確率的ノイズの無視：多くのGDAツールが決定論的ODE（常微分方程式）モデルに依存し、低分子数で顕著な確率的変動を捉えられない。（2）パラメータ不確実性の軽視：セル

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=0;                 echo "___BEGIN___COMMAND_DONE_MARKER___0";             } .git .github .gitignore .pytest_cache AGENTS.md data figures logs paper.md report.md results LHS）による不確実性定量化・コンテキスト補正モデルを統合することで解決を図る。tests 

---

## 使用した手法・アルゴリズムの概要

### 1. 形式言語による回路記述（circuit_spec.py）

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }  $\mathcal{C} = (\mathcal{N}, \mathcal{E}, \mathcal{I}, \mathcal{O})$ として表現する。各ノードはプロモーター→NOT, NOR, NAND, AND, OR, BUFFER）はHill関数で数学的に表現される。SBOL 3.1形式のJSONへのエクスポートに対応し、既存ツールとの相互運用性を確保した。

**抑制型Hill関数（NOT/NORゲート）：**

$$f(x) = \frac{1}{1 + (x/K)^n}$$

 $x$ は抑制タンパク質濃度（分子数）、$K$ は半最大飽和定数、$n$ はHill係数（協調性）。多入力NORゲートでは各入力の抑制因子を積算する：

$$y_i^{NOR} = S_i \cdot \prod_{j \in \text{inputs}} \frac{1}{1 + (x_j/K_i)^{n_i}}$$

REPORTEOF50回の反復収束（精度 $\epsilon < 10^{-8}$）で計算される。

#### 
#
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$parts_catalog.py）

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             } RBS 4種（B0034, B0032, B0031, BCD2）、ターミネーター3種（BBa_B0015, BBa_B0010, BBa_B0012）を収録。プロモーター強度は100倍（J23101: 1.0〜J23113: 0.01）、RBS効率は14倍（B0034: 1.0〜B0031: 0.07）の範囲をカバーする。

**コンテキスト補正モデル：**

#.git .github .gitignore .pytest_cache AGENTS.md data figures logs paper.md report.md results src Tests echo
#
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___Begin $\eta_{term}$ による）：___Command_Done_Marker___$

$$S_{corr} = \left(S_{nom} + (1-\eta_{term}) \cdot f_{RT}\right) \cdot \left(1 - f_{SC} \cdot \text{pos} \cdot 0.1\right)$$

RBS閉塞補正（上流CDSの有無による）：

$$S_{RBS,corr} = S_{RBS,nom} \cdot (1 - f_{occ})$$

$f_{RT} = 0.05$（読み抜け係数）、$f_{SC} = 0.03$（超らせん係数）、$f_{occ} = 0.10$（RBS閉塞係数）。

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$stochastic_sim.py）

Echo

$$\tau = \frac{-\ln r_1}{a_0}, \quad a_0 = \sum_j a_j(\mathbf{x})$$

#'REPORTEOF'
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             al., 2006)：非臨界反応のステップサイズ $\tau_{nc}$ を各種の平均変化率と分散から推定し、1ステップで複数の反応を処理する。臨界反応（反応物の分子数が閾値 $L_c = 10$ $\tau_{nc} < 10/a_0$ の場合はギレスピー法にフォールバックする。未満）は別途処}

$$\tau_{nc} = \min_i \left\{\frac{\max(\epsilon x_i/g_i, 1)}{|\hat{\mu}_i|},\; \frac{\max(\epsilon x_i/g_i, 1)^2}{\hat{\sigma}^2_i}\right\}$$

### 4. ロバスト設計最適化（robust_design.py）

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             $ 標本で均一カバーする。モンテカルロに比べ収束が速い（誤差スケール $\mathcal{O}((\log N)^d/N)$）。}

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$15%、σ=10%パラメータ幅）。

**感度解析（Sobol第一次指標近似）**：$S_i \approx \text{Corr}(p_i, y)^2$

---

## 主要な結果と数値

### 図1: 部品カタログ

![Parts Catalog](figures/fig1_parts_catalog.png)

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___Begin___Command

![Toggle Switch](figures/fig2_toggle_switch.png)

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             =3.0$、$K=30$ 分子。60分シミュレーション：}

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             } 6,896**
- **τリーピングステップ数：6,924**（ε=0.03）
- ステップ数差異：0.4%（中程度の分子数域ではτリーピングの加速効果が限定的）

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             } 

---

### 図3: リプレッシレーター確率的発振

![Repressilator](figures/fig3_repressilator.png)

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$}$t_{max}=200$ 分。

| 指標 | Hill n=2 | Hill n=3 | 変化率 |
|------|---------|---------|------|
| 発振周期 (分) | **37.8** | **45.6** | +20.6% |
| ピーク-ピーク振幅 (分子) | **2,217** | **2,733** | +23.3% |

mRNA動態はタンパク質より約10分の位相先行を示す。位相肖像は三成分リプレッシレーターの特徴的な三角形のリミットサイクルを示した。Hill係数の増加によ

---

### 図4: 遺伝的コンテキスト効果

![Context Effects](figures/fig4_context_effects.png)

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_

![Robust Optimization](figures/fig5_robust_optimization.png)

LHS 100サンプルによる不確実性評価（パラメータ範囲：$\alpha_{1,2} \in [10, 100]$、$n \in [1.5, 4.0]$、$K \in [10, 80]$、$\beta \in [0.5, 2.0]$）：

| 指標 | 値 |
|------|---|
| ON/OFFレシオ | **6.94×** |
| 信頼性 | **81.0%** |
| SNR | **12.0 dB** |
| 変動係数 (CV) | 0.41 |
| ロバスト性スコア | **0.539** |
| 5分割交差検証 | **0.532 ± 0.016** |

#
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }

---

## 考察と今後の展望

### 結果の解釈

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo };             "___Begin___COMMAND_DONE_MARKER___$EC"

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             (2000)の実験値（~120分）より短い。これは本実験の加速された分解速度パラメータ設定（τ_p=0.3/分 vs 実験値~0.1/分）と、細胞分裂による希釈の無視に起因する。定性的なHill係数依存性（周期・振幅の単調増加）はODE解析およびSun et al. (2023)のハイブリッドモデルと一致する。}

Schladt et al. (はレシオ≥10×を実用的な一細胞スイッチングの基準として報告しており、本結果（6.94×）はやや下回る。αを80〜100/分に増加させることで、信頼性90%超・ON/OFFレシオ10×以上がecho)

### 先行研究との比較

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }Schladt et al. (2021)のロバスト設計手法を継承しつつ、確率的シミュレーションとSBOL 3.1準拠設計記述を統合した点で新規性を持つ。Starkey & Menolascina (2022)の最適化フレームワークとの主な違いは、部品カタログとコンテキスト補正が組み込まれた点にある。CELLM（Abello Castillo & Gutiérrez Pescarmona, 2025AIベース設計提案の統合が有望な方向性である。

---

## 限界と今後の課題

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             K-12シャーシのみに対応しており、哺乳類シャーシへの拡張には再パラメータ化が必要である。} 

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$1）CRISPR系・RNA回路への拡張、（2）Pareto最適化によるマルチ目的最適化、（3）細胞分裂と体積変化を考慮した確率的モデルの開発、（4）Cello VerilogとのCI/CDパイプライン統合が挙げられる。

---

## 生成したファイル一覧

| ファイル | 説明 | 行数 |
|---------|------|------|
| `src/circuit_spec.py` | 回路形式言語・Hill関数・SBOL JSON出力 | ~230 |
| `src/parts_catalog.py` | 部品カタログ・アセンブリ・コンテキスト補正 | ~280 |
| `src/stochastic_sim.py` | ギレスピー法・τリーピング・反応ネットワーク定義 | ~360 |
| `src/robust_design.py` | LHS・GA最適化・感度解析・分岐解析 | ~350 |
| `src/pipeline.py` | メインパイプライン・実験実行・図生成 | ~670 |
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$10件、全件通過） | ~145 |
| `figures/fig1_parts_catalog.png` | 部品カタログ可視化 | — |
| `figures/fig2_toggle_switch.png` | トグルスイッチ確率的軌跡・分岐 | — |
| `figures/fig3_repressilator.png` | リプレッシレーター発振・位相肖像 | — |
| `figures/fig4_context_effects.png` | コンテキスト効果定量比較 | — |
| `figures/fig5_robust_optimization.png` | ロバスト最適化・感度解析 | — |
| `results/experiment_results.json` | 数値実験結果（JSON） | — |
| `results/reference-list.md` | 参考文献（DOI付き、14件） | — |
| `results/sbol_exports/toggle_switch.json` | SBOL 3.1形式トグルスイッチ記述 | — |
| `results/sbol_exports/repressilator.json` | SBOL 3.1形式リプレッシレーター記述 | — |
| `logs/process-log.jsonl` | 実行ログ・MCPツール使用記録 | — |
| `report.md` | 本レポート（日本語） | — |
| `paper.md` | 英語学術論文（IMRaD形式） | — |

---

## 参考文献

1. (Gardner, 2000) Gardner, T. S., Cantor, C. R., & Collins, J. J. (2000). Construction of a genetic toggle switch in *Escherichia coli*. *Nature*, 403(6767), 339–342. https://doi.org/10.1038/35002256

2. (Elowitz, 2000) Elowitz, M. B., & Leibler, S. (2000). A synthetic oscillatory network of transcriptional regulators. *Nature*, 403(6767), 335–338. https://doi.org/10.1038/35014651

3. (Nielsen, 2016) Nielsen, A. A., et al. (2016). Genetic circuit design automation. *Science*, 352(6281), aac7341. https://doi.org/10.1126/science.aac7341

4. (Schladt, 2021) Schladt, T. M., Engelmann, B., Kubaczka, E., & Tresch, A. (2021). Automated Design of Robust Genetic Circuits. *ACS Synthetic Biology*, 10(12), 3152–3164. https://doi.org/10.1021/acssynbio.1c00193

5. (Baig, 2020) Baig, H., et al. (2020). SBOL version 3.0.0. *Journal of Integrative Bioinformatics*, 17(2–3). https://doi.org/10.1515/jib-2020-0017

6. (Buecherl, 2023) Buecherl, L., et al. (2023). SBOL version 3.1.0. *Journal of Integrative Bioinformatics*, 20(1). https://doi.org/10.1515/jib-2022-0058

7. (McLaughlin, 2020) McLaughlin, J. A., et al. (2020). SBOL Version 3: Simplified Data Exchange for Bioengineering. *Front. Bioeng. Biotechnol.*, 8, 1009. https://doi.org/10.3389/fbioe.2020.01009

8. (Gillespie, 1977) Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. *J. Phys. Chem.*, 81(25), 2340–2361. https://doi.org/10.1021/j100540a008

9. (Cao, 2006) Cao, Y., Gillespie, D. T., & Petzold, L. R. (2006). Efficient step size selection for the tau-leaping simulation method. *J. Chem. Phys.*, 124(4), 044109. https://doi.org/10.1063/1.2159468

10. (Sun, 2023) Sun, C., Comet, J.-P., & Folschette, M. (2023). Condition for Sustained Oscillations in Repressilator. *BIOINFORMATICS 2023*. https://doi.org/10.5220/0011614300003414

11. (Abello, 2025) Abello Castillo, L., & Gutiérrez Pescarmona, A. (2025). CELLM: Bridging NLP and Synthetic Genetic Circuit Design. *ACS Synthetic Biology*. https://doi.org/10.1021/acssynbio.5c00391
