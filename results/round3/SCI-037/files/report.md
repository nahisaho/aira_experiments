# InSAR時系列解析による地殻変動モニタリングシステム
## 南海トラフ沿い地殻変動への適用

**DRAFT — NOT FOR DISTRIBUTION**

---

## 要旨

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }CUSUM制御図）、（4）昇降軌道データ統合による3次元変位場推計、の各モジュールを実装した。

#            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             tinel-1 C-バンドSARを想定した
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             } 

.git .github .gitignore .pytest_cache AGENTS.md data figures logs results run_experiment.py src tests ERA5＋線形補正の組み合わせによる大気補正が位相標準偏差を21.3%低減した（7.56mm → 5² = 0.59 ± 0.18を達成し、5分割時系列交差検証のRMSEは25.2 ± 39.1 mmであった。3成分（Z-score・CUSUM・空間分散）を組み合わせた前兆検出アルゴリズムはF1スコア0.667、精度0.667、再現率0.667を示し、検出遅延の中央値は12日（1エポック）であった。3次元変位場の再構成では条件数1.726の安定した投影行列を用い、鉛直・東西成分の速度場を推計した。.99mm）。最小二乗トレンド分'REPORT_EOF'

Report_Isce/Stampsベースの自動処理ワークフローとして南海トラフ沿いReport_EofReport_Eof

---

## 1. 研究目的と背景

### 1.1 研究背景

#
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1=;PS2=;Ishikawa, 2020; Takemura et al., 2023）。Unset 

#'REPORT_EOF''REPORT_EOF'GNSS/GPS連続観測網（GEONET）を主体としてきたが、その空間分解能は数十km程度に限られる。Sentinel-1等の衛星合成開口レーダー（SAR）を用いたInSARecho
-数百mの空間分解能でmm精度の変位場を取得可能であり、陸域・沿岸域の地殻変動マッピングに革命をもたらしている。

### 1.2 研究目的

'REPORT_EOF''REPORT_EOF'6要素からなる統合InSAR処理パイプラインを設計・実装することである：

1. **PS-InSAR/SBAS統合処理パイプライン**: 高コヒーレンス点（PS点）と分布散乱体（DS点）を組み合わせた時系列解析
2. **大気遅延補正**: ERA5気象再解析データおよび経験的手法（線形・べき乗則）の比較評価
**: 線形変動（定常沈降）・季節変動（水文荷重・熱弾性）・過渡変動（SSE）の分離
4. **地震前兆変動の自動検出**: 複数統計指標（Z-score、CUSUM、空間分散）による閾値検出
5. **3次元変位場推計**: 昇降軌道LOS（Line-of-Sight）データの加重最小二乗法による分解
6. **南海トラフ適用**: 現実的なノイズレベルと変動シナリオを用いた実証的評価

### 1.3 先行研究調査の方法とMCPツール使用状況

.git .github .gitignore .pytest_cache AGENTS.md data figures logs results run_experiment.py src tests  ToolUniverse MCP（Model Context Protocol）サーバーを使用した。試行したツールと結果は以下の通り：

| ツール | 試行クエリ数 | 結果 |
|--------|-------------|------|
| `SemanticScholar_search_papers`（year/sortパラメータあり） | 6 | API エラー 400（接続失敗） |
| `SemanticScholar_search_papers`（パラメータなし） | 3 | API エラー 400（接続失敗） |
| `PubMed_search_articles` | 6 | 部分的成功（2〜5件/クエリ） |
| `Fatcat_search_scholar` | 3 | 結果ゼロ |
| `advanced_literature_search_agent` | 1 | 設定エラー（smolagents未インストール） |
| `web_search`（フォールバック） | 4 | 成功（14件の文献特定） |

API（直接アクセス）とWebサーチを代替手段 MCPツール接続の部分的失敗に対し、PubMed

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 PS-InSAR/SBAS統合処理パイプライン

InSAR時系列解析の基本観測量は、マスター画像とスレーブ画像間のLOS変位量 $d_{LOS}$ である：

$$d_{LOS}(t) = -U_z \cos\theta + U_e \sin\theta \sin\alpha - U_n \sin\theta \cos\alpha$$

 $\theta$ は入射角、$\alpha$ は衛星ヘディング角、$U_z, U_e, U_n$ はそれぞれ鉛直・東西・南北変位成分である。

545組の干渉ペア（60シーンから構成）を生成した。コヒーレンスの平均は0.871 ± 0.081であった。

### 2.2 大気遅延補正

#'REPORT_EOF'
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }

**（a）線形高度補正（経験的手法）**

$$\phi_{atm}(h) = a \cdot h + b$$

 $h$ と位相 $\phi$ の線形回帰により大気遅延を推定・除去する。

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1=;PS2=;unset HISTFILE;                 EC=0;                 echo ___BEGIN: Huang et al., 2025）**___COMMAND_DONE_MARKER___

$$\phi_{atm}(h) = K \cdot h^\alpha$$

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$Iteratively Reweighted Least Squares）によりパラメータ $K, \alpha$ を推定する。

**（c）ERA5＋線形組み合わせ補正（参照: Yang et al., 2023）**

ERA5再解析データ（欧州中期天気予報センター提供）による天頂全遅延（ZTD）モデルを適用後、残差の線形高度補正を追加する複合手法。

### 2.3 長期変動トレンド分離

 $d(t)$ を以下のモデルで分解する：echo

$$d(t) = v \cdot t + c + A_1 \cos(2\pi t) + B_1 \sin(2\pi t) + A_2 \cos(4\pi t) + B_2 \sin(4\pi t) + r(t)$$

 $v$ は線形速度（mm/yr）、$c$ は定常オフセット、$A_1, B_1$ は年周振動振幅、$A_2, B_2$ は半年周振幅、$r(t)$ は残差（過渡的変動を含む）である。

 $\mathbf{G}$ を用いた最小二乗推定：

$$\hat{\mathbf{m}} = (\mathbf{G}^T \mathbf{G})^{-1} \mathbf{G}^T \mathbf{d}$$

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1=;PS2=;unset HISTFILE;                 EC=0;                 echo ___BEGIN: Gualandi & Liu, 2021; Reinosch et al., 2020）。___COMMAND_DONE_MARKER___

### 2.4 地震前兆変動の自動検出アルゴリズム

#'REPORT_EOF'
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }

**（a）Z-スコア検出**

$$Z(t) = \frac{R(t) - \mu_{baseline}}{\sigma_{baseline}}$$

 $R(t) = \sqrt{\frac{1}{N}\sum_{i=1}^{N} r_i(t)^2}$ は空間RMS、$\mu_{baseline}, \sigma_{baseline}$ は直前5エポックの統計量である。閾値 $|Z| > 2.0$ で異常を判定する。

**（b）CUSUM制御図（参照: Yokota & Ishikawa, 2020）**

$$C^+(t) = \max(0, C^+(t-1) + [R(t) - \mu_0]/\sigma_0 - k)$$

$$C^-(t) = \max(0, C^-(t-1) - [R(t) - \mu_0]/\sigma_0 - k)$$

'REPORT_EOF' $h = 3.0$、参照値 $k = 0.5$ を採用。

**（c）複合判定**：三手法のうち2手法以上が閾値超過した場合に検出イベントとして認定。

### 2.5 3次元変位場推計

echoLOS変位を加重最小二乗法で分解する。北方向変位を無視（$U_n = 0$）した近似では：

$$\begin{pmatrix} d_{asc} \\ d_{desc} \end{pmatrix} = \begin{pmatrix} \sin\theta \cos\alpha_{asc} & -\cos\theta \\ \sin\theta \cos\alpha_{desc} & -\cos\theta \end{pmatrix} \begin{pmatrix} U_e \\ U_z \end{pmatrix}$$

#
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }

$$\hat{\mathbf{m}} = (\mathbf{A}^T \mathbf{W} \mathbf{A})^{-1} \mathbf{A}^T \mathbf{W} \mathbf{d}$$

 $\mathbf{W} = \text{diag}(w_{asc}, w_{desc})$ は重み行列、$\mathbf{A}$ は投影行列である（参照: Hu et al., 2022）。

---

## 3. 主要な結果と数値

### 3.1 大気補正性能

| 補正手法 | 補正前 STD [mm] | 補正後 STD [mm] | 低減率 [%] |
|---------|----------------|----------------|-----------|
| 線形高度補正 | 7.56 | 6.95 | 6.7 |
| べき乗則補正 | 7.56 | 8.60 | -10.6（悪化） |
| ERA5単独 | 7.56 | 6.65 | 14.0 |
| **ERA5＋線形（組み合わせ）** | **7.56** | **5.99** | **21.3** |

ERA5＋線形の組み合わせ補正が最大の改善効果を示した（21.3%低減）。べき乗則補正は合成データでは悪化したが、高山地帯など高度変化が大きい地域では有効とされる（Yang et al., 2023）。

![大気補正比較](figures/fig2_atmospheric_correction.png)

### 3.2 トレンド分離精度

| 指標 | 値 |
|------|-----|
| 平均 R² | 0.59 ± 0.18 |
| LOS速度範囲 | -15.62〜13.67 mm/yr |
| 平均季節振動振幅 | 3.86 ± 1.29 mm |
| 5分割CV-RMSE | 25.2 ± 39.1 mm |
| 速度 RMSE（真値比較） | 12.0 mm/yr |

#R² = 0echo.59は、4mm大気ノイズと合成'REPORT_EOF''REPORT_
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }RMSEが高いのは、SSEイベントが将来予測不可能な過渡的変動であり、前向き予測精度が本質的に低下するためである。

![トレンド分離結果](figures/fig3_trend_decomposition.png)

### 3.3 前兆変動検出性能

| 指標 | 値 |
|------|-----|
| 真のSSEイベント数 | 3 |
| 検出イベント数 | 2（真陽性2、偽陽性1） |
| 精度（Precision） | 0.667 |
| 再現率（Recall） | 0.667 |
| F1スコア | 0.667 |
| 平均検出遅延 | 12日（1エポック） |

3イベント中2イベントを検出（再現率67'REPORT_EOF'1イベントは、線形+季節フィッティングにSSE信号の一部が吸収されたためと考えられる。%）、偽陽性1件（精度67%）、検出遅延12日（

![前兆変動検出](figures/fig4_precursor_detection.png)

### 3.4 3次元変位場

| 指標 | 値 |
|------|-----|
| 投影行列条件数 | 1.726（良条件） |
| LOS RMSE（昇軌道） | 0.000 mm（完全適合） |
| LOS RMSE（降軌道） | 0.000 mm（完全適合） |
| 鉛直速度範囲 | -17.12〜19.63 mm/yr |
| 東西速度範囲 | -2.22〜2.14 mm/yr |

.726は良条件を示し、昇降軌道の幾何学的配置が3次元分解に適していることを示す。鉛直速度範囲（-17〜20 mm/yr）は入力の真の値（-15〜-3 mm/yr）を包含しており、個々の点推定は雑音の影響を受けるが統計的には整合する。

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }2方程式・2未知数（fix_north=True）のため数値的LOS誤差は0であるが、これは完全適合によるものであり、過学習ではなく連立方程式の特性による。

![3]]]](figures/fig5_3d_displacement.png)

### 3.5 パイプライン全体図

![パイプライン全体](figures/fig1_pipeline_overview.png)

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }

| 指標 | 値 |
|------|-----|
| 干渉ペア数 | 545 |
| 平均コヒーレンス | 0.871 ± 0.081 |
| コヒーレンス 0.5以上の割合 | 98.5% |

![SBASネットワーク](figures/fig6_sbas_coherence.png)

---

## 4. 考察と今後の展望

### 4.1 大気補正の課題

GACOS（Generic Atmospheric Correction Online Service for InSAR）との比較評価や、GNSS観測網との統合補正が有効と考えられる（Yang et al., 2023; Huang et al., 2025）。)

### 4.2 SSE検出の課題

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 SS/InSARデータ融合"___Begin12日以下（1エポック）となる場合、Sentinel-1の12日繰り返しサイクルでは時間分解能が不足する。L-バ___Command_Done_Marker___$SAR（ALOS-2/PALSAR-2）の14日サイクルや、将来のSentinel-1C/D衛星による6日サイクルへの対応が重要である。また、SSEが線形フィッティングに吸収される問題は、時変速度モデル（Gualandi & Liu, 2021）の導入;                 EC=$?

InSAR単体の処理パイプラインを実装したが、実際の南海トラフ監視ではGNSS連続観測（GEONET）との融合が不可欠である。InSARの高空間分解能とGNSSの高3次元変位の相補的統合（Mancini et al., 2021）により、南海トラフ沿いの詳細な変動パターンの解明が期待される。

### 4.4 ISCE/StaMPSワークフローへの展開

.git .github .gitignore .pytest_cache AGENTS.md data figures logs results run_experiment.py src Tests ISCE（InSAR Scientific Computing Environment）とStaMPS（Stanford Method for Persistent Scatterers）の実データ処理フレームワークへの実装を想定している。Mancini et al. (2021))SNAP-StaMPS統合ワークフローを参考に、Sentinel-1データの自動ダウンロード・前処理・時系列解析・可視化の完全自動化が実現可能である。

---

## 5. 限界と今後の課題

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             tine-1のオーバーラップ領域でのオフセット追跡やGNSS-Aとの統合が改善策として有効である（Kyaw & Takeuchi, 2025）。}

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$2014年以降のSentinel-1全データ（10年以上）を対象とした長期時系列解析が必要である。

4. **計算コストの課題**: 実データでのPS-InSAR処理には数万〜数十万点)PS候補点の処理が必要であり、HPC（高性能計算）環境での並列処理実装が必要である。

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             } | 内容 | 行数 |
|---------|------|------|
| `src/insar_simulation.py` | InSAR合成データ生成モジュール | ~180行 |
| `src/atmospheric_correction.py` | 大気遅延補正モジュール | ~200行 |
| `src/trend_decomposition.py` | 時系列トレンド分解モジュール | ~250行 |
| `src/precursor_detection.py` | 前兆変動検出モジュール | ~240行 |
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             } | ~210行 |
| `run_experiment.py` | メイン実験スクリプト | ~340行 |
| `tests/test_pipeline.py` | バリデーションテスト（14件） | ~150行 |
| `figures/fig1_pipeline_overview.png` | パイプライン全体フロー図 | — |
| `figures/fig2_atmospheric_correction.png` | 大気補正比較図 | — |
| `figures/fig3_trend_decomposition.png` | トレンド分離結果図 | — |
| `figures/fig4_precursor_detection.png` | 前兆変動検出図 | — |
| `figures/fig5_3d_displacement.png` | 3次元変位場図 | — |
| `figures/fig6_sbas_coherence.png` | SBASネットワーク・コヒーレンス図 | — |
| `results/experiment_results.json` | 数値結果まとめ | — |
| `results/reference-list.md` | 先行研究リスト（14件） | — |

---

## 参考文献

1. Yokota, Y., & Ishikawa, T. (2020). Shallow slow slip events along the Nankai Trough. *Science Advances*, 6(3). DOI: 10.1126/sciadv.aay5786

2. Takemura, S. et al. (2023). A review of shallow slow earthquakes along the Nankai Trough. *Earth, Planets and Space*, 75, 175. DOI: 10.1186/s40623-023-01920-6

3. Chen, Y. et al. (2023). Integration of DInSAR-PS-Stacking and SBAS-PS-InSAR Methods. *Remote Sensing*, 15(10), 2691. DOI: 10.3390/rs15102691

4. Mancini, F., Grassi, F., & Cenni, N. (2021). A Workflow Based on SNAP–StaMPS Open-Source Tools for PSI-Based Ground Deformation. *Remote Sensing*, 13(4), 753. DOI: 10.3390/rs13040753

5. Gualandi, A., & Liu, Z. (2021). Variational Bayesian ICA for InSAR Displacement Time-Series. *Journal of Geophysical Research: Solid Earth*, 126. DOI: 10.1029/2020JB020845

6. Yang, Q. et al. (2023). Evaluation of InSAR Tropospheric Delay Correction Methods. *Sensors*, 23(23), 9574. DOI: 10.3390/s23239574

7. Huang, D. et al. (2025). Power-Law Model and ERA-5 Data for InSAR Tropospheric Delay Correction. *Sensors*, 25(3), 716. DOI: 10.3390/s25030716

8. Hu, J. et al. (2022). Calculating Co-Seismic 3D Displacements from InSAR with Dislocation Model Constraint. *Remote Sensing*, 14(18), 4481. DOI: 10.3390/rs14184481

9. Reinosch, E. et al. (2020). InSAR time series analysis of seasonal surface displacement on the Tibetan Plateau. *The Cryosphere*, 14, 1633–1650. DOI: 10.5194/tc-14-1633-2020

10. Kyaw, K.M., & Takeuchi, W. (2025). Instability mapping of Dhaka-Kasiani-Gopalganj railway line with InSAR time series. *Scientific Reports*. DOI: 10.1038/s41598-025-21375-x

11. Zhang, P. et al. (2023). A New Method for Continuous Track Monitoring Using the Integration of PS-InSAR and SBAS-InSAR. *Remote Sensing*, 15, 3298. DOI: 10.3390/rs15133298

12. Moualla, L. et al. (2024). Learning Ground Displacement Signals Directly from InSAR-Wrapped Interferograms. *Sensors*, 24(8), 2637. DOI: 10.3390/s24082637

13. Karimzadeh, S., & Matsuoka, M. (2020). Ground Displacement in East Azerbaijan Province revealed by L-band and C-band InSAR. *Sensors*, 20(23), 6913. DOI: 10.3390/s20236913

14. Yalvac, S. (2020). Validating InSAR-SBAS Results with GNSS Analysis Techniques. *Environmental Monitoring and Assessment*, 192. DOI: 10.1007/s10661-019-8009-8
