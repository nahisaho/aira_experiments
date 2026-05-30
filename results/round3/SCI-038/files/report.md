# ADR Mission Optimal Trajectory Design Report

## Abstract / 要旨
 では Active Debris Removal, ADR, ミッション の ための 統合 軌道設計 システム を Python で 実装 し, 再現可能 な 計算 実験 として 実行 した. 目的 は, デブリ 優先順位付け, 訪問 順序 最適化, 低推力 軌道遷移,  ランデブー, タンブリング 解析, 捕獲 手法 比較 を 単一 ワークフロー に 接続 し, 初期 設計 レベル で どのような トレード が 見えるか を 示す ことである. 実験 では 20 個 の 合成 LEO デブリ を 生成 し, リスク と 除去容易性 を 組み合わせた スコア により 10 個 の 重点 目標 を 選定 した. 単発 の シーケンス 最適化 では nearest-neighbor が 12.963 km/s, GA+2opt が 12.848 km/s で, 改善率 は 0.89% であった. ただし 6 seed 感度解析 では NN 平均 14.834 ± 4.614 km/s, GA+2opt 平均 14.807 ± 4.966 km/s, 平均 差 0.028 km/s, 95% CI [-0.438, 0.493], p = 0.885 となり, 進化計算 の 優位性 は この 合成 設定 では 統計 的 に は 強く 支持 されなかった. 低推力 区間 では 400 km, 51.6 deg から 約 800 km, 71 deg への 遷移 を 近似 Q-law で 計算 し, Δv 0.318 km/s, 時間 441.9 h, 推進剤 消費 10.76 kg を 得た. Hill ランデブー では 初期 距離 5001.2 m を 4532.5 s で 終端 に 近づけ, 合計 インパルス は 1262.1 mm/s であった. 姿勢 ダイナミクス では 初期 回転率 5.67 deg/s に対し, 卓越 周期 300.1 s を 推定 した. 捕獲 手法 比較 では harpoon が robotic arm より 平均 0.169 ± 0.157 高い 成 を 示した が, Bonferroni 補正 後 の p 値 は 0.092 であり, 強い 確証 とまでは 言えない. 以上 より, 本システム は 高忠実度 フライト ツール ではない が, ADR の 統合 概念設計 と 感度解析 の ベース近   として 有用 である.

## Introduction / 背景と目的
'MD''MD'echo 持続性 に対する 軌道 デブリ の 脅威 は 年々 大きく なっている. Murtaza et al. (2020) は 軌道 デブリ が 将来 の 宇宙 活動 と 保echo コスト に 与える 長期 的 な 影響 を 示し, 単なる 破片 数 の 問題 ではなく, 社会 的 な 持続性 問題 でもある と 指摘 した. Shan et al. (, ハープーン, ロボットアーム, ドラグ増強 など が それぞれ 異なる 対象 条件 と 制御 要求 を 持つ こと を 整理 した. Aglietti et al. (2019) の RemoveDEBRIS 実証 は, 技術 実証 レベル では あ, ADR が 純粋 な 理論 問題 ではなく, 宇宙空間 で 試験 可能 な 工学 問題 である こと を 示した. Papadopoulos et al. (2021) も また, 宇宙 ロボット 捕獲 では 姿勢 推定, 接触 力学, 遅延 制御 が 支配 的 な 難しさ である と 指摘 している.

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             } に対する Lyapunov 型 制御 の 方向性 を 与えている.

, 多く の 研究 は それぞれ 重要 な 一部分 を 深く 扱う 一方, 実際 の 概念設計 で 必要 な 統合 パイプライン を 直接 提供 していない. すなわち,  の 高い デブリ を どう 選ぶか, どの 順番 で 訪れるか, 低推力 による 遷移 を どう 見積もるか, 近傍 で どう 接近 するか, タンブリング 目標 に どの 捕獲 手法 が 向くか を 同時 に 評価 する ベースライン が 少ない. リス の 目的 は, それら の サブ問題 を Python モジュール に 分けて 実装 し, なおかつ 単一 の 実験 として 結合 する ことで, 初期 設計 に 必要 な 透明性 と 再現性 を 与える ことである. ここで の 主要 な 仮定 , 合成 カタログ を 用いても サブシステム 間 の 相互依存 を 観察 する には 十分 であり, 初期 段階 の 方法比較 には 意味 が ある という 点 に ある.

## Methods / 手法
 は `debris_catalog.py`, `orbital_mechanics.py`, `debris_dynamics.py`, `mission_optimizer.py`, `visualization.py`, `main.py` を 中核 とする. `tests/test_adr_system.py` では モジュール の 妥当性 を 最低限 検証 した. 依存 ライブラリ は numpy, scipy, matplotlib のみ とし, seed 42 を 基本 実験 の 既定値 とした. 文献 探索 は まず `SemanticScholar_search_papers` を 試行 した が HTTP 400 により 失敗 し, その後 `openalex_literature_search` と `Crossref_search_works` を 用いて 参考文献  DOI を 確認 した. この 手順 は, 文献収集 が 単なる 背景説明 ではなく, 実装 する モジュール の 選択 と 数式 の 定義 に 影響 を 与えた こと を 明示 する ため に 記録 している.

.git .github .gitignore .pytest_cache .venv AGENTS.md data figures logs report.md results src tests  優先度 は 衝突確率, 断面, 軌道 高度, 除去 難易度 を 正規化 して 線形 結合 した. 優先度 関数 は 次式 で 与える.

$$
P_i = \alpha \tilde c_i + \beta \tilde s_i + \gamma (1-\tilde h_i) + \delta (1-\tilde d_i)
$$

 $\tilde c_i$ は 正規化 衝突確率, $\tilde s_i$ は 正規化 断面積, $\tilde h_i$ は 正規化 高度, $\tilde d_i$ は 正規化 除去難易度 , 重み は $\alpha=0.4$, $\beta=0.3$, $\gamma=0.2$, $\delta=0.1$ とした. 補助 指標 として, 高衝突確率 かつ 大きな 断面積 を 強調 する ため, 次の リスク 指標 も 計算 した.で

$$
R_i = \tilde c_i \sqrt{\tilde s_i}
$$

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";                      最適化 では nearest-neighbor を ベースライン とし, それに 対して genetic algorithm と 2-opt 局所探索 を 適用 した. 各 レグ の 近似 Δv は Edelbaum 型 の 低推力 近似 と ランデブー オーバーヘッド の 和 として 評価 した.}

$$
\Delta v_{ij} \approx \sqrt{v_i^2 + v_j^2 - 2 v_i v_j \cos\left(\frac{\pi}{2}|\Delta i_{ij}|\right)} + \Delta v_{\mathrm{ops}}
$$

 $v_i=\sqrt{\mu/a_i}$, $v_j=\sqrt{\mu/a_j}$ である. 候補 手法 として は 厳密 な 混合整数 最適化 と Pointer Network 型 学習手法 も 考えられる が, 前者 は 軽量 な 再現 実装 には 重く, 後者 は 学習 データ と 訓練 コスト を 要する. そのため, 本研究 では  と 拡張性 を 優先 して heuristic + evolutionary search の 組合せ を 採用 した. この 判断 は, 概念設計 段階 では 「計算可能 で 説明 可能」な 手法 が 有利 という 工学 的 要請 に 基づく.

 区間 では Q-law に 着想 を 得た 近似 指標 を 用いた.

$$
Q = \left(\frac{a-a_t}{a_t}\right)^2 + \left(\frac{e-e_t}{1-e_t}\right)^2 + \left(\frac{i-i_t}{\pi/6}\right)^2
$$

 消費 は Tsiolkovsky 式 により 更新 した.

$$
\Delta m = m \left(1 - \exp\left(-\frac{\Delta v}{I_{sp} g_0}\right)\right)
$$

 ランデブー には Hill-Clohessy-Wiltshire 方程式 を 採用 し, 閉形式 の 状態遷移 行列 から 二回 インパルス 解 を 計算 した. 非線形 相対運動 モデル や J2 摂動 付 高忠実度 モデル も 候補 だが, 本研究 は 概念設計 が 目的 である ため, 線形 モデル の 可視化 と 解釈 性 を 優先 した. タンブリング 解析 では Euler 剛体 方程式 と quaternion kinematics を RK4 で 積分 し, 回転 周期 は periodogram により 推定 した. 統計 解析 では 6 seed の ペアデータ に 対して Shapiro-Wilk 正規性 検定 を 行った 後, paired t-test を 適用 した. 'MD' 比較 に対して は Bonferroni 補正 も 報告 した. 以上 の 選択 により, 本システム は 完全 性 より も モジュール 間 の 接続 性 と 再現可能 性 を 優先 する 設計 となっている.

## Results / 結果
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             } を 大きく 左右 する こと を 示している.

![Figure 1: デブリカタログ概観](figures/debris_catalog_overview.png)

#   順序 最適化 の 基本 実験 では NN が 12.963 km/s, GA+2opt が 12.848 km/s で, 単発 改善率 は 0.89% であった. しかし seed を 7, 11, 19, 23, 31, 42 に 変えた 感度解析 では, NN 平均 は 14.834 ± 4.614 km/s, 95% CI [9.992, 19.677], GA+2opt 平均 は 14.807 ± 4.966 km/s, 95% CI [9.596, 20.018] であった. 平均 差 は 0.028 ± 0.444 km/s, 95% CI [-0.438, 0.493], p = 0.885, Cohen's dz = 0.062 であり, 基本 実験 に 見られた 改善 は seed 横断 では 安定 的 な 優位 とまでは 言えない. この 結果 は, greedy baseline が すでに 強い 一         で, 進化計算 は 一部 の 幾何条件 では 有利 でも, 一様 に 優れる わけでは ない こと を 意味 する. 逆に 言えば, GA+2opt は ベースライン を 置き換える 魔法 の 手法 では なく, 既存 解 の 追加 改善 を 試みる 
 的 層 として 理解 する のが 適切 である.

![Figure 2: ミッション順序最適化](figures/mission_sequence_optimization.png)

 遷移 の 基本 実験 では Δv 0.318 km/s, 所要時間 441.9 h, 推進剤 消費 10.76 kg であり, 最終 傾斜角 は 70.97 deg に 到達 した. thrust を 0.18, 0.20, 0.22 N と 変化 させた とき, 時間 は 462.1 h, 441.9 h, 401.7 h と 単調 に 短縮 し, 最終 傾斜角 は 69.82–70.97 deg に 保たれた. この 振る舞い は, 近似 モデル で あっても パラメータ に 対する 応答 が 破綻 していない こと を 示す. また, 燃料 消費 が 約 10–11 kg に 収まる 点 は, 長期 の 電気 推進 遷移 を 模擬 する 初期 設計 として 妥当 な オーダー である. 高推力 側 では 時間短縮 の 利得 が 大きく, Δv と 質量 消費 は 大きく 崩れない ため, 簡略 モデル の 範囲 では thrust margin が 運用 余裕 に 変換 される こと も 読み取れる.

![Figure 3: 低推力軌道遷移](figures/low_thrust_trajectory.png)

Hill ランデブー 解析 では 初期  5001.2 m, 移動 時間 4532.5 s, 第1 インパルス 661.0 mm/s, 第2 インパルス 601.1 mm/s, 合計 1262.1 mm/s を 得た. 最終 距離 は 9.17e-13 m と 数値 的 に ほぼ 0 であり, 線形 モデル の 閉形式 解 が 初期 接近 設計 に 十分 利用 できる こと を 示す. 特に, 5 km 級 の 初期 オフセット に対して 1.3 m/s 未満 で 終端 条件 を 満たす という 結果 は, 近傍 運用 の 燃料 予算 を ミ   全体 の Δv から 切り分けて 考える 上で 有用 である. 線形 化 の 制約 は あるが, 早期 設計 に 必要 な 量級 把握 と 可視化 という 目的 には よく 合致 した.

![Figure 4: Hillランデブー解析](figures/hill_rendezvous.png)

 ダイナミクス では 初期 回転率 5.67 deg/s, 推定 周期 300.1 s, 最大 観測 回転率 5.72 deg/s であった. 図 5 の 相図 と Euler 角 応答 は, 非対慣性 を 持つ 目標 が 完全 に 一定 角速度 で 回る わけでは なく, 小さな 観測 ノイズ を 含めても 周波数 構造 が 読み取れる こと を 示す. この 情報 は 捕獲 手法 の 条件付き 選択 に 直結 する ため, mission-level の 序列化 と 切り離す べき では ない. すなわち, タンブリング 状態 は capture module の 入力 である と 同時に, 目標 価値 を 再評価 する ための mission planning input でもある.

![Figure 5: タンブリングデブリの姿勢応答](figures/tumbling_debris_dynamics.png)

  手法 比較 では 1–18 deg/s の 回転率 条件 で harpoon が 常に 最高 スコア を 得た. harpoon と robotic arm の 成功確率 差 は 平均 0.169 ± 0.157, 95% CI [0.004, 0.334], raw p = 0.046, Bonferroni 補正 後 p = 0.092, Cohen's dz = 1.076 であった. 効果量 は 大きい が, 補正 後 の 有意性 は 境界 的 である. したがって, 本結果 は 強い 結論 ではなく, 高スピン 条件 では standoff 距離 を 許容 する 手法 が 優位 に なる 可能性 を 示す 設計 指echo と 解釈 する のが 適切 である. また 18 deg/s 条件 では harpoon も 大きく 性能 を 落とす ため, 「どの 手法 でも 高回転 目標 は 難しい」という より 基本 的 な 事実 も 読み取れる.

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }](figures/capture_mechanism_analysis.png)

 ミッション 指標 として, 総 Δv は 12.848 km/s, 推進剤 要求量 は 353.8 kg, 10 目標 に対する 推定 ミッション 期間 は 526.5 日 であった. 全捕獲 成功率 の 連鎖 値 は 0.204% と 非常 に 低い が, これは 各 ターゲット の 成功確率 を 独立 に 乗算 した 保守 的 仮定 による. むしろ この 数値 は, 多目標 ADR が Δv 最適化 だけでなく, 失敗 許容 型 の 再計画 や 部分 完遂 指標 を 必要 とする こと を 示唆 する. 図 7 は Δv 予算, タイムライン, リスク 低減, 燃料, 成功率 を まとめた ダッシュボード であり, どの サブシステム が mission closure を 制約 しているか を 一目 で 確認 できる.

![Figure 7: ミッション要約ダッシュボード](figures/mission_summary_dashboard.png)

## Discussion / 考察
      の 重要 な 点 は, ADR の 設計 判断 が 単一 の サブ問題 では 決まらない こと を 実験 的 に 示した ところ に ある. シーケンス 最適化 だけ を 見れば, 基本 実験 では GA+2opt が NN を わずか に 上回った.  seed 横断 では 差 の 95% CI が 0 を またぎ, p 値 も 大きく, 「進化計算 が 常に 優位」と は 言えない. これは 実務 的 に 重要 であり, 初期 検討 では まず greedy heuristic で 妥当 な 解 を 作り, その後 に 限定 的 な 改良探索 を 加える 方針 が 合理 的 である こと を 示唆 する. Zhao et al. (2020), Medioni et al. (2022), Zona et al. (2023) が 示した ように, 多目標 ADR では 探索 手法 の 複雑 化 より, 問題 分解 と 評価 関数 の 設計 が 同じ くらい 重要 である.

, 低推力 遷移 と Hill ランデブー を 同じ パイプライン に 含めた ことで, 大域 的 な 軌道変換  運用 を 一続き に 見る こと が できた. Narayanaswamy & Damaren (2023) が 示した 低推力 誘導 の 文脈 と 比較 すると, 本モデル は かなり 単純 である. それでも, 低推力 区間 の Δv, 時間, .git .github .gitignore .pytest_cache .venv AGENTS.md data figures logs report.md results src tests  消費 が 明確 に 出力 され, ランデブー の mm/s オーダー の 操作量 と 同じ 結果 ファイル に 収まる 点 は, 概念設計 として 有用 である. すなわち, 「順序」, 「遷移」, 「接近」 を 同一 指標 系 の 中で 比較 . これは 将来, さらに 高忠実度 の 伝播 や 制約 条件 を 追加 しても, 基本 的 な ソフトウェア 構造 が そのまま 拡張 できる こと を 意味 する.できと 

  に関して は, Shan et al. (2016), Aglietti et al. (2019), Papadopoulos et al. (2021) が 指摘 した ように, 目標 の 回転 が 技術 選択 を 支配 する. 本モデル は 簡略 化 されている が, robotic  に なり, harpoon 系 が 相対 的 に 有利 と なる という 傾向 を 再現 した. したがって, 将来 の mission optimizer は Δv 最小化 だけでなく, target-specific capture feasibility を 組み込んだ expected mission value を 最適化 する 方向 へ 発展 させる べき である. Chen et al. (2024) や Guo et al. (2023) の ような 学習型 あるいは 部分 捕獲 の echo   は, まさに その 拡張 に 関連 する. 本結果 は 強い 実証 では ない が, どこに 次の モデル 拡張 を 入れるべきか という 設計 指針 を 与える.

 に, Simha et al. (2025) が 強調 した 政策 的 観点 から 見る と, 本研究 の ような 軽量 統合 モデル は 複数 シナリオ の 迅速 な 比較 に 向く. 高忠実度 の 飛行解析 だけ では, 政策 や 投資 判断 の 初期 段階 で 必要 な 多数 の what-if 分析 を こなす こと が . 一方, 本システム の ような 再現可能 ベースライン は, どの サブシステム が ボトルネック か を 可視化 し, 追加 投資 の 優先順位 を 議論 する 足場 を 提供 する. この 点で 本研究 は, 厳密 解 の 提示 も,  より      と 意思決定 を つなぐ 計算 基盤 を 提供 した と 位置づける のが 妥当 である.

## Limitations and Future Work / 限界と今後の課題
### Data Limitations
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             } TLE 誤差, 観測 欠損, 質量 推定 不確実性, 形状 情報 の 欠如, 断面積 の 時変性 を 反映 していない. 20 物体 という 規模 も 実運用 の カタログ に 比べて 小さい. そのため, 優先 スコア の 値 そのもの を 実運用 の ランキング と 見なす こと は できない. 本結果 は あくまで モジュール 統合 と 感度解析 の ベースライン である.

### Methodological Limitations
 を 重視 して 意図的 に 簡略 化 されている. 低推力 伝播 は averaged element update に 基づき, J2, eclipse, thrust steering saturation, mass-flow schedule を 明示 的 に 扱わない. Hill 方程式 は 円軌道 近似 と 線形 化 を 仮定 する ため, 長時間 あるいは 大規模 相対距離 では 誤差 が 増える. 捕獲 モデル も 接触 力学, 柔軟 構造, センサ ノイズ, 閉ループ 制御 遅延 を 含まない. この ため, 結果 は 「どの 条件 で 何が 効きそうか」を 示す が, flight qualification を 支える もの では ない.は 

### Evaluation Limitations
 指標 は 主に Δv, 時, 燃料, 成功確率 に 限られている. launch window, 通信 可視性, 故障 対応, 法制度, ターゲット 所有権, ミッション 保険 条件 は 組み込まれて いない. シーケンス 解析 の seed 感度 でも, 差 の 95% CI  0 を またいで おり, optimizer の 優位性 は 強く は 証明 されていない. External validation with independent real-world datasets is essential to confirm the generalizability of these findings beyond simulated conditions. この 文を そのまま 書く 理由 は, 本研究 が まさに 合成 条件 の 範囲 に ある から である. 現段階  現実 データ による 外部 妥当性 が 欠けている.

### Generalizability
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";              は LEO の 非協力 目標 を 想定 した ADR 概念検討 に 適している が, GEO servicing, cislunar cleanup, cooperative servicing に そのまま 適用 できる とは 限らない. 軌道 環境, 推進 方式, 姿勢 安定性, 観測 幾何 が 異なれば, 優先度 関数 と ダイナミクス モデル を 再設計 する 必要 が ある. 解析 も 単一 慣性 テンソル に 基づく ため, 大型 太陽電池 パドル や 非剛体 構造 を 持つ 目標 には 直ちに 一般化 できない.}

### Future Directions
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             } , J2 を 含む 軌道 遷移 モデル, partial capture policy, capture success の 条件付き expected value を 追加 する べき である. 中長期 的 には Chen et al. (2024) の 学習型 sequence generation や Guo et al. (2023) の 部分 捕獲 戦略 を 組み込み, Simha et al. (2025) の ような 政策 志向 指標 と 接続 する こと が 望ましい. また, 実験  harpoon 優位 の 傾向 が 示された が, これは 実験 装置 や 高忠実度 接触 シミュレーション で 検証 される 必要 が ある. 将来的 には 実データ, 高忠実度 伝播, 不確実性 伝播, 法制度 制約 を 一体 化 した multi-objective ADR platform へ 拡張 すべき である.

## References / 参考文献
1. Wijayatunga, M., Armellin, R., Holt, H., Pirovano, L., & Lidtke, A. (2023). Design and guidance of a multi-active debris removal mission. *Astrodynamics*. DOI: 10.1007/s42064-023-0159-3
2. Medioni, L., Gary, Y., Monclin, M., Oosterhof, C., Pierre, G., Semblanet, T., Comte, P., & Nocentini, K. (2022). Trajectory optimization for multi-target Active Debris Removal missions. *Advances in Space Research*. DOI: 10.1016/j.asr.2022.12.013
3. Zona, F., Zavoli, A., & Federici, L. (2023). Evolutionary Optimization for Active Debris Removal Mission Planning. *IEEE Access*. DOI: 10.1109/access.2023.3269305
4. Zhao, Z., Feng, F., & Yuan, J. (2020). A Novel Two-Level Optimization Strategy for Multi-Debris Active Removal Mission in LEO. *CMES*. DOI: 10.32604/cmes.2020.07504
5. Guo, Z., Pang, B., & Du, X. (2023). Optimal planning for a multi-debris active removal mission with a partial debris capture strategy. *Chinese Journal of Aeronautics*. DOI: 10.1016/j.cja.2023.03.013
6. Chen, S., Bai, X., & Zhao, Y. (2024). Rapid Sequence Generation for Active Debris Removal Mission Based on Attention Mechanism and Pointer Network. *IEEE Access*. DOI: 10.1109/access.2024.3425161
7. Papadopoulos, E., Aghili, F., Ma, O., & Lampariello, R. (2021). Robotic Manipulation and Capture in Space: A Survey. *Frontiers in Robotics and AI*. DOI: 10.3389/frobt.2021.686723
8. Aglietti, G. S., Taylor, B., Fellowes, S., et al. (2019). RemoveDEBRIS: An in-orbit demonstration of technologies for the removal of space debris. *The Aeronautical Journal*. DOI: 10.1017/aer.2019.136
9. Murtaza, A., Pirzada, S. J. H., Xu, T., & Liu, J. (2020). Orbital Debris Threat for Space Sustainability. *IEEE Access*. DOI: 10.1109/access.2020.2979505
10. Narayanaswamy, S., & Damaren, C. J. (2023). Equinoctial Lyapunov Control Law for Low-Thrust Rendezvous. *Journal of Guidance, Control, and Dynamics*. DOI: 10.2514/1.g006662
11. Simha, A., Servadio, S., & Lifson, M. (2025). Optimal Active Debris Removal mission planning to inform policy decisions. *Acta Astronautica*. DOI: 10.1016/j.actaastro.2024.11.050
12. Shan, M., Guo, J., & Gill, E. (2016). Review and Comparison of Active Space Debris Capturing and Removal Methods. *Progress in Aerospace Sciences*. DOI: 10.1016/j.paerosci.2015.11.001

## File Inventory / 生成ファイル一覧
- Source modules: `src/debris_catalog.py`, `src/orbital_mechanics.py`, `src/debris_dynamics.py`, `src/mission_optimizer.py`, `src/visualization.py`, `src/main.py`, `src/generate_documents.py`
- Tests: `tests/test_adr_system.py`
- Data: `data/debris_catalog.json`, `data/preprocessing-log.md`
- Results: `results/experiment_results.json`, `results/search-strategy.md`, `results/reference-list.md`, `results/statistical-summary.md`, `results/sensitivity-analysis.md`, `results/summary_metrics.json`
- Figures: `figures/debris_catalog_overview.png`, `figures/mission_sequence_optimization.png`, `figures/low_thrust_trajectory.png`, `figures/hill_rendezvous.png`, `figures/tumbling_debris_dynamics.png`, `figures/capture_mechanism_analysis.png`, `figures/mission_summary_dashboard.png`
- Logs: `logs/process-log.jsonl`
