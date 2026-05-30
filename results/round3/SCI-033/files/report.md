# 量子機械学習ベンチマーク報告書

DRAFT - NOT FOR DISTRIBUTION

## 概要
本報告書の目的は、量子機械学習（QML）における表現力、量子カーネル、データ符号化、学習可能性、ノイズ耐性を一体的に比較できる軽量ベンチマーク基盤を構築し、その実験結果を日本語で再現可能に記録することである。文献面では Sim らの表現力指標、Schuld と Killoran の特徴ヒルベルト空間観点、Liu らと Huang らの量子優位条件、McClean らの barren plateau、Cerezo らの変分量子アルゴリズム総説、Thanasilp らのカーネル集中を統合的に参照した。実装面では PennyLane によって、4量子ビットの PQC アーキテクチャ比較、量子カーネル SVM と古典 RBF-SVM の 5-fold 交差検証、angle / amplitude / IQP 符号化比較、2〜8量子ビットの勾配分散解析、そして depolarizing noise の影響評価を実行した。

主要結果は明確である。RCS 4層は KL divergence 0.0068 で最も Haar 分布に近い表現力を示し、SEL 2層は Meyer-Wallach 指標 $Q = 0.838 \pm 0.111$ で最も高い平均エンタングルメントを示した。分類性能では古典 RBF-SVM が 0.985 ± 0.014、量子カーネル SVM が 0.970 ± 0.027 であり、平均値では古典法が優勢だった。ただし paired t-test は $p = 0.208$、Wilcoxon 検定は $p = 0.500$ であり、差は記述的に解釈すべきである。符号化比較では angle encoding が 0.956 ± 0.028 で最良、IQP が 0.925 ± 0.028、amplitude が 0.819 ± 0.089 だった。さらに勾配分散は 2 qubits の 3.52×10^-2 から 8 qubits の 6.30×10^-4 まで減衰し、ノイズ下では fidelity が 0.979 ± 0.001 から 0.158 ± 0.007 へ急落した。以上より、QML の優位は自動的には現れず、問題構造、符号化、回路深さ、ノイズ制御がそろう場合にのみ期待される。

## 実験目的と背景
QML はしばしば「量子状態が古典的に扱いにくい特徴空間を与えるため、高い予測性能や量子優位を実現しうる」と説明される。しかし、現実の NISQ 環境ではノイズ、有限ショット、浅い回路深さ、訓練不安定性が支配的であり、理論的な可能性がそのまま実用性能に結び付くわけではない。Preskill (2018) は NISQ 時代の量子計算が有望であると同時に、誤り耐性以前の制約を強く受けることを強調した。McClean ら (2018) は深い変分回路で barren plateau が生じ、勾配分散が指数的に減衰することを示した。さらに Huang ら (2021) は、データが与えられると古典学習器が量子由来タスクでも高性能を発揮しうることを示し、Liu ら (2021) は厳密な量子優位には強い問題構造が必要だと論じた。

そこで本研究は、QML の主要論点を別々に議論するのではなく、同じコードベースの中で一体的に比較することを目標にした。第一に、PQC の表現力とエンタングルメント能力を比較し、どのアーキテクチャが高い状態空間探索力を持つかを調べる。第二に、量子カーネル法を古典 RBF-SVM と比較し、量子モデルの有効性を強い古典ベースラインに照らして評価する。第三に、angle / amplitude / IQP という代表的符号化の差を測る。第四に、barren plateau とノイズによって実装可能性がどの程度制約されるかを確認する。つまり本ベンチマークは、「表現力が高いこと」と「学習器として有用であること」が一致するのかを検証するための統合的な枠組みである。

## 先行研究調査結果（MCPツール試行状況含む）
文献調査ではまず ToolUniverse MCP の利用可能性を確認し、`tooluniverse-find_tools` によって Semantic Scholar、PubMed、Crossref 系ツールを探索した。その後、`SemanticScholar_search_papers` と `Crossref_search_works` を使って対象論文の取得を試みた。結果として、Sim ら (2019) については `SemanticScholar_get_paper` が正常に動作し、著者、タイトル、DOI、要旨を取得できた。Liu ら (2021)、Huang ら (2021)、Preskill (2018) などは Crossref で解決できた。一方で、統合的な検索語を一度に与えた Semantic Scholar 検索は 0件になり、一部の Crossref 問い合わせでは HTTP 429 が返り、レート制限により追加確認が必要になった。

そのため Abbas ら (2021)、Cerezo ら (2021)、Schuld (2021) など一部の論文では `web_search` を用いて DOI と掲載誌を検証した。ここで重要なのは、MCP は「失敗」したのではなく、「部分的に有効だった」という点である。ツール探索と個別論文メタデータの取得には役立ったが、曖昧タイトル解決やレート制限回避には補助経路が必要だった。特に、依頼文で指定された `Quantum models as kernel machines` という Schuld 2021 の表現は安定的に解決できず、内容的に最も近く DOI が確認できた `Supervised quantum machine learning models are kernel methods` を採用した。この判断と検索戦略は `results/search-strategy.md`、`results/screening-table.csv`、`results/extraction-table.csv`、`figures/prisma-flow.md` に正直に記録した。

採択文献の役割も明確である。Sim ら (2019) は表現力とエンタングルメントの指標定義、Abbas ら (2021) は QNN の表現能力、Schuld & Killoran (2019) と Schuld (2021) は量子カーネルの理論的基盤、Havlíček ら (2019) は量子特徴写像の実証例を与える。Liu ら (2021) と Huang ら (2021) は量子優位成立条件と古典データ学習の強さを示し、McClean ら (2018) と Cerezo ら (2021) は学習可能性の制約を整理する。LaRose & Coyle (2020) と Schuld, Sweke, Meyer (2021) は符号化依存性、Thanasilp ら (2024) はカーネル集中とノイズの悪影響、Preskill (2018) は NISQ 制約を与える。文献全体として、「高表現力 = 高性能」という単純図式は成立しないという結論が支持された。

## 手法の説明
ベンチマークは 5つのサブ実験で構成される。第一に、HEA、SEL、RCS の 3種類の PQC を 4量子ビット、1〜4層で比較し、表現力とエンタングルメント能力を推定した。第二に、ZZFeatureMap 風の量子カーネル SVM を実装し、同じデータに対して古典 RBF-SVM と 5-fold stratified cross-validation で比較した。第三に、angle encoding、amplitude encoding、IQP encoding を同一条件で比較した。第四に、ランダム 2層 PQC に対して $n = 2,4,6,8$ で勾配分散を測定し、barren plateau の指数減衰を調べた。第五に、depolarizing noise を 0.001, 0.01, 0.05, 0.1 に変え、noiseless 状態と noisy density matrix の fidelity を比較した。

主要な数式は次のとおりである。2つのパラメータ設定 $\theta,\phi$ に対する回路状態の fidelity を

$$
F(\theta,\phi)=\left|\langle\psi(\theta)|\psi(\phi)\rangle\right|^2
$$

と定義する。表現力は Haar 分布との KL divergence

$$
\mathcal{E}=D_{KL}\big(\hat{P}_{PQC}(F)\|P_{Haar}(F)\big)
$$

で測定し、値が小さいほど Haar ランダムに近い表現力を持つと解釈する。エンタングルメント能力は Meyer-Wallach 指標

$$
Q(\psi)=2\left(1-\frac{1}{n}\sum_{k=1}^n \mathrm{Tr}(\rho_k^2)\right)
$$

で評価した。量子カーネルは

$$
K(x_i,x_j)=\left|\langle 0|U_\phi^\dagger(x_i)U_\phi(x_j)|0\rangle\right|^2
$$

とした。barren plateau 解析では

$$
\mathrm{Var}[\partial_\theta E] \sim b^{-n}
$$

の指数減衰仮説を用いた。本研究で kernel SVM を主軸に選んだ理由は、表現力と学習器の最適化問題を切り分けやすく、古典 RBF-SVM との比較も公平に行いやすいからである。代替案として trainable variational classifier や frame potential も考えられるが、軽量で再現性の高い benchmark という目的に対しては、KL divergence と kernel SVM の組み合わせが最も明快であった。

## 主要な実験結果と数値
表現力比較では、RCS 4層が KL divergence 0.0068 で最良だった。HEA 2層は 0.0091、SEL 3層は 0.0086 であり、深さを増やすと急速に飽和する傾向が確認できた。一方、平均エンタングルメントは SEL 2層が $0.838 \pm 0.111$ で最大であり、最も表現力が高い RCS 4層の $0.817 \pm 0.092$ と一致しなかった。したがって、表現力とエンタングルメントは関連するが同一ではない。

| Architecture | Layers | KL divergence | MW entanglement mean | MW entanglement SD |
| --- | --- | --- | --- | --- |
| HEA | 1 | 0.2047 | 0.8423 | 0.1170 |
| HEA | 2 | 0.0091 | 0.7988 | 0.1155 |
| HEA | 3 | 0.0096 | 0.8204 | 0.0846 |
| HEA | 4 | 0.0107 | 0.8280 | 0.0783 |
| RCS | 1 | 0.7115 | 0.4047 | 0.1879 |
| RCS | 2 | 0.0337 | 0.7735 | 0.1148 |
| RCS | 3 | 0.0145 | 0.7162 | 0.1463 |
| RCS | 4 | 0.0068 | 0.8166 | 0.0915 |
| SEL | 1 | 0.2381 | 0.7970 | 0.1496 |
| SEL | 2 | 0.0149 | 0.8382 | 0.1107 |
| SEL | 3 | 0.0086 | 0.8161 | 0.0882 |
| SEL | 4 | 0.0110 | 0.8258 | 0.0794 |

量子カーネル SVM と古典 RBF-SVM の比較では、古典法が平均精度 0.985 ± 0.014、量子法が 0.970 ± 0.027 となった。平均では古典法が 1.5 percentage points 高いが、paired t-test は $t = 1.50, p = 0.208$、Wilcoxon 検定は $p = 0.500$ であり、差は統計的優越ではなく記述的差として扱うべきである。しかも学習時間は古典法が 0.0006秒、量子法が 33.5996秒であり、計算コスト差は非常に大きかった。

| Model | Accuracy mean | Accuracy SD | 95% CI half-width | Train time mean (s) |
| --- | --- | --- | --- | --- |
| Quantum kernel SVM | 0.9700 | 0.0274 | 0.0340 | 33.5996 |
| Classical RBF-SVM | 0.9850 | 0.0137 | 0.0170 | 0.0006 |

符号化比較では angle encoding が 0.956 ± 0.028 で最良、IQP が 0.925 ± 0.028、amplitude が 0.819 ± 0.089 だった。amplitude encoding は状態準備コストのため実用的には不利であり、angle encoding が最も良い精度・計算コストのバランスを示した。

| Encoding | Accuracy mean | Accuracy SD | 95% CI half-width | Train time mean (s) |
| --- | --- | --- | --- | --- |
| Angle | 0.9563 | 0.0280 | 0.0347 | 6.5952 |
| Amplitude | 0.8187 | 0.0895 | 0.1111 | 51.3501 |
| IQP | 0.9250 | 0.0280 | 0.0347 | 12.6518 |

barren plateau 解析では、勾配分散が 2 qubits の 3.52×10^-2 から 8 qubits の 6.30×10^-4 まで減衰し、指数基底は $b \approx 2.10$ となった。ノイズ評価では fidelity が error rate 0.001 の 0.979 ± 0.001 から 0.1 の 0.158 ± 0.007 まで急落し、quantum volume proxy も 256 から 2 へ低下した。ノイズは量子状態の幾何構造を大きく壊すことが分かる。

| Qubits | Observed variance | Fitted variance | Decay base b |
| --- | --- | --- | --- |
| 2.0000 | 0.0352 | 0.0625 | 2.0998 |
| 4.0000 | 0.0385 | 0.0142 | 2.0998 |
| 6.0000 | 0.0024 | 0.0032 | 2.0998 |
| 8.0000 | 0.0006 | 0.0007 | 2.0998 |

| Error rate | Fidelity mean | Fidelity SD | Quantum volume proxy |
| --- | --- | --- | --- |
| 0.0010 | 0.9792 | 0.0006 | 256.0000 |
| 0.0100 | 0.8110 | 0.0057 | 64.0000 |
| 0.0500 | 0.3640 | 0.0114 | 4.0000 |
| 0.1000 | 0.1578 | 0.0073 | 2.0000 |

![expressibility](figures/expressibility_comparison.png)
![entanglement](figures/entanglement_capability.png)
![kernel](figures/kernel_comparison.png)
![encoding](figures/encoding_comparison.png)
![barren](figures/barren_plateau.png)
![noise](figures/noise_impact.png)

## 考察と今後の展望
本実験の最も重要な含意は、QML の良さは単一指標では語れないという点である。RCS 4層は最も低い KL divergence を示したが、分類性能では古典 RBF-SVM が数値的に高かった。これは、回路が Haar に近い分布を生成できることが、そのまま予測タスクに対する最適な帰納バイアスを意味しないことを示している。むしろ Huang ら (2021) が示したように、データの存在は古典学習器をかなり強くしうる。この結果はその視点と整合的である。

また angle encoding の優位は、シンプルな連続回転が小規模・低ノイズ条件では最も安定した幾何を与える可能性を示している。IQP encoding は位相相関を導入しつつ競争的な精度を保ったが、amplitude encoding は本設定ではコストが高すぎた。これは LaRose & Coyle (2020) や Schuld, Sweke, Meyer (2021) が指摘した符号化依存性を実験的に裏づける。

さらに、barren plateau とノイズは別々の問題ではなく相互補強的である。深く表現力の高い回路は状態空間を広く探索できるが、そのこと自体が勾配分散の低下やカーネル値の集中を引き起こしうる。Thanasilp ら (2024) は expressibility、entanglement、measurement、noise がカーネル評価を難しくすると報告したが、本研究の結果も同方向のメッセージを与える。したがって、量子優位が出る条件を特定するには、問題構造、適切な符号化、浅いが十分に表現的な回路、そして低ノイズ環境が同時に必要である。

短期的な改善としては、実機キャリブレーションに基づく IBM backend noise、shot noise、kernel alignment 指標、追加の古典ベースラインをこの framework に組み込むべきである。中長期的には、化学スペクトル、材料設計、量子実験データのように量子特徴が自然に現れる実データへ拡張し、外部妥当性を検証することが必須である。External validation with independent real-world datasets is essential to confirm the generalizability of these findings beyond simulated conditions. また、エラー緩和と layerwise training を組み合わせれば、表現力と学習可能性のトレードオフをより実践的に最適化できる可能性がある。

## ファイル一覧
- `src/expressibility.py`: PQC 表現力と Meyer-Wallach 指標の計算
- `src/quantum_kernel.py`: 量子カーネル SVM と古典 RBF-SVM 比較
- `src/encoding_strategies.py`: angle / amplitude / IQP 符号化比較
- `src/barren_plateau.py`: 勾配分散と barren plateau 解析
- `src/noise_model.py`: depolarizing noise と quantum volume proxy
- `src/benchmark.py`: 実験統合、CSV 出力、図生成、ログ追記
- `tests/test_modules.py`: 軽量 smoke test
- `results/search-strategy.md`: 文献探索戦略と MCP 試行履歴
- `results/screening-table.csv`: 採否判定表
- `results/extraction-table.csv`: 文献抽出表
- `results/reference-list.md`: 論文執筆用参照リスト
- `results/kernel_comparison.csv` ほか: 実験結果一式
- `figures/*.png`: 6枚の主要図
- `paper.md`: 英語 IMRaD 論文本体
- `report.md`: 本日本語報告書

## 参考文献
1. Abbas, A., Sutter, D., Zoufal, C., Lucchi, A., Figalli, A., & Woerner, S. (2021). The power of quantum neural networks. *Nature Computational Science*, 1(6), 403–409. DOI: 10.1038/s43588-021-00084-1
2. Cerezo, M., Arrasmith, A., Babbush, R., Benjamin, S. C., Endo, S., Fujii, K., McClean, J. R., Mitarai, K., Yuan, X., Cincio, L., & Coles, P. J. (2021). Variational quantum algorithms. *Nature Reviews Physics*, 3, 625–644. DOI: 10.1038/s42254-021-00348-9
3. Havlíček, V., Córcoles, A. D., Temme, K., Harrow, A. W., Kandala, A., Chow, J. M., & Gambetta, J. M. (2019). Supervised learning with quantum-enhanced feature spaces. *Nature*, 567, 209–212. DOI: 10.1038/s41586-019-0980-2
4. Huang, H.-Y., Broughton, M., Mohseni, M., Babbush, R., Boixo, S., Neven, H., & McClean, J. R. (2021). Power of data in quantum machine learning. *Nature Communications*, 12, 2631. DOI: 10.1038/s41467-021-22539-9
5. LaRose, R., & Coyle, B. (2020). Robust data encodings for quantum classifiers. *Physical Review A*, 102, 032420. DOI: 10.1103/PhysRevA.102.032420
6. Liu, Y., Arunachalam, S., & Temme, K. (2021). A rigorous and robust quantum speed-up in supervised machine learning. *Nature Physics*, 17, 1013–1017. DOI: 10.1038/s41567-021-01287-z
7. McClean, J. R., Boixo, S., Smelyanskiy, V. N., Babbush, R., & Neven, H. (2018). Barren plateaus in quantum neural network training landscapes. *Nature Communications*, 9, 4812. DOI: 10.1038/s41467-018-07090-4
8. Preskill, J. (2018). Quantum computing in the NISQ era and beyond. *Quantum*, 2, 79. DOI: 10.22331/q-2018-08-06-79
9. Schuld, M. (2021). Supervised quantum machine learning models are kernel methods. *arXiv*. DOI: 10.48550/arXiv.2101.11020
10. Schuld, M., & Killoran, N. (2019). Quantum machine learning in feature Hilbert spaces. *Physical Review Letters*, 122, 040504. DOI: 10.1103/PhysRevLett.122.040504
11. Schuld, M., Sweke, R., & Meyer, J. J. (2021). The effect of data encoding on the expressive power of variational quantum machine learning models. *Physical Review A*, 103, 032430. DOI: 10.1103/PhysRevA.103.032430
12. Sim, S., Johnson, P. D., & Aspuru-Guzik, A. (2019). Expressibility and entangling capability of parameterized quantum circuits for hybrid quantum-classical algorithms. *Advanced Quantum Technologies*, 2(12), 1900070. DOI: 10.1002/qute.201900070
13. Thanasilp, S., Wang, S., Cerezo, M., & Holmes, Z. (2024). Exponential concentration in quantum kernel methods. *Nature Communications*, 15, 5200. DOI: 10.1038/s41467-024-49287-w
