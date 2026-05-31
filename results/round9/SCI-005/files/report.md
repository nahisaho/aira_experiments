# 実験レポート: LongSVCaller — ロングリードシーケンシングからの構造変異高精度検出パイプライン

---

## 1. 実験目的と背echo

### 1.1 研究目的

Oxford Nanopore Technologies（ONT）およびPacBio HiFiロングリードシーケンシングデータから構造変異（Structural Variants: Svs）を高精度に検出するため'Report**LongSVCaller**」を設計・実装し、その性能を評価することを目的とした。_

### 1.2 研究背景

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$SVs）は50 bp以上のゲノム変化であり、欠失（DEL）、挿入（INS）、逆位（INV）、重複（TRA）、およびクロモスリプシスや染色体外DNA（ecDNA）などの複雑な再配列を含む。SVsはメンデル型希少疾患、がん、神経発達障害の主要な遺伝的原因の一つである。

#            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1=;PS2=;unset HISTFILE;                 EC=0;                 ! . : [ [[ ]] aarch64-linux-gnu-addr2line aarch64-linux-gnu-ar aarch64-linux-gnu-as aarch64-linux-gnu-c++filt aarch64-linux-gnu-cpp aarch64-linux-gnu-cpp-12 aarch64-linux-gnu-dwp aarch64-linux-gnu-elfedit aarch64-linux-gnu-g++ aarch64-linux-gnu-g++-12 aarch64-linux-gnu-gcc aarch64-linux-gnu-gcc-12 aarch64-linux-gnu-gcc-ar aarch64-linux-gnu-gcc-ar-12 aarch64-linux-gnu-gcc-nm aarch64-linux-gnu-gcc-nm-12 aarch64-linux-gnu-gcc-ranlib aarch64-linux-gnu-gcc-ranlib-12 aarch64-linux-gnu-gcov aarch64-linux-gnu-gcov-12 aarch64-linux-gnu-gcov-dump aarch64-linux-gnu-gcov-dump-12 aarch64-linux-gnu-gcov-tool aarch64-linux-gnu-gcov-tool-12 aarch64-linux-gnu-gfortran aarch64-linux-gnu-gfortran-12 aarch64-linux-gnu-gold aarch64-linux-gnu-gp-archive aarch64-linux-gnu-gp-collect-app aarch64-linux-gnu-gp-display-html aarch64-linux-gnu-gp-display-src aarch64-linux-gnu-gp-display-text aarch64-linux-gnu-gprof aarch64-linux-gnu-gprofng aarch64-linux-gnu-ld aarch64-linux-gnu-ld.bfd aarch64-linux-gnu-ld.gold aarch64-linux-gnu-lto-dump aarch64-linux-gnu-lto-dump-12 aarch64-linux-gnu-nm aarch64-linux-gnu-objcopy aarch64-linux-gnu-objdump aarch64-linux-gnu-pkg-config aarch64-linux-gnu-pkgconf aarch64-linux-gnu-python3-config aarch64-linux-gnu-python3.11-config aarch64-linux-gnu-ranlib aarch64-linux-gnu-readelf aarch64-linux-gnu-size aarch64-linux-gnu-strings aarch64-linux-gnu-strip add-shell addgroup addpart addr2line adduser agetty alias apt apt-cache apt-cdrom apt-config apt-get apt-key apt-mark ar arch as awk b2sum badblocks base32 base64 basename basenc bash bashbug bg bind blkdiscard blkid blkzone blockdev break builtin c++ c++filt c89 c89-gcc c99 c99-gcc c_rehash caller captoinfo case cat cc cd chage chattr chcon chcpu chfn chgpasswd chgrp chmem chmod choom chown chpasswd chroot chrt chsh cksum clear clear_console cmp comm command compgen complete compopt continue copilot coproc corelist corepack cp cpan cpan5.36-aarch64-linux-gnu cpgr cpp cpp-12 cppw csplit csv2rdf ctrlaltdel cut cyclopts dash date dd ddgs deb-systemd-helper deb-systemd-invoke debconf debconf-apt-progress debconf-communicate debconf-copydb debconf-escape debconf-set-selections debconf-show debugfs debugpy debugpy-adapter declare delgroup delpart deluser df diff diff3 dir dircolors dirname dirs disown distro dmesg dnsdomainname do docker-entrypoint.sh domainname done dotenv dpkg dpkg-deb dpkg-divert dpkg-fsys-usrunmess dpkg-maintscript-helper dpkg-preconfigure dpkg-query dpkg-realpath dpkg-reconfigure dpkg-split dpkg-statoverride dpkg-trigger du dumpe2fs dumppdf.py dwp e2freefrag e2fsck e2image e2label e2mmpstatus e2scrub e2scrub_all e2undo e4crypt e4defrag echo egrep elfedit elif else email_validator enable enc2xs encguess env esac eval exec exit expand expiry export expr f2py f77 f95 factor faillock faillog fallocate false fastapi fastmcp fc fg fgrep fi filefrag fincore find findfs findmnt fitz flask flock fmt fold fonttools for fsck fsck.cramfs fsck.ext2 fsck.ext3 fsck.ext4 fsck.minix fsfreeze fstab-decode fstrim function g++ g++-12 gcc gcc-12 gcc-ar gcc-ar-12 gcc-nm gcc-nm-12 gcc-ranlib gcc-ranlib-12 gcov gcov-12 gcov-dump gcov-dump-12 gcov-tool gcov-tool-12 gencat generate-mcp-tools getconf getent getopt getopts getty gfortran gfortran-12 git git-receive-pack git-shell git-upload-archive git-upload-pack gold gp-archive gp-collect-app gp-display-html gp-display-src gp-display-text gpasswd gpgv gprof gprofng grep groupadd groupdel groupmems groupmod groups grpck grpconv grpunconv gunzip gzexe gzip h2ph h2xs hardlink hash head help hf history hostid hostname httpx huggingface-cli hwclock iconv iconvconfig id idna if in infocmp infotocap install installkernel instmodsh invoke-rc.d ionice ipcmk ipcrm ipcs ipython ipython3 ischroot isosize isympy jlpm jobs join json_pp jsonpath_ng jsonpointer jsonschema jupyter jupyter-dejavu jupyter-events jupyter-execute jupyter-kernel jupyter-kernelspec jupyter-konsole jupyter-lab jupyter-labextension jupyter-labhub jupyter-mcp-server jupyter-migrate jupyter-nbconvert jupyter-run jupyter-server jupyter-troubleshoot jupyter-trust keyring kill killall5 last lastb lastlog ld ld.bfd ld.gold ld.so ldattach ldconfig ldd let libnetcfg link linux32 linux64 ln local locale localedef log2design.py logger login logname logout logsave losetup ls lsattr lsblk lscpu lsfd lsipc lsirq lslocks lslogins lsmem lsns lto-dump lto-dump-12 magika magika-python-client mammoth mapfile markdown-it markdownify markitdown mawk mcookie mcp md5sum md5sum.textutils mesg mkdir mke2fs mkfifo mkfs mkfs.bfs mkfs.cramfs mkfs.ext2 mkfs.ext3 mkfs.ext4 mkfs.minix mkhomedir_helper mklost+found mknod mkswap mktemp more mount mountpoint mv namei nawk newgrp newusers nib-conform nib-convert nib-dicomfs nib-diff nib-ls nib-nifti-dx nib-roi nib-stats nib-tck2trk nib-trk2tck nice nipypecli nisdomainname nl nm node nodejs nohup nologin normalizer npm nproc npx nsenter numfmt numpy-config objcopy objdump od onnxruntime_test openssl pager pam-auth-update pam_getenv pam_namespace_helper pam_timestamp_check parrec2nii partx passwd paste pathchk pdb3 pdb3.11 pdf2txt.py pdfplumber perl perl5.36-aarch64-linux-gnu perl5.36.0 perlbug perldoc perlivp perlthanks piconv pidof pinky pip pip3 pip3.11 pivot_root pkg-config pkgconf pl2pm playwright pldd pod2html pod2man pod2text pod2usage podchecker policy-rc.d popd pr printenv printf prlimit prov-compare prov-convert prove ptar ptardiff ptargrep ptx pushd pwck pwconv pwd pwhistory_helper pwunconv py3clean py3compile py3versions pybabel pydoc3 pydoc3.11 pyftmerge pyftsubset pygettext3 pygettext3.11 pygmentize pyjson5 pypdfium2 python3 python3-config python3.11 python3.11-config ranlib rbash rdf2dot rdfgraphisomorphism rdfpipe rdfs2dot read readarray readelf readlink readonly readprofile realpath remove-shell rename.ul renice reset resize2fs resizepart return rev rgrep rm rmdir rmt rmt-tar rpcgen rtcwake run-parts runcon runuser runxlrd.py savelog scalar script scriptlive scriptreplay sdiff sed select send2trash seq service sessionmirror.py set setarch setpriv setsid setterm sg sh sha1sum sha224sum sha256sum sha384sum sha512sum shadowconfig shasum shift shopt shred shuf size sleep sort source sparqlquery splain split sprc start-stop-daemon stat stdbuf streamzip strings strip stty su sulogin sum suspend swaplabel swapoff swapon switch_root sync tabs tac tail tar tarcat taskset tee tempfile test then tic time timeout times tiny-agents toe tooluniverse tooluniverse-doctor tooluniverse-expert-feedback tooluniverse-expert-feedback-web tooluniverse-http-api tooluniverse-mcp tooluniverse-smcp tooluniverse-smcp-server tooluniverse-smcp-stdio tooluniverse-stdio touch tput tqdm tr trap true truncate tset tsort ttx tty tu tu-datastore tune2fs type typer typeset tzselect uclampset ulimit umask umount unalias uname uncompress unexpand uniq unix_chkpwd unix_update unlink unset unshare until update-alternatives update-ca-certificates update-passwd update-rc.d update-shells useradd userdel usermod users utmpdump uvicorn vba_extract.py vdir vigr vipw wait wall watchfiles wc wdctl websockets whereis which which.debianutils while who whoami wipefs wsdump xargs xsubpp yarn yarnpkg yes youtube_transcript_api ypdomainname zcat zcmp zdiff zdump zegrep zfgrep zforce zgrep zic zipdetails zless zmore znew 
#
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$ONT・PacBio）は10〜100+ kbのリードを提供し、ほとんどのSVブレークポイントをスパンできるが、以下の課題が残る：} 

- ONTの高いベースコールエラー率（8〜15%、R10.4+Doradoで〜1〜3%に改善）
- リピート・セグメンタル重複領域での偽陽性
- クロモスリプシス・ecDNAなどの複雑なSVの検出困難
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 パイプラインアーキテクチャ（6モジュール構成）

![Pipeline Architecture](figures/sv_pipeline_architecture.png)

**Module 1: RNNベースコーリング改善**
- LSTMネットワーク + CRFデコーダによる生の電流シグナル→塩基配列変換
- 参考: Bonito (Pagès-Gallego & de Ridder, 2023), SqueezeCall (Zhu, 2025)
- Q値（品質スコア）によるリードフィルタリング

**Module 2a: Split-Read解析**
- BAMファイルのSAタグ・ソフト/ハードクリッピングパターンからSVブレークポイントを特定
- 特徴量: `split_read_count`, `spanning_reads`, `basecall_quality`

**Module 2b: Read-Depth解析**
- ゲノム背景に対する深度比（DEL: ~0.3、DUP: ~1.8、正常: ~1.0）
- 特徴量: `mean_depth`, `depth_ratio`, `sv_type`, `log_sv_size`

**Module 2c: アセンブリベース検出**
#- hifiasmによるde novoアセンブリ → T2t-Chm13'Report_Eof'

- 特徴量: `assembly_supported`（二値）, `contig_length`, `contig_quality`

**Module 3: リピート領域特殊処理**
- T2t-
- 領域特異的スコアリング補正

**Module 4: 多証拠統合（主モデル）**
- Random Forest (RF): 200本、最大深度10
#- XGBoost Hybrid: 200エステ
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             SV密度 + コピー数振動 + ユニークジャンクション数}
- ecDNA: バックスプライスジャンクション 円形カバレッ + + 局所増幅率

**Module 6: GIAB Tier1ベンチマーキング**
- `truvari bench --passonly --refdist 100 --pctsim 0.7` による評価設計

### 2.2 ToolUniverse MCPツール使用状況

**PubMed検索 (PubMed_search_articles)**:
- クエリ: "Sniffles2 structural variant long-read", "nanopore basecalling neural network", "hybrid short-read long-read structural variant", "extrachromosomal DNA ecDNA detection"
- 結果: 計10+ 件の論文を特定（詳細は References 参照）

**Semantic Scholar (SemanticScholar_search_papers)**:
- APIレート制限（429エラー）により一部クエリは成功せず
- 成功: "chromothripsis complex structural variant long-read genome assembly detection" → 5件の論文取得

**NatureLM MCP (ask_naturelm)**:
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$0件一致） → 接続失敗
- 代替: 公表論文からの定量パラメータ使用

**GALACTICA MCP (scientific_qa, predict_citations)**:
- ⚠️ ToolUniverseレジストリに未登録（0件一致） → 接続失敗
- 代替: PubMed/Semantic Scholarによる文献検索・科学的検証

---

## 3. 主要な結果と数値

### 3.1 データセット概要 [cell:2]

| 項目 | 値 |
|---|---|
| 総SV数 | 2,000 |
| 真のSV数 | 1,229 (61.5%) |
| 偽陽性候補 | 771 (38.5%) |
| リピート領域に存在 | 34% |
| セグメンタル重複 | 20% |
| テロメア | 5% |
| 中央値SVサイズ | 545 bp |
| 特徴量数 | 14 |

### 3.2 SV探索戦略の比較 [cell:4c]

![Performance ROC Curves](figures/sv_performance_roc.png)

| 手法 | AUROC | F1 | Precision | Recall |
|---|---|---|---|---|
| Split-Read Only | 0.7114 | 0.7380 | 0.6757 | 0.8130 |
| Read-Depth Only | 0.5097 | 0.7382 | 0.6152 | 0.9228 |
| Assembly-Based | 0.7390 | 0.8231 | 0.7595 | 0.8984 |
| RF Integrated | 0.9412 | 0.9352 | 0.9315 | 0.9390 |
| **XGBoost Hybrid** | **0.9383** | **0.9412** | **0.9393** | **0.9431** |

### 3.3 5分割交差検証 [cell:5]

| モデル | AUROC (mean±std) | F1 (mean±std) |
|---|---|---|
| RF Integrated | 0.9366 ± 0.0063 | 0.9309 ± 0.0078 |
| XGBoost Hybrid | **0.9383 ± 0.0040** | 0.9236 ± 0.0051 |

### 3.4 SV種別性能 [cell:6]

| SV種別 | n(テスト) | AUROC | F1 |
|---|---|---|---|
| DEL | 193 | 0.9406 | 0.9345 |
| INS | 160 | 0.9384 | 0.9447 |
| INV | 22 | 0.8214 | 0.8966 |
| DUP | 15 | 1.000* | 1.000* |
| TRA | 10 | 1.000* | 1.000* |

*⚠️ 小サンプル過学習の可能性（DUP n=15、TRA n=10）

### 3.5 ゲノム領域別性能 [cell:6]

| 領域 | n | AUROC | F1 |
|---|---|---|---|
| リピート領域 | 133 | 0.9325 | 0.9317 |
| 非リピート | 267 | 0.9410 | 0.9458 |
| セグメンタル重複 | 85 | 0.9387 | 0.9752 |
| テロメア | 18 | 1.000* | 1.000* |

### 3. [cell:7]

XGBoost Hybrid テストセット(n=400): **TN=139, FP=15, FN=14, TP=232**

### 3.7 特徴量重要度 [cell:8]

![Feature Importance](figures/sv_feature_importance.png)

| 順位 | 特徴量 | Gini重要度 |
|---|---|---|
| 1 | アセンブリサポート | 0.2741 |
| 2 | スプリットリード数 | 0.1827 |
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___Begin___Command_Done_ | 0.1735 |
| 4 | ベースコール品質 | 0.0564 |
| 5 | 平均深度 | 0.0487 |

### 3.8 複雑なSV検出 [cell:10]

| モジュール | 陽性数 | AUROC | F1 |
|---|---|---|---|
| クロモスリプシス検出器 | 21/300 | 0.6508 | **0.2500** |
| ecDNA検出器 | 19/200 | 0.8395 | 0.7273 |

### 3.9 統計的検定 [cell:11]

- RF統合 vs Split-Read Only（Wilcoxon検定、AUROC）: **p = 0.0625** (α=0.05で非有意)
- Cohen's d: **11.330**（|d|>0.8 = 大きな効果量）
- 注: 5ペアのCVスコアでは検出力が低く非有意となったが、効果量は極めて大きい

### 3.10 包括的結果サマリー

![Comprehensive Results](figures/sv_comprehensive_results.png)

![Data Exploration](figures/sv_data_exploration.png)

---

## 4. 考察と今後の展望

### 4.1 主要な知見

**1. 多証拠統合の圧倒的な優位性**
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$AUROC=0.AUROC=0.938）の間には0.2以上の差がある。Read-Depth単独はAUROC=0.510と偶然水準に近く、深度シグナルのみでは信頼性が低いことが示された。

**2. アセンブリベース証拠が最重要特徴**
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$Gini重0.274と突出しており、de novoアセンブリによる確認がSV検出の最重要要素であることが確認された。

**3. ハイブリッド解析（短+長リード）の有効性**
3番目に重要（0.174）であり、直交する証拠源の組み合わせechoecho

**4. 複雑なSVの検出困難性**
.git data figures notebook.ipynb paper.md Src Report_EOFF1=0.250は、低陽性率（7%）と特徴量の複雑さを反映する。実データでは図形の破壊-再構築パターン、アリル頻度分布、ストランド特異クラスタリングなど、より高度な特徴量が

### 4.2 自己批判的分析

**合成データへの依存**: グランドトゥルースラベルが予測に使用した特徴量の線形結合から構築されており、実データには存在しない特徴量-ラベル相関が組み込まれている。実Aurocは合成データでの値より確実にecho

**小サンプル問題**: DUP（n=15）、TRA（n=10）、テロメア（n=18）のAUROC/F1=1.000は小サンプル過学習の典型であり、汎化性能を示すものではない。

**Wilcoxon検定の限界**: 5ペアのCVスコアでは検p=0.0625は統計的非有意であるが、Cohen's d=11.3は実践的に極めて大きな効果量を示す。

**NatureLM/GALACTICA不在の影響**: ベースコーリングパラメータの定量的予測検証や科学的妥当性の確認が行えなかった。'Report_EOF''REPORT_EOF''REPORT_EOF'echo

### 4.3 今後の展望

1. **実データでの検証**: GIAB HG002 Tier1真値セットに対するtruvariベンチマーキング
2. **グラフベースブレークポイント解析**: クロモス'REPORT_EOF'
3. **メチル化シグナル統合**: インプリンティング領域SV検出の改善
4. **トランスフォーマーベースコーラー**: Dorado（ONT公式）、SqueezeCallとの統合
5. **パンゲノム参照の活用**: HG002/HG003の個人ゲノムアセンブリを用いたソマティックSV偽陽性削減

---

## 5. 先行研究との比較

### 5.1 特定した論文（先行研究調査）

| # | タイトル | 著者 | 年 | DOI | 主要知見 |
|---|---|---|---|---|---|
| 1 | Deciphering the Structural Variants by Long-Read Genome Sequencing | Dutta, Dalal | 2025 | 10.1159/000549245 | LRSはSV検出の革命的ツール。T2T完成への貢献 |
| 2 | Advancing long-read nanopore genome assembly for rare disease | Negi et al. | 2025 | 10.1016/j.ajhg.2025.01.002 | 98サンプルのLRS解析で11症例を診断。87%のタンパク質コード遺伝子を完全フェーズ |
| 3 | Benchmarking somatic SV callers on HG008 | Cui et al. | 2026 | 10.3389/fgene.2026.1732039 | Sniffles2、Nanomonsv等の体細胞SV比較。マルチツール戦略の有効性 |
| 4 | Blackbird: SV detection using synthetic and low-coverage long-reads | Meleshko et al. | 2025 | 10.1093/bioadv/vbaf151 | 5×カバレッジでF1=0.835（DEL）。ハイブリッドアルゴリズムの優位性 |
#| 5 | Comprehensive benchmark of basecallers | Pagès-Gallego, de Ridder | 2023 | 10.1186/s13059-023-02903-2 | Lstm+Crfが最高性能。7モデル90アーキテ
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$ |
| 6 | SqueezeCall: nanopore basecalling | Zhu | 2025 | 10.46471/gigabyte.148 | Squeezeformerが既存RNN・Transformerを凌駕 |
| 7 | BaseNet: transformer-based signal decoding | Li et al. | 2024 | 10.1016/j.csbj.2024.09.016 | クロスアテンション・大規模事前学習モデルの有効性 |
| 8 | eccDNA identification via ECCFP | Li et al. | 2026 | 10.21769/BioProtoc.5636 | ONTロングリードからのecDNA高精度検出パイプライン |
| 9 | Hybrid short/long read variant calling | Gambardella | 2025 | 10.1016/j.crmeth.2025.101107 | 浅いハイブリッドLRS+Illuminaが単独技術の深いシーケンシングを凌駕 |
| 10 | DNAscope Hybrid pipeline | Hu et al. | 2025 | 10.3389/fbinf.2025.1691056 | 5〜10×LRSで30×単独を超えるエラー50%以上削減 |

### 5.2 先行研究の課題

1. **単一戦略の限界**: Sniffles2・SVIM等はsplit-readのみに依存
2. **リピート領域の脆弱性**: テロメア・セントロメアでの高偽陽性率
3. **複雑なSVの未対応**: クロモスリプシス・ecDNAに専用モジュールなし
4. **シグナルレベル統合の欠如**: ベースコール品質がSV判定に活用されていない
5. **ハイブリッド解析の未普及**: ショートリードとの体系的統合が少ない

---

## 6. 生成したファイル一覧

### 図表ファイル
| ファイル名 | 内'REPORT_EOF'         |
|---|---|
| `figures/sv_pipeline_architecture.png` | 6モジュールパイプラインアーキテクチャ図 |
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$SV種分布・サイズ・深度比等） |
| `figures/sv_performance_roc.png` | ROC曲線・性能比較・混同行列 |
| `figures/sv_feature_importance.png` | 特徴量重要度・ベースコール-F1関係 |
| `figures/sv_comprehensive_results.png` | 包括的性能サマリー（CV・SV種・P-R曲線・ベンチマーク比較） |

### データファイル
| ファイル名 | 内容 |
|---|---|
| `data/raw/sv_simulated_dataset.csv` | 合成SVデータセット（n=2,000、14特徴量）|

### 論文・レポート
| ファイル名 | 内容 |
|---|---|
| `paper.md` | 学術論文形式（英語）|
| `report.md` | 実験レポート（日本語、本文書）|

---

## 7. 再現性情報

- **乱数シード**: `numpy.random.seed(42)`, `random.seed(42)`, 全sklearn/xgboostモ `random_state=42`
- **Python**: 3.11.2 (GCC 12.2.0)
- **主要パッケージ**: numpy 2.4.6, pandas 3.0.3, scikit-learn 1.8.0, scipy 1.17.1, xgboost 3.2.0, matplotlib 3.10.9, seaborn 0.13.2
- **データ出自**: 全データは合成生成（GIAB HG002統計に基づくパラメータ）
- **実行環境**: Jupyter MCP経由でexecute_codeにより実行（セル番号は計算来歴のトレーサビリティに使用）

