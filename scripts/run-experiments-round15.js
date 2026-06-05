/**
 * AIRA Experiment Runner - Round 15
 * AIRA v3.4.8 — Fixed benchmark subset (n=30).
 * Selected 30 representative experiments from Round 8-12 analysis:
 *   - Domain diversity: general-science=14, genomics=6, molecular=4, protein=3, materials=3
 *   - High-variance experiments prioritized (gate changes, uncited claims, repair needs)
 *   - Enables direct cross-version comparison with reduced variance
 * Default 4 parallel workers (optimal stability/throughput from Round 12 analysis).
 * After each experiment, runs POST /validate; if gates fail, runs
 * POST /validate/repair → sends repair prompt → re-validates (max 3 iterations).
 * On 3-iteration failure, calls POST /validate/postmortem.
 * No timeout - waits until completion.
 * Uses direct HTTP API calls + WebSocket for chat trigger.
 * Target: http://192.168.1.15:3000 (AIRA)
 */

const fs = require('fs');
const path = require('path');
const http = require('http');
const WebSocket = require('ws');

const BASE_URL = process.env.AIRA_URL || 'http://192.168.1.15:3000';
const RESULTS_DIR = process.env.RESULTS_DIR_OVERRIDE || path.join(__dirname, '..', 'results', 'round15');
const ROUND9_DIR = path.join(__dirname, '..', '..', 'aira', 'co-scientist-optimization', 'round-9');
const TIMEOUT_MS = 24 * 60 * 60 * 1000; // 24 hours (effectively no timeout)
const POLL_INTERVAL_MS = 15000;
const NUM_WORKERS = parseInt(process.env.WORKERS || '4', 10);
const MAX_REPAIR_ITERATIONS = 3;

// Fixed benchmark subset (30 experiments selected from Round 8-12 analysis)
// Domain: general-science=14, genomics=6, molecular=4, protein=3, materials=3
const FIXED_SUBSET = new Set([
  'SCI-001','SCI-002','SCI-003','SCI-004','SCI-005','SCI-006','SCI-007',
  'SCI-011','SCI-015','SCI-017','SCI-018','SCI-020','SCI-021','SCI-022',
  'SCI-028','SCI-029','SCI-035','SCI-049','SCI-051','SCI-058','SCI-059',
  'SCI-060','SCI-065','SCI-073','SCI-082','SCI-083','SCI-084','SCI-090',
  'SCI-092','SCI-095'
]);

const CO_SCIENTIST_SKILL = '3776b03b-1800-4640-9f4a-8d5c1d418fa1';

// Global error handlers to prevent silent crashes
process.on('unhandledRejection', (reason, promise) => {
  console.error(`[UNHANDLED REJECTION] ${reason}`);
  console.error(reason?.stack || '');
});
process.on('uncaughtException', (err) => {
  console.error(`[UNCAUGHT EXCEPTION] ${err.message}`);
  console.error(err.stack || '');
});

function sleep(ms) { return new Promise(r => setTimeout(r, ms)); }

function log(workerId, msg) {
  const ts = new Date().toISOString().slice(11, 19);
  const line = `[${ts}][W${workerId}] ${msg}\n`;
  process.stdout.write(line);
}

// Direct HTTP helper (no CORS) with timeout
function httpRequest(method, urlPath, body, timeoutMs = 30000) {
  return new Promise((resolve, reject) => {
    const url = new URL(urlPath, BASE_URL);
    const options = {
      hostname: url.hostname,
      port: url.port,
      path: url.pathname + url.search,
      method,
      headers: {},
      timeout: timeoutMs
    };
    if (body) {
      const data = JSON.stringify(body);
      options.headers['Content-Type'] = 'application/json';
      options.headers['Content-Length'] = Buffer.byteLength(data);
    }

    const req = http.request(options, (res) => {
      let responseData = '';
      res.on('data', chunk => responseData += chunk);
      res.on('end', () => {
        try { resolve(JSON.parse(responseData)); }
        catch { resolve(responseData); }
      });
    });
    req.on('timeout', () => { req.destroy(); reject(new Error('HTTP request timeout')); });
    req.on('error', reject);
    if (body) req.write(JSON.stringify(body));
    req.end();
  });
}

async function apiGet(apiPath) {
  return httpRequest('GET', `/api${apiPath}`);
}

async function apiPost(apiPath, body) {
  return httpRequest('POST', `/api${apiPath}`, body);
}

async function apiDelete(apiPath) {
  return httpRequest('DELETE', `/api${apiPath}`);
}

// Get CSRF token
async function getCsrfToken() {
  const result = await apiGet('/csrf-token');
  return result.token;
}

// POST with CSRF
async function apiPostWithCsrf(apiPath, body) {
  const token = await getCsrfToken();
  return new Promise((resolve, reject) => {
    const url = new URL(`/api${apiPath}`, BASE_URL);
    const data = body ? JSON.stringify(body) : '';
    const options = {
      hostname: url.hostname,
      port: url.port,
      path: url.pathname + url.search,
      method: 'POST',
      timeout: 60000,
      headers: {
        'X-AIRA-Token': token,
        'Origin': `http://localhost:3000`,
      }
    };
    if (body) {
      options.headers['Content-Type'] = 'application/json';
      options.headers['Content-Length'] = Buffer.byteLength(data);
    }

    const req = http.request(options, (res) => {
      let responseData = '';
      res.on('data', chunk => responseData += chunk);
      res.on('end', () => {
        try { resolve(JSON.parse(responseData)); }
        catch { resolve(responseData); }
      });
    });
    req.on('timeout', () => { req.destroy(); reject(new Error('POST CSRF request timeout')); });
    req.on('error', reject);
    if (body) req.write(data);
    req.end();
  });
}

async function apiDeleteWithCsrf(apiPath) {
  const token = await getCsrfToken();
  return new Promise((resolve, reject) => {
    const url = new URL(`/api${apiPath}`, BASE_URL);
    const options = {
      hostname: url.hostname,
      port: url.port,
      path: url.pathname + url.search,
      method: 'DELETE',
      timeout: 30000,
      headers: {
        'X-AIRA-Token': token,
        'Origin': `http://localhost:3000`,
      }
    };

    const req = http.request(options, (res) => {
      let responseData = '';
      res.on('data', chunk => responseData += chunk);
      res.on('end', () => {
        try { resolve(JSON.parse(responseData)); }
        catch { resolve(responseData); }
      });
    });
    req.on('timeout', () => { req.destroy(); reject(new Error('DELETE CSRF request timeout')); });
    req.on('error', reject);
    req.end();
  });
}

/**
 * Open WebSocket to project. Returns WS instance that must stay open during run.
 * Includes automatic ping to keep connection alive.
 */
function openProjectWs(projectId, retries = 3) {
  return new Promise((resolve, reject) => {
    let attempt = 0;
    function tryConnect() {
      attempt++;
      const wsUrl = BASE_URL.replace('http://', 'ws://') + `/ws/projects/${projectId}/chat`;
      const ws = new WebSocket(wsUrl, {
        headers: { 'Origin': 'http://localhost:3000' }
      });
      let pingInterval;
      ws.on('open', () => {
        pingInterval = setInterval(() => {
          if (ws.readyState === WebSocket.OPEN) {
            ws.ping();
          }
        }, 30000);
        resolve(ws);
      });
      ws.on('close', () => { clearInterval(pingInterval); });
      ws.on('error', (err) => {
        clearInterval(pingInterval);
        if (attempt < retries) {
          setTimeout(tryConnect, 3000);
        } else {
          reject(err);
        }
      });
      setTimeout(() => {
        if (ws.readyState !== WebSocket.OPEN) {
          ws.terminate();
          if (attempt < retries) {
            setTimeout(tryConnect, 3000);
          } else {
            reject(new Error('WebSocket connection timeout'));
          }
        }
      }, 10000);
    }
    tryConnect();
  });
}

function sendChat(ws, content) {
  ws.send(JSON.stringify({ type: 'chat', content }));
}

/**
 * Classify experiment into domain category for targeted GALACTICA instructions.
 */
function classifyDomain(prompt) {
  const molecular = ['SMILES', '分子', '低分子', '薬物', '創薬', '化合物', 'ADMET', '溶解度', 'logP',
    '結合親和性', 'binding', 'レトロ合成', '合成経路', '有機合成', 'PROTAC', 'ADC', 'ペイロード', 'リンカー',
    'drug', 'molecule', 'リガンド'];
  const protein = ['タンパク質', 'タンパク', 'ペプチド', '抗体', 'protein', 'mRNA', 'ワクチン'];
  const materials = ['材料', '合金', 'ペロブスカイト', 'ポリマー', '触媒', '電解質', 'MOF', '金属有機',
    'CO2還元', 'リチウム', '電池', 'material', '生分解性', 'コポリマー', '自己組織化'];
  const genomics = ['CRISPR', 'ゲノム', 'RNA', 'DNA', 'エピジェネティ', 'メチル化',
    'ファーマコゲノミクス', '薬物応答', '代謝モデル', '代謝物'];

  if (molecular.some(k => prompt.includes(k))) return 'molecular';
  if (protein.some(k => prompt.includes(k))) return 'protein';
  if (materials.some(k => prompt.includes(k))) return 'materials';
  if (genomics.some(k => prompt.includes(k))) return 'genomics';

  const sciKeywords = ['モデル', '予測', 'シミュレーション', 'アルゴリズム', '最適化', 'ニューラル',
    '機械学習', 'ディープラーニング', '統計', 'パイプライン', 'フレームワーク'];
  if (sciKeywords.some(k => prompt.includes(k))) return 'general-science';
  return 'non-science';
}

function isGALACTICARelevant(prompt) {
  const domain = classifyDomain(prompt);
  return ['molecular', 'protein', 'materials', 'genomics'].includes(domain);
}

/**
 * Domain-specific NatureLM tool instructions (quantitative prediction)
 */
function getNatureLMInstructions(domain) {
  const base = `  - NatureLM MCP サーバーのツールを**直接呼び出す**こと（ToolUniverse を経由しない）`;
  switch (domain) {
    case 'molecular':
      return `${base}
  - \`generate_smiles\`: 目的の性質を持つ候補分子を複数生成し、探索空間を拡大する
  - \`predict_logp\`, \`predict_property\`: 生成した分子の物性を予測し、実験条件のベースラインを設定する
  - \`retrosynthesis\`: 合成可能性を検証し、実現可能な候補に絞り込む
  - \`ask_naturelm\`: 分子メカニズムに関する定量的パラメータ（結合エネルギー、IC50推定値、LogP等）を取得する`;
    case 'protein':
      return `${base}
  - \`generate_protein_sequence\`: 目的の機能を持つタンパク質配列を生成する
  - \`predict_property\`: タンパク質の物性予測を行い、設計パラメータを定量化する
  - \`ask_naturelm\`: 構造-活性相関、安定性条件、フォールディング特性に関する知見を取得する`;
    case 'materials':
      return `${base}
  - \`predict_material_composition\`: 目的の性質を持つ材料組成を予測する
  - \`predict_property\`: 候補材料の物性を予測し、スクリーニング基準を設定する
  - \`ask_naturelm\`: 材料の安定性、劣化メカニズム、界面特性に関する定量的知見を取得する`;
    case 'genomics':
      return `${base}
  - \`ask_naturelm\`: 生物学的メカニズムに関する定量的パラメータ（結合自由エネルギー、反応速度定数等）を取得する`;
    default:
      return `${base}
  - \`ask_naturelm\`: 研究テーマに関連する科学的知見・定量的パラメータを取得する`;
  }
}

/**
 * Domain-specific GALACTICA tool instructions (scientific validation & citations)
 */
function getGALACTICAInstructions(domain) {
  const base = `  - GALACTICA MCP サーバーのツールを**直接呼び出す**こと（ToolUniverse を経由しない）`;
  switch (domain) {
    case 'molecular':
      return `${base}
  - \`generate_molecule\`: 候補分子（SMILES）を複数生成し、NatureLMの予測と比較する
  - \`scientific_qa\`: 分子メカニズムの科学的妥当性を検証し、NatureLMの予測値の信頼性を評価する
  - \`predict_citations\`: 関連する先行研究の引用を予測し、文献調査を補完する
  - \`reasoning\`: 反応機構や合成経路の推論を行う`;
    case 'protein':
      return `${base}
  - \`predict_protein_annotations\`: 対象タンパク質のアミノ酸配列から機能予測・アノテーションを取得する
  - \`scientific_qa\`: タンパク質設計の科学的妥当性を検証する
  - \`predict_citations\`: 関連文献を予測し、文献レビューを補完する`;
    case 'materials':
      return `${base}
  - \`scientific_qa\`: 材料特性の科学的妥当性を検証し、NatureLMの予測と照合する
  - \`generate_molecule\`: 候補材料の分子構造を生成し、スクリーニング対象を拡大する
  - \`reasoning\`: 材料特性の物理的推論を行い、従来の探索範囲を超えた候補も検討する
  - \`generate_latex\`: 材料科学の数式（状態方程式、拡散方程式等）を生成する`;
    case 'genomics':
      return `${base}
  - \`scientific_qa\`: 生物学的メカニズムの科学的妥当性を検証する
  - \`predict_citations\`: ゲノミクス分野の関連文献を予測し、文献レビューを補完する`;
    default:
      return `${base}
  - \`scientific_qa\`: 研究テーマに関連する科学的知見を取得し、実験設計の妥当性を検証する
  - \`predict_citations\`: 関連文献を予測し、文献調査を補完する`;
  }
}

function getCombinedModelInstructions(domain) {
  const naturelmInstr = getNatureLMInstructions(domain);
  const galacticaInstr = getGALACTICAInstructions(domain);

  return `- **NatureLM MCP**（定量予測）を**直接呼び出して**以下のように活用すること（ToolUniverse を経由しない）：
${naturelmInstr}

- **GALACTICA MCP**（科学的検証・引用予測）を**直接呼び出して**以下のように活用すること（ToolUniverse を経由しない）：
${galacticaInstr}

- ⚠️ **両モデルの相互検証**: NatureLMの定量予測結果をGALACTICAのscientific_qaで科学的に検証し、矛盾がないか確認すること
- ⚠️ **直接呼び出しの徹底**: NatureLM/GALACTICA のツールは、それぞれの MCP サーバーから直接呼び出すこと。ToolUniverse MCP を経由した間接的な呼び出しは行わないこと
- 両モデルの予測結果はResultsセクションに定量的に記載し、一致・不一致を明示すること`;
}

/**
 * Generate round15 optimized prompt
 * Major optimizations from Round 12/13 analysis:
 *   - Golden rules at top for VM reduction
 *   - Citation ledger cell ([cell:results-summary]) workflow
 *   - Figure manifest cell ([cell:figure-manifest]) for orphan prevention
 *   - Streamlined steps (4→3 + pre-paper gate)
 *   - Direct NatureLM/GALACTICA calls (no ToolUniverse)
 *   - report.md downgraded to optional
 */
function generateRound15Prompt(originalPrompt) {
  const domain = classifyDomain(originalPrompt);
  const combinedInstructions = getCombinedModelInstructions(domain);

  const parts = originalPrompt.split('---\n');
  const body = parts.slice(1).join('---\n');

  const newPrefix = `以下の研究テーマについて、科学研究を実施し paper.md を作成してください。

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
🔴 **4つのゴールデンルール**（最優先 — 全ステップで遵守）
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

**Rule 1**: \`[cell:aira-env]\` と \`[cell:aira-seed]\` を**最初に実行**する。
**Rule 2**: paper.md の全ての定量的主張に \`[cell:<id>]\` 引用を付ける。引用する値は、そのセルの**最終出力**に含まれていなければならない。
**Rule 3**: paper.md を書く前に \`[cell:results-summary]\` を作成・実行し、全ての引用値を集約する（詳細はステップ3参照）。
**Rule 4**: paper.md を書く前に \`[cell:figure-manifest]\` を作成・実行し、全ての図の存在を確認する（詳細はステップ3参照）。

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

**ステップ1: 先行研究調査 + 実験計画**
- ToolUniverse MCP の学術検索ツール（Semantic Scholar / PubMed / Crossref 等）で先行研究を調査する
- 関連する最新論文（2020年以降）を5件以上特定し、タイトル・著者・DOI・主要知見をまとめる
- 先行研究の課題・限界を踏まえて実験計画を立てる
${combinedInstructions}
- ⚠️ NatureLM/GALACTICA 接続失敗時は、試行ツール名・エラー・代替手段を Methods に記録する

**ステップ2: Python実装・実行・検証（Jupyter MCP）**
- データ処理・分析・シミュレーションを Python で実装し、Jupyter MCP で実行する
- 推奨ライブラリ: \`pandas\`, \`numpy\`, \`scikit-learn\`, \`matplotlib\`, \`seaborn\`, \`scipy.stats\`, \`rdkit\`（分子系）
- \`random_state=42\` 等で乱数シードを固定する
- ⚠️ AUC/F1 等が 1.000 になったら過学習・データリークを疑い、交差検証の標準偏差付きで報告する
- ⚠️ コードで得た数値のみを論文に記載する（手計算・推測は不可）
- ⚠️ stderr エラーのあるセルの結果は引用しない — 修正して再実行すること
- ⚠️ 図は必ず \`plt.savefig('figures/xxx.png')\` で保存し、存在しない図を参照しない

**ステップ3: 引用台帳 + 図マニフェスト作成**（ゴールデンルール3）
paper.md 執筆前に、以下の2つのセルを作成・実行する：

\`[cell:results-summary]\` — 引用台帳セル：
\`\`\`python
# paper.md で引用する全ての数値をここにまとめる
results = {
    "metric_name": "exact_value_as_written_in_paper",
    # 例: "AUROC": "0.847 ± 0.023",
}
for k, v in results.items():
    print(f"{k}: {v}")
\`\`\`

\`[cell:figure-manifest]\` — 図マニフェストセル：
\`\`\`python
from pathlib import Path
figures = {
    "Figure 1": "figures/correlation_heatmap.png",
    # paper.md で参照する全ての図
}
for name, path in figures.items():
    exists = Path(path).exists()
    print(f"{name}: {path} -> {'✅' if exists else '❌ MISSING'}")
\`\`\`

⚠️ paper.md では \`[cell:results-summary]\` に含まれる値のみを引用する。
⚠️ paper.md では \`[cell:figure-manifest]\` に ✅ で存在確認済みの図のみを参照する。

📌 **引用ガイダンス**（Claims数の適正化）：
- 主要な結果（精度、有意差、p値、効果量など）を優先的に引用する
- 中間計算値やパラメータ設定値は本文に埋め込まず、Methodsでセル参照にとどめる
- 目安: 1セクションあたり定量的引用は5-8個が適切。過剰な数値の羅列を避ける
- Results セクションに集中させ、Abstract/Conclusion では主要指標のみ引用する

📌 **値の転記ルール**（VM削減のため厳守）：
- \`[cell:results-summary]\` の出力値を**そのまま**paper.mdに転記する（丸めない、単位を変えない）
- 例: セル出力が \`0.8734\` なら paper.md にも "0.8734" と記載する（"0.87" や "87.3%" にしない）
- 百分率変換や丸めが必要な場合は、**セル内で変換してから**出力し、変換後の値を引用する
- 表の数値もセル出力と完全一致させる

**ステップ4: paper.md 作成**（最重要 — 絶対に省略禁止）
⚠️⚠️⚠️ **paper.md は本タスクの最終目標です。** ⚠️⚠️⚠️
実験途中でも時間切れでも、**必ず paper.md を作成・保存すること**。

📄 **paper.md** — 以下のセクション構成：
  - **Abstract**: 200語以上。目的・手法・主要結果・意義
  - **Introduction**: 先行研究の位置づけと新規性
  - **Methods**: 手法・パラメータ・NatureLM/GALACTICA 使用状況（直接呼び出し）。Python コードを含む
  - **Results**: 定量結果を表形式で提示。交差検証の標準偏差付き。全数値に \`[cell:<id>]\` 引用。NatureLM予測とGALACTICA検証の結果を含む
  - **Discussion**: 結果の解釈・限界・先行研究比較。合成データの限界・実世界への一般化可能性を批判的に議論
  - **Conclusion**: 主要知見と今後の課題
  - **References**: DOI付き5件以上
  - **Reproducibility**: 乱数シード・Pythonバージョン・パッケージバージョン

📄 **report.md** — 任意。作成する場合は paper.md の要約とし、新しい未引用の数値を含めない

**補足ルール**:
- \`%\` と小数（\`83.0%\` = \`0.83\`）、指数表記（\`8.3e-1\`）は等価とみなされる
- DOI・年号・セクション番号・文献引用には \`[cell:<id>]\` 不要（自動除外）
- value_mismatches は repair-blocking gate ではないが、**意図的に残さない**こと。値がセル出力に含まれない場合は、セルを再実行するか記述を修正する
- データ（モックデータ含む）の生成方法を明記し、可能であれば \`data/raw/\` に保存する
- \`/validate\` 呼び出し前に paper.md の骨格（全セクション + 主要数値）が存在していること

`;

  return newPrefix + '---\n' + body;
}

function loadExperiments() {
  const dirs = fs.readdirSync(ROUND9_DIR)
    .filter(d => d.startsWith('SCI-'))
    .sort();

  return dirs.map(id => {
    const promptPath = path.join(ROUND9_DIR, id, 'input_prompt.txt');
    const originalPrompt = fs.readFileSync(promptPath, 'utf-8');
    const prompt = generateRound15Prompt(originalPrompt);
    const isRelevant = isGALACTICARelevant(originalPrompt);
    const domain = classifyDomain(originalPrompt);
    const afterSep = originalPrompt.split('---\n')[1] || originalPrompt;
    const titleMatch = afterSep.match(/^(.+?)(?:\n|。)/);
    const title = titleMatch ? titleMatch[1].trim().slice(0, 60) : id;
    return { id, title, prompt, galactica_relevant: isRelevant, domain };
  });
}

async function downloadFile(projectId, fileId) {
  return new Promise((resolve, reject) => {
    const url = new URL(`/api/projects/${projectId}/files/${fileId}/download`, BASE_URL);
    http.get({ hostname: url.hostname, port: url.port, path: url.pathname }, (res) => {
      const chunks = [];
      res.on('data', chunk => chunks.push(chunk));
      res.on('end', () => resolve(Buffer.concat(chunks)));
      res.on('error', reject);
    }).on('error', reject);
  });
}

/**
 * Wait for the latest run of a project to reach a terminal state.
 * prevRunCount: number of runs that existed before triggering the new run.
 *   If provided, waits until a NEW run appears (runs.length > prevRunCount)
 *   before checking for completion.
 * Returns { completed: boolean, status: string }
 */
async function waitForRunCompletion(projectId, workerId, startTime, label = '', prevRunCount = 0) {
  let lastStatus = '';
  let newRunDetected = (prevRunCount === 0); // if not tracking, skip detection
  for (let i = 0; i < Math.ceil(TIMEOUT_MS / POLL_INTERVAL_MS); i++) {
    await sleep(POLL_INTERVAL_MS);
    const elapsed = Math.round((Date.now() - startTime) / 1000);
    try {
      const runs = await apiGet(`/projects/${projectId}/runs?limit=10`);
      if (!Array.isArray(runs) || runs.length === 0) {
        if (elapsed > 30 && lastStatus === '') {
          log(workerId, `  ⏳ ${label}${elapsed}s [waiting for run to start]`);
          lastStatus = 'waiting';
        }
        continue;
      }

      // Wait for a new run to appear if prevRunCount is specified
      if (!newRunDetected) {
        if (runs.length > prevRunCount) {
          newRunDetected = true;
          log(workerId, `  ⏳ ${label}new run detected (${runs.length} runs)`);
        } else {
          if (elapsed > 30 && lastStatus === '') {
            log(workerId, `  ⏳ ${label}${elapsed}s [waiting for new run to appear]`);
            lastStatus = 'waiting-for-new';
          }
          continue;
        }
      }

      const latestRun = runs[0];
      if (latestRun.status === 'completed') {
        log(workerId, `  ✓ ${label}completed in ${elapsed}s`);
        return { completed: true, status: 'completed' };
      } else if (latestRun.status === 'failed') {
        log(workerId, `  ✗ ${label}failed after ${elapsed}s`);
        return { completed: false, status: 'failed' };
      }
      const newStatus = latestRun.status || 'unknown';
      if (newStatus !== lastStatus) {
        log(workerId, `  ⏳ ${label}${elapsed}s [${newStatus}]`);
        lastStatus = newStatus;
      } else if (elapsed % 120 < (POLL_INTERVAL_MS / 1000)) {
        log(workerId, `  ⏳ ${label}${elapsed}s...`);
      }
    } catch (e) {
      log(workerId, `  ⚠ poll error: ${e.message}`);
      // polling error, retry
    }
  }
  return { completed: false, status: 'timeout' };
}

/**
 * Run validate → repair loop (v3.4.8 — value match precision).
 * Returns { validation, repairIterations, repairHistory }
 */
async function runValidateRepairLoop(projectId, ws, workerId) {
  const repairHistory = [];
  let validation = null;
  let iteration = 0;

  for (iteration = 0; iteration <= MAX_REPAIR_ITERATIONS; iteration++) {
    // Run validation
    try {
      validation = await apiPostWithCsrf(`/projects/${projectId}/validate`);
    } catch (e) {
      log(workerId, `  ⚠ validation failed: ${e.message}`);
      return { validation: null, repairIterations: iteration, repairHistory };
    }

    if (!validation || typeof validation !== 'object' || !validation.gates) {
      log(workerId, `  ⚠ validation returned unexpected format`);
      return { validation, repairIterations: iteration, repairHistory };
    }

    const gates = validation.gates || [];
    const gatesPassed = gates.filter(g => g.passed).length;
    const totalGates = gates.length;
    const uncitedCount = (validation.uncited_claims || []).length;
    const totalClaims = (validation.claims || []).length;
    const valueMismatches = (validation.value_mismatches || []).length;
    const figureOrphans = (validation.figure_orphans || []).length;
    const reportThinness = validation.report_thinness || null;
    const gateResult = (name) => gates.find(g => g.name === name)?.passed ? '✓' : '✗';

    log(workerId, `  🔬 validate[${iteration}]: ${gatesPassed}/${totalGates} gates | seed=${gateResult('seed_presence')} env=${gateResult('env_capture')} no-err=${gateResult('no_error_in_cited')} cov=${gateResult('citation_coverage')} | claims=${totalClaims} uncited=${uncitedCount} val_mismatch=${valueMismatches} fig_orphan=${figureOrphans}${reportThinness ? ' ⚠thin=' + JSON.stringify(reportThinness) : ''}`);

    // If all gates pass, we're done
    if (gatesPassed === totalGates) {
      log(workerId, `  ✅ all gates passed on iteration ${iteration}`);
      break;
    }

    // If we've exhausted repair iterations, call postmortem and stop
    if (iteration >= MAX_REPAIR_ITERATIONS) {
      log(workerId, `  ⚠ max repair iterations (${MAX_REPAIR_ITERATIONS}) reached, ${gatesPassed}/${totalGates} gates passed`);
      // Auto-postmortem: v3.4.2 feature
      try {
        const postmortem = await apiPostWithCsrf(`/projects/${projectId}/validate/postmortem`);
        if (postmortem && postmortem.markdown_summary) {
          log(workerId, `  📋 postmortem generated (${postmortem.markdown_summary.length} chars)`);
          validation._postmortem = postmortem;
        }
      } catch (e) {
        log(workerId, `  ⚠ postmortem endpoint failed: ${e.message}`);
      }
      break;
    }

    // Call repair endpoint to get repair prompt
    let repairPayload = null;
    try {
      repairPayload = await apiPostWithCsrf(`/projects/${projectId}/validate/repair`);
    } catch (e) {
      log(workerId, `  ⚠ repair endpoint failed: ${e.message}`);
      break;
    }

    if (!repairPayload || typeof repairPayload !== 'object') {
      log(workerId, `  ⚠ repair returned unexpected format`);
      break;
    }

    // Extract repair prompt and send to agent
    const repairPrompt = repairPayload.prompt || repairPayload.repair_prompt || '';
    if (!repairPrompt) {
      log(workerId, `  ⚠ empty repair prompt — skipping repair`);
      break;
    }

    log(workerId, `  🔧 repair[${iteration}]: sending repair prompt (${repairPrompt.length} chars)`);
    repairHistory.push({
      iteration,
      gates_before: gatesPassed,
      total_gates: totalGates,
      uncited_before: uncitedCount,
      repair_prompt_length: repairPrompt.length
    });

    // Get current run count BEFORE sending chat (to detect new run later)
    let prevRunCount = 0;
    try {
      const currentRuns = await apiGet(`/projects/${projectId}/runs?limit=10`);
      prevRunCount = Array.isArray(currentRuns) ? currentRuns.length : 0;
    } catch (_) {}

    sendChat(ws, repairPrompt);

    // Wait a moment for AIRA to create the new run
    await sleep(3000);

    // Wait for repair run to complete (pass prevRunCount to detect new run)
    const repairStart = Date.now();
    const repairResult = await waitForRunCompletion(projectId, workerId, repairStart, `repair[${iteration}] `, prevRunCount);
    if (!repairResult.completed) {
      log(workerId, `  ⚠ repair run did not complete (${repairResult.status})`);
      break;
    }
  }

  return { validation, repairIterations: iteration, repairHistory };
}

async function runExperiment(experiment, workerId) {
  const expDir = path.join(RESULTS_DIR, experiment.id);
  fs.mkdirSync(expDir, { recursive: true });
  const prompt = experiment.prompt;

  const domainEmoji = { molecular: '🧪', protein: '🧬', materials: '⚗️', genomics: '🧬', 'general-science': '📊', 'non-science': '📝' };
  log(workerId, `▶ ${experiment.id}: ${experiment.title} ${domainEmoji[experiment.domain] || '📊'} [${experiment.domain}]`);
  fs.writeFileSync(path.join(expDir, 'input_prompt.txt'), prompt, 'utf-8');

  // Delete any existing project with this name
  const existing = await apiGet('/projects');
  if (Array.isArray(existing)) {
    for (const p of existing.filter(p => p.name === experiment.id)) {
      await apiDeleteWithCsrf(`/projects/${p.id}`);
      await sleep(1000);
    }
  }

  // Create project via API with skill
  const project = await apiPostWithCsrf('/projects', {
    name: experiment.id,
    description: experiment.title,
    skillSetId: CO_SCIENTIST_SKILL
  });
  if (!project || !project.id) {
    throw new Error(`Project creation failed: ${JSON.stringify(project)}`);
  }
  const projectId = project.id;
  log(workerId, `  project: ${projectId}`);

  // Add NatureLM MCP to the project (external registration)
  try {
    await apiPostWithCsrf(`/projects/${projectId}/mcp`, {
      name: 'naturelm',
      type: 'http',
      config: {
        url: 'http://172.17.0.1:3001/mcp',
        description: 'NatureLM MCP Server — IBM science foundation model for quantitative prediction. Tools: predict_logp, predict_property, generate_smiles, retrosynthesis, ask_naturelm, generate_protein_sequence, predict_material_composition.'
      }
    });
    log(workerId, `  naturelm MCP added`);
  } catch (e) {
    log(workerId, `  ⚠ naturelm MCP add failed: ${e.message}`);
  }

  // Add GALACTICA MCP to the project (external registration)
  try {
    await apiPostWithCsrf(`/projects/${projectId}/mcp`, {
      name: 'galactica',
      type: 'http',
      config: {
        url: 'http://172.17.0.1:3002/mcp',
        description: 'GALACTICA MCP Server — Meta AI science foundation model for scientific validation. Tools: generate_text, predict_citations, generate_latex, generate_molecule, predict_protein_annotations, scientific_qa, reasoning, summarize, model_info.'
      }
    });
    log(workerId, `  galactica MCP added`);
  } catch (e) {
    log(workerId, `  ⚠ galactica MCP add failed: ${e.message}`);
  }

  // Verify MCP tools are available
  try {
    const mcpTools = await apiGet(`/projects/${projectId}/mcp/tools`);
    const naturelmTools = Array.isArray(mcpTools) 
      ? mcpTools.filter(t => t.server === 'naturelm' || (t.name && t.name.includes('naturelm')))
      : [];
    const galacticaTools = Array.isArray(mcpTools) 
      ? mcpTools.filter(t => t.server === 'galactica' || (t.name && t.name.includes('galactica')))
      : [];
    log(workerId, `  MCP tools: ${Array.isArray(mcpTools) ? mcpTools.length : 0} total, ${naturelmTools.length} naturelm, ${galacticaTools.length} galactica`);
  } catch (e) {
    log(workerId, `  ⚠ MCP tools check failed: ${e.message}`);
  }

  // Open WebSocket and send chat message (keep open during run)
  const ws = await openProjectWs(projectId);
  sendChat(ws, prompt);
  log(workerId, `  prompt sent (${prompt.length} chars)`);

  // Wait for initial run completion
  const startTime = Date.now();
  const runResult = await waitForRunCompletion(projectId, workerId, startTime);

  // Collect messages
  const messages = await apiGet(`/projects/${projectId}/messages?limit=200&offset=0`);
  const assistantMessages = Array.isArray(messages)
    ? messages.filter(m => m.role === 'assistant')
    : [];
  const output = assistantMessages.map(m => m.content).join('\n\n---\n\n');
  fs.writeFileSync(path.join(expDir, 'output_response.md'), output, 'utf-8');

  // Collect files
  let files = await apiGet(`/projects/${projectId}/files`);
  let fileList = Array.isArray(files) ? files : [];
  log(workerId, `  📁 ${fileList.length} files`);

  // Check if paper.md exists; if not, send follow-up and wait
  let hasPaper = fileList.some(f => (f.file_path || f.filename || '').includes('paper.md'));
  if (runResult.completed && !hasPaper) {
    log(workerId, `  ⚠ paper.md missing — sending follow-up prompt`);
    sendChat(ws, 'paper.md がまだ作成されていません。今すぐ paper.md を作成してください。これは最重要の成果物です。Abstract, Introduction, Methods, Results, Discussion, Conclusion, References のセクションを含む学術論文形式で作成してください。');

    const retryStart = Date.now();
    await waitForRunCompletion(projectId, workerId, retryStart, 'retry ');

    // Re-collect files after retry
    files = await apiGet(`/projects/${projectId}/files`);
    fileList = Array.isArray(files) ? files : [];
    hasPaper = fileList.some(f => (f.file_path || f.filename || '').includes('paper.md'));
    log(workerId, `  📁 after retry: ${fileList.length} files, paper.md: ${hasPaper ? 'YES' : 'NO'}`);
  }

  // Download files
  if (fileList.length > 0) {
    const filesDir = path.join(expDir, 'files');
    fs.mkdirSync(filesDir, { recursive: true });

    for (const file of fileList) {
      try {
        const fileName = file.file_path || file.filename || file.id;
        const filePath = path.join(filesDir, fileName);
        fs.mkdirSync(path.dirname(filePath), { recursive: true });
        const fileBuffer = await downloadFile(projectId, file.id);
        if (fileBuffer && fileBuffer.length > 0) {
          fs.writeFileSync(filePath, fileBuffer);
        }
      } catch (e) {
        // skip failed file download
      }
    }
    fs.writeFileSync(path.join(expDir, 'file_list.json'), JSON.stringify(fileList, null, 2), 'utf-8');
  }

  // Check NatureLM and GALACTICA MCP usage in messages
  const allMessages = (Array.isArray(messages) ? messages : []).map(m => m.content || '').join(' ');
  const naturelmUsed = allMessages.includes('naturelm') ||
    allMessages.includes('NatureLM') ||
    allMessages.includes('predict_logp') ||
    allMessages.includes('predict_property') ||
    allMessages.includes('generate_smiles') ||
    allMessages.includes('retrosynthesis') ||
    allMessages.includes('ask_naturelm') ||
    allMessages.includes('generate_protein_sequence') ||
    allMessages.includes('predict_material_composition');
  const galacticaUsed = allMessages.includes('galactica') ||
    allMessages.includes('GALACTICA') ||
    allMessages.includes('Galactica') ||
    allMessages.includes('generate_molecule') ||
    allMessages.includes('scientific_qa') ||
    allMessages.includes('predict_protein_annotations') ||
    allMessages.includes('generate_text') ||
    allMessages.includes('predict_citations') ||
    allMessages.includes('generate_latex') ||
    allMessages.includes('reasoning') ||
    allMessages.includes('summarize');
  const jupyterUsed = allMessages.includes('jupyter') ||
    allMessages.includes('Jupyter') ||
    allMessages.includes('execute_code') ||
    allMessages.includes('run_cell') ||
    allMessages.includes('```python') ||
    allMessages.includes('In [') ||
    allMessages.includes('Out[');

  // Run validate → repair loop (v3.4.8)
  const { validation, repairIterations, repairHistory } = await runValidateRepairLoop(projectId, ws, workerId);

  // Save validation results
  if (validation) {
    fs.writeFileSync(path.join(expDir, 'validation.json'), JSON.stringify(validation, null, 2), 'utf-8');
  }
  if (repairHistory.length > 0) {
    fs.writeFileSync(path.join(expDir, 'repair-history.json'), JSON.stringify(repairHistory, null, 2), 'utf-8');
  }

  // Collect execution trace
  let trace = null;
  try {
    trace = await apiGet(`/projects/${projectId}/notebook/trace?latest=1`);
    if (trace && (Array.isArray(trace) ? trace.length > 0 : true)) {
      fs.writeFileSync(path.join(expDir, 'execution-trace.json'), JSON.stringify(trace, null, 2), 'utf-8');
      log(workerId, `  📋 execution trace saved`);
    }
  } catch (e) {
    // trace API may not be available — not critical
  }

  // Re-collect messages after repair loop (may have additional assistant messages)
  const finalMessages = await apiGet(`/projects/${projectId}/messages?limit=200&offset=0`);
  const finalAssistantMessages = Array.isArray(finalMessages)
    ? finalMessages.filter(m => m.role === 'assistant')
    : [];
  if (finalAssistantMessages.length > assistantMessages.length) {
    const finalOutput = finalAssistantMessages.map(m => m.content).join('\n\n---\n\n');
    fs.writeFileSync(path.join(expDir, 'output_response.md'), finalOutput, 'utf-8');
  }

  // Re-download files after repair (paper.md may have been updated)
  if (repairHistory.length > 0) {
    files = await apiGet(`/projects/${projectId}/files`);
    fileList = Array.isArray(files) ? files : [];
    if (fileList.length > 0) {
      const filesDir = path.join(expDir, 'files');
      for (const file of fileList) {
        try {
          const fileName = file.file_path || file.filename || file.id;
          const filePath = path.join(filesDir, fileName);
          fs.mkdirSync(path.dirname(filePath), { recursive: true });
          const fileBuffer = await downloadFile(projectId, file.id);
          if (fileBuffer && fileBuffer.length > 0) {
            fs.writeFileSync(filePath, fileBuffer);
          }
        } catch (e) { /* skip */ }
      }
    }
    hasPaper = fileList.some(f => (f.file_path || f.filename || '').includes('paper.md'));
  }

  // Save metadata
  const runs = await apiGet(`/projects/${projectId}/runs?limit=5`);
  const duration = Math.round((Date.now() - startTime) / 1000);
  const metadata = {
    experiment_id: experiment.id,
    title: experiment.title,
    project_id: projectId,
    prompt_length: prompt.length,
    output_length: output.length,
    file_count: fileList.length,
    duration_seconds: duration,
    status: runResult.completed ? (hasPaper ? 'success' : 'success-no-paper') : (Array.isArray(runs) && runs[0]?.status || 'unknown'),
    galactica_relevant: experiment.galactica_relevant,
    domain: experiment.domain,
    naturelm_used: naturelmUsed,
    galactica_used: galacticaUsed,
    jupyter_used: jupyterUsed,
    aira_version: 'v3.4.8',
    validation: validation && typeof validation === 'object' && validation.gates ? {
      seed_presence: validation.gates.find(g => g.name === 'seed_presence')?.passed || false,
      env_capture: validation.gates.find(g => g.name === 'env_capture')?.passed || false,
      no_error_in_cited: validation.gates.find(g => g.name === 'no_error_in_cited')?.passed || false,
      citation_coverage: validation.gates.find(g => g.name === 'citation_coverage')?.passed || false,
      total_claims: (validation.claims || []).length,
      uncited_claims: (validation.uncited_claims || []).length,
      unknown_citations: (validation.unknown_citations || []).length,
      value_mismatches: (validation.value_mismatches || []).length,
      figure_orphans: (validation.figure_orphans || []).length,
      report_thinness: validation.report_thinness || null,
      gates_passed: (validation.gates || []).filter(g => g.passed).length,
      postmortem: validation._postmortem ? validation._postmortem.markdown_summary : null
    } : null,
    repair: {
      iterations: repairIterations,
      history: repairHistory
    },
    timestamp: new Date().toISOString()
  };
  fs.writeFileSync(path.join(expDir, 'metadata.json'), JSON.stringify(metadata, null, 2), 'utf-8');

  // Save postmortem report if generated
  if (validation && validation._postmortem && validation._postmortem.markdown_summary) {
    fs.writeFileSync(path.join(expDir, 'postmortem.md'), validation._postmortem.markdown_summary, 'utf-8');
  }

  // Close WebSocket and delete project
  try { ws.close(); } catch (_) {}
  if (hasPaper) {
    await apiDeleteWithCsrf(`/projects/${projectId}`);
  } else {
    log(workerId, `  ⚠ keeping project ${projectId} for investigation (no paper.md)`);
  }
  log(workerId, `  ✓ done (${duration}s, ${fileList.length} files, NatureLM: ${naturelmUsed ? '✓' : '✗'}, GALACTICA: ${galacticaUsed ? '✓' : '✗'}, Jupyter: ${jupyterUsed ? '✓' : '✗'}, Gates: ${metadata.validation ? metadata.validation.gates_passed + '/4' : 'N/A'}, Repairs: ${repairIterations})`);

  return metadata;
}

async function worker(workerId, queue, results) {
  log(workerId, 'ready');

  while (true) {
    const experiment = queue.shift();
    if (!experiment) break;

    try {
      const meta = await runExperiment(experiment, workerId);
      results.push(meta);
    } catch (err) {
      log(workerId, `  ✗ ERROR: ${err.message}`);
      results.push({
        experiment_id: experiment.id,
        status: 'error',
        error: err.message,
        galactica_relevant: experiment.galactica_relevant,
        timestamp: new Date().toISOString()
      });
    }
    await sleep(3000);
  }

  log(workerId, 'finished');
}

async function main() {
  const allExperiments = loadExperiments();
  const relevantCount = allExperiments.filter(e => e.galactica_relevant).length;
  console.log(`📂 Loaded ${allExperiments.length} experiments (${relevantCount} GALACTICA-relevant)`);

  // Filter by CLI args, or run all
  const filterIds = process.argv.slice(2);
  let experiments;
  if (filterIds.length > 0) {
    if (filterIds.length === 1 && filterIds[0].includes(':')) {
      const [startId, endId] = filterIds[0].split(':');
      const startNum = parseInt(startId.replace('SCI-', ''));
      const endNum = parseInt(endId.replace('SCI-', ''));
      experiments = allExperiments.filter(e => {
        const n = parseInt(e.id.replace('SCI-', ''));
        return n >= startNum && n <= endNum;
      });
    } else {
      experiments = allExperiments.filter(e => filterIds.includes(e.id));
    }
  } else {
    // Default: use FIXED_SUBSET (30 benchmark experiments)
    experiments = allExperiments.filter(e => FIXED_SUBSET.has(e.id));
  }

  // Skip already completed experiments
  const completed = [];
  const pending = [];
  for (const exp of experiments) {
    const metaPath = path.join(RESULTS_DIR, exp.id, 'metadata.json');
    if (fs.existsSync(metaPath)) {
      try {
        const meta = JSON.parse(fs.readFileSync(metaPath, 'utf-8'));
        if (meta.status === 'success') {
          completed.push(exp.id);
          continue;
        }
      } catch (_) {}
    }
    pending.push(exp);
  }

  console.log(`\n🔬 AIRA Round15 Experiment Runner — Semantic Verification (v3.4.8) + NatureLM + GALACTICA + Jupyter (${NUM_WORKERS} workers)`);
  console.log(`📋 Total: ${experiments.length} | Already done: ${completed.length} | To run: ${pending.length}`);
  console.log(`🧪 GALACTICA-relevant: ${pending.filter(e => e.galactica_relevant).length}/${pending.length}`);
  console.log(`📂 Results → ${RESULTS_DIR}`);
  console.log(`🌐 AIRA: ${BASE_URL}`);
  console.log(`🔧 Max repair iterations: ${MAX_REPAIR_ITERATIONS}\n`);

  if (pending.length === 0) {
    console.log('All experiments already completed!');
    return;
  }

  fs.mkdirSync(RESULTS_DIR, { recursive: true });

  const queue = [...pending];
  const results = [];
  const totalStart = Date.now();

  // Launch workers with stagger delay to avoid overwhelming AIRA
  const workers = [];
  for (let i = 0; i < Math.min(NUM_WORKERS, pending.length); i++) {
    workers.push(worker(i + 1, queue, results));
    if (i < NUM_WORKERS - 1) await sleep(5000); // 5s stagger between workers
  }

  await Promise.all(workers);

  const totalDuration = ((Date.now() - totalStart) / 1000 / 60).toFixed(1);
  const successful = results.filter(r => r.status === 'success').length;
  const failed = results.filter(r => r.status !== 'success').length;
  const galacticaUsedCount = results.filter(r => r.galactica_used).length;
  const naturelmUsedCount = results.filter(r => r.naturelm_used).length;
  const jupyterUsedCount = results.filter(r => r.jupyter_used).length;
  const repairUsedCount = results.filter(r => r.repair && r.repair.iterations > 0).length;
  const avgRepairIters = results.filter(r => r.repair).length > 0
    ? (results.filter(r => r.repair).reduce((s, r) => s + (r.repair?.iterations || 0), 0) / results.filter(r => r.repair).length).toFixed(1)
    : '0';

  // Save summary
  const summary = {
    version: 'round15-semantic-verification',
    aira_version: 'v3.4.8',
    prompt_source: 'round-11 + semantic-verification + figure-provenance + postmortem + NatureLM MCP + GALACTICA MCP + Jupyter MCP',
    total_experiments: experiments.length,
    successful: successful + completed.length,
    failed,
    skipped: completed.length,
    galactica_relevant_count: experiments.filter(e => e.galactica_relevant).length,
    naturelm_actually_used: naturelmUsedCount,
    galactica_actually_used: galacticaUsedCount,
    jupyter_actually_used: jupyterUsedCount,
    repair_used_count: repairUsedCount,
    avg_repair_iterations: parseFloat(avgRepairIters),
    total_duration_minutes: parseFloat(totalDuration),
    workers: NUM_WORKERS,
    results
  };
  fs.writeFileSync(path.join(RESULTS_DIR, 'summary.json'), JSON.stringify(summary, null, 2), 'utf-8');

  console.log(`\n${'='.repeat(60)}`);
  console.log(`🏁 ALL DONE | Success: ${successful} | Failed: ${failed} | Skipped: ${completed.length}`);
  console.log(`🧪 NatureLM: ${naturelmUsedCount}/${successful} | GALACTICA: ${galacticaUsedCount}/${successful} | Jupyter: ${jupyterUsedCount}/${successful}`);
  console.log(`🔧 Repair used: ${repairUsedCount}/${successful} | Avg iterations: ${avgRepairIters}`);
  console.log(`   Duration: ${totalDuration} min | Workers: ${NUM_WORKERS}`);
  console.log(`${'='.repeat(60)}\n`);
}

main().catch(err => {
  console.error('Fatal:', err);
  process.exit(1);
});
