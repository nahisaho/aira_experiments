/**
 * AIRA Experiment Runner - Round 9
 * AIRA v3.2.1 — Computational Provenance edition.
 * Everything from Round 8 plus:
 *   - Pillar 1: Notebook execution trace auto-capture
 *   - Pillar 2: Numeric claim → [cell:<id>] citation linter
 *   - Pillar 3: Reproducibility gates (seed/env/no-error/citation-coverage)
 *   - Pillar 4: data/raw/ convention
 *   - Co-Scientist v4.7.0 with provenance instructions
 * After each experiment, runs POST /validate to collect gate results.
 * No timeout - waits until completion.
 * 5 Parallel Workers.
 * Uses direct HTTP API calls + WebSocket for chat trigger.
 * Target: http://localhost:3000 (AIRA)
 */

const fs = require('fs');
const path = require('path');
const http = require('http');
const WebSocket = require('ws');

const BASE_URL = process.env.AIRA_URL || 'http://192.168.1.15:3000';
const RESULTS_DIR = path.join(__dirname, '..', 'results', 'round9');
const ROUND9_DIR = path.join(__dirname, '..', '..', 'aira', 'co-scientist-optimization', 'round-9');
const TIMEOUT_MS = 24 * 60 * 60 * 1000; // 24 hours (effectively no timeout)
const POLL_INTERVAL_MS = 15000;
const NUM_WORKERS = parseInt(process.env.WORKERS || '5', 10);

const CO_SCIENTIST_SKILL = '3776b03b-1800-4640-9f4a-8d5c1d418fa1';

function sleep(ms) { return new Promise(r => setTimeout(r, ms)); }

function log(workerId, msg) {
  const ts = new Date().toISOString().slice(11, 19);
  console.log(`[${ts}][W${workerId}] ${msg}`);
}

// Direct HTTP helper (no CORS)
function httpRequest(method, urlPath, body) {
  return new Promise((resolve, reject) => {
    const url = new URL(urlPath, BASE_URL);
    const options = {
      hostname: url.hostname,
      port: url.port,
      path: url.pathname + url.search,
      method,
      headers: {}
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
    req.on('error', reject);
    req.end();
  });
}

/**
 * Open WebSocket to project. Returns WS instance that must stay open during run.
 * Includes automatic ping to keep connection alive.
 */
function openProjectWs(projectId) {
  return new Promise((resolve, reject) => {
    const wsUrl = BASE_URL.replace('http://', 'ws://') + `/ws/projects/${projectId}/chat`;
    const ws = new WebSocket(wsUrl, {
      headers: { 'Origin': 'http://localhost:3000' }
    });
    // Keep-alive ping every 30s
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
      reject(err);
    });
    setTimeout(() => {
      if (ws.readyState !== WebSocket.OPEN) {
        ws.terminate();
        reject(new Error('WebSocket connection timeout'));
      }
    }, 10000);
  });
}

function sendChat(ws, content) {
  ws.send(JSON.stringify({ type: 'chat', content }));
}

/**
 * Classify experiment into domain category for targeted GALACTICA instructions.
 * Returns: 'molecular', 'protein', 'materials', 'genomics', 'general-science', 'non-science'
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

  // General science keywords
  const sciKeywords = ['モデル', '予測', 'シミュレーション', 'アルゴリズム', '最適化', 'ニューラル',
    '機械学習', 'ディープラーニング', '統計', 'パイプライン', 'フレームワーク'];
  if (sciKeywords.some(k => prompt.includes(k))) return 'general-science';
  return 'non-science';
}

/**
 * Determine if an experiment topic is relevant to GALACTICA's capabilities
 */
function isGALACTICARelevant(prompt) {
  const domain = classifyDomain(prompt);
  return ['molecular', 'protein', 'materials', 'genomics'].includes(domain);
}

/**
 * Domain-specific NatureLM tool instructions (quantitative prediction)
 */
function getNatureLMInstructions(domain) {
  switch (domain) {
    case 'molecular':
      return `  - \`generate_smiles\`: 目的の性質を持つ候補分子を複数生成し、探索空間を拡大する
  - \`predict_logp\`, \`predict_property\`: 生成した分子の物性を予測し、実験条件のベースラインを設定する
  - \`retrosynthesis\`: 合成可能性を検証し、実現可能な候補に絞り込む
  - \`ask_naturelm\`: 分子メカニズムに関する定量的パラメータ（結合エネルギー、IC50推定値、LogP等）を取得する`;
    case 'protein':
      return `  - \`generate_protein_sequence\`: 目的の機能を持つタンパク質配列を生成する
  - \`predict_property\`: タンパク質の物性予測を行い、設計パラメータを定量化する
  - \`ask_naturelm\`: 構造-活性相関、安定性条件、フォールディング特性に関する知見を取得する`;
    case 'materials':
      return `  - \`predict_material_composition\`: 目的の性質を持つ材料組成を予測する
  - \`predict_property\`: 候補材料の物性を予測し、スクリーニング基準を設定する
  - \`ask_naturelm\`: 材料の安定性、劣化メカニズム、界面特性に関する定量的知見を取得する`;
    case 'genomics':
      return `  - \`ask_naturelm\`: 生物学的メカニズムに関する定量的パラメータ（結合自由エネルギー、反応速度定数等）を取得する`;
    default:
      return `  - \`ask_naturelm\`: 研究テーマに関連する科学的知見・定量的パラメータを取得する`;
  }
}

/**
 * Domain-specific GALACTICA tool instructions (scientific validation & citations)
 */
function getGALACTICAInstructions(domain) {
  switch (domain) {
    case 'molecular':
      return `  - \`generate_molecule\`: 候補分子（SMILES）を複数生成し、NatureLMの予測と比較する
  - \`scientific_qa\`: 分子メカニズムの科学的妥当性を検証し、NatureLMの予測値の信頼性を評価する
  - \`predict_citations\`: 関連する先行研究の引用を予測し、文献調査を補完する
  - \`reasoning\`: 反応機構や合成経路の推論を行う`;
    case 'protein':
      return `  - \`predict_protein_annotations\`: 対象タンパク質のアミノ酸配列から機能予測・アノテーションを取得する
  - \`scientific_qa\`: タンパク質設計の科学的妥当性を検証する
  - \`predict_citations\`: 関連文献を予測し、文献レビューを補完する`;
    case 'materials':
      return `  - \`scientific_qa\`: 材料特性の科学的妥当性を検証し、NatureLMの予測と照合する
  - \`generate_molecule\`: 候補材料の分子構造を生成し、スクリーニング対象を拡大する
  - \`reasoning\`: 材料特性の物理的推論を行い、従来の探索範囲を超えた候補も検討する
  - \`generate_latex\`: 材料科学の数式（状態方程式、拡散方程式等）を生成する`;
    case 'genomics':
      return `  - \`scientific_qa\`: 生物学的メカニズムの科学的妥当性を検証する
  - \`predict_citations\`: ゲノミクス分野の関連文献を予測し、文献レビューを補完する`;
    default:
      return `  - \`scientific_qa\`: 研究テーマに関連する科学的知見を取得し、実験設計の妥当性を検証する
  - \`predict_citations\`: 関連文献を予測し、文献調査を補完する`;
  }
}

/**
 * Combined NatureLM + GALACTICA instructions for each domain
 */
function getCombinedModelInstructions(domain) {
  const naturelmInstr = getNatureLMInstructions(domain);
  const galacticaInstr = getGALACTICAInstructions(domain);

  return `- **NatureLM MCP**（定量予測）を以下のように活用すること：
${naturelmInstr}

- **GALACTICA MCP**（科学的検証・引用予測）を以下のように活用すること：
${galacticaInstr}

- ⚠️ **両モデルの相互検証**: NatureLMの定量予測結果をGALACTICAのscientific_qaで科学的に検証し、矛盾がないか確認すること
- 両モデルの予測結果はResultsセクションに定量的に記載し、一致・不一致を明示すること`;
}

/**
 * Generate round9 optimized prompt
 * Round 8 base (NatureLM + GALACTICA + Jupyter) plus:
 *   - [cell:<id>] citation for every numeric claim
 *   - Reproducibility: seed fixing, env capture, data/raw/ convention
 *   - Provenance-aware reporting
 */
function generateRound9Prompt(originalPrompt) {
  const domain = classifyDomain(originalPrompt);
  const combinedInstructions = getCombinedModelInstructions(domain);

  // Split prompt: before and after the "---" separator
  const parts = originalPrompt.split('---\n');
  const body = parts.slice(1).join('---\n');

  const newPrefix = `以下の研究テーマについて、次のステップで進めてください：

**ステップ1: 先行研究調査（ToolUniverse MCP 使用）**
- ToolUniverse MCP の学術検索ツール（Semantic Scholar / PubMed / Crossref 等）を使って、この研究テーマに関連する先行研究を調査する
- 検索キーワードを複数設定し、関連する最新論文（特に2020年以降）を少なくとも5件特定する
- 各論文のタイトル、著者、年、DOI、主要な知見・手法をまとめる
- 先行研究の課題・限界を整理する

**ステップ2: 実験計画 + NatureLM定量予測 + GALACTICA科学的検証**
- 先行研究の知見を踏まえて実験計画を立てる
${combinedInstructions}
- ⚠️ NatureLM / GALACTICA MCP ツールへの接続に失敗した場合は、（1）試行したツール名、（2）エラー内容、（3）代替手段をMethodsセクションに記録すること。ツール接続の成否にかかわらず、試行の記録は科学的透明性として重要である

**ステップ3: Pythonコード実装と実行（Jupyter MCP 使用）**
- 実験で必要なデータ処理・分析・シミュレーション・モデル構築をPythonコードとして実装すること
- Jupyter MCP を使ってPythonコードを実際に実行し、実行結果（数値・統計量・図表等）を取得すること
- コードでは以下のライブラリを積極的に活用すること：
  - データ処理: \`pandas\`, \`numpy\`
  - 機械学習: \`scikit-learn\`, \`xgboost\`, \`lightgbm\`（必要に応じて）
  - 可視化: \`matplotlib\`, \`seaborn\`
  - 統計分析: \`scipy.stats\`
  - 化学・分子: \`rdkit\`（分子関連テーマの場合）
- ⚠️ **コードの実行結果をpaper.mdに反映すること**: コードで得られた数値・統計量・p値・信頼区間等を論文に引用する。手計算や推測ではなく、実行結果に基づく数値を報告すること
- ⚠️ **再現性を確保すること**: 乱数シード（\`random_state=42\` 等）を固定し、実験条件をコード内にコメントとして記録すること
- コードはステップ5でpaper.mdのAppendixまたはMethodsセクションに含めること

**ステップ3.5: 計算来歴（Computational Provenance）の確保** ⚠️ 重要
- ⚠️ **[cell:<id>] 引用**: paper.md / report.md に数値結果を記載する際、その数値を算出した Jupyter セルを \`[cell:<id>]\` 形式で引用すること（例: \`AUROC = 0.83 [cell:3]\`）。すべての定量的主張に引用を付けること
- ⚠️ **乱数シード固定**: \`np.random.seed(42)\`, \`random.seed(42)\`, \`torch.manual_seed(42)\` 等を実験コードの冒頭で設定すること
- ⚠️ **環境記録**: 実験ノートブックの最初または最後のセルで \`!pip freeze\` を実行し、依存パッケージのバージョンを記録すること
- ⚠️ **データ出自**: 使用するデータ（モックデータ含む）の生成方法・パラメータを明記し、可能であれば \`data/raw/\` に保存すること
- ⚠️ **エラーセルの引用禁止**: stderr にエラーが出力されたセルの結果を引用しないこと。エラーを修正してから再実行すること

**ステップ4: 自己批判的検証**
- 先行研究調査とNatureLM/GALACTICA予測の結果を踏まえて実験結果を検証する
- 実験中もNatureLM MCPツール（定量予測）とGALACTICA MCPツール（科学的検証）を適宜活用して仮説検証・データ補完を行う
- ⚠️ **現実的な結果を報告すること**: AUC/AUROC/F1等が1.000（完璧）になった場合は、過学習・データリーク・評価の不備を疑い、交差検証の標準偏差付きで報告する。合成データであっても現実的なノイズを含めること
- ⚠️ **自己批判的に結果を検証すること**: NatureLM/GALACTICAの予測と一致する結果であっても、以下の観点から自己の実験を批判的に評価する：
  - この結果は合成データ/シミュレーションの前提条件にどの程度依存しているか？
  - 実世界のデータに適用した場合、同等の性能が期待できるか？
  - 自分の実験設計に含まれるバイアスや限界は何か？
  - NatureLMの定量予測値とGALACTICAの科学的検証が矛盾する場合、どちらがより信頼できるか？

**ステップ5: 成果物作成（最重要 — 絶対に省略禁止）**
⚠️⚠️⚠️ **このステップは本タスクの最終目標です。ステップ1〜4はすべてこのステップのための準備です。** ⚠️⚠️⚠️
以下の2ファイルは**必ず**作成すること。**他のすべてのステップより優先度が高い**。
実験の途中であっても、時間がかかっていても、**必ず paper.md ファイルを作成して保存すること**。
paper.md を作成せずにタスクを終了することは許可されていません。

📄 **paper.md** （必須・最優先）— 以下のセクション構成に従う学術論文：
  - **Abstract**: 200語以上。研究目的・手法・主要結果・意義を含む
  - **Introduction**: 先行研究の位置づけと研究の新規性を明記
  - **Methods**: 使用した手法、パラメータ、NatureLM MCPツールとGALACTICA MCPツールの使用/試行状況を記載。両モデルの役割分担（定量予測 vs 科学的検証）を明記。**Jupyter MCPで実行したPythonコードを含める**
  - **Results**: 定量的な実験結果を表形式で提示。交差検証の標準偏差を含む。NatureLM予測結果とGALACTICA検証結果の両方を含める。**Jupyter実行結果から得られた数値を \`[cell:<id>]\` 形式で引用すること**
  - **Discussion**: 結果の解釈、限界、先行研究との比較。**NatureLMとGALACTICAの予測の一致・不一致を議論**。**自己の実験の限界・前提条件への依存・実世界への一般化可能性について批判的に議論すること**
  - **Conclusion**: 主要な知見と今後の課題
  - **References**: 先行研究調査で発見した文献をDOI付きで5件以上含める
  - **Reproducibility**: 乱数シード、Pythonバージョン、主要パッケージバージョンを記載

📄 **report.md** — 実験の全結果・手法・考察をまとめたレポート

🔴 **最終確認**: タスク完了前に、paper.md が保存されていることを確認すること。保存されていなければ、今すぐ作成すること。

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
    const prompt = generateRound9Prompt(originalPrompt);
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

  // Wait for completion
  const startTime = Date.now();
  let completed = false;
  let lastStatus = '';
  for (let i = 0; i < Math.ceil(TIMEOUT_MS / POLL_INTERVAL_MS); i++) {
    await sleep(POLL_INTERVAL_MS);
    const elapsed = Math.round((Date.now() - startTime) / 1000);
    try {
      const runs = await apiGet(`/projects/${projectId}/runs?limit=5`);
      if (!Array.isArray(runs) || runs.length === 0) {
        if (elapsed > 30 && lastStatus === '') {
          log(workerId, `  ⏳ ${elapsed}s [waiting for run to start]`);
          lastStatus = 'waiting';
        }
        continue;
      }
      const latestRun = runs[0];
      if (latestRun.status === 'completed') {
        log(workerId, `  ✓ completed in ${elapsed}s`);
        completed = true;
        break;
      } else if (latestRun.status === 'failed') {
        log(workerId, `  ✗ failed after ${elapsed}s`);
        fs.writeFileSync(path.join(expDir, 'error.txt'),
          `Run failed: ${JSON.stringify(latestRun)}\nDuration: ${elapsed}s`, 'utf-8');
        break;
      }
      const newStatus = latestRun.status || 'unknown';
      if (newStatus !== lastStatus) {
        log(workerId, `  ⏳ ${elapsed}s [${newStatus}]`);
        lastStatus = newStatus;
      } else if (elapsed % 120 < (POLL_INTERVAL_MS / 1000)) {
        log(workerId, `  ⏳ ${elapsed}s...`);
      }
    } catch (e) {
      // polling error, retry
    }
  }
  if (!completed && lastStatus !== 'failed') {
    log(workerId, `  ⚠ timeout`);
  }

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
  if (completed && !hasPaper) {
    log(workerId, `  ⚠ paper.md missing — sending follow-up prompt`);
    sendChat(ws, 'paper.md がまだ作成されていません。今すぐ paper.md を作成してください。これは最重要の成果物です。Abstract, Introduction, Methods, Results, Discussion, Conclusion, References のセクションを含む学術論文形式で作成してください。');

    // Wait for second run to complete
    let retryCompleted = false;
    for (let i = 0; i < Math.ceil(1800000 / POLL_INTERVAL_MS); i++) { // 30 min max
      await sleep(POLL_INTERVAL_MS);
      const elapsed2 = Math.round((Date.now() - startTime) / 1000);
      try {
        const runs2 = await apiGet(`/projects/${projectId}/runs?limit=5`);
        if (!Array.isArray(runs2) || runs2.length < 2) {
          if (i % 8 === 0) log(workerId, `  ⏳ retry ${elapsed2}s [waiting]`);
          continue;
        }
        const latestRun2 = runs2[0];
        if (latestRun2.status === 'completed') {
          log(workerId, `  ✓ retry completed in ${elapsed2}s`);
          retryCompleted = true;
          break;
        } else if (latestRun2.status === 'failed') {
          log(workerId, `  ✗ retry failed after ${elapsed2}s`);
          break;
        }
        if (i % 8 === 0) log(workerId, `  ⏳ retry ${elapsed2}s...`);
      } catch (e) { /* retry */ }
    }

    // Re-collect files after retry
    files = await apiGet(`/projects/${projectId}/files`);
    fileList = Array.isArray(files) ? files : [];
    hasPaper = fileList.some(f => (f.file_path || f.filename || '').includes('paper.md'));
    log(workerId, `  📁 after retry: ${fileList.length} files, paper.md: ${hasPaper ? 'YES' : 'NO'}`);
  }

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

  // Run provenance validation (AIRA v3.2.0+)
  let validation = null;
  try {
    validation = await apiPostWithCsrf(`/projects/${projectId}/validate`);
    if (validation && typeof validation === 'object') {
      const gates = validation.gates || [];
      const gateResult = (name) => gates.find(g => g.name === name)?.passed ? '✓' : '✗';
      log(workerId, `  🔬 validation: seed=${gateResult('seed_presence')} env=${gateResult('env_capture')} no-err=${gateResult('no_error_in_cited')} cov=${gateResult('citation_coverage')} claims=${(validation.claims||[]).length} uncited=${(validation.uncited_claims||[]).length}`);
      fs.writeFileSync(path.join(expDir, 'validation.json'), JSON.stringify(validation, null, 2), 'utf-8');
    }
  } catch (e) {
    log(workerId, `  ⚠ validation failed: ${e.message}`);
  }

  // Collect execution trace (AIRA v3.2.0+)
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
    status: completed ? (hasPaper ? 'success' : 'success-no-paper') : (Array.isArray(runs) && runs[0]?.status || 'unknown'),
    galactica_relevant: experiment.galactica_relevant,
    domain: experiment.domain,
    naturelm_used: naturelmUsed,
    galactica_used: galacticaUsed,
    jupyter_used: jupyterUsed,
    validation: validation && typeof validation === 'object' && validation.gates ? {
      seed_presence: validation.gates.find(g => g.name === 'seed_presence')?.passed || false,
      env_capture: validation.gates.find(g => g.name === 'env_capture')?.passed || false,
      no_error_in_cited: validation.gates.find(g => g.name === 'no_error_in_cited')?.passed || false,
      citation_coverage: validation.gates.find(g => g.name === 'citation_coverage')?.passed || false,
      total_claims: (validation.claims || []).length,
      uncited_claims: (validation.uncited_claims || []).length,
      unknown_citations: (validation.unknown_citations || []).length,
      gates_passed: (validation.gates || []).filter(g => g.passed).length
    } : null,
    timestamp: new Date().toISOString()
  };
  fs.writeFileSync(path.join(expDir, 'metadata.json'), JSON.stringify(metadata, null, 2), 'utf-8');

  // Close WebSocket and delete project (keep project if paper.md missing for investigation)
  try { ws.close(); } catch (_) {}
  if (hasPaper) {
    await apiDeleteWithCsrf(`/projects/${projectId}`);
  } else {
    log(workerId, `  ⚠ keeping project ${projectId} for investigation (no paper.md)`);
  }
  log(workerId, `  ✓ done (${duration}s, ${fileList.length} files, NatureLM: ${naturelmUsed ? '✓' : '✗'}, GALACTICA: ${galacticaUsed ? '✓' : '✗'}, Jupyter: ${jupyterUsed ? '✓' : '✗'}, Gates: ${metadata.validation ? metadata.validation.gates_passed + '/4' : 'N/A'})`);

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
    experiments = allExperiments;
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

  console.log(`\n🔬 AIRA Round9 Experiment Runner — Computational Provenance (v3.2.1) + NatureLM + GALACTICA + Jupyter (${NUM_WORKERS} workers)`);
  console.log(`📋 Total: ${experiments.length} | Already done: ${completed.length} | To run: ${pending.length}`);
  console.log(`🧪 GALACTICA-relevant: ${pending.filter(e => e.galactica_relevant).length}/${pending.length}`);
  console.log(`📂 Results → ${RESULTS_DIR}`);
  console.log(`🌐 AIRA: ${BASE_URL}\n`);

  if (pending.length === 0) {
    console.log('All experiments already completed!');
    return;
  }

  fs.mkdirSync(RESULTS_DIR, { recursive: true });

  const queue = [...pending];
  const results = [];
  const totalStart = Date.now();

  // Launch workers
  const workers = [];
  for (let i = 0; i < Math.min(NUM_WORKERS, pending.length); i++) {
    workers.push(worker(i + 1, queue, results));
  }

  await Promise.all(workers);

  const totalDuration = ((Date.now() - totalStart) / 1000 / 60).toFixed(1);
  const successful = results.filter(r => r.status === 'success').length;
  const failed = results.filter(r => r.status !== 'success').length;
  const galacticaUsedCount = results.filter(r => r.galactica_used).length;
  const naturelmUsedCount = results.filter(r => r.naturelm_used).length;
  const jupyterUsedCount = results.filter(r => r.jupyter_used).length;

  // Save summary
  const summary = {
    version: 'round8-naturelm-galactica-jupyter-selfcriticism',
    prompt_source: 'round-9 + NatureLM MCP + GALACTICA MCP + Jupyter MCP + quality prompt + self-criticism instructions',
    total_experiments: experiments.length,
    successful: successful + completed.length,
    failed,
    skipped: completed.length,
    galactica_relevant_count: experiments.filter(e => e.galactica_relevant).length,
    naturelm_actually_used: naturelmUsedCount,
    galactica_actually_used: galacticaUsedCount,
    jupyter_actually_used: jupyterUsedCount,
    total_duration_minutes: parseFloat(totalDuration),
    workers: NUM_WORKERS,
    results
  };
  fs.writeFileSync(path.join(RESULTS_DIR, 'summary.json'), JSON.stringify(summary, null, 2), 'utf-8');

  console.log(`\n${'='.repeat(60)}`);
  console.log(`🏁 ALL DONE | Success: ${successful} | Failed: ${failed} | Skipped: ${completed.length}`);
  console.log(`🧪 NatureLM used: ${naturelmUsedCount}/${successful} | GALACTICA used: ${galacticaUsedCount}/${successful} | Jupyter used: ${jupyterUsedCount}/${successful}`);
  console.log(`   Duration: ${totalDuration} min | Workers: ${NUM_WORKERS}`);
  console.log(`${'='.repeat(60)}\n`);
}

main().catch(err => {
  console.error('Fatal:', err);
  process.exit(1);
});

