/**
 * AIRA Experiment Runner v4.0 - Uses round-9 prompts (with 先行研究調査 prefix)
 * 5 Parallel Workers, reads prompts from co-scientist-optimization/round-9/
 */

const { chromium } = require('playwright');
const fs = require('fs');
const path = require('path');

const BASE_URL = 'http://localhost:5173';
const RESULTS_DIR = path.join(__dirname, 'results-v4.0');
const ROUND9_DIR = path.join(__dirname, '..', 'aira', 'co-scientist-optimization', 'round-9');
const TIMEOUT_MS = 60 * 60 * 1000;
const POLL_INTERVAL_MS = 10000;
const NUM_WORKERS = 5;

const CO_SCIENTIST_SKILL = '525b1100-7bdf-4cd1-9693-0e3079107206';

function sleep(ms) { return new Promise(r => setTimeout(r, ms)); }

function log(workerId, msg) {
  const ts = new Date().toISOString().slice(11, 19);
  console.log(`[${ts}][W${workerId}] ${msg}`);
}

function loadRound9Experiments() {
  const dirs = fs.readdirSync(ROUND9_DIR)
    .filter(d => d.startsWith('SCI-'))
    .sort();
  
  return dirs.map(id => {
    const promptPath = path.join(ROUND9_DIR, id, 'input_prompt.txt');
    const prompt = fs.readFileSync(promptPath, 'utf-8');
    // Extract title from the research theme (after "---" separator)
    const afterSep = prompt.split('---\n')[1] || prompt;
    const titleMatch = afterSep.match(/^(.+?)(?:\n|。)/);
    const title = titleMatch ? titleMatch[1].trim().slice(0, 60) : id;
    return { id, title, prompt };
  });
}

async function apiGet(page, apiPath) {
  return page.evaluate(async (p) => {
    const res = await fetch(`/api${p}`);
    return res.json();
  }, apiPath);
}

async function apiPost(page, apiPath, body) {
  return page.evaluate(async ({ p, b }) => {
    const tr = await fetch('/api/csrf-token');
    const { token } = await tr.json();
    const h = { 'X-AIRA-Token': token };
    if (b) h['Content-Type'] = 'application/json';
    const res = await fetch(`/api${p}`, {
      method: 'POST', headers: h,
      body: b ? JSON.stringify(b) : undefined
    });
    const text = await res.text();
    try { return JSON.parse(text); } catch { return text; }
  }, { p: apiPath, b: body });
}

async function apiDelete(page, apiPath) {
  return page.evaluate(async (p) => {
    const tr = await fetch('/api/csrf-token');
    const { token } = await tr.json();
    await fetch(`/api${p}`, { method: 'DELETE', headers: { 'X-AIRA-Token': token } });
  }, apiPath);
}

async function runExperiment(page, experiment, workerId) {
  const expDir = path.join(RESULTS_DIR, experiment.id);
  fs.mkdirSync(expDir, { recursive: true });
  // round-9 prompts already include report/paper instructions
  const prompt = experiment.prompt;

  log(workerId, `▶ ${experiment.id}: ${experiment.title}`);

  fs.writeFileSync(path.join(expDir, 'input_prompt.txt'), prompt, 'utf-8');

  // Delete any existing project with this name
  const existing = await apiGet(page, '/projects');
  for (const p of existing.filter(p => p.name === experiment.id)) {
    await apiDelete(page, `/projects/${p.id}`);
  }

  // Navigate to fresh state
  await page.goto(BASE_URL);
  await page.waitForLoadState('networkidle');
  await sleep(1000);

  // Create project via UI
  const newBtn = await page.$('button:has-text("+ New")') || await page.$('button:has-text("+ 新規")');
  if (!newBtn) throw new Error('Cannot find New Project button');
  await newBtn.click();
  await page.waitForTimeout(500);
  const nameInput = await page.$('input[placeholder="Enter project name"]') || await page.$('input[placeholder="プロジェクト名を入力"]');
  if (nameInput) await nameInput.fill(experiment.id);
  else await page.fill('input[type="text"]', experiment.id);
  const createBtn = await page.$('button:has-text("Create")') || await page.$('button:has-text("作成")');
  if (createBtn) await createBtn.click();
  else await page.click('button[type="submit"]');
  await page.waitForTimeout(2000);

  // Dismiss modals
  for (let attempt = 0; attempt < 3; attempt++) {
    const overlay = await page.$('div.fixed.inset-0');
    if (!overlay) break;
    const closeBtn = await page.$('div.fixed button:has-text("✕")')
      || await page.$('div.fixed button:has-text("×")')
      || await page.$('div.fixed button:has-text("Close")');
    if (closeBtn) await closeBtn.click();
    else await page.keyboard.press('Escape');
    await page.waitForTimeout(500);
  }

  // Get project ID
  const projects = await apiGet(page, '/projects');
  const project = projects.find(p => p.name === experiment.id);
  if (!project) throw new Error('Project creation failed');
  const projectId = project.id;
  log(workerId, `  project: ${projectId}`);

  // Assign skill via API
  await apiPost(page, `/projects/${projectId}/skills/${CO_SCIENTIST_SKILL}`, null);

  // Ensure no modal
  const modal = await page.$('div.fixed.inset-0');
  if (modal) {
    await page.keyboard.press('Escape');
    await page.waitForTimeout(500);
  }

  // Send prompt via UI
  const textarea = await page.$('textarea[placeholder="Type a message..."]') || await page.$('textarea[placeholder="メッセージを入力..."]') || await page.$('textarea');
  if (textarea) await textarea.fill(prompt);
  else throw new Error('Cannot find textarea');
  const sendBtn = await page.$('button:has-text("Send")') || await page.$('button:has-text("送信")');
  if (sendBtn) await sendBtn.click();
  else await page.keyboard.press('Enter');
  log(workerId, `  prompt sent`);

  // Wait for completion
  const startTime = Date.now();
  let completed = false;
  for (let i = 0; i < Math.ceil(TIMEOUT_MS / POLL_INTERVAL_MS); i++) {
    await sleep(POLL_INTERVAL_MS);
    const elapsed = Math.round((Date.now() - startTime) / 1000);
    try {
      const runs = await apiGet(page, `/projects/${projectId}/runs?limit=5`);
      const latestRun = runs[0];
      if (latestRun && latestRun.status === 'completed') {
        log(workerId, `  ✓ completed in ${elapsed}s`);
        completed = true;
        break;
      } else if (latestRun && latestRun.status === 'failed') {
        log(workerId, `  ✗ failed after ${elapsed}s`);
        fs.writeFileSync(path.join(expDir, 'error.txt'),
          `Run failed: ${latestRun.error_type}\nDuration: ${elapsed}s`, 'utf-8');
        break;
      }
      if (elapsed % 60 < (POLL_INTERVAL_MS / 1000)) {
        log(workerId, `  ⏳ ${elapsed}s...`);
      }
    } catch (e) {
      // polling error, retry
    }
  }
  if (!completed) {
    log(workerId, `  ⚠ timeout`);
  }

  // Collect messages
  const messages = await apiGet(page, `/projects/${projectId}/messages?limit=200&offset=0`);
  const assistantMessages = messages.filter(m => m.role === 'assistant');
  const output = assistantMessages.map(m => m.content).join('\n\n---\n\n');
  fs.writeFileSync(path.join(expDir, 'output_response.md'), output, 'utf-8');

  // Collect files
  const files = await apiGet(page, `/projects/${projectId}/files`);
  log(workerId, `  📁 ${files.length} files`);

  if (files.length > 0) {
    const filesDir = path.join(expDir, 'files');
    fs.mkdirSync(filesDir, { recursive: true });

    for (const file of files) {
      try {
        const fileName = file.file_path || file.filename || file.id;
        const filePath = path.join(filesDir, fileName);
        fs.mkdirSync(path.dirname(filePath), { recursive: true });
        const fileBuffer = await page.evaluate(async ({ projectId, fileId }) => {
          const res = await fetch(`/api/projects/${projectId}/files/${fileId}/download`);
          if (!res.ok) return null;
          const buf = await res.arrayBuffer();
          return Array.from(new Uint8Array(buf));
        }, { projectId, fileId: file.id });
        if (fileBuffer !== null) {
          fs.writeFileSync(filePath, Buffer.from(fileBuffer));
        }
      } catch (e) {
        // skip failed file
      }
    }
    fs.writeFileSync(path.join(expDir, 'file_list.json'), JSON.stringify(files, null, 2), 'utf-8');
  }

  // Save metadata
  const runs = await apiGet(page, `/projects/${projectId}/runs?limit=5`);
  const duration = Math.round((Date.now() - startTime) / 1000);
  const metadata = {
    experiment_id: experiment.id,
    title: experiment.title,
    project_id: projectId,
    prompt_length: prompt.length,
    output_length: output.length,
    file_count: files.length,
    duration_seconds: duration,
    status: completed ? 'success' : (runs[0]?.status || 'unknown'),
    timestamp: new Date().toISOString()
  };
  fs.writeFileSync(path.join(expDir, 'metadata.json'), JSON.stringify(metadata, null, 2), 'utf-8');

  // Delete project
  await apiDelete(page, `/projects/${projectId}`);
  log(workerId, `  ✓ done (${duration}s, ${files.length} files)`);

  return metadata;
}

async function worker(browser, workerId, queue, results) {
  const context = await browser.newContext();
  const page = await context.newPage();
  await page.goto(BASE_URL);
  await page.waitForLoadState('networkidle');
  log(workerId, 'ready');

  while (true) {
    const experiment = queue.shift();
    if (!experiment) break;

    try {
      const meta = await runExperiment(page, experiment, workerId);
      results.push(meta);
    } catch (err) {
      log(workerId, `  ✗ ERROR: ${err.message}`);
      results.push({
        experiment_id: experiment.id,
        status: 'error',
        error: err.message,
        timestamp: new Date().toISOString()
      });
      try {
        await page.goto(BASE_URL);
        await page.waitForLoadState('networkidle');
      } catch (_) {}
    }
    await sleep(2000);
  }

  await context.close();
  log(workerId, 'finished');
}

async function main() {
  const allExperiments = loadRound9Experiments();
  console.log(`📂 Loaded ${allExperiments.length} experiments from round-9`);

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

  console.log(`\n🔬 AIRA v4.0 Experiment Runner (round-9 prompts, ${NUM_WORKERS} workers)`);
  console.log(`📋 Total: ${experiments.length} | Already done: ${completed.length} | To run: ${pending.length}`);
  console.log(`📂 Results → ${RESULTS_DIR}\n`);

  if (pending.length === 0) {
    console.log('All experiments already completed!');
    return;
  }

  fs.mkdirSync(RESULTS_DIR, { recursive: true });

  const browser = await chromium.launch({ headless: true });
  const queue = [...pending];
  const results = [];
  const totalStart = Date.now();

  // Launch workers
  const workers = [];
  for (let i = 0; i < Math.min(NUM_WORKERS, pending.length); i++) {
    workers.push(worker(browser, i + 1, queue, results));
  }

  await Promise.all(workers);

  const totalDuration = ((Date.now() - totalStart) / 1000 / 60).toFixed(1);
  const successful = results.filter(r => r.status === 'success').length;
  const failed = results.filter(r => r.status !== 'success').length;

  // Save summary
  const summary = {
    version: 'v4.0',
    prompt_source: 'round-9 (with 先行研究調査 prefix)',
    total_experiments: experiments.length,
    successful: successful + completed.length,
    failed,
    skipped: completed.length,
    total_duration_minutes: parseFloat(totalDuration),
    workers: NUM_WORKERS,
    results
  };
  fs.writeFileSync(path.join(RESULTS_DIR, 'summary.json'), JSON.stringify(summary, null, 2), 'utf-8');

  console.log(`\n${'='.repeat(60)}`);
  console.log(`🏁 ALL DONE | Success: ${successful} | Failed: ${failed} | Skipped: ${completed.length}`);
  console.log(`   Duration: ${totalDuration} min | Workers: ${NUM_WORKERS}`);
  console.log(`${'='.repeat(60)}\n`);

  await browser.close();
}

main().catch(err => {
  console.error('Fatal:', err);
  process.exit(1);
});
