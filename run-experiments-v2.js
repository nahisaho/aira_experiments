/**
 * AIRA Experiment Runner v2
 * Uses Playwright UI automation (textarea + Send button) to run experiments
 * on AIRA at localhost:3000
 */

const { chromium } = require('playwright');
const fs = require('fs');
const path = require('path');

const BASE_URL = 'http://localhost:3000';
const RESULTS_DIR = path.join(__dirname, 'results');
const TIMEOUT_MS = 60 * 60 * 1000; // 60 min absolute timeout per experiment
const POLL_INTERVAL_MS = 10000; // 10 sec

// Additional instruction appended to every prompt
const REPORT_INSTRUCTION = `\n\n最後に、この実験の全結果・手法・考察をまとめた report.md を作成してください。report.md には以下を含めること：
- 実験目的と背景
- 使用した手法・アルゴリズムの概要
- 主要な結果と数値
- 考察と今後の展望
- 生成したファイル一覧`;

function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

async function apiGet(page, path) {
  return page.evaluate(async (path) => {
    const res = await fetch(`/api${path}`);
    return res.json();
  }, path);
}

async function apiPost(page, path, body) {
  return page.evaluate(async ({ path, body }) => {
    const tr = await fetch('/api/csrf-token');
    const { token } = await tr.json();
    const h = { 'X-AIRA-Token': token };
    if (body) h['Content-Type'] = 'application/json';
    const res = await fetch(`/api${path}`, {
      method: 'POST', headers: h,
      body: body ? JSON.stringify(body) : undefined
    });
    const text = await res.text();
    try { return JSON.parse(text); } catch { return text; }
  }, { path, body });
}

async function apiDelete(page, path) {
  return page.evaluate(async (path) => {
    const tr = await fetch('/api/csrf-token');
    const { token } = await tr.json();
    await fetch(`/api${path}`, { method: 'DELETE', headers: { 'X-AIRA-Token': token } });
  }, path);
}

async function runExperiment(page, experiment, index, total) {
  const expDir = path.join(RESULTS_DIR, experiment.id);
  fs.mkdirSync(expDir, { recursive: true });

  const prompt = experiment.prompt + REPORT_INSTRUCTION;

  console.log(`\n${'='.repeat(60)}`);
  console.log(`[${index + 1}/${total}] ${experiment.id}: ${experiment.title}`);
  console.log(`${'='.repeat(60)}`);

  // Save input prompt
  fs.writeFileSync(path.join(expDir, 'input_prompt.txt'), prompt, 'utf-8');
  console.log(`  ✓ Input prompt saved (${prompt.length} chars)`);

  // Delete any existing project with this name
  const existingProjects = await apiGet(page, '/projects');
  for (const p of existingProjects.filter(p => p.name === experiment.id)) {
    await apiDelete(page, `/projects/${p.id}`);
  }

  // 1. Create project via UI (click + New, fill name, click Create)
  await page.click('button:has-text("+ New")');
  await page.waitForTimeout(500);
  await page.fill('input[placeholder="Enter project name"]', experiment.id);
  await page.click('button:has-text("Create")');
  await page.waitForTimeout(2000);
  
  // Dismiss any modal that appears after creation
  for (let attempt = 0; attempt < 3; attempt++) {
    const overlay = await page.$('div.fixed.inset-0');
    if (!overlay) break;
    // Try clicking close buttons, X buttons, or pressing Escape
    const closeBtn = await page.$('div.fixed button:has-text("✕")') 
      || await page.$('div.fixed button:has-text("×")') 
      || await page.$('div.fixed button:has-text("Close")');
    if (closeBtn) {
      await closeBtn.click();
    } else {
      await page.keyboard.press('Escape');
    }
    await page.waitForTimeout(500);
  }
  
  // Get project ID from API
  const projects = await apiGet(page, '/projects');
  const project = projects.find(p => p.name === experiment.id);
  if (!project) throw new Error('Project creation failed');
  const projectId = project.id;
  console.log(`  ✓ Project created: ${projectId}`);

  // 2. Assign co-scientist skill via API
  const CO_SCIENTIST_SKILL = '525b1100-7bdf-4cd1-9693-0e3079107206';
  await apiPost(page, `/projects/${projectId}/skills/${CO_SCIENTIST_SKILL}`, null);
  console.log(`  ✓ Co-scientist skill assigned`);

  // Ensure no modal is open before sending
  const modal = await page.$('div.fixed.inset-0');
  if (modal) {
    await page.keyboard.press('Escape');
    await page.waitForTimeout(500);
  }

  // 3. Type prompt and send via UI
  await page.fill('textarea[placeholder="Type a message..."]', prompt);
  await page.click('button:has-text("Send")');
  console.log(`  ✓ Prompt sent via UI`);

  // 4. Wait for run to complete
  const startTime = Date.now();
  let completed = false;

  for (let i = 0; i < Math.ceil(TIMEOUT_MS / POLL_INTERVAL_MS); i++) {
    await sleep(POLL_INTERVAL_MS);
    const elapsed = Math.round((Date.now() - startTime) / 1000);

    const runs = await apiGet(page, `/projects/${projectId}/runs?limit=5`);
    const latestRun = runs[0];

    if (latestRun && latestRun.status === 'completed') {
      console.log(`  ✓ Run completed in ${elapsed}s`);
      completed = true;
      break;
    } else if (latestRun && latestRun.status === 'failed') {
      console.log(`  ✗ Run failed after ${elapsed}s (error: ${latestRun.error_type})`);
      fs.writeFileSync(path.join(expDir, 'error.txt'), 
        `Run failed: ${latestRun.error_type}\nDuration: ${elapsed}s`, 'utf-8');
      break;
    }

    // Print progress every 30s
    if (elapsed % 30 < (POLL_INTERVAL_MS / 1000)) {
      process.stdout.write(`  ⏳ ${elapsed}s elapsed...\r`);
    }
  }

  if (!completed) {
    const elapsed = Math.round((Date.now() - startTime) / 1000);
    console.log(`  ⚠ Timeout after ${elapsed}s`);
  }

  // 6. Collect messages (output)
  const messages = await apiGet(page, `/projects/${projectId}/messages?limit=200&offset=0`);
  const assistantMessages = messages.filter(m => m.role === 'assistant');
  const output = assistantMessages.map(m => m.content).join('\n\n---\n\n');
  fs.writeFileSync(path.join(expDir, 'output_response.md'), output, 'utf-8');
  console.log(`  ✓ Output saved (${output.length} chars, ${assistantMessages.length} messages)`);

  // 7. Collect files
  const files = await apiGet(page, `/projects/${projectId}/files`);
  console.log(`  📁 Generated files: ${files.length}`);

  if (files.length > 0) {
    const filesDir = path.join(expDir, 'files');
    fs.mkdirSync(filesDir, { recursive: true });

    for (const file of files) {
      try {
        const fileName = file.file_path || file.filename || file.id;
        const filePath = path.join(filesDir, fileName);
        fs.mkdirSync(path.dirname(filePath), { recursive: true });
        
        const fileResponse = await page.evaluate(async ({ projectId, fileId }) => {
          const res = await fetch(`/api/projects/${projectId}/files/${fileId}/view`);
          if (!res.ok) return null;
          return res.json();
        }, { projectId, fileId: file.id });
        
        if (fileResponse !== null) {
          const content = fileResponse.content || '';
          fs.writeFileSync(filePath, content, 'utf-8');
        }
      } catch (e) {
        console.log(`  ⚠ Could not download file: ${file.filename || file.id} (${e.message})`);
      }
    }
    console.log(`  ✓ Downloaded ${files.length} files`);

    // Save file list
    fs.writeFileSync(path.join(expDir, 'file_list.json'), JSON.stringify(files, null, 2), 'utf-8');

    // Try zip download  
    try {
      const zipResult = await page.evaluate(async (projectId) => {
        const res = await fetch(`/api/projects/${projectId}/files/zip`);
        if (!res.ok) return null;
        const blob = await res.blob();
        const ab = await blob.arrayBuffer();
        return Array.from(new Uint8Array(ab));
      }, projectId);
      if (zipResult) {
        fs.writeFileSync(path.join(expDir, 'files.zip'), Buffer.from(zipResult));
        console.log(`  ✓ files.zip saved`);
      }
    } catch (e) { /* zip not available */ }
  }

  // 8. Save metadata
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
    runs: runs,
    timestamp: new Date().toISOString()
  };
  fs.writeFileSync(path.join(expDir, 'metadata.json'), JSON.stringify(metadata, null, 2), 'utf-8');
  console.log(`  ✓ Metadata saved`);

  // 9. Delete project
  await apiDelete(page, `/projects/${projectId}`);
  console.log(`  ✓ Project deleted`);

  return metadata;
}

async function main() {
  const experimentsPath = path.join(__dirname, 'experiments', 'experiment_prompts.json');
  let experiments = JSON.parse(fs.readFileSync(experimentsPath, 'utf-8'));

  // Filter by experiment ID if provided as CLI argument
  const filterIds = process.argv.slice(2);
  if (filterIds.length > 0) {
    experiments = experiments.filter(e => filterIds.includes(e.id));
  }

  console.log(`\n🔬 AIRA Experiment Runner v2`);
  console.log(`📋 ${experiments.length} experiments to run`);
  console.log(`🎯 Each prompt includes report.md generation instruction`);
  console.log(`📂 Results → ${RESULTS_DIR}\n`);

  fs.mkdirSync(RESULTS_DIR, { recursive: true });

  const browser = await chromium.launch({ headless: true });
  const page = await browser.newPage();
  await page.goto(BASE_URL);
  await page.waitForLoadState('networkidle');
  console.log(`✓ Connected to AIRA at ${BASE_URL}\n`);

  const results = [];
  const totalStart = Date.now();

  for (let i = 0; i < experiments.length; i++) {
    try {
      const metadata = await runExperiment(page, experiments[i], i, experiments.length);
      results.push(metadata);
    } catch (err) {
      console.log(`  ✗ FATAL ERROR: ${err.message}`);
      results.push({ 
        experiment_id: experiments[i].id, 
        status: 'error', 
        error: err.message,
        timestamp: new Date().toISOString()
      });
    }

    // Brief pause between experiments
    if (i < experiments.length - 1) {
      await sleep(5000);
      // Refresh page to reset state
      await page.goto(BASE_URL);
      await page.waitForLoadState('networkidle');
    }
  }

  const totalDuration = ((Date.now() - totalStart) / 1000 / 60).toFixed(1);

  // Save summary
  const summary = {
    total_experiments: experiments.length,
    successful: results.filter(r => r.status === 'success').length,
    failed: results.filter(r => r.status !== 'success').length,
    total_duration_minutes: parseFloat(totalDuration),
    results
  };
  fs.writeFileSync(path.join(RESULTS_DIR, 'summary.json'), JSON.stringify(summary, null, 2), 'utf-8');

  console.log(`\n${'='.repeat(60)}`);
  console.log(`🏁 ALL EXPERIMENTS COMPLETE`);
  console.log(`   Total: ${experiments.length} | Success: ${summary.successful} | Failed: ${summary.failed}`);
  console.log(`   Duration: ${totalDuration} minutes`);
  console.log(`   Results: ${RESULTS_DIR}`);
  console.log(`${'='.repeat(60)}\n`);

  await browser.close();
}

main().catch(err => {
  console.error('Fatal error:', err);
  process.exit(1);
});
