/**
 * AIRA Experiment Runner
 * Playwright automation to run scientific experiments on AIRA (localhost:3000)
 * - Creates project per experiment
 * - Assigns co-scientist skill
 * - Sends prompt with report.md generation instruction
 * - Waits for completion
 * - Downloads generated files as zip
 * - Saves prompt and output per experiment directory
 */

const { chromium } = require('playwright');
const fs = require('fs');
const path = require('path');

const BASE_URL = 'http://localhost:3000';
const CO_SCIENTIST_SKILL_ID = '525b1100-7bdf-4cd1-9693-0e3079107206';
const RESULTS_DIR = path.join(__dirname, 'results');
const TIMEOUT_MS = 60 * 60 * 1000; // 60 min absolute timeout
const IDLE_TIMEOUT_MS = 30 * 60 * 1000; // 30 min idle timeout
const POLL_INTERVAL_MS = 10000; // 10 sec polling

// Additional instruction to append to every prompt
const REPORT_INSTRUCTION = `\n\n最後に、この実験の全結果・手法・考察をまとめた report.md を作成してください。report.md には以下を含めること：
- 実験目的と背景
- 使用した手法・アルゴリズムの概要
- 主要な結果と数値
- 考察と今後の展望
- 生成したファイル一覧`;

async function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

async function apiCall(page, method, path, body = null) {
  return await page.evaluate(async ({ method, path, body }) => {
    const tokenRes = await fetch('/api/csrf-token');
    const { token } = await tokenRes.json();
    const headers = { 'X-AIRA-Token': token };
    if (body) headers['Content-Type'] = 'application/json';
    const res = await fetch(`/api${path}`, {
      method,
      headers,
      body: body ? JSON.stringify(body) : undefined
    });
    const text = await res.text();
    try { return JSON.parse(text); } catch { return text; }
  }, { method, path, body });
}

async function createProject(page, name) {
  return await apiCall(page, 'POST', '/projects', { name });
}

async function assignSkill(page, projectId, skillId) {
  await apiCall(page, 'POST', `/projects/${projectId}/skills/${skillId}`);
}

async function sendMessage(page, projectId, content) {
  // First create the message via API
  const msg = await apiCall(page, 'POST', `/projects/${projectId}/messages`, { content });
  
  // Then trigger the run via WebSocket
  await page.evaluate(async ({ projectId, content, messageId }) => {
    return new Promise((resolve, reject) => {
      const ws = new WebSocket(`ws://localhost:3000/ws/projects/${projectId}/chat`);
      ws.onopen = () => {
        ws.send(JSON.stringify({ type: 'chat', content, messageId }));
        // Give it a moment then close our connection (the server keeps the run going)
        setTimeout(() => { ws.close(); resolve(); }, 1000);
      };
      ws.onerror = (e) => reject(new Error('WebSocket error'));
      setTimeout(() => reject(new Error('WebSocket timeout')), 10000);
    });
  }, { projectId, content, messageId: msg.id });
  
  return msg;
}

async function waitForCompletion(page, projectId, timeoutMs = TIMEOUT_MS) {
  const startTime = Date.now();

  // First wait for a run to start (status becomes 'running')
  let runStarted = false;
  while (!runStarted) {
    if (Date.now() - startTime > 60000) {
      // If no run starts within 60s, check if there's a completed run
      const runs = await page.evaluate(async (pid) => {
        const res = await fetch(`/api/projects/${pid}/runs?limit=5`);
        return res.json();
      }, projectId);
      if (runs.length > 0 && runs[0].status === 'completed') return runs[0];
      throw new Error('Run did not start within 60s');
    }

    const current = await page.evaluate(async (pid) => {
      const res = await fetch(`/api/projects/${pid}/runs/current`);
      return res.json();
    }, projectId);

    if (current.status === 'running') {
      runStarted = true;
      break;
    }

    // Also check runs list for newly started run
    const runs = await page.evaluate(async (pid) => {
      const res = await fetch(`/api/projects/${pid}/runs?limit=5`);
      return res.json();
    }, projectId);
    if (runs.length > 0 && (runs[0].status === 'running' || runs[0].status === 'completed')) {
      runStarted = true;
      if (runs[0].status === 'completed') return runs[0];
      break;
    }

    await sleep(2000);
  }

  // Now wait for completion
  while (true) {
    const elapsed = Date.now() - startTime;
    if (elapsed > timeoutMs) {
      throw new Error(`Absolute timeout (${timeoutMs / 60000} min) exceeded`);
    }

    const current = await page.evaluate(async (pid) => {
      const res = await fetch(`/api/projects/${pid}/runs/current`);
      return res.json();
    }, projectId);

    if (current.status === 'idle') {
      // Double-check that the run actually completed
      const runs = await page.evaluate(async (pid) => {
        const res = await fetch(`/api/projects/${pid}/runs?limit=5`);
        return res.json();
      }, projectId);
      if (runs.length > 0 && runs[0].status === 'completed') {
        return runs[0];
      }
    }

    await sleep(POLL_INTERVAL_MS);
  }
}

async function getMessages(page, projectId) {
  const messages = await page.evaluate(async (projectId) => {
    const res = await fetch(`/api/projects/${projectId}/messages?limit=200&offset=0`);
    return res.json();
  }, projectId);
  return messages;
}

async function getFiles(page, projectId) {
  const files = await page.evaluate(async (projectId) => {
    const res = await fetch(`/api/projects/${projectId}/files`);
    return res.json();
  }, projectId);
  return files;
}

async function downloadFile(page, projectId, fileId) {
  const content = await page.evaluate(async ({ projectId, fileId }) => {
    const res = await fetch(`/api/projects/${projectId}/files/${fileId}/view`);
    if (!res.ok) return null;
    return res.text();
  }, { projectId, fileId });
  return content;
}

async function downloadFilesAsZip(page, projectId) {
  // Try zip endpoint
  const result = await page.evaluate(async (projectId) => {
    const res = await fetch(`/api/projects/${projectId}/files/zip`);
    if (!res.ok) return { ok: false, status: res.status };
    const blob = await res.blob();
    const arrayBuffer = await blob.arrayBuffer();
    return { ok: true, data: Array.from(new Uint8Array(arrayBuffer)) };
  }, projectId);
  if (result.ok) return result.data;
  return null;
}

async function getRuns(page, projectId) {
  const runs = await page.evaluate(async (projectId) => {
    const res = await fetch(`/api/projects/${projectId}/runs?limit=20`);
    return res.json();
  }, projectId);
  return runs;
}

async function deleteProject(page, projectId) {
  await apiCall(page, 'DELETE', `/projects/${projectId}`);
}

async function runExperiment(page, experiment, index) {
  const expDir = path.join(RESULTS_DIR, experiment.id);
  fs.mkdirSync(expDir, { recursive: true });

  const prompt = experiment.prompt + REPORT_INSTRUCTION;

  console.log(`\n${'='.repeat(60)}`);
  console.log(`[${index + 1}/7] ${experiment.id}: ${experiment.title}`);
  console.log(`${'='.repeat(60)}`);

  // Save input prompt
  fs.writeFileSync(path.join(expDir, 'input_prompt.txt'), prompt, 'utf-8');
  console.log(`  ✓ Saved input prompt`);

  // Create project
  const project = await createProject(page, experiment.id);
  const projectId = project.id;
  console.log(`  ✓ Created project: ${projectId}`);

  // Assign co-scientist skill
  await assignSkill(page, projectId, CO_SCIENTIST_SKILL_ID);
  console.log(`  ✓ Assigned co-scientist skill`);

  // Send prompt
  const startTime = Date.now();
  await sendMessage(page, projectId, prompt);
  console.log(`  ✓ Sent prompt (${prompt.length} chars)`);
  console.log(`  ⏳ Waiting for completion...`);

  // Wait for completion
  try {
    await waitForCompletion(page, projectId);
    const duration = ((Date.now() - startTime) / 1000).toFixed(1);
    console.log(`  ✓ Completed in ${duration}s`);
  } catch (err) {
    console.log(`  ✗ ${err.message}`);
    const duration = ((Date.now() - startTime) / 1000).toFixed(1);
    fs.writeFileSync(path.join(expDir, 'error.txt'), `${err.message}\nDuration: ${duration}s`, 'utf-8');
  }

  // Get messages (output)
  const messages = await getMessages(page, projectId);
  const assistantMessages = messages.filter(m => m.role === 'assistant');
  const output = assistantMessages.map(m => m.content).join('\n\n---\n\n');
  fs.writeFileSync(path.join(expDir, 'output_response.md'), output, 'utf-8');
  console.log(`  ✓ Saved output (${output.length} chars, ${assistantMessages.length} messages)`);

  // Get and download files
  const files = await getFiles(page, projectId);
  console.log(`  📁 Generated files: ${files.length}`);

  if (files.length > 0) {
    // Try zip download first
    const zipData = await downloadFilesAsZip(page, projectId);
    if (zipData) {
      fs.writeFileSync(path.join(expDir, 'files.zip'), Buffer.from(zipData));
      console.log(`  ✓ Downloaded files.zip`);
    } else {
      // Download files individually
      const filesDir = path.join(expDir, 'files');
      fs.mkdirSync(filesDir, { recursive: true });
      for (const file of files) {
        const fileId = file.id || file.path || file.name;
        const content = await downloadFile(page, projectId, fileId);
        if (content !== null) {
          const fileName = file.name || file.path || fileId;
          const filePath = path.join(filesDir, fileName);
          fs.mkdirSync(path.dirname(filePath), { recursive: true });
          fs.writeFileSync(filePath, content, 'utf-8');
        }
      }
      console.log(`  ✓ Downloaded ${files.length} files individually`);
    }

    // Save file list
    fs.writeFileSync(
      path.join(expDir, 'file_list.json'),
      JSON.stringify(files, null, 2),
      'utf-8'
    );
  }

  // Save run metadata
  const runs = await getRuns(page, projectId);
  const metadata = {
    experiment_id: experiment.id,
    title: experiment.title,
    project_id: projectId,
    prompt_length: prompt.length,
    output_length: output.length,
    file_count: files.length,
    runs: runs,
    timestamp: new Date().toISOString()
  };
  fs.writeFileSync(path.join(expDir, 'metadata.json'), JSON.stringify(metadata, null, 2), 'utf-8');
  console.log(`  ✓ Saved metadata`);

  // Delete project to clean up
  await deleteProject(page, projectId);
  console.log(`  ✓ Deleted project`);

  return metadata;
}

async function main() {
  // Load experiments
  const experimentsPath = path.join(__dirname, 'experiments', 'experiment_prompts.json');
  const experiments = JSON.parse(fs.readFileSync(experimentsPath, 'utf-8'));

  console.log(`\n🔬 AIRA Experiment Runner`);
  console.log(`📋 ${experiments.length} experiments to run`);
  console.log(`🎯 Each prompt includes report.md generation instruction`);
  console.log(`📂 Results will be saved to: ${RESULTS_DIR}\n`);

  fs.mkdirSync(RESULTS_DIR, { recursive: true });

  const browser = await chromium.launch({ headless: true });
  const page = await browser.newPage();

  // Navigate to AIRA and establish session
  await page.goto(BASE_URL);
  await page.waitForLoadState('networkidle');
  console.log(`✓ Connected to AIRA at ${BASE_URL}\n`);

  const results = [];
  const startTime = Date.now();

  for (let i = 0; i < experiments.length; i++) {
    try {
      const metadata = await runExperiment(page, experiments[i], i);
      results.push({ ...metadata, status: 'success' });
    } catch (err) {
      console.log(`  ✗ EXPERIMENT FAILED: ${err.message}`);
      results.push({ experiment_id: experiments[i].id, status: 'failed', error: err.message });
    }

    // Brief pause between experiments
    if (i < experiments.length - 1) {
      await sleep(3000);
    }
  }

  const totalDuration = ((Date.now() - startTime) / 1000 / 60).toFixed(1);

  // Save summary
  const summary = {
    total_experiments: experiments.length,
    successful: results.filter(r => r.status === 'success').length,
    failed: results.filter(r => r.status === 'failed').length,
    total_duration_minutes: totalDuration,
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
