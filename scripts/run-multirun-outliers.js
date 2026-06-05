#!/usr/bin/env node
/**
 * Multi-run experiment for VM outlier analysis
 * Runs 4 outlier experiments (SCI-004, SCI-018, SCI-059, SCI-092) × 3 times each
 * to measure LLM non-determinism impact on VM.
 * Results stored in results/round14-multirun/run-{1,2,3}/SCI-XXX/
 * 
 * Reuses run-experiments-round14.js by setting RESULTS_DIR per run.
 */

const { execSync, spawn } = require('child_process');
const fs = require('fs');
const path = require('path');

const OUTLIERS = ['SCI-004', 'SCI-018', 'SCI-059', 'SCI-092'];
const NUM_RUNS = 3;
const SCRIPT = path.join(__dirname, 'run-experiments-round14.js');
const BASE_RESULTS = path.join(__dirname, '..', 'results', 'round14-multirun');

async function runBatch(runNum) {
  const resultsDir = path.join(BASE_RESULTS, `run-${runNum}`);
  fs.mkdirSync(resultsDir, { recursive: true });

  console.log(`\n${'='.repeat(60)}`);
  console.log(`🔄 Run ${runNum}/${NUM_RUNS} — ${OUTLIERS.join(', ')}`);
  console.log(`📂 Results → ${resultsDir}`);
  console.log(`${'='.repeat(60)}\n`);

  return new Promise((resolve, reject) => {
    const child = spawn('node', [SCRIPT, ...OUTLIERS], {
      env: {
        ...process.env,
        RESULTS_DIR_OVERRIDE: resultsDir,
        WORKERS: '2',  // 2 workers for 4 experiments
      },
      stdio: 'inherit',
    });

    child.on('close', (code) => {
      if (code === 0) resolve();
      else reject(new Error(`Run ${runNum} exited with code ${code}`));
    });
  });
}

async function analyzeResults() {
  console.log(`\n${'='.repeat(60)}`);
  console.log('📊 MULTI-RUN ANALYSIS');
  console.log(`${'='.repeat(60)}\n`);

  const table = [];
  for (const id of OUTLIERS) {
    const row = { id, vms: [], unciteds: [], figOrps: [], claims: [], durations: [] };
    for (let r = 1; r <= NUM_RUNS; r++) {
      try {
        const meta = JSON.parse(fs.readFileSync(
          path.join(BASE_RESULTS, `run-${r}`, id, 'metadata.json'), 'utf8'
        ));
        const v = meta.validation || {};
        row.vms.push(v.value_mismatches || 0);
        row.unciteds.push(v.uncited_claims || 0);
        row.figOrps.push(v.figure_orphans || 0);
        row.claims.push(v.total_claims || 0);
        row.durations.push(meta.duration_seconds || 0);
      } catch (e) {
        row.vms.push('?');
      }
    }
    table.push(row);
  }

  const median = arr => {
    const sorted = [...arr].filter(v => typeof v === 'number').sort((a, b) => a - b);
    const mid = Math.floor(sorted.length / 2);
    return sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
  };
  const stddev = arr => {
    const nums = arr.filter(v => typeof v === 'number');
    const avg = nums.reduce((s, v) => s + v, 0) / nums.length;
    return Math.sqrt(nums.reduce((s, v) => s + (v - avg) ** 2, 0) / nums.length);
  };

  console.log('ID        | Run1 VM | Run2 VM | Run3 VM | Median | StdDev | R14-orig | Claims(runs)');
  console.log('----------|---------|---------|---------|--------|--------|----------|-------------');
  for (const row of table) {
    // Get original R14 VM
    let origVM = '?';
    try {
      const meta = JSON.parse(fs.readFileSync(
        path.join(__dirname, '..', 'results', 'round14', row.id, 'metadata.json'), 'utf8'
      ));
      origVM = meta.validation?.value_mismatches || 0;
    } catch (e) {}

    const med = median(row.vms);
    const sd = stddev(row.vms);
    console.log(
      `${row.id.padEnd(10)}| ${String(row.vms[0]).padEnd(8)}| ${String(row.vms[1]).padEnd(8)}| ${String(row.vms[2]).padEnd(8)}| ${String(med).padEnd(7)}| ${sd.toFixed(1).padEnd(7)}| ${String(origVM).padEnd(9)}| ${row.claims.join(',')}`
    );
  }

  // Save analysis
  const analysis = { outliers: OUTLIERS, num_runs: NUM_RUNS, results: table };
  fs.writeFileSync(
    path.join(BASE_RESULTS, 'analysis.json'),
    JSON.stringify(analysis, null, 2), 'utf8'
  );
  console.log(`\n📝 Analysis saved to ${path.join(BASE_RESULTS, 'analysis.json')}`);
}

async function main() {
  console.log('🔬 Multi-run VM outlier analysis');
  console.log(`   Experiments: ${OUTLIERS.join(', ')}`);
  console.log(`   Runs: ${NUM_RUNS}`);
  console.log(`   Total: ${OUTLIERS.length * NUM_RUNS} experiment executions\n`);

  for (let r = 1; r <= NUM_RUNS; r++) {
    await runBatch(r);
  }

  await analyzeResults();
}

main().catch(err => {
  console.error('Fatal:', err);
  process.exit(1);
});
