#!/usr/bin/env node
/**
 * Extract all experiment prompts from the Qiita article
 * Fetches all pages and parses SCI-XXX prompts
 */
const fs = require('fs');
const path = require('path');

async function fetchPage(startIndex = 0) {
  const url = `https://qiita.com/hisaho/items/ad00285df52fd84494f0`;
  // Use the Qiita API to get the raw markdown
  const res = await fetch(`https://qiita.com/api/v2/items/ad00285df52fd84494f0`, {
    headers: { 'Accept': 'application/json' }
  });
  const data = await res.json();
  return data.body; // raw markdown
}

async function main() {
  console.log('Fetching article from Qiita API...');
  const markdown = await fetchPage();
  console.log(`Article length: ${markdown.length} chars`);
  
  // Save raw markdown for reference
  fs.writeFileSync('/tmp/qiita_article_raw.md', markdown, 'utf-8');
  
  // Parse all SCI-XXX experiments
  // Pattern: ## SCI-XXX: Title (with optional ✅)
  // Then find the prompt in ``` code blocks after "プロンプト"
  
  const experiments = [];
  const lines = markdown.split('\n');
  
  let i = 0;
  while (i < lines.length) {
    const line = lines[i];
    // Match experiment headers like: ## SCI-001: CRISPR-Cas9オフターゲット予測モデル ✅
    const headerMatch = line.match(/^##\s+(SCI-(\d{3}))[:：]\s*(.+?)(?:\s*[✅❌⚠️])?\s*$/);
    
    if (headerMatch) {
      const id = headerMatch[1];
      const num = parseInt(headerMatch[2]);
      const title = headerMatch[3].trim();
      
      // Search for the prompt code block
      let prompt = '';
      let j = i + 1;
      let foundPromptSection = false;
      
      while (j < lines.length) {
        // Check if we hit the next experiment header
        if (lines[j].match(/^##\s+SCI-\d{3}/)) break;
        
        // Look for "プロンプト" indicator
        if (lines[j].includes('プロンプト') || lines[j].includes('prompt')) {
          foundPromptSection = true;
        }
        
        // If we found the prompt section, look for the code block
        if (foundPromptSection && lines[j].trim() === '```') {
          // Start of code block - collect until closing ```
          j++;
          const promptLines = [];
          while (j < lines.length && lines[j].trim() !== '```') {
            promptLines.push(lines[j]);
            j++;
          }
          if (promptLines.length > 0) {
            prompt = promptLines.join('\n').trim();
            break;
          }
        }
        j++;
      }
      
      if (prompt) {
        experiments.push({ id, title, prompt });
        console.log(`  ✓ ${id}: ${title} (${prompt.length} chars)`);
      } else {
        console.log(`  ⚠ ${id}: ${title} - NO PROMPT FOUND`);
      }
    }
    i++;
  }
  
  console.log(`\nTotal experiments found: ${experiments.length}`);
  
  // Save to experiment_prompts.json
  const outPath = path.join(__dirname, 'experiments', 'experiment_prompts.json');
  fs.writeFileSync(outPath, JSON.stringify(experiments, null, 2), 'utf-8');
  console.log(`Saved to ${outPath}`);
}

main().catch(console.error);
