const { chromium } = require('playwright');
(async () => {
  const browser = await chromium.launch({ headless: true });
  const page = await browser.newPage();
  await page.goto('http://localhost:3000');
  await page.waitForLoadState('networkidle');
  
  // Create project
  await page.click('button:has-text("+ New")');
  await page.waitForTimeout(500);
  await page.fill('input[placeholder="Enter project name"]', 'TestExp');
  await page.click('button:has-text("Create")');
  await page.waitForTimeout(2000);
  
  // Check what's on screen
  const overlays = await page.$$('div.fixed.inset-0');
  console.log('Overlays:', overlays.length);
  
  const textarea = await page.$('textarea');
  console.log('Textarea found:', !!textarea);
  
  const buttons = await page.$$eval('button', els => els.map(e => e.textContent.trim().substring(0, 40)));
  console.log('Buttons:', buttons);
  
  if (overlays.length > 0) {
    const overlayText = await overlays[0].evaluate(el => el.textContent.substring(0, 500));
    console.log('Overlay text:', overlayText);
  }

  // Press Escape
  await page.keyboard.press('Escape');
  await page.waitForTimeout(500);
  
  const overlays2 = await page.$$('div.fixed.inset-0');
  console.log('After Escape, overlays:', overlays2.length);
  
  const textarea2 = await page.$('textarea');
  console.log('After Escape, textarea:', !!textarea2);

  if (!textarea2) {
    // Maybe need to click on the project in sidebar
    const projectLinks = await page.$$('a[href*="project"], button:has-text("TestExp")');
    console.log('Project links:', projectLinks.length);
    
    // Check URL
    console.log('URL:', page.url());
    
    // Take screenshot
    await page.screenshot({ path: '/tmp/debug-modal.png' });
  }
  
  // Clean up
  const projects = await page.evaluate(async () => {
    const res = await fetch('/api/projects');
    return res.json();
  });
  for (const p of projects.filter(p => p.name === 'TestExp')) {
    await page.evaluate(async (pid) => {
      const tr = await fetch('/api/csrf-token');
      const { token } = await tr.json();
      await fetch('/api/projects/' + pid, { method: 'DELETE', headers: { 'X-AIRA-Token': token } });
    }, p.id);
  }
  
  await browser.close();
})();
