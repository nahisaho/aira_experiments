Now let me implement the full simulation framework:Now let me create the main experiment runner:Let me fix the noise model implementation to produce distinct results:Now let me check the figures and then write the report:Now let me save the literature results and write the report and paper:Now write the academic paper:Now add basic tests and verify everything:All tests pass. Now verify word counts and clean up:Report needs expanding - let me add more detail:958 words now. Final cleanup:Both files meet requirements. Let me append the final log entry:---

## Experiment Complete: Surface Code Simulation Framework

### Key Scientific Findings

1. **Circuit-level threshold** — $p_{th} = 0.40\% \pm 0.13\%$ (R² = 0.944) for rotated surface code with depolarizing noise and MWPM decoding, consistent with literature range of 0.3–0.7%.

2. **Noise model sensitivity** — At $p=1\%$, $d=5$: phase damping yields **14× lower LER** (0.57%) and amplitude damping **5.5× lower LER** (1.47%) vs. symmetric depolarizing (8.13%), highlighting the importance of hardware-matched decoder design.

3. **Leakage is catastrophic** — Leakage noise with $f_{\rm leak}=0.2$ increases LER by **5.8×** (to 47.1% at $p=1\%$), underscoring the necessity of leakage reduction units in physical devices.

4. **Error suppression onset** — At $p=0.004 \approx p_{th}$, LER decreases from 0.99% ($d=3$) to 0.43% ($d=7$), indicating the onset of distance-dependent exponential error suppression.

5. **Framework speed** — All 4 experiments completed in **under 5 seconds** on single-core hardware using Stim 1.16.0 + PyMatching 2.4.0.

### Most Important Figure
![Noise Model Comparison](figures/noise_comparison.png)

### Deliverables
| Type | Count | Details |
|------|-------|---------|
| Source modules | 5 | noise_models, surface_code, decoders, simulation, visualization |
| Main script | 1 | main.py (~200 lines) |
| Tests | 5 | All passing |
| Figures | 7 | PNG, 150 DPI |
| Results | 10 | JSON data files |
| Papers | 2 | report.md (958 words), paper.md (3,464 words) |

### Literature Survey (ArXiv via ToolUniverse MCP)
12 papers identified; SemanticScholar API returned 400 errors (fallback to ArXiv successful). Key sources: Gidney 2021 (Stim), Higgott & Gidney 2023 (PyMatching 2), Huang et al. 2020 (Union-Find), Chang et al. 2024 (leakage), Lin & Lai 2025 (UIUF).

### Limitations
- Small code distances ($d \leq 7$) limit threshold estimation accuracy
- UF decoder approximated via PyMatching (genuine UF not available)
- Lattice surgery uses single-patch proxy, not full 2-patch merge