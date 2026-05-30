Now let me build the full framework in parallel while doing more targeted literature searches:The workflow ran successfully. Now let me create the test file and then write the papers:Let me fix that test — at Γ the gap depends only on M, not A1/A2. Need to test at finite k:All 14 tests pass. Now let me load the results and write the reports:report.md needs expansion to 850+ words. Let me expand it:Still need a bit more in report.md:Both files pass the word-count gates (report: 949, paper: 3568). Now let's do final cleanup:All files are in order. Here is the complete summary:

---

## Experiment Complete: Topological Insulator Design Framework

### Key Scientific Findings

1. **SOC-induced topological gap**: Band gap opens from 0.048 eV (without SOC) to **0.419 eV** with SOC (8.8× enhancement) in the Bi₂Se₃ k·p model; negative Dirac mass M = −0.28 eV confirms Z₂ = (1;1,1,0) Strong TI classification.

2. **Surface state Dirac cone**: Slab calculation (18 layers) resolves topologically protected surface states with **v_D = 11.4 eV·Å** (≈ 1.7 × 10⁵ m/s) and Dirac point at E_D = 0.348 eV — consistent with experimental ARPES order-of-magnitude.

3. **Phase diagram**: SOC scaling maps the topological phase boundary; the Bi₂Se₃ parameterisation lies firmly in the topological phase at all λ > 0, with gap monotonically increasing with SOC strength.

4. **Material screening**: Of 20 Bi₂Se₃-class candidates, **16 pass all 6 TI criteria** (score 6/6). Top targets: **TlBiTe₂** and **SnBi₂Te₄** (gap = 0.25 eV), and **BiSbTeSe₂** (gap = 0.26 eV) — optimal for room-temperature applications.

5. **Z₂ invariants**: Fu-Kane parity criterion correctly classifies Z₂ = (1;1,1,0); Wilson-loop 3D calculation requires all 6 TRI planes (current 2D implementation gives Z₂ = 0 on single kz=0 plane — documented limitation).

### Most Important Figure
![Surface State Dirac Cone](figures/fig2_surface_states.png)

### Deliverables
- **Source code**: 6 modules (1,595 lines total)
- **Tests**: 14 test cases, 14 passed
- **Report**: report.md (949 words)
- **Paper**: paper.md (3,568 words, full IMRaD structure)
- **Figures**: 5 (band structure, surface states, phase diagram, Wilson loop, screening)
- **Results**: screening_results.csv (20 materials), results_summary.json

### Limitations and Future Work
- Wilson-loop Z₂ requires all 6 TRI planes with ≥50×50 k-mesh for 3D accuracy
- k·p model overestimates gap (~30%); GW many-body corrections needed for quantitative agreement
- Magnetic TIs (MnBi₂Te₄) require time-reversal-breaking invariant methods (Chern number / axion angle)
- Full pipeline integration with Quantum ESPRESSO + Wannier90 + Z2Pack required for ab-initio validation

### MCP Tool Notes
ArXiv API timed out twice (network error); literature search succeeded via OpenAlex, CORE, and Crossref. 8 key references (2020+) identified and cited with DOIs.