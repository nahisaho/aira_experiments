# Literature Survey: Surface Code Simulation Framework

## MCP Tool Search Record
- **Tools attempted**: SemanticScholar_search_papers (error 400), ArXiv_search_papers (success)
- **SemanticScholar status**: API returned 400 errors for queries with multiple parameters; fallback to ArXiv
- **ArXiv status**: Successfully retrieved papers

## Key Papers Found

### 1. Stim: a fast stabilizer circuit simulator (Gidney, 2021)
- **Authors**: Craig Gidney
- **Year**: 2021
- **ArXiv**: 2103.02202
- **Key findings**: Fast tableau-based stabilizer simulator using SIMD; linear-time reference samples + Pauli frame propagation; can analyze distance-100 surface code (20k qubits) in 15s
- **Relevance**: Core simulator used in this work

### 2. Sparse Blossom: correcting a million errors per core second (Higgott & Gidney, 2023)
- **Authors**: Oscar Higgott, Craig Gidney
- **Year**: 2023
- **ArXiv**: 2303.15933
- **Key findings**: Sparse blossom MWPM algorithm; avoids all-to-all Dijkstra; processes d=17 circuit at <1μs/round at 0.1% noise; open-sourced in PyMatching 2
- **Relevance**: MWPM decoder used in this work

### 3. Performance enhancement of surface codes via recursive MWPM decoding (deMarti iOlius et al., 2022)
- **Authors**: Antonio deMarti iOlius, Josu Etxezarreta Martinez, Patricio Fuentes, Pedro M. Crespo
- **Year**: 2022
- **ArXiv**: 2212.11632
- **Key findings**: Recursive MWPM decoder improves threshold by 18% under depolarizing noise; 105.5% improvement under i.ni.d. noise; considers correlated X/Y/Z errors
- **Relevance**: Advanced MWPM optimization baseline

### 4. Fault-Tolerant Weighted Union-Find Decoding on the Toric Code (Huang et al., 2020)
- **Authors**: Shilin Huang, Michael Newman, Kenneth R. Brown
- **Year**: 2020
- **ArXiv**: 2004.04693
- **Key findings**: Weighted UF increases threshold from 0.38% to 0.62% (toric code, circuit-level); near-linear complexity; benchmarked vs matching decoder
- **Relevance**: Union-Find decoder baseline

### 5. Union-find quantum decoding without union-find (Griffiths & Browne, 2023)
- **Authors**: Sam J. Griffiths, Dan E. Browne
- **Year**: 2023
- **ArXiv**: 2306.09767
- **Key findings**: Analysis showing UF decoder underutilizes disjoint-set structure; linear worst-case complexity; improved architectural designs
- **Relevance**: Union-Find decoder analysis

### 6. Decoder Performance in Hybrid CV-Discrete Surface-Code Threshold Estimation (Wayo et al., 2026)
- **Authors**: Dennis Delali Kwesi Wayo et al.
- **Year**: 2026
- **ArXiv**: 2603.06730
- **Key findings**: At d=5 Pauli baseline, MWPM reduces mean LER from 0.384 to 0.260 vs UF; threshold crossing median ~0.053
- **Relevance**: MWPM vs UF comparison baseline

### 7. Surface Code with Imperfect Erasure Checks (Chang et al., 2024)
- **Authors**: Kathleen Chang et al.
- **Year**: 2024
- **ArXiv**: 2408.00842
- **Key findings**: Non-Pauli leakage noise; threshold still >2x that of Pauli noise with imperfect checks; effective error distance degradation
- **Relevance**: Non-Pauli/leakage noise effects

### 8. Spatially parallel decoding for multi-qubit lattice surgery (Lin et al., 2024)
- **Authors**: Sophia Fuhui Lin et al.
- **Year**: 2024
- **ArXiv**: 2403.01353
- **Key findings**: Spatially parallel decoding for lattice surgery; compatible with hardware accelerators; maintains fault-tolerance during merges
- **Relevance**: Lattice surgery simulation

### 9. Minimum-Weight Parity Factor Decoder for Quantum Error Correction (Wu et al., 2025)
- **Authors**: Yue Wu et al.
- **Year**: 2025
- **ArXiv**: 2508.04969
- **Key findings**: HyperBlossom/Hyperion: 4.8x lower LER vs MWPM at d=11; almost-linear average runtime; unifies graph-based decoders
- **Relevance**: Advanced decoder comparison

### 10. Lattice surgery with Bell measurements (Haug et al., 2025)
- **Authors**: Trond Hjerpekjøn Haug et al.
- **Year**: 2025
- **ArXiv**: 2510.13541
- **Key findings**: Bell-measurement lattice surgery; 40% entanglement resource saving; stronger logical error suppression for given entanglement rate
- **Relevance**: Lattice surgery protocol analysis

### 11. Union-Intersection Union-Find for Decoding Depolarizing Errors (Lin & Lai, 2025)
- **Authors**: Tzu-Hao Lin, Ching-Yi Lai
- **Year**: 2025
- **ArXiv**: 2506.14745
- **Key findings**: UIUF algorithm reduces LER by >1 order of magnitude at ~1e-5; outperforms MWPM on rotated surface codes; linear-time complexity
- **Relevance**: Advanced UF decoder variant

### 12. High threshold error correction for the surface code (Wootton & Loss, 2012)
- **Authors**: James R. Wootton, Daniel Loss
- **Year**: 2012
- **ArXiv**: 1202.4316
- **Key findings**: Algorithm correcting depolarizing noise up to 18.5% threshold; polynomial time complexity
- **Relevance**: Foundational surface code decoding result
