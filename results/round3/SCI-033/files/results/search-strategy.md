# Search Strategy

## Objective
Identify literature relevant to quantum machine learning expressibility, kernel methods, data encoding, barren plateaus, and noisy execution for a benchmark-oriented comparison with classical baselines.

## Databases and tools
- ToolUniverse discovery: `tooluniverse-find_tools` for Semantic Scholar, PubMed, and Crossref capability discovery.
- ToolUniverse execution: `SemanticScholar_search_papers`, `SemanticScholar_get_paper`, and `Crossref_search_works`.
- Web fallback: `web_search` for citation disambiguation when Crossref rate limiting or Semantic Scholar query failures occurred.

## Search strings
1. `quantum machine learning expressibility entangling capability parameterized quantum circuits quantum kernels barren plateaus`
2. `Quantum Machine Learning in Feature Hilbert Spaces Schuld Killoran 2019`
3. `A rigorous and robust quantum speed-up in supervised machine learning Liu 2021`
4. `Variational Quantum Algorithms Cerezo 2021 Nature Reviews Physics`
5. `Supervised learning with quantum-enhanced feature spaces Havlicek 2019 Nature`
6. `Power of data in quantum machine learning Huang 2021`
7. `The power of quantum neural networks Abbas 2021`
8. `Exponential concentration in quantum kernel methods Thanasilp`

## Inclusion criteria
- Peer-reviewed journal articles or clearly labeled preprints directly relevant to QML benchmarking.
- Papers addressing at least one of: expressibility, entanglement, kernel methods, encoding, barren plateaus, or NISQ/noise constraints.
- English-language sources published from 2018 onward, plus foundational context where necessary.

## Exclusion criteria
- Articles focused only on chemistry or optimization with no machine-learning or benchmarking relevance.
- Duplicated metadata entries.
- Non-specific overview pages lacking bibliographic traceability.

## Honest MCP outcome log
- `tooluniverse-find_tools` succeeded and listed Semantic Scholar, PubMed, and Crossref tools.
- Broad `SemanticScholar_search_papers` query returned zero items for the combined benchmark query.
- `SemanticScholar_get_paper` succeeded for Sim et al. (2019) and confirmed DOI metadata.
- `Crossref_search_works` resolved most target titles, including Liu et al. (2021), Huang et al. (2021), and Preskill (2018).
- One Crossref request returned HTTP 429 during title disambiguation, so `web_search` was used to verify Abbas et al. (2021), Cerezo et al. (2021), and the kernel-method reference associated with Schuld (2021).
- Because the exact title requested as `Quantum models as kernel machines` did not resolve cleanly through MCP, the verified closely related source `Supervised quantum machine learning models are kernel methods` was included and explicitly labeled.

## Screening limitation
This was a targeted, single-reviewer screening exercise rather than a full PRISMA two-reviewer review. The limitation is recorded in the report and paper.
