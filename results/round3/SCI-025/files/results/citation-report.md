# Citation Report

## Verification summary
- Survey references checked against `results/reference-list.md`: 12/12 accounted for in `paper.md`.
- In-text citation mapping: 12/12 references are cited in the manuscript body.
- DOI verification via Crossref API: 12/12 DOIs returned valid records.
- Unverified reference rate: 0%.
- Bulk numeric citation patterns: none detected.

## Semantic alignment notes
- Introduction and Related Work correctly use Min et al. (2020), Fransen et al. (2023), and Baldera-Moreno et al. (2022) to frame the literature gap.
- Methods cites Malashin et al. (2025), Knott et al. (2020), Tournier et al. (2023), Orlando et al. (2023), and Penas et al. (2024) in method-specific contexts rather than bulk unsupported clusters.
- Results remain primarily data-driven, which is appropriate because claims there are supported by generated outputs rather than literature alone.

## Fixes applied
- Corrected a misspelling in `report.md`: `Fransenn et al.` -> `Fransen et al.`

## Remaining risks
- `report.md` uses mixed Japanese and English citation punctuation in a few places (e.g., full-width parentheses). This is stylistic rather than a metadata error.
- The manuscript is internally consistent, but experimental validation is still absent because the framework is based on synthetic and simulation-derived data.
