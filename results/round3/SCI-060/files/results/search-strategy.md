# Search Strategy for VE Methodology Framework

## Research Objective
Identify methodological literature on real-world vaccine effectiveness (VE) estimation, covering: Test-Negative Design (TND), waning immunity, variant-specific VE, healthy-vaccinee bias, booster causal effects, and mRNA hospitalization VE.

## MCP Tool Attempts

| Tool | Status | Notes |
|------|--------|-------|
| `PubMed_search_articles` | ✅ Success | Used for all primary searches |
| `SemanticScholar_search_papers` | ❌ Error 400 | API returned HTTP 400; query likely too long |
| `PMC_search_papers` | Fallback | Not required; PubMed sufficient |

## Search Queries Executed

1. `"test-negative design vaccine effectiveness real-world estimation statistical methods"` → 8 results (PubMed)
2. `"waning COVID-19 vaccine effectiveness Omicron variant mRNA booster"` → 8 results (PubMed)
3. `"healthy vaccinee bias frailty bias observational vaccine study"` → 4 results (PubMed)
4. `"mRNA vaccine effectiveness hospitalization COVID-19 real-world test-negative design 2022 2023"` → 6 results (PubMed)
5. `"variant-specific vaccine effectiveness SARS-CoV-2 Delta Omicron statistical estimation"` → 1 result (PubMed)
6. `"causal inference booster dose vaccine effectiveness propensity score"` → 0 results (PubMed)

## Inclusion Criteria
- Year: 2017–2026
- Language: English
- Study type: Original research, systematic review, or methodological paper
- Topic: VE estimation, TND methodology, waning, bias correction, causal inference for vaccines

## Exclusion Criteria
- Animal studies
- Immunological-only studies (no VE estimate)
- Non-peer-reviewed preprints as primary evidence (labeled ⚠️ if used)

## PRISMA Summary
- Records identified: 27
- Screened: 27
- Included after full-text review: 14 (listed in reference-list.md)
- Excluded: 13 (topic mismatch or duplicate)
