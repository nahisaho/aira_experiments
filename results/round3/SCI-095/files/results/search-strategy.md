# Search Strategy Documentation

## Research Topic
Quantitative analysis of the impact of Open Access (OA) and Open Data on the research community, covering: OA citation advantage, data sharing patterns, preprint server roles, FAIR compliance, citizen science, and life sciences case studies.

## Databases Searched
| Database | Tool Used | Status |
|----------|-----------|--------|
| Semantic Scholar | `SemanticScholar_search_papers` | ⚠️ Rate-limited (HTTP 429 / 400) |
| PubMed | `PubMed_search_articles` | ✅ Success |
| Crossref | `Crossref_search_works` | ✅ Success |

**MCP Connection Log:**
- `SemanticScholar_search_papers` attempted 4 times with queries: "open access citation advantage bibliometrics causal inference", "FAIR principles open data research impact altmetrics", "preprint bioRxiv peer review", "open access citation advantage". All returned HTTP 429 (rate-limit) or 400 errors. No Semantic Scholar papers retrieved via MCP.
- Fallback: PubMed and Crossref MCP tools succeeded; 20+ candidate papers retrieved.
- Additional citation data recovered from Crossref JSON output files.

## Search Queries

### Cluster 1: OA Citation Advantage (OACA)
- "open access citation advantage" (Crossref, 2018–)
- "open access journal article citation count comparison meta-analysis" (Crossref, 2020–)
- "open access citation advantage OA effect quasi-experiment difference-in-differences" (Crossref, 2020–)

### Cluster 2: FAIR Data Principles & Open Data
- "FAIR data principles open science research reproducibility" (PubMed)
- "Wilkinson FAIR data principles guiding findable accessible interoperable reusable" (PubMed)
- "open access data sharing research impact" (Crossref, 2019–)

### Cluster 3: Preprints & Peer Review
- "preprint server bioRxiv medRxiv peer review scientific publication" (PubMed)
- "open access scientific publication impact altmetrics social media" (PubMed)

### Cluster 4: Citizen Science
- "citizen science volunteer research contribution impact ecology biodiversity" (PubMed)

## Inclusion Criteria
- Peer-reviewed journal articles or high-quality reviews/systematic reviews
- Published 2016–2026 (with emphasis on 2020+)
- English language
- Directly relevant to OA, open data, FAIR principles, preprints, or citizen science impact

## Exclusion Criteria
- Non-English publications
- Conference abstracts only
- Tangentially related health/clinical studies unrelated to OA/open data
- Papers without DOIs
