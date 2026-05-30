# Literature Search Strategy

## Objective
Build a literature-backed microkinetic modeling framework for heterogeneous catalysis with a Fischer-Tropsch synthesis (FTS) case study.

## Databases and tools
1. ToolUniverse MCP: `Crossref_search_works`, `Crossref_get_work`, `PubMed_search_articles`, `SemanticScholar_search_papers`.
2. Fallback web search for DOI/citation verification when MCP returned tool-name errors, rate limits, or metadata mismatches.

## Queries executed
- "microkinetic modeling heterogeneous catalysis DFT"
- "Fischer-Tropsch synthesis microkinetics"
- "transition state theory surface reactions coverage effects"
- "lateral interactions mean field approximation catalysis"
- "CatMAP microkinetics framework"
- Follow-up exact-title searches for the rate-control, kinetic Monte Carlo, and Fischer-Tropsch mechanism references.

## Inclusion criteria
- Peer-reviewed journal article, review, or book chapter with identifiable DOI.
- Direct relevance to microkinetic modeling, DFT-derived kinetics, rate-control theory, lateral interactions, surface kMC, or Fischer-Tropsch catalysis.
- Sufficient bibliographic metadata to support downstream citation and synthesis.

## Exclusion criteria
- Irrelevant catalysis subfields with no microkinetic content.
- Duplicate records across databases.
- Records with unverified or conflicting bibliographic metadata unless corrected by a second source.

## MCP tool status
- Initial attempts using `SemanticScholar_search` and `PubMed_search` failed because those exact tool names were unavailable in the MCP catalog.
- `SemanticScholar_search_papers` later returned a 429 rate-limit response for one query, so the search was not relied upon as the sole source.
- `Crossref_get_work` failed for some user-supplied DOIs because several targets contained bibliographic mismatches. These were corrected by exact-title searches and web verification.

## Corrected metadata notes
- Stegelmann et al. rate-control paper: correct DOI is `10.1021/ja9000097`.
- Nørskov et al. “Universality in Heterogeneous Catalysis”: correct DOI is `10.1006/jcat.2002.3615`.
- Campbell (2017) “The Degree of Rate Control: A Powerful Tool for Catalysis Research”: correct DOI is `10.1021/acscatal.7b00115`.
- van Santen et al. “Mechanistic Issues in Fischer-Tropsch Catalysis” is indexed as an *Advances in Catalysis* chapter with DOI `10.1016/B978-0-12-387772-7.00003-4`.
- Andersen et al. “A Practical Guide to Surface Kinetic Monte Carlo Simulations” is indexed in *Frontiers in Chemistry* with DOI `10.3389/fchem.2019.00202`.

## Screening summary
- Records identified: 23
- Duplicates removed: 8
- Records screened: 15
- Excluded after title/abstract screening: 2
- Full-text / metadata assessment: 13
- Included studies: 13

## Limitations
This was a single-reviewer screening workflow, so selection bias cannot be excluded. Where a claim depended on a single source, it is marked as such in downstream synthesis.
