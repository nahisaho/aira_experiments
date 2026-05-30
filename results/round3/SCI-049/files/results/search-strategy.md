# Search Strategy

## Objective
Identify recent, real, DOI-verifiable literature for automated quality control and anomaly detection in large-scale scientific data streams, emphasizing changepoint detection, concept drift, unsupervised anomaly detection, explainability, and CERN/LIGO-style monitoring.

## Databases and tools
- ToolUniverse MCP discovery: `tooluniverse-find_tools`
- ToolUniverse MCP literature APIs: `Crossref_search_works`, `Crossref_get_work`, `SemanticScholar_search_papers`, `SemanticScholar_get_paper`, `PubMed_search_articles`, `InspireHEP_search_papers`
- Web fallback: `web_search`

## Core queries
- "change point detection time series"
- "anomaly detection scientific data"
- "concept drift detection streaming"
- "isolation forest anomaly detection"
- "PELT changepoint detection"
- "Deep SVDD anomaly detection"
- "anomaly detection CERN detector monitoring"
- "explainable anomaly detection SHAP"

## Inclusion criteria
1. Direct relevance to changepoint detection, anomaly detection, concept drift, explainability, or scientific data quality monitoring.
2. Real publication metadata with a verifiable DOI.
3. Preference for 2020 or later publications; foundational pre-2020 methods retained when they are methodologically essential.
4. Sufficient detail to support method selection and discussion.

## Exclusion criteria
1. Irrelevant search hits caused by noisy keyword matching.
2. Peer-review reports, decision letters, or metadata records without the target article.
3. Results without trustworthy citation metadata or without a DOI.
4. Domain-mismatched results that did not contribute to the pipeline design.

## MCP outcome
- `tooluniverse-find_tools`: succeeded.
- `tooluniverse-execute_tool` with Crossref, PubMed, and InspireHEP: mixed success; Crossref and InspireHEP returned useful records, PubMed and some broad Semantic Scholar queries returned sparse results.
- `SemanticScholar_get_paper`: rate-limited (HTTP 429) for at least one DOI lookup.
- `Crossref_get_work` failed for one arXiv DOI lookup (`10.48550/arXiv.2501.13789`), so DOI verification for that item was not used in the final source-of-truth list.
- `web_search`: succeeded and supplied fallback metadata for 2020+ papers, especially explainability and concept-drift references.

## Screening note
Single-reviewer screening was performed in this run; this is a PRISMA limitation and may increase selection bias.
