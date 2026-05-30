# Search Strategy

## Objective
To identify recent (2020 onward) studies relevant to whole-brain connectome pipelines integrating structural connectivity, functional connectivity, graph-theoretic analysis, disease biomarkers, and reliability assessment.

## Databases and tools
- PubMed via ToolUniverse MCP (`PubMed_search_articles`)
- Crossref via ToolUniverse MCP (`Crossref_search_works`, `Crossref_get_work`)
- Semantic Scholar via ToolUniverse MCP (`SemanticScholar_search_papers`)

## Queries executed
1. "whole brain connectome fMRI dMRI pipeline" (Semantic Scholar; 0 results)
2. "dynamic functional connectivity brain network" (Semantic Scholar; API error 400)
3. "structural functional connectivity graph theory brain" (Semantic Scholar; API error 400)
4. "schizophrenia connectome biomarker resting state fMRI" (PubMed)
5. "Alzheimer disease white matter tractography connectome" (PubMed)
6. "test-retest reliability functional connectivity" (PubMed)
7. "small world network brain modularity hub connectome" (Crossref)

## Inclusion criteria
- Published in 2020 or later
- Human neuroimaging, connectomics, graph/network neuroscience, or connectome biomarker studies
- Relevant to at least one of: structural connectivity, resting-state functional connectivity, dynamic FC, graph measures, disease classification, or reliability
- Peer-reviewed journal article preferred

## Exclusion criteria
- Preprints when peer-reviewed alternatives were available
- Non-brain or non-connectome studies
- Papers lacking neuroimaging methodology relevance to the target pipeline
- Duplicates across searches

## Screening notes
- Single-reviewer screening was used; this is a limitation versus formal two-reviewer PRISMA practice.
- PubMed provided the most relevant disease-specific neuroimaging studies.
- Crossref was used to validate DOI metadata for selected included articles.
- Semantic Scholar MCP partially failed for these queries (400 errors on two searches), so the final synthesis relied mainly on PubMed plus Crossref verification. Because PubMed and Crossref returned sufficient recent peer-reviewed papers, no additional REST fallback was required.

## PRISMA counts used in this project
- Records identified: 18
- Duplicates removed: 4
- Records screened: 14
- Records excluded at title/abstract stage: 3
- Full-text/metadata assessed: 11
- Excluded with reasons: 1 (preprint / lower direct relevance)
- Studies included: 10
