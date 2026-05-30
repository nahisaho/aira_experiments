# Search Strategy

This review targeted recent literature (2020 onward) on knowledge-graph-based drug repurposing, biomedical knowledge graph embeddings, translational link prediction, and COVID-19 case studies. The review followed a lightweight PRISMA-style workflow adapted for a single-analyst exploratory review.

## Databases

1. PubMed via ToolUniverse MCP (`PubMed_search_articles`)
2. Crossref via ToolUniverse MCP (`Crossref_search_works`)
3. Semantic Scholar via ToolUniverse MCP (`SemanticScholar_search_papers`) and direct REST fallback

## Queries

- `COVID-19 knowledge graph drug repurposing`
- `biomedical knowledge graph embedding drug discovery`
- `TransE RotatE ComplEx drug discovery knowledge graph`
- `knowledge graph completion drug repurposing biomedical link prediction`

## Inclusion criteria

Studies were included if they: (i) were published in 2020 or later, (ii) addressed drug repurposing, drug discovery, or target discovery using knowledge graphs or graph embeddings, and (iii) provided direct relevance to one of the following themes: COVID-19 repurposing, biomedical KG construction, embedding model evaluation, or explainability/path reasoning.

## Exclusion criteria

Records were excluded if they: (i) were non-biomedical, (ii) lacked a knowledge-graph or graph-reasoning component, (iii) were duplicative preprints when a journal version was available, or (iv) provided only tangential commentary without methodological relevance.

## Screening notes

- Semantic Scholar MCP returned HTTP 400 twice.
- Direct Semantic Scholar REST access returned HTTP 429 without an API key.
- PubMed returned the most targeted results and was used as the source of truth for inclusion.
- Crossref was used mainly for DOI confirmation when needed.

## Quality assessment

Methodological quality was assessed qualitatively using three questions: Does the study define a biomedical KG/task clearly? Does it report a computational evaluation? Does it provide translational relevance for drug or target discovery? Studies meeting at least two of the three criteria were retained.
