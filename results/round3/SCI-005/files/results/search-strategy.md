# Search Strategy

## Objective
Identify recent studies relevant to high-accuracy structural variant detection from Oxford Nanopore and PacBio long-read data, with emphasis on multi-strategy calling, repeat-region handling, complex rearrangements, and hybrid short/long-read validation.

## Databases and tools
1. ToolUniverse discovery for PubMed and Semantic Scholar tools.
2. MCP execution attempts with `PubMed_search_articles`, `SemanticScholar_search_papers`, and `Crossref_search_works`.
3. Fallback to Python `requests` for Crossref lookups after Semantic Scholar MCP returned HTTP 400 and Crossref rate-limited some broad requests.

## Queries
- long-read structural variant detection nanopore pacbio
- PBSV Sniffles SVABA structural variant caller benchmark
- repeat region SV detection telomere centromere
- chromothripsis extrachromosomal DNA detection
- hybrid short long read SV detection

## Inclusion criteria
- Published 2018 or later, preferably 2020+
- Direct relevance to long-read SV detection, benchmarking, repeat/low-complexity analysis, complex SV analysis, or hybrid evidence integration
- DOI available or verifiable
- Method, benchmark, systems, or high-value review papers with concrete findings

## Exclusion criteria
- General sequencing papers without a clear SV focus
- Duplicate preprint/journal pairs when the journal article was identifiable
- Papers with insufficient metadata or unclear relevance

## Screening note
Two-database coverage was achieved through PubMed and Crossref-backed retrieval. Semantic Scholar MCP attempts were documented but failed for the supplied queries. Screening was performed by a single reviewer, which is a PRISMA limitation and should be considered when interpreting study selection.
