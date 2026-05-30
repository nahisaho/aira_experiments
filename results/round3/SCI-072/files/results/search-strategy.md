# Search Strategy

## Objective
Identify efficient algorithms and benchmark studies for Multi-Agent Path Finding (MAPF), with emphasis on optimal, bounded-suboptimal, and fast scalable methods relevant to warehouse robotics.

## Databases and tools
1. ToolUniverse MCP attempt: `SemanticScholar_search_papers` (rate-limited / 400 errors on several MAPF queries)
2. ToolUniverse MCP fallback: `PubMed_search_articles` (very limited relevance for robotics/planning MAPF queries)
3. ToolUniverse MCP success: `Crossref_search_works`
4. Python `urllib` fallback to Crossref REST API for title-based retrieval
5. Python `urllib` fallback to Semantic Scholar Graph API for title-based retrieval (received HTTP 429 rate limits during this run)

## Search strings
- "multi-agent path finding CBS conflict based search"
- "EECBS bounded-suboptimal multi-agent path finding"
- "LaCAM quick multi-agent pathfinding"
- "lifelong multi-agent path finding warehouse"
- "multi-agent pathfinding continuous time"
- "priority inheritance with backtracking multi-agent path finding"
- "anytime multi-agent path finding"

## Inclusion criteria
- Directly addresses MAPF algorithms, variants, benchmarks, or warehouse/lifelong MAPF.
- Published in peer-reviewed venues or archival conference proceedings.
- Provides algorithmic insight, benchmark results, or domain-specific MAPF extensions.
- English language.
- Priority given to 2020+ papers, while retaining seminal pre-2020 anchors.

## Exclusion criteria
- Generic motion-planning papers without MAPF formulation.
- Non-archival abstracts or tutorial-only records.
- Duplicates of the same work across extended abstracts and journal versions unless needed for provenance.

## Date range
Primary emphasis: 2020-2025. Anchor papers before 2020 were retained for foundational context.

## Screening notes
This review was completed by a single reviewer, which is a PRISMA limitation. To mitigate that limitation, key claims were cross-validated across foundational MAPF surveys, algorithm papers, and application-oriented warehouse studies.
