# Literature Search Strategy — ADR Mission Trajectory Design

## Research Topic
Active Debris Removal (ADR) mission optimal trajectory design system:
- Multi-target debris removal trajectory optimization
- Low-thrust orbit transfers
- Rendezvous & proximity operations (Hill/CW equations)
- Tumbling debris attitude estimation
- Capture mechanism dynamics
- Mission sequence optimization

## Databases Searched
| Database | Tool Used | Status | Notes |
|----------|-----------|--------|-------|
| Semantic Scholar | `SemanticScholar_search_papers` | ❌ HTTP 400 error | Year-filtered queries failed |
| OpenAlex | `openalex_literature_search` | ✅ Success | Primary source |
| Crossref | `Crossref_search_works` | ✅ Success | Supplementary source |

## Search Queries

### Query 1: ADR Trajectory Optimization
- **Tool**: openalex_literature_search  
- **Query**: "active debris removal orbital trajectory optimization"  
- **Filters**: year_from=2019, max_results=6  
- **Papers found**: 6 (including Wijayatunga 2023, Medioni 2022, Narayanaswamy 2023, Huang 2022)

### Query 2: Capture Mechanisms / Robotic Systems
- **Tool**: openalex_literature_search  
- **Query**: "tumbling debris attitude estimation rotation capture active debris removal"  
- **Filters**: year_from=2019, max_results=5  
- **Papers found**: 5 (Papadopoulos 2021, Aglietti 2019, Ellery 2019, etc.)

### Query 3: Debris Risk & Sustainability
- **Tool**: openalex_literature_search  
- **Query**: "debris collision risk scoring removal prioritization LEO space sustainability"  
- **Filters**: year_from=2020, max_results=5  
- **Papers found**: 5 (Murtaza 2020, Banken 2022, etc.)

### Query 4: Mission Planning (Crossref)
- **Tool**: Crossref_search_works  
- **Query**: "active debris removal mission sequence optimization delta-v minimization"  
- **Filters**: from-pub-date:2019-01-01, type:journal-article  
- **Papers found**: 6 (Zona 2023, Zhao 2020, Guo 2023, Chen 2024, Simha 2025, Borelli 2026)

## Inclusion Criteria
- Published 2016–2026
- English language
- Subject: orbital mechanics, spacecraft trajectory, ADR, space debris
- Peer-reviewed or IEEE/AIAA conference proceedings

## Exclusion Criteria
- Papers not related to orbital mechanics or debris removal
- Non-English papers
- Papers without DOI

## PRISMA Flow
- Records identified: ~22
- Records screened: 22
- Records excluded (off-topic): 10
- Full-text assessed: 12
- Included in synthesis: 12
