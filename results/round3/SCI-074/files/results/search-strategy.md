# Search Strategy

## Research Topic
GPS-denied autonomous UAV navigation using Visual-SLAM + obstacle avoidance for indoor warehouse inspection

## Databases Queried
1. Crossref Works API (primary, successful)
2. Semantic Scholar API (rate-limited 429; 3 queries attempted, 0 successful)
3. PubMed via ToolUniverse MCP (returned empty for robotics topic)

## MCP Tool Status
| Tool | Status | Notes |
|------|--------|-------|
| `SemanticScholar_search_papers` | ❌ Rate-limited (HTTP 429) | Attempted 2× in parallel, then waited; all failed |
| `Crossref_search_works` | ✅ Success | Primary source; 8 queries executed |
| `PubMed_search_articles` | ⚠️ Partial | Returned 0 results (robotics not indexed in PubMed) |

## Search Queries

### Crossref
1. `visual inertial odometry SLAM UAV` (filter: from-pub-date:2020)
2. `VSLAM autonomous drone indoor navigation` (filter: from-pub-date:2020)
3. `OctoMap 3D mapping ROS robot` (filter: from-pub-date:2020)
4. `EGO-Planner trajectory UAV` (filter: from-pub-date:2020)
5. `ORB-SLAM visual SLAM` (filter: from-pub-date:2020)
6. `VINS-Mono visual inertial state estimation` (filter: from-pub-date:2019)
7. `ORB-SLAM3 visual inertial SLAM` (filter: from-pub-date:2019)
8. `dynamic object detection tracking UAV flight`
9. `ROS2 PX4 autonomous quadrotor`
10. `warehouse inventory drone autonomous inspection`
11. `FASTER agile trajectory planning UAV Tordesillas`

## Inclusion Criteria
- Published 2018–2026
- Topic: visual SLAM / VIO / obstacle avoidance / UAV path planning / warehouse drone
- Peer-reviewed (conference or journal)

## Exclusion Criteria
- Medical/biological papers
- Not directly relevant to UAV navigation

## Results
- Retrieved: ~50 candidate records
- Included after screening: 12 (see reference-list.md)
- Excluded: off-topic, duplicate DOIs, or insufficient metadata
