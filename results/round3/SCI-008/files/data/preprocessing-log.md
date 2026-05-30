# Preprocessing Log

- No private or patient-level data were used.
- Literature metadata were normalized into a common table with citation key, year, DOI, and thematic tags.
- For the synthetic biomedical knowledge graph, entity identifiers will be generated programmatically and linked through typed relations.
- Random seeds fixed at 42 for `random`, `numpy`, and `torch` (when available).
- COVID-19 seed entities will be manually curated to ensure that candidate-drug ranking remains interpretable while still containing realistic label noise.
