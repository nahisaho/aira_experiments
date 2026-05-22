# Preprocessing Log

- Seed fixed to 42 for `random` and `numpy`.
- Simulated 100 lots across 5 supply-chain stages.
- Generated 500 stage-level transaction records.
- Serialized certifications and alerts for CSV export.
- Batched transactions into 20 mined blocks plus one genesis block.
- Verified blockchain link integrity and Merkle inclusion before saving outputs.
