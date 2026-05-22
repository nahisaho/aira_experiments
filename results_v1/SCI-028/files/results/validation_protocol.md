
# Validation Protocol: Tokamak Disruption Prediction AI

## 1. Dataset Requirements

### JET (primary training + test)
- Total shots: ~90,000 (1983–2022)
- Disruption rate: ~15–20% of high-performance shots
- Required signals: Mirnov array (32 coils), Thomson scattering, ECE, magnetics
- Recommended: JET Disruption Database (JDD) — contact EUROfusion
- Shot range for test: 90000–99000 (most recent campaigns)

### KSTAR (cross-device validation)
- Total shots: ~30,000 (2008–2023)
- Disruption rate: ~10% (higher H-mode fraction)
- Required signals: Same as JET (mapped via signal name registry)
- Data access: KSTAR data portal or IMAS MDSplus server

### ASDEX Upgrade (optional supplementary)
- Total shots: ~40,000
- Useful for: impurity-induced disruption scenarios

## 2. Preprocessing Checklist
- [ ] Synchronise signal timestamps to 0.1 ms grid
- [ ] Apply per-device normalisation (IP/IP_max, BT/BT_nom)
- [ ] Remove shots with incomplete data (<80% signal coverage)
- [ ] Label disruption time as last timestamp before dIp/dt < −0.5 MA/ms
- [ ] Assign disruption cause from Mirnov/radiation pattern (automated + manual audit)

## 3. Split Rules
- NEVER mix shot segments from the same discharge across train/test
- Use temporal split (shot ID ascending) to prevent future leakage
- Hold out entire JET campaigns for final test (e.g., C38, C40)
- KSTAR test set: all shots, zero-shot transfer evaluation

## 4. Evaluation Protocol
1. Train on JET train set (shots before 2018)
2. Validate hyperparameters on JET val set (2018–2020)
3. Report final metrics on JET test set (2020+)
4. Report zero-shot metrics on KSTAR (no fine-tuning)
5. Report few-shot metrics on KSTAR (10/50/100 shots fine-tuning)
6. Compare against:
   - APODIS baseline (Versace et al., 2010)
   - SVM baseline (Rattá et al., 2010)
   - LSTM baseline (Kates-Harbeck et al., 2019 — FRNN)

## 5. Minimum Acceptance Thresholds
| Metric              | Minimum | Target |
|---------------------|---------|--------|
| AUC-ROC (JET)       | 0.92    | 0.97   |
| TPR @ FPR=0.05      | 0.90    | 0.95   |
| Avg. Warning Time   | 30 ms   | 100 ms |
| KSTAR AUC (0-shot)  | 0.80    | 0.90   |
| Inference Latency   | <30 ms  | <20 ms |

## 6. Statistical Testing
- Compare models using DeLong test for AUC-ROC differences
- Bootstrap 95% CI (n=1000) for all point estimates
- McNemar's test for shot-level TPR/FPR comparison
- Bonferroni correction for multiple device comparisons (k=3)

## 7. Prospective Validation (future)
- Submit model to JET / KSTAR control room for shadow-mode validation
- Run in parallel with existing system for ≥100 disruptive shots
- Report prospective TPR/FPR without post-hoc threshold adjustment
