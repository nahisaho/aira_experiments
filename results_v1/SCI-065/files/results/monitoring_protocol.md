# Brain Organoid Bioreactor Monitoring Protocol

## Scope
This protocol integrates online, at-line, and offline biomarker monitoring for maturation control in brain organoid bioreactors.

## Monitoring cadence
- **Online (continuous, 30 min acquisition):** pH, dissolved oxygen (DO), glucose, lactate
- **At-line (daily):** LDH release, multiplex cytokine panel
- **Offline (weekly):** qPCR and immunostaining for OCT4, NANOG, PAX6, SOX2, NESTIN, TBR1, CTIP2, SATB2, MAP2, SYN1, GFAP; weekly electrophysiology benchmarking

## Release criteria by phase
1. **Neural induction (day 0-6):** OCT4 and NANOG falling, pH 7.25-7.40, DO > 55%
2. **Patterning (day 6-25):** PAX6/SOX2/NESTIN peak trajectory, glucose consumption increasing without LDH surge
3. **Cortical differentiation (day 25-50):** TBR1 and CTIP2 increasing, lactate rise matched by stable DO
4. **Maturation (day 50+):** MAP2/SYN1/GFAP increase, low LDH, sustained electrophysiology score

## Process-control actions
- Shewhart rule violation: verify sensor health, then inspect perfusion and aeration hardware within 2 h
- CUSUM trigger without Shewhart breach: check slow drift in feed composition, calibration, and gas blending
- Multivariate anomaly score > 4.5: hold automated feed escalation and perform targeted microscopy/qPCR confirmation

## Escalation logic
- DO low + lactate high -> increase gas transfer and inspect aggregate density
- Glucose low + LDH high -> reduce residence time, refresh medium, and assess necrotic cores
- Cytokine elevation without sensor deviation -> inspect contamination or inflammatory glial overgrowth
