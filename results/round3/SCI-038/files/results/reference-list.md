# Reference List — ADR Mission Trajectory Design

## Included Studies (Post-Screening)

1. (Wijayatunga, 2023) Wijayatunga, M., Armellin, R., Holt, H., Pirovano, L., & Lidtke, A. (2023). Design and guidance of a multi-active debris removal mission. *Astrodynamics*, 7(4), 415–432. https://doi.org/10.1007/s42064-023-0159-3

2. (Medioni, 2022) Medioni, L., Gary, Y., Monclin, M., Oosterhof, C., Pierre, G., Semblanet, T., Comte, P., & Nocentini, K. (2022). Trajectory optimization for multi-target Active Debris Removal missions. *Advances in Space Research*. https://doi.org/10.1016/j.asr.2022.12.013

3. (Zona, 2023) Zona, F., Zavoli, A., & Federici, L. (2023). Evolutionary Optimization for Active Debris Removal Mission Planning. *IEEE Access*, 11, 39812–39825. https://doi.org/10.1109/access.2023.3269305

4. (Zhao, 2020) Zhao, Z., Feng, F., & Yuan, J. (2020). A Novel Two-Level Optimization Strategy for Multi-Debris Active Removal Mission in LEO. *Computer Modeling in Engineering & Sciences*, 122(1), 151–175. https://doi.org/10.32604/cmes.2020.07504

5. (Guo, 2023) Guo, Z., Pang, B., & Du, X. (2023). Optimal planning for a multi-debris active removal mission with a partial debris capture strategy. *Chinese Journal of Aeronautics*, 36(10), 308–322. https://doi.org/10.1016/j.cja.2023.03.013

6. (Chen, 2024) Chen, S., Bai, X., & Zhao, Y. (2024). Rapid Sequence Generation for Active Debris Removal Mission Based on Attention Mechanism and Pointer Network. *IEEE Access*, 12. https://doi.org/10.1109/access.2024.3425161

7. (Papadopoulos, 2021) Papadopoulos, E., Aghili, F., Ma, O., & Lampariello, R. (2021). Robotic Manipulation and Capture in Space: A Survey. *Frontiers in Robotics and AI*, 8, 686723. https://doi.org/10.3389/frobt.2021.686723

8. (Aglietti, 2019) Aglietti, G. S., Taylor, B., Fellowes, S., et al. (2019). RemoveDEBRIS: An in-orbit demonstration of technologies for the removal of space debris. *The Aeronautical Journal*, 124(1271), 1–23. https://doi.org/10.1017/aer.2019.136

9. (Murtaza, 2020) Murtaza, A., Pirzada, S. J. H., Xu, T., & Liu, J. (2020). Orbital Debris Threat for Space Sustainability and Way Forward. *IEEE Access*, 8, 44119–44169. https://doi.org/10.1109/access.2020.2979505

10. (Narayanaswamy, 2023) Narayanaswamy, S., & Damaren, C. J. (2023). Equinoctial Lyapunov Control Law for Low-Thrust Rendezvous. *Journal of Guidance, Control, and Dynamics*, 46(5), 998–1006. https://doi.org/10.2514/1.g006662

11. (Simha, 2025) Simha, A., Servadio, S., & Lifson, M. (2025). Optimal Active Debris Removal mission planning to inform policy decisions. *Acta Astronautica*, 227, 1–12. https://doi.org/10.1016/j.actaastro.2024.11.050

12. (Shan, 2016) Shan, M., Guo, J., & Gill, E. (2016). Review and Comparison of Active Space Debris Capturing and Removal Methods. *Progress in Aerospace Sciences*, 80, 18–32. https://doi.org/10.1016/j.paerosci.2015.11.001

---

## Key Themes Synthesized

### Theme 1: Multi-Target Mission Sequence Optimization
(Wijayatunga, 2023), (Medioni, 2022), (Zona, 2023), (Zhao, 2020), (Guo, 2023), (Chen, 2024), (Simha, 2025)
- Primary challenge: combinatorial explosion of sequence ordering
- Methods: evolutionary algorithms, Bayesian optimization, attention networks, two-level methods

### Theme 2: Low-Thrust Trajectory Computation
(Narayanaswamy, 2023), (Wijayatunga, 2023)
- Q-law, Lyapunov guidance laws, Edelbaum's formula
- Electric propulsion (Isp 1000–3000 s)

### Theme 3: Capture Mechanism Technologies
(Papadopoulos, 2021), (Aglietti, 2019), (Shan, 2016)
- Net, harpoon, robotic arm — each with different rotation-rate tolerances
- RemoveDEBRIS mission as benchmark

### Theme 4: Space Debris Risk & Sustainability
(Murtaza, 2020)
- LEO congestion at 800–1000 km altitude
- Cascade effect mitigation requires removing large objects first

## Research Gaps Identified
⚠️ Gap 1: Integrated frameworks combining all sub-problems (scoring, trajectory, capture) in a single pipeline are lacking.  
⚠️ Gap 2: Real-time adaptation to tumbling rate estimation and method selection is underexplored.  
⚠️ Gap 3: Policy-informed optimization (Simha, 2025) is a growing area not yet linked to operational trajectory planners.
