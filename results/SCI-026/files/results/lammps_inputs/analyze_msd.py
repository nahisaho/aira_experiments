#!/usr/bin/env python3
"""
Post-process LAMMPS MSD output to extract Li-ion diffusivity.
D = lim_{t→∞} MSD(t) / 6t   (3D)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

data = np.loadtxt("msd_Li.txt", comments="#")
t_ps  = data[:, 0] * 0.002  # timestep 2 fs → ps
msd   = data[:, 4]           # total MSD (Å²)

# Linear fit to last 60% of trajectory
fit_start = int(len(t_ps) * 0.4)
coeffs = np.polyfit(t_ps[fit_start:], msd[fit_start:], 1)
D_Aps  = coeffs[0] / 6.0                       # Å²/ps
D_cm2s = D_Aps * 1e-16 / 1e-12                 # cm²/s

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(t_ps, msd, "b-", lw=1.5, label="Li MSD")
ax.plot(t_ps[fit_start:], np.polyval(coeffs, t_ps[fit_start:]),
        "r--", lw=2, label=f"Linear fit  D={D_cm2s:.2e} cm²/s")
ax.set_xlabel("Time (ps)")
ax.set_ylabel("MSD (Å²)")
ax.set_title("Li-ion Mean Squared Displacement at Interface")
ax.legend()
plt.tight_layout()
plt.savefig("figures/msd_Li_interface.png", dpi=300, bbox_inches="tight")
print(f"D_Li = {D_cm2s:.3e} cm²/s")
