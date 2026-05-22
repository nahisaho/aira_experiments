"""Active Debris Removal mission design toolkit."""

from pathlib import Path

GM_EARTH = 3.986e14  # m^3/s^2
R_EARTH = 6371e3  # m
G0 = 9.80665  # m/s^2
SEED = 42

WORKSPACE = Path(__file__).resolve().parent.parent
DATA_DIR = WORKSPACE / "data"
RESULTS_DIR = WORKSPACE / "results"
FIGURES_DIR = WORKSPACE / "figures"
LOGS_DIR = WORKSPACE / "logs"

for _path in (DATA_DIR, RESULTS_DIR, FIGURES_DIR, LOGS_DIR):
    _path.mkdir(parents=True, exist_ok=True)
