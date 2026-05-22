# Urban Traffic Microsimulation + Real-Time Control Optimization
# SUMO/Flow/RLlib Framework

## Installation

```bash
pip install -r requirements.txt
```

### SUMO Installation (required for full simulation)
```bash
sudo add-apt-repository ppa:sumo/stable
sudo apt-get update
sudo apt-get install sumo sumo-tools sumo-doc
export SUMO_HOME="/usr/share/sumo"
```

## Quick Start

```bash
# Run synthetic evaluation (no SUMO required)
python -m src.main_orchestrator

# Run individual component tests
python -m src.models.idm_model
python -m src.agents.marl_signal_control
python -m src.models.demand_estimation
python -m src.models.dynamic_routing
```

## Architecture

```
src/
├── main_orchestrator.py       # Pipeline integration
├── models/
│   ├── idm_model.py           # IDM + MOBIL car-following
│   ├── demand_estimation.py   # Kalman filter OD estimation
│   └── dynamic_routing.py     # A* rerouting + incident detection
├── agents/
│   └── marl_signal_control.py # MAPPO signal controller
├── network/
│   └── sumo_environment.py    # SUMO/Flow interface
└── utils/
configs/
└── simulation_config.yaml     # Full configuration
```
