# README — pyssem-v2

## Overview
This repository implements a **stochastic compartment model** for Low Earth Orbit (LEO) populations — **satellites, derelicts, and debris**.  

Two primary solvers are included:  
- **Euler–Maruyama (SDE, “EM”) solver**: stochastic approximation of the Markov Jump Process.  
- **Discrete-Event Simulation (DES)**: event-driven simulation of collisions, drag, launches, and post-mission disposal.  

In both approaches, the **ODE drift** provides the deterministic backbone and is embedded as the first stage of the solver.  

Theoretical derivations and mid-stage results are summarized in the accompanying files:  
- *LEO_slides_20250826.pdf*  
- *Report_LEO_Draft2.pdf*  

---

## Repository Layout (actual contents)

```
pyssem-v2 3/
├─ Simulations_EM.ipynb                 # Main notebook: runs EM (SDE) simulations
├─ Simultations_DiscreteEvent.ipynb     # Main notebook: runs DES (helper funcs in notebook)
├─ model.py                             # Model facade: builds ScenarioProperties, plotting helpers
├─ example-purdue.json                  # Example scenario config (time horizon, shells, species, etc.)
├─ scenario-properties.pkl              # Cached scenario setup (optional)
├─ scenario-properties-baseline.pkl     # Cached baseline setup (optional)
├─ utils/
│  ├─ simulation/
│  │  ├─ scen_properties.py             # ScenarioProperties: core state, EM integrator, ODE drift
│  │  ├─ species.py                     # Species definitions, symbols/vectors
│  │  ├─ species_pair_class.py          # Collision pair structure
│  │  └─ __init__.py
│  ├─ collisions/
│  │  └─ collisions.py                  # Collision kernels, fragment sizing, pair creation
│  ├─ drag/
│  │  └─ drag.py                        # Atmospheric density models (static/JB2008), drag terms
│  ├─ launch/
│  │  ├─ launch.py                      # Launch models (constant/ADEPT-like), init helpers
│  │  └─ data/
│  │     ├─ launch_36.csv
│  │     ├─ x0_20.csv
│  │     └─ x0_36.csv
│  ├─ pmd/
│  │  └─ pmd.py                         # Post-Mission Disposal terms
│  └─ handlers/
│     └─ handlers.py                    # I/O helpers (downloads, etc.)
└─ .idea/ …                             # IDE metadata (not required to run)
```

---

## Core Modules (by category)

### Main Notebooks
- **Simulations_EM.ipynb** — Run Euler–Maruyama (SDE) simulations.  
- **Simultations_DiscreteEvent.ipynb** — Run Discrete-Event Simulations (DES). Event functions are implemented within the notebook.  

### Core Scripts
- **model.py** — Entry point wrapper. Loads scenario configuration, builds `ScenarioProperties`, and provides plotting utilities.  
- **example-purdue.json** — Example scenario configuration (time horizon, shell setup, parameters).  
- **scenario-properties.pkl / scenario-properties-baseline.pkl** — Saved scenario objects (optional).  

### `utils/simulation/`
- **scen_properties.py** — Defines `ScenarioProperties`.  
  - Core state container for time, shells, and species.  
  - Implements ODE drift, Euler–Maruyama integration, and scenario-level simulation control.  
- **species.py** — Defines parameters for each species (satellite, derelict, debris), including radii, avoidance/disable ratios, drag/PMD flags, and initial states.  
- **species_pair_class.py** — Encodes species-pair interactions (e.g., SS, SD, SN). Stores impact parameter, relative velocity, and collision type metadata.  

### `utils/collisions/`
- **collisions.py** — Collision kernel. Computes interaction rates, splits lethal/disable components, and applies fragment generation formulas for catastrophic vs. non-catastrophic events.  

### `utils/drag/`
- **drag.py** — Atmospheric drag model. Provides simple exponential or JB2008 density models. Calculates outflow (to lower shells) and inflow (from upper shells).  

### `utils/launch/`
- **launch.py** — Launch processes. Defines constant and scenario-based launch rates. Includes helpers for initial conditions (`x0`) and CSV loaders.  
- **data/** — Example initialization and launch CSV files.  

### `utils/pmd/`
- **pmd.py** — Post-Mission Disposal (PMD). Defines removal rates for satellites and derelicts, modeled as delayed exponential processes.  

### `utils/handlers/`
- **handlers.py** — I/O helpers for scenario loading, downloads, and file management.  

---

## How to Use

1. **Environment**  
   Python ≥ 3.9, with standard scientific libraries (`numpy`, `scipy`, `pandas`, `matplotlib`, `sympy`).  

2. **Configure a scenario**  
   Adjust `example-purdue.json` to set duration, shell size, species parameters, and functions (launch, drag, PMD).  

3. **Run simulations**  
   - Open **Simulations_EM.ipynb** to generate stochastic (SDE) sample paths.  
   - Open **Simultations_DiscreteEvent.ipynb** to run event-driven simulations.  

4. **View results**  
   - Plots are produced with `model.py` utilities.  
   - Outputs include species evolution by shell and solver comparisons (ODE vs. SDE vs. DES).  

---

## Notes
- The EM (SDE) solver and DES are the **main focus** of this code.  
- The ODE solver is embedded within both solvers as the deterministic first step.  
- DES helper functions (`run_des`, event rate updates) are currently defined inside the notebook, not as a separate module.  
- Launch rates are parameterized; mapping real-world capacity milestones to time-dependent rates is still under study.  
