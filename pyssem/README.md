# README — PODAS-pyssem

## Overview
This repository implements a **stochastic compartment model** for Low Earth Orbit (LEO) populations  
**satellites(S), derelicts(N_223kg), and debris(N_0.64kg)**.  

Two primary solvers are included:  
- **Euler–Maruyama (SDE, “EM”) solver**: stochastic approximation of the Markov Jump Process.  
- **Discrete-Event Simulation (DES)**: event-driven simulation of collisions, drag, launches, and post-mission disposal.  

In both approaches, the **ODE drift** provides the deterministic backbone and is embedded as the first stage of the solver.  

Theoretical derivations and mid-stage results are summarized in the accompanying files:  
- *LEO_slides_20250826.pdf*  
- *Report_LEO_Draft2.pdf*  

---

## Repository Layout

```
pyssem/
├─ Simulations_EM.ipynb                 # Main notebook: runs EM (SDE) simulations
├─ Simultations_DiscreteEvent.ipynb     # Main notebook: runs DES (helper funcs in notebook)
├─ model.py                             # Model facade: builds ScenarioProperties, plotting helpers
├─ example-purdue.json                  # Example scenario config (time horizon, shells, species, etc.)
├─ figures                              # Resulting Figures
├─ frames                               # Contains frames for making figures
├─ utils/
│  ├─ simulation/
│  │  ├─ scen_properties.py             # Implements ODE drift, Euler–Maruyama
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
│  ├─ pmd/
│  │  └─ pmd.py                         # Post-Mission Disposal terms
│  └─ handlers/
│     └─ handlers.py                    
└─ .idea/ …                             
```

---

## Core Modules

### Main Notebooks
- **Simulations_EM.ipynb** — Runs Euler–Maruyama (SDE) simulations.  
- **Simultations_DiscreteEvent.ipynb** — Runs Discrete-Event Simulations (DES). 

### Core Scripts
- **model.py** — Loads scenario configuration, builds `ScenarioProperties`, and provides plotting utilities.  
- **example-purdue.json** — Example scenario configuration (time horizon, shell setup, parameters).  

### `utils/simulation/`
- **scen_properties.py** — Defines `ScenarioProperties`.  
  - Core state container for time, shells, and species.  
  - Implements ODE drift, Euler–Maruyama integration, and scenario-level simulation control.  
- **species.py** — Defines parameters for each species, initial states.  
- **species_pair_class.py** — Encodes species-pair interactions (e.g., SS, SD, SN). Stores impact parameter, relative velocity, and collision type metadata.  

### `utils/collisions/`
- **collisions.py** — Collision kernel. Computes interaction rates, splits lethal/disable components, and applies fragment generation formulas for catastrophic vs. non-catastrophic events.  

### `utils/drag/`
- **drag.py** — Atmospheric drag model. Calculates outflow (to lower shells) and inflow (from upper shells).  

### `utils/launch/`
- **launch.py** — Launch processes. Defines constant and scenario-based **launch rates** and **initial populations**
### `utils/pmd/`
- **pmd.py** — Post-Mission Disposal (PMD). Defines removal rates for satellites and derelicts, modeled as delayed exponential processes.  

---

## How to Use

1. **Configure a scenario**  
   Adjust `example-purdue.json` to set duration, shell size, species parameters, and functions (launch, drag, PMD).  

2. **Run simulations**  
   - Open **Simulations_EM.ipynb** to generate stochastic (SDE) sample paths.  
   - Open **Simultations_DiscreteEvent.ipynb** to run event-driven simulations.  

3. **View results**  
   - Plots are stored in `figures/`

---

## Notes
- The EM (SDE) solver and DES are the main focus of this code.  
- The ODE solver is embedded within both solvers as the deterministic first step.  
- Launch rates are parameterized; time-dependent rates is still under work.  
