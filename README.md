# README — RPO Phase Simulation (Asteroid Capture Mission)

This repository simulates the **Rendezvous & Proximity Operations (RPO)** portion of an asteroid capture mission. The main entry point is `RPO_Main_code.py`, which successively runs the different steps of the RPO campaign and outputs per‑step **time of flight (TOF)**, **ΔV**, and **propellant mass**.

---

## 1) Context & Scope
This code is one piece of a larger **integrated modeling & simulation framework** for the asteroid capture missions development described in the attached report (Multidisciplinary Framework for the Conceptual Analysis of
Asteroid Capture Missions). That framework couples:

- **Low‑thrust trajectory design** (outbound to the asteroid and inbound to cislunar space),
- **Rendezvous & Proximity Operations (RPO)** including approach, spin‑axis alignment, spin‑rate matching, capture, and post‑capture detumbling, and
- **Spacecraft subsystem sizing**, including capture system mass estimation.

Within that end‑to‑end context, this repository focuses on the **RPO phase only**, producing the RPO timeline and propellant usage that the broader framework consumes.

---

## 2) Problem Description
During RPO, a chaser spacecraft must: (1) approach and position itself along the asteroid’s spin axis; (2) **match the asteroid’s spin** so that the capture system is properly oriented; (3) perform a close approach and capture; (4) **detumble** the combined spacecraft‑asteroid composite to a safe state. Each action has guidance/controls logic and a corresponding propellant cost. The goal is to estimate **TOF and propellant** per sub‑phase and the **capture system (CS) mass** as a function of asteroid size.

---

## 3) Methods:
- **Step 1 — Observation & Approach**: here is no dedicated code module for this phase. Its propellant consumption is assumed negligible, and its time of flight is fixed based on literature values. These activities are modeled only by constants in the main code, as the detailed dynamics and operations are beyond the present scope.
- **Step 2 — Spin‑axis rendezvous**: A low‑thrust transfer inside the asteroid‑fixed RTN frame places the chaser at a target point on the asteroid’s spin axis (e.g., ~10 m above surface). The model computes thrust history, ΔV, TOF, and propellant, with thrust capped by available power & Isp.
- **Step 3 — Spin‑rate matching**: A sliding‑mode attitude controller commands body torques to track the asteroid spin while keeping the capture system aligned. Torques are mapped to RCS thruster firings via a configuration pseudo‑inverse; propellant is integrated.
- **Step 4 — Close approach**: A short final approach with a fixed ΔV assumption (literature‑based) is converted to propellant via the rocket equation.
- **Step 5 — Capture system sizing**: A parametric **bag capture system** estimate (inflatable exoskeleton + airbags + robotic arms + actuators) scales with asteroid diameter and basic material properties.
- **Step 6 — Detumbling (post‑capture)**: Another sliding‑mode controller acts on the **composite** inertia (spacecraft + asteroid) to bring angular rates near zero while bounding max thruster levels; propellant is integrated from commanded thrust.

---

## 4) Repository Map

- **`RPO_Main_code.py`** — Orchestrates RPO steps 1–6; returns final mass, total TOF, and per‑step arrays for ΔV, TOF, propellant, and capture system mass.
- **`step_2.py`** — Low‑thrust guidance from initial RTN state to a target point on the spin axis; enforces thrust available from power & Isp; outputs TOF, ΔV, propellant.
- **`Step3_code.py`** — Spin‑rate matching & pointing control (regulation + tracking phases). Sliding‑mode control; RCS configuration and pseudo‑inverse mapping; outputs propellant, TOF, ΔV; optional plots.
- **`step5_CS_mass.py`** (and **`BagCS_mass.py`**) — Bag capture system mass model vs. asteroid diameter & uncertainty; material/geometry scalings and actuator assumptions.
- **`Step6_Python.py`** — Detumbling of the **captured composite**, including parallel‑axis inertia combination, RCS mapping, and propellant/TOF estimation; optional plots.
- **`step_2_README.txt`** — In‑code reference & I/O description for Step 2.
- **`ACM_SciTech_Draft.pdf`** — Report providing the overall mission framework and RPO context.

---

## 5) Inputs & Outputs

### Common physical/mission inputs
- Asteroid: mean diameter, mass, assumed spin axis & spin rate
- Spacecraft: initial mass at RPO start, Isp, available power, efficiency
- Geometry & materials for capture system (internal in Step 5)

### RPO_Main outputs
- **`mass_sc_f`** — Spacecraft mass after RPO
- **`tof_total`** — Total RPO time of flight
- **`prop_consumption_vals`** — Array of propellant mass per step (1…6)
- **`tof_vals`** — Array of TOF per step
- **`deltaV_vals`** — Array of ΔV per step
- **`mass_CS`** — Capture system mass estimate

Each submodule may also show optional plots (attitude, torques/commands, quaternions, trajectory).

---

## 6) Quick Start
1. **Install**: Python 3.10+, `numpy`, `matplotlib`, (and `scipy` for Step 3).
2. **Open** `RPO_Main_code.py` and scroll to the **“EXAMPLE FOR TESTING”** block. Adjust:
   - Asteroid: `mass_ast`, `diameter_ast`, `spin_axis_ast`, `spin_rate_ast`
   - Spacecraft / propulsion: `mass_sc_i`, `Isp`, `power`, `efficiency`
   - RTN initial & target states: `x0` and `x_t`
3. **Run** the script. It will print the final mass, per‑step propellant/TOF/ΔV, and CS mass.

> Note: Set plotting flags in `step_2.py`, `Step3_code.py`, and `Step6_Python.py` to visualize transfers and attitude behavior.

---

## 7) Step‑by‑Step Explanation (what each file does)

### Step 2 — Maneuver to the spin axis (`step_2.py`)
- Computes coefficients for RTN‑frame guidance to the spin‑axis target.
- Limits thrust by power and Isp (SEP‑style); integrates thrust usage into ΔV and propellant.
- Returns **total thrust on‑time (N·s)**, **TOF (s)**, **propellant (kg)**, **ΔV (m/s)**.

### Step 3 — Match spin rate & maintain pointing (`Step3_code.py`)
- Two phases: **regulation** (align CS with asteroid spin axis) then **tracking** (spin‑rate following).
- Sliding‑mode gains tune convergence; error quaternion thresholds define success.
- **RCS mapping** converts torque commands to thruster firings using a reduced configuration matrix (pseudo‑inverse); optional PWPF shaping.
- Returns **propellant**, **TOF**, **ΔV**; optional plots of quaternions, torque commands, body‑frame target axis.

### Step 4 — Close approach
- Uses a fixed ΔV and short TOF assumption to compute propellant with the rocket equation.

### Step 5 — Capture system mass (`step5_CS_mass.py` / `BagCS_mass.py`)
- Parametric **bag capture system** model:
  - Inflatable exoskeleton & airbags (Kevlar‑like areal density),
  - Robotic arms (carbon fiber volume scaling),
  - Arm actuators (mass per unit).
- Geometric scalings derive surfaces/volumes from asteroid diameter and an uncertainty term; returns **total CS mass**.

### Step 6 — Detumble the composite (`Step6_Python.py`)
- Builds composite inertia (spacecraft + asteroid) via parallel‑axis theorem and estimated COM shift.
- Sliding‑mode regulation drives rates to zero; **RCS mapping** enforces max thrust constraints through iteration on gains.
- Integrates commanded thrust to propellant; returns **propellant**, **TOF**, **ΔV**; optional diagnostic plots.

---

## 8) Assumptions & Limitations
- Asteroid is approximated as a sphere with uniform density; spin axis & rate are assumed/estimated and held constant.
- Step 4 (close approach) uses fixed ΔV/TOF from literature; Step 1 (observation) uses placeholder prop/time.
- RCS configuration is simplified (paired thrusters, zero net translation) and mapped via a pseudo‑inverse; saturation/logic are idealized.
- SEP thrust in Step 2 is simplified from detailed low‑thrust optimal control; thrust is bounded by power & efficiency.
- Capture system mass is a first‑order parametric estimate from geometry and material densities (no detailed FEM or packaging).

---

## 9) Contributors
- Akira Jinkoji, Graduate Research Assistant, Aerospace System Design Laboratory (ASDL), Georgia Institude of Tehnology
- Salma Benhissoune, Graduate Research Assistant, Aerospace System Design Laboratory (ASDL), Georgia Institude of Tehnology
