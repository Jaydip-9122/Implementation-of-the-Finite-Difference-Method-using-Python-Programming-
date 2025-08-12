# Finite Difference Method (FDM) Wave Equation Solver

This program solves a 1D wave equation using various finite difference methods (FDM) schemes. It supports user-defined initial conditions and CFL numbers for simulation. The implemented FDM schemes include FTFS, FTCS, FTBS, LW, BW, FR, and BTCS.

---

## Features
1. **Finite Difference Schemes:**
   - Forward-Time Forward-Space (FTFS)
   - Forward-Time Central-Space (FTCS)
   - Forward-Time Backward-Space (FTBS)
   - Lax-Wendroff (LW)
   - Beam-Warming (BW)
   - Fromm’s Scheme (FR)
   - Backward-Time Central-Space (BTCS)

2. **Initial Conditions:**
   - Five pre-defined initial conditions (FIRST, SECOND, THIRD, FOURTH, FIFTH).

3. **Customizable Parameters:**
   - CFL number
   - Boundary conditions

4. **Visualization:**
   - Plots initial and final states of the wave.

---

## How to Run the Code

### Prerequisites
Ensure you have the following Python libraries installed:
- `numpy`
- `matplotlib`

You can install them using:
```bash
pip install numpy matplotlib
```

### Steps to Run
1. Save the script to a file, e.g., `24M0004_EXPLICIT.py` or `24M0004_IMPLICIT.py`.
2. Open a terminal or IDE and navigate to the script's directory.
3. Run the script:
   ```bash
   python 24M0004_EXXPLICIT.py
   ```
4. Follow the on-screen prompts to select:
   - FDM scheme (e.g., FTFS, FTCS, FTBS, LW, BW, FR, BTCS).
   - Initial condition (e.g., FIRST, SECOND, etc.).
   - CFL number.

### Example Input
When prompted:
- **FDM Scheme:** `FTFS`
- **Initial Condition:** `FIRST`
- **CFL Number:** `0.5`

### Output
- A plot showing the initial condition (dashed blue line) and the final state of the wave (magenta line).

---

## Code Walkthrough

### Parameters
- `L`: Domain length (default: 1.0).
- `N`: Number of grid points (default: 101).
- `a`: Wave speed (default: 1).
- `T`: Total simulation time (default: 0.35).
- `CFL`: Courant-Friedrichs-Lewy number.

### Functions
- **`FTFS`**: Implements the Forward-Time Forward-Space scheme.
- **`FTCS`**: Implements the Forward-Time Central-Space scheme.
- **`FTBS`**: Implements the Forward-Time Backward-Space scheme.
- **`LW`**: Implements the Lax-Wendroff scheme.
- **`BW`**: Implements the Beam-Warming scheme.
- **`FR`**: Implements Fromm’s scheme.
- **`BTCS`**: Implements the Backward-Time Central-Space scheme.

### Initial Conditions
- **FIRST**: Step function.
- **SECOND, THIRD, FOURTH**: Sine wave packets with increasing frequency.
- **FIFTH**: Gaussian distribution.

### Example Outputs
#### FTFS Scheme with FIRST Initial Condition and CFL = 0.5
```text
FDM Scheme: FTFS
Initial Condition: FIRST
CFL Number: 0.5
```
**Output:** A plot displaying the wave evolution.

---

## Contact
For further assistance, reach out to **Jaydip Patel** at [24m0004@iitb.ac.in](mailto:24m0004@iitb.ac.in).



