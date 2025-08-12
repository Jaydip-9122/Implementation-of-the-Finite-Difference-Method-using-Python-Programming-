# Lid-Driven Cavity Flow Solver (Vorticity-Streamfunction Formulation)

This assignment implements a 2D incompressible lid-driven cavity flow solver using the vorticity-streamfunction formulation. The method uses finite difference techniques and solves the vorticity transport equation coupled with a Poisson equation for the streamfunction.

## Features

- Solves 2D steady lid-driven cavity problem using:
  - Vorticity transport equation (explicit time integration)
  - Poisson equation for streamfunction (Jacobi iteration)
- Uses upwind differencing for convection terms
- Variable time step based on convection and diffusion stability criteria
- Visualizes streamfunction contours and streamlines
- Compares centerline velocity profiles with benchmark data
- Tracks RMS error convergence for u and v velocity fields

## Grid and Domain

- Domain: Square cavity \[0,1\] x \[0,1\]
- Grid size: 31x31 points
- Grid spacing computed as:
  - `delx = LenX / (NX - 1)`
  - `dely = LenY / (NY - 1)`

## Parameters

- Reynolds number: `Re = 100`
- Kinematic viscosity: `ν = 1 / Re`
- Time-stepping:
  - CFL condition-based for convection (`σc = 0.4`)
  - Diffusion-based condition (`σd = 0.6`)
- Lid velocity: `u = 1.0` on top wall

## Files Required

Make sure to include the following data files for comparison:
- `mid verticle.txt` – Vertical centerline (u vs y)
- `mid horizontal.txt` – Horizontal centerline (v vs x)

The files should be in the format:

## Output Plots

- Streamfunction contour (`ψ`)
- Streamlines using `streamplot()`
- Comparison of simulation vs benchmark data:
  - Vertical line (u vs y)
  - Horizontal line (v vs x)
- RMS convergence of velocity fields

## How to Run

Simply execute the Python script using:

```bash
python lid_cavity.py
you need text file for reference plot.

Title
Lid-Driven Cavity Flow Solver (Vorticity-Streamfunction)

Description
A 2D incompressible lid-driven cavity flow simulator using the vorticity-streamfunction formulation. Solves steady-state flow using finite difference methods and explicit time stepping, with visualization and benchmark comparisons.

Features
- Solves vorticity-streamfunction equations for 2D incompressible flow
- Explicit time integration for vorticity transport
- Jacobi iteration for streamfunction (Poisson solver)
- Adaptive time step based on convection and diffusion stability
- Visualization of streamfunction and streamlines
- Velocity profile comparisons with benchmark data (Ghia et al.)
- RMS convergence tracking

Installation
1. Clone the repository or download the script
2. Ensure required data files are in the same directory:
   - `mid verticle.txt`
   - `mid horizontal.txt`
3. Install dependencies:
   - `pip install numpy matplotlib`

Usage
Run the Python script:
python lid_cavity.py

Make sure the benchmark `.txt` files are present for centerline comparison plots.

Input Files
- `mid verticle.txt`: Vertical centerline data (u vs y)
- `mid horizontal.txt`: Horizontal centerline data (v vs x)
Each file should have two columns: coordinate and reference velocity.

Output
- Contour plot of streamfunction (ψ)
- Streamline plot of the velocity field
- Comparison plots:
  - u velocity along vertical centerline
  - v velocity along horizontal centerline
- Log-scale RMS convergence history for u and v

Parameters
- Domain: 1.0 x 1.0 square
- Grid: 31 x 31 uniform points
- Reynolds number: Re = 100
- Lid velocity: u = 1.0 (top wall)
- CFL numbers: σc = 0.4, σd = 0.6
- Tolerance: 1e-8

Author
Jaydip Patel
24M0004@iitb.ac.in



