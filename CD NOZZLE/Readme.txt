1D Quasi-One-Dimensional Nozzle Flow Simulation with Shock Handling
===================================================================

Overview:
---------
This Python script simulates compressible, quasi-one-dimensional flow through a converging-diverging nozzle using the finite volume method. The simulation accounts for shock waves and applies Van Leer flux vector splitting for flux computation. It includes both numerical and analytical (exact) solutions for comparison.

Main Features:
--------------
- Conservative form of 1D Euler equations
- Van Leer flux vector splitting scheme
- Time marching with dynamic CFL-based time step
- Boundary conditions: subsonic inlet & fixed pressure outlet
- Shock location estimation using exact area-Mach relations and normal shock equations
- Smoothing around shock region to improve numerical stability
- Comparison with analytical solution using area-Mach relation
- Visualization of:
  - Pressure distribution (P/P₀)
  - Mach number along the nozzle length

Key Parameters:
---------------
- γ = 1.4 (Ratio of specific heats)
- R = 287 J/kg·K (Gas constant for air)
- P₀ = 101330 Pa (Total inlet pressure)
- T₀ = 300 K (Total inlet temperature)
- Pe/P₀ = 0.585 (Exit-to-inlet pressure ratio)
- Nozzle length = 2 m
- Grid points = 101
- CFL = 0.95
- Tolerance = 1e-4
- Maximum Iterations = 50,000

Nozzle Geometry:
----------------
Area distribution is defined as:

    A(x) = 1.0 + 2.0 * (x - 1.0)^2

This gives a symmetric converging-diverging nozzle with a throat at x = 1.0 m.

How It Works:
-------------
1. Initializes the flow field with subsonic conditions.
2. Iteratively updates conservative variables using flux vector splitting.
3. Estimates shock location by matching exit pressure ratio.
4. Applies numerical smoothing around the detected shock.
5. Computes exact solution using area-Mach and normal shock relations.
6. Compares numerical and exact solutions visually.

Output:
-------
- Console output: Iteration count, maximum Mach number, convergence error
- Plots:
  - Pressure ratio (P/P₀) vs nozzle length
  - Mach number vs nozzle length

Dependencies:
-------------
- numpy
- matplotlib
- scipy

Author: Patel Jaydip



