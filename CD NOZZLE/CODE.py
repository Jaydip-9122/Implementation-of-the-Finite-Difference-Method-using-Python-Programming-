import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Constants
GAMMA = 1.4  # Specific heat ratio
R = 287.0    # Specific gas constant for air (J/(kg*K))
P0 = 1.0133e5  # Initial pressure (Pa)
T0 = 300.0   # Initial temperature (K)
Pe_P0 = 0.585   # Exit pressure ratio
TOL = 1e-3
Pe = Pe_P0 * P0  # Exit pressure (Pa)
CFL = 0.95   # Conservative CFL for stability
ITER_MAX = 20000    # Maximum number of iterations

# Nozzle geometry
NozzleLength = 2.0  # Length of the nozzle (m)
Imax = 101   # Number of grid points
x = np.linspace(0, NozzleLength, Imax)  # Grid points
dx = NozzleLength / (Imax - 1)  # Grid spacing
A_XLOC = 1.0 + 2.0 * (x - 1.0)**2    # Area distribution

# Initialization
MACHinitial = 0.15  # Initial Mach Number
TEMPERATURE = T0 / (1 + (GAMMA - 1)/2 * MACHinitial**2) * np.ones(Imax)
PRESSURE = P0 / (1 + (GAMMA - 1)/2 * MACHinitial**2)**(GAMMA/(GAMMA-1)) * np.ones(Imax)
DENSITY = PRESSURE / (R * TEMPERATURE)
VELOCITY = MACHinitial * np.sqrt(GAMMA * PRESSURE / DENSITY)

# Conservative or State Vector
U = np.zeros((3, Imax))
U[0] = DENSITY * A_XLOC  # Mass density
U[1] = DENSITY * VELOCITY * A_XLOC  # Momentum density
ENERGY = PRESSURE / (GAMMA - 1) + 0.5 * DENSITY * VELOCITY**2
U[2] = ENERGY * A_XLOC  # Energy density

# Initialize flux arrays
F_plus = np.zeros((3, Imax))
F_minus = np.zeros((3, Imax))

# Function to calculate shock location 
def SHOCK_LOCATION():
    def MACH_AREA_RELATION(M, A_ratio):
        return (1/M) * ((2/(GAMMA+1))*(1 + (GAMMA-1)/2 * M**2))**((GAMMA+1)/(2*(GAMMA-1))) - A_ratio
    
    def NSW(M1):
        return np.sqrt((1 + ((GAMMA - 1) / 2) * M1**2) / (GAMMA * M1**2 - (GAMMA - 1) / 2))
    
    X_VALUES = np.linspace(1.74, 1.78, 1000)
    BEST_XS, BST_ERROR = None, float('inf')
    
    for xs in X_VALUES:
        A_s = A_XLOC[np.argmin(np.abs(x - xs))]
        A_ratio = A_s / np.min(A_XLOC)
        
        try:
            M1 = fsolve(MACH_AREA_RELATION, 2.5, args=(A_ratio))[0]
            M2 = NSW(M1)
        except:
            continue
        
        # Pressure ratio before shock
        P1_P0 = (1 + ((GAMMA - 1) / 2) * M1**2)**(-GAMMA / (GAMMA - 1))
        P2_P1 = (2 * GAMMA * M1**2 - (GAMMA - 1)) / (GAMMA + 1)
        P2_P0 = P2_P1 * P1_P0
        
        # Exit conditions
        A_EXIT = A_XLOC[-1]
        A_EXIT_A2 = A_EXIT / A_s
        try:
            M_EXIT = fsolve(MACH_AREA_RELATION, 0.3, args=(A_EXIT_A2))[0]
        except:
            continue
        
        P_EXIT_P2 = (1 + ((GAMMA - 1) / 2) * M_EXIT**2)**(-GAMMA / (GAMMA - 1))
        PE_P0 = P_EXIT_P2 * P2_P0
        
        error = abs(PE_P0 - Pe_P0)
        if error < BST_ERROR:
            BEST_XS, BST_ERROR = xs, error
            best_M1, best_M2, best_pe_p0 = M1, M2, PE_P0
    
    return BEST_XS, best_M1, best_M2

# Calculate shock location
SHOCK_X, M1, M2 = SHOCK_LOCATION()

# Main loop
Iter_count = 0
ERROR = 1.0
MACH = []

while Iter_count < ITER_MAX and ERROR > TOL:    # loop for Van Leer method
    U_old = np.copy(U)
    P_Old = np.copy(PRESSURE)
    
    # Time step calculation
    a = np.sqrt(GAMMA * PRESSURE / DENSITY)
    lambda_max = np.max(np.abs(VELOCITY) + a)
    dt = CFL * dx / lambda_max
    
    # Boundary conditions
    VELOCITY[0] = VELOCITY[1]  # Extrapolate velocity at inlet

    PRESSURE[-1] = Pe  # Fixed pressure at exit
    DENSITY[-1] = DENSITY[-2]  # Extrapolate density at exit
    VELOCITY[-1] = VELOCITY[-2]  # Extrapolate velocity at exit
    
    # Update conservative variables at boundaries
    ENERGY = PRESSURE / (GAMMA - 1) + 0.5 * DENSITY * VELOCITY**2
    U[0, [0, -1]] = DENSITY[[0, -1]] * A_XLOC[[0, -1]]  # Mass conservation
    U[1, [0, -1]] = DENSITY[[0, -1]] * VELOCITY[[0, -1]] * A_XLOC[[0, -1]]  # Momentum conservation
    U[2, [0, -1]] = ENERGY[[0, -1]] * A_XLOC[[0, -1]]   # Energy conservation

    # Calculate fluxes for all points
    for i in range(Imax):
        density = U[0, i] / A_XLOC[i]
        velocity = U[1, i] / U[0, i]
        pressure = (GAMMA - 1) * (U[2, i]/A_XLOC[i] - 0.5 * density * velocity**2)
        sound_speed = np.sqrt(GAMMA * pressure / density)
        M = velocity / sound_speed
        
        # Calculate flux vector
        F = np.array([
            U[1, i],
            (density * velocity**2 + pressure) * A_XLOC[i],
            velocity * (U[2, i] + pressure * A_XLOC[i])
            ])
                
        # Van Leer splitting
        if M <= -1:
            F_plus[:, i] = np.zeros(3)
            F_minus[:, i] = F
        elif -1 < M < 1:
            alpha = 0.25 * density * sound_speed * (M + 1)**2 * A_XLOC[i]
            F_plus[:, i] = np.array([
                alpha,
                alpha * (2 * sound_speed / GAMMA) * (1 + (GAMMA - 1)/2 * M),
                alpha * (2 * sound_speed**2 / (GAMMA**2 - 1)) * (1 + (GAMMA - 1)/2 * M)**2
            ])
            F_minus[:, i] = F - F_plus[:, i]
            
        else:
            F_plus[:, i] = F
            F_minus[:, i] = np.zeros(3)

    # Update interior points with improved shock handling
    for i in range(1, Imax - 1):
        p_i = DENSITY[i] * a[i]**2 / GAMMA
        S = np.array([0, p_i * (A_XLOC[i+1] - A_XLOC[i-1])/(2*dx), 0])
        
        # Main iteration step
        U[:, i] = U[:, i] - dt/dx * (F_plus[:, i] - F_plus[:, i-1]) - dt/dx * (F_minus[:, i+1] - F_minus[:, i]) + dt * S
        
    
    # Update primitive variables
    DENSITY = U[0, :] / A_XLOC
    VELOCITY = U[1, :] / U[0, :]
    PRESSURE = (GAMMA - 1) * (U[2, :]/A_XLOC - 0.5 * DENSITY * VELOCITY**2)
    
    # Calculate Mach number properly
    M = VELOCITY / np.sqrt(GAMMA * PRESSURE / DENSITY)
    MACH.append(np.max(M))
    
    ERROR = np.max(np.abs(PRESSURE - P_Old))
    Iter_count += 1
    
    if Iter_count % 100 == 0:
        print(f"Iteration: {Iter_count}, Max Mach: {np.max(M):.3f}")

# Exact solution calculation
def area_mach_relation(M, A_ratio):
    return (1/M) * ((2/(GAMMA+1))*(1 + (GAMMA-1)/2 * M**2))**((GAMMA+1)/(2*(GAMMA-1))) - A_ratio

# Throat location
A_throat = np.min(A_XLOC)
throat_idx = np.argmin(A_XLOC)

# Shock index
shock_idx = np.argmin(np.abs(x - SHOCK_X))

# Primitive variables 
DENSITY = U[0, :] / A_XLOC
VELOCITY = U[1, :] / U[0, :]
PRESSURE = (GAMMA - 1) * (U[2, :]/A_XLOC - 0.5 * DENSITY * VELOCITY**2)
M = VELOCITY / np.sqrt(GAMMA * PRESSURE / DENSITY)
PRESSURE_P = PRESSURE / P0


M_exact = np.ones(Imax)
p_exact = np.ones(Imax) * P0

# Exact solution up to shock
for i in range(Imax):
    if x[i] < x[throat_idx]:
        # Subsonic region before throat
        A_ratio = A_XLOC[i] / A_throat
        try:
            M_exact[i] = fsolve(area_mach_relation, 0.3, args=(A_ratio))[0]
        except:
            M_exact[i] = 0.3
        p_exact[i] = P0 / (1 + (GAMMA-1)/2 * M_exact[i]**2)**(GAMMA/(GAMMA-1))
    elif x[i] < SHOCK_X:
        # Supersonic region after throat but before shock
        A_ratio = A_XLOC[i] / A_throat
        try:
            M_exact[i] = fsolve(area_mach_relation, 1.5, args=(A_ratio))[0]
        except:
            M_exact[i] = 1.5
        p_exact[i] = P0 / (1 + (GAMMA-1)/2 * M_exact[i]**2)**(GAMMA/(GAMMA-1))
    else:
        pass

if shock_idx < Imax - 1:
    # Calculate normal shock relations at the shock point
    p1 = P0 / (1 + (GAMMA-1)/2 * M1**2)**(GAMMA/(GAMMA-1))
    p2_p1 = (2 * GAMMA * M1**2 - (GAMMA - 1)) / (GAMMA + 1)
    p_shock = p1 * p2_p1
    
    # Set the shock point exactly
    p_exact[shock_idx] = p_shock
    M_exact[shock_idx] = M2
    
    p_ratio = p_exact[shock_idx] / PRESSURE[shock_idx]
    m_ratio = M_exact[shock_idx] / M[shock_idx]
    
    for i in range(shock_idx+1, Imax):
        p_exact[i] = PRESSURE[i] * p_ratio
        M_exact[i] = M[i] * m_ratio

# Final normalization of pressure
P_exact = p_exact / P0

# Plotting functions
def plot_grid(X, Y, exact1, color1, color2, xtitle, ytitle):
    plt.figure(figsize=(10 ,8))
    plt.plot(X, Y, color=color1, linewidth=3, label="Van Leer")
    plt.plot(X, exact1, color=color2, linewidth=3, linestyle="dotted", label="GD THEORY")
    plt.xlabel(xtitle, fontsize=12)
    plt.ylabel(ytitle, fontsize=12)
    plt.xlim(min(X), max(X))
    plt.grid(True, linestyle='-.', alpha=1)
    plt.legend(fontsize=10, loc='best')
    plt.tight_layout()
    plt.show()

#  Plots
plot_grid(x, PRESSURE_P, P_exact, 'magenta', 'black', "Nozzle Length (m)", "P/P0")
plot_grid(x, M, M_exact, 'silver', 'purple', "Nozzle Length (m)", "Mach Number")


