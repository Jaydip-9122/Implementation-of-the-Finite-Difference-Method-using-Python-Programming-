import numpy as np
import matplotlib.pyplot as plt

# Domain and Grid
LenX = 1.0  #   Length in x-direction
LenY = 1.0  #   Length in y-direction
NX = 31 #   Number of grid points in x-direction
NY = 31 #   Number of grid points in y-direction
delx = LenX / (NX - 1)  #   Grid spacing in x-direction
dely = LenY / (NY - 1)  #   Grid spacing in y-direction
# Physical Parameters
ReNo = 100.0
kinematicVisc = 1.0 / ReNo  #   Kinematic viscosity
# Stability Parameters
sigmaC = 0.4;sigmaD = 0.6
# Convergence
TOL = 1.0e-8;maxIter = 100000
UrmsHist = []  # RMS history for u velocity
VrmsHist = []  # RMS history for v velocity
# Initial Conditions
si = np.full((NX, NY), 100.0)
OMEGA = np.full((NX, NY),1*0.0)
U_velocity = np.full((NX, NY),1*0.0)
V_velocity = np.full((NX, NY),1*0.0)

for iter in range(maxIter):
    U = U_velocity.copy()
    V = V_velocity.copy()
    
    for i in range(NX):
    #   Top boundary
        OMEGA[i, NY-1] = -2 * (si[i, NY-2] - si[i, NY-1]) / dely**2 - 2 * 1.0 / dely
    #   Bottom boundary
        OMEGA[i, 0] = -2 * (si[i, 1] - si[i, 0]) / dely**2
    for j in range(NY):
    #   Left boundary 
        OMEGA[0, j] = -2 * (si[1, j] - si[0, j]) / delx**2
    #   Right boundary  
        OMEGA[NX-1, j] = -2 * (si[NX-2, j] - si[NX-1, j]) / delx**2

    # Time Step
    umax = np.max(np.abs(U_velocity))
    vmax = np.max(np.abs(V_velocity))
    dTC = sigmaC * delx * dely / (umax * dely + vmax * delx)    #   Convection Time Step
    dTD = sigmaD / (2 * kinematicVisc * (1/delx**2 + 1/dely**2))    #   Diffusion Time Step
    dT = min(dTC, dTD)  #   Time step size

    # Vorticity Transport
    Omega = OMEGA.copy()    
    for i in range(1, NX-1):
        for j in range(1, NY-1):
        # Second-order upwind in x-direction except special case @ i=0,1,NX-1,NX-2
            if U_velocity[i, j] > 0:
                if i>=2:
                    dwdx = U_velocity[i, j]*(3 * OMEGA[i, j] - 4 * OMEGA[i-1, j] + OMEGA[i-2, j]) / (2 * delx)
                else:
                    dwdx = U_velocity[i, j]*(OMEGA[i+1, j] - OMEGA[i-1, j]) / (2 * delx)
            else:
                if i<=NX-3: 
                    dwdx = U_velocity[i, j]*(-3 * OMEGA[i, j] + 4 * OMEGA[i+1, j] - OMEGA[i+2, j]) / (2 * delx)
                else:
                    dwdx = U_velocity[i, j]*(OMEGA[i+1, j] - OMEGA[i-1, j]) / (2 * delx)
        # Second-order upwind in y-direction except special case @ j=0,1,NY-1,NY-2
            if V_velocity[i, j] > 0:
                if j>=2:
                    dwdy = V_velocity[i, j]*(3 * OMEGA[i, j] - 4 * OMEGA[i, j-1] + OMEGA[i, j-2]) / (2 * dely)
                else:
                    dwdy = V_velocity[i, j]*(OMEGA[i, j+1] - OMEGA[i, j-1]) / (2 * dely)
            else:
                if j<=NY-3:
                    dwdy = V_velocity[i, j]*(-3 * OMEGA[i, j] + 4 * OMEGA[i, j+1] - OMEGA[i, j+2]) / (2 * dely)
                else:   
                    dwdy = V_velocity[i, j]*(OMEGA[i, j+1] - OMEGA[i, j-1]) / (2 * dely)
            diffusion = kinematicVisc * (
                (OMEGA[i, j+1] - 2*OMEGA[i, j] + OMEGA[i, j-1]) / dely**2 +
                (OMEGA[i+1, j] - 2*OMEGA[i, j] + OMEGA[i-1, j]) / delx**2)  #   Diffusion term
            
            Omega[i, j] = OMEGA[i, j] + dT * (-dwdx - dwdy + diffusion) #   Update vorticity
    OMEGA = Omega
    # Streamfunction Solution
    for i in range(maxIter):
        si_old = si.copy()
        si[1:-1, 1:-1] = 0.25 * (
            si[2:, 1:-1] + si[:-2, 1:-1] + si[1:-1, 2:] + si[1:-1, :-2] + delx**2 * OMEGA[1:-1, 1:-1])  
        if np.sqrt(np.mean((si - si_old)**2)) < 1e-2:
            break
    # Update velocities
    for i in range(1, NX-1):
        for j in range(1, si.shape[1] - 1):
            U_velocity[i][j] = (si[i][j + 1] - si[i][j - 1]) / (2 * dely)
            V_velocity[i][j] = -(si[i + 1][j] - si[i - 1][j]) / (2 * delx)
# Lid-driven top wall condition
    for i in range(NX):
        U_velocity[i][-1] = 1.0  # Lid velocity
        V_velocity[i][-1] = 0.0  # No penetration at the lid

    # Convergence check
    RMS_u = np.sqrt(np.mean((U_velocity - U)**2))
    RMS_v = np.sqrt(np.mean((V_velocity - V)**2))

    if iter % 1000 == 0:
        print(f"Iteration: {iter}, RMS_u: {RMS_u:.2e}, RMS_v: {RMS_v:.2e}")
    
    UrmsHist.append(RMS_u)
    VrmsHist.append(RMS_v)

    if RMS_u < TOL and RMS_v < TOL:
        print(f"Converged at iteration {iter}")
        break

import numpy as np
import matplotlib.pyplot as plt

def plot_contourf_streamfunction(psi, title='Stream Function Contour'):
    plt.figure(figsize=(8, 6))
    contour = plt.contourf(psi.T, levels=50, cmap='gist_rainbow')
    plt.colorbar(contour, label='Stream Function (ψ)')
    plt.xlabel('i')
    plt.ylabel('j')
    plt.title(title)
    plt.tight_layout()
    plt.show()

def plot_streamlines(U, V, LenX, LenY, NX, NY):
    x = np.linspace(0, LenX, NX)
    y = np.linspace(0, LenY, NY)
    X, Y = np.meshgrid(x, y)
    
    plt.figure(figsize=(8, 6))
    plt.streamplot(X, Y, U.T, V.T, color='green', density=2, linewidth=1, arrowsize=1)
    plt.title("Streamlines of Lid-Driven Cavity Flow")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.show()

def plot_vertical_centerline_velocity(U, y, center_i, filename='mid verticle.txt'):
    data = np.loadtxt(filename, skiprows=1)
    y_ref, u_ref = data[:, 0], data[:, 1]

    plt.figure(figsize=(5, 4))
    plt.plot(y, U[center_i, :], '--', label='Simulation')
    plt.plot(y_ref, u_ref, 'green', label='Reference Data')
    plt.xlabel("u velocity")
    plt.ylabel("y")
    plt.title("Vertical Centerline Velocity (u vs y)")
    plt.grid(True)
    plt.legend()
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()

def plot_horizontal_centerline_velocity(V, x, center_j, filename='mid horizontal.txt'):
    data = np.loadtxt(filename, skiprows=1)
    x_ref, v_ref = data[:, 0], data[:, 1]

    plt.figure(figsize=(6, 4))
    plt.plot(x, V[:, center_j], '--', label='Simulation')
    plt.plot(x_ref, v_ref, 'green', label='Reference Data')
    plt.xlabel("x")
    plt.ylabel("v velocity")
    plt.title("Horizontal Centerline Velocity (v vs x)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_rms_convergence(Urms, Vrms):
    plt.figure(figsize=(8, 5))
    plt.semilogy(Urms, label='RMS_u')
    plt.semilogy(Vrms, label='RMS_v')
    plt.xlabel("Iteration")
    plt.ylabel("RMS Value (log scale)")
    plt.title("RMS Convergence History")
    plt.legend()
    plt.tight_layout()
    plt.show()

x = np.linspace(0, LenX, NX)
y = np.linspace(0, LenY, NY)
center_i = NX // 2
center_j = NY // 2
plot_contourf_streamfunction(si)
plot_streamlines(U_velocity, V_velocity, LenX, LenY, NX, NY)
plot_vertical_centerline_velocity(U_velocity, y, center_i)
plot_horizontal_centerline_velocity(V_velocity, x, center_j)
plot_rms_convergence(UrmsHist, VrmsHist)
