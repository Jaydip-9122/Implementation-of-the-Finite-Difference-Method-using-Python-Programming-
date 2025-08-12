import numpy as np 
import matplotlib.pyplot as plt

L = 1.0     #   DOMAIN LENGTH
N = 101     #   NO. OF GRID POINTS
h = L/(N-1) #   STEP LENGTH
a = 1       #   WAVE SPEED
T = 0.35  #   TOTAL SIMULATION TIME
x = np.linspace(0,L,N)  #   DOMAIN
u = np.zeros(N)     #   INITIALIZING
CFL_NUMBER = float(input("Please Enter the value of CFL Number you want to go for: "))

def TDMA(lower, main, upper, N):
    
    AMATRIX = np.zeros((N, N))
    np.fill_diagonal(AMATRIX, main)
    np.fill_diagonal(AMATRIX[:-1, 1:], upper)  # Upper diagonal,exception of last row and first column
    np.fill_diagonal(AMATRIX[1:, :-1], lower)  # Lower diagonal,exception of last column and first row
    return AMATRIX

lower = -CFL_NUMBER / 2
main = 1 
upper = CFL_NUMBER / 2
A = TDMA(lower, main, upper, N)

def BTCS(A,u,CFL,NO_t,BOUNDARY_CONDITION):        
        for n in range(NO_t):            
            u_current = u.copy()
              # Boundary conditions
        
            u_current[0] = BOUNDARY_CONDITION[0]
            u_current[-1] =BOUNDARY_CONDITION[-1]

            u_NEW = np.linalg.solve(A, u_current)  # Solve the linear system
             
            u = u_NEW     
    
        return u

INITIAL_CONDITION = input("Please Enter the INITIAL CONDITION form giving list [FIRST,SECOND,THIRD,FOURTH,FIFTH]: ")
if INITIAL_CONDITION is None:
    raise ValueError(f"Invalid Initial Condition '{INITIAL_CONDITION}'. Please choose from [FIRST,SECOND,THIRD,FOURTH,FIFTH].")


def plot_scheme(x, u1,u2, T, CFL, scheme_name, initial_condition, xlabel="x", ylabel="u", xlim=(0, 1), ylim=(-1.5, 1.5)):
    plt.figure(figsize=(10, 4))
    plt.plot(x, u1, label=f'Time = {T:.2f}', color='magenta')  # Magenta line color
    plt.plot(x, u2, label="Initial condition,T = 0", color='blue',linestyle = "--")
    # Generate the title
    title = f"{scheme_name} Scheme, CFL = {CFL_NUMBER}, {INITIAL_CONDITION} INITIAL CONDITION"
    plt.title(title);plt.xlabel(xlabel);plt.ylabel(ylabel);plt.xlim(*xlim);plt.ylim(*ylim);plt.legend();
     
    # # Save the plot
    # filename = title.replace(" ", "_").replace(",", "").replace("=", "_") + ".png"
    # plt.savefig(filename, dpi=300, bbox_inches='tight') 
    plt.show()

if INITIAL_CONDITION == "FIRST":
    for i in range(N):
        if x[i] >= 0.2:
            u[i] = 1
        else:
            u[i] = 0
    
    BOUNDARY_CONDITION = (0,1)

elif INITIAL_CONDITION == "SECOND":
    for i in range(N):
        if 0 <= x[i] <= 0.05:
            u[i] = 0
        elif 0.05 <= x[i] <= 0.35:
            u[i] = np.sin(4 * np.pi * ((x[i] - 0.05) / 0.3))
        elif 0.35 <= x[i] <= 1:
            u[i] = 0
    BOUNDARY_CONDITION = (0,0)

elif INITIAL_CONDITION == "THIRD":
    for i in range(N):
        if 0 <= x[i] <= 0.05:
            u[i] = 0
        elif 0.05 <= x[i] <= 0.35:
            u[i] = np.sin(8 * np.pi * ((x[i] - 0.05) / 0.3))
        elif 0.35 <= x[i] <= 1:
            u[i] = 0
    BOUNDARY_CONDITION = (0,0,0,0)

elif INITIAL_CONDITION == "FOURTH":
    for i in range(N):
        if 0 <= x[i] <= 0.05:
            u[i] = 0
        elif 0.05 <= x[i] <= 0.35:
            u[i] = np.sin(12 * np.pi * ((x[i] - 0.05) / 0.3))
        elif 0.35 <= x[i] <= 1:
            u[i] = 0
    BOUNDARY_CONDITION = (0,0,0,0)

elif INITIAL_CONDITION == "FIFTH":
    for i in range(N):
        u[i]= np.exp((-50/(0.4**2)*((x[i]-0.2)**2)))

    BOUNDARY_CONDITION = (0,0,0,0)

del_T = CFL_NUMBER * h / a     #   DELTA T IN CFL NUMBER
NO_t = int(T / del_T)
u_BTCS = BTCS(A,u, CFL_NUMBER, NO_t,BOUNDARY_CONDITION)
plot_scheme(x, u_BTCS, u, T, CFL_NUMBER, "BTCS", INITIAL_CONDITION)