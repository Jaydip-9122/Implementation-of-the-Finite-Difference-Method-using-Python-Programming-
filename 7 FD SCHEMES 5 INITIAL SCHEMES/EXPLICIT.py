import numpy as np 
import matplotlib.pyplot as plt

L = 1.0     #   DOMAIN LENGTH
N = 101     #   NO. OF GRID POINTS
h = (L/(N-1)) #   STEP LENGTH
a = 1       #   WAVE SPEED
T = 0.35  #   TOTAL SIMULATION TIME
x = np.linspace(0,L,N)  #   DOMAIN
u = np.zeros(N)     #   INITIALIZING

def FTFS(u, CFL, T, h, a,BOUNDARY_CONDITION):
    u_current = np.copy(u)

    for n in range(NO_t):
        u_new = np.copy(u_current)
        for i in range(1, N - 1):
            u_new[i] = u_current[i] - CFL * (u_current[i + 1] - u_current[i])
        u_current = u_new

        # Boundary Conditions
        u_current[0] = BOUNDARY_CONDITION[0]
        u_current[-1] = BOUNDARY_CONDITION[-1]

    return u_current
 
def FTCS(u, CFL, T, h, a,BOUNDARY_CONDITION):
    u_current = np.copy(u)
    for n in range(NO_t):
        u_new = np.copy(u_current)
        for i in range(1, N - 1):
            u_new[i] = u_current[i] - (CFL/2) * (u_current[i + 1] - u_current[i])
        u_current = u_new

        # Boundary Conditions
        u_current[0] = BOUNDARY_CONDITION[0]
        u_current[-1] = BOUNDARY_CONDITION[-1]

    return u_current

def FTBS(u, CFL, T, h, a,BOUNDARY_CONDITION):
    u_current = np.copy(u)

    for n in range(NO_t):
        u_new = np.copy(u_current)
        for i in range(1, N - 1):
            u_new[i] = u_current[i] - CFL * (u_current[i] - u_current[i-1])
        u_current = u_new

        # Boundary Conditions
        u_current[0] = BOUNDARY_CONDITION[0]
        u_current[-1] = BOUNDARY_CONDITION[-1]

    return u_current

def LW(u, CFL, T, h, a,BOUNDARY_CONDITION):
    u_current = np.copy(u)

    for n in range(NO_t):
        u_new = np.copy(u_current)
        for i in range(1, N - 1):
            u_new[i] = u_current[i] - (CFL/2) * (u_current[i+1] - u_current[i-1]) +  ((CFL**2)/2) * (u_current[i+1] - 2*u_current[i]+ u_current[i-1])
        u_current = u_new

        # Boundary Conditions
        u_current[0] = BOUNDARY_CONDITION[0]
        u_current[-1] = BOUNDARY_CONDITION[-1]
    return u_current

def BW(u, CFL, T, h, a,BOUNDARY_CONDITION):
    u_current = np.copy(u)

    for n in range(NO_t):
        u_new = np.copy(u_current)
        for i in range(2, N - 2):
            u_new[i] = u_current[i] - (CFL/2) * (3*u_current[i] - 4*u_current[i-1] +u_current[i-2]) +  ((CFL**2)/2) * (u_current[i] - 2*u_current[i-1]+ u_current[i-2])
        u_current = u_new

        # Boundary Conditions
        u_current[0] = BOUNDARY_CONDITION[0]
        u_current[1] = BOUNDARY_CONDITION[1]
        u_current[-1] = BOUNDARY_CONDITION[-1]
        u_current[-2] = BOUNDARY_CONDITION[-2]
    return u_current

def FR(u, CFL, T, h, a,BOUNDARY_CONDITION):
    u_current = np.copy(u)

    for n in range(NO_t):
        u_new = np.copy(u_current)
        for i in range(2, N - 2):
            u_new[i] = ((u_current[i] - (CFL/2) * (3*u_current[i] - 4*u_current[i-1] +u_current[i-2]) +  ((CFL**2)/2) * (u_current[i] - 2*u_current[i-1]+ u_current[i-2]))+(u_current[i] - (CFL/2) * (u_current[i+1] - u_current[i-1]) +  ((CFL**2)/2) * (u_current[i+1] - 2*u_current[i]+ u_current[i-1])))/2
        u_current = u_new

        # Boundary Conditions
        u_current[0] = BOUNDARY_CONDITION[0]
        u_current[1] = BOUNDARY_CONDITION[1]
        u_current[-1] = BOUNDARY_CONDITION[-1]
        u_current[-2] = BOUNDARY_CONDITION[-2]
    return u_current

def plot_scheme(x, u,u2, T, CFL, scheme_name, initial_condition, xlabel="x", ylabel="u", xlim=(0, 1), ylim=(-1.5, 1.5)):
    plt.figure(figsize=(10, 4))
    plt.plot(x, u, label=f'Time = {T:.2f}', color='magenta')  # Magenta line color
    plt.plot(x, u2, label="Initial condition,T = 0", color='blue',linestyle = "--")
    # Generate the title
    title = f"{(str(SCHEME))} SCHEME, CFL = {CFL_NUMBER}, {INITIAL_CONDITION} INITIAL_CONDITION"
    plt.title(title);plt.xlabel(xlabel);plt.ylabel(ylabel);plt.xlim(*xlim);plt.ylim(*ylim);plt.legend()
    
    # # Save the plot
    # filename = title.replace(" ", "_").replace(",", "").replace("=", "_") + ".png"
    # plt.savefig(filename, dpi=300, bbox_inches='tight') 
    plt.show()

SCHEME_DICTIONARY = {
    "FTFS": FTFS,
    "FTCS": FTCS,
    "FTBS": FTBS,
    "LW": LW,
    "BW": BW,
    "FR": FR
}
SCHEME = input("Please Enter the FDM scheme form giving list [FTFS,FTCS,FTBS,LW,BW,FR]: ")
SCHEME_NAME = SCHEME_DICTIONARY.get(SCHEME)
if SCHEME_NAME is None:
    raise ValueError(f"Invalid scheme name '{SCHEME}'. Please choose from [FTFS, FTCS, FTBS, LW, BW, FR].")

INITIAL_CONDITION = input("Please Enter the INITIAL CONDITION form giving list [FIRST,SECOND,THIRD,FOURTH,FIFTH]: ")
if INITIAL_CONDITION is None:
    raise ValueError(f"Invalid Initial Condition '{INITIAL_CONDITION}'. Please choose from [FIRST,SECOND,THIRD,FOURTH,FIFTH].")
CFL_NUMBER = float(input("Please Enter the value of CFL Number you want to go for: "))

if INITIAL_CONDITION == "FIRST":
    for i in range(N):
        if x[i] >= 0.2:
            u[i] = 1
        else:
            u[i] = 0
    if  SCHEME == "FTFS" or"FTCS" or "FTBS" or "LW":
        BOUNDARY_CONDITION = (0,1)
    elif SCHEME == "BR" or"FR":
        BOUNDARY_CONDITION = (0,0,1,1)

elif INITIAL_CONDITION == "SECOND":
    for i in range(N):
        if 0 <= x[i] <= 0.05:
            u[i] = 0
        elif 0.05 <= x[i] <= 0.35:
            u[i] = np.sin(4 * np.pi * ((x[i] - 0.05) / 0.3))
        elif 0.35 <= x[i] <= 1:
            u[i] = 0
    if  SCHEME == "FTFS" or"FTCS" or "FTBS" or "LW":
        BOUNDARY_CONDITION = (0,0)
    elif SCHEME == "BR" or"FR":
        BOUNDARY_CONDITION = (0,0,0,0)

elif INITIAL_CONDITION == "THIRD":
    for i in range(N):
        if 0 <= x[i] <= 0.05:
            u[i] = 0
        elif 0.05 <= x[i] <= 0.35:
            u[i] = np.sin(8 * np.pi * ((x[i] - 0.05) / 0.3))
        elif 0.35 <= x[i] <= 1:
            u[i] = 0
    if SCHEME == "FTFS" or"FTCS" or "FTBS" or "LW":
        BOUNDARY_CONDITION = (0,0)
    elif SCHEME == "BR" or"FR":
        BOUNDARY_CONDITION = (0,0,0,0)

elif INITIAL_CONDITION == "FOURTH":
    for i in range(N):
        if 0 <= x[i] <= 0.05:
            u[i] = 0
        elif 0.05 <= x[i] <= 0.35:
            u[i] = np.sin(12 * np.pi * ((x[i] - 0.05) / 0.3))
        elif 0.35 <= x[i] <= 1:
            u[i] = 0
    if SCHEME == "FTFS" or"FTCS" or "FTBS" or "LW":
        BOUNDARY_CONDITION = (0,0)
    elif SCHEME == "BR" or"FR":
        BOUNDARY_CONDITION = (0,0,0,0)

elif INITIAL_CONDITION == "FIFTH":
    for i in range(N):
        u[i]= np.exp((-50/(0.4**2)*((x[i]-0.2)**2)))

    if SCHEME == "FTFS" or"FTCS" or "FTBS" or "LW":
        BOUNDARY_CONDITION = (0,0)
    elif SCHEME == "BR" or"FR":
        BOUNDARY_CONDITION = (0,0,0,0)

del_T = CFL_NUMBER * h / a     #   DELTA T IN CFL NUMBER
NO_t = int(T / del_T)
u_FINAL = SCHEME_NAME(u, CFL_NUMBER, T, h, a, BOUNDARY_CONDITION)
plot_scheme(x, u_FINAL,u, T, CFL_NUMBER, "FTFS", INITIAL_CONDITION)