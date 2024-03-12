import numpy as np

a = 0
b = 1
N = 100
h = (b - a)/N
x_node = np.linspace(a, b, N + 1)
x_cell = np.linspace(a + h/2, b - h/2, N)
x_disc = (a + b)/2

T = 1
time = 0
CFL = 0.3
gamma = 1.4

U_cons = np.zeros((N, 4))
U_flux = np.zeros((N + 1, 4))

#U_cons[i] == U_flux[i + 1/2]
#U = (rho, u, P, E)


#%% функции

def Riemann(U_L, U_R):  #U = (rho, u, P, E),  U.shape = 4
    
    rho_L = U_L[0]
    rho_R = U_R[0]
    u_L   = U_L[1]
    u_R   = U_R[1]
    P_L   = U_L[2]
    P_R   = U_R[2]
    E_L   = U_L[3]
    E_R   = U_R[3]
    
    H_L = E_L + P_L/rho_L
    H_R = E_R + P_R/rho_R
    
    u_RL = (rho_R**0.5*u_R + rho_L**0.5*u_L) / (rho_L**0.5 + rho_R**0.5)
    H_RL = (rho_R**0.5*H_R + rho_L**0.5*H_L) / (rho_L**0.5 + rho_R**0.5)
    rho_RL = (rho_L*rho_R)**0.5
    c_RL = ((gamma - 1)*(H_RL - 0.5*u_RL**2))**0.5
    
    lambda1 = u_RL
    lambda2 = u_RL + c_RL
    lambda3 = u_RL - c_RL
    
    dv1 = rho_R - rho_L - (P_R - P_L)/c_RL**2
    dv2 = u_R - u_L + (P_R - P_L)/(c_RL*rho_RL)
    dv3 = u_R - u_L - (P_R - P_L)/(c_RL*rho_RL)
    
    r1 = np.array([1,  u_RL,  0.5*u_RL**2])
    r2 = np.array([1,  u_RL + c_RL,  H_RL + c_RL*u_RL])* rho_RL/(2*c_RL)
    r3 = np.array([1,  u_RL - c_RL,  H_RL - c_RL*u_RL])*-rho_RL/(2*c_RL)

    if 0 <= lambda3:
        U = U_L.copy()
    elif lambda3 <= 0 <= lambda1:
        rho = rho_L - rho_RL/(2*c_RL)*dv3
        u = 1/rho*(rho_L*u_L - rho_RL/(2*c_RL)*(u_RL - c_RL)*dv3)
        E = 1/rho*(rho_L*E_L - rho_RL/(2*c_RL)*(H_RL - c_RL*u_RL)*dv3)
        P = (gamma - 1)*rho*(E - 0.5*u**2)
        U = np.array([rho, u, P, E])
    elif lambda1 <= 0 <= lambda2:
        rho = rho_R - rho_RL/(2*c_RL)*dv2
        u = 1/rho*(rho_R*u_R - rho_RL/(2*c_RL)*(u_RL + c_RL)*dv2)
        E = 1/rho*(rho_R*E_R - rho_RL/(2*c_RL)*(H_RL + c_RL*u_RL)*dv2)
        P = (gamma - 1)*rho*(E - 0.5*u**2)
        U = np.array([rho, u, P, E])
    else:
        U = U_R.copy()
    
    return U


def tau_update(U_cons):  #U[i] = (rho, u, P, E), U.shape = (N, 4)
    rho_cons = U_cons[:, 0]
    u_cons   = U_cons[:, 1]
    P_cons   = U_cons[:, 2]
    c_cons = np.sqrt(gamma*P_cons/rho_cons)
    return CFL*h/np.max(np.abs(u_cons) + c_cons)


def flux_generator(U_cons, U_L, U_R):  #U[i] = (rho, u, P, E), U.shape = (N, 4)
    U_flux = np.zeros((N+1, 4))
    for i in range(1, N):
        U_flux[i] = Riemann(U_cons[i-1], U_cons[i])
    U_flux[0] = U_L
    U_flux[N] = U_R
    return U_flux


def corrector(U_cons, U_flux):  #U[i] = (rho, u, P, E), U.shape = (N+1, 4)
    rho_cons = U_cons[:, 0]
    u_cons   = U_cons[:, 1]
    P_cons   = U_cons[:, 2]
    E_cons   = U_cons[:, 3]
    rho_flux = U_flux[:, 0]
    u_flux   = U_flux[:, 1]
    P_flux   = U_flux[:, 2]
    E_flux   = U_flux[:, 3]
    
    for i in range(N):
        rho_new = rho_cons[i] - tau/h*(rho_flux[i+1]*u_flux[i+1] - rho_flux[i]*u_flux[i])
        
        u_cons[i] = (rho_cons[i]*u_cons[i] - tau/h*(rho_flux[i+1]*u_flux[i+1]**2 \
                                                  - rho_flux[ i ]*u_flux[ i ]**2 \
                        + P_flux[i+1] - P_flux[i]))/rho_new
        
        E_cons[i] = (rho_cons[i]*E_cons[i] - tau/h*(rho_flux[i+1]*E_flux[i+1]*u_flux[i+1] \
                                                  - rho_flux[ i ]*E_flux[ i ]*u_flux[ i ] \
                        + P_flux[i+1]*u_flux[i+1] - P_flux[i]*u_flux[i]))/rho_new
        
        rho_cons[i] = rho_new
        
        P_cons[i] = (gamma - 1)*rho_cons[i]*(E_cons[i] - 0.5*u_cons[i]**2)
    

def export(U_cons, step, time):
    import os
    dir_name = 'out' 
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    file_name = dir_name + '\\out_' + str(step) + '.dat'
    with open(file_name, 'w') as f:
        
        # rho_cons
        f.write("TITLE=\"OUT_CONSVT\"\n")
        f.write("VARIABLES = \"X\", \"rhoc\"\n")
        f.write("ZONE\n")
        f.write("T = CabLeftNoSonic_FluxCorr_rhoc, I = {}\n".format(N))
        f.write("ZONETYPE = ORDERED, DATAPACKING = POINT, C = BLUE\n")
        f.write("STRANDID = 1\n")
        f.write("SOLUTIONTIME = {}\n".format(time))
        f.write("\n")
        for i in range(N):
            f.write(f"{x_cell[i]} {U_cons[i][0]}\n")
        f.write("\n")
        
        # u_cons
        f.write("TITLE=\"OUT_CONSVT\"\n")
        f.write("VARIABLES = \"X\", \"uc\"\n")
        f.write("ZONE\n")
        f.write("T = CabLeftNoSonic_FluxCorr_uc, I = {}\n".format(N))
        f.write("ZONETYPE = ORDERED, DATAPACKING = POINT, C = BLUE\n")
        f.write("STRANDID = 2\n")
        f.write("SOLUTIONTIME = {}\n".format(time))
        f.write("\n")
        for i in range(N):
            f.write(f"{x_cell[i]} {U_cons[i][1]}\n")
        f.write("\n")
        
        # P_cons
        f.write("TITLE=\"OUT_CONSVT\"\n")
        f.write("VARIABLES = \"X\", \"Pc\"\n")
        f.write("ZONE\n")
        f.write("T = CabLeftNoSonic_FluxCorr_Pc, I = {}\n".format(N))
        f.write("ZONETYPE = ORDERED, DATAPACKING = POINT, C = BLUE\n")
        f.write("STRANDID = 3\n")
        f.write("SOLUTIONTIME = {}\n".format(time))
        f.write("\n")
        for i in range(N):
            f.write(f"{x_cell[i]} {U_cons[i][2]}\n")
        f.write("\n")


#%% начальные данные

rho_L = 1
u_L   = 0
P_L   = 2
E_L = 0.5*u_L**2 + P_L/(rho_L*(gamma-1))

rho_R = 1
u_R   = 0
P_R   = 1
E_R   = 0.5*u_R**2 + P_R/(rho_R*(gamma-1))

U_L = np.array([rho_L, u_L, P_L, E_L])
U_R = np.array([rho_R, u_R, P_R, E_R])

for i in range(N):
    if x_cell[i] < x_disc:
        U_cons[i][0] = rho_L
        U_cons[i][1] = u_L
        U_cons[i][2] = P_L
    else:
        U_cons[i][0] = rho_R
        U_cons[i][1] = u_R
        U_cons[i][2] = P_R
    
#E_cons = 0.5*u**2 + P/(rho*(gamma-1))
U_cons[:, 3] = 0.5*U_cons[:,1]**2 + U_cons[:,2]/(U_cons[:,0]*(gamma-1))


#%% цикл по времени

time = 0
step = 0

export(U_cons, step, time)

#while time < T:
for _ in range(120):
    tau = tau_update(U_cons)
    U_flux = flux_generator(U_cons, U_L, U_R)
    corrector(U_cons, U_flux)
    time += tau
    step += 1
    export(U_cons, step, time)
    

#%% графики консервативных

import matplotlib.pyplot as plt

plt.grid(True)
plt.plot(x_cell, U_cons[:, 0], label='rho_cons')
plt.plot(x_cell, U_cons[:, 1], label='u_cons')
plt.plot(x_cell, U_cons[:, 2], label='P_cons')
plt.plot(x_cell, U_cons[:, 3], label='E_cons')
plt.legend()

#%% графики потоковых

import matplotlib.pyplot as plt

plt.grid(True)
plt.plot(x_node, U_flux[:, 0], label='rho_flux')
plt.plot(x_node, U_flux[:, 1], label='u_flux')
plt.plot(x_node, U_flux[:, 2], label='P_flux')
plt.plot(x_node, U_flux[:, 3], label='E_flux')
plt.legend()

