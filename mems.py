import numpy as np
import matplotlib.pyplot as plt
from  scipy.integrate import solve_ivp
from scipy.interpolate import splrep, BSpline
from scipy.ndimage import uniform_filter1d

L = 1E-3
N = 300
temp_left = 70
temp_right = 1273

alpha = 7.18E-5
material = "AlN"
timeRes = 2000
p = 13
sigma = 5.68E-8
epsilon = 0.9
Cs = 760
rho = 3270
Ac = 250E-6 * 250E-6
q = p /Ac 
endTime = 0.06
k = alpha * rho * Cs 
dt = endTime/(timeRes-1)
dx = L / (N + 1)
x = np.linspace(dx, L-dx, N+1)



# Adjusted construction of A
A = np.zeros((N+1, N+1))
A[0, 0] = 0  # Dirichlet condition at x=0, ensuring it stays constant
#A[N,N] = 0
for i in range(1, N):
    A[i, i-1] = A[i, i+1] = alpha
    A[i, i] = -2 * alpha
A[N, N-1] = 1 #* alpha  # Neumann condition approximation
A[N, N] = -1 #* alpha
A /= dx**2  # Apply the spatial discretization scaling



def dTdt(t, T):
    T_mod = np.copy(T)
    T_mod[0] = temp_left  # Apply Dirichlet condition at x=0
    ret = A.dot(T_mod)
    point = int(N-1)
    if t < endTime/2:
        ret[point] = alpha/dx**2 * (1*(q*dx/k) - T_mod[point] + T_mod[point-1])#heat in the middle
    ret[-1] = alpha/dx**2 * (( (-1*sigma * epsilon * T_mod[-1]**4) * dx/k) - T_mod[-1] + T_mod[-2]) #radiate at the end
    
    return ret

# Initial condition adjustment
T_init = np.zeros_like(x)
#T_init[0] = temp_left  # Ensure the left boundary condition is applied
#T_init[-1] = temp_left
for i in range(0, len(T_init)):
   T_init[i] = 70# temp_left + ((temp_right - temp_left) / N * i)
#T_init[-1] = temp_left

t_span = (0,endTime)
sol = solve_ivp(dTdt, t_span, T_init, method='Radau', t_eval=np.linspace(*t_span, timeRes))
plt.plot(x,T_init)
plt.title("Starting Condition")
plt.xlabel("Distance along rod [m]")
plt.ylabel("Temperature [K]")
plt.show()

T = np.zeros_like(sol.t)
dT_dt = np.zeros_like(sol.t)


for i in range(len(sol.t)-1):
    T[i] = sol.y[-1,i]
    dT_dt[i] = (sol.y[-1,i] - sol.y[-1,i+1]) / dt * Cs * rho * Ac * dx


#plt.plot(T[1:-2], dT_dt[1:-2])
#plt.title("Counter Heat Power")
#plt.xlabel(r'$T$')
#plt.ylabel(r'$\frac{dT}{dt}$ [W]')
#plt.show()
#plt.figure()




for i in range(0, len(sol.t)):
    plt.plot(x, sol.y[:,i], 'black',label="TipTemp: " + str(np.round(sol.y[-1,i],0)) + "K" )
    plt.legend()
    if sol.t[i] < 0.07:
        plt.title("time: " + str(np.round(sol.t[i],4)) + "s   Power: " + str(p) + "W")
    else:
        plt.title("time: " + str(np.round(sol.t[i],4))+ "s   Power: " + str(0) + "W")
    #plt.ylim([np.min(T_init)-5,np.max(T_init)+5])
    plt.pause(.0001)
    plt.cla()
    plt.ylabel("Temperature (K)")
    plt.xlabel("Distance Along Rod (m)")
    

plt.figure()
for i in range(20):
    plt.plot(x,sol.y[:,int(timeRes/20*i)], label=str(np.round(sol.t[int(timeRes/20*i)],3))+"s  "+ str(np.round(sol.y[-1,int(timeRes/20*i)],1))+ "K")
plt.plot(x,sol.y[:,-1],label=str(np.round(sol.t[-1],3))+"s "+ str(np.round(sol.y[-1,-1]))+"K")
plt.legend()
plt.title(material)
plt.ylabel("Temperature (K)")
plt.xlabel("Distance Along Rod (m)")
plt.show()
