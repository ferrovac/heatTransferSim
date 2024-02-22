# heatTransferSim
A 1D heat transfer solver written in python
# Background 
The goal is to solve the heat equation 
$$\frac{\partial T}{ \partial t} = \alpha  \nabla ^{2} T$$
we simplify and only consider the 1D case:
$$\frac{\partial T}{ \partial t} = \alpha  \frac{\partial ^2 T}{\partial x²}$$
To solve the equation numerically we convert the PDE to system of ODE using the finite differences methode. The system of ODEs is then integrated using scypy.ivp_solve module.
# Discretisation
First we define the grid. We use $N$ equally spaced datapoints to represent the rod, where $T_{N}$ denotes the temperature at point $N$.
Using the FD methode we get the following approximation for the second spacial derivative:
$$\frac{\partial T}{\partial t} \approx \alpha \frac{T_{n+1} -2 T_n + T_{n-1}}{\Delta x²} $$
where $\Delta x$ is the spacial difference between two datapoints.
We can use this approximation to descripe the PDE as system of ODEs as matrix equation in the form:
```math
\frac{d}{dt} \underbrace{ \begin{bmatrix} T_0 \\ .\\.\\.\\ T_{N-1}\\T_{N} \end{bmatrix}}_{=:\vec{T}} = \alpha
\underbrace{
\begin{bmatrix} 
   0 & \cdot & \cdot & \cdot & \cdot & \cdot & 0 & 0\\
   1 & -2 & 1 & 0 & \cdot & \cdot & \cdot  & 0 \\
   0 & 1 & -2 & 1 & 0 & \cdot & \cdot & 0 \\
   \cdot & \cdot  & \cdot & \cdot & \cdot & \cdot &   & \cdot\\
   \cdot &   & \cdot & \cdot & \cdot & \cdot &  \cdot & \cdot\\
   \cdot &   &   & \cdot & \cdot & \cdot & \cdot  & 0\\
   0 & \cdot  & \cdot & \cdot & 0 & 1 & -2 & 1\\
   0 & \cdot  & \cdot & \cdot & \cdot & \cdot & 1 & -1
   \end{bmatrix}
}_{=:\boldsymbol{dTdt}}
\frac{1}{\Delta x²}
```
The first and last row of the matrix already contain a boundary condition. Not all boundary condition can be directly encoded in the matrix.
# Boundary Conditions
## Constant Temperature (Dirichlet)
The condition states:
$$T(x) = T_B$$
for $x = 0$ we get:
$$T_{0} = T_B$$
to incorporate the condition into the ode system we need to consider this in terms of $\frac{dT}{dt}$ i.e. $\frac{dT_{0}}{dt} = 0$
this leads to the first row of the matrix $\boldsymbol{dTdt}$ being all zeros.
## Heating (Neuman)
We can use Fourier's Law to model a heating point:
$$q = -k \nabla T$$
where $q$ is the local heat flux in $Wm²$, $k$ is the conductivity in $\frac{W}{mK}$.
note: 
$$q = \frac{P}{A}$$
where $P$ is the heating power and $A$ is the heating area.
We apply discretisation and consider the end point $x=L$:
$$q \approx -k \frac{T_{N+1}-T_{N}}{\Delta x}$$
note we used the forward methode. This leads to the point $T_{N+1}$ being present, which we do not know, because it is not part or our domain.
This is called a ghost point and we can get ridd of it with the heat equation itself. 
Since we need this boundary condition in terms of the time derivative as well we can solve the condition for $T_{N+1}$ and plug it into the heat equation. This also eliminates the ghost point:
$$T_{N+1} = T_N - \frac{q \Delta x}{k}$$
substituted into the heat equation:
```math
\left. \frac{\partial T(x,t)}{\partial t} \right|_{x=L} = \alpha \left.\frac{\partial² T(x,t)}{\partial x²} \right|_{x=L} \approx \alpha \frac{\left( T_N - \frac{q \Delta x}{k} \right) -2T_N + T_{N-1} }{\Delta x²}
```
this simplifies to
```math
\left. \frac{\partial T(x,t)}{\partial t} \right|_{x=L} \approx \alpha \frac{T_{N-1} -T_N -\frac{q \Delta x}{k}}{\Delta x²}
```
This boundary condition can not be directly implemented in the matrix. But we can implement it in the dTdt function we pass to scipy.ivp_solve:
```python
def dTdt(t, T):
    T_mod = np.copy(T)
    T_mod[0] = temp_left  # constant temp at left end
    ret = A.dot(T_mod)
    ret[N] = alpha/dx**2 * (1*(q*dx/k) - T_mod[N] + T_mod[N-1]) #heat source at right end 
    return ret
```
We can also use this to model thermal radiation losses because it's basically just a negative heat source that is proportional to temperature.
The Stefan-bolzmann law states:
$$q = \epsilon \sigma \left( T⁴ - T_{\text{env}}⁴ \right)$$
We can substitue this in the expression above and get radiation in our model for free.
# Model Validation
We can use analytic solutions for the steady state case to validate the model.
## Conduction
The rate of heat transfer through a homogenious rod is given by:
$$Q= -kA \frac{\Delta T}{L}$$
## Radiation
After we have validated Conduction we can use it to validate Radiation.
We plant a heat source at $n=N/2$ and the radiating boundary at $n=N$
In the steady state there will be a temperature difference between $T_{N/2}$ and $T_{N}$
using the conduction formula above we get the rate of heat transfer. This will have to match 
the Steffan-Bolzmann law:
$$Q = A \epsilon \sigma \left( T⁴ - T_{\text{env}}⁴ \right)$$
