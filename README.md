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
   0 & \cdot  & \cdot & \cdot & \cdot & \cdot & 1 & 1
   \end{bmatrix}
}_{=:\boldsymbol{dTdt}}
\frac{1}{\Delta x²}
```

