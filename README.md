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

