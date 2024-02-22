# heatTransferSim
A 1D heat transfer solver written in python
# Background 
The goal is to solve the heat equation 
$$\frac{\partial T}{ \partial t} = \alpha  \nabla ^{2} T$$
we simlify and only consider the 1d case:
$$\frac{\partial T}{ \partial t} = \alpha  \frac{\partial ^2 T}{\partial x²} ^{2} $$
To solve the equation numerically we use the finite differences methode:
