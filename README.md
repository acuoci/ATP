# ATP: Advanced Transport Phenomena
Collection of numerical codes presented in the Advanced Transport Phenomena (ATP) class at Politecnico di Milano.

Description
-----------
The numerical codes are organized in the following folders:

1. `lectures`: codes presented and discussed in the ATP lectures

2. `practicals`: codes and problems presented and discussed in the ATP practical sessions


## 1. Advection-diffusion-reaction equation in 1D with the Finite Difference (FD) method
The advection-diffusion equation is solved on a 1D domain using the finite-difference method. Constant, uniform velocity and diffusion coefficients are assumed. Spatial derivatives are discretized using 2nd-order, centered schemes. Time integration is carried out with different methods. Application to a metallic bar with fixed temperatures at the boundaries and possible heat exchange with the external environment.
* Time integration via ode45 solver: [bar_1d_ode45.m](lectures/FDM1D/bar_1d_ode45.m)
* Time integration via explicit (forward) Euler method: [bar_1d_explicit_euler.m](lectures/FDM1D/bar_1d_explicit_euler.m)
* Time integration via implicit (backward) Euler method: [bar_1d_implicit_euler.m](lectures/FDM1D/bar_1d_implicit_euler.m)

## 2. Advection-diffusion-reaction equation in 1D with the Finite Volume (FV) method
The advection-diffusion equation is solved on a 1D domain using the finite-volume method. Constant, uniform velocity and diffusion coefficients are assumed. Spatial derivatives are discretized using 2nd-order, centered schemes. Time integration is carried out with ode15s solver. Application to a metallic bar with fixed temperatures at the boundaries and possible heat exchange with the external environment.
* Time integration via ode15s solver: [bar_1d_ode15s.m](lectures/FVM1D/bar_1d_ode15s.m)
