# ATP: Advanced Transport Phenomena
Collection of numerical codes presented in the Advanced Transport Phenomena (ATP) class at Politecnico di Milano.

Description
-----------
The numerical codes are organized in the following folders:

1. `lectures`: codes presented and discussed in the ATP lectures

2. `practicals`: codes and problems presented and discussed in the ATP practical sessions


## 1. Advection-diffusion-reaction equation in 1D with the Finite Difference (FD) method
The advection-diffusion-reaction equation is solved on a 1D domain using the finite-difference method. Constant, uniform velocity and diffusion coefficients are assumed. Spatial derivatives are discretized using 2nd-order, centered schemes. Time integration is carried out with different methods. Application to a metallic bar with fixed temperatures at the boundaries and possible heat exchange with the external environment.
* Time integration via ode45 solver: [bar_1d_ode45.m](lectures/FDM1D/bar_1d_ode45.m)
* Time integration via explicit (forward) Euler method: [bar_1d_explicit_euler.m](lectures/FDM1D/bar_1d_explicit_euler.m)
* Time integration via implicit (backward) Euler method: [bar_1d_implicit_euler.m](lectures/FDM1D/bar_1d_implicit_euler.m)

## 2. Advection-diffusion equation in 1D with the Finite Difference (FD) method
The advection-diffusion equation is solved on a 1D domain using the finite-difference method. Constant, uniform velocity and diffusion coefficients are assumed. The forward (or explicit) Euler method is adopted for the time discretization, while spatial derivatives are discretized using 2nd-order, centered schemes.
* Matlab script: [advection_diffusion_1d.m](lectures/FDM1D//advection_diffusion_1d.m)
* Matlab live script: [advection_diffusion_1d_live.mlx](lectures/FDM1D/advection_diffusion_1d_live.mlx)

## 3. Advection-diffusion-reaction equation in 1D with the Finite Volume (FV) method
The advection-diffusion-reaction equation is solved on a 1D domain using the finite-volume method. Constant, uniform velocity and diffusion coefficients are assumed. Spatial derivatives are discretized using 2nd-order, centered schemes. Time integration is carried out with ode15s solver. Application to a metallic bar with fixed temperatures at the boundaries and possible heat exchange with the external environment.
* Time integration via ode15s solver: [bar_1d_ode15s.m](lectures/FVM1D/bar_1d_ode15s.m)

## 4. Advection-diffusion-reaction equation in 2D with the Finite Difference (FD) method
The advection-diffusion-reaction equation is solved on a 2D rectangular domain using the finite-difference method. Analyically prescribed velocity fields are assumed. Constant and uniform diffusion coefficients are assumed. Spatial derivatives are discretized using 2nd-order, centered schemes.
* Evaporating pool: [evaporating_pool_2d_ode15s.m](lectures/FDM2D/evaporating_pool_2d_ode15s.m)
* Tubular reactor: [tubular_reactor_2d_ode15s.m](lectures/FDM2D/tubular_reactor_2d_ode15s.m)
* Tubular reactor (multiple reactions): [tubular_reactor_2d_multiple_reactions_ode15s.m](lectures/FDM2D/tubular_reactor_2d_multiple_reactions_ode15s.m)

## 5. Advection-diffusion equation in 2D with the Finite Difference (FD) method
The advection-diffusion equation is solved on a 2D rectangular domain using the finite-difference method. Constant, uniform velocity components and diffusion coefficients are assumed. The forward (or explicit) Euler method is adopted for the time discretization, while spatial derivatives are discretized using 2nd-order, centered schemes.
* Matlab script: [advection_diffusion_2d.m](lectures/FDM2D/advection_diffusion_2d.m)
* Matlab live script: [advection_diffusion_2d_live.mlx](lectures/FDM2D/advection_diffusion_2d_live.mlx)

## 6. Poisson equation in 2D
The Poisson equation is solved on a 2D rectangular domain using the finite-difference method. A constant source term is initially adopted. Spatial derivatives are discretized using 2nd-order, centered schemes. Different methods are adopted for solving the equation: the Jacobi method, the Gauss-Siedler method, and the Successive Over-Relaxation (SOR) method
* Matlab script: [poisson_2d.m](lectures/FDM2D/poisson_2d.m)
* Matlab live script: [poisson_2d_live.mlx](lectures/FDM2D/poisson_2d_live.mlx)

The same Poisson equation is solved by explicitly assembling the sparse matrix corresponding to the linear system arising after the spatial discretization
* Matlab script: [poisson_2d_matrix.m](lectures/FDM2D/poisson_2d_matrix.m)

## 7. Lorenz equations
The Lorenz system is a system of ordinary differential equations first studied by mathematician and meteorologist Edward Lorenz. It is notable for having chaotic solutions for certain parameter values and initial conditions. See also: [Lorenz system](https://en.wikipedia.org/wiki/Lorenz_system)
* Matlab script: [lorenz_equations.m](lectures/turbulence/lorenz_equations.m)

## 8. Population Balance Equations (PBE)
Numerical solution of a diffusion controlled growth 1D population balance equation using the Discrete Sectional Method (DSM) or the Quadrature Method of  Moments (QMOM).
* Matlab script (DSM): [dsm_diffusion_controlled_growth.m](lectures/PBE/dsm_diffusion_controlled_growth.m)
* Matlab script (QMOM): [qmom_diffusion_controlled_growth.m](lectures/PBE/qmom_diffusion_controlled_growth.m)
* Matlab live script (DSM): [dsm_diffusion_controlled_growth_live.m](lectures/PBE/dsm_diffusion_controlled_growth_live.m)
* Matlab live script (QMOM): [qmom_diffusion_controlled_growth_live.m](lectures/PBE/qmom_diffusion_controlled_growth_live.m)
