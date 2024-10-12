# ATP: Practical Sessions

List of the practical sessions discussed during the Advanced Transport Phenomena class.

## [PS1](ps1): Discretization of the Time Derivative

Solution of the ODE system of equations:

$$
\begin{cases}
  \dfrac{d x_1}{dt} = -x_2 \\
  \dfrac{d x_2}{dt} = x_1
\end{cases}
$$

* **Analytical Solution:** in parametric form $x_1 = \cos(t)$, $x_2 = \sin(t)$.
* **Metodology:** 1<sup>st</sup> order forward Euler approach for the discretization of the time derivative.
* **Key Test:** comparison between the numerical and analytical solution, and evaluation of the convergence rate at time $t = 2\pi$.

<p align="middle" >
  <img src="doc/ps1-plots.png" width="40%" /> 
  <img src="doc/ps1-errors.png" width="40%" />
</p>

## [PS2](ps2): Solution of the Heat Equation

Solution of the 1D Heat Equation with source term, for the evolution of the temperature field $T$ along a metal slab.

$$
  \dfrac{\partial T}{\partial t} = \alpha \nabla^2 T + \dfrac{\dot{q}}{\rho Cp}
$$

where $\alpha = \lambda/\rho/Cp$ is the thermal diffusion coefficient, assuming that the density $\rho$, heat capacity $Cp$, and thermal conductivity $\lambda$ are constant. The equation is integrated for a sufficiently long time until reaching the steady state. We assume an initial temperature of the metal slab, and Dirichlet boundary conditions for the left (hot) and right (cold) boundaries.

* **Metodology:** 1<sup>st</sup> order forward Euler approach for the discretization of the time derivative; 2<sup>nd</sup> order centered approach for the discretization of the diffusion term. The system is approximated both using both the [Finite Difference](ps2/diffusioneq_1D_FDM_explicit.m) and the [Finite Volume](ps2/diffusioneq_1D_FVM_explicit.m) discretization methods.
* **Ket Test:** differences between FDM and FVM, including the implementation of the boundary conditions. Effect of the introduction of a source term on the steady state solution. Importance of the choice of the time step for the stability of the system.
