# ATP: Practical Sessions

List of the practical sessions discussed during the Advanced Transport Phenomena class.

## [Discretization of the Time Derivative](ps1)

Solution of the ODE system of equations:

$$
\begin{cases}
  \dfrac{d x_1}{dt} = -x_2 \\
  \dfrac{d x_2}{dt} = x_1
\end{cases}
$$

* **Analytical Solution:** in parametric form $x_1 = \cos(t)$, $x_2 = \sin(t)$.
* **Methodology:** 1<sup>st</sup> order forward Euler approach for the discretization of the time derivative.
* **Key Tests:** comparison between the numerical and analytical solution, and evaluation of the convergence rate at time $t = 2\pi$.

<p align="middle" >
  <img src="doc/ps1-plots.png" width="40%" /> 
  <img src="doc/ps1-errors.png" width="40%" />
</p>

## [Solution of the Heat Equation](ps2)

Solution of the 1D Heat Equation with source term, for the evolution of the temperature field $T$ along a metal slab.

$$
  \dfrac{\partial T}{\partial t} = \alpha \nabla^2 T + \dfrac{\dot{q}}{\rho Cp}
$$

where $\alpha = \lambda/\rho/Cp$ is the thermal diffusion coefficient, assuming that the density $\rho$, heat capacity $Cp$, and thermal conductivity $\lambda$ are constant. The equation is integrated for a sufficiently long time until reaching the steady state. We assume an initial temperature of the metal slab, and Dirichlet boundary conditions for the left (hot) and right (cold) boundaries.

* **Methodology:** 1<sup>st</sup> order forward Euler approach for the discretization of the time derivative; 2<sup>nd</sup> order centered approach for the discretization of the diffusion term. The system is approximated both using both the [Finite Difference](ps2/diffusioneq_1D_FDM_explicit.m) and the [Finite Volume](ps2/diffusioneq_1D_FVM_explicit.m) discretization methods.
* **Key Tests:** differences between FDM and FVM, including the implementation of the boundary conditions. Effect of the introduction of a source term on the steady state solution. Importance of the choice of the time step for the stability of the system.

<p align="middle" >
  <img src="doc/ps2-heatequation.gif" width="49%" /> 
  <img src="doc/ps2-heatequation-sources.gif" width="49%" />
</p>

## [Solution of a Benchmark Advection-Diffusion Equation](ps3)

Solution of the 1D Advection-Diffusion equation for the transport of a sine wave.

$$
  \dfrac{\partial f}{\partial t} + \mathbf{u}\cdot\nabla f = \Gamma \nabla^2 f
$$

where the velocity field $\mathbf{u}$ and the diffusion coefficient $\Gamma$ are constant and uniform. The problem is initialized to the analytical solution, and periodic boundary conditions are used.

* **Analytical Solution:** $f(x,t) = A \sin(2\pi k (x - \mathbf{u}t)) e^{-4\pi^2k^2\Gamma t}$, where $A$ is the amplitude, while $k$ is the wave number.
* **Methodology:** 1<sup>st</sup> order forward Euler approach for the discretization of the time derivative. Comparison of 1<sup>st</sup> order and 2<sup>nd</sup> order approaches for the discretization of the convective term: BDS and CDS approximations for the finite difference case, Upwind and Centered approximations for the finite volume case. 2<sup>nd</sup> order centered approach for the discretization of the diffusion term.
* **Key Tests:** implementation of the problem using FDM and FVM, including the implementation of periodic boundary conditions. Testing the stability of the time discretization, and the behavior of the system at different velocity and diffusivity, with particular focus on pure convective conditions ($\Gamma = 0$).

<p align="middle" >
  <img src="doc/ps3-plots.gif" width="49%" /> 
</p>

## [Transport of a Discontinuity](ps4/advection_discontinuity_1D_FVM.m)

Solution of the 1D Advection equation:

$$
  \dfrac{\partial f}{\partial t} + \mathbf{u}\cdot\nabla f = 0
$$

for a step function $f$. In these conditions, the centered scheme leads to strong oscillations, while the upwind scheme leads to unphysical diffusion. We implement a possible solution to this problem.

* **Methodology:** 1<sup>st</sup> order forward Euler approach for the discretization of the time derivative. Comparison of 1<sup>st</sup> order upwind and 2<sup>nd</sup> order centered approaches for the discretization of the convective terms, in finite volumes. Implementation of flux limiters to switch between upwind and centered schemes depending on the local gradient of the function.
* **Key Tests:** comparison between centered, upwind, and flux limiter cases.

<p align="middle" >
  <img src="doc/ps4-discontinuity.gif" width="49%" /> 
</p>

## [Time implicit solution of Advection-Diffusion equations](ps4/implicit_advection_diffusion_1D_FVM.m)

Solution of the 1D Advection equation:

$$
  \dfrac{\partial f}{\partial t} + \mathbf{u}\cdot\nabla f = \Gamma\nabla^2 f
$$

intialized as a Gaussian function:

$$
  f(x) = \dfrac{1}{\sigma \sqrt{2\pi}} \exp\left(-\dfrac{(x - \mu)^2}{2\sigma^2}\right)
$$

using a Dirichlet boundary condition for the right side of the domain $f(x=0,t) = 0$ and a Neumann boundary condition for the right side $\nabla f|_{x=L} = 0$.

* **Methodology:** 1<sup>st</sup> order backward Euler approach for the discretization of the time derivative. 1<sup>st</sup> order upwind approach for the discretization of the convective terms, 2<sup>nd</sup> order centered approach for diffusion term, using the FVM.
* **Key Tests:** verify that the solution is stable regardless the time step, and analyze the different solutions at varying Peclet number.

<p align="middle" >
  <img src="doc/ps4-implicit.gif" width="49%" /> 
</p>

## [Iterative solution of Poisson Equations](ps5)

Solution of the 1D Poisson equation:

$$
  \nabla\cdot\left(\nabla f\right) = S
$$

where the the source term is constant and equal to 1, the left and right boundaries of the domain are prescribed with Dirichlet boundary conditions ($f_{left} = 0$, and $f_{right} = 1$).

* **Methodology:** Jacobi method, Gauss-Seidel, and Successive Over-Relaxation for the iterative solution of the problem. Discretization using the Finite Difference Method, using a centered 2<sup>nd</sup> order approach for the Laplacian operator.
* **Key Tests:** observe the different number of iterations required to reach convergence using Jacobi, Gauss-Seidel and the SOR methods (at different over-relaxation parameters $\beta$).

<div align="center">

| Method  | $\boldsymbol{\beta}$ | iterations |
| ------- | ------- |----------- |
| Jacobi  | -  | 4043 |
| Gauss-Seidel |  - | 1982 |
| Gauss-Seidel + SOR | $\beta = 1.5$  | 667 |
| Gauss-Seidel + SOR | $\beta = 1.7$  | 355 |
| Gauss-Seidel + SOR | $\beta = 1.9$  | 128 |
| Gauss-Seidel + SOR | $\beta = 1.95$ | 244 |

</div>
