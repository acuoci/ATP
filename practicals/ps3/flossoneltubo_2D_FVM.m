% ----------------------------------------------------------------------- %
%   __  __       _______ _               ____  _  _       _______ _____   %
%  |  \/  |   /\|__   __| |        /\   |  _ \| || |   /\|__   __|  __ \  %
%  | \  / |  /  \  | |  | |       /  \  | |_) | || |_ /  \  | |  | |__) | %
%  | |\/| | / /\ \ | |  | |      / /\ \ |  _ <|__   _/ /\ \ | |  |  ___/  %
%  | |  | |/ ____ \| |  | |____ / ____ \| |_) |  | |/ ____ \| |  | |      %
%  |_|  |_/_/    \_|_|  |______/_/    \_|____/   |_/_/    \_|_|  |_|      %
%                                                                         %
% ----------------------------------------------------------------------- %
%                                                                         %
%   Authors: Alberto Cuoci <alberto.cuoci@polimi.it>                      %
%            Edoardo Cipriano <edoardo.cipriano@polimi.it>                %
%   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            %
%   Department of Chemistry, Materials and Chemical Engineering           %
%   Politecnico di Milano                                                 %
%   P.zza Leonardo da Vinci 32, 20133 Milano                              %
%                                                                         %
% ----------------------------------------------------------------------- %
%                                                                         %
%   This file is part of Matlab4ATP framework.                            %
%                                                                         %
%   License                                                               %
%                                                                         %
%   Copyright(C) 2022 Alberto Cuoci                                       %
%   Matlab4ATP is free software: you can redistribute it and/or           %
%   modify it under the terms of the GNU General Public License as        %
%   published by the Free Software Foundation, either version 3 of the    %
%   License, or (at your option) any later version.                       %
%                                                                         %
%   Matlab4CFDofRF is distributed in the hope that it will be useful,     %
%   but WITHOUT ANY WARRANTY; without even the implied warranty of        %
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         %
%   GNU General Public License for more details.                          %
%                                                                         %
%   You should have received a copy of the GNU General Public License     %
%   along with Matlab4ATP. If not, see <http://www.gnu.org/licenses/>.    %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
%  Code: 2D Projection method for solving the flow in a tube problem      %
%        using the FVM and explicit in time discretization on a           %
%        staggered grid.                                                  %
%                                                                         %
%-------------------------------------------------------------------------%

clc; clear; close all;

% ----------------------------------------------------------------------- %
% Problem Data                                                            %
% ----------------------------------------------------------------------- %

% Problem Data
Lx  = 5.4;          % length of the tube along x [m]
Ly  = 0.9;          % length of the tube along y [m]
h   = 0.05/2;       % cell size
nu  = 1e-2;         % Kinematic viscosity of the fluid [m2/s]
tau = 20;           % Total simulation time (long enough for steady state)

% Boundary conditions
un = 0; us = 0; ve = 0; vw = 0;     % Tangential wall velocities [m/s]
uin = 0.5;                          % Inlet velocity [m/s]

% Parameters for SOR
maxiter = 10000;          % maximum number of iterations
beta = 1.9;               % SOR coefficient
tolerance = 1e-6;         % tolerance for convergence of Poisson equation

% ----------------------------------------------------------------------- %
% Pre-Processing                                                          %
% ----------------------------------------------------------------------- %

% Build mesh
nx = Lx/h;          % Number of cells along x direction
ny = Ly/h;          % Number of cells along y direction
x = 0:h:Lx;         % Coordinates of the points/nodes along x direction (for post-processing) 
y = 0:h:Ly;         % Coordinates of the points/nodes along y direction (for post-processing)

% Check even number of cells
if (mod(nx,2)~=0 || mod(ny,2)~=0)
    error('Only even number of cells can be accepted (for graphical purposes only)');
end

% [INOUT] Inlet section (west side)
% Define cells belonging to the inlet section of the tube: y=2nd cell to y = nyth cell + 1 
nin_start = 2;      % first cell index 
nin_end = ny+1;     % last cell index

% [INOUT] Outlet section (east side)
% Define cells belonging to the outlet section of the tube: same as inlet
nout_start = 2;     % first cell index 
nout_end = ny+1;    % last cell index

% [INOUT] Inlet/Outlet section areas
% Define the area of the inlet and outlet section of the tube
% Area = (number of cells in the section of the tube) * height of a cell
Ain  = h*(nin_end-nin_start+1);      % inlet section area [m]
Aout = h*(nout_end-nout_start+1);    % outlet section area [m]

% [INOUT] Estimated max velocity considering also inlet and outlet velocity
% Needed in order to find the most limiting time step for the simulation
umax=max([abs(un),abs(us),abs(ve),abs(vw),uin,uin*Ain/Aout]);

% Build Time
sigma = 0.95;                       % safety factor for time step (stability)
dt_diff = h^2/4/nu;                 % maximum time step for stability of the diffusion
dt_conv = 4*nu/umax^2;              % maximum time step for stability of the convection
dt = sigma*min(dt_diff, dt_conv);   % time step of the simulation

% Print time step informations
fprintf("Maximum time step for convection = %f\n", dt_conv);
fprintf("Maximum time step for diffusion  = %f\n", dt_diff);
fprintf("Selected time step               = %f\n", dt);

% Number of time steps to run
nsteps = tau/dt;                      % number of steps

% ----------------------------------------------------------------------- %
% Memory allocation
% ----------------------------------------------------------------------- %

% Main fields (velocities and pressure)
u = zeros(nx+1,ny+2);   % u is defined on the uStaggered grid
v = zeros(nx+2,ny+1);   % v is defined on the vStaggered grid
p = zeros(nx+2,ny+2);   % p is defined on the main grid with dim [nx+2ghost cells]x[ny+2ghost cells]

% Gamma coefficient is used in order to set the boundary conditions for the
% pressure equation. Using the staggered grid, the Poisson equation for
% pressure can be solved in every cell of the domain as usual just by changing
% the gamma coefficient depending on the position of the cell:
%   if the cell is internal:            gamma = 1/4
%   if the cell is close to a boundary: gamma = 1/3
%   if the cell is close to an edge:    gamma = 1/2
%
% pij = gamma*beta*(pi+1,j + pi-1,j + pi,j+1 + pi,j-1 - S*h^2) + (1-beta)*pij
%
gamma = zeros(nx+2, ny+2);      % Create the null field on the main grid

% 1. Internal cells
gamma(3:nx,3:ny) = 1/4;         % Internal cells

% 2. Boundary cells
gamma(3:nx,2) = 1/3;            % South boundary
gamma(3:nx,ny+1) = 1/3;         % North boundary
gamma(nx+1,3:ny) = 1/3;         % East boundary
gamma(2,3:ny) = 1/3;            % West boundary

% 3. Edges cells
gamma(2,2) = 1/2;               % South-west edge cell
gamma(2,ny+1) = 1/2;            % North-west edge cell
gamma(nx+1,2) = 1/2;            % South-east edge cell
gamma(nx+1,ny+1) = 1/2;         % North-east edge cell 

% [INOUT] Correction of gamma coefficients for taking into account inlet and outlet sections
gamma(nx+1,nout_start:nout_end) = 1/4;      % Outlet section is treated as an internal cell
gamma(2,nin_start:nin_end) = 1/3;           % Inlet section is treated like a boundary

% [INOUT] Correction of gamma coeffiecient for edges
gamma(2,nin_start) = 1/2;   % Correct edges
gamma(2,nin_end) = 1/2;

% Initial conditions: set reasonable initial velocity value instead of
% initializing everything to zero
u(:,:) = 0.5;
ut = u; vt = v;

% ----------------------------------------------------------------------- %
% Solution over time
% ----------------------------------------------------------------------- %
t=0.0;
for is=1:nsteps
    
    % Boundary conditions (for tangential wall velocities, they're all 0)
    u(1:nx+1,1)    = 2*us - u(1:nx+1,2);                % south wall
    u(1:nx+1,ny+2) = 2*un - u(1:nx+1,ny+1);             % north wall
    v(1,1:ny+1)    = 2*vw - v(2,1:ny+1);                % west wall
    v(nx+2,1:ny+1) = 2*ve - v(nx+1,1:ny+1);             % east wall
    
    % [INOUT] Over-writing inlet conditions
    u(1,nin_start:nin_end) = uin;               % fixed inlet velocity
    
    % [INOUT] Over-writing outlet conditions
    u(nx+1,nout_start:nout_end) = u(nx,nout_start:nout_end);    % zero-gradient outlet velocity
    v(nx+2,nout_start:nout_end) = v(nx+1,nout_start:nout_end);  % zero-gradient outlet velocity
    
    % ------------------------------------------------------------------- %
    % 1. Advection-diffusion equation (predictor)
    % ------------------------------------------------------------------- %

    % u velocity
    for i=2:nx
        for j=2:ny+1

            % Advection
            ue = 0.25*(u(i+1,j) + u(i,j))^2;
            uw = 0.25*(u(i,j) + u(i-1,j))^2;
            un = 0.25*(u(i,j+1) + u(i,j))*(v(i+1,j) + v(i,j));
            us = 0.25*(u(i,j) + u(i,j-1))*(v(i+1,j-1) + v(i,j-1));

            Aij = (ue - uw + un - us)/h;
            
            % Diffusion
            Dij = nu/h^2*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4*u(i,j));

            ut(i,j) = u(i,j) + dt*(-Aij+Dij);
        end
    end

    % v velocity
    for i=2:nx+1
        for j=2:ny
            
            % Advection
            vn = 0.25*(v(i,j+1) + v(i,j))^2;
            vs = 0.25*(v(i,j) + v(i,j-1))^2;
            ve = 0.25*(u(i,j+1) + u(i,j))*(v(i+1,j) + v(i,j));
            vw = 0.25*(u(i-1,j+1) + u(i-1,j))*(v(i,j) + v(i-1,j));
            
            Aij = (ve - vw + vn - vs)/h;

            % Diffusion
            Dij = nu/h^2*(v(i+1,j) + v(i-1,j) + v(i,j+1) + v(i,j-1) - 4*v(i,j));

            vt(i,j) = v(i,j) + dt*(-Aij+Dij);
        end
    end

    % [INOUT] Update boundary conditions for temporary velocity
    ut(1,nin_start:nin_end) = u(1,nin_start:nin_end);            % fixed inlet velocity
    ut(nx+1,nout_start:nout_end) = u(nx+1,nout_start:nout_end);  % zero-gradient outlet velocity
    vt(nx+2,nout_start:nout_end) = v(nx+2,nout_start:nout_end);  % zero-gradient outlet velocity

    % ------------------------------------------------------------------- %
    % 2. Pressure equation (Projection)
    % ------------------------------------------------------------------- %
    for iter=1:maxiter
        
        for i=2:nx+1
            for j=2:ny+1
                delta = p(i+1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1);
                S = h/dt*(ut(i,j) - ut(i-1,j) + vt(i,j) - vt(i,j-1));
                p(i,j) = beta*gamma(i,j)*(delta - S) + (1-beta)*p(i,j);
            end
        end

        res = 0.;
        for i=2:nx+1
            for j=2:ny+1
                delta = p(i+1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1);
                S = h/dt*(ut(i,j) - ut(i-1,j) + vt(i,j) - vt(i,j-1));
                res = res + abs( p(i,j) - gamma(i,j)*(delta-S) );
            end
        end
        res = res/(nx*ny);

        if (res <= tolerance)
            break;
        end

    end
    
    % ------------------------------------------------------------------- %
    % 3. Correction: Correct the velocity (only outside the obstacle)
    % ------------------------------------------------------------------- %
    for i=2:nx
        for j=2:ny+1
            u(i,j) = ut(i,j) - dt*(p(i+1,j)-p(i,j))/h;
        end
    end
    for i=2:nx+1
        for j=2:ny
            v(i,j) = vt(i,j) - dt*(p(i,j+1)-p(i,j))/h;
        end
    end

    % [INOUT] Correct the velocity on the outlet patch
    u(nx+1,nout_start:nout_end) = ut(nx+1,nout_start:nout_end) - ...
        dt/h*(p(nx+2,nout_start:nout_end) - p(nx+1,nout_start:nout_end));

    % [INOUT] Because of numerical errors in the solution of the equations,
    %         the overall continuity equation, i.e. the conservation of
    %         mass cannot be guaranteed. It is better to correct the outlet
    %         velocity in order to force conservation of mass
    Qin = mean(u(1,nin_start:nin_end))*Ain;         % inlet volumetric flowrate
    Qout = mean(u(nx+1,nout_start:nout_end))*Aout;  % outlet volumetric flowrate  
    if (abs(Qout)>1.e-6)    % Only if some fluid is flowing outside the tube
        % Correct outlet velocity to ensure conservation of mass
        u(nx+1,nout_start:nout_end) = u(nx+1,nout_start:nout_end)*abs(Qin/Qout);
    end
    
    % Print on the screen
    if (mod(is,20)==1)
        fprintf( 'Step: %d - Time: %f - Poisson iterations: %d\n', is, t, iter );
    end
    
    % Advance in time
    t=t+dt;
    
end

% ----------------------------------------------------------------------- %
% Final Post-processing                                                   %
% ----------------------------------------------------------------------- %

% Interpolate on the nodes: values are stored at different locations of the
% grid. However, we want to plot every quantity in the same location which
% correspond to the set of coordinates x and y. To do this, we need to
% interpolate the computed quantities from cells to nodes using linear
% interpolation:
%    field_to_plot = node_interp(field, gridPosition, nx, ny);
uu = node_interp(u,"u",nx,ny);
vv = node_interp(v,"v",nx,ny);
pp = node_interp(p,"p",nx,ny);

% Surface map: u-velocity
subplot(311);
surf(x,y,uu'); view(2); axis tight; colormap jet; colorbar; shading interp;
title('u'); xlabel('x'); ylabel('y');

% Surface map: v-velocity
subplot(312);
surf(x,y,vv'); view(2); axis tight; colormap jet; colorbar; shading interp;
title('v'); xlabel('x'); ylabel('y');

% Surface map: pressure
subplot(313);
surf(x,y,pp'); view(2); axis tight; colormap jet; colorbar; shading interp;
title('p'); xlabel('x'); ylabel('y');

% ----------------------------------------------------------------------- %
% Functions                                                               %
% ----------------------------------------------------------------------- %
function [fnode] = node_interp(f, grid, nx, ny)

    fnode = zeros(nx+1,ny+1);

    if (grid == "u")
        for i=1:nx+1
            for j=1:ny+1
                fnode(i,j) = 0.50*(f(i,j) + f(i,j+1));
            end
        end

    elseif (grid == "v")
        for i=1:nx+1
            for j=1:ny+1
                fnode(i,j) = 0.50*(f(i,j) + f(i+1,j));
            end
        end

    elseif (grid == "p")
        for i=1:nx+1
            for j=1:ny+1
                fnode(i,j) = 0.25*(f(i,j) + f(i+1,j) + f(i,j+1) + f(i+1,j+1));
            end
        end
    end
end
