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
%  Code: 1D Projection method for solving the flow in a tube problem      %
%        using the FVM and explicit in time discretization on a           %
%        staggered grid including a drag force term in the momentum       %
%        equation in order to get a pressure drop.                        %
%                                                                         %
% ----------------------------------------------------------------------- %

clc; clear; close all;

% ----------------------------------------------------------------------- %
% Problem Data                                                            %
% ----------------------------------------------------------------------- %

% Problem Data
Lx = 50;
nx = 50;
hx  = Lx/nx;
Ain  = hx^2*pi/4;
Aout = hx^2*pi/4;
nu  = 1e-2;
tau = 20;

% Boundary conditions
uin = 1;
pout = 0;

% Parameters for SOR
maxiter = 1000;           % maximum number of iterations
beta = 1.9;               % SOR coefficient
tolerance = 1e-12;         % tolerance for convergence of Poisson equation

% ----------------------------------------------------------------------- %
% Pre-Processing                                                          %
% ----------------------------------------------------------------------- %

% Build mesh
x = 0:hx:Lx;         % Coordinates of the points/nodes along x direction (for post-processing) 
 
% [INOUT] Estimated max velocity considering also inlet and outlet velocity
% Needed in order to find the most limiting time step for the simulation
umax=max([abs(uin),abs(uin)*Ain/Aout]);

% Build Time
sigma = 0.1;                       % safety factor for time step (stability)
dt_diff = hx^2/4/nu;                 % maximum time step for stability of the diffusion
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
u = zeros(1,nx+1);
p = zeros(1,nx+2);

% Gamma coefficient is used in order to set the boundary conditions for the
% pressure equation. Using the staggered grid, the Poisson equation for
% pressure can be solved in every cell of the domain as usual just by changing
% the gamma coefficient depending on the position of the cell:
gamma = zeros(1,nx+2);      % Create the null field on the main grid

% 1. Internal cells
gamma(:) = 1/2;

% 2. Boundary cells (inlet and outlet section)
gamma(2) = 1;           % inlet treated like a boundary cell
gamma(nx+1) = 1/2;      % outlet treated like an internal cell

% 3. Ghost cells
gamma(1) = 0;
gamma(nx+2) = 0;

% Initial conditions: set reasonable initial velocity value instead of
% initializing everything to zero
u(:) = 1;
ut = u;

% Parameters for the drag force
Re = uin*hx/nu;

% ----------------------------------------------------------------------- %
% Solution over time
% ----------------------------------------------------------------------- %
t=0.0;
for is=1:nsteps

    % Boundary Conditions
    u(1) = uin;             % fixedValue, inlet boundary condition
    u(nx+1) = u(nx);        % zeroGradient, outlet boundary condition

    % ------------------------------------------------------------------- %
    % 1. Advection-diffusion equation (predictor)
    % ------------------------------------------------------------------- %

    for i=2:nx
        
        ue = 0.5*(u(i+1) + u(i));
        uw = 0.5*(u(i) + u(i-1));

        he = hx; hw = hx; dv = hx^2;

        Ai = ue^2*he - uw^2*hw;
        Di = nu*(u(i+1) - u(i)) - nu*(u(i) - u(i-1));

        f = 16/Re;                   % [DRAG]
        Fdrag = 0.5*u(i)^2*f*(2*hx); % [DRAG]

        ut(i) = u(i) + dt*(-Ai + Di)/dv - dt*Fdrag; %[DRAG]
    end

    % Boundary Conditions for temporary velocity
    ut(1) = uin;               % fixedValue, inlet boundary condition
    ut(nx+1) = u(nx+1);        % zeroGradient, outlet boundary condition

    % ------------------------------------------------------------------- %
    % 2. Pressure equation (Projection)
    % ------------------------------------------------------------------- %

    for iter=1:maxiter

        for i=2:nx+1
            delta = p(i+1) + p(i-1);
            S = hx/dt*(ut(i) - ut(i-1));
            p(i) = beta*gamma(i)*(delta - S) + (1-beta)*p(i);
        end

        res = 0.;
        for i=2:nx+1
            delta = p(i+1) + p(i-1);
            S = hx/dt*(ut(i) - ut(i-1));
            res = res + abs( p(i) - gamma(i)*(delta-S) );
        end
        res = res/nx;

        if (res <= tolerance)
            break;
        end

    end

    % ------------------------------------------------------------------- %
    % 3. Correction: Correct the velocity (only outside the obstacle)
    % ------------------------------------------------------------------- %

    for i=2:nx
        u(i) = ut(i) - dt*(p(i+1)-p(i))/hx;
    end

    % Print on the screen
    if (mod(is,20)==1)
        fprintf( 'Step: %d - Time: %f - Poisson iterations: %d\n', is, t, iter );
    end
    
    % Advance in time
    t=t+dt;

end

pplot = zeros(size(x));
p(1) = p(2);
for i=1:nx+1
    pplot(i) = 0.5*(p(i+1) + p(i));
end

% Plot velocity
figure;
plot(x,u,"LineWidth",1.8); xlabel("length [m]"); ylabel("Velocity [m/s]"); % plot velocity
ylim([0.95 1.05])

% Plot pressure
figure;
plot(x,pplot,"LineWidth",1.8); xlabel("length [m]"); ylabel("Pressure [Pa]");
