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
%  Code: 1D diffusion equation by the FV method using implicit            %
%        in time discretization method                                    %
%                                                                         %
% ----------------------------------------------------------------------- %

clc; close all; clear;

% ----------------------------------------------------------------------- %
% Pre-Processing
% ----------------------------------------------------------------------- %

% Data
alpha = 0.01;       % Diffusivity coefficient: alpha=k/rho/cp [m2/s]
L = 1;              % length of the domain [m]

% Build Mesh
ncells = 20;       % number of cells that discretize the 1D domain
h = L/ncells;      % distance between two consecutive points

% create a vector with the position of the points (face coordinates)
% the vector is composed of "npoints" evenly spaced from x=0 to x=L
x = linspace(0,L,ncells+1);

% Time
tau = 50;                 % total simulation time (high enough to reach steady state) [s]
dt_diff = 0.5*h^2/alpha;  % Maximum time step that accounts for the diffusion phenomena (from Di=0.5)
sigma = 3;                % Safety factor to avoid to work exactly at the minimum stability conditions
dt = sigma*dt_diff;       % Choice of the most limiting time step

% Print the computed minimum delta t. %f tells a "floating-point" number
% has to be printed. \n goes to the next line
fprintf("Maximum time step for diffusion  = %f\n", dt_diff);
fprintf("Selected time step               = %f\n", dt);

nsteps = tau/dt;    % Number of time steps to run

% Initial Conditions & Boundary Conditions
Tleft = 500;        % BC for temperature on the left side of the domain  [K]
Tright = 300;       % BC for temperature on the right side of the domain [K]
Tinit = 300;        % IC for temperature at time=0 [K].

% Memory allocations
T = ones(ncells+2)*Tinit; % Number of cells + 2 ghost cells
Tplot = zeros(size(x));

% ----------------------------------------------------------------------- %
% Solution loop
% ----------------------------------------------------------------------- %

% Re-set index for convenience
n = ncells + 2;

% Set known terms (known at time t)
b  = zeros(n,1);

% Build diagonals of the global tridiagonal matrix
% Ap = central diagonal
% Ae = upper diagonal
% Aw = lower diagonal
Ap = -(1 + 2*alpha*dt/h^2);
Ae = alpha*dt/h^2;
Aw = alpha*dt/h^2;

% Create global matrix using sparse in order to avoid to store the "0"s
A = sparse(n,n);

% Fill the global matrix using the three diagonals
A(1,1)=1; A(1,2)=1;                             % From the boundary conditions
for i=2:n-1, A(i,i-1) = Aw; end
for i=2:n-1, A(i,i)   = Ap; end
for i=2:n-1, A(i,i+1) = Ae; end
A(n,n)=1; A(n,n-1)=1;   % From the boundary conditions


% loop over all the time-steps
for t=1:nsteps

    To = T;  % Store temperature at time t

    % Update RHS term using info from the old time step
    for i=2:n-1
        b(i) = -To(i);
    end
    b(1) = 2*Tleft;
    b(n) = 2*Tright;

    % Solve the linear system of equations
    T = A\b;

    % On-The-Fly Post Processing
    if (mod(t,50)==1)   % => Every 50 time steps

        % Linear interpolation, in order to interpolate cell-centered
        % temperature values to the face-centered position.
        for i=1:ncells+1
            Tplot(i) = 0.5*(T(i+1) + T(i));
        end

        plot(x,Tplot, "LineWidth", 1.8);      % Plot results
        grid on;                          % Show a grid
        xlabel("length [m]")              % Name of the x axis
        ylabel("Temperature [K]")         % Name of the y axis
        drawnow;                          % Show the plot
    end

end

