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
%  Code: 1D advection diffusion equation using FVM and a time implicit    %
%        discretization of the time derivative.                           %
%                                                                         %
% ----------------------------------------------------------------------- %

clc; close all; clear;

%-------------------------------------------------------------------------%
% User-defined data
%-------------------------------------------------------------------------%

ncells = 500;           % number of cells
np = ncells + 1;        % number of grid points
L = 5;                  % domain length [m]
u = 1;                  % velocity [m/s]
D = 0.;                 % diffusion coefficient [m2/s]
tau = 1.5;              % total simulation time [s]
fleft = 2.;             % left side boundary value
fright = 1.;            % right side boundary value

%-------------------------------------------------------------------------%
% Pre-processing of user-defined data
%-------------------------------------------------------------------------%

% Grid step calculation
h = L/ncells;           % grid step [m]
x = linspace(0,L,np);

% Check the stability conditions on time step
dt_max = 1*h/u;
sigma = 2;
dt = sigma*dt_max;
nstep = tau/dt;

% Memory allocation
f  = zeros (ncells+2, 1);
fo = zeros (ncells+2, 1);
fp = zeros (ncells+1, 1);

% Initial solution
for i=1:ncells+2
    xi = (h*(i-2) + 0.5*h);
    if (xi > 0.5*L)
        f(i) = 1.;
    else
        f(i) = 2.;
    end
end

% Create Matrix A
n = ncells + 2;

% Known terms vector
b = zeros (ncells+2, 1);

% Build diagonals of the global tridiagonal matrix
alpha = u*dt/h;
beta  = D*dt/h^2;
Ap = -2*beta - alpha - 1;
Ae = beta;
Aw = beta + alpha;

% Create global matrix using sparse in order to avoid to store the "0"s
A = sparse (ncells+2, ncells+2);

% Fill the global matrix using the three diagonals
A(1,1)=1; A(1,2)=1;                         % From the boundary conditions
for i=2:ncells+1, A(i,i-1) = Aw; end
for i=2:ncells+1, A(i,i) = Ap;   end
for i=2:ncells+1, A(i,i+1) = Ae; end
A(n,n)=1; A(n,n-1)=1;                       % From the boundary conditions

% loop over all the time steps
for t=1:nstep
    
    % Store solution at the old time step
    fo = f;

    % Update RHS term using info from the old time step
    for i=2:ncells+1
        b(i) = -fo(i);
    end
    b(1) = 2*fleft;
    b(ncells+2) = 2*fright;

    % Solve the linear system of equations
    f = A\b;

    % Plot results
    if (mod(t,1) == 0 || t == 1)
        
        % Linear interpolation, in order to interpolate cell-centered
        % temperature values to the face-centered position.
        for i=1:ncells+1
            fp(i) = 0.5*(f(i+1) + f(i));
        end

        hold off;
        plot (x, fp, "LineWidth", 2);
        grid on;
        xlabel ("length [m]");
        ylabel ("Solution");
        drawnow;
    end
end
