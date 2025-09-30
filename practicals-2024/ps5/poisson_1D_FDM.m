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
%  Code: 1D Poisson equation using different methods to accelerate        %
%        the convergence: Jacobi, Gauss-Seidel, SOR. Using FDM.           %
%                                                                         %
% ----------------------------------------------------------------------- %

clc; close all; clear;

% Data
L = 1;          % Length of  the domain
fleft = 0;      % Boundary condition f(x=0) = 0
fright = 1;     % Boundary Condition f(x=1) = 1
S = 1;          % source term of the Poisson equation: laplacian(f) = S = 0

% Parameters for the iterative solution
maxiter = 10000;        % Maximum number of iterations to run
tolerance = 1e-3;       % Tolerance used to decide if the iterative solution converged
beta = 1.9;             % Over Relaxation parameter (for SOR method)

% Mesh
nx = 50;                % Number of points that discretize the 1D domain
h = L/(nx-1);           % Distance between two consecutive points
x = linspace(0,L,nx);   % Vector with coordinates of the points used to discretized x

% Create fields
f = zeros(1,nx);

% Set Boundary Conditions
% Since the FDM is used, the boundary conditions are directly imposed on
% the function f at the point belonging to the boundaries of the domain1
f(1) = fleft;
f(end) = fright;

% Iterative solution loop
% Loop from iteration 1 to the user-provided maximum number of iterations
for iter=1:maxiter
    
    % Store f at the previous iteration (just for Jacobi Method)
    fo = f;

    % Loop over all the internal points
    for i=2:nx-1
        f(i) = 0.5*(fo(i+1)+fo(i-1)-S*h^2);                     % Jacobi Method
        %f(i) = 0.5*(f(i+1)+f(i-1)-S*h^2);                       % Gauss-Seidel
        %f(i) = 0.5*beta*(f(i+1)+f(i-1)-S*h^2) + (1-beta)*f(i);  % SOR
    end

    % Compute residual:
    % Residual = |laplacian(f) - S|
    %          = | (f(i+1) + f(i-1) -2*f(i))/h^2 - S | <= tolerance
    res = 0;
    for i=2:nx-1
        res = res + abs( (f(i+1)+f(i-1)-2*f(i))/h^2 - S );
    end

    % Normalize residual with the number of points
    % in order to find the residual in the single point
    res = res/(nx-2);

    % If the residual is smaller than the user-defined tolerance
    if (res <= tolerance)
        % Leave the for loop
        break;
    end

    % If the maximum number of iteration was reached print this message
    % which tells you that the system did not converge
    if (iter==maxiter-1)
        fprintf("Maximum number of iterations reached\n");
    end
end

% Plot the solution: function f in every point x
plot(x, f, "LineWidth", 2);
xlabel("length [m]");
ylabel("scalar field f");
grid on;

% Print the number of iterations that were used to convergence
fprintf("Number of iterations = %d\n", iter);

