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
%   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       %
%           Edoardo Cipriano <edoardo.cipriano@polimi.it>                 %
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
%  Code: Forward Euler discretization of an ODE system                    %
%        of equations.                                                    %
%        The problem setup was inspired by                                %
%        https://antoonvanhooft.nl/min_ns/chapt2                          %
%                                                                         %
% ----------------------------------------------------------------------- %

clc; close all; clear;

% We create a vector with the number of time steps of the simulation
Nstepsvec(1) = 10;
while (Nstepsvec(end) < 10000.)
    Nstepsvec(end+1) = Nstepsvec(end)*2;
end

% We allocate a vector with dimensions equal to the number of time steps,
% and we will use this vector to store the error at each time step value.
errorsvec = zeros(size(Nstepsvec));

% We call the function run, at a specific time step Nstepsvec(i), and we
% save the error returned by the function.
for i=1:length(Nstepsvec)
    errorsvec(i) = run (Nstepsvec(i));
end

% The errors are plotted in logarithmic scale, for both x and y axis, and
% compared with the 1^st and 2^nd order trends. It's important to set hold
% on after the first 'loglog' command. Otherwise, the log scale is not
% applied.
figure;
loglog (Nstepsvec, errorsvec, 'x', LineWidth=2);
hold on;
loglog (Nstepsvec, 20*Nstepsvec.^-1, LineWidth=2);
loglog (Nstepsvec, 400*Nstepsvec.^-2, LineWidth=2);
xlabel ("N");
ylabel ("L_2");
legend ("Results", "1^{st} order", "2^{nd} Order");
title ("Convergence rate of the time derivative");
axis("square");
grid ("on");

% We write a function that solves a simple ODE system of equations using a
% segregated approach and a Forward Euler approximation of the time
% derivative:
%   dx1/dt = -x2
%   dx2/dt = x1
%   x1(t = 0) = 1.
%   x2(t = 0) = 0.
%
% The analytical solution reads:
%   x1(t) = cos(t)
%   x2(t) = sin(t)
%
function [error] = run(Nsteps)
    % We set the initial values of the variables x1 and x2 being solved
    x1 = 1.;            % Initial value of x1
    x2 = 0.;            % Initial value of x2
    tend = 2.*pi;       % Total simulation time
    dt = tend/Nsteps;   % Time step
    
    x1vec(1) = x1;
    x2vec(1) = x2;
    
    % We solve the system of ODE
    for i=1:Nsteps
        x1old = x1;
        x2old = x2;
    
        x1 = x1old - dt*x2old;  % Advance solution of x1 at the new time step
        x2 = x2old + dt*x1old;  % Advance solution of x2 at the new time step
    
        x1vec(i+1) = x1;
        x2vec(i+1) = x2;
    
        fprintf ("x1 = %g - x2 = %g\n", x1, x2);   % Print current solution
    end
    
    % The error is computed at the final time
    error1 = x1 - 1.;
    error2 = x2 - 0.;
    L2norm = (error1*error1 + error2*error2)^0.5;
    error = L2norm;
    
    % We plot the numerical and analytical solutions
    figure(1);
    t = linspace (0., 2.*pi, 200);
    plot (cos(t), sin(t), LineWidth=2);
    hold on;
    plot (x1vec, x2vec, 'x', LineWidth=2);
    xlim ([-1.5 1.5]);
    ylim ([-1.5 1.5]);
    title ("Comparison Between Theretical and Numerical Solution")
    legend ("Analytic", "Numerical")
    grid on;
    axis ("square");
end
