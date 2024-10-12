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
%  Code: 1D advection-diffusion equation by the FD method                 %
%        solution using explicit in time discretization method            %
%                                                                         %
% ----------------------------------------------------------------------- %

clc; close all; clear;

% ----------------------------------------------------------------------- %
% Pre-Processing
% ----------------------------------------------------------------------- %

global u D;

% Data
L = 2.0;        % domain length [m]
u = 1;          % velocity [m/s]
D = 0.05;       % diffusion coefficient [m2/s]
tau = 2*pi;     % total integration time [s] 

% Analytical/Exact solution
function f = exact (x, t)
    global u D;
    
    A = 0.5;        % amplitude of initial solution
    k = 1;          % wave number [1/m]
    
    f = A*exp(-4*pi^2*k^2*D*t)*sin(2*pi*k*(x - u*t));
end

% Build mesh
npoints = 51;               % number of grid points    
h = L/(npoints-1);          % grid step [m]
x = linspace(0,L,npoints);  % point coordinates (for post-processing)

% Build Time
dt_conv = 2*D/u^2;                  % Maximum time step for convection
dt_diff = h^2/2/D;                  % Maximum time step for diffusion
dt_max = min (dt_conv, dt_diff);    % Maximum allowed time step
sigma = 0.9;                        % Regulates dt with respect to dt_max
dt = sigma*dt_max;                  % Final simulation time step
nsteps = tau/dt;

fprintf('dt_conv=%g, dt_diff=%g, dt=%g\n', dt_conv, dt_diff, dt);

% Memory allocation
fo = zeros(npoints,1);      % temporary numerical solution
f  = zeros(npoints,1);      % current numerical solution
a  = zeros(npoints,1);      % exact solution

for i=1:npoints
    f(i) = exact (x(i), 0.);
end

% ----------------------------------------------------------------------- %
% Solution loop
% ----------------------------------------------------------------------- %

t = 0.;
for m=1:nsteps
    
    % New time step
    t = t+dt;       

    % Update analytical solution
    for i=1:npoints
        a(i) = exact (x(i), t);
    end
    
    % Forward Euler method
    fo=f;
    for i=2:npoints-1
		f(i) = fo(i)-(u*dt/2/h)*(fo(i+1)-fo(i-1))+...      % advection CDS
			   D*(dt/h^2)*(fo(i+1)-2*fo(i)+fo(i-1));       % diffusion
    end 
    
    % Periodic boundary condition
    f(npoints) = fo(npoints)-(u*dt/2/h)*(fo(2)-fo(npoints-1))+...   % advection CDS
            D*(dt/h^2)*(fo(2)-2*fo(npoints)+fo(npoints-1)); 
    f(1)  = f(npoints);

    if (mod(m,5) == 1)
        % Graphical output
        hold off; plot(0:h:L,f,'linewidth',2); axis([0 L -1, 1]);
        hold on; plot(0:h:L,a,'r','linewidth',2);
        hold on; legend('numerical', 'exact');
        xlabel('spatial coordinate [m]');
        ylabel('f');
        title('Evolution of the solution in space and time');
        grid on;
        drawnow;
    end
end

