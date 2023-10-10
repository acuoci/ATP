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
%  Code: 1D advection equation for a discontinuous field using FVM        %
%        testing different advection schemes and flux limiters.           %
%                                                                         %
% ----------------------------------------------------------------------- %

clc; close all; clear;

%-------------------------------------------------------------------------%
% User-defined data
%-------------------------------------------------------------------------%

ncells = 200;         % number of cells
np = ncells + 1;     % number of grid points
dt = 0.001;           % time step [s]
nstep = 500;         % number of time steps
L = 2.0;             % domain length [m]
u = 1;               % velocity [m/s]
D = 0.;            % [Discontinuity] diffusion coefficient [m2/s]
A = 0.5;             % amplitude of initial solution
k = 1;               % wave number [1/m]

%-- [Discontinuity]
fleft = 2.;
fright = 1.;

%-------------------------------------------------------------------------%
% Pre-processing of user-defined data
%-------------------------------------------------------------------------%

% Grid step calculation
h = L/ncells;         % grid step [m]
x = linspace(0,L,np);

% Memory allocation
fo = zeros(ncells+2,1);      % temporary numerical solution
f  = zeros(ncells+2,1);      % current numerical solution
fp = zeros(ncells+1,1);      % solution interpolated on the nodes (for plot)
a  = zeros(ncells+2,1);      % exact solution
ap = zeros(ncells+1,1);

% Initial solution
for i=1:ncells+2
    xi = (h*(i-2) + 0.5*h);
    if (xi > 0.5*L)
        f(i) = 1;
    else
        f(i) = 2;
    end
end

% % Check the stability conditions on time step
Co = u*dt/h;                        % Courant number
dt_max = 1*h/u;
dt = 0.1*dt_max;

%-------------------------------------------------------------------------%
% Advancing in time
%-------------------------------------------------------------------------%

t = 0.;
for m=1:nstep

    % Interpolate f to nodes
    for i=1:np
        fp(i) = 0.5*(f(i+1) + f(i));
    end

    % Graphical output
    if (mod(m,100)==0 || m==1)
        plot(0:h:L,fp,'linewidth',2); axis([0 L 0, 3]); % plot num. 
        hold on;
        xlabel('spatial coordinate [m]');
        ylabel('solution');
        grid on;
        drawnow;
    end

    % Store old field
    fo = f;

    % [Discontinuity] Impose the west and east wall values using ghost cells
    f(1) = 2.*fleft - f(2);
    f(ncells+2) = 2.*fright - f(ncells+1);

    % Advance the solution in the internal cells except the last one
    for i=2:ncells+1 % [Discontinuity] (+1) to include last ghost cell
        % Ai = u/2/h*(fo(i+1) - fo(i-1));             % centered
        % Ai = u/h*(fo(i) - fo(i-1));                 % upwind

        % Advection term using flux limiters
        rE = r (fo, i);
        rW = r (fo, i-1);
        psiE = psi_superbee (rE);
        psiW = psi_superbee (rW);
        fE = fo(i) + 0.5*psiE*(fo(i) - fo(i-1));
        if (i-1 == 1)
            fW = fo(i);
        else
            fW = fo(i-1) + 0.5*psiW*(fo(i-1) - fo(i-2));
        end
        Ai = u/h*(fE - fW);

        Di = D/h^2*(fo(i+1) + fo(i-1) - 2*fo(i));
        f(i) = fo(i) + dt*(-Ai+Di);
    end

    % New time step
    t=t+dt;       
end

% Helper functions

function ri = r(f, i)
    if (i == 1)
        ri = 1;
    else
        ri = (f(i+1) - f(i))/(f(i) - f(i-1));
    end
end

% Flux limiters: see pag. 96 Ferziger-Peric for more

function psi = psi_minmod (r)
    psi = max (0, min (r, 1));
end

function psi = psi_vanleer (r)
    psi = (r + abs(r))/(1 + abs(r));
end

function psi = psi_muscl (r)
    psi = max (0., min ([2, 2*r, 0.5*(1+r)]));
end

function psi = psi_superbee (r)
    psi = max ([0, min(2*r, 1), min(r, 2)]);
end
