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
%  Code: Solution of Advection-Diffusion-Reaction equations in 1D using   %
%        FVM discretization. The reaction terms are treated implicitly    %
%        and they refer to the reactions: A->B->C                         %
%                                                                         %
%        RA = -k1*CA^2;                                                   %
%        RB =  k1*CA^2 - k2*CB;                                           %
%        RC =  k2*CB;                                                     %
%                                                                         %
%-------------------------------------------------------------------------%

clc; clear; close all;

%-------------------------------------------------------------------------%
% Problem Data
%-------------------------------------------------------------------------%

L = 2;
tau = 5;
ncells = 100;
u = 1;
Dmix = 1.e-2;
k1 = 10;
k2 = 1;

CAin = 1.;
CBin = 0.;
CCin = 0.;

nsteps = 10000;

%-------------------------------------------------------------------------%
% Pre-Processing
%-------------------------------------------------------------------------%

% Build Mesh
h = L/ncells;
x = linspace(0, L, ncells+1);

% Choose dt for stability
dt = tau/nsteps;

%-------------------------------------------------------------------------%
% Memory Allocations
%-------------------------------------------------------------------------%

ZERO = zeros(1, ncells+2);

CA = ZERO; CB = ZERO; CC = ZERO;
RA = ZERO; RB = ZERO; RC = ZERO;
JA = ZERO; JB = ZERO; JC = ZERO;

CAo = CA; CBo = CB; CCo = CC;

SAimp = ZERO; SBimp = ZERO; SCimp = ZERO;
SAexp = ZERO; SBexp = ZERO; SCexp = ZERO;

%-------------------------------------------------------------------------%
% Solution loop
%-------------------------------------------------------------------------%

t = 0.;
for is=1:nsteps

    %---------------------------------------------------------------------%
    % Set the Boundary Conditions
    %---------------------------------------------------------------------%
    
    % Inlet Section
    CA(1) = 2*CAin - CA(2);
    CB(1) = 2*CBin - CB(2);
    CC(1) = 2*CCin - CC(2);

    % Outlet Section
    CA(ncells+2) = CA(ncells+1);
    CB(ncells+2) = CB(ncells+1);
    CC(ncells+2) = CC(ncells+1);

    %---------------------------------------------------------------------%
    % Compute Reaction Terms
    %---------------------------------------------------------------------%
    
    for i=2:ncells+1
        % Reaction rates at time t
        RA(i) = -k1*CA(i);
        RB(i) =  k1*CA(i) - k2*CB(i);
        RC(i) =  k2*CB(i);

        % Jacobian at time t
        JA(i) = -2*k1*CA(i);
        JB(i) =  2*k1*CA(i) - k2;
        JC(i) =  k2;

        % Explicit source terms
        SAexp(i) = RA(i) - JA(i)*CA(i);
        SBexp(i) = RB(i) - JB(i)*CB(i);
        SCexp(i) = RC(i) - JC(i)*CC(i);

        % Implicit source terms
        SAimp(i) = -JA(i);
        SBimp(i) = -JB(i);
        SCimp(i) = -JC(i);
    end

    %---------------------------------------------------------------------%
    % Solve Advection-Diffusion-Reaction Equations
    %---------------------------------------------------------------------%

    [MA,bA] = AdvectionDiffusionMatrix1D ...
        (CA, CAin, SAexp, SAimp, u, Dmix, dt, h, ncells);
    [MB,bB] = AdvectionDiffusionMatrix1D ...
        (CB, CBin, SBexp, SBimp, u, Dmix, dt, h, ncells);
    [MC,bC] = AdvectionDiffusionMatrix1D ...
        (CC, CCin, SCexp, SCimp, u, Dmix, dt, h, ncells);

    CA = MA\bA;
    CB = MB\bB;
    CC = MC\bC;

    %---------------------------------------------------------------------%
    % Post-Processing
    %---------------------------------------------------------------------%
    
    CAp = CellToFaceInterpolation (CA, ncells);
    CBp = CellToFaceInterpolation (CB, ncells);
    CCp = CellToFaceInterpolation (CC, ncells);

    if (mod(is,20)==1)
        fprintf ("Iter: %d - Time: %f\n", is, t);
        hold off;
        plot (x, CAp, "LineWidth", 1.8); hold on;
        plot (x, CBp, "LineWidth", 1.8);
        plot (x, CCp, "LineWidth", 1.8);
        xlabel ("lenght [m]"); ylabel ("Concentration [kmol/m3]");
        legend ("CA", "CB", "CC");
        xlim([0 L]); ylim([0 1]);
        drawnow;
    end

    % Advance the simulation time
    t = t + dt;
end

%-------------------------------------------------------------------------%
% Useful functions
%-------------------------------------------------------------------------%

function [A,b] = AdvectionDiffusionMatrix1D (C, Cin, expS, impS, u, Dmix, dt, h, ncells)

    % Build diagonals of the global tridiagonal matrix
    % Ap = central diagonal
    % Ae = upper diagonal
    % Aw = lower diagonal
    delta = u*dt/2/h;
    gamma = Dmix*dt/h^2;

    Ap =  1 + 2*gamma;
    Ae =  delta - gamma;
    Aw = -delta - gamma;

    n = ncells + 2;
    A = sparse(n,n);

    A(1,1) = 1; A(1,2) = 1;
    for i=2:n-1, A(i,i-1) = Aw;              end
    for i=2:n-1, A(i,i)   = Ap + impS(i)*dt; end
    for i=2:n-1, A(i,i+1) = Ae;              end
    A(n,n) = 1;  A(n,n-1) = -1;

    b = zeros(n);
    for i=2:ncells+1
        b(i) = C(i) + expS(i)*dt;
    end
    b(1) = 2*Cin;
    b(ncells+2) = 0.;
end

function Cface = CellToFaceInterpolation (Ccell, ncells)

    Cface = zeros(1,ncells+1);
    for i=1:ncells+1
        Cface(i) = 0.5*(Ccell(i+1) + Ccell(i));
    end

end