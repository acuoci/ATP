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
%  Code: Population Balance Equation (PBE) with growth term and           %
%        breakage sink and source terms. Particles with radius            %
%        greater than r0 are destroyed forming two particles with         %
%        half the radius of the destroyed particle.                       %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% Initial distribution: f=a*r^2*exp(-b*r)
% r = particle radius (mum)
% a and b: distribution parameters
a = 0.108;  % (1/mum/cm3)
b = 0.60;   % (1/mum)

% Domain of integration
% rMax = 30;  % maximum radius (mum)
% M = 500;    % number of intervals
rMax = 120;   % maximum radius (mum)
M = 1000;     % number of intervals
tf = 20;      % maximum time (s)
r0 = 12;      % radius at which particles are broken
kG = 1.e-1;   % rate of breakage

% Initial distribution (#/cm3/mum)
r = 0:rMax/M:rMax;
fIn = fInitial(r,a,b);

% Number of particles (per unit of volume) in each interval (#/cm3)
NIn = zeros(M,1);
for i=1:M
    NIn(i) = fIn(i+1)*(r(i+1)-r(i));
end

[t, N] = ode15s(@ODESystem, 0:0.25:tf, NIn, [], r, kG, r0);
ntimes = length(t);

% Reconstruction of density function (#/cm3/mum)
f = zeros(ntimes, M+1);
for i=1:ntimes
    f(i, 2:end) = N(i,:)./(r(2:end)-r(1:end-1));
end

% Dynamic evolution of density function
figure;
for k=1:length(t)
     hold off;
     plot(r, f(k,:), 'b');
     hold on;
     xlabel('r (\mum)'); ylabel('f (#/micron/cm3)');
     titlestring = strcat('time= ', num2str(t(k), '%.2f'), ' s'); title (titlestring);
     xlim([0 20]); ylim([0 0.6]);
     xline(r0);
     legend('numerical');
     frame = getframe(gcf);
end

% Equations 
function dN = ODESystem(~, N, r, kG, r0)

    M = length(N);
    Ldot = 0.68;

    f = zeros(M+1,1);
    for i=1:M
        f(i+1) = N(i)/(r(i+1)-r(i));
    end

    integ = 0;
    for i=1:M
        rh = 0.5*(r(i+1) + r(i));
        integ = integ + 2*kG*heaviside(rh-r0)*N(i);
    end

    dN = zeros(M,1);
    dN(1) = -Ldot*f(2) + integ*dirac(0.5*r0, r(2), r(1)) - kG*heaviside(rh-r0)*N(1);
    for i=2:M-1
        rh = 0.5*(r(i+1) + r(i));
        % dN(i) = 0 ...
        dN(i) = -Ldot*f(i+1) + Ldot*f(i) ...
            + integ*dirac(0.5*r0, r(i+1), r(i)) ...
            - kG*heaviside(rh-r0)*N(i);
    end
    dN(M) = Ldot*f(M) + integ*dirac(0.5*r0, r(M), r(M-1)) - kG*heaviside(rh-r0)*N(M);
end

function delta = dirac(r0, rr, rl)
    if (r0 >= rl && r0 <= rr)
        delta = 1;
    else
        delta = 0;
    end
end

% Initial solution
function f = fInitial(r, a, b)
    f = a*(r.^2).*exp(-b*r);
end

