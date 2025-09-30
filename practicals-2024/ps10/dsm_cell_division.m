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
%        breakage sink term which destroys particles with age             %
%        greater than tau0. When we destroy a particle we form            %
%        two new particles with age (tau) zero.                           %
%                                                                         %
% ----------------------------------------------------------------------- %

close all;
clear variables;

% Initial distribution: f=a*r^2*exp(-b*r)
% a and b: distribution parameters
a = 0.108;  % (1/mum/cm3)
b = 0.60;   % (1/mum)

% Domain of integration
tauMax = 120;  % maximum radius (mum)
M = 1000;      % number of intervals
tf = 200;      % maximum time (s)
tau0 = 70;     % age at which cells are broken
kG = 1.e-0;    % rate of breakage

% Initial distribution (#/cm3/mum)
tau = 0:tauMax/M:tauMax;
fIn = fInitial(tau,a,b);

% Number of particles (per unit of volume) in each interval (#/cm3)
NIn = zeros(M,1);
for i=1:M
    NIn(i) = fIn(i+1)*(tau(i+1)-tau(i));
end

[t, N] = ode15s(@ODESystem, 0:0.25:tf, NIn, [], tau, kG, tau0);
ntimes = length(t);

% Reconstruction of density function (#/cm3/mum)
f = zeros(ntimes, M+1);
for i=1:ntimes
    f(i, 2:end) = N(i,:)./(tau(2:end)-tau(1:end-1));
end

% Dynamic evolution of density function
figure;
for k=1:length(t)
     hold off;
     plot(tau, f(k,:), 'b');
     hold on;
     xlabel('r (\mum)'); ylabel('f (#/micron/cm3)');
     titlestring = strcat('time= ', num2str(t(k), '%.2f'), ' s'); title (titlestring);
     xlim([0 120]); ylim([0 0.6]);
     xline(tau0);
     legend('numerical');
     frame = getframe(gcf);
end

% Equations 
function dN = ODESystem(~, N, r, kG, r0)

    M = length(N);

    f = zeros(M+1,1);
    for i=1:M
        f(i+1) = N(i)/(r(i+1)-r(i));
    end

    fBc = 0.;
    for i=1:M
        rh = 0.5*(r(i+1) + r(i));
        fBc = fBc + 2*kG*heaviside(rh-r0)*N(i);
    end

    dN = zeros(M,1);
    dN(1) = -(f(2) - fBc);
    for i=2:M-1
        rh = 0.5*(r(i+1) + r(i));
        dN(i) = -f(i+1) + f(i) ...
            - kG*heaviside(rh-r0)*N(i);
    end

    rh = 0.5*(r(M+1) + r(M));
    dN(M) = f(M) - kG*heaviside(rh-r0)*N(M);
end

% Initial solution
function f = fInitial(r, a, b)
    f = a*(r.^2).*exp(-b*r);
end

