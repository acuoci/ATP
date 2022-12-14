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
%  Code: Population Balance Equation (PBE) with growth term only          %
%        Solved using the Discrete Sectional Method                       %
%        R. McGraw (1997) Description of Aerosol Dynamics by the          %
%        Quadrature Method of Moments, Aerosol Science and Technology,    %
%        27:2, 255-265 (1997), DOI: 10.1080/02786829708965471             %
%                                                                         %
% ----------------------------------------------------------------------- %

close all;
clear variables;

% Initial distribution: f=a*r^2*exp(-b*r)
% r = particle radius (mum)
% a and b: distribution parameters
a = 0.108;  % (1/mum/cm3)
b = 0.60;   % (1/mum)

% Growth rate: phi(r) = kG/r
kG = 0.78;  % growth rate constant (mum2/s)

% Domain of integration
rMax = 30;  % maximum radius (mum)
M = 500;    % number of intervals
tf = 20;    % maximum time (s)

% Initial distribution (#/cm3/mum)
r = 0:rMax/M:rMax;
fIn = fInitial(r,a,b);

% Number of particles (per unit of volume) in each interval (#/cm3)
NIn = zeros(M,1);
for i=1:M
    NIn(i) = fIn(i+1)*(r(i+1)-r(i));
end

[t, N] = ode15s(@ODESystem, 0:0.25:tf, NIn, [], r, kG);
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
     legend('numerical');
     frame = getframe(gcf);
end

% ODE System of Equations 
function dN = ODESystem(~, N, r, kG)

    M = length(N);

    f = zeros(M+1,1);
    for i=1:M
        f(i+1) = N(i)/(r(i+1)-r(i));
    end

    dN = zeros(M,1);
    dN(1) = -f(2)*Ldot(kG, r(2));
    for i=2:M-1
        dN(i) = -f(i+1)*Ldot(kG, r(i+1)) + f(i)*Ldot(kG, r(i));
    end
    dN(M) = f(M)*Ldot(kG, r(M));
end

% Growth rate function
function Lprime = Ldot(kG, r)
    Lprime = kG/r;
end

% Initial solution
function f = fInitial(r, a, b)
    f = a*(r.^2).*exp(-b*r);
end
