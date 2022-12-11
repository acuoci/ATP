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
tf = 10;    % maximum time (s)

% Initial distribution (#/cm3/mum)
r = 0:rMax/M:rMax;
fIn = fInitial(r,a,b);

% Number of particles (per unit of volume) in each interval (#/cm3)
NIn = zeros(M,1);
for i=1:M
    NIn(i) = fIn(i+1)*(r(i+1)-r(i));
end

% Continous phase variables
Y0 = 0;

% Setup initial conditions for ODE system
Y = Y0;
N0 = NIn;
N0(end+1) = Y0;
neq = length(N0);

[t, N] = ode15s(@ODESystem, 0:0.25:tf, N0, [], r, kG);
ntimes = length(t);

% Reconstruction of density function (#/cm3/mum)
f = zeros (ntimes, neq);
for i=1:ntimes
    for j=2:(length(N)-1)
        f(i+1,j) = N(i,j)/(r(i+1)-r(i));
    end
end

% Dynamic evolution of density function
figure;
for k=1:length(t)
     hold off;
     plot(r, f(k,:), 'b');
     hold on;
     xlabel('r (\mum)'); ylabel('f (#/micron/cm3)'); title('time=20 s'); 
     xlim([0 20]); ylim([0 0.6]);
     legend('numerical');
     frame = getframe(gcf);
end

% Dynamic of the continuous variable
figure;
plot (t,N(:,end));
xlabel('r (\mum)'); ylabel('Y'); title('time=20 s');

% Equations 
function dN = ODESystem(~, N, r, kG)

    neq = length(N);
    M = length(N) - 1;
    Y = N(end);

    f = zeros(M+1,1);
    for i=1:M
        f(i+1) = N(i)/(r(i+1)-r(i));
    end

    dN = zeros(neq,1);
    dN(1) = -f(2)*Ldot(kG, r(2), Y);
    for i=2:M-1
        dN(i) = -f(i+1)*Ldot(kG, r(i+1), Y) + f(i)*Ldot(kG, r(i), Y);
    end
    dN(M) = f(M)*Ldot(kG, r(M), Y);

    % Equation for the continuous variable
    sumNint = 0;
    for i=1:M
        sumNint = sumNint + N(i);
    end
    dN(M+1) = Ldot(kG, r(i+1), Y)*sumNint;
end

% Growth rate function
function Lprime = Ldot(kG, r, Y)
    % Lprime = kG/r;
    Lprime = kG*(5. - Y);
end

% Initial solution
function f = fInitial(r, a, b)

    f = a*(r.^2).*exp(-b*r);

end

