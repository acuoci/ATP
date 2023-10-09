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
%        Solved using the QMOM (McGraw, 1997)                             %
%        R. McGraw (1997) Description of Aerosol Dynamics by the          %
%        Quadrature Method of Moments, Aerosol Science and Technology,    %
%        27:2, 255-265 (1997), DOI: 10.1080/02786829708965471             %
%                                                                         %
% ----------------------------------------------------------------------- %

close all;
clear variables;

% Initial distribution: f=a*r^2*exp(-b*r) <=> f=a*V^(2/3)*exp(-b*V^(1/3))
% r = particle radius (mum)
% V = particle volume (mum^3)
% a and b: distribution parameters
a = 0.108;  % (1/mum3/cm3)
b = 0.60;   % (1/mum)

% Growth rate: phi(r) = kG/r <=> phi(V) = kG/V^(1/3)
kG = 0.78;  % growth rate constant (mum2/s)

% Domain of integration
rMax = 50;  % maximum radius (mum)
M = 1000;   % number of intervals
tf = 20;    % maximum time (s)

% Initial distribution
Vmax = rMax^3;
r = 0:rMax/M:rMax;
V = r.^3;
fIn = a*V.^(2/3).*exp(-b*V.^(1/3));

NIn = fIn(1:end-1).*(V(2:end)-V(1:end-1));


% Solution of equation moments
[t, N] = ode15s(@ODESystem, 0:1:tf, NIn, [], V, kG);
ntimes = length(t);

f = zeros(ntimes, M+1);
for i=1:ntimes
    f(i, 1:end-1) = N(i,:)./(V(2:end)-V(1:end-1));
end

figure; hold on;
plot(V(1:end).^(1/3), fIn);
plot(V(1:end).^(1/3), f(end,:));
hold off;

% Calculation of moments
m0 = Moment(V, f(end,:), 0);


function m = Moment(r, f, order)

    np = length(r);
    m = 0;
    for i=1:np-1
        deltar = r(i+1)-r(i);
        F = 0.50*(f(i+1)*r(i+1)^order+f(i)*r(i)^order);
        m = m + deltar*F;
    end

end

%% Equations 
function dN = ODESystem(~, N, V, kG)

    M = length(N);

    f = zeros(M,1);
    for i=1:M
        f(i) = N(i)/(V(i+1)-V(i));
    end

    dN = zeros(M,1);
    dN(1) = -f(1)*Vdot(kG, V(1));
    for i=2:M-1
        dN(i) = f(i-1)*Vdot(kG, V(i-1)) - f(i)*Vdot(kG, V(i));
    end
    dN(M) = f(M-1)*Vdot(kG, V(M-1)) - f(M)*Vdot(kG, V(M));

end


%% Growth rate function
function Vprime = Vdot(kG, V)

     Vprime = 3*kG*V^(1/3);

end
