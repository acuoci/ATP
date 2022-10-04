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
%  Code: Lorenz's equations                                               %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

global sigma beta rho

%% User-defined data
%-------------------------------------------------------------------------%
sigma = 10;
beta = 8/3;
rho = 28;


%% Numerical solution
%-------------------------------------------------------------------------%

% Initial conditions: case 1
psiInitial = [0.1 0.1 0.1]';
[t1,psi1] = ode15s(@LorenzEquations, 0:0.1:100, psiInitial);

% Initial conditions: case 2 (perturbed)
psiInitial = [0.1000001 0.1 0.1]';
[t2,psi2] = ode15s(@LorenzEquations, 0:0.1:100, psiInitial);


%% Plotting the results
%-------------------------------------------------------------------------%

% Plotting x vs time
subplot(3,1,1);
plot(t1, psi1(:,1), 'b');
title ('Case 1'); ylabel('x');
subplot(3,1,2);
plot(t2, psi2(:,1), 'r');
title('Case 2'); ylabel('x');
subplot(3,1,3)
plot(t1, psi1(:,1)-psi2(:,1), 'g');
title('Difference'); ylabel('x_1-x_2');
xlabel('time');


%% Lorenz Equations
%-------------------------------------------------------------------------%
function dpsidt = LorenzEquations(t,psi)

    global sigma beta rho

    dpsidt = zeros(3,1);

    dpsidt(1) = sigma*(psi(2)-psi(1));
    dpsidt(2) = rho*psi(1)-psi(2)-psi(1)*psi(3);
    dpsidt(3) = -beta*psi(3)+psi(1)*psi(2);

end