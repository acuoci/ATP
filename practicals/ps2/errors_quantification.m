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
%  Code: automatically run the script advection diffusion 1D FVM explicit %
%        and store the error at different mesh size in order to study     %
%        the spacial convergence of the discretization schemes.           %
%        Before running this code comment line 45 and 51 in the script    %
%        advection_diffusion_1D_FVM_explicit.m                            %
%                                                                         %
% ----------------------------------------------------------------------- %
clc; close all; clear;

ncvec = [10, 20, 30, 40, 50, 60];
sigma = 1;
Evec  = zeros(length(ncvec),1);
for runi=1:length(ncvec)
    ncells = ncvec(runi);
    run advection_diffusion_1D_FVM_explicit.m;
    Evec(runi) = E;
end

hvec = zeros(size(ncvec));
order1 = zeros(size(ncvec));
order2 = zeros(size(ncvec));
const = 2;
for runi=1:length(ncvec)
    hvec(runi) = 1/(ncvec(runi)-1);
    order1(runi) = const*(hvec(runi))^1;
    order2(runi) = const*(hvec(runi))^2;
end

figure; hold on; grid on;
plot(log10(1./hvec), log10(Evec), "o", "MarkerSize", 8, "LineWidth", 1.8);
plot(log10(1./hvec), log10(order1), "LineWidth", 1.8);
plot(log10(1./hvec), log10(order2), "LineWidth", 1.8);
legend ("Results", "1^{st} order", "2^{nd} order");
title("Space Convergence")
xlabel("log_{10} (1/h)"); ylabel("log_{10} (Error)")
