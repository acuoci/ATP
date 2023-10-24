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
%  Code: 2D Poisson equation convergence rate using Gauss-Seidel and      %
%        SOR methods in a FV discretization.                              %
%                                                                         %
% ----------------------------------------------------------------------- %

clc; close all; clear;

% Run at different levels of refinement
% ncells = 2^level for each dimension
levels = [4, 5, 6, 7];
errors = zeros(size(levels));

for i=1:length(levels)
    errors(i) = run(levels(i));
end

figure;
loglog (2.^levels, errors, 'x', "LineWidth", 2); hold on;
loglog (2.^levels, 20.*(2.^levels).^-2, "LineWidth", 2); hold on;
xlabel ("N");
ylabel ("L_2");
xlim ([2^(levels(1)-1), 2^(levels(end)+1)]);
legend ("Results", "2^{nd} Order");
title ("Convergence rate of the spatial derivative");
axis("square");
grid ("on");

function err = run (level)
    
    % Data
    L = 10;
    fleft = 0;
    
    % Parameters for the iterative solution
    maxiter = 100000;
    tolerance = 1e-9;
    beta = 1.9;
    
    % Grid
    ncells = (2^level);
    h = L/ncells;
    xp = linspace (-0.5*L, 0.5*L, ncells+1);
    xc = linspace (-0.5*L-0.5*h, 0.5*L+0.5*h, ncells+2);
    
    % Create fields
    f  = zeros (ncells+2);
    fo = zeros (ncells+2);
    fp = zeros (ncells+1);
    
    % Loop over total number of iterations
    for iter=1:maxiter
        
        % Update boundary conditions
        f(1,:) = f(2,:);
        f(ncells+2,:) = f(ncells+1,:);
        f(:,1) = 2*fleft - f(:,2);
        f(:,ncells+2) = f(:,ncells+1);
    
        for i=2:ncells+1
            for j=2:ncells+1
                xi = xc(i); yi = xc(j);
    
                f(i,j) = beta/4.*((f(i+1,j) + f(i-1,j) + f(i,j+1) + f(i,j-1)) - source (xi,yi)*h^2) ...
                    + (1 - beta)*f(i,j);
            end
        end
    
        % Compute residuals
        res = 0;
        for i=2:ncells+1
            for j=2:ncells+1
                xi = xc(i); yi = xc(j);
    
                res = res + abs( (f(i+1,j) + f(i-1,j) + f(i,j+1) + f(i,j-1) ...
                    - 4.*f(i,j))/h^2 - source (xi,yi));
            end
        end
        
        % Normalize residual with the number of points
        % in order to find the residual in the single point
        res = res/((ncells)*(ncells));
    
        % If the residual is smaller than the user-defined tolerance
        if (res <= tolerance)
            % Leave the for loop
            break;
        end
    
        % If the maximum number of iteration was reached print this message
        % which tells you that the system did not converge
        if (iter==maxiter-1)
            fprintf("Maximum number of iterations reached\n");
        end
    end
    
    % Linear interpolation
    for i=1:ncells+1
        for j=1:ncells+1
            fp(i,j) = 0.25*(f(i,j) + f(i+1,j) + f(i,j+1) + f(i+1,j+1));
        end
    end
    
    % Plot results as a 2D surface
    surf (xp, xp, fp');
    axis square;
    colormap jet;
    shading interp;
    view (2);
    xlabel ("x");
    ylabel ("y");
    colorbar;
    
    % Compute the error (L2-norm)
    err = 0;
    for i=2:ncells+1
        for j=2:ncells+1
            xi = xc(i); yi = xc(j);
    
            el = abs (f(i,j) - exact(xi, yi));
            err = err + el^2;
        end
    end
    err = h*sqrt (err); 

end

%-- Helper functions

function s = source (x, y)
    s = (4*(x^2 + y^2) - 4)*exp(-x^2 - y^2);
end

function e = exact (x, y)
    e = exp(-x^2 - y^2);
end
