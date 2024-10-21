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
%  Code: 1D advection diffusion equation using FVM and a time implicit    %
%        discretization of the time derivative.                           %
%                                                                         %
% ----------------------------------------------------------------------- %

clc; close all; clear;

%-------------------------------------------------------------------------%
% User-defined data
%-------------------------------------------------------------------------%

ncells = 500;           % number of cells
L = 50;                 % domain length [m]
u = 1.;                 % velocity [m/s]
D = 0.1;                % diffusion coefficient [m2/s]
tau = 25;               % total simulation time [s]
fleft = 0.;             % left side boundary value

%-------------------------------------------------------------------------%
% Pre-processing of user-defined data
%-------------------------------------------------------------------------%

% Grid step calculation
h = L/ncells;           % grid step [m]
x = linspace (0, L, ncells+1);
xc = linspace (-0.5*h, L+0.5*h, ncells+2);

% Check the stability conditions on time step
dt_max = 1*h/u;
sigma = 2;
dt = sigma*dt_max;
nstep = tau/dt;

% Memory allocation
f  = zeros (ncells+2, 1);
fo = zeros (ncells+2, 1);
fp = zeros (ncells+1, 1);

% Initial solution
for i=1:ncells+2
    sigma2 = 0.2;
    sigma = sqrt(sigma2);
    mu = 0.25*L;
    f(i) = 1/(sigma*sqrt(2*pi))*exp(-(xc(i) - mu)^2/(2*sigma2));
end

% Create Matrix A
n = ncells + 2;

% Known terms vector
b = zeros (ncells+2, 1);

% Build diagonals of the global tridiagonal matrix
alpha = u*dt/h;
beta  = D*dt/h^2;
Ap = -2*beta - alpha - 1;
Ae = beta;
Aw = beta + alpha;

% Create global matrix using sparse in order to avoid to store the "0"s
A = sparse (ncells+2, ncells+2);

% Fill the global matrix using the three diagonals
A(1,1)=1; A(1,2)=1;                         % From the boundary conditions
for i=2:ncells+1, A(i,i-1) = Aw; end
for i=2:ncells+1, A(i,i) = Ap;   end
for i=2:ncells+1, A(i,i+1) = Ae; end
A(n,n)=1; A(n,n-1)=-1;                      % From the boundary conditions

% loop over all the time steps
for m=1:nstep
    
    % Store solution at the old time step
    fo = f;

    % Update RHS term using info from the old time step
    for i=2:ncells+1
        b(i) = -fo(i);
    end
    b(1) = 2*fleft;
    b(ncells+2) = 0;

    % Solve the linear system of equations
    f = A\b;

    % Plot results
    if (mod(m,1) == 0 || m == 1)
        
        % Linear interpolation, in order to interpolate cell-centered
        % temperature values to the face-centered position.
        for i=1:ncells+1
            fp(i) = 0.5*(f(i+1) + f(i));
        end

        hold off;
        plot (x, fp, "LineWidth", 2);
        ylim([0 1]);
        grid on;
        xlabel ("length [m]");
        ylabel ("Solution");
        title ("Evolution of f using time implicit discretization");
        drawnow;

        % Assembly plots in a gif
        if(0)
            filename = 'ps4-implicit.gif';
            message = sprintf('time=%g s\n', m);
            time = annotation ('textbox', [0.8 0.6 0.2 0.2], 'String', message, 'EdgeColor', 'none');
            set(gcf, 'Color', 'w');
            frame = getframe(gcf);
            img = frame2im(frame);
            [imind, cm] = rgb2ind(img, 256);
            annotation ('arrow', [0.4 0.6], [0.6 0.6]);
            text (0.47*L, 0.63, 'u');
            
            % Write to the GIF file
            if m == 1
                imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1);
            else
                imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1);
            end
        
            delete(time);
        end
    end
    
end
