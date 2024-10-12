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
%  Code: 1D advection-diffusion equation by the FV method                 %
%        solution using explicit in time discretization method            %
%                                                                         %
% ----------------------------------------------------------------------- %

clc; close all; clear;

% ----------------------------------------------------------------------- %
% Pre-Processing
% ----------------------------------------------------------------------- %

global u D;

% Data
L = 2.0;        % domain length [m]
u = 1;          % velocity [m/s]
D = 0.05;       % diffusion coefficient [m2/s]
tau = 2*pi;     % total integration time [s] 

% Analytical/Exact solution
function f = exact (x, t)
    global u D;
    
    A = 0.5;        % amplitude of initial solution
    k = 1;          % wave number [1/m]
    
    f = A*exp(-4*pi^2*k^2*D*t)*sin(2*pi*k*(x - u*t));
end

% Build mesh
ncells = 50;                 % number of grid points    
h = L/ncells;                % grid step [m]
x = linspace(0,L,ncells+1);  % point coordinates (for post-processing)
xc = -0.5*h:h:L+0.5*h;       % cell-centers coordinates (including ghosts)

% Build Time
dt_conv = 2*D/u^2;                  % Maximum time step for convection
dt_diff = h^2/2/D;                  % Maximum time step for diffusion
dt_max = min (dt_conv, dt_diff);    % Maximum allowed time step
sigma = 0.9;                        % Regulates dt with respect to dt_max
dt = sigma*dt_max;                  % Final simulation time step
nsteps = tau/dt;

fprintf('dt_conv=%g, dt_diff=%g, dt=%g\n', dt_conv, dt_diff, dt);

% Memory allocation
fo = zeros(ncells+2,1);      % temporary numerical solution
f  = zeros(ncells+2,1);      % current numerical solution
fp = zeros(ncells+1,1);      % solution interpolated on the nodes (for plot)
a  = zeros(ncells+2,1);      % exact solution
ap = zeros(ncells+1,1);      % exact solution (for plot)

% Initial solution
for i=1:ncells+2
    f(i) = exact (xc(i), 0.);
end

% ----------------------------------------------------------------------- %
% Solution loop
% ----------------------------------------------------------------------- %

t = 0.;
for m=1:nsteps
    
    % New time step
    t = t+dt;       

    % Update the analytical solution
    for i=1:ncells+2
        a(i) = exact (xc(i), t);
    end

    % Store old field
    fo = f;

    % Advance the solution in the internal cells except the last one
    for i=2:ncells
        Ai = u/2/h*(fo(i+1) - fo(i-1));
        Di = D/h^2*(fo(i+1) + fo(i-1) - 2*fo(i));
        f(i) = fo(i) + dt*(-Ai+Di);
    end

    % Transport in the last internal cell (@ncells+1) considering cell (@2)
    % as neighboring cell
    i = ncells+1; ip1 = 2; im1 = ncells;
    Ai = u/2/h*(fo(ip1) - fo(im1));
    Di = D/h^2*(fo(ip1) + fo(im1) - 2*fo(i));
    f(i) = fo(i) + dt*(-Ai+Di);

    % Find west wall value interpolating (@ncells+1) and (@2)
    fwall = 0.5*(f(ncells+1) + f(2));

    % Impose the west wall value at the ghost cells for the next iteration
    f(ncells+2) = 2*fwall - f(ncells+1);
    f(1)        = 2*fwall - f(2);

    if (mod(m,5) == 1)
    
        % Interpolate f to nodes
        for i=1:ncells+1
            fp(i) = 0.5*(f(i+1) + f(i));
            ap(i) = 0.5*(a(i+1) + a(i));
        end

        % Graphical output
        hold off; plot(0:h:L,fp,'linewidth',2); axis([0 L -1, 1]);
	    hold on; plot(0:h:L,ap,'r','linewidth',2);
        hold on; legend('numerical', 'exact');
        xlabel('spatial coordinate [m]');
        ylabel('f');
        title('Evolution of the solution in space and time');
        grid on;
        drawnow;
    
        if(0)
            filename = 'ps3-plots.gif';    
            message = sprintf('time=%g s\n', t);
            time = annotation ('textbox', [0.7 0.05 0.2 0.2], 'String', message, 'EdgeColor', 'none');
            set(gcf, 'Color', 'w');
            frame = getframe(gcf);
            img = frame2im(frame);
            [imind, cm] = rgb2ind(img, 256);
        
            % Write to the GIF file
            if m == 1
                imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
            else
                imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
            end
        
            delete(time);
        end
    end


end
