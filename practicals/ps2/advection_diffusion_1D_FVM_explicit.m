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

%-------------------------------------------------------------------------%
% User-defined data
%-------------------------------------------------------------------------%

ncells = 60;         % number of cells
np = ncells + 1;      % number of grid points
dt = 0.01;           % time step [s]
nstep = 100;         % number of time steps
%dt = 0.5*0.004;
%tottime = 0.5;        % total simulation seconds  
%nstep = ceil(tottime/dt);   % number of time steps
L = 2.0;              % domain length [m]
u = 1;                % velocity [m/s]
D = 0.05;             % diffusion coefficient [m2/s]
A = 0.5;              % amplitude of initial solution
k = 1;                % wave number [1/m]

%-------------------------------------------------------------------------%
% Pre-processing of user-defined data
%-------------------------------------------------------------------------%

% Grid step calculation
% h = L/(np-1);         % grid step [m]
h = L/ncells;         % grid step [m]
x = linspace(0,L,np);

% Memory allocation
fo = zeros(ncells+2,1);      % temporary numerical solution
f  = zeros(ncells+2,1);      % current numerical solution
fp = zeros(ncells+1,1);      % solution interpolated on the nodes (for plot)
a  = zeros(ncells+2,1);      % exact solution
ap = zeros(ncells+1,1);

% Error fields
eps  = zeros(ncells+2,1);     % Relative error
L1   = 0.;                    % L1 = volume_integral(eps dV) / volume_integral(dV)
Linf = 0.;                    % Linf = max (eps)

% Initial solution
for i=1:ncells+2
    xi = (h*(i-2) + 0.5*h);
    f(i) = A*sin(2*pi*k*xi);
end

% Check the stability conditions on time step
Co = u*dt/h;                        % Courant number
Di = D*dt/h^2;                      % Diffusion number
dt_max = min(1*h/u, 0.5*h*h/D);     % Maximum allowed time step
fprintf('Co=%f, Di=%f, dt=%f, dt(max)=%f\n', Co, Di, dt, dt_max);

dt = 0.8*dt_max;

%-------------------------------------------------------------------------%
% Video setup
%-------------------------------------------------------------------------%
% video_name = 'advection_diffusion_1d.mp4';
% videompg4 = VideoWriter(video_name, 'MPEG-4');
% open(videompg4);

%-------------------------------------------------------------------------%
% Advancing in time
%-------------------------------------------------------------------------%

t = 0.;
for m=1:nstep
    
    % Update the analytical solution
    for i=1:ncells+2
        xi = (h*(i-2) + 0.5*h);
        a(i) = A*exp(-4*pi*pi*k*k*D*t)*sin(2*pi*k*(xi-u*t));
    end

    % Interpolate f to nodes
    for i=1:np
        fp(i) = 0.5*(f(i+1) + f(i));
        ap(i) = 0.5*(a(i+1) + a(i));
    end

    % Squared areas below the analytical and numerical solutions
    a2_int = 0.;
    f2_int = 0.;
    for i=1:np-1
         a2_int = a2_int + h/2*(ap(i)^2+ap(i+1)^2);
         f2_int = f2_int + h/2*(fp(i)^2+fp(i+1)^2);
    end  
    
    % Graphical output
    message = sprintf('time=%d\na^2(int)=%d\ny^2(int)=%d', t, a2_int, f2_int);
	hold off; plot(0:h:L,fp,'linewidth',2); axis([0 L -1, 1]); % plot num. 
	hold on; plot(0:h:L,ap,'r','linewidth',2);                 % plot exact
    hold on; legend('numerical', 'exact');                    % legend
    xlabel('spatial coordinate [m]');
    ylabel('solution');    
    time = annotation('textbox',[0.15 0.8 0.1 0.1],'String',message,'EdgeColor','none');
    frame = getframe(gcf);
    % writeVideo(videompg4,frame);
    delete(time);

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

    % Interpolate f to nodes
    for i=1:np
        fp(i) = 0.5*(f(i+1) + f(i));
        ap(i) = 0.5*(a(i+1) + a(i));
    end

    % Update the error between numerical and analytical solution
    E = 0;
    for i=1:np
        E = E + (fp(i)-ap(i))^2;
    end
    E = h*sqrt(E);

    % Error quantification
    for i=2:ncells+1
        eps(i) = abs(f(i) - a(i));
    end

    L1 = 0; Vol = 0;
    for i=2:ncells+1
        L1 = L1 + eps(i)*h^2;
        Vol = Vol + h^2;
    end
    L1 = L1/Vol;
    Linf = max(eps);

    % New time step
    t=t+dt;
       
    % Print the current time (every 25 steps)
    if (mod(m,25)==1), fprintf('time=%d E=%d L1=%d Linf=%d\n', t, E, L1, Linf); end
end

fprintf('time=%d E=%e\n', t, E);


% Closing the video stream
% close(videompg4);
