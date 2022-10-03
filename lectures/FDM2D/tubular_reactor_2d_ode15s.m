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
%  Code: 2D advection-diffusion-reaction by the FD method                 %
%        solution via the ode15s solver                                   %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% Global variables (meaning reported below)
global Nx Ny hx hy u v gamma kappa CinMean f A


%-------------------------------------------------------------------------%
% Data
%-------------------------------------------------------------------------%
Lx = 1.00;          % reactor length (m)
Ly = 0.10;          % reactor width (m)
mu = 1e-3;          % viscosity (kg/m/s)
gamma = 0.242e-4;   % diffusion coefficient (m2/s)
kappa = 0.1;        % kinetic constant (1/s)
uIn = 0.1;          % axial velocity (average) (m/s)
tf = 100;           % total time of simulation (s)
Nx = 50;            % number of grid points along x (-)
Ny = 20;            % numebr of grid points along y (-)

% Sinusoidal inlet concentration
CinMean = 1;        % inlet concentration (kmol/m3)
f = 0.025;          % frequency (1/s)    
A = 0.50;             % amplitude (-)

%-------------------------------------------------------------------------%
% Preprocessing
%-------------------------------------------------------------------------%
hx = Lx/(Nx-1);     % grid spacing along x (m)
hy = Ly/(Ny-1);     % grid spacing along y (m)

% Assembling the computational grid (mesh)
x = 0:hx:Lx;                % x axis (m)
y = 0:hy:Ly;                % y axis (m)
[X,Y] = meshgrid(x,y);      % mesh object (for graphical purposes)

% Velocity Field
[u, v] = VelocityField(Nx,Ny,Ly,y,uIn,mu);

% Mass matrix for defining the algebraic equations along the boundaries
Mdiag = ones(Nx,Ny);
Mdiag(:,1) = 0;         % south side
Mdiag(:,Ny) = 0;        % north side
Mdiag(1,:) = 0;         % west side
Mdiag(Nx,:) = 0;        % east side
options = odeset('Mass',diag(Mdiag(:)));

%-------------------------------------------------------------------------%
% Solution via ode15s solver
%-------------------------------------------------------------------------%

% Initial solution
C = zeros(Nx,Ny);
C(:,1) = Cin(0)*ones(Nx,1);

% DAE solver
[t, C] = ode15s(@ODESystem, 0:1:tf, C(:), options);


%-------------------------------------------------------------------------%
% Video setup
%-------------------------------------------------------------------------%
video_name = 'tubular_reactor_2d_ode15s.mp4';
videompg4 = VideoWriter(video_name, 'MPEG-4');
open(videompg4); figure;

for k=1:length(t)
     solution = reshape( C(k,:), Nx,Ny)';
     hold off;
     surf(X,Y, solution)
     colormap(jet); shading interp; colorbar; view(2);
     clim([0 CinMean*(1+A)]);
     hold on;
     xlabel('x [m]'); ylabel('y [m]');    
     message = sprintf('time=%f', t(k));
     time = annotation('textbox',[0.15 0.8 0.1 0.1],'String',message,'EdgeColor','none', 'Color', 'white');
     frame = getframe(gcf);
     writeVideo(videompg4,frame);
     delete(time);
end

% Closing the video stream
close(videompg4);



function dC_over_dt = ODESystem(t, C)

    global Nx Ny hx hy u v gamma kappa

    % Allocate memory
    dC_over_dt = zeros(Nx,Ny);

    % Reshape C vector into a matrix
    C = reshape(C,Nx,Ny);

    % Boundaries
    for i=1:Nx
        dC_over_dt(i,1)  = C(i,1) - C(i,2);             % south (Dirichlet)
        dC_over_dt(i,Ny) = C(i,Ny) - C(i,Ny-1);         % north (Neumann)
    end
    for j=2:Ny-1
        dC_over_dt(1,j)  = u(1,j)*(C(1,j)-Cin(t)) + ...
                           -gamma*(C(2,j)-C(1,j))/hx;   % west (Danckwerts)
        dC_over_dt(Nx,j) = C(Nx,j) - C(Nx-1,j);         % east (Neumann)
    end

    % Internal points
    for i=2:Nx-1

        for j=2:Ny-1

            dC_over_dx = (C(i,j)-C(i-1,j))/hx;
            d2C_over_dx2 = (C(i+1,j)-2*C(i,j)+C(i-1,j))/hx^2;

            dC_over_dy = (C(i,j)-C(i,j-1))/hy;
            d2C_over_dy2 = (C(i,j+1)-2*C(i,j)+C(i,j-1))/hy^2;        

            dC_over_dt(i,j) =   -u(i,j)*dC_over_dx -v(i,j)*dC_over_dy + ...
                                gamma*(d2C_over_dx2 + d2C_over_dy2) + ...
                                -kappa*C(i,j);
        end

    end

    % Unrolling (flattening) the matrix
    dC_over_dt = dC_over_dt(:);

end


% Velocity field
% Analytical solution for laminar flow between parallel plates
function [u, v] = VelocityField(Nx,Ny,Ly,y,uIn,mu)

    u = zeros(Nx,Ny);
    v = zeros(Nx,Ny);
        
    G = 12*mu*uIn/Ly^2;

    for i=1:Nx
        for j=1:Ny
            u(i,j) = G/(2*mu)*y(j)*(Ly-y(j));
        end
    end

end

% Dynamic inlet concentration
function CinDynamic = Cin(t)

    global CinMean f A

    CinDynamic = CinMean*(1+A*sin(2*pi*f*t));

end