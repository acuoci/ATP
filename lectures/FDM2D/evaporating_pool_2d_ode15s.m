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
global Nx Ny hx hy u v gamma kappa Cvp


%-------------------------------------------------------------------------%
% Data
%-------------------------------------------------------------------------%
Lx = 0.20;          % horizontal length (m)
Ly = 0.10;          % vertical length (m)
Pv = 2300;          % water vapor pressure at 298K (Pa)
T = 298;            % temperature (K)
gamma = 0.242e-4;   % diffusion coefficient of water in air (m2/s)
nu = 1e-5;          % air kinematic viscosity (m2/s)
kappa = 0.;         % kinetic constant (1/s)
uInf = 0.03;         % free stream horizontal velocity (m/s)
tf = 100;           % total time of simulation (s)
Lplate = 0.10;      % plate length (before water) (m)
Nx = 20;            % number of grid points along x (-)
Ny = 30;            % numebr of grid points along y (-)


%-------------------------------------------------------------------------%
% Preprocessing
%-------------------------------------------------------------------------%
hx = Lx/(Nx-1);     % grid spacing along x (m)
hy = Ly/(Ny-1);     % grid spacing along y (m)
Cvp = Pv/8314/T;    % water concentration on the surface (kmol/m3)

% Assembling the computational grid (mesh)
x = 0:hx:Lx;                % x axis (m)
y = 0:hy:Ly;                % y axis (m)
[X,Y] = meshgrid(x,y);      % mesh object (for graphical purposes)

% Velocity Field
[u, v] = VelocityField(Nx,Ny,x,y,uInf,Lplate,nu);

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
C(:,1) = Cvp;

% DAE solver
[t, C] = ode15s(@ODESystem, 0:1:tf, C(:), options);


%-------------------------------------------------------------------------%
% Video setup
%-------------------------------------------------------------------------%
video_name = 'evaporating_pool_2d_ode15s.mp4';
videompg4 = VideoWriter(video_name, 'MPEG-4');
open(videompg4); figure;

for k=1:length(t)
     solution = reshape( C(k,:), Nx,Ny)';
     hold off;
     surf(X,Y, solution)
     colormap(jet); shading interp; colorbar; view(2);
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



function dC_over_dt = ODESystem(~, C)

    global Nx Ny hx hy u v gamma kappa Cvp

    % Allocate memory
    dC_over_dt = zeros(Nx,Ny);

    % Reshape C vector into a matrix
    C = reshape(C,Nx,Ny);

    % Boundaries
    for i=1:Nx
        dC_over_dt(i,1)  = C(i,1) - Cvp;                % south (Dirichlet)
        dC_over_dt(i,Ny) = C(i,Ny) - C(i,Ny-1);         % north (Neumann)
    end
    for j=2:Ny-1
        dC_over_dt(1,j)  = u(1,j)*(C(1,j)-0) + ...
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


% Velocity field along the horizontal direction (Blasius velocity profile)
% See: Savas, An approximate compact analytical expression for the 
%      Blasius velocity profile, Communications in Nonlinear Science and 
%      Numerical Simulation, 17(10), p. 3772-3775 (2012)

function [u, v] = VelocityField(Nx,Ny,x,y,uInf,Lplate,nu)

    a = 0.33206;        % parameter in the analytical correlation
    n = 1.50;           % parameter in the analytical correlation
    u = zeros(Nx,Ny);
    v = zeros(Nx,Ny);
    
    [etaVec, fVec] = ode45(@BlasiusFunction, 0:0.1:1000, 0, [], a,n);
    
    for i=1:Nx
        xc = x(i)+Lplate;
        for j=1:Ny
            yc = y(j);
            eta = sqrt(uInf/(nu*xc))*yc;
            fprime = (tanh((a*eta)^n))^(1/n);
            f = interp1(etaVec,fVec, eta);
            u(i,j) = uInf*fprime;
            v(i,j) = sqrt(nu*uInf/xc)* 0.50*(eta*fprime-f);
        end
    end

end

function df_over_deta = BlasiusFunction(eta,~, a,n)

    df_over_deta = (tanh((a*eta)^n))^(1/n);

end