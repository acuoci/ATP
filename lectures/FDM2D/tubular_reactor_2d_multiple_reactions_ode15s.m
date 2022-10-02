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
global Nx Ny hx hy u v gamma kappa Cin

% A->B      r1 = k1*CA
% A+B->C    r2 = k2*CA*CB

%-------------------------------------------------------------------------%
% Data
%-------------------------------------------------------------------------%
Lx = 1.00;              % reactor length (m)
Ly = 0.10;              % reactor width (m)
mu = 1e-3;              % viscosity (kg/m/s)
gamma = 0.242e-4;       % diffusion coefficient (m2/s)
kappa = [ 0.1, 0.05];   % kinetic constants (1/s, m3/kmol/s)
uIn = 0.1;              % axial velocity (average) (m/s)
tf = 25;                % total time of simulation (s)
Nx = 24;                % number of grid points along x (-)
Ny = 16;                % numebr of grid points along y (-)
Cin = [1 0 0];          % inlet concentrations (kmol/m3)


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
M = ones(Nx,Ny);
M(:,1) = 0;                 % south side
M(:,Ny) = 0;                % north side
M(1,:) = 0;                 % west side
M(Nx,:) = 0;                % east side
Mdiag = diag([M(:);M(:);M(:)]);
options = odeset('Mass',Mdiag);

%-------------------------------------------------------------------------%
% Solution via ode15s solver
%-------------------------------------------------------------------------%

% Initial solution
CA = zeros(Nx,Ny);
CB = zeros(Nx,Ny);
CC = zeros(Nx,Ny);
CA(:,1) = Cin(1)*ones(Nx,1);
CB(:,1) = Cin(2)*ones(Nx,1);
CC(:,1) = Cin(3)*ones(Nx,1);
C = [CA(:); CB(:); CC(:)];

% DAE solver
[t, C] = ode15s(@ODESystem, 0:0.5:tf, C(:), options);

%-------------------------------------------------------------------------%
% Video setup
%-------------------------------------------------------------------------%
video_name = 'tubular_reactor_2d_ode15s.mp4';
videompg4 = VideoWriter(video_name, 'MPEG-4');
open(videompg4); figure;

for k=1:length(t)
     CA = C(k, 1:Nx*Ny);
     CB = C(k, Nx*Ny+1:2*Nx*Ny);
     CC = C(k, 2*Nx*Ny+1:3*Nx*Ny);

     
     hold off;

     subplot(311);
     solution = reshape( CA, Nx,Ny)';
     surf(X,Y, solution)
     colormap(jet); shading interp; colorbar; view(2);
     clim([0 max(max(C(:,1:Nx*Ny)))]);
     hold on;
     xlabel('x [m]'); ylabel('y [m]');    
     message = sprintf('time=%f', t(k));
     time = annotation('textbox',[0.15 0.8 0.1 0.1],'String',message,'EdgeColor','none', 'Color', 'white');

     subplot(312);
     solution = reshape( CB, Nx,Ny)';
     surf(X,Y, solution)
     colormap(jet); shading interp; colorbar; view(2);
     clim([0 max(max(C(:,1*Nx*Ny+1:2*Nx*Ny)))]);
     hold on;
     xlabel('x [m]'); ylabel('y [m]');    
     
     subplot(313);
     solution = reshape( CC, Nx,Ny)';
     surf(X,Y, solution)
     colormap(jet); shading interp; colorbar; view(2);
     clim([0 max(max(C(:,2*Nx*Ny+1:3*Nx*Ny)))]);
     hold on;
     xlabel('x [m]'); ylabel('y [m]');    
     
     
     frame = getframe(gcf);
     writeVideo(videompg4,frame);
     delete(time);
end

% Closing the video stream
close(videompg4);



function dC_over_dt = ODESystem(~, C)

    global Nx Ny hx hy u v gamma kappa Cin

    % Allocate memory
    dCA_over_dt = zeros(Nx,Ny);
    dCB_over_dt = zeros(Nx,Ny);
    dCC_over_dt = zeros(Nx,Ny);

    % Reshape C vector into a matrix
    CA = reshape(C(1:Nx*Ny),Nx,Ny);
    CB = reshape(C(Nx*Ny+1:2*Nx*Ny),Nx,Ny);
    CC = reshape(C(2*Nx*Ny+1:3*Nx*Ny),Nx,Ny);

    % Boundaries
    for i=1:Nx
        dCA_over_dt(i,1)  = CA(i,1) - CA(i,2);             % south (Neumann)
        dCA_over_dt(i,Ny) = CA(i,Ny) - CA(i,Ny-1);         % north (Neumann)
        dCB_over_dt(i,1)  = CB(i,1) - CB(i,2);             % south (Neumann)
        dCB_over_dt(i,Ny) = CB(i,Ny) - CB(i,Ny-1);         % north (Neumann)
        dCC_over_dt(i,1)  = CC(i,1) - CC(i,2);             % south (Neumann)
        dCC_over_dt(i,Ny) = CC(i,Ny) - CC(i,Ny-1);         % north (Neumann)
    end
    for j=2:Ny-1
        dCA_over_dt(1,j)  = u(1,j)*(CA(1,j)-Cin(1)) + ...
                           -gamma*(CA(2,j)-CA(1,j))/hx;     % west (Danckwerts)
        dCA_over_dt(Nx,j) = CA(Nx,j) - CA(Nx-1,j);          % east (Neumann)
        dCB_over_dt(1,j)  = u(1,j)*(CB(1,j)-Cin(2)) + ...
                           -gamma*(CB(2,j)-CB(1,j))/hx;     % west (Danckwerts)
        dCB_over_dt(Nx,j) = CB(Nx,j) - CB(Nx-1,j);          % east (Neumann)
        dCC_over_dt(1,j)  = u(1,j)*(CC(1,j)-Cin(3)) + ...
                           -gamma*(CC(2,j)-CC(1,j))/hx;     % west (Danckwerts)
        dCC_over_dt(Nx,j) = CC(Nx,j) - CC(Nx-1,j);          % east (Neumann)        
    end

    % Internal points
    for i=2:Nx-1

        for j=2:Ny-1

            r = [kappa(1)*CA(i,j), kappa(2)*CA(i,j)*CB(i,j)];
            R = [-r(1)-r(2), r(1)-r(2), r(2)];

            dCA_over_dx = (CA(i,j)-CA(i-1,j))/hx;
            d2CA_over_dx2 = (CA(i+1,j)-2*CA(i,j)+CA(i-1,j))/hx^2;
            dCB_over_dx = (CB(i,j)-CB(i-1,j))/hx;
            d2CB_over_dx2 = (CB(i+1,j)-2*CB(i,j)+CB(i-1,j))/hx^2;
            dCC_over_dx = (CC(i,j)-CC(i-1,j))/hx;
            d2CC_over_dx2 = (CC(i+1,j)-2*CC(i,j)+CC(i-1,j))/hx^2;

            dCA_over_dy = (CA(i,j)-CA(i,j-1))/hy;
            d2CA_over_dy2 = (CA(i,j+1)-2*CA(i,j)+CA(i,j-1))/hy^2; 
            dCB_over_dy = (CB(i,j)-CB(i,j-1))/hy;
            d2CB_over_dy2 = (CB(i,j+1)-2*CB(i,j)+CB(i,j-1))/hy^2;                 
            dCC_over_dy = (CC(i,j)-CC(i,j-1))/hy;
            d2CC_over_dy2 = (CC(i,j+1)-2*CC(i,j)+CC(i,j-1))/hy^2;     


            dCA_over_dt(i,j) =   -u(i,j)*dCA_over_dx -v(i,j)*dCA_over_dy + ...
                                gamma*(d2CA_over_dx2 + d2CA_over_dy2) + ...
                                R(1);

            dCB_over_dt(i,j) =   -u(i,j)*dCB_over_dx -v(i,j)*dCB_over_dy + ...
                                gamma*(d2CB_over_dx2 + d2CB_over_dy2) + ...
                                R(2);

            dCC_over_dt(i,j) =   -u(i,j)*dCC_over_dx -v(i,j)*dCC_over_dy + ...
                                gamma*(d2CC_over_dx2 + d2CC_over_dy2) + ...
                                R(3);            
        end

    end

    % Unrolling (flattening) the matrix
    dC_over_dt = [ dCA_over_dt(:); dCB_over_dt(:); dCC_over_dt(:)];

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
