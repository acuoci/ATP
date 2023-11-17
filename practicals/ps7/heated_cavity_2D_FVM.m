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
%  Code: Heated cavity test case - solution of incompressible             %
%        Navier-Stokes equations, with flow driven by natural             %
%        convective fluxes, described using the Bousinesq                 %
%        approximation. Problem setup from Ferziger Peric book            %
%        4th ed., pag. 264.                                               %
%                                                                         %
% ----------------------------------------------------------------------- %

clc; close all; clear;

% Pre-Processing
L = 1.;             % Length of the domain [m]
nu = 1.e-3;         % [HEATED] Kinematic viscosity [m2/s]
tau = 40;           % [HEATED] Total simulation time [s]
level = 5;          % Maximum level of refinement

% [HEATED] Pre-Processing
Tinit = 6.;         % Initial temperature value
Thot = 10.;         % Temperature of the bottom wall
Tcold = 0.;         % temperature of the cold wall
rho = 1.;           % baseline/initial density
g = 10;             % gravitational acceleration
betaexp = 0.01;     % thermal expansion coefficient
Pr = 0.1;           % Prandtl number
alpha = nu*rho/Pr;  % thermal diffusivity

% Grid setup
Lx = L;
Ly = L;
nx = 2^level;
hx = Lx/nx;
hy = hx;
ny = Ly/hy;
x = linspace (0., L, nx+1);
h = hx;

% Dirichlet boundary conditions everywhere
unwall = 0.;      % [HEATED]
uswall = 0.;
vewall = 0.;
vwwall = 0.;

% Poisson solver settings
maxiter = 10000;
beta = 1.6;
tolerance = 1.e-6;

% Time step setup
% sigma = 0.5;
% dt_diff = h^2/4/nu;
% dt_conv = 4*nu/unwall^2;
% dt = sigma*min (dt_diff, dt_conv);
dt = 0.005;
nsteps = tau/dt;
Re = unwall*L/nu;

% Print initial info
fprintf ("Time step = %f - Re = %f\n", dt, Re);

% Memory allocations
u = zeros (nx+1, ny+2);
v = zeros (nx+2, ny+1);
p = zeros (nx+2, ny+2);

ut = zeros (nx+1, ny+2);
vt = zeros (nx+2, ny+1);

up = zeros (nx+1, ny+1);
vp = zeros (nx+1, ny+1);
pp = zeros (nx+1, ny+1);

% [HEATED] Memory allocations
T  = zeros (nx+2, ny+2) + Tinit;
To = zeros (nx+2, ny+2) + Tinit;
Tp = zeros (nx+1, ny+1);

% Coefficient for pressure equation
gamma = zeros (nx+2, ny+2);
gamma(:,:) = 1/4;
gamma(2,:) = 1/3;
gamma(:,2) = 1/3;
gamma(end-1,:) = 1/3;
gamma(:,end-1) = 1/3;
gamma(2,2) = 1/2;
gamma(2,end-1) = 1/2;
gamma(end-1,2) = 1/2;
gamma(end-1,end-1) = 1/2;

% Time solution loop
time = 0.;
for m=1:nsteps
    time = time + dt;

    %-- Set boundary conditions

    u(:,1)   = 2*uswall - u(:,2);        % south wall
    u(:,end) = 2*unwall - u(:,end-1);    % north wall
    v(1,:)   = 2*vwwall - v(2,:);        % west wall
    v(end,:) = 2*vewall - v(end-1,:);    % east wall

    %-- [HEATED] Set boundary conditions for temperature
    
    T(:,1) = 2.*Thot - T(:,2);              % south wall
    T(:,end) = 2.*Tcold - T(:,end-1);       % north wall
    T(1,:) = T(2,:);                        % west wall
    T(end,:) = T(end-1,:);                  % east wall

    %-- [HEATED] Solution of the temperature equation using u at time t

    To = T;

    for i=2:nx+1
        for j=2:ny+1
            Te = 0.5*(To(i+1,j) + To(i,j));
            Tw = 0.5*(To(i,j) + To(i-1,j));
            Tn = 0.5*(To(i,j+1) + To(i,j));
            Ts = 0.5*(To(i,j) + To(i,j-1));

            ue = u(i,j);
            uw = u(i-1,j);
            vn = v(i,j);
            vs = v(i,j-1);

            Aij = (Te*ue - Tw*uw + Tn*vn - Ts*vs)/h;
            Dij = alpha/h^2*(To(i+1,j) + To(i-1,j) + To(i,j+1) + To(i,j-1) - 4.*To(i,j));

            T(i,j) = To(i,j) + dt*(-Aij + Dij);
        end
    end

    %-- Prediction: find temporary velocity

    for i=2:nx
        for j=2:ny+1
            ue = 0.5*(u(i+1,j) + u(i,j));
            uw = 0.5*(u(i,j) + u(i-1,j));
            un = 0.5*(u(i,j+1) + u(i,j));
            us = 0.5*(u(i,j) + u(i,j-1));
            vn = 0.5*(v(i+1,j) + v(i,j));
            vs = 0.5*(v(i+1,j-1) + v(i,j-1));

            Aij = (ue^2 - uw^2 + un*vn - us*vs)/h;
            Dij = (nu/h^2)*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4.*u(i,j));

            ut(i,j) = u(i,j) + dt*(-Aij + Dij);
        end
    end

    for i=2:nx+1
        for j=2:ny
            vn = 0.5*(v(i,j+1) + v(i,j));
            vs = 0.5*(v(i,j) + v(i,j-1));
            ve = 0.5*(v(i+1,j) + v(i,j));
            vw = 0.5*(v(i,j) + v(i-1,j));
            ue = 0.5*(u(i,j+1) + u(i,j));
            uw = 0.5*(u(i-1,j+1) + u(i-1,j));

            Aij = (ve*ue - vw*uw + vn^2 - vs^2)/h;
            Dij = (nu/h^2)*(v(i+1,j) + v(i-1,j) + v(i,j+1) + v(i,j-1) - 4.*v(i,j));

            vt(i,j) = v(i,j) + dt*(-Aij + Dij);
        end
    end

    %-- Acceleration: we add the gravity acceleration term (along y)

    for i=2:nx
        for j=2:ny+1
            acc = betaexp*rho*g*(T(i,j) - Tinit);
            ut(i,j) = ut(i,j) + dt*acc;
        end
    end

    %-- Projection: find the pressure that statisfies the continuity

    for iter=1:maxiter

        % Update pressure using Jacobi/Gauss-Seidel/SOR
        for i=2:nx+1
            for j=2:ny+1
                delta = p(i+1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1);
                S = h/dt*(ut(i,j) - ut(i-1,j) + vt(i,j) - vt(i,j-1));
                p(i,j) = beta*gamma(i,j)*(delta - S) + (1 - beta)*p(i,j);
            end
        end

        % Compute residuals
        res = 0.;
        for i=2:nx+1
            for j=2:ny+1
                delta = p(i+1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1);
                S = h/dt*(ut(i,j) - ut(i-1,j) + vt(i,j) - vt(i,j-1));
                res = res + abs ( p(i,j) - gamma(i,j)*(delta - S) );
            end
        end
        res = res/(nx*ny);

        % Check residuals
        if (res <= tolerance)
            break;
        end

        % Check maxiter
        if (iter == maxiter-1)
            fprintf ("WARNING: Maximum number of iteration reached.\n");
        end
    end

    if (mod(m, 20) == 0)
        fprintf ("time = %f - Poisson iterations = %d\n", time, iter);
    end

    %-- Correction: update velocity according to the new pressure gradient

    for i=2:nx
        for j=2:ny+1
            u(i,j) = ut(i,j) - dt/h*(p(i+1,j) - p(i,j));
        end
    end

    for i=2:nx+1
        for j=2:ny
            v(i,j) = vt(i,j) - dt/h*(p(i,j+1) - p(i,j));
        end
    end

    if (mod(m, 50) == 0)
        for i=1:nx+1
            for j=1:ny+1
                Tp(i,j) = 0.25*(T(i,j) + T(i+1,j) + T(i,j+1) + T(i+1,j+1));
            end
        end

        figure(1);
        surf (x, x, Tp');
        view (2);
        xlabel ("x");
        ylabel ("y");
        title ("T");
        axis square;
        colormap jet;
        colorbar;
        shading interp;
    end
end

%-- Linear interpolations
for i=1:nx+1
    for j=1:ny+1
        pp(i,j) = 0.25*(p(i,j) + p(i+1,j) + p(i,j+1) + p(i+1,j+1));
        % [HEATED]
        Tp(i,j) = 0.25*(T(i,j) + T(i+1,j) + T(i,j+1) + T(i+1,j+1));
    end
end

for i=1:nx+1
    for j=1:ny+1
        up(i,j) = 0.5*(u(i,j+1) + u(i,j));
    end
end

for i=1:nx+1
    for j=1:ny+1
        vp(i,j) = 0.5*(v(i+1,j) + v(i,j));
    end
end

figure(2);
subplot (2,3,1);
surf (x, x, up');
view (2);
xlabel ("x");
ylabel ("y");
title ("u");
axis square;
colormap jet;
colorbar;
shading interp;

subplot (2,3,2);
surf (x, x, vp');
view (2);
xlabel ("x");
ylabel ("y");
title ("v");
axis square;
colormap jet;
colorbar;
shading interp;

subplot (2,3,3);
surf (x, x, Tp');
view (2);
xlabel ("x");
ylabel ("y");
title ("T");
axis square;
colormap jet;
colorbar;
shading interp;

% Read reference solution data from file
ref_u_along_y = dlmread('benchmark-solution', '', 1, 0);

u_profile = up(nx/2+1,:);
v_profile = vp(:,ny/2+1);
p_vertical_profile = pp(nx/2+1,:);
p_horizontal_profile = pp(:,ny/2+1);

subplot (2,3,4);
plot (ref_u_along_y(:,2), ref_u_along_y(:,1), "o", "LineWidth", 2);
hold on;
plot (x, u_profile, "LineWidth", 2);
title('u along y (centerline)');
xlabel('y');
ylabel('u');
xlim([0 L]);
axis square;
grid on;

subplot (2,3,5);
sx = 0:2*h:L;
sy = 0:2*h:L;
[X, Y] = meshgrid (x, x);
streamline (X, Y, up', vp', sx, sy);
title('streamlines');
xlabel('x');
ylabel('y');
xlim([0 L]);
axis square;
grid on;

subplot (2,3,6);
contour(x, x, Tp', 20, 'black', "LineWidth", 1.2);
title('isolines');
xlabel('x');
ylabel('y');
xlim([0 L]);
axis square;
grid on;

