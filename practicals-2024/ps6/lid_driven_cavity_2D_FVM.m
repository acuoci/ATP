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
%  Code: Solution of 2D incompressible Navier-Stokes equations            %
%        using a staggered FVM, and time-explicit discretization.         %
%                                                                         %
% ----------------------------------------------------------------------- %

clc; close all; clear;

% Pre-Processing
L = 1.;             % Length of the domain [m]
nu = 1.e-2;         % Kinematic viscosity [m2/s]
tau = 20.;          % Total simulation time [s]
level = 5;          % Maximum level of refinement

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
unwall = 1.;            % north side velocity [m/s]
uswall = 0.;
vewall = 0.;
vwwall = 0.;

% Poisson solver settings
maxiter = 10000;
beta = 1.9;
tolerance = 1.e-6;

% Time step setup
sigma = 0.5;
dt_diff = h^2/4/nu;
dt_conv = 4*nu/unwall^2;
dt = sigma*min (dt_diff, dt_conv);
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
end

%-- Linear interpolations
for i=1:nx+1
    for j=1:ny+1
        pp(i,j) = 0.25*(p(i,j) + p(i+1,j) + p(i,j+1) + p(i+1,j+1));
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

subplot (2,3,1);
surf (x, x, up');
view (2);
title('u.x');
axis square;
colormap jet;
colorbar;
shading interp;

subplot (2,3,2);
surf (x, x, vp');
view (2);
title('u.y');
axis square;
colormap jet;
colorbar;
shading interp;

subplot (2,3,3);
surf (x, x, pp');
view (2);
title('p');
axis square;
colormap jet;
colorbar;
shading interp;

%-- Compare with exp data (available only for Re=100, 400, and 1000)

% Read experimental data from file
exp_u_along_y = readmatrix ("exp_u_along_y");
exp_v_along_x = readmatrix ("exp_v_along_x");

u_profile = up(nx/2+1,:);
v_profile = vp(:,ny/2+1);
p_vertical_profile = pp(nx/2+1,:);
p_horizontal_profile = pp(:,ny/2+1);

% Be careful: cols 1,2 for Re=100, 3,4 for Re=400, 5,6 for Re=1000
subplot (2,3,4);
plot (exp_u_along_y(:,1), exp_u_along_y(:,2), "o", "LineWidth", 2);
hold on;
plot (x, u_profile, "LineWidth", 2);
title('u along y (centerline)');
legend ("Experiments", "Results");
xlabel('y');
ylabel('u');
xlim([0 L]);
axis square;
grid on;

subplot (2,3,5);
plot (exp_v_along_x(:,1), exp_v_along_x(:,2), "o", "LineWidth", 2);
hold on;
plot (x, v_profile, "LineWidth", 2);
title('v along x (centerline)');
legend ("Experiments", "Results");
xlabel('y');
ylabel('v');
xlim([0 L]);
axis square;
grid on;

subplot (2,3,6);
sx = 0:4*h:L;
sy = 0:4*h:L;
[X, Y] = meshgrid (x, x);
streamline (X, Y, up', vp', sx, sy);
title('streamlines');
xlabel('x');
ylabel('y');
xlim([0 L]);
axis square;
grid on;

if (0)
    set(gcf, 'Color', 'w');

    x0=10;
    y0=10;
    width=1800;
    height=1200;
    set(gcf,'position',[x0,y0,width,height]);

    img = getframe(gcf);
    imwrite(img.cdata, 'lid.png');
end
