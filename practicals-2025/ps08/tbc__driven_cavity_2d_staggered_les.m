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
%   Copyright(C) 2025 Alberto Cuoci                                       %
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
%  Code: 2D driven-cavity problem based on a staggered grid               %
%        The code solves the LES equations for turbulent conditions       %
%        Turbulence is modelled via Smagorinsky model                     %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% ----------------------------------------------------------------------- %
% User data
% ----------------------------------------------------------------------- %

% Only even numbers of cells are acceptable
nx=32;              % number of (physical) cells along x
ny=nx;              % number of (physical) cells along y
L=1;                % length [m]
nu=0.01;            % laminar kinematic viscosity [m2/s] 
tau=20;             % total time of simulation [s]

% Boundary conditions (velocities)
un=100;     % north wall velocity [m/s]
us=0;       % south wall velocity [m/s]
ve=0;       % east wall velocity [m/s]
vw=0;       % west wall velocity [m/s]

% Parameters for SOR
max_iterations=10000;   % maximum number of iterations
beta=1.5;               % SOR coefficient
max_error=1e-3;         % error for convergence

% ----------------------------------------------------------------------- %
% Data processing
% ----------------------------------------------------------------------- %
if (mod(nx,2)~=0 || mod(ny,2)~=0)
    error('Only even number of cells can be accepted (for graphical purposes only)');
end

% Grid step
h=L/nx;             % grid step [m]

% Time step
sigma = 0.5;                        % safety factor for time step (stability)
u2=(un^2+us^2+ve^2+vw^2);           % velocity measure [m2/s2]
dt_diff=h^2/4/nu;                   % time step (diffusion stability)
dt_conv=4*nu/u2;                    % time step (convection stability)
dt=sigma*min(dt_diff, dt_conv);     % time step (stability)
nsteps=tau/dt;                      % number of steps
Re = un*L/nu;                       % Reynolds' number

fprintf('Time step: %f\n', dt);
fprintf(' - Diffusion:  %f\n', dt_diff);
fprintf(' - Convection: %f\n', dt_conv);
fprintf('Reynolds number: %f\n', Re);

% Grid construction
x=0:h:L;                         % grid coordinates (x axis)
y=0:h:L;                         % grid coordinates (y axis)
[X,Y] = meshgrid(x,y);           % MATLAB grid

% ----------------------------------------------------------------------- %
% Memory allocation
% ----------------------------------------------------------------------- %

% Main fields (velocities and pressure)
u=zeros(nx+1,ny+2);
v=zeros(nx+2,ny+1);
p=zeros(nx+2,ny+2);
nut=zeros(nx+2,nx+2)+nu;
kappa=zeros(nx+2,nx+2);

% Temporary velocity fields
ut=zeros(nx+1,ny+2);
vt=zeros(nx+2,ny+1);

% Fields used only for graphical post-processing purposes
uu=zeros(nx+1,ny+1);
vv=zeros(nx+1,ny+1);

% Coefficient for pressure equation
c=zeros(nx+2,ny+2)+1/4;
c(2,3:ny)=1/3;c(nx+1,3:ny)=1/3;c(3:nx,2)=1/3;c(3:nx,ny+1)=1/3;
c(2,2)=1/2;c(2,ny+1)=1/2;c(nx+1,2)=1/2;c(nx+1,ny+1)=1/2;

% ----------------------------------------------------------------------- %
% Solution over time
% ----------------------------------------------------------------------- %
t=0.0;
for is=1:nsteps
    
    % [LES] Update the turbulent viscosity (Smagorinsky model)
    nut = TurbulentViscosity( u, v, nx, ny, h );
                
    % Boundary conditions (velocity)
    u(1:nx+1,1)=2*us-u(1:nx+1,2);
    u(1:nx+1,ny+2)=2*un-u(1:nx+1,ny+1);
    v(1,1:ny+1)=2*vw-v(2,1:ny+1);
    v(nx+2,1:ny+1)=2*ve-v(nx+1,1:ny+1);
    
    % Advection-diffusion equation (predictor)
    [ut, vt] = AdvectionDiffusion2D( ut, vt, u, v, nx, ny, h, dt, nu, nut);
    
    % Pressure equation (Poisson)
    [p, iter] = Poisson2D( p, ut, vt, c, nx, ny, h, dt, ...
                           beta, max_iterations, max_error);
    
    % Correct the velocity
    u(2:nx,2:ny+1)=ut(2:nx,2:ny+1)-(dt/h)*(p(3:nx+1,2:ny+1)-p(2:nx,2:ny+1));
    v(2:nx+1,2:ny)=vt(2:nx+1,2:ny)-(dt/h)*(p(2:nx+1,3:ny+1)-p(2:nx+1,2:ny));
    
    % Print time step on the screen
    if (mod(is,50)==1)
        nut_mean = mean(mean(nut(2:nx+1, 2:ny+1)));
        fprintf('Step: %d - Time: %f - Poisson iterations: %d - Mean nut: %f \n', ...
                is, t, iter, nut_mean);
    end
  
    % Advance time
    t=t+dt;
 
    % ----------------------------------------------------------------------- %
    % Update graphical output
    % ----------------------------------------------------------------------- %
    if (mod(is,200)==1)
        
        % Post-processing (reconstructi on the corners of pressure cells)
        uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
        vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
        nutnut(1:nx+1,1:ny+1)=0.25*(nut(1:nx+1,1:ny+1)+nut(2:nx+2,1:ny+1) + ...
                                    nut(1:nx+1,2:ny+2)+nut(2:nx+2,2:ny+2)); 
        kk(1:nx+1,1:ny+1)=0.25*(kappa(1:nx+1,1:ny+1)+kappa(2:nx+2,1:ny+1) + ...
                                kappa(1:nx+1,2:ny+2)+kappa(2:nx+2,2:ny+2));

        % Vorticity: omega = dv/dx - du/dy
        omega=zeros(nx+1,ny+1);
        for i=2:nx
            for j=2:ny
                dv_dx = (v(i+1,j) - v(i,j))/h;
                du_dy = (u(i,j+1) - u(i,j))/h;
                omega(i,j) = dv_dx - du_dy;
            end
        end
           
        % Contour: x-velocity
        subplot(231);
        contourf(X,Y,uu',20,'LineColor','none');
        axis('square'); title('u velocity (m/s)');
        colormap('jet'); colorbar;

        % Contour: y-velocity
        subplot(232);
        contourf(X,Y,vv',20,'LineColor','none');
        axis('square'); title('v velocity (m/s)');
        colormap('jet'); colorbar;

        % Plot: velocity components along the horizontal middle axis
        subplot(233);
        plot(x,uu(:, round(ny/2)));
        hold on;
        plot(x,vv(:, round(ny/2)));
        axis('square');
        title('velocity along x-axis');
        legend('x-velocity', 'y-velocity');
        hold off;

        % Colormap vorticity
        subplot(234);
        contourf(X,Y,omega',20,'LineColor','none');
        axis('square'); title('vorticity (m/s^2)');
        colormap('jet'); colorbar;
        clim([-max(abs(omega(:))), max(abs(omega(:)))]);  % Symmetric color scale
        
        % Colormap turbulent viscosity 
        subplot(235);
        contourf(X,Y,nutnut',20,'LineColor','none');
        axis('square'); title('turbulent viscosity (m^2/s)');
        colormap('jet'); colorbar;

        % Streamlines
        subplot(236); cla;
        streamslice(X,Y,uu',vv',1);  % The '2' controls density of streamlines
        axis('square'); title('Streamlines');
        xlim([0 L]); ylim([0 L]);
        hold off;

        pause(0.01)
    
    end
    
end


% --------------------------------------------------------------------------------------
% Poisson equation solver
% --------------------------------------------------------------------------------------
function [p, iter] = Poisson2D( p, ut, vt, gamma, nx, ny, h, dt, ...
                                beta, max_iterations, max_error)

    for iter=1:max_iterations
        
        for i=2:nx+1
            for j=2:ny+1
                
                delta = p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1);
                S = (h/dt)*(ut(i,j)-ut(i-1,j)+vt(i,j)-vt(i,j-1));
                p(i,j)=beta*gamma(i,j)*( delta-S )+(1-beta)*p(i,j);
                
            end
        end
        
        % Estimate the error
        epsilon=0.0; 
        for i=2:nx+1
            for j=2:ny+1
                delta = p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1);
                S = (h/dt)*(ut(i,j)-ut(i-1,j)+vt(i,j)-vt(i,j-1));              
                epsilon=epsilon+abs( p(i,j) - gamma(i,j)*( delta-S ) );
            end
        end
        epsilon = epsilon / (nx*ny);
        
        % Check the error
        if (epsilon <= max_error) % stop if converged
            break;
        end 
        
    end

end

% --------------------------------------------------------------------------------------
% Advection-diffusion equation
% --------------------------------------------------------------------------------------
function [ut, vt] = AdvectionDiffusion2D( ut, vt, u, v, nx, ny, h, dt, nu, nut)
                            
    % Temporary u-velocity
    for i=2:nx
        for j=2:ny+1 
            
            nutot = nu + nut(i,j);
            
            ue2 = 0.25*( u(i+1,j)+u(i,j) )^2;
            uw2 = 0.25*( u(i,j)+u(i-1,j) )^2;
            unv = 0.25*( u(i,j+1)+u(i,j) )*( v(i+1,j)+v(i,j) );
            usv = 0.25*( u(i,j)+u(i,j-1) )*( v(i+1,j-1)+v(i,j-1) );
            
            A = (ue2-uw2+unv-usv)/h;
            D = (nutot/h^2)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j));
            
            ut(i,j)=u(i,j)+dt*(-A+D);
            
        end
    end
    
    % Temporary v-velocity
    for i=2:nx+1
        for j=2:ny 
            
            nutot = nu + nut(i,j);
            
            vn2 = 0.25*( v(i,j+1)+v(i,j) )^2;
            vs2 = 0.25*( v(i,j)+v(i,j-1) )^2;
            veu = 0.25*( u(i,j+1)+u(i,j) )*( v(i+1,j)+v(i,j) );
            vwu = 0.25*( u(i-1,j+1)+u(i-1,j) )*( v(i,j)+v(i-1,j) );
            A = (vn2 - vs2 + veu - vwu)/h;
            D = (nutot/h^2)*(v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1)-4*v(i,j));
            
            vt(i,j)=v(i,j)+dt*(-A+D);
            
        end
    end
    
end

% --------------------------------------------------------------------------------------
% [LES] Turbulent viscosity
% --------------------------------------------------------------------------------------
function nut = TurbulentViscosity( u, v, nx, ny, h )

    % Update the turbulent viscosity (Smagorinsky LES model)
    Cs = 0.1;   % Smagorinsky constant (typically 0.1-0.2)
    Delta = h;  % Filter width (grid size)
    
    nut = zeros(nx+2,ny+2);
    for i=2:nx+1
        for j=2:ny+1

            % Calculate strain rate tensor magnitude
            
            ue = u(i,j);
            uw = u(i-1,j);
            vn = v(i,j);
            vs = v(i,j-1);
            
            uen = (u(i,j)+u(i,j+1))/2;
            ues = (u(i,j)+u(i,j-1))/2;
            uwn = (u(i-1,j)+u(i-1,j+1))/2;
            uws = (u(i-1,j)+u(i-1,j-1))/2;
            un = (uen + uwn)/2;
            us = (ues + uws)/2;
            
            vne = (v(i,j)+v(i+1,j))/2;
            vse = (v(i,j-1)+v(i+1,j-1))/2;
            vnw = (v(i,j)+v(i-1,j))/2;
            vsw = (v(i,j-1)+v(i-1,j-1))/2;
            ve = (vne + vse)/2;
            vw = (vnw + vsw)/2;
            
            % Strain rate components
            % [TBC]
            
            % Strain rate tensor magnitude: |S| = sqrt(2*Sij*Sij)
            % [TBC]
            
            % Smagorinsky eddy viscosity
            nut(i,j) = % [TBC]

        end
    end

end