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
%  Code: 1D advection-diffusion equation but neglecting the
%        diffusion term. Study how the different discretization
%        of the convective term changes the solution                      %
%                                                                         %
% ----------------------------------------------------------------------- %

clc; close all; clear all;

% ----------------------------------------------------------------------- %
% Pre-Processing
% ----------------------------------------------------------------------- %

% Data
Gamma = 0;                 % Diffusion coefficient [m2/s]
L = 1;                     % Length of the domain [m]
tau = 7;                   % Total simulation time [s]
ncells = 128;               % Total number of cells along one direction: 64 or 128

% Build Mesh
nx = ncells;                % Number of cells along the x dimension
ny = nx;                    % Number of cells along the y dimension
h = L/nx;                   % Size of a mesh cells (squared cells) hx=hy=h
x = 0:h:L;                  % (POSTPROC) vector with the coordinates of x points/nodes, dimensions = ncells+1
y = 0:h:L;                  % (POSTPROC) vector with the coordinates of y points/nodes, dimensions = ncells+1

% Build Time
%dt_Co = 4*Gamma/0.05^2;      % Maximum time step that accounts for the advection phenomena
%dt_Di = h^2/4/Gamma;         % Maximum time step that accounts for the diffusion phenomena
dt_Co = h/0.5;
sigma = 0.1;                % Safety factor
dt = sigma*dt_Co;

% Print the computed minimum delta t. %f tells a "floating-point" number
% has to be printed. \n goes to the next line
fprintf("Maximum time step for convection = %f\n", dt_Co);
%fprintf("Maximum time step for diffusion  = %f\n", dt_Di);
fprintf("Selected time step               = %f\n",dt);

% Number of time steps to run
nsteps = tau/dt;

% Memory allocations
% Define fields on a finite volume grid, which comprises the number of
% cells and the ghost cells (+2) used to set the boundary conditions
alpha = zeros([nx+2, ny+2]);    % Passive scalar field
u     = zeros([nx+2, ny+2]);    % Velocity along the x direction field
v     = zeros([nx+2, ny+2]);    % Velocity along the y direction field


% Set initial condition for the passive scalar alpha reading the cell 
% values stored in the filed "blob-64.txt" or "blob-128.txt" depending
% on the selected resolution.
% At time zero alpha = 1 in a circle of radius 0.15 at the coordinate (0.5,
% 0.75) of the domain
alpha = readVofi(alpha, nx);
alphacds = alpha;

% Set velocity field
[u, v] = updatevelocity_diagonal(u,v,nx);      % Make a constant vel field directed toward south-east
%[u, v] = updatevelocity_circular(u,v,nx);      % Make a constant circular velocity field

% ----------------------------------------------------------------------- %
% Solution loop
% ----------------------------------------------------------------------- %

time = 0;
for t=1:nsteps

    % Update alpha boundary condition
    alphaWall = 0.;
    alpha(1:nx+2,1)    = 2*alphaWall-alpha(1:nx+2,2);          % south
    alpha(1:nx+2,nx+2) = 2*alphaWall-alpha(1:nx+2,nx+1);       % north
    alpha(1,1:nx+2)    = 2*alphaWall-alpha(2,1:nx+2);          % west
    alpha(nx+2,1:ny+2) = 2*alphaWall-alpha(nx+1,1:ny+2);       % east

    alphaWall = 0.;
    alphacds(1:nx+2,1)    = 2*alphaWall-alphacds(1:nx+2,2);          % south
    alphacds(1:nx+2,nx+2) = 2*alphaWall-alphacds(1:nx+2,nx+1);       % north
    alphacds(1,1:nx+2)    = 2*alphaWall-alphacds(2,1:nx+2);          % west
    alphacds(nx+2,1:ny+2) = 2*alphaWall-alphacds(nx+1,1:ny+2);       % east

    % Explicit forward in time integration of the advection-diffusion eq
    % for the passive scalar alpha
    
    alphao = alpha;     % Store alpha at current time t
    alphaocds = alphacds;

    % Loop over all the internal cells (no ghost cells). 
    % The boundaries where already solved setting the boundary conditions
    for i=2:nx+1
        for j=2:ny+1

            % Advection term, discretized using a 2nd order CDS scheme
            
            % *-- CDS --*
            Aijcds = u(i,j)*(alphaocds(i+1,j)-alphaocds(i-1,j))/2/h + ...
                     v(i,j)*(alphaocds(i,j+1)-alphaocds(i,j-1))/2/h;

            % *-- FDS --*
            %Aij = u(i,j)*(alphao(i+1,j)-alphao(i,j))/h + ...
            %      v(i,j)*(alphao(i,j+1)-alphao(i,j))/h;

            % *-- BDS --*
            Aij = u(i,j)*(alphao(i,j)-alphao(i-1,j))/h + ...
                  v(i,j)*(alphao(i,j)-alphao(i,j-1))/h;

            % Find the passive scalar at time t+1 using an explicit
            % (1st order, forward in time) scheme. This scheme require
            % strict time-step limits for stability. This time-step is
            % computed at the beginning of this code (Build Time Section)
            alpha(i,j) = alphao(i,j) + dt*(-Aij);
            alphacds(i,j) = alphaocds(i,j) + dt*(-Aijcds);


        end
    end

    % On-The-Fly Post Processing
    if (mod(t,25)==1)       % => Every 25 time-steps
        fprintf("Time = %f [s]\n", time);   % Print current simulation time
        % The alpha field has dimensions ncells+2 because it comprises the ghost
        % cells. However, to visualize the results we are not interested in
        % the absolute value of the ghost cells. What we want to visualize is
        % the value of alpha in the domain and at the boundaries.
        % To do so, it is necessary to interpolate the cells values of alpha on
        % the nodes. We can simply use a linear interpolation. In this way
        % we can plot also the boundaries but not the ghost cells
        alphanode = nodeInterpolation(alpha, nx);
        alphanodecds = nodeInterpolation(alphacds, nx);

        subplot(1,2,1);
        surf(x, y, alphanode');     % Plot results
        axis("square");             % Adjust axis proportions
        colormap("jet");            % Set of colors used
        colorbar;                   % Show the colorbar
        caxis([0,1]);              % Fix the color interval in the range [0,1]
        shading interp;             % Shade the colors intead of visualizing the mesh
        %zlim([0,1]);                % Set limits to the z coordinate
        view(2);                    % Visualize the plot in 2 dimensions
        xlabel("x");                % Set a name to the x axis
        ylabel("y");                % Set a name to the y axis
        title("Pure Advection: UPWIND")     % Set title of the plot
                
        subplot(1,2,2);
        surf(x, y, alphanodecds');     % Plot results
        axis("square");             % Adjust axis proportions
        colormap("jet");            % Set of colors used
        colorbar;                   % Show the colorbar
        caxis([0,1]);              % Fix the color interval in the range [0,1]
        shading interp;             % Shade the colors intead of visualizing the mesh
        %zlim([0,1]);                % Set limits to the z coordinate
        view(2);                    % Visualize the plot in 2 dimensions
        xlabel("x");                % Set a name to the x axis
        ylabel("y");                % Set a name to the y axis
        title("Pure Advection: LINEAR")     % Set title of the plot

        drawnow;                    % Show plot (needed for plots inside loops)
    end

    time = time + dt;               % Advance value of the current time
end


% ~~~~~~~~ Functions

function [u,v] = updatevelocity_circular(u,v,ncells)
    % Compute a circular velocity field
    % u = 0.5 - y
    % v = x - 0.5
    step = 1/ncells;
    x = step/2:step:ncells+1;
    y = x;

    for i=1:ncells+2
        for j=1:ncells+2
            u(i,j) = 0.5 - y(j);
            v(i,j) = x(i) - 0.5;
        end
    end

end


function [u,v] = updatevelocity_diagonal(u,v,ncells)
    % Compute constant diagonal velocity field
    % u = +0.02;
    % v = -0.05;
    for i=1:ncells+2
        for j=1:ncells+2
            u(i,j) = 0.05;
            v(i,j) = 0.0;
        end
    end

end


function nodeField = nodeInterpolation(field, ncells)
    % Interpolate a field defined at the center of the cells to the nodes.
    % For post-processing purposes
    npoints = ncells+1;
    nodeField = zeros(npoints);

    for i=1:npoints
        for j=1:npoints
            nodeField(i,j) = 0.25*(field(i,j) + field(i,j+1) + field(i+1,j) + field(i+1,j+1));
        end
    end

end


function initfield = readVofi(field, ncells)
    % Read the content of "blob-64.txt" or "blob-128.txt" file, which contains
    % the internal cell values of the initial condition of the alpha field
    % and store these initial values in the variable initfield
    %namefile = append("blob-",num2str(ncells),".txt");
    namefile = append("centered-",num2str(ncells),".txt");
    fileID = fopen(namefile,'r');
    A = fscanf(fileID,"%f");

    % Store A into initfield
    initfield = field;
    counter = 1;
    for i=2:ncells+1
        for j=2:ncells+1
            initfield(i,j) = A(counter);
            counter = counter+1;
        end
    end
end
