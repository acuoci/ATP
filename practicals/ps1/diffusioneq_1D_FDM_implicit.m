
% 1D Diffusion Equations: using FDM and implicit time discretization
clc; close all; clear;

% ----------------------------------------------------------------------- %
% Pre-Processing
% ----------------------------------------------------------------------- %

% Data
alpha = 0.01;       % Diffusivity coefficient: alpha=k/rho/cp [m2/s]
L = 1;              % length of the domain [m]

% Build Mesh
npoints = 20;       % number of points that discretize the 1D domain
h = L/(npoints-1);  % distance between two consecutive points

% create a vector with the coordinates of the points
% this vector is composed of "npoints" evenly spaced from x=0 to x=L
x = linspace(0,L,npoints);

% Time
tau = 50;                 % total simulation time (high enough to reach steady state) [s]
dt_diff = 0.5*h^2/alpha;  % Maximum time step that accounts for the diffusion phenomena (from Di=0.5)
sigma = 2;                % Safety factor to avoid to work exactly at the minimum stability conditions
dt = sigma*dt_diff;       % Choice of the most limiting time step

% Print the computed minimum delta t. %f tells a "floating-point" number
% has to be printed. \n goes to the next line
fprintf("Maximum time step for diffusion  = %f\n", dt_diff);
fprintf("Selected time step               = %f\n", dt);

nsteps = tau/dt;    % Number of time steps to run

% Initial Conditions & Boundary Conditions
Tleft = 500;        % BC for temperature on the left side of the domain  [K]
Tright = 300;       % BC for temperature on the right side of the domain [K]
Tinit = 300;        % IC for temperature at time=0 [K].

% Memory allocations
T = ones(size(x))*Tinit;  % Create the temperature fields with dimension = number of points

% ----------------------------------------------------------------------- %
% Solution loop
% ----------------------------------------------------------------------- %

% Set known terms (known at time t)
b  = zeros(npoints,1);

% Build diagonals of the global tridiagonal matrix
% Ap = central diagonal
% Ae = upper diagonal
% Aw = lower diagonal
Ap = -(1 + 2*alpha*dt/h^2);
Ae = alpha*dt/h^2;
Aw = alpha*dt/h^2;

% Create global matrix using sparse in order to avoid to store the "0"s
A = sparse(npoints,npoints);

% Fill the global matrix using the three diagonals
A(1,1)=1; A(1,2)=0;                         % From the boundary conditions
for i=2:npoints-1, A(i,i-1) = Aw; end
for i=2:npoints-1, A(i,i)   = Ap; end
for i=2:npoints-1, A(i,i+1) = Ae; end
A(npoints,npoints)=1; A(npoints,npoints-1)=0;   % From the boundary conditions

% loop over all the time-steps
for t=1:nsteps

    To = T;  % Store temperature at time t

    % Update RHS term using info from the old time step
    for i=1:npoints
        b(i) = -To(i);
    end
    b(1) = Tleft;
    b(npoints) = Tright;

    % Solve the linear system of equations
    T = A\b;

    % On-The-Fly Post Processing
    if (mod(t,50)==1)   % => Every 50 time steps
        plot(x,T, "LineWidth", 1.8);      % Plot results
        grid on;                          % Show a grid
        xlabel("length [m]")              % Name of the x axis
        ylabel("Temperature [K]")         % Name of the y axis
        drawnow;                          % Show the plot
    end

end
