clc; close all; clear;

% We create a vector with the number of points at which we want to run the
% simulation.
NP = 10;
while (NP(end) < 10000.)
    NP(end+1) = NP(end)*2;
end

% We allocate two vectors in memory, with dimensions equal to the number of
% points vector.
err_cds_vec = zeros(size(NP));
err_fds_vec = zeros(size(NP));

% We fill the two vectors with the errors, from the run function.
for i=1:length(NP)
    [err_cds, err_fds] = run (NP(i));
    err_cds_vec(i) = err_cds;
    err_fds_vec(i) = err_fds;
end

% The errors are plotted in logarithmic scale, for both x and y axis, and
% compared with the 1^st and 2^nd order trends. It's important to set hold
% on after the first 'loglog' command. Otherwise, the log scale is not
% applied.
figure;
loglog (NP, err_cds_vec, 'x', LineWidth=2);
hold on;
loglog (NP, err_fds_vec, 'x', LineWidth=2);
loglog (NP, 1*NP.^-1, LineWidth=2);
loglog (NP, 0.5*NP.^-2, LineWidth=2);
set (gca, 'YScale', 'log');
xlabel ("N");
ylabel ("L_2");
legend ("CDS", "FDS", "1^{st} order", "2^{nd} Order");
title ("Convergence rate of the space derivative");
axis("square");
grid ("on");

% We write a function that, given a function 'fun', computes the numerical
% spatial derivatives using the FDM, with the CDS and FDS approaches. The
% goal is to compare the errors obtained with the two different methods.
function [err_cds, err_fds] = run(NP)
    np = NP;                        % Total number of points
    x0 = 0;                         % Initial spatial coordinate
    xL = 1;                         % Final spatial coordinate
    h = (xL - x0)/np;               % Grid step size
    x = linspace (x0, xL, np);      % Vector with the point coordinates
    fun = exp(-x.^2);               % Given function
    d_an = -2.*x.*exp(-x.^2);       % Analytic derivative of 'fun'
    
    % Compute the numerical derivative
    % Two vectors with dimensions equal to the size of 'x' are allocated
    d_cds = zeros(size(x));
    d_fds = zeros(size(x));
    
    % loop over the internal points, and compute the derivative
    for i=2:np-1
        d_cds(i) = (fun(i+1) - fun(i-1))/(2.*h);    % CDS scheme
        d_fds(i) = (fun(i+1) - fun(i))/h;           % FDS scheme
    end
    
    % Plot the results
    % figure;
    % hold on;
    % plot (x, fun, LineWidth=1.5);
    % plot (x, d_an, LineWidth=1.5);
    % plot (x, d_cds, LineWidth=1.5);
    % plot (x, d_fds, LineWidth=1.5);
    % legend ("Function", "Derivative", "CDS", "FDS");

    % Compute the errors in each internal points (np-2 points).
    err_cds_vec = zeros (np-2);
    err_fds_vec = zeros (np-2);
    for i=2:np-1
        err_cds_vec(i) = abs (d_cds(i) - d_an(i));
        err_fds_vec(i) = abs (d_fds(i) - d_an(i));
    end

    % Compute the total errors summing the squares of the contributions of
    % each point.
    err_cds = 0.;
    err_fds = 0.;
    for i=2:np-1
        err_cds = err_cds + sqrt(err_cds_vec(i)^2)*h^2;
        err_fds = err_fds + sqrt(err_fds_vec(i)^2)*h;
    end
end

