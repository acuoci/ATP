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
%  Code: 1D advection equation for a discontinuous field using FVM        %
%        testing different advection schemes and flux limiters.           %
%                                                                         %
% ----------------------------------------------------------------------- %

clc; close all; clear;

%-------------------------------------------------------------------------%
% User-defined data
%-------------------------------------------------------------------------%

global UNDEFINED;

ncells = 200;           % number of cells
L = 5.0;                % domain length [m]
u = 1;                  % velocity [m/s]
D = 0.;                 % diffusion coefficient [m2/s]
tau = 1.5;              % total simulation time [s]
fleft = 2.;             % left side boundary value
fright = 1.;            % right side boundary value
UNDEFINED = 1e30;       % Large value

%-------------------------------------------------------------------------%
% Pre-processing of user-defined data
%-------------------------------------------------------------------------%

% Grid step calculation
h = L/ncells;           % grid step [m]
x = linspace (0, L, ncells+1);
xc = linspace (-0.5*h, L+0.5*h, ncells+2);

% Memory allocation
cds.f  = zeros (ncells+2, 1);
cds.fo = zeros (ncells+2, 1);
cds.fp = zeros (ncells+1, 1);

uds.f  = zeros (ncells+2, 1);
uds.fo = zeros (ncells+2, 1);
uds.fp = zeros (ncells+1, 1);

flm.f  = zeros (ncells+2, 1);
flm.fo = zeros (ncells+2, 1);
flm.fp = zeros (ncells+1, 1);

% Initial solution
for i=1:ncells+2
    if (xc(i) > 0.5*L)
        cds.f(i) = 1.;
        uds.f(i) = 1.;
        flm.f(i) = 1.;
    else
        cds.f(i) = 2.;
        uds.f(i) = 2.;
        flm.f(i) = 2.;
    end
end

% Check the stability conditions on time step
dt_max = 1*h/u;
sigma = 0.01;

dt = sigma*dt_max;
nstep = tau/dt;

%-------------------------------------------------------------------------%
% Advancing in time
%-------------------------------------------------------------------------%

t = 0.;
for m=1:nstep

    % Advance time
    t = t + dt;       

    % Interpolate f to nodes
    for i=1:ncells+1
        cds.fp(i) = 0.5*(cds.f(i+1) + cds.f(i));
        uds.fp(i) = 0.5*(uds.f(i+1) + uds.f(i));
        flm.fp(i) = 0.5*(flm.f(i+1) + flm.f(i));
    end

    % Plot results every 100 iterations and at the first iteration
    if (mod(m,100)==0 || m==1)
        hold off;
        plot (x, cds.fp, 'b', "LineWidth", 2); hold on;
        plot (x, uds.fp, 'r', "LineWidth", 2); hold on;
        plot (x, flm.fp, '--black', "LineWidth", 2);
        legend (["centered", "upwind", "flux limiter"]);
        title ('Impact of the Advection Scheme on the transport of a discontinuity');
        xlabel('spatial coordinate [m]');
        ylabel('solution');
        grid on;
        ylim([1 2.4]);
        drawnow;

        % Assembly plots in a gif
        if(0)
            filename = 'ps4-discontinuity.gif';
            message = sprintf('time=%g s\n', t);
            time = annotation ('textbox', [0.2 0.05 0.2 0.2], 'String', message, 'EdgeColor', 'none');
            set(gcf, 'Color', 'w');
            frame = getframe(gcf);
            img = frame2im(frame);
            [imind, cm] = rgb2ind(img, 256);
        
            % Write to the GIF file
            if m == 1
                imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
            else
                imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
            end
        
            delete(time);
        end
    end

    % Store old field
    cds.fo = cds.f;
    uds.fo = uds.f;
    flm.fo = flm.f;

    % Impose the west and east wall values using ghost cells
    cds.f(1) = 2.*fleft - cds.f(2);
    uds.f(1) = 2.*fleft - uds.f(2);
    flm.f(1) = 2.*fleft - flm.f(2);

    cds.f(ncells+2) = 2.*fright - cds.f(ncells+1);
    uds.f(ncells+2) = 2.*fright - uds.f(ncells+1);
    flm.f(ncells+2) = 2.*fright - flm.f(ncells+1);

    % Advance the solution in the internal cells except the last one
    for i=2:ncells+1
        cds.Ai = u/2/h*(cds.fo(i+1) - cds.fo(i-1));             % centered
        uds.Ai = u/h*(uds.fo(i) - uds.fo(i-1));                 % upwind

        % Advection term using flux limiters
        rE = r (flm.fo, i);
        rW = r (flm.fo, i-1);
        psiE = psi_muscl (rE);
        psiW = psi_muscl (rW);
        fE = flm.fo(i) + 0.5*psiE*(flm.fo(i) - flm.fo(i-1));

        if (i-1 == 1)
            fW = 0.5*(flm.fo(i) + flm.fo(i-1));
        else
            fW = flm.fo(i-1) + 0.5*psiW*(flm.fo(i-1) - flm.fo(i-2));
        end

        if (rE == UNDEFINED)
            fE = 0.5*(flm.fo(i+1) + flm.fo(i));
        end

        if (rW == UNDEFINED)
            fW = 0.5*(flm.fo(i) + flm.fo(i-1));
        end
        flm.Ai = u/h*(fE - fW);

        cds.Di = D/h^2*(cds.fo(i+1) + cds.fo(i-1) - 2*cds.fo(i));
        uds.Di = D/h^2*(uds.fo(i+1) + uds.fo(i-1) - 2*uds.fo(i));
        flm.Di = D/h^2*(flm.fo(i+1) + flm.fo(i-1) - 2*flm.fo(i));

        cds.f(i) = cds.fo(i) + dt*(-cds.Ai + cds.Di);
        uds.f(i) = uds.fo(i) + dt*(-uds.Ai + uds.Di);
        flm.f(i) = flm.fo(i) + dt*(-flm.Ai + flm.Di);
    end
end

%-- Helper functions

function ri = r(f, i)
    global UNDEFINED;

    if (i == 1)
        ri = UNDEFINED;
    else
        if ((f(i) - f(i-1)) == 0)
            ri = UNDEFINED;
        else
            ri = (f(i+1) - f(i))/(f(i) - f(i-1));
        end
    end
end

% Flux limiters: see pag. 96 Ferziger-Peric for more

function psi = psi_minmod (r)
    psi = max (0, min (r, 1));
end

function psi = psi_vanleer (r)
    psi = (r + abs(r))/(1 + abs(r));
end

function psi = psi_muscl (r)
    psi = max (0., min ([2, 2*r, 0.5*(1+r)]));
end

function psi = psi_superbee (r)
    psi = max ([0, min(2*r, 1), min(r, 2)]);
end

function psi = psi_upwind (r)
    psi = 0;
end

function psi = psi_centered (r)
    psi = r;
end