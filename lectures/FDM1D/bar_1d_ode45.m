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
%  Code: 1D advection-diffusion-reaction by the FD methos                 %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% Global variables (meaning reported below)
global u gamma alpha beta N h


%-------------------------------------------------------------------------%
% Bar with fixed-temperature boundaries
%-------------------------------------------------------------------------%
L = 0.10;       % length (m)
D = 0.01;       % diameter (m)
rho = 7750;     % density (kg/m3)
Cp = 466;       % constant pressure specific heat (J/kg/K)
lambda = 45.;   % thermal conductivity (W/m/K)
u = 0;          % velocity (m/s)
Tleft = 373;    % hot-boundary (left) temperature (K)
Tright = 273;   % cold-boundary (right) temperature (K)
Tex = 298;      % external temperature (K)
U = 10;         % heat exchange coefficient (W/m2/K)
Tin = 273;      % initial temperature (K)
N = 50;         % number of grid points (-)
tf = 600;       % total time (s)


%-------------------------------------------------------------------------%
% Preprocessing
%-------------------------------------------------------------------------%
gamma = lambda/rho/Cp;      % thermal diffusivity (m2/s)
beta = -U/(rho*Cp)*(4/D);   % source term (1/s)
alpha = -beta*Tex;          % source term (K/s)
h = L/(N-1);                % grid spacing (m)


%-------------------------------------------------------------------------%
% Solution via ode45 solver
%-------------------------------------------------------------------------%
T = [Tleft; ones(N-2,1)*Tin; Tright];
[t, T] = ode45(@ODESystem, 0:1:tf, T);


%-------------------------------------------------------------------------%
% Video setup
%-------------------------------------------------------------------------%
video_name = 'bar.mp4';
videompg4 = VideoWriter(video_name, 'MPEG-4');
open(videompg4);

for k=1:size(T,1)
    hold off;
    plot(0:h:L, T(k,:), 'linewidth',2);
    hold on;
    xlabel('axial length [m]');
    ylabel('temperature');    
    message = sprintf('time=%f', t(k));
    time = annotation('textbox',[0.15 0.8 0.1 0.1],'String',message,'EdgeColor','none');
    frame = getframe(gcf);
    writeVideo(videompg4,frame);
    delete(time);
end

% Closing the video stream
close(videompg4);


%-------------------------------------------------------------------------%
% ODE system
%-------------------------------------------------------------------------%
function dT_over_dt = ODESystem(~,T)

    global u gamma alpha beta N h

    dT_over_dt = zeros(N,1);

    % Boundary @ x=0
    dT_over_dt(1) = 0.;

    % Internal points
    for i=2:N-1

        dT_over_dx = (T(i+1)-T(i-1))/(2*h);
        d2T_over_dx2 = (T(i+1)-2.*T(i)+T(i-1))/h^2;

        dT_over_dt(i) = -u*dT_over_dx + ...
                         gamma*d2T_over_dx2 + ...
                         alpha + beta*T(i);

    end

    % Boundary @ x=L
    dT_over_dt(N) = 0.;

end
