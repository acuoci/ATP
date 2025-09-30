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
%  Code: 1D advection-diffusion-reaction                                  %
%        Comparison between explicit and implicit Euler methods           %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% User-defined data
%-------------------------------------------------------------------------%
np=21;              % number of grid points                
nstep=100;          % number of time steps
L=1.0;              % domain length [m]
dt=0.025;           % time step [s]
u=0.5;              % velocity [m/s]
D=0.01;             % diffusion coefficient [m2/s]
kappa=1;            % kinetic constant [1/s]
C0=0.;              % initial concentration [kmol/m3]
Cin=1.;             % inlet concentration [kmol/m3]

% Pre-processing of user-defined data
%-------------------------------------------------------------------------%
% Grid step calculation
h=L/(np-1);             % grid step [m]

% Memory allocation
Co_exp=zeros(np,1);     % previous numerical solution (explicit Euler)
C_exp=zeros(np,1);      % current numerical solution (explicit Euler)
Co_imp=zeros(np,1);     % previous numerical solution (implicit Euler)
C_imp=zeros(np,1);      % current numerical solution (implicit Euler)
A=zeros(np,np);         % linear system matrix
b=zeros(np);            % linear system RHS

% Initial solution
C_exp(1)=Cin;
C_imp(1)=Cin;
for i=2:np
	C_exp(i)=C0;
    C_imp(i)=C0;
end

% Check the stability conditions on time step (explicit Euler only)
Co = u*dt/h;                        % Courant number
Di = D*dt/h^2;                      % Diffusion number
dt_max = min(1*h/u, 0.5*h*h/D);     % Maximum allowed time step
fprintf('Co=%f, Di=%f, dt=%f, dt(max)=%f\n', Co, Di, dt, dt_max);

% Video setup
%-------------------------------------------------------------------------%
video_name = 'advection_diffusion_reaction_1d.mp4';
videompg4 = VideoWriter(video_name, 'MPEG-4');
open(videompg4);

% Advancing in time
%-------------------------------------------------------------------------%
t = 0.;
for m=1:nstep
    
    % Graphical output
    message = sprintf('time=%d\n', t); 
    hold off; 
    plot(0:h:L,C_exp,'linewidth',2);
    hold on; 
    plot(0:h:L,C_imp,'linewidth',2);
    axis([0 L 0 1]);
    legend('explicit','implicit');
    xlabel('spatial coordinate [m]');
    ylabel('concentration [kmol/m3]');    
    time = annotation('textbox',[0.15 0.8 0.1 0.1],'String',message,'EdgeColor','none');
    frame = getframe(gcf);
    writeVideo(videompg4,frame);
    delete(time);

    % Current solution
    Co_exp=C_exp;
    Co_imp=C_imp;


    % Explicit (forward) Euler
    %---------------------------------------------------------------------%
    C_exp(1) = Cin;
    for i=2:np-1 
        dC_over_dx = (u/2/h)*(Co_exp(i+1)-Co_exp(i-1));
        d2C_over_dx2 = D/h^2*(Co_exp(i+1)-2*Co_exp(i)+Co_exp(i-1));
        C_exp(i) =  Co_exp(i) - dt*dC_over_dx + ...   % advection
			        dt*d2C_over_dx2 + ...             % diffusion
                    -kappa*dt*C_exp(i);               % reaction
    end
    C_exp(np) = Co_exp(np)-(u*dt/h)*(Co_exp(np)-Co_exp(np-1))+...  % advection
                -kappa*dt*C_exp(np);                               % reaction



    % Implicit (backward) Euler
    %---------------------------------------------------------------------%

    % Matrix coefficients (for internal points only)
    AP = 1. + 2*D*dt/h^2 +kappa*dt;
    AE =  u*dt/2/h - D*dt/h^2;
    AW = -u*dt/2/h - D*dt/h^2;
    
    % Inlet boundary: prescribed inlet concentration
    A(1,1) = 1.;
    b(1) = Cin;

    % Internal points: complete convection-diffusion-reaction equation
    for i=2:np-1 
        A(i,i+1) = AE;
		A(i,i)   = AP;
        A(i,i-1) = AW;
        b(i) = Co_imp(i);
    end 

    % Outlet boundary (negligible diffusion)
    A(np,np-1) = -u*dt/h; 
    A(np,np) = 1. + u*dt/h + kappa*dt; 
    b(np) = Co_imp(np);
    
    % Solve the linear system of equations
    C_imp = A\b;



    % New time step
    t=t+dt; 
    
    % Print the current time (every 25 steps)
    if (mod(m,25)==1), fprintf('time=%d\n', t); end
end

% Closing the video stream
close(videompg4);
