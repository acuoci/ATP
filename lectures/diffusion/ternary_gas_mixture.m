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
%   Copyright(C) 2023 Alberto Cuoci                                       %
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
%  Code: Diffusion in an ideal ternary gas mixture                        %
%        Duncan and Toor (1962)                                           %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% Global variables (meaning reported below)
global D L A V1 V2 Ctot

%-------------------------------------------------------------------------%
% Input data
%-------------------------------------------------------------------------%
L = 8.59e-2;        % length (m)
d = 2.08e-3;        % diameter (m)
V1 = 78.63e-6;      % volume (m3)
V2 = 78.63e-6;      % volume (m3)
P = 1e5;            % pressure (Pa)
T = 35+273;         % temperature (K)
tFinal = 20*3600;   % final time (s)

% Initial mole fractions (1=H2, 2=N2, 3=CO2)
x1in = [0    0.50 0.50];
x2in = [0.50 0.50 0.  ];

% Binary diffusion coefficients
D = [ 0       8.33e-5 6.80e-5; ...
      8.33e-5 0       1.68e-5; ...
      6.80e-5 1.68e-5 0        ];

% Experimental data
expdata = dlmread("ExpDuncanToor.out",'', 1);

%-------------------------------------------------------------------------%
% Pre-processing
%-------------------------------------------------------------------------%

% Cross section area (m2)
A = pi*d^2/4;

% Total concentration (kmol/m3)
Ctot = P/8314./T;


%-------------------------------------------------------------------------%
% Classical Fick's model (without flux correction)
%-------------------------------------------------------------------------%

[tFick,xFick]=ode45(@FickDiffusion,[0 tFinal], [x1in x2in]');  

xFick1 = xFick(:,1:3);      % mole fractions in vessel 1
xFick2 = xFick(:,4:6);      % mole fractions in vessel 2

PlotSingleSpeciesWithExp(tFick,xFick1(:,1),xFick2(:,1),expdata(:,1),expdata(:,2),expdata(:,3),"Hydrogen H_2 (Fick's model)");
PlotSingleSpeciesWithExp(tFick,xFick1(:,2),xFick2(:,2),expdata(:,1),expdata(:,4),expdata(:,5),"Nitrogen N_2 (Fick's model)");
PlotSingleSpeciesWithExp(tFick,xFick1(:,3),xFick2(:,3),expdata(:,1),expdata(:,6),expdata(:,7),"Carbon dioxide CO_2 (Fick's model)");
PlotSingleVessel(tFick,xFick1,"Vessel 1 (Fick's model)");
PlotSingleVessel(tFick,xFick2,"Vessel 2 (Fick's model)");


%-------------------------------------------------------------------------%
% Corrected Fick's model (i.e., including the diffusion flux correction)
%-------------------------------------------------------------------------%

[tFickCor,xFickCor]=ode45(@FickCorrectedDiffusion,[0 tFinal], [x1in x2in]');  

xFickCor1 = xFickCor(:,1:3);      % mole fractions in vessel 1
xFickCor2 = xFickCor(:,4:6);      % mole fractions in vessel 2

PlotSingleSpeciesWithExp(tFickCor,xFickCor1(:,1),xFickCor2(:,1),expdata(:,1),expdata(:,2),expdata(:,3),"Hydrogen H_2 (Corrected Fick's model)");
PlotSingleSpeciesWithExp(tFickCor,xFickCor1(:,2),xFickCor2(:,2),expdata(:,1),expdata(:,4),expdata(:,5),"Nitrogen N_2 (Corrected Fick's model)");
PlotSingleSpeciesWithExp(tFickCor,xFickCor1(:,3),xFickCor2(:,3),expdata(:,1),expdata(:,6),expdata(:,7),"Carbon dioxide CO_2 (Corrected Fick's model)");
PlotSingleVessel(tFickCor,xFickCor1,"Vessel 1 (Corrected Fick's model)");
PlotSingleVessel(tFickCor,xFickCor2,"Vessel 2 (Corrected Fick's model)");

%-------------------------------------------------------------------------%
% Maxwell-Stefan model (i.e., including the diffusion flux correction)
%-------------------------------------------------------------------------%

[tMS,xMS]=ode45(@MaxwellStefanDiffusion,[0 tFinal], [x1in x2in]');  

xMS1 = xMS(:,1:3);      % mole fractions in vessel 1
xMS2 = xMS(:,4:6);      % mole fractions in vessel 2

PlotSingleSpeciesWithExp(tMS,xMS1(:,1),xMS2(:,1),expdata(:,1),expdata(:,2),expdata(:,3),"Hydrogen H_2 (Maxwell-Stefan's model)");
PlotSingleSpeciesWithExp(tMS,xMS1(:,2),xMS2(:,2),expdata(:,1),expdata(:,4),expdata(:,5),"Nitrogen N_2 (Maxwell-Stefan's model)");
PlotSingleSpeciesWithExp(tMS,xMS1(:,3),xMS2(:,3),expdata(:,1),expdata(:,6),expdata(:,7),"Carbon dioxide CO_2 (Maxwell-Stefan's model)");
PlotSingleVessel(tMS,xMS1,"Vessel 1 (Maxwell-Stefan's model)");
PlotSingleVessel(tMS,xMS2,"Vessel 2 (Maxwell-Stefan's model)");

%-------------------------------------------------------------------------%
% Generalized Fick's model
%-------------------------------------------------------------------------%

[tFickGen,xFickGen]=ode45(@FickGeneralizedDiffusion,[0 tFinal], [x1in x2in]');  

xFickGen1 = xFickGen(:,1:3);      % mole fractions in vessel 1
xFickGen2 = xFickGen(:,4:6);      % mole fractions in vessel 2

PlotSingleSpeciesWithExp(tFickGen,xFickGen1(:,1),xFickGen2(:,1),expdata(:,1),expdata(:,2),expdata(:,3),"Hydrogen H_2 (Generalized Fick's model)");
PlotSingleSpeciesWithExp(tFickGen,xFickGen1(:,2),xFickGen2(:,2),expdata(:,1),expdata(:,4),expdata(:,5),"Nitrogen N_2 (Generalized Fick's model)");
PlotSingleSpeciesWithExp(tFickGen,xFickGen1(:,3),xFickGen2(:,3),expdata(:,1),expdata(:,6),expdata(:,7),"Carbon dioxide CO_2 (Generalized Fick's model)");
PlotSingleVessel(tFickGen,xFickGen1,"Vessel 1 (Generalized Fick's model)");
PlotSingleVessel(tFickGen,xFickGen2,"Vessel 2 (Generalized Fick's model)");


%-------------------------------------------------------------------------%
% Generalized Fick's model
%-------------------------------------------------------------------------%

PlotComparisonSingleSpeciesWithExp(tFick,tFickCor,tMS, xFick1(:,1),xFickCor1(:,1),xMS1(:,1), ...
                    xFick2(:,1),xFickCor2(:,1),xMS2(:,1), expdata(:,1),expdata(:,2),expdata(:,3),"Hydrogen H_2");

PlotComparisonSingleSpeciesWithExp(tFick,tFickCor,tMS, xFick1(:,2),xFickCor1(:,2),xMS1(:,2), ...
                    xFick2(:,2),xFickCor2(:,2),xMS2(:,2), expdata(:,1),expdata(:,4),expdata(:,5),"Nitrogen N_2");

PlotComparisonSingleSpeciesWithExp(tFick,tFickCor,tMS, xFick1(:,3),xFickCor1(:,3),xMS1(:,3), ...
                    xFick2(:,3),xFickCor2(:,3),xMS2(:,3), expdata(:,1),expdata(:,6),expdata(:,7),"Carbon dioxide CO_2");


%-------------------------------------------------------------------------%
% Classical Fick's model (without flux correction)
%-------------------------------------------------------------------------%

function dx_over_dt = FickDiffusion(~,x)

    global D L A V1 V2 Ctot

    x1 = x(1:3);                % mole fractions in vessel 1
    x2 = x(4:6);                % mole fractions in vessel 2

    % Average mole fractions
    X(1) = 0.50*(x1(1)+x2(1));    % average mole fraction of H2
    X(2) = 0.50*(x1(2)+x2(2));    % average mole fraction of N2
    X(3) = 0.50*(x1(3)+x2(3));    % average mole fraction of CO2

    Dmix = FickDiffusionCoefficients(X,D);

    J = -Ctot*Dmix.*(x2-x1)/L;  % diffusion flux from 1 to 2 (kmol/m2/s)

    dx1_over_dt = -J/Ctot*A/V1; % governing equations for vessel 1 
    dx2_over_dt =  J/Ctot*A/V2; % governing equations for vessel 2

    dx_over_dt = [dx1_over_dt ; dx2_over_dt];

end

%-------------------------------------------------------------------------%
% Corrected Fick's model (without flux correction)
%-------------------------------------------------------------------------%

function dx_over_dt = FickCorrectedDiffusion(~,x)

    global D L A V1 V2 Ctot

    x1 = x(1:3);                % mole fractions in vessel 1
    x2 = x(4:6);                % mole fractions in vessel 2

    X(1) = 0.50*(x1(1)+x2(1));    % average mole fraction of H2
    X(2) = 0.50*(x1(2)+x2(2));    % average mole fraction of N2
    X(3) = 0.50*(x1(3)+x2(3));    % average mole fraction of CO2

    Dmix = FickDiffusionCoefficients(X,D);

    J = -Ctot*Dmix.*(x2-x1)/L;  % diffusion flux from 1 to 2 (kmol/m2/s)

    % Correction flux
    Jc = -sum(J);
    J = J+0.50*(x1+x2)*Jc;

    dx1_over_dt = -J/Ctot*A/V1; % governing equations for vessel 1 
    dx2_over_dt =  J/Ctot*A/V2; %

    dx_over_dt = [dx1_over_dt ; dx2_over_dt];

end

%-------------------------------------------------------------------------%
% Maxwell-Stefan's model 
%-------------------------------------------------------------------------%

function dx_over_dt = MaxwellStefanDiffusion(~,x)

    global D L A V1 V2 Ctot

    x1 = x(1:3);        % mole fractions in vessel 1
    x2 = x(4:6);        % mole fractions in vessel 2

    X(1) = 0.50*(x1(1)+x2(1));
    X(2) = 0.50*(x1(2)+x2(2));
    X(3) = 0.50*(x1(3)+x2(3));

    dx_over_dz = (x2-x1)/L;  % mole fraction gradient 

    N = length(X);
    M = zeros(N,N);
    b = zeros(N,1);

    % Assembling matrix M
    M(N,:) = 1;     % component N as reference
    for i=1:N-1
        for j=1:N
            if (i~=j)   % off-diagonal terms
                M(i,j) = -X(i)/Ctot/D(i,j);
            else        % diagonal terms
                for k=1:N
                    if (k~=i)
                        M(i,j) = M(i,j) + X(k)/Ctot/D(i,k);
                    end
                end
            end
        end
    end

    b(1) = -dx_over_dz(1);
    b(2) = -dx_over_dz(2);
    b(3) = 0;

    % Flux inversion
    J = M\b;

    dx1_over_dt = -J/Ctot*A/V1; % governing equations for vessel 1 
    dx2_over_dt =  J/Ctot*A/V2; % governing equations for vessel 2

    dx_over_dt = [dx1_over_dt ; dx2_over_dt];

end

%-------------------------------------------------------------------------%
% Generalized Fick's model
%-------------------------------------------------------------------------%

function dx_over_dt = FickGeneralizedDiffusion(~,x)

    global D L A V1 V2 Ctot

    x1 = x(1:3);        % mole fractions in vessel 1
    x2 = x(4:6);        % mole fractions in vessel 2

    dx_over_dz = (x2-x1)/L;  % mole fraction gradient 

    % Average mole fractions
    X(1) = 0.50*(x1(1)+x2(1));
    X(2) = 0.50*(x1(2)+x2(2));
    X(3) = 0.50*(x1(3)+x2(3));

    % Assembling the B matrix
    N = length(X);
    B = zeros(N-1,N-1);
    for i=1:N-1
        for j=1:N-1
            if (i~=j)
                B(i,j) = -X(i)*(1/D(i,j)-1/D(i,N));
            else
                B(i,j) = X(i)/D(i,N);
                for k=1:N
                    if (k~=i)
                        B(i,j) = B(i,j) + X(k)/D(i,k);
                    end
                end
            end
        end
    end
    
    % Ideal mixture
    Gamma = diag(ones(N-1,1));

    % Fick diffusion matrix
    Dmix = B\Gamma;
    
    % Fluxes
    J = zeros(N,1);
    J(1) = -Ctot*(Dmix(1,1)*dx_over_dz(1)+Dmix(1,2)*dx_over_dz(2));
    J(2) = -Ctot*(Dmix(2,1)*dx_over_dz(1)+Dmix(2,2)*dx_over_dz(2));
    J(3) = -J(1)-J(2);

    % Mass balances
    dx1_over_dt = -J/Ctot*A/V1; % governing equations for vessel 1 
    dx2_over_dt =  J/Ctot*A/V2; %

    dx_over_dt = [dx1_over_dt ; dx2_over_dt];

end


%-------------------------------------------------------------------------%
% Fick's diffusion coefficients in the mixture
%-------------------------------------------------------------------------%
function Dmix = FickDiffusionCoefficients(X,D)

    N = length(X);
    Dmix = zeros(N,1);
    for i=1:N
        sum = 0;
        for k=1:N
            if (i~=k)
                sum = sum+X(k)/D(i,k);
            end
        end
        Dmix(i) = (1-X(i))/sum;
    end

end



function PlotSingleSpeciesWithExp(t,x1,x2,expdatax,expdata1,expdata2,title_fig)

    figure; hold on;
    plot(t/3600, x1,'LineWidth',2);
    plot(t/3600, x2,'LineWidth',2);
    scatter(expdatax, expdata1, 60, 'r');
    scatter(expdatax, expdata2, 60, 'b');
    title(title_fig)
    xlabel('time (hr)'); ylabel('mole fraction');
    legend('Vessel 1', 'Vessel 2');
    set(gca,'fontsize', 14);
    hold off;

end

function PlotSingleVessel(t,x,title_fig)

    figure; hold on;
    plot(t/3600, x(:,1),'LineWidth',2);
    plot(t/3600, x(:,2),'LineWidth',2);
    plot(t/3600, x(:,3),'LineWidth',2);
    plot(t/3600, sum(x,2),'LineWidth',2);
    title(title_fig)
    xlabel('time (hr)'); ylabel('mole fraction');
    legend('H2', 'N2', 'CO2', 'Tot');
    set(gca,'fontsize', 14);
    hold off;

end

function PlotComparisonSingleSpeciesWithExp(ta,tb,tc,x1a,x1b,x1c,x2a,x2b,x2c, expdatax,expdata1,expdata2, title_fig)

    figure; hold on;
    plot(ta/3600, x1a,'b--','LineWidth',2);
    plot(tb/3600, x1b,'b.','LineWidth',2);
    plot(tc/3600, x1c,'b','LineWidth',2);
    plot(ta/3600, x2a,'r--','LineWidth',2);
    plot(tb/3600, x2b,'r.','LineWidth',2);
    plot(tc/3600, x2c,'r','LineWidth',2);
    scatter(expdatax, expdata1, 60, 'r');
    scatter(expdatax, expdata2, 60, 'b');
    title(title_fig)
    xlabel('time (hr)'); ylabel('mole fraction');
    legend('Fick (1)', 'Corrected Fick (1)', 'Maxwell-Stefan (1)', 'Fick (2)', 'Corrected Fick (2)', 'Maxwell-Stefan (2)');
    set(gca,'fontsize', 14);
    hold off;

end
