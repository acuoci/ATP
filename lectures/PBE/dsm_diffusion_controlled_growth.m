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
%  Code: Population Balance Equation (PBE) with growth term only          %
%        Solved using the Discrete Sectional Method                       %
%        R. McGraw (1997) Description of Aerosol Dynamics by the          %
%        Quadrature Method of Moments, Aerosol Science and Technology,    %
%        27:2, 255-265 (1997), DOI: 10.1080/02786829708965471             %
%                                                                         %
% ----------------------------------------------------------------------- %

close all;
clear variables;

% Initial distribution: f=a*r^2*exp(-b*r)
% r = particle radius (mum)
% a and b: distribution parameters
a = 0.108;  % (1/mum/cm3)
b = 0.60;   % (1/mum)

% Growth rate: phi(r) = kG/r
kG = 0.78;  % growth rate constant (mum2/s)

% Domain of integration
rMax = 30;  % maximum radius (mum)
M = 500;    % number of intervals
tf = 20;    % maximum time (s)

% Initial distribution (#/cm3/mum)
r = 0:rMax/M:rMax;
fIn = fInitial(r,a,b);

% Number of particles (per unit of volume) in each interval (#/cm3)
NIn = zeros(M,1);
for i=1:M
    NIn(i) = fIn(i+1)*(r(i+1)-r(i));
end

% Solution of discrete section equations
[t, N] = ode15s(@ODESystem, 0:0.25:tf, NIn, [], r, kG);
ntimes = length(t);

% Reconstruction of density function (#/cm3/mum)
f = zeros(ntimes, M+1);
for i=1:ntimes
    f(i, 2:end) = N(i,:)./(r(2:end)-r(1:end-1));
end


% Plot solution
figure; hold on;
plot(r, fIn, 'r');
plot(r, f(end,:), 'b');
plot(r, fAnalytical(r, tf, kG, a, b), 'g');
xlabel('r (micron)');
ylabel('f (#/micron/cm^3)');
legend('initial', 'numerical', 'analytical');
hold off;


% Calculation of moments
m0(1) = Moment(r, fIn, 0);
m0(2) = Moment(r, fIn, 1);
mf(1) = Moment(r, f(end,:), 0);
mf(2) = Moment(r, f(end,:), 1);

fprintf('Time=0 s:  Ntot(#/cm3)=%f Rm(micron)=%f\n', m0(1), m0(2));
fprintf('Time=20 s: Ntot(#/cm3)=%f Rm(micron)=%f\n', mf(1), mf(2));


% Evolution of mean radius

%  Analytical results
muAnalytical1 = zeros(length(t),1);
muAnalytical2 = zeros(length(t),1);
for k=1:length(t)
    muAnalytical1(k) = AnalyticalMoments(1, r, t(k), kG, a, b);
    muAnalytical2(k) = AnalyticalMoments(2, r, t(k), kG, a, b);
end

% Numerical results
m1 = zeros(length(t),1);
m2 = zeros(length(t),1);
for i=1:length(t)
    m1(i) = Moment(r, f(i,:), 1);
    m2(i) = Moment(r, f(i,:), 2);
    
end

figure; hold on;
yyaxis left; ylabel('mean radius (micron)');
plot(t, m1); 
plot(t, muAnalytical1, 'b--');
yyaxis right; ylabel('std deviation (micron)');
plot(t, sqrt(m2-m1.^2)); 
plot(t, sqrt(muAnalytical2-muAnalytical1.^2), 'r--');
xlabel('time (s)'); hold off;
legend('mean radius', 'mean radius (analytical)', 'std deviation', 'std deviation (analytical)');


%% Dynamic evolution of density function
video_name = 'dsm.mp4';
videompg4 = VideoWriter(video_name, 'MPEG-4');
open(videompg4);

figure;
for k=1:length(t)
     hold off;
     plot(r, fAnalytical(r, t(k), kG, a, b), 'g');
     hold on;
     plot(r, f(k,:), 'b');
     hold on;
     xlabel('r (\mum)'); ylabel('f (#/micron/cm3'); title('time=20 s'); 
     xlim([0 20]); ylim([0 0.6]);
     legend('numerical','analytical');
     frame = getframe(gcf);
     writeVideo(videompg4, frame);
end
close(videompg4);



%% Equations 
function dN = ODESystem(~, N, r, kG)

    M = length(N);

    f = zeros(M+1,1);
    for i=1:M
        f(i+1) = N(i)/(r(i+1)-r(i));
    end

    dN = zeros(M,1);
    dN(1) = -f(2)*Ldot(kG, r(2));
    for i=2:M-1
        dN(i) = -f(i+1)*Ldot(kG, r(i+1)) + f(i)*Ldot(kG, r(i));
    end
    dN(M) = f(M)*Ldot(kG, r(M));

end


%% Evaluation of moments
function m = Moment(r, f, order)

    np = length(r);
    m = 0;
    for i=1:np-1
        deltar = r(i+1)-r(i);
        I = 0.50*(f(i+1)*r(i+1)^order+f(i)*r(i)^order);
        m = m + deltar*I;
    end

end


%% Growth rate function
function Lprime = Ldot(kG, r)

    Lprime = kG/r;

end


%% Initial solution
function f = fInitial(r, a, b)

    f = a*(r.^2).*exp(-b*r);

end


%% Analytical solution
function f = fAnalytical(r, t, kG, a, b)

    f = zeros(length(r),1);
    for i=1:length(r)
        arg = r(i)^2-2*kG*t;
        if (arg > 0)
            f(i) = r(i)/sqrt(arg).*fInitial(sqrt(arg), a, b);
        end
    end

end


%% Moments of the analytical solution
function muk = AnalyticalMoments(k, r, t, kG, a, b)

    muk = 0.;
    for j=1:length(r)-1
        deltar = r(j+1)-r(j);
        Ic = 0.50*(r(j+1)^(k)*fAnalytical(r(j+1), t, kG, a, b)+r(j)^(k)*fAnalytical(r(j), t, kG, a, b));
        muk = muk + Ic*deltar;
    end

end