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
%        Solved using the QMOM (McGraw, 1997)                             %
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

% Number of quadrature points
N = 3;

% Domain of integration
rMax = 30;  % maximum radius (mum)
tf = 20;    % maximum time (s)

% Initial distribution (#/cm3/mum)
r = 0:rMax/1000:rMax;
fIn = fInitial(r, a, b);

% Initial moments
Neq = 2*N;  % number of moments/equations
muIn = zeros(Neq,1);
for j=1:length(r)-1
    deltar = r(j+1)-r(j);
    for i=0:Neq-1
        Ic = 0.50*(r(j+1)^(i)*fIn(j+1)+r(j)^(i)*fIn(j));
        muIn(i+1) = muIn(i+1) + Ic*deltar;
    end
end

% Normalization
muNormIn = muIn/muIn(1);

% Solution of equation moments
[t, muNorm] = ode15s(@ODESystem, 0:1:tf, muNormIn, [], kG);

% Denormalization
mu = muNorm*muIn(1);

% Plot temporal evolution of moments
figure;
tiledlayout(2,3);
nexttile; plot(t, muNorm(:,1)); xlabel('time (s)'); ylabel('\mu_0^{norm}');
nexttile; plot(t, muNorm(:,2)); xlabel('time (s)'); ylabel('\mu_1^{norm}');
nexttile; plot(t, muNorm(:,3)); xlabel('time (s)'); ylabel('\mu_2^{norm}');
nexttile; plot(t, muNorm(:,4)); xlabel('time (s)'); ylabel('\mu_3^{norm}');
nexttile; plot(t, muNorm(:,5)); xlabel('time (s)'); ylabel('\mu_4^{norm}');
nexttile; plot(t, muNorm(:,6)); xlabel('time (s)'); ylabel('\mu_5^{norm}');

%% Evolution of mean radius
%  Analytical results
for k=1:length(t)
    muAnalytical1(k) = AnalyticalMoments(1, r, t(k), kG, a, b);
    muAnalytical2(k) = AnalyticalMoments(2, r, t(k), kG, a, b);
end
figure; hold on;
yyaxis left; ylabel('mean radius (micron)');
plot(t, muNorm(:,2)); 
plot(t, muAnalytical1/muIn(1), 'b--');
yyaxis right; ylabel('std deviation (micron)');
plot(t, sqrt(muNorm(:,3)-muNorm(:,2).^2)); 
plot(t, sqrt(muAnalytical2-muAnalytical1.^2), 'r--');
xlabel('time (s)'); hold off;
legend('mean radius', 'mean radius (analytical)', 'std deviation', 'std deviation (analytical)');


% Plot the distribution function
figure;
tiledlayout(1,2);

% Initial time
nexttile; hold on; xlabel('r (\mum)'); ylabel('f (#/micron/cm3'); title('time=0'); xlim([0 20]);
plot(r, fIn, 'r');
[w, L] = MomentInversion(muNorm(1, :));
for i=1:N
    line([L(i),L(i)], [0 w(i)*muIn(1)]);
end
legend('initial', 'numerical');
hold off;

% Final time
nexttile; hold on; xlabel('r (\mum)'); ylabel('f (#/micron/cm3'); title('time=20 s'); xlim([0 20]);
plot(r, fIn, 'r');
[w, L] = MomentInversion(muNorm(end, :));
for i=1:N
    line([L(i),L(i)], [0 w(i)*muIn(1)]);
end
plot(r, fAnalytical(r, tf, kG, a, b), 'g');
legend('initial', 'numerical','analytical');
hold off;


[w, L] = MomentInversion(muNorm(1, :));
fprintf('Time=0 s: rm(micron)=%f sigma(micron)=%f\n', Lm(w,L), StdDev(w,L));

[w, L] = MomentInversion(muNorm(end, :));
fprintf('Time=20 s: rm(micron)=%f sigma(micron)=%f\n', Lm(w,L), StdDev(w,L));


%% Dynamic evolution of density function
video_name = 'qmom.mp4';
videompg4 = VideoWriter(video_name, 'MPEG-4');
open(videompg4);

figure;
for k=1:length(t)
     
     [w, L] = MomentInversion(muNorm(k, :));

     hold off;
     plot(r, fAnalytical(r, t(k), kG, a, b), 'g');
     hold on;
     for i=1:N
         line([L(i),L(i)], [0 w(i)*muIn(1)]);
     end
     hold on;
     xlabel('r (\mum)'); ylabel('f (#/micron/cm3'); title('time=20 s');
     legend('numerical','analytical');
     xlim([0 20]); ylim([0 0.6]);
     frame = getframe(gcf);
     writeVideo(videompg4, frame);
end
close(videompg4);


%% Moment equations
function dmuNorm = ODESystem(~, muNorm, kG)

    Neq = length(muNorm);
    N = Neq/2;

    [w,L] = MomentInversion(muNorm);

    dmuNorm = zeros(Neq, 1);

    for i=1:Neq
        k = i-1;
        sum = 0;
        for j=1:N
            sum = sum + L(j)^(k-1)*Ldot(kG, L(j))*w(j);
        end
        dmuNorm(i) = k*sum;
    end

end


%% Growth rate function
function Lprime = Ldot(kG, r)

    Lprime = kG/r;

end


%% Moment inversion function (PD algorithm, Gordon 1968)
function [w, L] = MomentInversion(muNorm)

    N = length(muNorm)/2;

    P = zeros(2*N+1, 2*N+1);
    P(1,1) = 1;
    P(1,2) = 1;
    for i=2:2*N
        P(i,2) = (-1)^(i-1)*muNorm(i-1+1);
    end
    
    for j=3:2*N+1
        for i=1:2*N+2-j
            P(i,j) = P(1,j-1)*P(i+1,j-2)-P(1,j-2)*P(i+1,j-1);
        end
    end

    alpha = zeros(2*N,1);
    alpha(1) = 0;
    for i=2:2*N
        alpha(i) = P(1,i+1)/P(1,i)/P(1,i-1);
    end
    
    a = zeros(N,1);
    for i=1:N
        a(i) = alpha(2*i)+alpha(2*i-1);
    end
    
    b = zeros(N-1,1);
    for i=1:N-1
        b(i) = sqrt(alpha(2*i+1)*alpha(2*i));
    end

    A = diag(a);
    for i=1:N-1
        A(i,i+1) = b(i);
        A(i+1,i) = b(i);
    end

    [V,csi] = eig(A);
    L = diag(csi);

    w = zeros(N,1);
    for i=1:N
        w(i) = V(1,i)^2;
    end

end


%% Mean value
function m = Lm(w,L)

    m = dot(w,L);

end


%% Standard deviation
function sigma = StdDev(w,L)

    m = Lm(w,L);
    sigma = sqrt( dot(w, L.^2) - m^2);

end

%% Initial density function
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
