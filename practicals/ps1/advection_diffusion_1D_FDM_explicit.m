
% 1D Advection-Diffusionn Equations:
% using FDM and explicit time discretization

clc; close all; clear;

%-------------------------------------------------------------------------%
% User-defined data
%-------------------------------------------------------------------------%

np = 21;             % number of grid points    
dt = 0.05;          % time step [s]
nstep = 100;         % number of time steps
%tottime = 0.5;        % total simulation seconds  
%nstep = tottime/dt;   % number of time steps
L = 2.0;              % domain length [m]
u = 1;                % velocity [m/s]
D = 0.05;             % diffusion coefficient [m2/s]
A = 0.5;              % amplitude of initial solution
k = 1;                % wave number [1/m]

%-------------------------------------------------------------------------%
% Pre-processing of user-defined data
%-------------------------------------------------------------------------%

% Grid step calculation
h = L/(np-1);         % grid step [m]
x = linspace(0,L,np);

% Memory allocation
fo = zeros(np,1);     % temporary numerical solution
f  = zeros(np,1);      % current numerical solution
a  = zeros(np,1);      % exact solution

% Initial solution
for i=1:np
	f(i)=A*sin(2*pi*k*h*(i-1)); 
end

% Check the stability conditions on time step
Co = u*dt/h;                        % Courant number
Di = D*dt/h^2;                      % Diffusion number
dt_max = min(1*h/u, 0.5*h*h/D);     % Maximum allowed time step
fprintf('Co=%f, Di=%f, dt=%f, dt(max)=%f\n', Co, Di, dt, dt_max);

%-------------------------------------------------------------------------%
% Video setup
%-------------------------------------------------------------------------%
% video_name = 'advection_diffusion_1d.mp4';
% videompg4 = VideoWriter(video_name, 'MPEG-4');
% open(videompg4);

%-------------------------------------------------------------------------%
% Advancing in time
%-------------------------------------------------------------------------%
t = 0.;
for m=1:nstep
    
    % Update the analytical solution
    for i=1:np 
		% a(i) = A*exp(-4*pi*pi*k*k*D*t)*sin(2*pi*k*(h*(i-1)-u*t)); 
        a(i) = A*exp(-4*pi*pi*k*k*D*t)*sin(2*pi*k*(x(i)-u*t)); 
    end  
    
    % Squared areas below the analytical and numerical solutions
    a2_int = 0.;
    f2_int = 0.;
    for i=1:np-1
         a2_int = a2_int + h/2*(a(i)^2+a(i+1)^2);
         f2_int = f2_int + h/2*(f(i)^2+f(i+1)^2);
    end  
    
    % Graphical output
    message = sprintf('time=%d\na^2(int)=%d\ny^2(int)=%d', t, a2_int, f2_int);
	hold off; plot(0:h:L,f,'linewidth',2); axis([0 L -1, 1]); % plot num. 
	hold on; plot(0:h:L,a,'r','linewidth',2);                 % plot exact
    hold on; legend('numerical', 'exact');                    % legend
    xlabel('spatial coordinate [m]');
    ylabel('solution');    
    time = annotation('textbox',[0.15 0.8 0.1 0.1],'String',message,'EdgeColor','none');
    frame = getframe(gcf);
    % writeVideo(videompg4,frame);
    delete(time);
    
    % Forward Euler method
    fo=f;   
    for i=2:np-1
		f(i) = fo(i)-(u*dt/2/h)*(fo(i+1)-fo(i-1))+...      % advection
			   D*(dt/h^2)*(fo(i+1)-2*fo(i)+fo(i-1));       % diffusion
    end 
    
    % Periodic boundary condition
    f(np) = fo(np)-(u*dt/2/h)*(fo(2)-fo(np-1))+...
            D*(dt/h^2)*(fo(2)-2*fo(np)+fo(np-1)); 
    f(1)  = f(np);
    
    % Update the error between numerical and analytical solution
    E = 0;
    for i=1:np 
        E = E + (f(i)-a(i))^2;
    end
    E = h*sqrt(E);

    % New time step
    t=t+dt;
       
    % Print the current time (every 25 steps)
    if (mod(m,25)==1), fprintf('time=%d E=%e\n', t, E); end
end

fprintf('time=%d E=%e\n', t, E);


% Closing the video stream
% close(videompg4);
