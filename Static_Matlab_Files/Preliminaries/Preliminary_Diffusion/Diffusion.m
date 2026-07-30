%% Spatial ODE with space varying Wnt Signal
close; clear all; clc

tic 

% Discretize time: 
N_t   = 2e2;
t_0   = 0;
t_end = 2;
t     = linspace(t_0,t_end,N_t);

% Discretize space: 
N_x   = 1e2; 
x_0   = 0; 
x_end = 1; 
x     = linspace(x_0,x_end,N_x)';
dx    = x(2) - x(1);

% Parameters:
k = .1;

% Differentiation Matrix: 
Dxx            = toeplitz([-2 1 zeros(1,N_x-2)]);
Dxx(1,2)       = 2;
Dxx(end,end-1) = 2;
Dxx            = k*Dxx/dx^2;

% Initial condition: 
f = @(x) 2*x.*cos(2*pi*x)+1;

U0 = f(x);

% Define the Right hand side: 
dUdt = @(t,U) Dxx*U + .09*U.*(1-U/1.6);

% Solve the Heat equation:
tic
options = odeset('Stats','on');
[t,U] = ode23s(dUdt,t,U0,options);
toc

% Animation: 
for i = 1:N_t
   
    figure(1)
    plot(x,U(i,:),'k','linewidth',3); hold off; 
    if i == 1
        pause
    end
    grid on
    grid minor
    ylim([0 2])

end


% Surface Plot: 
[T,X] = meshgrid(t,x);

figure(2)
surf(T,X,U')
shading interp
hold on

for i = 1:6:N_x

    plot3(t,x(i)*ones(size(t)),U(:,i),'k','linewidth',2)
    if i == 1
        plot3(zeros(size(x)),x,U(1,:),'k','linewidth',5)
    end
    hold on
end