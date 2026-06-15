%% Diffusion 2D Nonlinear
close all; clear; clc

% Discretize time: 
N_t = 2e2;
t_0 = 0;
t_end = 5;
t = linspace(t_0,t_end,N_t);

% Discretize x space:
N_x = 3e1;
x_0 = 0;
x_end = 1;
x = linspace(x_0, x_end, N_x)';
dx = x(2) - x(1);

% Discretize y space:
N_y = 3e1;
y_0 = 0;
y_end = 1;
y = linspace(y_0, y_end, N_y)';
dy = y(2) - y(1);

%Create mesh grid:

[Y,X] = meshgrid(y,x);

% Parameters

alpha = 0.005;

% Differentiation Matrix
Dx = diag(1/2*ones(N_x-1,1),1)+diag(-1/2*ones(N_x-1,1),-1);
Dx(1,2) = 0;
Dx(end, end-1) = 0;
Dx = Dx/dx;

% Differentiation Matrix:
Dy = diag(1/2*ones(N_y-1,1),1)+diag(-1/2*ones(N_y-1,1),-1);
Dy(1,2) = 0;
Dy(end, end-1) = 0;
Dy = Dy/dy;

% Initial Condition
%f = @(x,y) ones(size(x))+2*ones(size(x)).*(x>.25).*(x<.75).*(y>.25).*(y<.75);
f = @(x,y) 2*exp(-(x-.5).^2*5 - (y - .5).^2*5);
U0 = f(X,Y);
U0 = U0(:);

% Define the right hand side
dUdt = @(t, U) nonlin_diff(t, U, N_x, N_y, Dx, Dy, alpha);

% Solve the Heat equation:
tic
options = odeset('Stats','on');
[t, U] = ode23s(dUdt, t, U0, options);
toc

%% Animation:
for i = 1: N_t

    U_plot = reshape(U(i,:), N_x, N_y);

    figure(1)
    surf(X,Y,U_plot); hold off;
    if i == 1
        pause
    else
        drawnow
    end
    
    
    grid on
    grid minor
end
