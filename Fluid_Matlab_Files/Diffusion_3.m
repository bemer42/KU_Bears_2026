%% Diffusion 3
% Discretize time: 
N_t = 1e3;
t_0 = 0;
t_end = 5;
t = linspace(t_0,t_end,N_t);

% Discretize space:
N_x = 100;
x_0 = 0;
x_end = 1;
x = linspace(x_0, x_end, N_x)';
dx = x(2) - x(1);

% Parameters
k = .1;
alpha = 0.25;
m = 5;

% Differentiation Matrix:
Dxx = toeplitz([-2 1 zeros(1,N_x-2)]);
Dxx(1,2) = 2;
Dxx(end, end-1) = 2;
Dxx = k*Dxx/dx^2;



% Initial Condition
f = @(x) cos(x*2*pi) + 1;
U0 = f(x);

% Define the Right hand side:

dUdt = @(t,U) Dxx * U - alpha * U .* (1 - U/m);

% Solve the heat equation:
tic
options = odeset('Stats', 'on');
[t,U] = ode23s(dUdt, t, U0, options);
toc


% Animation:

for i = 1: N_t

    figure(1)
    plot(x,U(i,:), 'w', 'LineWidth',3); hold off;
    if i == 1
        pause
    end
    grid on
    grid minor
    ylim([0,2])

end