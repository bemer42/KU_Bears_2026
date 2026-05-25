%% Growth_Decay

close all; clear; clc

% Discretize time: 
N = 1e4;
t_0 = 0;
t_end = 20;
t = linspace(t_0,t_end,N);

% Parameters
k = 0.1;

% Initial Condition
u0 = 1;

% Right hand side:

dudt = @(t,u)k*u;

% Numerically solve the IVP
[t,u] = ode45(dudt, t, u0);

%Plot
figure(1)
plot(t,u,'k','lineWidth',3)
title('Exponential Growth', 'FontSize',10)
xlabel('time(t)', 'FontSize', 8)
ylabel('u(t)', 'FontSize', 8)
set(gca, 'fontsize', 6)
grid on
grid minor

