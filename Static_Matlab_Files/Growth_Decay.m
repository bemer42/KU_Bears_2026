%% Exponential Growth/Decay
close all; clear; clc

% Discretize time
N     = 1e4; 
t_0   = 0;
t_end = 20;
t     = linspace(t_0, t_end, N);

% Parameters
k = 0.1;

% Initial condition
u0 = 1;

% Right hand side
dudt = @(t,u) k*u;

% Numerically solve the IVP
[t, u] = ode45(dudt, t, u0);

% Plot
figure(1)
plot(t,u,'k','linewidth',3)
title('Exponential Growth','fontsize',20)
xlabel('time (t)','fontsize',16)
ylabel('u(t)','fontsize',16)
set(gca, 'fontsize',10)
grid on
grid minor



