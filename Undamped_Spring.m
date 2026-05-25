%% Undamped Spring System

close all; clear; clc

% Discretize time

N = 1e4;
t_0 = 0;
t_end = 20;
t = linspace(t_0,t_end,N);

% Parameters

k = 0.5;
m = 50;

u0 = 1;
v0 = 2;
x0 = [u0,v0];


% Right hand side

dudt = @(t, x)[x(2); (-k/m)*x(1)];

% Solve numerically using ode45
    
[t, x] = ode45(dudt, t, x0);
    
% get u and v from matrix
u = x(:, 1);
v = x(:, 2); 
    
figure;
subplot(2,1,1);
plot(t, u, 'LineWidth', 1.5); hold on;
plot(t, v, 'LineWidth', 1.5);
title('Undamped Mass-on-a-Spring System Dynamics');
xlabel('Time (t)');
ylabel('State Values');
grid on;
grid minor;
    
   
subplot(2,1,2);
plot(u, v, 'g-', 'LineWidth', 2); hold on;
xlabel('u0');
ylabel('v0');
grid on;
grid minor;

