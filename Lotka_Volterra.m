%% Lotka_Volterra

close all; clear; clc

% Discretize time: 
N = 1e4;
t_0 = 0;
t_end = 20;
t = linspace(t_0,t_end,N);

u0 = linspace(1,139,20);
v0 = linspace(1,139,20);

u = 1;
v = 2;




% Parameters
a = .5;
b = 9;
y = 10;
s = 20;




% Right hand side:

dudt = @(t,u)(a*u-b*u*v);
dvdt = @(t,v)(-y*v+s*u*v);

% Numerically solve the IVP
[t,u] = ode45(dudt, t, u0);
[t,v] = ode45(dudt, t, v0);

% Define Equilibrium:
% u_star = us + 2*sin(t);

%Plot
figure;
subplot(2,1,1);
plot(t, u, 'LineWidth', 1.5); hold on;
plot(t, v, 'LineWidth', 1.5);
title('Lotka-Volterra');
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
