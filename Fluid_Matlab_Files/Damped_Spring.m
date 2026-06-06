%% Damped Spring

close all; clear; clc

% Discretize time

N = 1e4;
t_0 = 0;
t_end = 200;
t = linspace(t_0,t_end,N);

% Parameters

k = 0.5;
b = 0;
m = 5;

u0 = 1;
v0 = 2;

U0 = [u0;v0];


% Right hand side
dudt = @(t,u,v) v;
dvdt = @(t,u,v) -k/m*u - b*v/m;

dUdt = @(t, U)[dudt(t,U(1),U(2));
               dvdt(t,U(1),U(2))];

% Solve numerically using ode45

[t, U] = ode45(dUdt, t, U0);

% get u and v from matrix
u = U(:, 1);
v = U(:, 2); 

%Plot
figure(1);
plot(t, u, 'LineWidth', 1.5); hold on;
plot(t, v, 'LineWidth', 1.5);
title('Damped Mass-on-a-Spring System Dynamics');
xlabel('Time (t)');
ylabel('State Values');
legend('u(t','v(t)');
grid on;
grid minor;

A = [0 1;
    -k/m -b/m];
eigs(A)

%%Phase Plane
N_mesh = 20;
u_vec = linspace(min(u),max(u),N_mesh);
v_vec = linspace(min(v),max(v),N_mesh);

[U_mesh,V_mesh] = meshgrid(u_vec,v_vec);

Fu_vec = @(u,v) v;
Fv_vec = @(u,v) -k/m*u -b/m*v;

figure(2)
quiver(U_mesh,V_mesh, Fu_vec(U_mesh,V_mesh),Fv_vec(U_mesh,V_mesh));
hold on;
plot(u, v, 'k-', 'LineWidth', 2); hold on;
plot(u0,v0, 'go','LineWidth', 4)
plot(u(end),v(end), 'ro', 'LineWidth',3)
xlabel('u');
ylabel('v');
grid on;
grid minor;