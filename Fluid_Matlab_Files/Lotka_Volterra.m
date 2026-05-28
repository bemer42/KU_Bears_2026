%% Lotka_Volterra

close all; clear; clc

% Discretize time: 
N = 1e4;
t_0 = 0;
t_end = 100;
t = linspace(t_0,t_end,N);

u0 = 1;
v0 = 2;

U0 = [u0,v0];




% Parameters
a = .07;
b = .02;
y = .02;
d = .1;


% Right hand side:


dUdt = @(t, U)[a*U(1)-b*U(1)*U(2);
                -y*U(2)+d*U(1)*U(2)];

% Numerically solve the IVP
[t,U] = ode45(dUdt, t, U0);

u = U(:,1);
v = U(:,2);




%Plot
figure(1);
plot(t, u, 'LineWidth', 1.5); hold on;
plot(t, v, 'LineWidth', 1.5);
title('Lotka-Volterra');
xlabel('Time (t)');
ylabel('Population');
legend('Prey','Predator');
grid on;
grid minor;
    
%%Phase Plane
N_mesh = 20;
u_vec = linspace(min(u),max(u),N_mesh);
v_vec = linspace(min(v),max(v),N_mesh);

[U_mesh,V_mesh] = meshgrid(u_vec,v_vec);

Fu_vec = @(u,v) a*u-b*u*v;
Fv_vec = @(u,v) -y*v+d*u*v;

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
