%% APC_Axin_Beta_Catenin Model

close all; clear; clc

% Discretize time: 
N = 1e4;
t_0 = 0;
t_end = 20;
t = linspace(t_0,t_end,N);

u0 = 1;
v0 = 2;

U0 = [u0,v0];




% Parameters
k1 = 0.182;
k2 = 1.82e-2;
k3 = 5.0e-2;
k4 = 0.267;
k5 = 0.133;
k6 = 9.09e-2;
k_6 = 0.909;
k9 = 206.0;
k10 = 206.0;
k11 = 0.417;
k13 = 2.57e-4;
k19 = 8.33e-4;
K7 = 50.0;
K8 = 120.0;
K16 = 30.0;
K17 = 1;
K20 = 1;
K21 = 1;
Kt = 39.2115;
Kb = 34.0445;
GSK0 = 50.0;
W = 1;
DSH0 = 100.0;
TCF0 = 15.0;
v18 = 0.1774;


% Right hand side:

f = @(W) (k1*k4*k6*k9*K21/k5*K7*K8)*((W+k2/k1)/k1*(Dsh0k3+k_6)*W+k2*k_6);
g = @(A, B, X , W) ((K21+X)/GSK0)+((k9+k10)*A*X/k9*k10*GSK0)*(((k4+k5)*K8*k10/k4*(k9+k10)+B))*f;
F = @F(A, B) (v18/(1+(TCF0*B/Kt*(K6+B)+B/Kb)-k19*A);
G = @G(A, B, X, W);
H = @H(A, B, X);


dUdt = @(t, U)[k1*U(1)*(1-U(1)/m1) - d1*U(1)*U(2);
    k2*U(2)*(1-U(2)/m2) - d2*U(1)*U(2)];

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
legend('Species 1','Species 2');
grid on;
grid minor;
