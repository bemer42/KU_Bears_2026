%% APC_Axin_Beta_Catenin Model

close all; clear; clc

% Discretize time: 
N = 1e5;
t_0 = 0;
t_end = 5e6;
t = linspace(t_0,t_end,N);

A0 = 90;
B0 = 200;
X0 = 2e-3;

U0 = [A0,B0,X0];

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
k15 = 0.333;
k19 = 8.33e-4;
K7 = 50.0;
K8 = 120.0;
K16 = 30.0;
K17 = 1200;
K20 = 1;
K21 = 1;
Kt = 39.2115;
Kb = 34.0445;
GSK0 = 50.0;
Dsh0 = 100.0;
TCF0 = 15.0;
v18 = 0.1774;
v12 = 0.423;
v14 = 8.22e-5;
Km = 98;

W = 1;

% Right hand side:

f = @(W) (k1*k4*k6*k9*K21/ (k5*K7*K8))*((W+k2/k1)./(k1*(Dsh0*k3+k_6)*W+k2*k_6));
g = @(A, B, X , W) ((K21+X)./GSK0)+(((k9+k10).*A.*X)./(k9*k10*GSK0)).*(((k4+k5)*K8*k10)./(k4*(k9+k10))+B).*f(W);
F = @(A, B) (v18./(1+(TCF0.*B)./(Kt.*(K16+B))+B./Kb))-k19.*A;
G = @(A, B, X, W) v12 - (k13 + A.*X.*(f(W)./g(A,B,X, W))+F(A,B)./K17).*B ;
H = @(A, B, X)  v14 - ((k15.*A)./(Km+A)+F(A,B)./K7).*X;


dUdt = @(t, U)[F(U(1),U(2)) ; 
        
      ((1 + U(1)/K7 + U(2)/K20 + K21./((K21+U(3)).*g(U(1),U(2),U(3),W))).*G(U(1),U(2),U(3),W) - (U(2)/K20).*H(U(1),U(2),U(3)))./ ...
      ((1+ U(1)/K7 + U(3)/K20).*(1+ U(1)/K7 + U(2)/K20 + K21./((K21+U(3)).*g(U(1),U(2),U(3),W)))-(U(2).*U(3)/(K20^2)));

      ((-U(3)/K20).*G(U(1),U(2),U(3),W)+(1+U(1)/K7+U(3)/K20).*H(U(1),U(2),U(3)))./ ...
      ((1+U(1)/K7+U(3)/K20).*(1+ U(1)/K7+ U(2)/K20 + K21./((K21+U(3)).*g(U(1),U(2),U(3),W)))-(U(2).*U(3)/(K20^2))) ];

% Numerically solve the IVP
[t,U] = ode23s(dUdt, t, U0);

A = U(:,1);
B = U(:,2);
X = U(:,3);

sigma = (B .* sqrt(TCF0)) ./ (A .* sqrt(K16 + B));
dsigmadA = -(B .* sqrt(TCF0)) ./ (A.^2 .* sqrt(K16 + B));
dsigmadB = (sqrt(TCF0) .* (2*K16 + B)) ./ (2 .* A .* (K16 + B).^(1.5));



%Plot
figure(1);
plot(t,A,'k','linewidth',3); hold on;
plot(t,B,'k','linewidth',3);
plot(t,X,'k','linewidth',3);
title('APC-Axin-Beta-Catenin Model');

grid on;
grid minor;

%%Phase Plane
N_mesh = 5;
vec = linspace(0,1,N_mesh);


[Am,Bm,Xm] = meshgrid(vec,vec,vec);

% Fa = @(A,B,X) F(A,B);
% Fb = @(A,B,X) ((1 + A/K7 + B/K20 + K21./((K21+X).*g(A,B,X,W))).*G(A,B,X,W) - (B/K20).*H(A,B,X))./ ...
%       ((1+ A/K7 + X/K20).*(1+ A/K7 + B/K20 + K21./((K21+X).*g(A,B,X,W)))-(B.*X/(K20^2)));
% 
% Fx = @(A,B,X) ((-X/K20).*G(A,B,X,W)+(1+A/K7+X/K20).*H(A,B,X))./ ...
%       ((1+A/K7+X/K20).*(1+ A/K7+ B/K20 + K21./((K21+X).*g(A,B,X,W)))-(B.*X/(K20^2)));

figure(2)
plot3(A, B, X, 'k-', 'LineWidth', 2); hold on;
plot3(A0,B0,X0, 'go','LineWidth', 4)
plot3(A(end),B(end),X(end), 'ro', 'LineWidth',3)
plot3(A,zeros(size(t)),X,'k', 'LineWidth',.5);
plot3(A,B, zeros(size(t)),'k', 'LineWidth',.5);
plot3(zeros(size(t)),B,X,'k', 'LineWidth',.5);
xlabel('APC');
ylabel('Beta Catenin');
zlabel('Axin');
xlim([0, 100]);
ylim([0, 500]);
zlim([0, 4e-3]);
grid on;
grid minor;

