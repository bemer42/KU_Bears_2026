%% Classic SIR

close all; clear; clc

% Discretize time: 
N = 1e4;
t_0 = 0;
t_end = 100;
t = linspace(t_0,t_end,N);

% Initial conditions
S0 = .99;
I0 = .01;
R0 = 0;

U0 = [S0; I0; R0];


% Parameters
alpha = .1;
beta = 1.4;
gamma = .4;
mu = .1;


% Right hand side:


dUdt = @(t, U)[alpha-beta*U(1).*U(2) - mu*U(1);
                beta*U(1).*U(2)-gamma*U(2) - mu*U(2);
                gamma*U(2) - mu*U(3)];

% Numerically solve the IVP
[t,U] = ode45(dUdt, t, U0);

S = U(:,1);
I = U(:,2);
R = U(:,3);




%Plot
figure(1);
plot(t, S, 'LineWidth', 1.5); hold on;
plot(t, I, 'LineWidth', 1.5);
plot(t, R, 'LineWidth', 1.5);
title('Classic SIR');
xlabel('Time (t)');
ylabel('Population');
legend('S','I', 'R');
grid on;
grid minor;

%%Phase Plane
N_mesh = 8;
vec = linspace(0,1,N_mesh);


[Sm,Im,Rm] = meshgrid(vec,vec,vec);

Fs = @(S,I,R) alpha * beta*S.*I-mu*S;
Fi = @(S,I,R) beta *S.*I - gamma*I - mu*I;
Fr = @(S,I,R) gamma*I - mu*R;

figure(2)
%quiver3(Sm,Im, Rm, Fs(Sm,Im,Rm),Fi(Sm,Im,Rm),Fr(Sm,Im,Rm));
%hold on;
plot3(S, I, R, 'k-', 'LineWidth', 2); hold on;
plot3(S0,I0,R0, 'go','LineWidth', 4)
plot3(S(end),I(end),R(end), 'ro', 'LineWidth',3)
plot3(S,zeros(size(t)),R,'w', 'LineWidth',.5);
plot3(S,I, zeros(size(t)),'w', 'LineWidth',.5);
plot3(zeros(size(t)),I,R,'w', 'LineWidth',.5);
xlabel('S');
ylabel('I');
zlabel('R');
xlim([0,1]);
ylim([0,1]);
zlim([0,1]);
grid on;
grid minor;