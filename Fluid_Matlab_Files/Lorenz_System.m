%% Lorenz System

close all; clear; clc

% Discretize time: 
N = 1e4;
t_0 = 0;
t_end = 100;
t = linspace(t_0,t_end,N);

% Initial conditions
x0 = .99;
y0 = .01;
z0 = 0;

U0 = [x0; y0; z0];


% Parameters
sigma = 10;
beta = 8/3;
ro = 28;


% Right hand side:


dUdt = @(t, U)[ sigma.*(U(2)-U(1));
                U(1).*(ro-U(3)) - U(2);
                U(1).*U(2) - beta*U(3) ];

% Numerically solve the IVP
[t,U] = ode45(dUdt, t, U0);

x = U(:,1);
y = U(:,2);
z = U(:,3);




%Plot
figure(1);
plot(t, x, 'LineWidth', 1.5); hold on;
plot(t, y, 'LineWidth', 1.5);
plot(t, z, 'LineWidth', 1.5);
title('Lorenz System');
xlabel('Time (t)');
ylabel('Population');
legend('x','y', 'z');
grid on;
grid minor;

%%Phase Plane
N_mesh = 8;
vec = linspace(0,1,N_mesh);


[Sm,Im,Rm] = meshgrid(vec,vec,vec);

Fs = @(S,I,R) sigma.*(y-x);
Fi = @(S,I,R)  x.*(ro-z) - y;
Fr = @(S,I,R) x.*y - beta*z;

figure(2)
%quiver3(Sm,Im, Rm, Fs(Sm,Im,Rm),Fi(Sm,Im,Rm),Fr(Sm,Im,Rm));
%hold on;
plot3(x, y, z, 'k-', 'LineWidth', 2); hold on;
plot3(x0,y0,z0, 'go','LineWidth', 4)
plot3(x(end),y(end),z(end), 'ro', 'LineWidth',3)
plot3(x,zeros(size(t)),z,'w', 'LineWidth',.5);
plot3(x,y, zeros(size(t)),'w', 'LineWidth',.5);
plot3(zeros(size(t)),y,z,'w', 'LineWidth',.5);
xlabel('x');
ylabel('y');
zlabel('z');
xlim([0,1]);
ylim([0,1]);
zlim([0,1]);
grid on;
grid minor;