%% Competing Species
close all; clear; clc

% Discretize time
N     = 1e4; 
t_0   = 0;
t_end = 30;
t     = linspace(t_0, t_end, N);

% Parameters
k1 = 1; 
k2 = 1.4;
m1 = 50;
m2 = 80;
d1 = .001;
d2 = .007;

% Initial condition
u0 = 2;
v0 = 8;

U0 = [u0; v0];

% Right hand side
dudt = @(t,u,v) k1*u.*(1-u/m1) - d1*u.*v;
dvdt = @(t,u,v) k2*v.*(1-v/m2) - d2*u.*v;

dUdt = @(t,U) [dudt(t,U(1),U(2));
               dvdt(t,U(1),U(2))];

% Numerically solve the IVP
[t, U] = ode45(dUdt, t, U0);
u = U(:,1);
v = U(:,2);

% Coexistence criteria: 
d1*m1/k1 - m1/m2
d2*m2/k2 - m2/m1

d1*m1/k1*d2*m2/k2

% Plot the time series
figure(1)
plot(t,u,'b','linewidth',5)
hold on
plot(t,v,'r','linewidth',5)
set(gca, 'fontsize',12)
title('Competing Species -- Time-Series','fontsize',20)
xlabel('time (t)','fontsize',16)
ylabel('u(t), v(t)','fontsize',16)
legend('u(t)', 'v(t)', 'fontsize', 12)
grid on
grid minor

% Phase Plane with vector field
N = 30; 
buff = .3;
u_vec = linspace(0, max(u)+buff*abs(max(u)), N);
v_vec = linspace(0, max(v)+buff*abs(max(v)), N);
[Umesh,Vmesh] = meshgrid(u_vec,v_vec);

Norm = sqrt(dudt(0,Umesh,Vmesh).^2 + dvdt(0,Umesh,Vmesh).^2);

figure(2)
quiver(Umesh, Vmesh, dudt(0,Umesh,Vmesh)./Norm,dvdt(0,Umesh,Vmesh)./Norm,.5,...
    'linewidth',2,'color',[.75 .75 .75],'ShowArrowHead','on')
hold on
plot(u,v,'k','linewidth',3)
plot(u0,v0,'go','linewidth',5)
plot(u(end),v(end),'ro','linewidth',4)
set(gca, 'fontsize',12)
title('Competing Species -- Phase Plane','fontsize',20)
xlabel('u(t)','fontsize',16)
ylabel('v(t)','fontsize',16)
xlim([u_vec(1) u_vec(end)])
ylim([v_vec(1) v_vec(end)])
grid on
grid minor

