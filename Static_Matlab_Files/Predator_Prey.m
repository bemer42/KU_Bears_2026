%% Predator-Prey
close all; clear; clc

% Discretize time
N     = 1e4; 
t_0   = 0;
t_end = 30;
t     = linspace(t_0, t_end, N);

% Parameters
a = 1;
b = .1;
c = 2;
d = .08;

% Initial condition
u0 = 10;
v0 = 8;

U0 = [u0; v0];

% Right hand side
dudt = @(t,u,v) a*u - b*u.*v;
dvdt = @(t,u,v) -c*v + d*u.*v; 

dUdt = @(t,U) [dudt(t,U(1),U(2));
               dvdt(t,U(1),U(2))];

% Numerically solve the IVP
[t, U] = ode45(dUdt, t, U0);
u = U(:,1);
v = U(:,2);

% Plot the time series
figure(1)
plot(t,u,'b','linewidth',5)
hold on
plot(t,v,'r','linewidth',5)
set(gca, 'fontsize',12)
title('Predator-Prey -- Time-Series','fontsize',20)
xlabel('time (t)','fontsize',16)
ylabel('u(t), v(t)','fontsize',16)
legend('u(t)', 'v(t)', 'fontsize', 12)
grid on
grid minor

% Phase Plane with vector field
N = 30; 
buff = .3;
u_vec = linspace(min(u)-buff*abs(min(u)), max(u)+buff*abs(max(u)), N);
v_vec = linspace(min(v)-buff*abs(min(v)), max(v)+buff*abs(max(v)), N);
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
title('Predator-Prey -- Phase Plane','fontsize',20)
xlabel('u(t)','fontsize',16)
ylabel('v(t)','fontsize',16)
xlim([u_vec(1) u_vec(end)])
ylim([v_vec(1) v_vec(end)])
grid on
grid minor

