%% Lorenz System
close all; clear; clc

% Time discretization
N     = 1e4;
t_0   = 0;
t_end = 50;
t     = linspace(t_0,t_end,N);

% Parameters
sigma = 10;
beta  = 8/3;
rho   = 28;

% Initial conditions
x0 = 0.5;
y0 = 0.7;
z0 = 1;

U0 = [x0; y0; z0];

% Right Hand Side functions
dxdt = @(t,x,y,z) sigma.*(y-x);
dydt = @(t,x,y,z) rho.*x-y-x.*z;
dzdt = @(t,x,y,z) x.*y-beta.*z;

dUdt = @(t,U) [dxdt(t,U(1),U(2),U(3));
    dydt(t,U(1),U(2),U(3));
    dzdt(t,U(1),U(2),U(3))];

% Solve the system
[t,U] = ode15s(dUdt,t,U0);

x = U(:,1);
y = U(:,2);
z = U(:,3);

% Time-Series
figure(1)
plot(t,x,'b','linewidth',5)
hold on
plot(t,y,'r','linewidth',5)
plot(t,z,'g','linewidth',5)
set(gca, 'fontsize',12)
title('Lorenz System -- Time-Series','fontsize',20)
xlabel('time (t)','fontsize',16)
ylabel('x(t), y(t), z(t)','fontsize',16)
legend('x(t)', 'y(t)', 'z(t)', 'fontsize', 12)
grid on
grid minor

% Phase Portait
dt = 4;
projections = 1;

figure(2)
plot3(x,y,z,'k','linewidth',1)
hold on
plot3(x(1:dt:end),y(1:dt:end),z(1:dt:end),'ko','linewidth',2)
if projections > 0
    plot3(min(x)*ones(size(t)),y,z,'k','linewidth',.5)
    plot3(x,min(y)*ones(size(t)),z,'k','linewidth',.5)
    plot3(x,y,min(z)*ones(size(t)),'k','linewidth',.5)
end
set(gca,'fontsize',20)
title('Lorenz System -- Phase Portrait','fontsize',35)
xlabel('x','fontsize',25)
ylabel('y','fontsize',25)
zlabel('z','fontsize',25)
grid on
grid minor
box on
view([1 1 1])
xlim([min(x) max(x)])
ylim([min(y) max(y)])
zlim([min(z) max(z)])

%% Movie Plot

tail = 300;
dt = 10;
figure(3)


for i = 1:dt:length(t)

    start = max(1,i-tail);

    figure(3)
    plot3(x,y,z,'b','linewidth',.2)
    hold on
    plot3(x(i),y(i),z(i),'ko','linewidth',5)
    hold on
    plot3(x(start:i),y(start:i),z(start:i),'k','linewidth',2)
    set(gca,'fontsize',20)
    title('Lorenz System -- Phase Portrait','fontsize',35)
    xlabel('x','fontsize',25)
    ylabel('y','fontsize',25)
    zlabel('z','fontsize',25)
    set(gca,'fontsize',20)
    title('Lorenz System -- Phase Portrait','fontsize',35)
    xlabel('x','fontsize',25)
    ylabel('y','fontsize',25)
    zlabel('z','fontsize',25)
    view([1 1 1])
    xlim([min(x) max(x)])
    ylim([min(y) max(y)])
    zlim([min(z) max(z)])
    grid on
    grid minor
    box on
    hold off
    if i == 1
        pause
    end
end













