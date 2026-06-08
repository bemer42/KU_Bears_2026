%% Spatial ODE with space varying Wnt Signal

close all; clear; clc

tic

% Discretize time: 
N_t = 1e3;
t_0 = 0;
t_end = 5e3;
t = linspace(t_0,t_end,N_t);

% Discretize space:
N_z = 83;
z_0 = 0;
z_end = 82;
z = linspace(z_0, z_end, N_z);

% Define a Wnt signal:

W = @(z) exp(-z/10);

% Define initial conditon:

A0 = @(z) 25 + (114-25)/82*z;
B0 = @(z) 289 + (23-289)/82*z;
X0 = @(z) 9e-4 + (5e-4 - 9e-4)/82*z;

% Initialize matrices:

A = zeros(N_t,N_z); B = A; X = A;


% Solve for APC, Beta-Cat, and Axin for every z value

for i = 1: N_z
    
    U0 = [A0(z(i)), B0(z(i)), X0(z(i))];
    [t,U] = ABX_fun(t,U0 , W(z(i)));

    A(:, i) = U(:, 1);
    B(:, i) = U(:, 2);
    X(:, i) = U(:, 3);


end

toc

[T,Z] = meshgrid(t,z);


%%  APC Plot: 
figure(1)
surf(T,Z,A'); 
hold on;
shading interp
title('APC','fontsize',20)
xlabel('time (t)','fontsize',16)
ylabel('crypt level (z)','fontsize',16)
set(gca,'fontsize',14)

% Plot of trajectories
for i = 1:6:N_z

    figure(1)
    plot3(t,z(i)*ones(size(t)),A(:,i),'k','linewidth',3); hold on;

end

%%  Beta Catenin Plot: 
figure(2)
surf(T,Z,B'); 
hold on;
shading interp
title('Beta-Catenin','fontsize',20)
xlabel('time (t)','fontsize',16)
ylabel('crypt level (z)','fontsize',16)
set(gca,'fontsize',14)

% Plot of trajectories
for i = 1:6:N_z

    figure(2)
    plot3(t,z(i)*ones(size(t)),B(:,i),'k','linewidth',3); hold on;

end



%%  Axin Plot: 
figure(3)
surf(T,Z,X'); 
hold on;
shading interp
title('Axin','fontsize',20)
xlabel('time (t)','fontsize',16)
ylabel('crypt level (z)','fontsize',16)
set(gca,'fontsize',14)

% Plot of trajectories
for i = 1:6:N_z

    figure(3)
    plot3(t,z(i)*ones(size(t)),X(:,i),'k','linewidth',3); hold on;

end

%% Animation

for i=1: N_t

    figure(4)
    plot(z,A(i, :), 'w', 'LineWidth',3); hold on;
    if i == 1
        pause
    end

end