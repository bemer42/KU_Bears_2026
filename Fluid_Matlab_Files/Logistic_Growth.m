%% Logistic Growth

close all; clear; clc

% Discretize time

N = 1e4;
t_0 = 0;
t_end = 20;
t = linspace(t_0,t_end,N);

% Parameters

k = 0.5;
M = 100;

u0_vector = [10, 50, 100, 150];

% Right hand side

dudt = @(t, u) k * u * (1 - u / M);

for i = 1:length(u0_vector)
    u0 = u0_vector(i);

    % Solve numerically using ode45
    [t, u] = ode45(dudt, t, u0);

    plot(t,u,'k','lineWidth',3)
    hold on
    plot(t, u, 'r:', 'LineWidth',3)
    title('Logistic Growth', 'FontSize',10)
    xlabel('time(t)', 'FontSize', 8)
    ylabel('u(t)', 'FontSize', 8)
    set(gca, 'fontsize', 6)
    grid on
    grid minor

end

    

