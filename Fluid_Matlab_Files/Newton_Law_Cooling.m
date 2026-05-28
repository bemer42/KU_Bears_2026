%% Newton's Law of Cooling

close all; clear; clc

% Discretize time: 
N = 1e4;
t_0 = 0;
t_end = 20;
t = linspace(t_0,t_end,N);

u0_vec = linspace(1,139,20);

for i = 1: length(u0_vec)

    % Parameters
    k = .5;
    us = 70;
    
    % Initial Condition
    u0 = u0_vec(i);
    
    % Right hand side:
    
    dudt = @(t,u)-k*(u-(us + 2*sin(t)));
    
    % Numerically solve the IVP
    [t,u] = ode45(dudt, t, u0);
    
    % Define Equilibrium:
    u_star = us + 2*sin(t);
    
    %Plot
    figure(1)
    plot(t,u,'k','lineWidth',3)
    hold on
    plot(t, u_star, 'r:', 'LineWidth',3)
    title('Exponential Growth', 'FontSize',10)
    xlabel('time(t)', 'FontSize', 8)
    ylabel('u(t)', 'FontSize', 8)
    set(gca, 'fontsize', 6)
    grid on
    grid minor

end

