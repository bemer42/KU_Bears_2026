%% Newton's Law of Cooling
close all; clear; clc

% Discretize time
N     = 1e4; 
t_0   = 0;
t_end = 20;
t     = linspace(t_0, t_end, N);

% For loop over initial condition
u0_vec = linspace(1,139,20);
for i = 1:length(u0_vec)

    % Parameters
    k  = .5;
    us = 70;
    
    % Initial condition
    u0 = u0_vec(i);
    
    % Right hand side
    dudt = @(t,u) -k*(u - (us + 2*sin(t)));
    
    % Numerically solve the IVP
    [t, u] = ode45(dudt, t, u0);
    
    % Define equlibrium:
    u_star = us + 2*sin(t);
    
    % Plot
    figure(1)
    plot(t,u,'k','linewidth',3)
    hold on
    plot(t,u_star,'r:','linewidth',3)
    set(gca, 'fontsize',12)
    title('Newton Law of Cooling','fontsize',20)
    xlabel('time (t)','fontsize',16)
    ylabel('u(t)','fontsize',16)
    grid on
    grid minor

end


