% Cell movement PDE in one dimension along crypt shape
clear; close all; clc

% Discretize time: 
N_t = 1e3;
t_0 = 0;
t_end = 2e5;
t = linspace(t_0,t_end,N_t);

% Discretize space:
N_z = 2e2;
z_0 = 0;
z_end = 78;
z = linspace(z_0, z_end, N_z)';
dz = z(2) - z(1);

% Define the differentiation matrix:
Dz = diag(1/2*ones(N_z-1,1),1) - diag(1/2*ones(N_z-1,1),-1);
Dz(1,1:3) = [-3/2 2 -1/2];
Dz(end, end-2:end) = [1/2 -2 3/2];
Dz = Dz/dz;


% Define parameters

r_b = 41/2/pi;
r_t = 10/pi;   
a = 0.3;   
alpha = 0.05;
k = 1e-6;


% Crypt width
x = @(z) r_b*(1 - exp(-a * z)) + r_t * exp(a * (z - z_end));
xp = @(z) -a*r_b*exp(-a*z) + a*r_t*exp(a*z-z_end);

% Initial condition:
%q_0 = 1 + ones(size(z)).*(z<30);
q_0 = ones(size(z));
q_0_int = q_0(2:end-1);



%Define right hand side function:
dQdt = @(t,Q) dqdt_snipsnap(t, Q, z, xp(z), Dz, alpha, k);

%Solve the system of ODEs that represents the PDE:
tic
options = odeset('Stats', 'on');
[t,Q] = ode23s(dQdt, t, q_0_int, options);
toc


% Plot an animation:

dt = round(.01*N_t);

for i = 1: dt : N_t

    Q_l = (-2*Q(i,1) + 1/2*Q(i,2))/(-3/2);
    Q_r = 1;
    Q_full = [Q_l Q(i, :) Q_r];

    % Average Value
    avg = trapz(z, Q_full) / (z_end-z_0);

    figure(1)
    plot(z,Q_full, 'k', 'LineWidth',5);
    title(sprintf('Time Step: %d    Average Density: %.4f', i, avg));
    grid on
    grid minor
    
    if i == 1
        pause
    else
        drawnow
    end
    ylim([0 2])
    hold off;
   
    

end

%% Plot on crypt domain
figure(2)
for i = 1: dt : N_t
    Q_l = (-2*Q(i,1) + 1/2*Q(i,2))/(-3/2);
    Q_r = 1;
    Q_full = [Q_l, Q(i, :), Q_r];
    scatter3(x(z), zeros(size(z)),z, 20, Q_full)
    colormap("turbo")
    colorbar
    xlim([-30 30])
    ylim([-30 30])
    zlim([0 80])
    clim([0 2])
    grid on; 

    drawnow
end
