%% Dynamic_Cell_1D
% Cell movement PDE in one dimension along crypt shape with evolving r(z,t)
clear; close all; clc

% Grid/Time Points:
Nt    = 1e3;
t_end = 1e3;
Nz    = 1e2;

% Predefine structures:
geom    = struct();
diffmat = struct();
par     = struct();

% Define BCs for q:
bcType  = "NbDt";

% Define parameters:
alpha_z = 1.6406e3;
zp      = 27; 
Tc      = 15.1954;
eta     = 1e-2;
epsQ    = .01;
alpha_s = 1e-4;
mu      = 1.2e-2;
beta    = 1.2e-2;
epsS    = 0.01;
s_th    = 0.5;
q_th    = 1.06;

% Define geometry parameters:
r_b   = 41/2/pi;
r_t   = 10/pi;
a     = 0.3;

% Discretize time:
t_0   = 0;
tspan = linspace(t_0, t_end, Nt);

% Discretize space:
L     = 78.8783;
z     = linspace(0, L, Nz)';
dz    = z(2) - z(1);

% Define the differentiation matrix for q and r:
Dz = diag(1/2*ones(Nz-1,1),1) - diag(1/2*ones(Nz-1,1),-1);
Dz(1,1:3) = [-3/2 2 -1/2];
Dz(end, end-2:end) = [1/2 -2 3/2];
Dz = Dz/dz;

% Upwind derivative and diffusion matrices for s:
DzF = spdiags([-ones(Nz,1) ones(Nz,1)], [0 1], Nz, Nz)/dz;   
DzB = spdiags([-ones(Nz,1) ones(Nz,1)], [-1 0], Nz, Nz)/dz;  
DzF(end,end-1:end) = [-1 1]/dz;  
DzB(1,1:2)         = [-1 1]/dz;   

% Initial condition:
q0     = ones(size(z));
q0_int = q0(2:end-1);

% Initial Condition for crypt radius:
r0 = r_b*(1 - exp(-a*z)) + r_t*exp(a*(z - L));

% Initial Condition for signal:
s0     = zeros(size(z));
s0_int = s0(2:end-1);

% Coupled initial condition:
y0 = [q0_int; r0; s0_int];

% Build structures for solver:
diffmat.Dz  = Dz;
diffmat.DzF = DzF;
diffmat.DzB = DzB;

geom.z   = z;
geom.L   = L;
geom.r_t = r_t;
geom.r_b = r_b;
geom.a   = a;
geom.r0  = r0;

par.alpha_z = alpha_z;
par.zp      = zp;
par.Tc      = Tc;
par.eta     = eta;
par.epsQ    = epsQ;
par.alpha_s = alpha_s;
par.mu      = mu;
par.beta    = beta;
par.epsS    = epsS;
par.s_th  = s_th;
par.q_th     = q_th;

% Define right hand side function:
dYdt = @(t,y) dqdt_dyn_1D_snipsnap(t, y, diffmat, geom, par, bcType);

% Solve the system:
tic
options = odeset('Stats','on');
[t,Y] = ode15s(dYdt, tspan, y0, options);
toc

% Unpack solution:
Nq = Nz - 2;
Nr = Nz;
Ns = Nz - 2;

Q_int_sol  = Y(:,1:Nq);                 % Nt x (Nz-2)
R_full_sol = Y(:,Nq+1:Nq+Nr).';         % Nz x Nt
S_int_sol  = Y(:,Nq+Nr+1:end);          % Nt x (Nz-2)  (optional)

% Define Storages:
Q_full = zeros(Nz, numel(t));
V_full = Q_full; R_full = Q_full; Rz_full = Q_full;
G_full = Q_full; S_full = Q_full;

for i = 1:numel(t)

    % Radius from solution:
    r  = R_full_sol(:,i);
    rz = Dz*r;
    g  = sqrt(1 + rz.^2);

    % Apply boundary conditions to reconstruct q_full
    q_int = Q_int_sol(i,:)';

    switch bcType
        case "NbDt"
            Q_l = -(Dz(1,2)*q_int(1) + Dz(1,3)*q_int(2))/Dz(1,1);
            Q_r = 1;
        case "NbNt"
            Q_l = -(Dz(1,2)*q_int(1) + Dz(1,3)*q_int(2))/Dz(1,1);
            Q_r = -(Dz(end,end-1)*q_int(end) + Dz(end,end-2)*q_int(end-1))/Dz(end,end);
        case "DbNt"
            Q_l = 1;
            Q_r = -(Dz(end,end-1)*q_int(end) + Dz(end,end-2)*q_int(end-1))/Dz(end,end);
        case "DbDt"
            Q_l = 1;
            Q_r = 1;
        otherwise
            error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
    end

    q_full = [Q_l; q_int; Q_r];

    % Reconstruct s_full from s_int (Neumann BCs: s_z = 0 at both ends)
    s_int = S_int_sol(i,:)';

    S_l = -(Dz(1,2)*s_int(1) + Dz(1,3)*s_int(2))/Dz(1,1);
    S_r = -(Dz(end,end-1)*s_int(end) + Dz(end,end-2)*s_int(end-1))/Dz(end,end);

    s_full = [S_l; s_int; S_r];

    Q_full(:,i)  = q_full;
    S_full(:,i)  = s_full;

    R_full(:,i)  = r;
    Rz_full(:,i) = rz;
    G_full(:,i)  = g;

    % Velocity using time-dependent g:
    V_full(:,i) = -(1./g).*par.alpha_z./q_full.^3.*(Dz*q_full);
end

% Steady State Analysis
Q_ss = Q_full(:,end);
V_ss = V_full(:,end);
R_ss = R_full(:,end);
g_ss = G_full(:,end);

% Total number of cells (using r(z,t_end) and g(z,t_end))
Total_Cells  = 2*pi*trapz(z, R_ss .* Q_ss .* g_ss)
Prolif_Cells = 2*pi*trapz(z, R_ss .* Q_ss .* (z<par.zp) .* g_ss)

% Crypt renewal time
Arc_Length = cumtrapz(z, g_ss);
pos = find(Arc_Length >= 1);
Crypt_Renewal_Time = trapz(z(pos), g_ss(pos)./V_ss(pos))/24

% Plots
figure(1)
plot(z, Q_ss,'k', 'LineWidth', 5); hold on
plot(z, V_ss,'k:','LineWidth', 5)
set(gca,'fontsize',24)
title('Cell Density and Velocity','interpreter','latex')
xlabel('$z$','interpreter','latex')
ylabel('$q_s(z)$ and $v_s(z)$','interpreter','latex')
legend('$q_s(z)$','$v_s(z)$','interpreter','latex')
grid on; grid minor
xlim([0 L])

figure(2)
plot(z, R_ss,'k','LineWidth',5)
set(gca,'fontsize',24)
title('Crypt Radius at Final Time','interpreter','latex')
xlabel('$z$','interpreter','latex')
ylabel('$r(z,t_{end})$','interpreter','latex')
grid on; grid minor
xlim([0 L])

%% Animation:
dt = ceil(.01*Nt);
for i = 1:dt:Nt

    Q  = Q_full(:,i);
    S  = S_full(:,i);
    V  = V_full(:,i);
    R  = R_full(:,i);
    Rz = Rz_full(:,i);
    G  = G_full(:,i);

    figure(3)
    subplot(2,2,[1 3])
    plot(z,Q, 'k', 'LineWidth',5);
    hold on
    plot(z,V,'k--','linewidth',5)
    plot(z,q_th*ones(size(z)),'k:','linewidth',2)
    set(gca,'fontsize',42)
    title('Cell Density and Velocity','fontsize',55,'interpreter','latex')
    xlabel('$z$','fontsize',50,'interpreter','latex')
    ylabel('$q(z,t)$ and $v(z,t)$','fontsize',50,'interpreter','latex')
    legend('$q(z,t)$','$v(z,t)$','$q_{th}$','fontsize',50,'interpreter','latex')
    grid on
    grid minor
    xlim([0 L])
    ylim([min(min([Q_ss; V_ss])) max(max([Q_ss; V_ss]))+.5])
    hold off
    subplot(2,2,2)
    plot(z,R,'k','linewidth',5)
    set(gca,'fontsize',42)
    title('Crypt Radius','fontsize',55,'interpreter','latex')
    xlabel('$z$','fontsize',50,'interpreter','latex')
    grid on
    grid minor
    xlim([0 L])
    ylim([0 max(2*r0)])
    hold off
    subplot(2,2,4)
    plot(z,S,'k','linewidth',5)
    hold on
    plot(z,s_th*ones(size(z)),'k:','linewidth',5)
    set(gca,'fontsize',42)
    title('Signal','fontsize',55,'interpreter','latex')
    xlabel('$z$','fontsize',50,'interpreter','latex')
    grid on
    grid minor
    xlim([0 L])
    ylim([0 beta/mu])
    hold off
    if i == 1
        pause
    end
end

%% Animation on crypt domain
dt = ceil(.01*Nt);
theta = linspace(0,2*pi,Nz);

for k = 1:dt:Nt

    [T,R] = meshgrid(theta,R_full(:,k));

    X_plot = R.*cos(T);
    Y_plot = R.*sin(T);
    Z_plot = repmat(z,1,Nz);

    Q_dens = repmat(Q_full(:,k)',Nz,1);

    figure(3)
    surf(X_plot,Y_plot,Z_plot,Q_dens')
    set(gca,'fontsize',42)
    title('Cell Density on Crypt','fontsize',55,'interpreter','latex')
    xlabel('$x$','fontsize',50,'interpreter','latex')
    ylabel('$y$','fontsize',50,'interpreter','latex')
    zlabel('$z$','fontsize',50,'interpreter','latex')
    grid on
    grid minor
    box on
    shading interp
    colormap turbo
    colorbar
    lighting gouraud
    xlim([-40 40])
    ylim([-40 40])
    zlim([0 80])
    caxis([1 1.1])

    if k == 1
        pause
    end
end
