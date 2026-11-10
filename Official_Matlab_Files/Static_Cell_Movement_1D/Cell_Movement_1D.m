%% Cell movement PDE in one dimension along crypt shape
clear; close all; clc; ssp = 1; uniform = 0;

% -------------------------
% Grid/Time Points
% -------------------------
N_t   = 7e3;
t_end = 4e5;
Nz    = 1e3;

% -------------------------
% Predefine structures and boundary type:
% -------------------------
geom    = struct();
diffmat = struct();
par     = struct();
bcType  = "NbDt";

% -------------------------
% Define parameters:
% -------------------------
alpha = 908.29;
zp    = 23.4087;
Tc    = 15.1954;

% -------------------------
% Define geometry parameters:
% -------------------------
r_b   = 18/pi;
r_t   = 10/pi;
a     = 0.3;

% -------------------------
% Discretize time:
% -------------------------
t_0   = 0;
t     = linspace(t_0,t_end,N_t);

% -------------------------
% Discretize space:
% -------------------------
z_0   = 0;
L     = 79.3515;
if uniform == 1
    % --- Uniform Grid ---
    z    = linspace(z_0,L,Nz)';
else
    % --- Non-Uniform Grid ---
    alp = .75;
    bet = 4;
    xi  = linspace(0,1,Nz)';
    z   = L*(tanh(bet*(xi-alp))-tanh(-bet*alp))./(tanh(bet*(1-alp))-tanh(-bet*alp));
end

% -------------------------
% Crypt geometry factors:
% -------------------------
r  = r_b*(1 - exp(-a*z)) + r_t * exp(a*(z-L));
rz = a*r_b*exp(-a*z) + a*r_t*exp(a*(z-L));
g  = sqrt(1 + rz.^2);

% -------------------------
% Define the differentiation matrix:
% -------------------------
Dz = diffmat_nonuniform(z);

% -------------------------
% Initial condition:
% -------------------------
q_0     = ones(size(z));
q_0_int = q_0(2:end-1);

% -------------------------
% Build structures for solver:
% -------------------------
diffmat.Dz = Dz;

geom.z  = z;
geom.rz = rz;
geom.g  = g;

par.alpha = alpha;
par.zp    = zp;
par.Tc    = Tc;

% -------------------------
% Define right hand side function:
% -------------------------
dQdt = @(t,q_int) dqdt_1D_snipsnap(t, q_int, diffmat, geom, par, bcType);

% -------------------------
%Solve the system of ODEs that represents the PDE:
% -------------------------
tic
options = odeset('Stats', 'on');
[t,Q_int] = ode15s(dQdt, t, q_0_int, options);
toc

% -------------------------
% Extend to full Q:
% -------------------------
Q_full = zeros(Nz,N_t);
V_full = zeros(Nz,N_t);
wL     = full(Dz(1,1:3));
wR     = full(Dz(end,end-2:end));
for i = 1:N_t
    qi = Q_int(i,:).';   % interior as column

    switch bcType
        case "NbDt"
            Q_l = -(wL(2)*qi(1) + wL(3)*qi(2))/wL(1);
            Q_r = 1;

        case "NbNt"
            Q_l = -(wL(2)*qi(1) + wL(3)*qi(2))/wL(1);
            Q_r = -(wR(1)*qi(end-1) + wR(2)*qi(end))/wR(3);

        case "DbNt"
            Q_l = 1;
            Q_r = -(wR(1)*qi(end-1) + wR(2)*qi(end))/wR(3);

        case "DbDt"
            Q_l = 1;
            Q_r = 1;

        otherwise
            error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
    end

    Q_full(:,i) = [Q_l; qi; Q_r];
    V_full(:,i) = -(1./g).*par.alpha./Q_full(:,i).^3.*(Dz*Q_full(:,i));
end

% -------------------------
% Steady State Analysis:
% -------------------------
Q_ss = Q_full(:,end);
V_ss = V_full(:,end);

% -------------------------
% Total number of cells:
% -------------------------
Total_Cells  = 2*pi*trapz(z,r.*Q_ss.*g)
Prolif_Cells = 2*pi*trapz(z,r.*Q_ss.*(z<zp).*g)

% -------------------------
% Crypt renewal time:
% -------------------------
Arc_Length = cumtrapz(z,g);
pos = find(Arc_Length>=1);
Crypt_Renewal_Time = trapz(z(pos),g(pos)./V_ss(pos))/24

% -------------------------
% Steady state plot:
% -------------------------
while ssp == 1

    fig = figure(1);
    fig.WindowState = 'maximized';
    
    figure(1)
    subplot(1,5,1:3)
    plot(z,Q_ss,'k', 'LineWidth',10);
    hold on
    plot(z,V_ss,'k:','linewidth',10)
    set(gca,'fontsize',42)
    title('Cell Density and Velocity','fontsize',55,'interpreter','latex')
    xlabel('$z$','fontsize',50,'interpreter','latex')
    ylabel('$q_s(z)$ and $v_s(z)$','fontsize',50,'interpreter','latex')
    legend('$q_s(z)$','$v_s(z)$','fontsize',50,'interpreter','latex')
    grid on
    grid minor
    xlim([0 L])
    ylim([min(min([Q_ss; V_ss])) max(max([Q_ss; V_ss]))+.5])

    theta = linspace(0,2*pi,Nz);
    [T,R] = meshgrid(theta,r);

    X_plot = R.*cos(T);
    Y_plot = R.*sin(T);
    Z_plot = repmat(z,1,Nz);

    Q_dens = repmat(Q_full(:,end)',Nz,1);

    subplot(1,5,4:5)
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
    hold on
    mesh(1.2*R.*cos(T),1.2*R.*sin(T),zp*ones(size(T)),'edgecolor',[.5 .5 .5])
    legend('$q_s(z)$','fontsize',50,'interpreter','latex')
    xlim([-20 20])
    ylim([-20 20])
    zlim([0 80])
    caxis([1 max(max(Q_ss))])
    text(.725,.8,0,['$\mathcal{C} = $ ' num2str(ceil(Total_Cells))],...
        'fontsize',30,'Units','normalized','interpreter','latex')
    text(.725,.75,0,['$\mathcal{P} = $ ' num2str(ceil(Prolif_Cells))],...
        'fontsize',30,'Units','normalized','interpreter','latex')    
    text(.725,.7,0,['$\mathcal{R} = $ ' num2str(round(Crypt_Renewal_Time,2))],...
        'fontsize',30,'Units','normalized','interpreter','latex')

    ssp = 0;
end
%% Animation:
dt = round(.01*N_t);
for i = 1:dt:N_t

    Q = Q_full(:,i);
    V = V_full(:,i);

    figure(2)
    plot(z,Q, 'k', 'LineWidth',5);
    hold on
    plot(z,V,'k:','linewidth',5)
    set(gca,'fontsize',42)
    title('Cell Density and Velocity','fontsize',55,'interpreter','latex')
    xlabel('$z$','fontsize',50,'interpreter','latex')
    ylabel('$q(z,t)$ and $v(z,t)$','fontsize',50,'interpreter','latex')
    legend('$q(z,t)$','$v(z,t)$','fontsize',50,'interpreter','latex')
    grid on
    grid minor
    xlim([0 L])
    ylim([min(min([Q_ss; V_ss]))-.5 max(max([Q_ss; V_ss]))+.5])
    if i == 1
        pause
    end
    hold off;
end

%% Animation on crypt domain
dt = round(.1*N_t);
theta = linspace(0,2*pi,Nz);
[T,R] = meshgrid(theta,r);

X_plot = R.*cos(T);
Y_plot = R.*sin(T);
Z_plot = repmat(z,1,Nz);

for k = 1:dt:N_t

    Q_dens = repmat(Q_full(:,k)',Nz,1);

    figure(3)
    surf(X_plot,Y_plot,Z_plot,Q_dens')
    set(gca,'fontsize',42)
    title('Cell Density on Crypt','fontsize',55,'interpreter','latex')
    xlabel('$x$','fontsize',50,'interpreter','latex')
    ylabel('$y$','fontsize',50,'interpreter','latex')
    zlabel('$z$','fontsize',50,'interpreter','latex')
    legend('$q_s(z)$','fontsize',50,'interpreter','latex')
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
    caxis([1 max(max(Q_full))])
    if k == 1
        pause
    end

end

