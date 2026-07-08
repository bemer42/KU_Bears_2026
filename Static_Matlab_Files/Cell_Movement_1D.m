%% Cell movement PDE in one dimension along crypt shape
clear; close all; clc

% Predefine structures:
geom    = struct();
diffmat = struct();
par     = struct();

% Define BCs and dimesional pde:
bcType  = "NbDt";

% Define parameters:
alpha = 1.6406e3;
zp    = 27; 
Tc    = 15.1954;

% Define geometry parameters:
r_b   = 41/2/pi;
r_t   = 10/pi;
a     = 0.3;

% Discretize time:
Nt    = 7e3;
t_0   = 0;
t_end = 5e3;
t     = linspace(t_0,t_end,Nt);

% Discretize space:
Nz    = 1e3;
z_0   = 0;
L     = 78.8783;
z     = linspace(z_0, L, Nz)';
dz    = z(2) - z(1);

% Define geometry:
r  = r_b*(1 - exp(-a*z)) + r_t * exp(a*(z-L));
rz = a*r_b*exp(-a*z) + a*r_t*exp(a*(z-L));
g  = sqrt(1 + rz.^2);

% Define the differentiation matrix:
Dz = diag(1/2*ones(Nz-1,1),1) - diag(1/2*ones(Nz-1,1),-1);
Dz(1,1:3) = [-3/2 2 -1/2];
Dz(end, end-2:end) = [1/2 -2 3/2];
Dz = Dz/dz;

% Initial condition:
q_0     = ones(size(z));
q_0_int = q_0(2:end-1);

% Build structures for solver: 
diffmat.Dz = Dz;

geom.z  = z;
geom.rz = rz;
geom.g  = g; 

par.alpha = alpha;
par.zp    = zp; 
par.Tc    = Tc;

%Define right hand side function:
dQdt = @(t,q_int) dqdt_1D_snipsnap(t, q_int, diffmat, geom, par, bcType);

%Solve the system of ODEs that represents the PDE:
tic
options = odeset('Stats', 'on');
[t,Q_int] = ode15s(dQdt, t, q_0_int, options);
toc

% Extend to full Q:
Q_full = zeros(Nz,Nt);
V_full = zeros(Nz,Nt);
for i = 1:Nt

    switch bcType
        case "NbDt"
            % bottom Neumann
            Q_l = -(Dz(1,2)*Q_int(i,1) + Dz(1,3)*Q_int(i,2))/Dz(1,1);
            % top Dirichlet
            Q_r = 1;
        case "NbNt"
            % bottom Neumann
            Q_l = -(Dz(1,2)*Q_int(i,1) + Dz(1,3)*Q_int(i,2))/Dz(1,1);
            % top Neumann
            Q_r = -(Dz(end,end-1)*Q_int(i,end) + Dz(end,end-2)*Q_int(i,end-1))/Dz(end,end);
        case "DbNt"
            % bottom Dirichlet
            Q_l = 1;
            % top Neumann
            Q_r = -(Dz(end,end-1)*Q_int(i,end) + Dz(end,end-2)*Q_int(i,end-1))/Dz(end,end);
        case "DbDt"
            % bottom Dirichlet
            Q_l = 1;
            % top Dirichlet
            Q_r = 1;
        otherwise
            error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
    end

    Q_full(:,i) = [Q_l; Q_int(i, :)'; Q_r];
    V_full(:,i) = -(1./g).*alpha./Q_full(:,i).^3.*(Dz*Q_full(:,i));

end

% Steady State Analysis:
Q_ss = Q_full(:,end);
V_ss = V_full(:,end);

% Total number of cells:
Total_Cells = 2*pi*trapz(z,r.*Q_ss.*g)

% Crypt renewal time:
Arc_Length = cumtrapz(z,g);
pos = find(Arc_Length>=1);
Crypt_Renewal_Time = trapz(z(pos),g(pos)./V_ss(pos))/24

% Steady state plot:
figure(1)
plot(z,Q_ss,'k', 'LineWidth',5);
hold on
plot(z,V_ss,'k:','linewidth',5)
set(gca,'fontsize',16)
title('Cell Density and Velocity','fontsize',25,'interpreter','latex')
xlabel('$z$','fontsize',20,'interpreter','latex')
ylabel('$q(z,t)$ and $v(z,t)$','fontsize',20,'interpreter','latex')
legend('$q(z,t)$','fontsize',18,'interpreter','latex')
grid on
grid minor
xlim([0 L])
ylim([min(min([Q_ss; V_ss]))-.5 max(max([Q_ss; V_ss]))+.5])


%% Animation:
dt = round(.01*Nt);
for i = 1:dt:Nt

    Q = Q_full(:,i);
    V = V_full(:,i);

    figure(2)
    plot(z,Q, 'k', 'LineWidth',5);
    hold on
    plot(z,V,'k:','linewidth',5)
    set(gca,'fontsize',16)
    title('Cell Density and Velocity','fontsize',25)
    xlabel('z','fontsize',20)
    ylabel('q(z,t) and v(z,t)','fontsize',20)
    legend('q(z,t)','v(z,t)','fontsize',18)
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
dt = round(.1*Nt);
theta = linspace(0,2*pi,Nz);
r_plot = r(z);
[T,R] = meshgrid(theta,r_plot);

X_plot = R.*cos(T);
Y_plot = R.*sin(T);
Z_plot = repmat(z,1,Nz);

for k = 1:dt:Nt

    Q_dens = repmat(Q_full(:,k)',Nz,1);

    figure(3)
    surf(X_plot,Y_plot,Z_plot,Q_dens')
    shading interp
    colormap turbo
    colorbar
    xlim([-30 30])
    ylim([-30 30])
    zlim([0 80])
    caxis([1 1.5])
    if k == 1
        pause
    end

end

