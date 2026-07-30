%% Cell Movement 2D Dynamic
close all; clear; clc; ssp = 1;

% Grid/Time Points:
N_t   = 1e3;
t_end = 5e2;
Nz    = 3e1;
Nth   = 3e1;

% Predefine structures:
geom    = struct();
diffmat = struct();
par     = struct();

% Boundary condition set:
bcType = "NbDt";

% Define cell density parameters
alpha_z  = 1.3e3;
alpha_th = 1.3e1;
Tc       = 15.1954;
zp       = 27;

% Define signal parameters:
alpha_s = 1e-4;
beta    = 1e-2;
mu      = 1e-2;
q_th    = 1.04;
epsQ    = 0.01;

% Define radius parameters:
eta   = 1e-2;
s_th  = 0.5;
epsS  = 0.001;

% Define geometry parameters
r_b = 41/2/pi;
r_t = 10/pi;
a   = 0.3;

% Discretize time
t_0   = 0;
t     = linspace(t_0,t_end,N_t);

% Discretize theta space
th_0    = 0;
th_end  = 2*pi;
th      = linspace(0,2*pi,Nth+1)';
th(end) = [];
dth     = th(2)-th(1);

% Discretize z space (nonuniform)
z_0  = 0;
L    = 78.8783;
b    = 2;
xi   = linspace(0,1,Nz)';
z    = L/2*(1+tanh(b*(2*xi-1))/tanh(b));

% Create mesh grid:
[Z,Th] = meshgrid(z,th);

% Differentiation Matrix in theta
Dth        = spdiags([-1/2*ones(Nth,1) 1/2*ones(Nth,1)],[-1 1],Nth,Nth);
Dth(1,end) = -1/2;
Dth(end,1) =  1/2;
Dth        = Dth/dth;

% Differentiation Matrix in z
Dz = diffmat_nonuniform(z);

% FV data for nonuniform z-grid
zv  = z(:);
dzf = diff(zv);
zf  = [zv(1); 0.5*(zv(1:end-1)+zv(2:end)); zv(end)];
dzc = diff(zf);

% Initial radius r0 and target radius r_tar
r0z = r_b*(1-exp(-a*z)) + r_t*exp(a*(z-L));
R0  = repmat(r0z.', Nth, 1);
bulge = exp(-((Th - 3*pi/2).^2)/0.2) .* exp(-((Z - 30).^2)/50);
Rtar = R0 .* (1 + 0.5*bulge);
Rtar = R0 .* (2);

% Initial conditions
f = @(th,z) 1 + exp(-(z - zp).^2/.1-(th-3*pi/2).^2/.2);
f = @(th,z) ones(size(th));
Q0 = f(Th,Z);
S0 = zeros(Nth,Nz);

Q0_int  = Q0(:,2:end-1);
S0_int  = S0(:,2:end-1);
R0_full = R0;

y0 = [Q0_int(:); S0_int(:); R0_full(:)];

% =========================
% Store structures for PDE
diffmat.Dth  = Dth;
diffmat.Dz   = Dz;

geom.Th   = Th;
geom.Z    = Z;
geom.z    = z;
geom.th   = th;
geom.dth  = dth;
geom.dzf  = dzf;
geom.dzc  = dzc;
geom.L    = L;
geom.R0   = R0;
geom.Rtar = Rtar;

par.alpha_z  = alpha_z;
par.alpha_th = alpha_th;
par.Tc       = Tc;
par.zp       = zp;
par.alpha_s  = alpha_s;
par.beta     = beta;
par.mu       = mu;
par.q_th     = q_th;
par.epsQ     = epsQ;
par.eta      = eta;
par.s_th     = s_th;
par.epsS     = epsS;
par.r_b      = r_b;
par.r_t      = r_t;
par.a        = a;
% =========================

% Define RHS and solve
dYdt = @(t, y) dqdt_dyn_2D_snipsnap(t, y, diffmat, geom, par, bcType);

tic
options = odeset('Stats','on','MaxStep',inf);
[t, y] = ode15s(dYdt, t, y0, options);
toc

% Expand interior into full snapshots for each time point:
nq = Nth*(Nz-2);
ns = Nth*(Nz-2);
nr = Nth*Nz;

Q_full   = zeros(N_t, Nth*Nz);
S_full   = zeros(N_t, Nth*Nz);
R_full   = zeros(N_t, Nth*Nz);
Vth_full = zeros(N_t, Nth*Nz);
Vz_full  = zeros(N_t, Nth*Nz);

% One-sided z-derivative stencil coefficients from Dz
b0 = Dz(1,1);       b1 = Dz(1,2);       b2 = Dz(1,3);
t2 = Dz(end,end-2); t1 = Dz(end,end-1); t0 = Dz(end,end);

for i = 1:N_t

    % Unpack time i:
    q_int = y(i, 1:nq).';
    s_int = y(i, nq+(1:ns)).';
    r_vec = y(i, nq+ns+(1:nr)).';

    Q_int = reshape(q_int, Nth, Nz-2);
    S_int = reshape(s_int, Nth, Nz-2);
    R     = reshape(r_vec, Nth, Nz);

    % Geometry factors for this snapshot
    Rth = Dth*R;
    Rz  = R*Dz';
    G   = sqrt(Rth.^2 + R.^2.*(Rz.^2 + 1));

    % Reassemble Q and S full
    Q = zeros(Nth, Nz);
    S = zeros(Nth, Nz);
    Q(:,2:Nz-1) = Q_int;
    S(:,2:Nz-1) = S_int;

    % ---- Apply BCs for q based on bcType (same as your original code)
    switch bcType
        case "NbDt"
            % Bottom Neumann (Flux_z=0): solve Q(:,1)
            Cz  = spdiags(par.alpha_z*(Rth(:,1).^2 + R(:,1).^2), 0, Nth, Nth);
            Cth = spdiags(par.alpha_z*(Rth(:,1).*Rz(:,1)),      0, Nth, Nth);
            A = b0*Cz - Cth*Dth;
            b = -Cz*(b1*Q(:,2) + b2*Q(:,3));
            Q(:,1) = A \ b;

            % Top Dirichlet
            Q(:,end) = 1;

        case "NbNt"
            % Bottom Neumann
            Cz  = spdiags(par.alpha_z*(Rth(:,1).^2 + R(:,1).^2), 0, Nth, Nth);
            Cth = spdiags(par.alpha_z*(Rth(:,1).*Rz(:,1)),      0, Nth, Nth);
            A = b0*Cz - Cth*Dth;
            b = -Cz*(b1*Q(:,2) + b2*Q(:,3));
            Q(:,1) = A \ b;

            % Top Neumann
            Cz  = spdiags(par.alpha_z*(Rth(:,Nz).^2 + R(:,Nz).^2), 0, Nth, Nth);
            Cth = spdiags(par.alpha_z*(Rth(:,Nz).*Rz(:,Nz)),      0, Nth, Nth);
            A = t0*Cz - Cth*Dth;
            b = -Cz*(t1*Q(:,end-1) + t2*Q(:,end-2));
            Q(:,end) = A \ b;

        case "DbNt"
            % Bottom Dirichlet
            Q(:,1) = 1;

            % Top Neumann
            Cz  = spdiags(par.alpha_z*(Rth(:,Nz).^2 + R(:,Nz).^2), 0, Nth, Nth);
            Cth = spdiags(par.alpha_z*(Rth(:,Nz).*Rz(:,Nz)),      0, Nth, Nth);
            A = t0*Cz - Cth*Dth;
            b = -Cz*(t1*Q(:,end-1) + t2*Q(:,end-2));
            Q(:,end) = A \ b;

        case "DbDt"
            % Bottom Dirichlet
            Q(:,1) = 1;

            % Top Dirichlet
            Q(:,end) = 1;

        otherwise
            error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
    end

    % ---- Apply BCs for s: Neumann in z (s_z=0) using one-sided Dz rows
    S(:,1)   = -(b1*S(:,2) + b2*S(:,3))/b0;
    S(:,end) = -(t1*S(:,end-1) + t2*S(:,end-2))/t0;

    % ---- Store full snapshots
    Q_full(i,:) = Q(:);
    S_full(i,:) = S(:);
    R_full(i,:) = R(:);

    % ---- Reassemble velocity fields (FV-consistent with your dqdt code)
    geom_i     = geom;
    geom_i.R   = R;
    geom_i.Rth = Rth;
    geom_i.Rz  = Rz;
    geom_i.G   = G;

    [Vth, Vz] = velocityField_FVconsistent_2D(Q, diffmat, geom_i, par, bcType);

    Vth_full(i,:) = Vth(:);
    Vz_full(i,:)  = Vz(:);

end

% Steady State Analysis
Q_ss   = reshape(Q_full(end,:), Nth, Nz);
S_ss   = reshape(S_full(end,:), Nth, Nz);
R_ss   = reshape(R_full(end,:), Nth, Nz);
Vth_ss = reshape(Vth_full(end,:), Nth, Nz);
Vz_ss  = reshape(Vz_full(end,:), Nth, Nz);

% Geometry factor at final time
Rth = Dth*R_ss;
Rz  = R_ss*Dz';
G   = sqrt(Rth.^2 + R_ss.^2.*(Rz.^2 + 1));

% Total number of cells:
Total_Cells  = trapz(z,trapz(th,Q_ss.*G))
Prolif_Cells = trapz(z,trapz(th,Q_ss.*(Z<par.zp).*G))

% Crypt Renewal Time (reuse your functions)
% r_val   = @(z) r_b*(1-exp(-a*z))+r_t*exp(a*(z-L));
% z_val = fzero(@(z) 2*pi*r_val(z)-1,.4);
% ind = find(z>z_val);
%
% [Rmean, Rall, theta0] = avgTravelTimeFromVectorField(th, z, Vth_ss, Vz_ss, ...
%     'theta0', linspace(pi/2,3*pi/2,60), 'z0', z(ind(1)));
% Crypt_Renewal_Time = Rmean/24

% Steady state plot:
while ssp == 1
    ex = @(x) [x; x(1,:)];
    dm = 1;
    figure(1)
    surf([Th; 2*pi*ones(1,Nz)],ex(Z),ex(Q_ss));
    hold on
    quiver3(Th(2:dm:end,2:dm:end-dm),Z(2:dm:end,2:dm:end-dm),Q_ss(2:dm:end,2:dm:end-dm),...
        Vth_ss(2:dm:end,2:dm:end-dm),Vz_ss(2:dm:end,2:dm:end-dm),zeros(size(Th(2:dm:end,2:dm:end-dm))),...
        1.0,'color',[0.15 0.15 0.15],'Linewidth',2.5,'MaxHeadSize',.008);
    set(gca,'fontsize',42)
    title('Cell Density and Velocity','fontsize',55,'interpreter','latex')
    xlabel('$\theta$','fontsize',50,'interpreter','latex')
    ylabel('$z$','fontsize',50,'interpreter','latex')
    zlabel('$q(\theta,z)$ and $\vec{v}(\theta,z)$','fontsize',50,'interpreter','latex')
    legend('$q(\theta,z)$','$\vec{v}(\theta,z)$','fontsize',50,'interpreter','latex')
    grid on
    grid minor
    box on
    colormap turbo
    shading interp
    colorbar
    caxis([1 max(max(Q_ss))])
    xlim([0 2*pi])
    ylim([0 L])
    view([0 0 1])
    hold off

    % Intial, Time End, and Target Domain:
    X0_plot   = R0.*cos(Th);
    Y0_plot   = R0.*sin(Th);
    Xend_plot = R_ss.*cos(Th);
    Yend_plot = R_ss.*sin(Th);
    Xtar_plot = Rtar.*cos(Th);
    Ytar_plot = Rtar.*sin(Th);    
    
    figure(2)
    subplot(1,3,1)
    mesh(ex(X0_plot),ex(Y0_plot),ex(Z),'EdgeColor',[.5 .5 .5])
    set(gca,'fontsize',42)
    title('Intial Crypt State','fontsize',55,'interpreter','latex')
    xlabel('$x$','fontsize',50,'interpreter','latex')
    ylabel('$y$','fontsize',50,'interpreter','latex')
    zlabel('$z$','fontsize',50,'interpreter','latex')
    grid on
    grid minor
    box on
    lighting gouraud
    xlim([-40 40])
    ylim([-40 40])
    zlim([0 80])
    subplot(1,3,2)
    mesh(ex(Xend_plot),ex(Yend_plot),ex(Z),'EdgeColor',[.5 .5 .5])
    set(gca,'fontsize',42)
    title('Crypt State at Time End','fontsize',55,'interpreter','latex')
    xlabel('$x$','fontsize',50,'interpreter','latex')
    ylabel('$y$','fontsize',50,'interpreter','latex')
    zlabel('$z$','fontsize',50,'interpreter','latex')
    grid on
    grid minor
    box on
    lighting gouraud
    xlim([-40 40])
    ylim([-40 40])
    zlim([0 80])
    subplot(1,3,3)
    mesh(ex(Xtar_plot),ex(Ytar_plot),ex(Z),'EdgeColor',[.5 .5 .5])
    set(gca,'fontsize',42)
    title('Target Crypt State','fontsize',55,'interpreter','latex')
    xlabel('$x$','fontsize',50,'interpreter','latex')
    ylabel('$y$','fontsize',50,'interpreter','latex')
    zlabel('$z$','fontsize',50,'interpreter','latex')
    grid on
    grid minor
    box on
    lighting gouraud
    xlim([-40 40])
    ylim([-40 40])
    zlim([0 80])

    ssp = 0;
end


%% Animation
dm = 1;
dt  = ceil(.01*N_t);
for i = 1:dt:N_t

    Q  = reshape(Q_full(i,:),  Nth, Nz);
    R  = reshape(R_full(i,:), Nth, Nz);
    S  = reshape(S_full(i,:), Nth, Nz);
    Vth = reshape(Vth_full(i,:),Nth, Nz);
    Vz  = reshape(Vz_full(i,:), Nth, Nz);

    figure(3)
    subplot(2,2,[1 3])
    surf([Th; 2*pi*ones(1,Nz)],ex(Z),ex(Q));
    hold on
    quiver3(Th(2:dm:end,2:dm:end-dm),Z(2:dm:end,2:dm:end-dm),Q(2:dm:end,2:dm:end-dm),...
        Vth(2:dm:end,2:dm:end-dm),Vz(2:dm:end,2:dm:end-dm),zeros(size(Th(2:dm:end,2:dm:end-dm))),...
        1.0,'color',[0.15 0.15 0.15],'Linewidth',2.5,'MaxHeadSize',.008);
    set(gca,'fontsize',42)
    title('Cell Density and Velocity','fontsize',55,'interpreter','latex')
    xlabel('$\theta$','fontsize',50,'interpreter','latex')
    ylabel('$z$','fontsize',50,'interpreter','latex')
    zlabel('$q(\theta,z)$ and $\vec{v}(\theta,z)$','fontsize',50,'interpreter','latex')
    legend('$q(\theta,z)$','$\vec{v}(\theta,z)$','fontsize',50,...
           'interpreter','latex','location','northwest')
    grid on
    grid minor
    box on
    colormap turbo
    shading interp
    colorbar
    caxis([1 max(max(Q_ss))])
    xlim([0 2*pi])
    ylim([0 L])
    view([0 0 1])
    hold off
    subplot(2,2,2)
    mesh(ex(R).*cos([Th; 2*pi*ones(1,Nz)]),...
         ex(R).*sin([Th; 2*pi*ones(1,Nz)]),ex(Z),'edgecolor',[0.5 0.5 0.5])
    set(gca,'fontsize',42)
    title('Crypt State','fontsize',55,'interpreter','latex')
    xlabel('$x$','fontsize',50,'interpreter','latex')
    ylabel('$y$','fontsize',50,'interpreter','latex')
    zlabel('$z$','fontsize',50,'interpreter','latex')
    legend('$\vec{X}(\theta,z,t)$','fontsize',50,'interpreter','latex')
    xlim([-40 40])
    ylim([-40 40])
    zlim([0 80])
    hold off
    subplot(2,2,4)
    surf(Th,Z,S)
    hold on
    mesh([Th; 2*pi*ones(1,Nz)],ex(Z),...
        s_th*ones(size(ex(Z))),'edgecolor',[0.5 0.5 0.5])
    set(gca,'fontsize',42)
    title('Signal','fontsize',55,'interpreter','latex')
    xlabel('$\theta$','fontsize',50,'interpreter','latex')
    ylabel('$z$','fontsize',50,'interpreter','latex')
    zlabel('$s(\theta,z,t)$','fontsize',50,'interpreter','latex')
    legend('$s(\theta,z,t)$','$s_{th}$','fontsize',50,'interpreter','latex')
    zlim([0 beta/mu])
    hold off
    if i == 1
        pause
    end
end

%% Animation plot on crypt domain
dt = ceil(.01*N_t);

Z_plot = Z;

for k = 1:dt:N_t

    Q = reshape(Q_full(k,:), Nth, Nz);
    R = reshape(R_full(k,:), Nth, Nz);

    X_plot = ex(R).*cos([Th; 2*pi*ones(1,Nz)]);
    Y_plot = ex(R).*sin([Th; 2*pi*ones(1,Nz)]);

    figure(4)
    surf(X_plot,Y_plot,ex(Z_plot),ex(Q))
    set(gca,'fontsize',42)
    title('Cell Density on Crypt','fontsize',55,'interpreter','latex')
    xlabel('$x$','fontsize',50,'interpreter','latex')
    ylabel('$y$','fontsize',50,'interpreter','latex')
    zlabel('$z$','fontsize',50,'interpreter','latex')
    legend('$q(\theta,z)$','fontsize',50,'interpreter','latex')
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
    caxis([1 1.05])
    if k == 1
        pause
    end

end