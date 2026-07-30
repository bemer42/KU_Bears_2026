%% Cell Movement 2D Nonlinear
close all; clear; clc

% Predefine structures:
geom    = struct();
diffmat = struct();
par     = struct();

% Boundary condition set:
bcType = "NbDt";

% Define parameters:
alpha_z  = 1.6406e3;
alpha_th = alpha_z/10;
Tc       = 15.1954;
zp       = 27;

% Define geometry parameters:
r_b   = 41/2/pi;
r_t   = 10/pi;
a     = 0.3;

% Discretize time:
N_t = 2e3;
t_0 = 0;
t_end = 8.5e-1;
t = linspace(t_0,t_end,N_t);

% Discretize theta space:
Nth    = 4e1;
th_0   = 0;
th_end = 2*pi;
th     = linspace(th_0,th_end,Nth)';
dth    = th(2) - th(1);

% Discretize z space:
Nz   = 7e1;
z_0  = 0;
L    = 78.8783;
% --- Uniform Grid ---
z    = linspace(z_0,L,Nz)';
dz   = z(2) - z(1);
% --- Non-Uniform Grid ---
beta = .4;
xi   = linspace(0,1,Nz)';
z    = L/2*(1+tanh(beta*(2*xi-1))/tanh(beta));

% Create mesh grid:
[Z,Th] = meshgrid(z,th);

% Differentiation Matrix in theta
Dth = spdiags([-1/2*ones(Nth,1) 1/2*ones(Nth,1)],[-1 1],Nth,Nth);
Dth(1,end) = -1/2;
Dth(end,1) =  1/2;
Dth = Dth/dth;

% Define the differentiation matrix in z
% --- Uniform Grid ---
Dz = spdiags([-1/2*ones(Nz,1) 1/2*ones(Nz,1)],[-1 1],Nz,Nz);
Dz(1,1:3) = [-3/2 2 -1/2];
Dz(end, end-2:end) = [1/2 -2 3/2];
Dz = Dz/dz;
% --- Non-Uniform Grid ---
Dz = diffmat_nonuniform(z);

R   = r_b*(1-exp(-a*Z))+r_t*exp(a*(Z-L));
Rth = 0*Th;
Rz  = a*r_b*exp(-a*Z) + a*r_t*exp(a*(Z-L));
G   = sqrt(Rth.^2+R.^2.*(Rz.^2 + 1));

% Initial Condition:
f = @(th,z) 1 + exp(-(z - zp).^2/.1-(th-pi/2).^2/.2);
% f = @(th,z) ones(size(th));
Q0 = f(Th,Z);
Q0_int = Q0(:,2:end-1);
q0_int = Q0_int(:);

% Construct structures for pde:
diffmat.Dth = Dth;
diffmat.Dz  = Dz;

geom.Th  = Th;
geom.Z   = Z;
geom.R   = R;
geom.Rth = Rth;
geom.Rz  = Rz;
geom.G   = G;

par.alpha_th = alpha_th;
par.alpha_z  = alpha_z;
par.Tc       = Tc;
par.zp       = zp;

% Define the right hand side
dQdt = @(t, q_int) dqdt_2D_snipsnap(t, q_int, diffmat, geom, par, bcType);

% Solve the PDE:
tic
% dtdiff = 0.95*(min(diff(z))^2/par.alpha_z);
dtdiff = inf;
options = odeset('Stats','on','MaxStep',dtdiff);
[t, q_int] = ode15s(dQdt, t, q0_int, options);
toc

% Expand interior into full for each time point:
Q_full = zeros(N_t,Nth*Nz);

% One-sided z-derivative stencil coefficients from Dz
b0 = Dz(1,1);       b1 = Dz(1,2);       b2 = Dz(1,3);
t2 = Dz(end,end-2); t1 = Dz(end,end-1); t0 = Dz(end,end);

for i = 1:N_t

    % Reshape q_int into matrix:
    Q_int = reshape(q_int(i,:), Nth, Nz-2);

    % Initialize full snapshot and insert interior
    Q = zeros(Nth, Nz);
    Q(:,2:Nz-1) = Q_int;

    % Apply BCs based on bcType
    switch bcType
        case "NbDt"
            % Bottom Neumann
            Cz  = spdiags(par.alpha_z*(Rth(:,1).^2 + R(:,1).^2), 0, Nth, Nth);
            Cth = spdiags(par.alpha_th*(Rth(:,1).*Rz(:,1)),      0, Nth, Nth);
            A = b0*Cz - Cth*Dth;
            b = -Cz*(b1*Q(:,2) + b2*Q(:,3));
            Q(:,1) = A \ b;
            % Top Dirichlet
            Q(:,end) = 1;
        case "NbNt"
            % Bottom Neumann
            Cz  = spdiags(par.alpha_z*(Rth(:,1).^2 + R(:,1).^2), 0, Nth, Nth);
            Cth = spdiags(par.alpha_th*(Rth(:,1).*Rz(:,1)),      0, Nth, Nth);
            A = b0*Cz - Cth*Dth;
            b = -Cz*(b1*Q(:,2) + b2*Q(:,3));
            Q(:,1) = A \ b;
            % Top Neumann
            Cz  = spdiags(par.alpha_z*(Rth(:,Nz).^2 + R(:,Nz).^2), 0, Nth, Nth);
            Cth = spdiags(par.alpha_th*(Rth(:,Nz).*Rz(:,Nz)),      0, Nth, Nth);
            A = t0*Cz - Cth*Dth;
            b = -Cz*(t1*Q(:,end-1) + t2*Q(:,end-2));
            Q(:,end) = A \ b;
        case "DbNt"
            % Bottom Dirichlet
            Q(:,1) = 1;
            % Top Neumann (Flux_z = 0): solve Q(:,end)
            Cz  = spdiags(par.alpha_z*(Rth(:,Nz).^2 + R(:,Nz).^2), 0, Nth, Nth);
            Cth = spdiags(par.alpha_th*(Rth(:,Nz).*Rz(:,Nz)),      0, Nth, Nth);
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
    % Save flattened snapshot
    Q_full(i,:) = Q(:);
end

% Steady State Analysis
Q_ss = reshape(Q_full(end,:),Nth,Nz);

% Total number of cells:
Total_Cells = trapz(z,trapz(th,Q_ss.*G))
Prolif_Cells = trapz(z,trapz(th,Q_ss.*(Z<zp).*G))

[Vth, Vz] = velocityField_FVconsistent_2D(Q_ss, diffmat, geom, par, bcType);

% Crypt Renewal Time
[Rmean, Rall, theta0] = avgTravelTimeFromVectorField(th, z, Vth, Vz, ...
    'theta0', linspace(pi/2,3*pi/2,60), 'z0', z(2));
Crypt_Renewal_Time = Rmean/24


% Steady state plot:
dm = 5;
figure(1)
surf(Th,Z,Q_ss);
hold on
quiver3(Th(2:dm:end,2:dm:end-dm),Z(2:dm:end,2:dm:end-dm),Q_ss(2:dm:end,2:dm:end-dm),...
    Vth(2:dm:end,2:dm:end-dm),Vz(2:dm:end,2:dm:end-dm),zeros(size(Th(2:dm:end,2:dm:end-dm))),...
    1.0,'color',[0.15 0.15 0.15],'Linewidth',2.5,'MaxHeadSize',.008);
set(gca,'fontsize',42)
title('Cell Density and Velocity','fontsize',55,'interpreter','latex')
xlabel('$\theta$','fontsize',50,'interpreter','latex')
ylabel('$z$','fontsize',50,'interpreter','latex')
zlabel('$q_s(\theta,z)$ and $\vec{v}_s(\theta,z)$','fontsize',50,'interpreter','latex')
legend('$q_s(\theta,z)$','$\vec{v}_s(\theta,z)$','fontsize',50,'interpreter','latex')
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


%% Animation
dm = 5;
dt  = round(.05*N_t);

for i = 1:dt:N_t

    Q = reshape(Q_full(i,:), Nth, Nz);

    [Vth, Vz] = velocityField_FVconsistent_2D(Q, diffmat, geom, par, bcType);

    figure(2)
    surf(Th,Z,Q);
    hold on
    quiver3(Th(2:dm:end,2:dm:end-dm),Z(2:dm:end,2:dm:end-dm),Q(2:dm:end,2:dm:end-dm),...
        Vth(2:dm:end,2:dm:end-dm),Vz(2:dm:end,2:dm:end-dm),zeros(size(Th(2:dm:end,2:dm:end-dm))),...
        1.0,'color',[0.15 0.15 0.15],'Linewidth',2.5,'MaxHeadSize',.008);    set(gca,'fontsize',42)
    title('Cell Density and Velocity','fontsize',55,'interpreter','latex')
    xlabel('$\theta$','fontsize',50,'interpreter','latex')
    ylabel('$z$','fontsize',50,'interpreter','latex')
    zlabel('$q_s(\theta,z)$ and $\vec{v}_s(\theta,z)$','fontsize',50,'interpreter','latex')
    legend('$q_s(\theta,z)$','$\vec{v}_s(\theta,z)$','fontsize',50,'interpreter','latex')
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
    if i == 1
        pause
    end
    hold off
end

%% Plot on crypt domain
dt = round(.1*N_t);


X_plot = R.*cos(Th);
Y_plot = R.*sin(Th);
Z_plot = Z;

for k = 1:dt:N_t

    Q = reshape(Q_full(k,:),Nth,Nz);

    figure(2)
    surf(X_plot,Y_plot,Z_plot,Q)
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
    caxis([1 1.05])
    if k == 1
        pause
    end

end


