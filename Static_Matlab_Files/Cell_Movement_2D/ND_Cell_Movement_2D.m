%% Nondimensional Cell Movement 2D
close all; clear; clc

% Predefine structures:
geom_hat    = struct();
diffmat_hat = struct();
par_hat     = struct();

% Boundary condition set:
bcType = "NbDt";

% -------------------------
%  Dimensional parameters
% -------------------------
alpha_z  = 1.6406e3;
alpha_th = alpha_z/10;
Tc       = 15.1954;
zp       = 27;
L        = 78.8783;
ell      = 1;

% -------------------------
% Geometry parameters:
% -------------------------
r_b   = 41/2/pi;
r_t   = 10/pi;
a     = 0.3;

% -------------------------
%  Nondimensional groups
% -------------------------
epsilon  = r_b/L;
gamma_th = log(2)*r_b^2/(alpha_th*ell^2*Tc);
gamma_z  = log(2)*L^2  /(alpha_z *ell^2*Tc);
Zp_hat   = zp/L;

% -------------------------
%  Time discretization
% -------------------------
Nt        = 7e3;
t_hat_0   = 0;
t_hat_end = 5e1*log(2)./Tc;
t_hat     = linspace(t_hat_0,t_hat_end,Nt);

% -------------------------
%  Theta grid
% -------------------------
Nth    = 5e1;
th_0   = 0;
th_end = 2*pi;
th     = linspace(th_0, th_end, Nth)';
dth    = th(2) - th(1);

% -------------------------
%  z-hat grid in [0,1]
% -------------------------
Nz     = 5e1;
beta   = 3;
xi     = linspace(0,1,Nz)';
z_hat  = 0.5*(1 + tanh(beta*(2*xi-1))/tanh(beta));
z      = L*z_hat;
dz_hat = z_hat(2) - z_hat(1);

% Create mesh grids
[Z_hat, Th] = meshgrid(z_hat, th);
Z           = L*Z_hat;

% -------------------------
%  Differentiation in theta
% -------------------------
Dth = spdiags([-1/2*ones(Nth,1) 1/2*ones(Nth,1)],[-1 1],Nth,Nth);
Dth(1,end) = -1/2;
Dth(end,1) =  1/2;
Dth = Dth/dth;

% -------------------------
%  Differentiation in z
% -------------------------
Dz = diffmat_nonuniform(z_hat);

% -------------------------
%  Nondimensional geometry on (theta, zhat)
% -------------------------
R_hat   = (1 - exp(-a*L*Z_hat)) + (r_t/r_b)*exp(a*L*(Z_hat - 1));
Rth_hat = 0*Th;
Rz_hat  = a*L*exp(-a*L*Z_hat) + (r_t/r_b)*a*L*exp(a*L*(Z_hat - 1));
S_hat   = sqrt(Rth_hat.^2 + R_hat.^2 .* (1 + epsilon^2*(Rz_hat.^2)));

% physical radius for plotting on crypt:
Rphys = r_b * R_hat;
Gphys = sqrt( (r_b*Rth_hat).^2 + (Rphys).^2 .* (1 + (epsilon*Rz_hat).^2) );

% -------------------------
%  Structures for RHS:
% -------------------------
diffmat_hat.Dth  = Dth;
diffmat_hat.Dz   = Dz;

geom_hat.Th      = Th;
geom_hat.Z_hat   = Z_hat;
geom_hat.Zphys   = Z;
geom_hat.R_hat   = R_hat;
geom_hat.Rth_hat = Rth_hat;
geom_hat.Rz_hat  = Rz_hat;
geom_hat.S_hat   = S_hat;
geom_hat.Rphys   = Rphys;
geom_hat.Gphys   = Gphys;

par_hat.gamma_th = gamma_th;
par_hat.gamma_z  = gamma_z;
par_hat.epsilon  = epsilon;
par_hat.ell      = ell;
par_hat.Zp_hat   = Zp_hat;

% -------------------------
%  Initial condition
% % -------------------------
% f = @(th, z_hat) 1 + exp(-(z_hat - Zp_hat).^2/.1 - (th - pi/2).^2/.2);
f = @(th,z_hat) ones(size(Th));
Q_hat_0 = f(Th, Z_hat);
Q_hat_0_int = Q_hat_0(:,2:end-1);
q_hat_0_int = Q_hat_0_int(:);

% -------------------------
%  ND RHS handle and solve
% -------------------------
dQhatdt = @(t_hat, q_hat_int) ND_dqdt_2D_snipsnap(t_hat, q_hat_int, diffmat_hat, geom_hat, par_hat, bcType);

tic
options = odeset('Stats','on','MaxStep',inf);
[t_hat, q_hat_int] = ode15s(dQhatdt, t_hat, q_hat_0_int, options);
toc

% -------------------------
%  Expand interior back to full Q for each time point 
% -------------------------
Qhat_full = zeros(Nt, Nth*Nz);

% One-sided zhat-derivative stencil coefficients from Dz
b0 = Dz(1,1);       b1 = Dz(1,2);       b2 = Dz(1,3);
t2 = Dz(end,end-2); t1 = Dz(end,end-1); t0 = Dz(end,end);

nz = Nz - 2;

for i = 1:Nt

    Q_int = reshape(q_hat_int(i,:), Nth, nz);

    Q = zeros(Nth, Nz);
    Q(:,2:Nz-1) = Q_int;

    % Apply BCs (ND IBVP: usually bottom Neumann, top Dirichlet qhat=ell)
    switch bcType
        case "NbDt"
            Q(:,1)   = solveFluxNeumannAtBoundary(1,  b0, b1, b2, Q(:,2),Q(:,3),diffmat_hat, geom_hat);
            Q(:,end) = ell;
        case "NbNt"
            Q(:,1)   = solveFluxNeumannAtBoundary(1,  b0, b1, b2, Q(:,2),Q(:,3),diffmat_hat, geom_hat);
            Q(:,end) = solveFluxNeumannAtBoundary(Nz, t0, t1, t2, Q(:,end-1), Q(:,end-2),diffmat_hat, geom_hat);
        case "DbNt"
            Q(:,1)   = ell;  % if you really want Dirichlet at the base in ND units
            Q(:,end) = solveFluxNeumannAtBoundary(Nz, t0, t1, t2, Q(:,end-1), Q(:,end-2),diffmat_hat, geom_hat);
        case "DbDt"
            Q(:,1)   = ell;
            Q(:,end) = ell;
        otherwise
            error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
    end

    Qhat_full(i,:) = Q(:);
end

% -------------------------
%  Steady state analysis 
% -------------------------
Qhat_ss = reshape(Qhat_full(end,:), Nth, Nz);

Total_Cells  = (r_b*L/ell)*trapz(z_hat, trapz(th, Qhat_ss .* S_hat))
Prolif_Cells = (r_b*L/ell)*trapz(z_hat, trapz(th, Qhat_ss .* (Z_hat < Zp_hat) .* S_hat))

% -------------------------
%  Velocity field (ND) for plotting arrows
% -------------------------
[Vth_ss, Vz_hat_ss] = ND_velocityField_FVconsistent_2D(Qhat_ss, diffmat_hat, geom_hat, par_hat, bcType);

% Crypt Renewal Time
[Rmean, Rall, theta0] = avgTravelTimeFromVectorField(th, z_hat, Vth_ss, Vz_hat_ss, ...
    'theta0', linspace(pi/2,3*pi/2,60), 'z0', z(2));
Crypt_Renewal_Time = (Tc/log(2))*Rmean/24

% -------------------------
%  Steady state plot
% -------------------------
dm = 5;
figure(1)
surf(Th, Z_hat, Qhat_ss);
hold on
quiver3(Th(2:dm:end,2:dm:end-dm), Z_hat(2:dm:end,2:dm:end-dm), Qhat_ss(2:dm:end,2:dm:end-dm), ...
        Vth_ss(2:dm:end,2:dm:end-dm), Vz_hat_ss(2:dm:end,2:dm:end-dm), zeros(size(Th(2:dm:end,2:dm:end-dm))), ...
        1.0, 'color',[0.15 0.15 0.15], 'Linewidth',2.5, 'MaxHeadSize',.008);
set(gca,'fontsize',42)
title('ND Cell Density and Velocity','fontsize',55,'interpreter','latex')
xlabel('$\theta$','fontsize',50,'interpreter','latex')
ylabel('$\hat{z}$','fontsize',50,'interpreter','latex')
zlabel('$\hat q(\theta,\hat{z})$ and $\vec{v}(\theta,\hat{z})$','fontsize',50,'interpreter','latex')
grid on; grid minor; box on
colormap turbo
shading interp
colorbar
caxis([1 max(Qhat_ss(:))])
xlim([0 2*pi])
ylim([0 1])
view([0 0 1])
hold off

%%  Animation 
dtPlot = round(.05*Nt);
dm = 5;

for i = 1:dtPlot:Nt

    Q = reshape(Qhat_full(i,:), Nth, Nz);
    [Vth, Vz] = ND_velocityField_fromFlux(Q, diffmat_hat, geom_hat, par_hat);

    figure(2)
    surf(Th, Z, Q);
    hold on
    quiver3(Th(2:dm:end,2:dm:end-dm), Z(2:dm:end,2:dm:end-dm), Q(2:dm:end,2:dm:end-dm), ...
        Vth(2:dm:end,2:dm:end-dm), Vz(2:dm:end,2:dm:end-dm), zeros(size(Th(2:dm:end,2:dm:end-dm))), ...
        1.0, 'color',[0.15 0.15 0.15], 'Linewidth',2.5, 'MaxHeadSize',.008);

    set(gca,'fontsize',42)
    title('ND Cell Density and Velocity','fontsize',55,'interpreter','latex')
    xlabel('$\theta$','fontsize',50,'interpreter','latex')
    ylabel('$z$','fontsize',50,'interpreter','latex')
    zlabel('$\hat q(\theta,z)$ and $\vec{v}(\theta,z)$','fontsize',50,'interpreter','latex')
    grid on; grid minor; box on
    colormap turbo
    shading interp
    colorbar
    xlim([0 2*pi])
    ylim([0 L])
    view([0 0 1])

    if i == 1
        pause
    end
    hold off
end

%% -------------------------
%  Plot on crypt domain (physical xyz, same as your last block)
% -------------------------
X_plot = Rphys.*cos(Th);
Y_plot = Rphys.*sin(Th);
Z_plot = Z;

dtPlot2 = round(.1*Nt);

for k = 1:dtPlot2:Nt

    Q = reshape(Qhat_full(k,:),Nth,Nz);

    figure(3)
    surf(X_plot, Y_plot, Z_plot, Q)
    set(gca,'fontsize',42)
    title('ND Cell Density on Crypt','fontsize',55,'interpreter','latex')
    xlabel('$x$','fontsize',50,'interpreter','latex')
    ylabel('$y$','fontsize',50,'interpreter','latex')
    zlabel('$z$','fontsize',50,'interpreter','latex')
    grid on; grid minor; box on
    shading interp
    colormap turbo
    colorbar
    lighting gouraud
    xlim([-40 40])
    ylim([-40 40])
    zlim([0 80])

    if k == 1
        pause
    end
end

%% Helpers

% Bottom/top Neumann (Flux_z = 0) for boundary column
function qB = solveFluxNeumannAtBoundary(jB, a0, a1, a2, q1, q2, diffmat_hat, geom_hat)
R_hat   = geom_hat.R_hat; 
Rth_hat = geom_hat.Rth_hat; 
Rz_hat  = geom_hat.Rz_hat; 
Dth     = diffmat_hat.Dth;
Nth     = length(Dth); 
Cz  = spdiags((Rth_hat(:,jB).^2 + R_hat(:,jB).^2), 0, Nth, Nth);
Cth = spdiags((Rth_hat(:,jB).*Rz_hat(:,jB)),    0, Nth, Nth);
A = a0*Cz - Cth*Dth;
b = -Cz*(a1*q1 + a2*q2);
qB = A \ b;
end