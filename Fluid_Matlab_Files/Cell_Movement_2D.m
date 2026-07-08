%% Cell Movement 2D Nonlinear
close all; clear; clc

% Boundary condition set:
bcType = "NbDt";
 
% Parameters
par = struct();
par.alpha_th = 1e-1;
par.alpha_z  = 1e0;
par.k        = 1e-4;
par.zu       = 27;
par.zb       = 0;
par.thc      = 0;
par.thw      = 3*pi;

% Discretize time:
N_t = 2e1;
t_0 = 0;
t_end = 1e5;
t = linspace(t_0,t_end,N_t);

% Discretize theta space:
Nth    = 5e1;
th_0   = 0;
th_end = 2*pi;
th     = linspace(th_0,th_end,Nth)';
dth    = th(2) - th(1);

% Discretize z space:
Nz   = 5e1;
z_0  = 0;
L    = 78;
% --- Uniform Grid ---
z    = linspace(z_0,L,Nz)';
dz   = z(2) - z(1);
% --- Non-Uniform Grid --- 
beta = 3;
xi   = linspace(0,1,Nz)';
z    = L/2*(1+tanh(beta*(2*xi-1))/tanh(beta));

% Create mesh grid:
[Z,Th] = meshgrid(z,th);

% Differentiation Matrices
Dth = spdiags([-1/2*ones(Nth,1) 1/2*ones(Nth,1)],[-1 1],Nth,Nth);
Dth(1,end) = -1/2;
Dth(end,1) =  1/2;
Dth = Dth/dth;

% Define the differentiation matrix:
% --- Uniform Grid --- 
Dz = spdiags([-1/2*ones(Nz,1) 1/2*ones(Nz,1)],[-1 1],Nz,Nz);
Dz(1,1:3) = [-3/2 2 -1/2];
Dz(end, end-2:end) = [1/2 -2 3/2];
Dz = Dz/dz;
% --- Non-Uniform Grid --- 
Dz = diffmat_nonuniform(z);

% Radius function:
r_b   = 41/2/pi;
r_t   = 10/pi;
a     = 0.3;

r   = @(theta,z) r_b*(1 - exp(-a * z)) + r_t * exp(a * (z - L));
rth = @(theta,z) 0*theta;
rz  = @(theta,z) -a*r_b*exp(-a*z) + a*r_t*exp(a*z-L);

% Initial Condition
gau = @(th,z,thc,zc,a,s) a*exp(-(th-thc).^2/s-(z-zc).^2/s/5);
% f = @(th,z) 1 + exp(-(th-pi).^2*5 - (z - 20).^2/15);
% f = @(th,z) 1 + ones(size(th)).*(abs(z-35)<10).*(abs(th-pi)<pi/4) + exp(-(th-pi).^2*5 - (z - 20).^2/15);
% f = @(th,z) 1 + gau(th,z,pi/2,50,2,10) + gau(th,z,pi,20,2,1) + gau(th,z,3*pi/2,70,3,10);
f = @(th,z) ones(size(th));
Q0 = f(Th,Z);
Q0_int = Q0(:,2:end-1);
q0_int = Q0_int(:);

% Define the right hand side
dQdt = @(t, q_int) dqdt_2D_snipsnap(t, q_int, Th, Z, r(Th,Z), rth(Th,Z), rz(Th,Z), Dth, Dz, par, bcType);

% Solve the PDE:
tic
% dtdiff = 0.95*(min(diff(z))^2/par.alpha_z);
dtdiff = inf;
options = odeset('Stats','on','MaxStep',dtdiff);
[t, q_int] = ode15s(dQdt, t, q0_int, options);
toc

% Expand interior into full for each time point:
Q_full = zeros(N_t,Nth*Nz);

% Precompute geometry arrays:
R   = r(Th,Z);
RTH = rth(Th,Z);
RZ  = rz(Th,Z);

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
            Cz  = spdiags(par.alpha_z*(RTH(:,1).^2 + R(:,1).^2), 0, Nth, Nth);
            Cth = spdiags(par.alpha_th*(RTH(:,1).*RZ(:,1)),      0, Nth, Nth);
            A = b0*Cz - Cth*Dth;
            b = -Cz*(b1*Q(:,2) + b2*Q(:,3));
            Q(:,1) = A \ b;
            % Top Dirichlet
            Q(:,end) = 1;
        case "NbNt"
            % Bottom Neumann
            Cz  = spdiags(par.alpha_z*(RTH(:,1).^2 + R(:,1).^2), 0, Nth, Nth);
            Cth = spdiags(par.alpha_th*(RTH(:,1).*RZ(:,1)),      0, Nth, Nth);
            A = b0*Cz - Cth*Dth;
            b = -Cz*(b1*Q(:,2) + b2*Q(:,3));
            Q(:,1) = A \ b;
            % Top Neumann
            Cz  = spdiags(par.alpha_z*(RTH(:,Nz).^2 + R(:,Nz).^2), 0, Nth, Nth);
            Cth = spdiags(par.alpha_th*(RTH(:,Nz).*RZ(:,Nz)),      0, Nth, Nth);
            A = t0*Cz - Cth*Dth;
            b = -Cz*(t1*Q(:,end-1) + t2*Q(:,end-2));
            Q(:,end) = A \ b;
        case "DbNt"
            % Bottom Dirichlet
            Q(:,1) = 1;
            % Top Neumann (Flux_z = 0): solve Q(:,end)
            Cz  = spdiags(par.alpha_z*(RTH(:,Nz).^2 + R(:,Nz).^2), 0, Nth, Nth);
            Cth = spdiags(par.alpha_th*(RTH(:,Nz).*RZ(:,Nz)),      0, Nth, Nth);
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

%% Animation
dm = 4;
Thm = Th(1:dm:end, 1:dm:end);
Zm  = Z( 1:dm:end, 1:dm:end);
dt  = round(.1*N_t);
 
for i = 1:dt:N_t

    Q = reshape(Q_full(i,:), Nth, Nz);
    Qm = Q(1:dm:end, 1:dm:end);

    figure(1)
    surf(Th, Z, Q);
    %     surf(Thm,Zm,Qm);
%     shading interp
    colormap summer
    lighting gouraud
    grid on
    set(gca,'fontsize',16)
    xlabel('\theta','fontsize',20)
    ylabel('z','fontsize',20)
    zlabel('q','fontsize',20)
    title('Cell Movement PDE on (\theta,z)','fontsize',25)
    xlim([0 2*pi])
    ylim([0 78])
    zlim([0.9 1.2])
    caxis([min(Q_full(:)) max(Q_full(:))])
    view(45,30)
    if i == 1
        pause
    else
        drawnow
    end
    hold off
end

%% Plot on crypt domain
dt = round(.1*N_t);

X_plot = r(Th,Z).*cos(Th);
Y_plot = r(Th,Z).*sin(Th);
Z_plot = Z;

for k = 1:dt:N_t

    Q = reshape(Q_full(k,:),Nth,Nz);

    figure(2)
    surf(X_plot,Y_plot,Z_plot,Q)
    %shading interp
    colormap turbo
    colorbar
    xlim([-30 30])
    ylim([-30 30])
    zlim([0 80])
    caxis([1 1.25])
    if k == 1
        pause
    else
        drawnow
    end

end

%% Steady State Analysis

X_plot = r(Th,Z).*cos(Th);
Y_plot = r(Th,Z).*sin(Th);
Z_plot = Z;

Q_ss = reshape(Q_full(end,:),Nth,Nz);

figure(3)
surf(X_plot,Y_plot,Z_plot,Q_ss)
shading interp
colormap turbo
colorbar
xlim([-30 30])
ylim([-30 30])
zlim([0 80])
caxis([1 1.05])

