%% Diffusion 2D Nonlinear
close all; clear; clc

% Boundary condition set: 
bcType = "NbDt";

% Parameters
par = struct();
par.alpha_th = 1e-1;
par.alpha_z  = 1e0;
par.k        = 2e-4;
par.zu       = 50;
par.zb       = 5; 
par.thc      = pi/2;
par.thw      = pi/2;

% Discretize time:
N_t = 2e3;
t_0 = 0;
t_end = 2e3;
t = linspace(t_0,t_end,N_t);

% Discretize theta space:
Nth    = 4e1;
th_0   = 0;
th_end = 2*pi;
th     = linspace(th_0,th_end,Nth)';
dth    = th(2) - th(1);

% Discretize z space:
Nz    = 4e1;
z_0   = 0;
z_end = 78;
z     = linspace(z_0,z_end,Nz)';
dz    = z(2) - z(1);

%Create mesh grid:
[Z,Th] = meshgrid(z,th);

% Differentiation Matrices
Dth = spdiags([-1/2*ones(Nth,1) 1/2*ones(Nth,1)],[-1 1],Nth,Nth);
Dth(1,end) = -1/2;
Dth(end,1) =  1/2;
Dth = Dth/dth;

% Define the differentiation matrix:
Dz = spdiags([-1/2*ones(Nz,1) 1/2*ones(Nz,1)],[-1 1],Nz,Nz);
Dz(1,1:3) = [-3/2 2 -1/2];
Dz(end, end-2:end) = [1/2 -2 3/2];
Dz = Dz/dz;

% Radius function:
r_b   = 41/2/pi;
r_t   = 10/pi;
a     = 0.3;

r   = @(theta,z) r_b*(1 - exp(-a * z)) + r_t * exp(a * (z - z_end));
rth = @(theta,z) 0;
rz  = @(theta,z) -a*r_b*exp(-a*z) + a*r_t*exp(a*z-z_end);

% Initial Condition
gau = @(th,z,thc,zc,a,s) a*exp(-(th-thc).^2/s-(z-zc).^2/s/5);
f = @(th,z) 1 + 3*exp(-(th-pi).^2*5 - (z - 20).^2/15);
f = @(th,z) 1 + ones(size(th)).*(abs(z-15)<10).*(abs(th-pi)<.5) + exp(-(th-pi).^2*5 - (z - 20).^2/15);
% f = @(th,z) 1 + gau(th,z,pi/2,50,2,10) + gau(th,z,pi,20,2,1) + gau(th,z,3*pi/2,70,3,10);
Q0 = f(Th,Z);
Q0_int = Q0(:,2:end-1);
q0_int = Q0_int(:);

% Define the right hand side
dQdt = @(t, q_int) dqdt_2D_snipsnap(t, q_int, Th, Z, r(Th,Z), rth(Th,Z), rz(Th,Z), Dth, Dz, par, bcType);

% Solve the PDE:
tic
options = odeset('Stats','on');
[t, q_int] = ode23s(dQdt, t, q0_int, options);
toc

% Expand interior into full for each time point:
Q_full = zeros(N_t, Nth*Nz);
rth_0  = rth(Th,Z); rth_0 = rth_0(:,1);
rz_0   = rz(Th,Z);  rz_0  = rz_0(:,1);
r_0    = r(Th,Z);   r_0   = r_0(:,1);



for i = 1:N_t

    
    Q_int = reshape(q_int(i,:), Nth, Nz-2);

    % Initialize full snapshot:
    Q = zeros(Nth, Nz);

    % Put interior z-columns back in:
    Q(:,2:Nz-1) = Q_int;

    % Dirichlet at z = z_end:
    Q(:,Nz) = 1;

    % No-flux at z = 0 (Neumann q_z = 0), 2nd-order one-sided:
    % Z Interior Boundary Conditions:
    % -- no flux at z = z_0 (left side of matrix):
    % Need to solve for q_0 = Q(:,1);
    Cz   = spdiags(par.alpha_z*(rth_0.^2+r_0.^2),0,Nth,Nth);
    Cth  = spdiags(par.alpha_th*(rth_0.*rz_0),0,Nth,Nth);

    % Build A and b from:
    %  Cz*((-3/2*q0 + 2*q1 - 1/2*q2)/dz) - Cth*(Dth*q0) = 0
    A = Dz(1,1)*Cz - Cth*Dth;
    b = -Cz*(Dz(1,2)*Q_int(:,1) + Dz(1,3)*Q_int(:,2));
    % Solve for boundary values Q(:,1)
    Q(:,1) = A \ b;
    
    % Save flattened snapshot
    Q_full(i,:) = Q(:);
end

%% Animation 
dm = 4;                        
Thm = Th(1:dm:end, 1:dm:end);
Zm  = Z( 1:dm:end, 1:dm:end);
dt  = round(.01*N_t);

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
    zlim([0 5])
    clim([min(Q_full(:)) max(Q_full(:))])
    view(45,30)
    if i == 1
        pause 
    else
        drawnow
    end
    hold off
end

%% Plot on crypt domain
dt = round(.01*N_t);

X_plot = r(Th,Z).*cos(Th);
Y_plot = r(Th,Z).*sin(Th);
Z_plot = Z;

for k = 1:dt:N_t

    Q = reshape(Q_full(k,:),Nth,Nz);

    figure(2)
    surf(X_plot,Y_plot,Z_plot,Q)
%     shading interp
    colormap turbo
    colorbar
    xlim([-30 30])
    ylim([-30 30])
    zlim([0 80])
    clim([1 1.25])
    if k == 1
        pause
    else
        drawnow
    end

end


% Surface integral
R_grid   = r(Th, Z);
Rth_grid = rth(Th, Z);
Rz_grid  = rz(Th, Z);


dA_factor = sqrt(Rth_grid.^2 + R_grid.^2 .* (Rz_grid.^2 + 1));

total_population = zeros(N_t, 1);

for i = 1:N_t
    Q_frame = reshape(Q_full(i,:), Nth, Nz);

    
    integrand = Q_frame .* dA_factor;

    
    total_population(i) = trapz(z, trapz(th, integrand, 1));
end
figure(3)
plot(t, total_population, 'LineWidth', 3)
grid on; grid minor;
