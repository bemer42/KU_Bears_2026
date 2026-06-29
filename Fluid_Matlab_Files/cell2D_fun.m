
function total_cells = cell2D_fun(k, zu)



% Boundary condition set:
bcType = "DbDt";

% Parameters
par = struct();
par.alpha_th = 1e-1;
par.alpha_z  = 1e0;
par.k        = k;
par.zu       = zu;
par.zb       = 0;
par.thc      = 0;
par.thw      = 3*pi;

% Discretize time:
N_t = 2e3;
t_0 = 0;
t_end = 2e5;
t = linspace(t_0,t_end,N_t);

% Discretize theta space:
Nth    = 5e1;
th_0   = 0;
th_end = 2*pi;
th     = linspace(th_0,th_end,Nth)';
dth    = th(2) - th(1);

% Discretize z space:
Nz    = 5e1;
z_0   = 0;
z_end = 78;
z     = linspace(z_0,z_end,Nz)';
dz    = z(2) - z(1);

% Create mesh grid:
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
rth = @(theta,z) 0*theta;
rz  = @(theta,z) -a*r_b*exp(-a*z) + a*r_t*exp(a*z-z_end);

% Initial Condition
f = @(th,z) ones(size(th));
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


%% Surface integral
R_grid   = r(Th, Z);
Rth_grid = rth(Th, Z);
Rz_grid  = rz(Th, Z);

Q_ss = reshape(Q_full(i,:), Nth, Nz);


dA_factor = sqrt(Rth_grid.^2 + R_grid.^2 .* (Rz_grid.^2 + 1));


integrand = Q_ss .* dA_factor;


total_cells = trapz(z, trapz(th, integrand, 1));


end