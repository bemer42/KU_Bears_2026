function [C, P, Renewal_Time] = ND_2D_Function(z_0, gamma_th, gamma_z, zp_hat, bcType)
%ND_2D_Function  Solve ND 2D model and return C, P, and renewal time.
%
% Inputs
%   gamma_th   : nondimensional group (theta) 
%   gamma_z    : nondimensional group (z)
%   Zp_hat_in  : either
%                 - scalar in (0,1) meaning constant proliferative height, or
%                 - function handle Zp_hat(th,t_hat) returning size(th)
%   bcType     : "NbDt", "NbNt", "DbNt", or "DbDt"
%
% Outputs
%   C            : total cells (from compute_CP)
%   P            : proliferative cells (from compute_CP)
%   Renewal_Time : renewal time (output of compute_ND_Renewal_Time)
%
% Requirements (must be on MATLAB path):
%   diffmat_nonuniform, ND_dqdt_2D_snipsnap, solveFluxNeumannAtBoundary,
%   compute_CP, ND_velocityField_FVconsistent_2D, compute_ND_Renewal_Time

% -------------------------
% User-set / run controls
% -------------------------
uniform = 0;

% -------------------------
% Grid/Time points:
% -------------------------
Nt     = 5e3;
t_end  = 1e5;
Nz     = 6e1;
Nth    = 6e1;

% -------------------------
% Predefine structures and boundary type:
% -------------------------
geom_hat    = struct();
diffmat_hat = struct();
par_hat     = struct();

% -------------------------
%  Dimensional (fixed) parameters:
% -------------------------
Tc       = 15.1954;
L        = 79.3515;
ell      = 1;

% -------------------------
% Parameter functions:
% -------------------------
Zp_hat   = @(th,t_hat) zp_hat*ones(size(th));
% Zp_hat   = @(th,t_hat) (27 - 21*exp(-(th-pi/2).^2*3) - 21.*exp(-(th-3*pi/2).^2*3))/L;

% -------------------------
% Geometry parameters:
% -------------------------
r_b   = 18/pi;
r_t   = 10/pi;
a     = 0.3;

% -------------------------
% Boundary condition set:
% -------------------------
bcType = string(bcType);
mustBeMember(bcType, ["NbDt","NbNt","DbNt","DbDt"]);

% -------------------------
% Nondimensional groups
% -------------------------
epsilon = r_b/L;

% -------------------------
% Time discretization (ND time)
% -------------------------
t_hat_0   = 0;
t_hat_end = t_end*log(2)./Tc;
t_hat     = linspace(t_hat_0, t_hat_end, Nt);

% -------------------------
% Theta grid
% -------------------------
th_0    = 0;
th_end  = 2*pi;
th      = linspace(th_0, th_end, Nth+1).';
th(end) = [];
dth     = th(2) - th(1);

% -------------------------
% z-hat grid in [0,1]
% -------------------------
z0 = 0;
if uniform == 1
    z_hat  = linspace(z0,1,Nz).';
else
    alpha = .75;
    beta  = 4;
    xi    = linspace(0,1,Nz).';
    z_hat = (tanh(beta*(xi-alpha))-tanh(-beta*alpha))./(tanh(beta*(1-alpha))-tanh(-beta*alpha));
end

% -------------------------
% Mesh grids (Nth x Nz)
% -------------------------
[Z_hat, Th] = meshgrid(z_hat, th);

% -------------------------
% Differentiation in theta
% -------------------------
Dth = spdiags([-0.5*ones(Nth,1) 0.5*ones(Nth,1)], [-1 1], Nth, Nth);
Dth(1,end) = -0.5;
Dth(end,1) =  0.5;
Dth = Dth / dth;

% -------------------------
% Differentiation in z (nonuniform)
% -------------------------
Dz = diffmat_nonuniform(z_hat);

% -------------------------
% ND geometry on (theta, zhat)
% -------------------------
R_hat   = (1 - exp(-a*L*Z_hat)) + (r_t/r_b)*exp(a*L*(Z_hat - 1));
Rz_hat  = R_hat*Dz';
Rth_hat = Dth*R_hat;
S_hat   = sqrt(Rth_hat.^2 + R_hat.^2 .* (1 + epsilon^2*(Rz_hat.^2)));

% -------------------------
% Structures for RHS
% -------------------------
diffmat_hat.Dth   = Dth;
diffmat_hat.Dz    = Dz;

geom_hat.th       = th;
geom_hat.z_hat    = z_hat;
geom_hat.Th       = Th;
geom_hat.Z_hat    = Z_hat;
geom_hat.R_hat    = R_hat;
geom_hat.Rth_hat  = Rth_hat;
geom_hat.Rz_hat   = Rz_hat;
geom_hat.S_hat    = S_hat;

par_hat.gamma_th  = gamma_th;
par_hat.gamma_z   = gamma_z;
par_hat.epsilon   = epsilon;
par_hat.ell       = ell;
par_hat.Zp_hat    = @(th, t_hat) Zp_hat(th, t_hat);

% -------------------------
% Initial condition
% -------------------------
Q_hat_0     = ones(size(Th));        
Q_hat_0_int = Q_hat_0(:,2:end-1);
q_hat_0_int = Q_hat_0_int(:);

% -------------------------
% ND RHS handle and solve
% -------------------------
dQhatdt = @(t_hat, qint) ND_dqdt_2D_snipsnap(t_hat, qint, diffmat_hat, geom_hat, par_hat, bcType);

options = odeset('Stats','off','MaxStep',inf);
[t_hat, q_hat_int] = ode15s(dQhatdt, t_hat, q_hat_0_int, options);

% -------------------------
%  Expand interior back to full Q:
% -------------------------
Qhat_full = zeros(Nt, Nth*Nz);

% -------------------------
% One-sided zhat-derivative stencil coefficients from Dz:
% -------------------------
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
            Q(:,1)   = ell;  
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
%  Steady state analysis and crypt metrics:
% -------------------------
Qhat_ss = reshape(Qhat_full(end,:), Nth, Nz);

[Vth_ss, Vz_hat_ss] = ND_compute_Velocity_Field(Qhat_ss, diffmat_hat, geom_hat, par_hat, bcType);

[C, P] = compute_Cell_Totals(Qhat_ss, S_hat, th, z_hat, Z_hat < Zp_hat(th,t_hat(end)), r_b*L/ell);

th_0_vec = linspace(th(2), th(end-1), 10);
Renewal_Time = ND_compute_Renewal_Time(z_0, th_0_vec, Vth_ss, Vz_hat_ss, Tc, geom_hat);

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