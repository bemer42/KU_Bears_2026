function [Vth, Vz_hat] = ND_compute_Velocity_Field(Q_hat, diffmat_hat, geom_hat, par_hat, bcType)

% -------------------------
% Collect nondimensional parameters:
% -------------------------
epsilon  = par_hat.epsilon;
gamma_th = par_hat.gamma_th;
gamma_z  = par_hat.gamma_z;

% -------------------------
% Collect geometry: 
% -------------------------
r     = geom_hat.R_hat;
rth   = geom_hat.Rth_hat;
rz    = geom_hat.Rz_hat;
S_hat = geom_hat.S_hat;

% -------------------------
% Collect differentiation matrices: 
% -------------------------
Dz  = diffmat_hat.Dz;
Dth = diffmat_hat.Dth; 

% -------------------------
% Define partial derivatives:
% -------------------------
Qth = Dth*Q_hat;
Qz  = Q_hat*Dz';

% -------------------------
% Define grid sizes:
% -------------------------
[Nth, Nz] = size(Q_hat);

% -------------------------
% Flux in theta:
% -------------------------
Flux_th = (1 ./ (Q_hat.^2) ./ S_hat) .* ( (1/gamma_th) * ((1+epsilon^2*rz.^2).*Qth - epsilon^2*(rth.*rz).*Qz) );
Vth = -(1 ./ Q_hat) .* Flux_th;

% -------------------------------------------------------------------------
% FV z-divergence on NONUNIFORM grid:
% -------------------------------------------------------------------------

% Extract 1D z-grid from your mesh
zh = geom_hat.Z_hat(1,:).';
dzf = diff(zh);

% Face averages of geometry + solution
rf   = 0.5*(r(:,1:end-1)   + r(:,2:end));
rthf = 0.5*(rth(:,1:end-1) + rth(:,2:end));
rzf  = 0.5*(rz(:,1:end-1)  + rz(:,2:end));
Qf   = 0.5*(Q_hat(:,1:end-1)   + Q_hat(:,2:end));
Qf   = max(Qf, 1e-10);

sF = sqrt(rthf.^2 + rf.^2 .* (1 + epsilon^2*rzf.^2));

% z-gradient at faces
qz_face = (Q_hat(:,2:end) - Q_hat(:,1:end-1)) ./ (dzf.' );

% theta-derivative at faces
DthQf   = Dth * Qf;

% Build face flux Flux_z at faces using the SAME numerator as your nodal Flux_z
num_face = -(rthf.*rzf).*DthQf + (rthf.^2 + rf.^2).*qz_face;

Flux_z_face = (1 ./ (Qf.^2) ./ sF) .* (1/gamma_z) .* num_face;

% enforce Neumann as face flux = 0 (same as RHS)
if bcType == "NbDt" || bcType == "NbNt"
    Flux_z_face(:,1) = 0;
end
if bcType == "NbNt" || bcType == "DbNt"
    Flux_z_face(:,end) = 0;
end

% face velocity (consistent with your definition V = -J/Q)
Vz_face = -(1 ./ Qf) .* Flux_z_face;   

% -------------------------------------------------------------------------

% -------------------------
% Convert face velocity to nodal velocity (simple adjacent-face average):
% -------------------------
Vz_hat           = zeros(Nth, Nz);
Vz_hat(:,2:Nz-1) = 0.5*(Vz_face(:,1:end-1) + Vz_face(:,2:end));
Vz_hat(:,1)      = Vz_face(:,1);
Vz_hat(:,Nz)     = Vz_face(:,end);
end