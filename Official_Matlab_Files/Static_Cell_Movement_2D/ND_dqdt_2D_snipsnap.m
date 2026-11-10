function dqhatdt = ND_dqdt_2D_snipsnap(t, q_hat_int, diffmat_hat, geom_hat, par_hat, bcType)

% -------------------------
% Collect physical parameters:
% -------------------------
gamma_th = par_hat.gamma_th;
gamma_z  = par_hat.gamma_z;
epsilon  = par_hat.epsilon;
ell      = par_hat.ell;
Zp_hat   = @(th,t) par_hat.Zp_hat(th,t);

% -------------------------
% Collect geometry: 
% -------------------------
Th      = geom_hat.Th;
Z_hat   = geom_hat.Z_hat;
R_hat   = geom_hat.R_hat;
Rth_hat = geom_hat.Rth_hat;
Rz_hat  = geom_hat.Rz_hat;   
S_hat   = geom_hat.S_hat;

% -------------------------
% Collect differentiation matrices: 
% -------------------------
Dz  = diffmat_hat.Dz;
Dth = diffmat_hat.Dth;

% -------------------------
% Define grid sizes:
% -------------------------
Nth = size(Dth,1);
Nz  = size(Dz,1);
nz  = Nz - 2;

% -------------------------
% Reshape interior state into matrix (Nth x (Nz-2)):
% -------------------------
Q_hat_int = reshape(q_hat_int, Nth, nz);

% -------------------------
% Initialize full matrix and insert interior z-columns:
% -------------------------
Q_hat = zeros(Nth, Nz);
Q_hat(:,2:Nz-1) = Q_hat_int;

% -------------------------
% Boundary derivative stencil coefficients from Dz:
% -------------------------
b0 = Dz(1,1);       b1 = Dz(1,2);       b2 = Dz(1,3);
t2 = Dz(end,end-2); t1 = Dz(end,end-1); t0 = Dz(end,end);

% -------------------------
% Helper that solves Flux_z(:,boundary)=0 for the boundary column qB = Q(:,jB)
% using one-sided stencil: a0*qB + a1*q1 + a2*q2 for q_z at boundary
% -------------------------    
    function qB = solveFluxNeumannAtBoundary(jB, a0, a1, a2, q1, q2)
        % Numerator for z-flux:  (Rth_hat^2+R_hat^2)*Q_zhat - (Rth_hat*Rz_hat)*Q_theta  = 0
        Cz  = spdiags((Rth_hat(:,jB).^2 + R_hat(:,jB).^2), 0, Nth, Nth);
        Cth = spdiags((Rth_hat(:,jB).*Rz_hat(:,jB)),       0, Nth, Nth);

        A = a0*Cz - Cth*Dth;
        b = -Cz*(a1*q1 + a2*q2);
        qB = A \ b;
    end

% -------------------------    
% Apply boundary conditions
% -------------------------    
switch bcType
    case "NbDt"
        Q_hat(:,1)   = solveFluxNeumannAtBoundary(1,  b0, b1, b2, Q_hat(:,2),     Q_hat(:,3));
        Q_hat(:,end) = ell;
    case "NbNt"
        Q_hat(:,1)   = solveFluxNeumannAtBoundary(1,  b0, b1, b2, Q_hat(:,2),     Q_hat(:,3));
        Q_hat(:,end) = solveFluxNeumannAtBoundary(Nz, t0, t1, t2, Q_hat(:,end-1), Q_hat(:,end-2));
    case "DbNt"
        Q_hat(:,1)   = ell;
        Q_hat(:,end) = solveFluxNeumannAtBoundary(Nz, t0, t1, t2, Q_hat(:,end-1), Q_hat(:,end-2));
    case "DbDt"
        Q_hat(:,1)   = ell;
        Q_hat(:,end) = ell;
    otherwise
        error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
end
Q_hat = max(Q_hat, 1e-10);

% -------------------------
% Flux in theta:
% -------------------------
Qth = Dth*Q_hat;
Qz  = Q_hat*Dz';

Jth = (1 ./ (Q_hat.^2) ./ S_hat) .* ...
      ( (1/gamma_th) * ( (1 + epsilon^2*Rz_hat.^2).*Qth - epsilon^2*(Rth_hat.*Rz_hat).*Qz ) );

% -------------------------------------------------------------------------
% FV z-divergence on NONUNIFORM grid:
% -------------------------------------------------------------------------

% Extract 1D z-grid from your mesh
zh = Z_hat(1,:).';
dzf = diff(zh);

% Control-volume "cell widths" around nodes, using midpoints
zf  = [zh(1); 0.5*(zh(1:end-1)+zh(2:end)); zh(end)];
dzc = diff(zf);

% Face averages
rf   = 0.5*(R_hat(:,1:end-1)   + R_hat(:,2:end));
rthf = 0.5*(Rth_hat(:,1:end-1) + Rth_hat(:,2:end));
rzf  = 0.5*(Rz_hat(:,1:end-1)  + Rz_hat(:,2:end));
Qf   = 0.5*(Q_hat(:,1:end-1)   + Q_hat(:,2:end));
Qf = max(Qf, 1e-10);

sHat_f = sqrt(rthf.^2 + rf.^2 .* (1 + epsilon^2*rzf.^2));

% Gradients at faces
qz_face  = (Q_hat(:,2:end) - Q_hat(:,1:end-1)) ./ (dzf.' );
DthQf    = Dth * Qf;

% Face z-flux J_z at faces
num_face = -(rthf.*rzf).*DthQf + (rthf.^2 + rf.^2).*qz_face;

Jz_face = (1 ./ (Qf.^2) ./ sHat_f) .* (1/gamma_z) .* num_face;

% Enforce Neumann as FACE FLUX = 0 (FV-consistent)
if bcType == "NbDt" || bcType == "NbNt"
    Jz_face(:,1) = 0;
end
if bcType == "NbNt" || bcType == "DbNt"
    Jz_face(:,end) = 0;
end

% Divergence in zhat at nodes (only interior used)
divJz = zeros(Nth, Nz);
divJz(:,2:Nz-1) = (Jz_face(:,2:Nz-1) - Jz_face(:,1:Nz-2)) ./ (dzc(2:Nz-1).');
% -------------------------------------------------------------------------

% -------------------------
% Cell proliferation:
% -------------------------
P = Q_hat .* (Z_hat < Zp_hat(Th,t));

% -------------------------
% PDE:
% -------------------------
dQhatdt_full = (1 ./ S_hat) .* (Dth*Jth) + (1 ./ S_hat) .* divJz + P;

% -------------------------
% Return only interior z-columns:
% -------------------------
dQhatdt = dQhatdt_full(:,2:end-1);
dqhatdt = dQhatdt(:);

end