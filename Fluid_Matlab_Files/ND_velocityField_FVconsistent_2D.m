function [Vth, Vz_hat] = ND_velocityField_FVconsistent_2D(Q, diffmat_hat, geom_hat, par_hat, bcType)

Dth = diffmat_hat.Dth;

% geometry
r     = geom_hat.R_hat;
rth   = geom_hat.Rth_hat;
rz    = geom_hat.Rz_hat;
S     = geom_hat.S_hat;

eps   = par_hat.epsilon;
gth   = par_hat.gamma_th;
gz    = par_hat.gamma_z;

% ---- Vtheta from nodal Jtheta (this part is fine)
Dz = diffmat_hat.Dz;
Qth = Dth*Q;
Qz  = Q*Dz';

Jth = (1 ./ (Q.^2) ./ S) .* ( (1/gth) * ((1+eps^2*rz.^2).*Qth - eps^2*(rth.*rz).*Qz) );
Vth = -Jth ./ Q;

% ---- Vz_hat from FV face flux (consistent with RHS)
zh = geom_hat.Z_hat(1,:).';
dzf = diff(zh);

rf   = 0.5*(r(:,1:end-1)   + r(:,2:end));
rthf = 0.5*(rth(:,1:end-1) + rth(:,2:end));
rzf  = 0.5*(rz(:,1:end-1)  + rz(:,2:end));
Qf   = 0.5*(Q(:,1:end-1)   + Q(:,2:end));
Qf   = max(Qf, 1e-10);

sF = sqrt(rthf.^2 + rf.^2 .* (1 + eps^2*rzf.^2));

qz_face = (Q(:,2:end) - Q(:,1:end-1)) ./ (dzf.' );
DthQf   = Dth * Qf;

num_face = -(rthf.*rzf).*DthQf + (rthf.^2 + rf.^2).*qz_face;
Jz_face = (1 ./ (Qf.^2) ./ sF) .* (1/gz) .* num_face;

% enforce Neumann as face flux = 0 (same as RHS)
if bcType == "NbDt" || bcType == "NbNt"
    Jz_face(:,1) = 0;
end
if bcType == "NbNt" || bcType == "DbNt"
    Jz_face(:,end) = 0;
end

% map face fluxes to nodes by averaging adjacent faces
[Nth, Nz] = size(Q);
Jz_node = zeros(Nth, Nz);
Jz_node(:,2:Nz-1) = 0.5*(Jz_face(:,1:end-1) + Jz_face(:,2:end));

Vz_hat = -Jz_node ./ Q;
end