function [Vth_node, Vz_node, Flux_th, Flux_z_face, Vz_face] = velocityField_FVconsistent_2D(Q, diffmat, geom, par, bcType)
% velocityField_FVconsistent_2D
% Build velocity fields (Vth,Vz) from fluxes in a way that is consistent with
% dqdt_2D_snipsnap:
%   - theta: uses the same nodal Flux_th expression
%   - z    : uses the same FV face flux Flux_z_face and converts to nodal Vz
%
% INPUTS
%   Q      : Nth x Nz full snapshot (with BC columns already filled)
%   diffmat.Dth, diffmat.Dz
%   geom.Z, geom.R, geom.Rth, geom.Rz, geom.G
%   par.alpha_th, par.alpha_z
%   bcType : "NbDt","NbNt","DbNt","DbDt"
%
% OUTPUTS
%   Vth_node   : Nth x Nz nodal v^theta
%   Vz_node    : Nth x Nz nodal v^z (from FV face flux averaged to nodes)
%   Flux_th    : Nth x Nz nodal flux in theta used by PDE
%   Flux_z_face: Nth x (Nz-1) face flux in z used by PDE
%   Vz_face    : Nth x (Nz-1) face velocity in z

alpha_th = par.alpha_th;
alpha_z  = par.alpha_z;

Z   = geom.Z;
R   = geom.R;
Rth = geom.Rth;
Rz  = geom.Rz;
G   = geom.G;

Dth = diffmat.Dth;

[Nth, Nz] = size(Q);

% ---------------- theta flux (same as dqdt_2D_snipsnap) ----------------
Flux_th = (1 ./ (Q.^2) ./ G) .* ...
          (alpha_th*(Rz.^2+1).*(Dth*Q) - alpha_th*(Rth.*Rz).*(Q*diffmat.Dz') );

% Velocity relation used in your driver:
%   Vth = -(1/Q) * Flux_th
Vth_node = -(1 ./ Q) .* Flux_th;

% ---------------- z flux via FV faces (same as dqdt_2D_snipsnap) ----------------
% Extract 1D z-grid
zv = Z(1,:).';
dzf = diff(zv);                 % 1 x (Nz-1) face spacings

% Face averages
rf   = 0.5*(R(:,1:end-1)   + R(:,2:end));
rthf = 0.5*(Rth(:,1:end-1) + Rth(:,2:end));
rzf  = 0.5*(Rz(:,1:end-1)  + Rz(:,2:end));
Qf   = 0.5*(Q(:,1:end-1)   + Q(:,2:end));

det_gf = rthf.^2 + rf.^2 .* (rzf.^2 + 1);

% z-gradient at faces
qz_face = (Q(:,2:end) - Q(:,1:end-1)) ./ (dzf.' );   % (Nth x (Nz-1))

% theta-derivative at faces
DthQf = Dth * Qf;

% Face flux numerator (exactly as dqdt_2D_snipsnap)
num_face = -alpha_z*(rthf.*rzf).*DthQf + alpha_z*(rthf.^2 + rf.^2).*qz_face;

Flux_z_face = (1 ./ (Qf.^2) ./ sqrt(det_gf)) .* num_face;

% Enforce Neumann as FACE FLUX = 0 (FV-consistent BC handling)
if bcType == "NbDt" || bcType == "NbNt"
    Flux_z_face(:,1) = 0;      % bottom face
end
if bcType == "NbNt" || bcType == "DbNt"
    Flux_z_face(:,end) = 0;    % top face
end

% Face velocity consistent with your driver relation V = -(1/Q)*Flux
Vz_face = -(1 ./ Qf) .* Flux_z_face;

% Convert face velocity to nodal velocity (simple adjacent-face average)
Vz_node = zeros(Nth, Nz);
Vz_node(:,2:Nz-1) = 0.5*(Vz_face(:,1:Nz-2) + Vz_face(:,2:Nz-1));
Vz_node(:,1)      = Vz_face(:,1);
Vz_node(:,Nz)     = Vz_face(:,end);

end