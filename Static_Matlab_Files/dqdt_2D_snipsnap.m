function dqdt = dqdt_2D_snipsnap(t, q_int, diffmat, geom, par, bcType)

% Collect all parameters:
alpha_th = par.alpha_th;
alpha_z  = par.alpha_z;
Tc       = par.Tc;
zp       = par.zp;

% Collect geometry: 
Z   = geom.Z;
R   = geom.R;
Rth = geom.Rth; 
Rz  = geom.Rz; 
G   = geom.G;

% Collect differentiation matrices: 
Dz  = diffmat.Dz;
Dth = diffmat.Dth; 

% Sizes
Nth = size(Dth,1);
Nz  = size(Dz,1);
nz  = Nz - 2;

% Reshape interior state into matrix (Nth x (Nz-2))
Q_int = reshape(q_int, Nth, nz);

% Initialize full matrix and insert interior z-columns
Q = zeros(Nth, Nz);
Q(:,2:Nz-1) = Q_int;

% Convenience: boundary derivative stencil coefficients from Dz
b0 = Dz(1,1);       b1 = Dz(1,2);       b2 = Dz(1,3);
t2 = Dz(end,end-2); t1 = Dz(end,end-1); t0 = Dz(end,end);

% Helper that solves Flux_z(:,boundary)=0 for the boundary column qB = Q(:,jB)
% using one-sided stencil: a0*qB + a1*q1 + a2*q2 for q_z at boundary
    function qB = solveFluxNeumannAtBoundary(jB, a0, a1, a2, q1, q2)

        % Coefficients at that boundary column
        Cz  = spdiags(alpha_z*(Rth(:,jB).^2 + R(:,jB).^2), 0, Nth, Nth);
        Cth = spdiags(alpha_z*(Rth(:,jB).*Rz(:,jB)),      0, Nth, Nth);

        % Enforce (numerator of Flux_z)=0:
        %   Cz*(a0*qB + a1*q1 + a2*q2) - Cth*(Dth*qB) = 0
        A = a0*Cz - Cth*Dth;
        b = -Cz*(a1*q1 + a2*q2);
        qB = A \ b;
    end

% Apply boundary conditions
switch bcType
    case "NbDt"
        % bottom Neumann
        Q(:,1) = solveFluxNeumannAtBoundary(1, b0, b1, b2, Q(:,2), Q(:,3));
        % top Dirichlet
        Q(:,end) = 1;
    case "NbNt"
        % bottom Neumann
        Q(:,1) = solveFluxNeumannAtBoundary(1, b0, b1, b2, Q(:,2), Q(:,3));
        % top Neumann
        Q(:,end) = solveFluxNeumannAtBoundary(Nz, t0, t1, t2, Q(:,end-1), Q(:,end-2));
    case "DbNt"
        % bottom Dirichlet
        Q(:,1) = 1;
        % top Neumann
        Q(:,end) = solveFluxNeumannAtBoundary(Nz, t0, t1, t2, Q(:,end-1), Q(:,end-2));
    case "DbDt"
        % bottom Dirichlet
        Q(:,1) = 1;
        % top Dirichlet
        Q(:,end) = 1;
    otherwise
        error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
end

% --- Flux_th
Flux_th = (1 ./ (Q.^2) ./ G) .* ...
          (alpha_th*(Rz.^2+1).*(Dth*Q) - alpha_th*(Rth.*Rz).*(Q*Dz') );

% FV z-divergence on NONUNIFORM grid

% Extract 1D z-grid from your mesh
zv = Z(1,:).';
dzf = diff(zv);

% Control-volume "cell widths" around nodes, using midpoints
zf  = [zv(1); 0.5*(zv(1:end-1)+zv(2:end)); zv(end)];
dzc = diff(zf);

% Face averages of geometry + solution
rf   = 0.5*(R(:,1:end-1)   + R(:,2:end));
rthf = 0.5*(Rth(:,1:end-1) + Rth(:,2:end));
rzf  = 0.5*(Rz(:,1:end-1)  + Rz(:,2:end));
Qf   = 0.5*(Q(:,1:end-1)   + Q(:,2:end));

det_gf = rthf.^2 + rf.^2 .* (rzf.^2 + 1);

% z-gradient at faces
qz_face = (Q(:,2:end) - Q(:,1:end-1)) ./ (dzf.' );

% theta-derivative at faces
DthQf = Dth * Qf;

% Build face flux Flux_z at faces using the SAME numerator as your nodal Flux_z
num_face = -alpha_z*(rthf.*rzf).*DthQf + alpha_z*(rthf.^2 + rf.^2).*qz_face;

Flux_z_face = (1 ./ (Qf.^2) ./ sqrt(det_gf)) .* num_face;

% Enforce Neumann as FACE FLUX = 0 (this is the FV-consistent way)
if bcType == "NbDt" || bcType == "NbNt"
    Flux_z_face(:,1) = 0;          % bottom face
end
if bcType == "NbNt" || bcType == "DbNt"
    Flux_z_face(:,end) = 0;        % top face
end

% Divergence at nodes (Nth x Nz), only interior needed
divFz = zeros(Nth, Nz);
divFz(:,2:Nz-1) = (Flux_z_face(:,2:Nz-1) - Flux_z_face(:,1:Nz-2)) ./ (dzc(2:Nz-1).');

% --- Cell proliferation
P = log(2)/Tc * Q .* (Z < zp);

% --- PDE
dQdt_full = (1 ./ G) .* (Dth*Flux_th) + ...
    (1 ./ G) .* divFz + P;

% Return only interior z-columns
dQdt = dQdt_full(:,2:end-1);
dqdt = dQdt(:);

end