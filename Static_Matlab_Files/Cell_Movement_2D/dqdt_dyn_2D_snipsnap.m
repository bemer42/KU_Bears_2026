function dYdt = dqdt_dyn_2D_snipsnap(t, Y, diffmat, geom, par, bcType)

% Collect all parameters
alpha_th   = par.alpha_th;
alpha_z    = par.alpha_z;
Tc         = par.Tc;
zp         = par.zp;
alpha_s    = par.alpha_s;
beta       = par.beta;
mu         = par.mu;
q_th       = par.q_th;
epsQ       = par.epsQ;
eta        = par.eta;
s_th       = par.s_th;
epsS       = par.epsS;

% Collect geometry + grids
Z    = geom.Z;
Th   = geom.Th;
dzf  = geom.dzf;
dzc  = geom.dzc;
dth  = geom.dth;
R0   = geom.R0;
Rtar = geom.Rtar;

% Collect differentiation matrices
Dz   = diffmat.Dz;
Dth  = diffmat.Dth;

% Sizes
Nth = size(Dth,1);
Nz  = size(Dz,1);
nz  = Nz - 2;

% Unpack state
nq = Nth*(Nz-2);
ns = Nth*(Nz-2);
nr = Nth*Nz;

q_int = Y(1:nq);
s_int = Y(nq+(1:ns));
r_vec = Y(nq+ns+(1:nr));

Q_int = reshape(q_int, Nth, nz);
S_int = reshape(s_int, Nth, nz);
R     = reshape(r_vec, Nth, Nz);

% Geometry from current radius
Rth = Dth*R;
Rz  = R*Dz';
G   = sqrt(Rth.^2 + R.^2.*(Rz.^2 + 1));   % sqrt(|g|)

% Build full Q and S and apply BCs
Q = zeros(Nth, Nz);
Q(:,2:Nz-1) = Q_int;

S = zeros(Nth, Nz);
S(:,2:Nz-1) = S_int;

% One-sided z-derivative stencil coefficients from Dz
b0 = Dz(1,1);       b1 = Dz(1,2);       b2 = Dz(1,3);
t2 = Dz(end,end-2); t1 = Dz(end,end-1); t0 = Dz(end,end);

% Nested helper that solves Flux_z(:,boundary)=0 for Q(:,jB)
    function qB = solveFluxNeumannAtBoundary(jB, a0, a1, a2, q1, q2)
        Cz  = spdiags(alpha_z*(Rth(:,jB).^2 + R(:,jB).^2), 0, Nth, Nth);
        Cth = spdiags(alpha_z*(Rth(:,jB).*Rz(:,jB)),      0, Nth, Nth);
        A = a0*Cz - Cth*Dth;
        b = -Cz*(a1*q1 + a2*q2);
        qB = A \ b;
    end

% Apply q BCs
switch bcType
    case "NbDt"
        Q(:,1)   = solveFluxNeumannAtBoundary(1,  b0, b1, b2, Q(:,2), Q(:,3));
        Q(:,end) = 1;
    case "NbNt"
        Q(:,1)   = solveFluxNeumannAtBoundary(1,  b0, b1, b2, Q(:,2), Q(:,3));
        Q(:,end) = solveFluxNeumannAtBoundary(Nz, t0, t1, t2, Q(:,end-1), Q(:,end-2));
    case "DbNt"
        Q(:,1)   = 1;
        Q(:,end) = solveFluxNeumannAtBoundary(Nz, t0, t1, t2, Q(:,end-1), Q(:,end-2));
    case "DbDt"
        Q(:,1)   = 1;
        Q(:,end) = 1;
    otherwise
        error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
end

% Apply s BCs: Neumann s_z=0 at both ends
S(:,1)   = -(b1*S(:,2) + b2*S(:,3))/b0;
S(:,end) = -(t1*S(:,end-1) + t2*S(:,end-2))/t0;

% r_t: signal-gated relaxation toward target radius
Hs   = 0.5*(1 + tanh((S - s_th)/epsS));
Reff = R0 + Hs.*(Rtar - R0);
Rt   = eta*(Reff - R);

% Dilution term for q: -(G_t/G) q
Rt_th = Dth*Rt;
Rt_z  = Rt*Dz';

detg = G.^2;
Gt_over_G = (Rth.*Rt_th + R.*Rt.*(1+Rz.^2) + (R.^2).*Rz.*Rt_z) ./ detg;

% q theta-flux 
Flux_th = (1 ./ (Q.^2) ./ G) .* ...
          (alpha_th*(Rz.^2+1).*(Dth*Q) - alpha_th*(Rth.*Rz).*(Q*Dz') );

% q z-divergence 
% Face averages
rf   = 0.5*(R(:,1:end-1)   + R(:,2:end));
rthf = 0.5*(Rth(:,1:end-1) + Rth(:,2:end));
rzf  = 0.5*(Rz(:,1:end-1)  + Rz(:,2:end));
Qf   = 0.5*(Q(:,1:end-1)   + Q(:,2:end));

det_gf = rthf.^2 + rf.^2 .* (rzf.^2 + 1);

% z-gradient at faces
qz_face = (Q(:,2:end) - Q(:,1:end-1)) ./ (dzf.' );

% theta-derivative at faces
DthQf = Dth * Qf;

% Face z-flux
num_face = -alpha_z*(rthf.*rzf).*DthQf + alpha_z*(rthf.^2 + rf.^2).*qz_face;
Flux_z_face = (1 ./ (Qf.^2) ./ sqrt(det_gf)) .* num_face;

% Neumann as FACE FLUX = 0 (FV-consistent)
if bcType == "NbDt" || bcType == "NbNt"
    Flux_z_face(:,1) = 0;
end
if bcType == "NbNt" || bcType == "DbNt"
    Flux_z_face(:,end) = 0;
end

divFz = zeros(Nth, Nz);
divFz(:,2:Nz-1) = (Flux_z_face(:,2:Nz-1) - Flux_z_face(:,1:Nz-2)) ./ (dzc(2:Nz-1).');

% Proliferation term
th0 = 3*pi/2;  w = 0.5;  epsTh = 0.05;
Wth = 0.5*(1 - tanh((abs(Th-th0)-w)/epsTh));  
% P   = log(2)/Tc * Q .* (Z<(zp+8*cos(2*Th)-8));
P   = log(2)/Tc * Q .* (Z<zp);

% Dynamic q PDE 
dQdt_full = (1 ./ G) .* (Dth*Flux_th) + (1 ./ G) .* divFz + P - Gt_over_G .* Q;

% Velocity fields (use your FV-consistent function)
geom_local     = geom;
geom_local.R   = R;
geom_local.Rth = Rth;
geom_local.Rz  = Rz;
geom_local.G   = G;

[Vth, Vz] = velocityField_FVconsistent_2D(Q, diffmat, geom_local, par, bcType);

% Signal source switch H(q-q_th)
Hq = 0.5*(1 + tanh((Q - q_th)/epsQ));

% Signal advection (upwind in theta and z)
% theta upwind (periodic)
Sth_fwd = (circshift(S,-1,1) - S)/dth;
Sth_bwd = (S - circshift(S, 1,1))/dth;
adv_th  = -( max(Vth,0).*Sth_bwd + min(Vth,0).*Sth_fwd );

% z upwind (nonuniform)
Sz_fwd = zeros(Nth,Nz);
Sz_bwd = zeros(Nth,Nz);

Sz_fwd(:,1:Nz-1) = (S(:,2:Nz) - S(:,1:Nz-1)) ./ (dzf.');
Sz_fwd(:,Nz)     = Sz_fwd(:,Nz-1);

Sz_bwd(:,2:Nz) = (S(:,2:Nz) - S(:,1:Nz-1)) ./ (dzf.');
Sz_bwd(:,1)    = Sz_bwd(:,2);

adv_z = -( max(Vz,0).*Sz_bwd + min(Vz,0).*Sz_fwd );

adv = adv_th + adv_z;

% Signal diffusion (theta + FV in z) + source/decay
% --- Laplace–Beltrami diffusion: (1/G)( dth(Jth) + dz(Jz) )
Sth = Dth*S;
Sz  = S*Dz';

% theta flux Jth at nodes
Jth_s = alpha_s * ( (Rz.^2+1).*Sth - (Rth.*Rz).*Sz ) ./ G;

% z flux Jz on faces (FV-consistent, analogous to q)
Sf   = 0.5*(S(:,1:end-1) + S(:,2:end));
DthSf = Dth*Sf;
sz_face = (S(:,2:end) - S(:,1:end-1)) ./ (dzf.');

% face-averaged geometry (same as q block)
rf   = 0.5*(R(:,1:end-1)   + R(:,2:end));
rthf = 0.5*(Rth(:,1:end-1) + Rth(:,2:end));
rzf  = 0.5*(Rz(:,1:end-1)  + Rz(:,2:end));
det_gf = rthf.^2 + rf.^2 .* (rzf.^2 + 1);

num_face_s = -(rthf.*rzf).*DthSf + (rthf.^2 + rf.^2).*sz_face;
Jz_s_face  = alpha_s * (1 ./ sqrt(det_gf)) .* num_face_s;

% Neumann (no diffusive flux) in z:
Jz_s_face(:,1)   = 0;
Jz_s_face(:,end) = 0;

divJz = zeros(Nth,Nz);
divJz(:,2:Nz-1) = (Jz_s_face(:,2:Nz-1) - Jz_s_face(:,1:Nz-2)) ./ (dzc(2:Nz-1).');

LBs = (1 ./ G) .* (Dth*Jth_s) + (1 ./ G) .* divJz;

dSdt_full = adv + LBs + beta*Hq - mu*S;

% Return only interior z-columns for q and s; full for r
dQdt = dQdt_full(:,2:end-1);
dSdt = dSdt_full(:,2:end-1);

dYdt = [dQdt(:); dSdt(:); Rt(:)];

end