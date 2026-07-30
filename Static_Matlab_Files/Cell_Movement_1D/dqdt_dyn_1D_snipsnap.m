function dydt = dqdt_dyn_1D_snipsnap(t, y, diffmat, geom, par, bcType)

% Collect parameters:
alpha_z = par.alpha_z;
Tc      = par.Tc;
zp      = par.zp;
eta     = par.eta;
q_th    = par.q_th;
epsQ    = par.epsQ;
alpha_s = par.alpha_s;    
mu      = par.mu;      
beta    = par.beta;    
s_th    = par.s_th; 
epsS    = par.epsS;  

% Collect geometry grid:
z  = geom.z;
Nz = numel(z);

% Collect differentiation matrices:
Dz  = diffmat.Dz;
DzF = diffmat.DzF;
DzB = diffmat.DzB;

% --- Unpack state ---
Nq = Nz - 2;
Nr = Nz;
Ns = Nz - 2;

q_int  = y(1:Nq);
r_full = y(Nq+1:Nq+Nr);
s_int  = y(Nq+Nr+1:end);

% Apply boundary conditions for q (UNCHANGED)
switch bcType
    case "NbDt"
        q_l = -(Dz(1,2)*q_int(1) + Dz(1,3)*q_int(2))/Dz(1,1);
        q_r = 1;
    case "NbNt"
        q_l = -(Dz(1,2)*q_int(1) + Dz(1,3)*q_int(2))/Dz(1,1);
        q_r = -(Dz(end,end-1)*q_int(end) + Dz(end,end-2)*q_int(end-1))/Dz(end,end);
    case "DbNt"
        q_l = 1;
        q_r = -(Dz(end,end-1)*q_int(end) + Dz(end,end-2)*q_int(end-1))/Dz(end,end);
    case "DbDt"
        q_l = 1;
        q_r = 1;
    otherwise
        error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
end

% Extend to full q:
q_full = [q_l; q_int; q_r];

% Neumann Boundary conditions for s
s_l = -(Dz(1,2)*s_int(1) + Dz(1,3)*s_int(2))/Dz(1,1);
s_r = -(Dz(end,end-1)*s_int(end) + Dz(end,end-2)*s_int(end-1))/Dz(end,end);
s_full = [s_l; s_int; s_r];

% Geometry from evolving r
rz = Dz * r_full;
g  = sqrt(1 + rz.^2);

% Compute cell velocity and use it as wave speed
v_full = -(1./g).*alpha_z./(q_full.^3).*(Dz*q_full);

% Signal source turns on when q>1
Hq = 0.5*(1 + tanh((q_full - q_th)/epsQ));

% Signal PDE
vP = max(v_full, 0);   
vM = min(v_full, 0);   

dsdt_full = -(vP.*(DzB*s_full) + vM.*(DzF*s_full)) +...
            alpha_s*(1./g).*(Dz*((1./g).*(Dz*s_full))) - mu*s_full +...
            beta * Hq;

% Radius swelling drivin by signal s:
Hs = 0.5*(1 + tanh((s_full - s_th)/epsS));
r_tar = geom.r0 .* (1 + Hs);
drdt  = eta * (r_tar - r_full);

% Need rzt = (r_t)_z for dilution term
rzt = Dz * drdt;

% Full PDE for q (same as yours, but uses q_safe in division)
dqdt_full  = (1./g).*(Dz*((1./g).*alpha_z./q_full.^2 .* (Dz*q_full)))  + ...
             log(2)/Tc .* q_full .*(z<zp) - (rz.*rzt)./(g.^2).*q_full;

% Return stacked RHS
dydt = [dqdt_full(2:end-1);   
        drdt;              
        dsdt_full(2:end-1)]; 
end