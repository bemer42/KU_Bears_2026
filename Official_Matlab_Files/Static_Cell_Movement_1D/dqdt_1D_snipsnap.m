function dqdt = dqdt_1D_snipsnap(t, q_int, diffmat, geom, par, bcType)

% Note that input column vector q is only interior values
% -------------------------
% Collect parameters:
% -------------------------
alpha = par.alpha;
Tc    = par.Tc;
zp    = par.zp;

% -------------------------
% Collect geometry:
% -------------------------
z  = geom.z;
rz = geom.rz;
g  = geom.g;

% -------------------------
% Collect differentiation matrix:
% -------------------------
Dz = diffmat.Dz;

% -------------------------
% Boundary derivative weights from Dz:
% -------------------------
wL = full(Dz(1,1:3));                  
wR = full(Dz(end,end-2:end));          

switch bcType
    case "NbDt"
        q_l = -(wL(2)*q_int(1) + wL(3)*q_int(2))/wL(1);
        q_r = 1;

    case "NbNt"
        q_l = -(wL(2)*q_int(1) + wL(3)*q_int(2))/wL(1);
        q_r = -(wR(1)*q_int(end-1) + wR(2)*q_int(end))/wR(3);

    case "DbNt"
        q_l = 1;
        q_r = -(wR(1)*q_int(end-1) + wR(2)*q_int(end))/wR(3);

    case "DbDt"
        q_l = 1;
        q_r = 1;

    otherwise
        error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
end

% -------------------------
% Extend to full q:
% -------------------------
q_full = [q_l; q_int; q_r];

% -------------------------
% Apply PDE to q_full:
% -------------------------
dqdt  = (1./g).*(Dz*((1./g).*alpha./q_full.^2 .*(Dz*q_full)))  + ...
    log(2)/Tc .* q_full .*(z<zp);

% -------------------------
% Chop to interior q:
% -------------------------
dqdt = dqdt(2:end-1);

end