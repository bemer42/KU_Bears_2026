
function dqdt = ND_dqdt_1D_snipsnap(t_hat, q_hat_int, diffmat_hat, geom_hat, par_hat, bcType)


% Collect parameters:
gamma  = par_hat.gamma;
ell    = par_hat.ell;
zp_hat = par_hat.zp_hat;

% Collect geometry parameters:
z_hat  = geom_hat.z_hat;
g_hat  = geom_hat.g_hat;

% Collect differentiation matrix:
Dz_hat = diffmat_hat.Dz_hat;

% Apply boundary conditions
switch bcType
    case "NbDt"
        % bottom Neumann
        q_l = -(Dz_hat(1,2)*q_hat_int(1) + Dz_hat(1,3)*q_hat_int(2))/Dz_hat(1,1);
        % top Dirichlet
        q_r = ell;
    case "NbNt"
        % bottom Neumann
        q_l = -(Dz_hat(1,2)*q_hat_int(1) + Dz_hat(1,3)*q_hat_int(2))/Dz_hat(1,1);
        % top Neumann
        q_r = -(Dz_hat(end,end-1)*q_hat_int(end) + Dz_hat(end,end-2)*q_hat_int(end-1))/Dz_hat(end,end);
    case "DbNt"
        % bottom Dirichlet
        q_l = ell;
        % top Neumann
        q_r = -(Dz_hat(end,end-1)*q_hat_int(end) + Dz_hat(end,end-2)*q_hat_int(end-1))/Dz_hat(end,end);
    case "DbDt"
        % bottom Dirichlet
        q_l = ell;
        % top Dirichlet
        q_r = ell;
    otherwise
        error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
end

%Extend to full q:
q_hat_full = [q_l; q_hat_int; q_r];

%Apply PDE to q_full:
dqdt  = (1./gamma).*(1./g_hat).*(Dz_hat*((1./g_hat).*(1./q_hat_full.^2).*(Dz_hat*q_hat_full)))  + ...
         q_hat_full .*(z_hat<zp_hat);

%Chop to interior q:
dqdt = dqdt(2:end-1);



end