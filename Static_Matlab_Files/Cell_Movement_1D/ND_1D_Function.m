function [C,P,R] = ND_1D_Function(gamma,zp_hat,bcType)

% Predefine structures:
geom_hat    = struct();
diffmat_hat = struct();
par_hat     = struct();

% Define parameters:
Tc     = 15.1594; 
L      = 78.8783; 
ell    = 1; 

% Define geometry parameters:
r_b     = 41/2/pi;
r_t     = 10/pi;
a       = 0.3;
epsilon = r_b/L; 

% Discretize time:
Nt        = 1e3;
t_hat_0   = 0;
t_hat_end = 5e5*log(2)./Tc;
t_hat     = linspace(t_hat_0,t_hat_end,Nt);

% Discretize space:
Nz        = 1e3;
z_hat_0   = 0;
z_hat_end = 1;
z_hat     = linspace(z_hat_0, z_hat_end, Nz)';
dz_hat    = z_hat(2) - z_hat(1);

% Define geometry:
r_hat  = (1 - exp(-a*L*z_hat)) + r_t/r_b * exp(a*L*(z_hat-1));
rz_hat = a*L*exp(-a*L*z_hat) + a*L*r_t/r_b*exp(a*L*(z_hat-1));
g_hat  = sqrt(1+(epsilon*rz_hat).^2);

% Define the differentiation matrix:
Dz = diag(1/2*ones(Nz-1,1),1) - diag(1/2*ones(Nz-1,1),-1);
Dz(1,1:3) = [-3/2 2 -1/2];
Dz(end, end-2:end) = [1/2 -2 3/2];
Dz_hat = Dz/dz_hat;

% Initial condition:
q_hat_0     = ones(size(z_hat));
q_hat_0_int = q_hat_0(2:end-1);

% Build structures for solver: 
diffmat_hat.Dz_hat = Dz_hat;

geom_hat.z_hat   = z_hat;
geom_hat.g_hat   = g_hat;
geom_hat.epsilon = epsilon;

par_hat.gamma  = gamma;
par_hat.zp_hat = zp_hat;
par_hat.ell    = ell;

%Define right hand side function:
dQhatdt = @(t,q_hat_int) ND_dqdt_1D_snipsnap(t, q_hat_int, diffmat_hat, geom_hat, par_hat, bcType);

%Solve the system of ODEs that represents the PDE:
options = odeset('Stats', 'off');
[t_hat,Qhat_int] = ode15s(dQhatdt, t_hat, q_hat_0_int, options);

% Extend to full Q:
Qhat_full = zeros(Nz,Nt);
Vhat_full = zeros(Nz,Nt);
for i = 1:Nt

    switch bcType
        case "NbDt"
            % bottom Neumann
            Q_l = -(Dz(1,2)*Qhat_int(i,1) + Dz(1,3)*Qhat_int(i,2))/Dz(1,1);
            % top Dirichlet
            Q_r = ell;
        case "NbNt"
            % bottom Neumann
            Q_l = -(Dz(1,2)*Qhat_int(i,1) + Dz(1,3)*Qhat_int(i,2))/Dz(1,1);
            % top Neumann
            Q_r = -(Dz(end,end-1)*Qhat_int(i,end) + Dz(end,end-2)*Qhat_int(i,end-1))/Dz(end,end);
        case "DbNt"
            % bottom Dirichlet
            Q_l = ell;
            % top Neumann
            Q_r = -(Dz(end,end-1)*Qhat_int(i,end) + Dz(end,end-2)*Qhat_int(i,end-1))/Dz(end,end);
        case "DbDt"
            % bottom Dirichlet
            Q_l = ell;
            % top Dirichlet
            Q_r = ell;
        otherwise
            error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
    end

    Qhat_full(:,i) = [Q_l; Qhat_int(i, :)'; Q_r];
    Vhat_full(:,i) = -(1./par_hat.gamma).*(1./g_hat).*(1./Qhat_full(:,i).^3).*(Dz_hat*Qhat_full(:,i));

end

% Steady State Analysis: 
Qhat_ss = Qhat_full(:,end);
Vhat_ss = Vhat_full(:,end);

% Total number of cells:
C = (L.^2*epsilon./ell).*2*pi*trapz(z_hat,r_hat.*Qhat_ss.*g_hat);
P = (L.^2*epsilon./ell).*2*pi*trapz(z_hat,r_hat.*Qhat_ss.*(z_hat<zp_hat).*g_hat);

% Crypt renewal time:
pos = find(z_hat>=1/L);
R = (Tc/log(2)).*trapz(z_hat(pos),g_hat(pos)./Vhat_ss(pos))/24;

end



