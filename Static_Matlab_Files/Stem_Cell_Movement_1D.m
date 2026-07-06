%% Cell Movement 1D with Stemness
close all; clear; clc

bcType = "NbDt";

% % Parameters
par = struct();
par.alpha_z  = 1e0;
par.alpha_s  = 1e-3;
par.k        = 1e-4;     
par.s_stem   = 1.8;
par.s_ta     = .8;
par.limiter  = "vanleer";  % "vanleer", "mc", "minmod"

% Time grid
N_t   = 2e3;
t_0   = 0;
t_end = 1e4;
tspan = linspace(t_0,t_end,N_t);

% Discretize s
Ns  = 8e1;
s_0 = 0;
L_s = 3;
s   = linspace(s_0,L_s,Ns)';
ds  = s(2)-s(1);

% Discretize z
Nz  = 8e1;
z_0 = 0;
L_z = 78;
z   = linspace(z_0,L_z,Nz)';
dz  = z(2)-z(1);

% Mesh
[Z,S] = meshgrid(z,s);

% Convective Velocity in Stem dimension:
load('Average_Values')
S_avg = B_avg.*sqrt(15)./A_avg./sqrt(B_avg + 30);
s_star  = interp1(z_old,S_avg,z,'pchip','extrap');
s_star = L_s * (s_star - s_0) / (L_s - s_0 + eps);    
v0 = 1e-3;
w  = .5;
vs_fun = @(s,z) -v0 * tanh( (s - repmat(s_star.', size(s,1), 1)) / w);
VS = vs_fun(S,Z);         

% Differentiation matrices
Ds = spdiags([-1/2*ones(Ns,1) 1/2*ones(Ns,1)],[-1 1],Ns,Ns);
Ds(1,1:3) = [-3/2 2 -1/2];
Ds(end, end-2:end) = [1/2 -2 3/2];
Ds = Ds/ds;

Dz = spdiags([-1/2*ones(Nz,1) 1/2*ones(Nz,1)],[-1 1],Nz,Nz);
Dz(1,1:3) = [-3/2 2 -1/2];
Dz(end, end-2:end) = [1/2 -2 3/2];
Dz = Dz/dz;

% Curvature of crypt:
r_b = 41/2/pi;
r_t = 10/pi;
a   = 0.3;
rz_fun = @(zz) -a*r_b*exp(-a*zz) + a*r_t*exp(a*zz - L_z);

RZ = rz_fun(Z);           

% Initial condition
sigS = 10; sigZ = 8;
P0 = (1/L_s) * (1 + 0.05*exp(-((S-30).^2)/(2*sigS^2)).*exp(-((Z-15).^2)/(2*sigZ^2)));
P0(P0<=0) = 1e-12;

p0_int = P0(2:end-1,2:end-1);
p0_int = p0_int(:);

% RHS
dPdt = @(t, p_int) dpdt_1D_snipsnap(t, p_int, S, Z, RZ, VS, Ds, Dz, par, bcType);

% Solve
options = odeset('Stats','on','MaxStep',inf);

% Strongly recommended for this model (prevents tiny negative p causing q<=0):
options = odeset(options,'NonNegative',1:numel(p0_int));

tic
[t, p_int] = ode15s(dPdt, tspan, p0_int, options);
toc

% Reconstruction (q-based BCs to match RHS)
wts = ones(Ns,1); wts(1)=0.5; wts(end)=0.5;
qcol = @(pcol) ds * (wts.' * pcol);

b0 = Dz(1,1); b1 = Dz(1,2); b2 = Dz(1,3);
t2 = Dz(end,end-2); t1 = Dz(end,end-1); t0 = Dz(end,end);

a0L = Ds(1,1); a1L = Ds(1,2); a2L = Ds(1,3);
a0R = Ds(end,end); a1R = Ds(end,end-1); a2R = Ds(end,end-2);

alpha_s = par.alpha_s;

P_full = zeros(N_t, Ns*Nz);

for i = 1:N_t
    Pint = reshape(p_int(i,:), Ns-2, Nz-2);
    P = zeros(Ns, Nz);
    P(2:Ns-1, 2:Nz-1) = Pint;

    % s no-flux
    P = applyNoFlux_s(P, 2:Nz-1, alpha_s, VS, a0L,a1L,a2L,a0R,a1R,a2R);

    % z BCs (q-based)
    switch bcType
        case "NbDt"
            q2 = qcol(P(:,2));
            q3 = qcol(P(:,3));
            q1 = -(b1*q2 + b2*q3)/b0;
            P(:,1) = (q1/q2) * P(:,2);

            qNm1 = qcol(P(:,end-1));
            P(:,end) = P(:,end-1)/qNm1;

        case "NbNt"
            q2 = qcol(P(:,2));
            q3 = qcol(P(:,3));
            q1 = -(b1*q2 + b2*q3)/b0;
            P(:,1) = (q1/q2) * P(:,2);

            qNm1 = qcol(P(:,end-1));
            qNm2 = qcol(P(:,end-2));
            qN   = -(t1*qNm1 + t2*qNm2)/t0;
            P(:,end) = (qN/qNm1) * P(:,end-1);

        case "DbNt"
            q2 = qcol(P(:,2));
            P(:,1) = P(:,2)/q2;

            qNm1 = qcol(P(:,end-1));
            qNm2 = qcol(P(:,end-2));
            qN   = -(t1*qNm1 + t2*qNm2)/t0;
            P(:,end) = (qN/qNm1) * P(:,end-1);

        case "DbDt"
            q2 = qcol(P(:,2));
            P(:,1) = P(:,2)/q2;

            qNm1 = qcol(P(:,end-1));
            P(:,end) = P(:,end-1)/qNm1;

        otherwise
            error('Unknown bcType.');
    end

    % s no-flux on z-boundary columns too
    P = applyNoFlux_s(P, [1 Nz], alpha_s, VS, a0L,a1L,a2L,a0R,a1R,a2R);

    P_full(i,:) = P(:);
end

%% Animation
dtPlot = round(.01*N_t);

for i = 1:dtPlot:N_t
    P = reshape(P_full(i,:), Ns, Nz); 
    figure(2); clf
    surf(Z,S,P);
    colormap summer; grid on; grid minor;
    set(gca,'fontsize',16);
    xlabel('z','fontsize',18); 
    ylabel('s','fontsize',18); 
    zlabel('p','fontsize',18);
    title(sprintf('p(s,z,t), t = %.3g', t(i)),'fontsize',25);
    xlim([0 L_z]); ylim([0 L_s]); zlim([min(min(P_full)) max(max(P_full))])
    xaxislocation = 'origin'; 
    yaxislocation = 'origin';
    set(gca,'XDir','normal','YDir','normal')
    view([0 0 1]);
    if i == 1
        pause
    end
end

%% Helpers
function Pout = applyNoFlux_s(Pin, zCols, alpha_s, VS, a0L,a1L,a2L,a0R,a1R,a2R)
    Pout = Pin;

    vsL = VS(1,   zCols);
    vsR = VS(end, zCols);

    if alpha_s == 0
        maskL = abs(vsL) > 0;
        maskR = abs(vsR) > 0;
        Pout(1,   zCols(maskL)) = 0;
        Pout(end, zCols(maskR)) = 0;
        Pout(1,   zCols(~maskL)) = Pout(2,     zCols(~maskL));
        Pout(end, zCols(~maskR)) = Pout(end-1, zCols(~maskR));
        return
    end

    denomL = (-vsL + alpha_s*a0L);
    denomR = (-vsR + alpha_s*a0R);

    if any(abs(denomL) < 1e-14) || any(abs(denomR) < 1e-14)
        error('s no-flux enforcement ill-conditioned (driver).');
    end

    Pout(1, zCols)   = -alpha_s*(a1L*Pout(2,zCols)     + a2L*Pout(3,zCols))     ./ denomL;
    Pout(end, zCols) = -alpha_s*(a1R*Pout(end-1,zCols) + a2R*Pout(end-2,zCols)) ./ denomR;
end