%% Crypt Geometry 2-Patch (North/South) -- Geometry only
% Two patches: rho_n on [0,pi], rho_s on [pi,2pi]
% Interior-in-z state, free-ring Neumann in z, seam coupling in theta.
close all; clear; clc

% Predefine structures:
geom    = struct();
diffmat = struct();
par     = struct();

% Boundary condition set (for rho in z):
bcType = "NbNt";   % Neumann bottom + Neumann top: rho_z=0 at z=0,L

% Parameters (match your q-code values for geometry)
par.r_b = 41/2/pi;
par.r_t = 10/pi;
par.a   = 0.3;
par.L   = 78.8783;

% Motion law parameters (tune)
par.beta   = 1e0;     % density forcing (set 0 while testing geometry)
par.eta    = 5e-1;    % pull-to-terminal
par.kappa  = 5e-2;    % curvature smoothing
par.rhoMin = 1e-6;

% Time discretization (same spirit)
N_t   = 2e2;
t_0   = 0;
t_end = 1e1;                 % geometry timescale (adjust)
tspan = linspace(t_0,t_end,N_t);

% Theta grids (two patches)
Nth = 3e1;                   % like your code
th_n = linspace(0, pi,   Nth)';      % includes endpoints
th_s = linspace(pi, 2*pi, Nth)';     % includes endpoints
dth  = th_n(2) - th_n(1);            % should match on both

% Z grid (same nonuniform mapping)
Nz   = 4e1;
L    = par.L;
beta = 3;
xi   = linspace(0,1,Nz)';
z    = L/2*(1+tanh(beta*(2*xi-1))/tanh(beta));  % Nz x 1

% Placeholder density field q_i(th,z,t) (replace later with your q-solver)
amp = 0.2; zp = 27; sig = 10;
par.qfun_n = @(th,z,t) (1.0 + .5*(L-z)/L);
par.qfun_s = @(th,z,t) (1.0 + .5*(L-z)/L);

% Meshes
[Z_n,Th_n] = meshgrid(z, th_n);
[Z_s,Th_s] = meshgrid(z, th_s);

% Differentiation matrix in z (same as your code)
Dz = diffmat_nonuniform(z);

% Pack diffmat
diffmat.Dz  = Dz;
diffmat.dth = dth;

% Pack geom
geom.Th_n = Th_n;  geom.Z_n = Z_n;
geom.Th_s = Th_s;  geom.Z_s = Z_s;

% Initial condition: rho_i(θ,z,0)=r(z) on BOTH patches
r_n = par.r_b*(1-exp(-par.a*Z_n)) + par.r_t*exp(par.a*(Z_n-L));
r_s = par.r_b*(1-exp(-par.a*Z_s)) + par.r_t*exp(par.a*(Z_s-L));

RHO0_n = r_n;
RHO0_s = r_s;

% interior-in-z state vectors
RHO0_n_int = RHO0_n(:,2:end-1);
RHO0_s_int = RHO0_s(:,2:end-1);

y0 = [RHO0_n_int(:); RHO0_s_int(:)];

% Solve with ode15s
rhs = @(t,y) drhodt_2D_snipsnap(t, y, diffmat, geom, par, bcType);

options = odeset('Stats','on','RelTol',1e-5,'AbsTol',1e-7,'MaxStep',inf);

tic
[t, ysol] = ode23s(rhs, tspan, y0, options);
toc

% Expand final snapshot for plotting (apply z-BC + seam coupling)
Nth = numel(th_n); Nz = numel(z); nz = Nz-2;

yend = ysol(end,:).';

rho_n_int = reshape(yend(1:Nth*nz), Nth, nz);
rho_s_int = reshape(yend(Nth*nz+1:end), Nth, nz);

rho_n = zeros(Nth,Nz); rho_n(:,2:end-1)=rho_n_int;
rho_s = zeros(Nth,Nz); rho_s(:,2:end-1)=rho_s_int;

% Apply z-Neumann and seam coupling (same routines as RHS uses)
rho_n = applyZNeumannBC(rho_n, Dz, bcType);
rho_s = applyZNeumannBC(rho_s, Dz, bcType);
[rho_n, rho_s] = enforceSeamsValueContinuity(rho_n, rho_s);

rho_n = max(rho_n, par.rhoMin);
rho_s = max(rho_s, par.rhoMin);

%% Animation: both patches on the same 3D axes
dtFrame = round(0.05*N_t);

figure(10); clf
ax = axes('Parent',gcf);
hold(ax,'on');
axis(ax,'equal');
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
view(ax,45,25);
colormap(ax,turbo);
cb = colorbar(ax);

% initialize first frame
k = 1;
yK = ysol(k,:).';

Nth = numel(th_n); Nz = numel(z); nz = Nz-2;
Dz  = diffmat.Dz;

rho_n_int = reshape(yK(1:Nth*nz), Nth, nz);
rho_s_int = reshape(yK(Nth*nz+1:end), Nth, nz);

rho_n = zeros(Nth,Nz); rho_n(:,2:end-1)=rho_n_int;
rho_s = zeros(Nth,Nz); rho_s(:,2:end-1)=rho_s_int;

rho_n = applyZNeumannBC(rho_n, Dz, bcType);
rho_s = applyZNeumannBC(rho_s, Dz, bcType);
[rho_n, rho_s] = enforceSeamsValueContinuity(rho_n, rho_s);

rho_n = max(rho_n, par.rhoMin);
rho_s = max(rho_s, par.rhoMin);

Xn = rho_n.*cos(Th_n);  Yn = rho_n.*sin(Th_n);  Zn = Z_n;
Xs = rho_s.*cos(Th_s);  Ys = rho_s.*sin(Th_s);  Zs = Z_s;

hn = surf(ax, Xn, Yn, Zn, rho_n, 'EdgeColor','none', 'FaceAlpha',1.0);
hs = surf(ax, Xs, Ys, Zs, rho_s, 'EdgeColor','none', 'FaceAlpha',1.0);

caxis(ax,[min([rho_n(:);rho_s(:)]), max([rho_n(:);rho_s(:)])]);
title(ax,sprintf('t = %.2f', t(k)));
camlight(ax); lighting(ax,'gouraud');
grid on 
grid minor
box on

for k = 1:dtFrame:N_t
    yK = ysol(k,:).';

    rho_n_int = reshape(yK(1:Nth*nz), Nth, nz);
    rho_s_int = reshape(yK(Nth*nz+1:end), Nth, nz);

    rho_n = zeros(Nth,Nz); rho_n(:,2:end-1)=rho_n_int;
    rho_s = zeros(Nth,Nz); rho_s(:,2:end-1)=rho_s_int;

    rho_n = applyZNeumannBC(rho_n, Dz, bcType);
    rho_s = applyZNeumannBC(rho_s, Dz, bcType);
    [rho_n, rho_s] = enforceSeamsValueContinuity(rho_n, rho_s);

    rho_n = max(rho_n, par.rhoMin);
    rho_s = max(rho_s, par.rhoMin);

    Xn = rho_n.*cos(Th_n);  Yn = rho_n.*sin(Th_n);
    Xs = rho_s.*cos(Th_s);  Ys = rho_s.*sin(Th_s);

    set(hn, 'XData', Xn, 'YData', Yn, 'ZData', Zn, 'CData', rho_n);
    set(hs, 'XData', Xs, 'YData', Ys, 'ZData', Zs, 'CData', rho_s);

    caxis(ax,[min([rho_n(:);rho_s(:)]), max([rho_n(:);rho_s(:)])]);
    title(ax,sprintf('t = %.2f', t(k)));

    drawnow;

    if k==1
        pause
    end
end


function rho = applyZNeumannBC(rho, Dz, bcType)
% Enforce rho_z=0 at z=0 and z=L using one-sided stencil coefficients from Dz

Nz = size(rho,2);

b0 = Dz(1,1);       b1 = Dz(1,2);       b2 = Dz(1,3);
t2 = Dz(end,end-2); t1 = Dz(end,end-1); t0 = Dz(end,end);

switch bcType
    case "NbNt"
        rho(:,1)   = -(b1*rho(:,2) + b2*rho(:,3))/b0;
        rho(:,end) = -(t1*rho(:,end-1) + t2*rho(:,end-2))/t0;
    otherwise
        error('applyZNeumannBC: only bcType="NbNt" implemented here.');
end
end

function [rho_n, rho_s] = enforceSeamsValueContinuity(rho_n, rho_s)
% Enforce rho continuity at:
%   theta=pi: north(end,:) == south(1,:)
%   wrap:     north(1,:)   == south(end,:)

% theta = pi seam
avgPi = 0.5*(rho_n(end,:) + rho_s(1,:));
rho_n(end,:) = avgPi;
rho_s(1,:)   = avgPi;

% wrap seam (0 ~ 2pi)
avgWrap = 0.5*(rho_n(1,:) + rho_s(end,:));
rho_n(1,:)   = avgWrap;
rho_s(end,:) = avgWrap;
end
