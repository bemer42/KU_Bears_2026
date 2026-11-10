%% ND_2D_Analysis_IsoSurfaces.m
% 3-parameter scan of the ND 2D model over:
%   X = log10(gamma_theta), Y = z_p_hat, Z = log10(gamma_z)
% and visualization via isosurfaces for C, P, R plus feasible-region mask.
%
% REQUIREMENT:
%   ND_2D_Function(z0_hat, gamma_z, gamma_theta, z_p_hat, bcType)
% returns [C, P, R]
%
% IMPORTANT:
%   This version avoids the "Input grid is not a valid MESHGRID" error by:
%     - using x,y,z as VECTORS (not full ndgrid arrays) for isosurface/isonormals
%     - permuting metric arrays into MATLAB's meshgrid ordering: V(y,x,z)

clear; close all; clc

% -------------------------
% Fixed parameters
% -------------------------
L      = 79.3515;
z0_hat = 2.11/L;
bcType = "NbDt";

% -------------------------
% Parameter vectors (set these)
% -------------------------
N_gth = 3;     % gamma_theta samples
N_zp  = 3;     % z_p_hat samples
N_gz  = 3;     % gamma_z samples

gth_vec   = logspace(-0.6, -0.2, N_gth);      % gamma_theta
zphat_vec = linspace(20, 30, N_zp)/L;         % z_p_hat
gz_vec    = logspace(-0.8, -0.2, N_gz);        % gamma_z

% For plotting axes:
x = log10(gth_vec);     % X-axis: log10(gamma_theta)
y = zphat_vec;          % Y-axis: z_p_hat
z = log10(gz_vec);      % Z-axis: log10(gamma_z)

% -------------------------
% Preallocate metric arrays:
% C(i,j,k) = metric at (gth_i, zphat_j, gz_k)
% -------------------------
C = nan(N_gth, N_zp, N_gz);
P = nan(N_gth, N_zp, N_gz);
R = nan(N_gth, N_zp, N_gz);

% -------------------------
% Parameter scan
% -------------------------
tStart = tic;
for k = 1:N_gz
    gamma_z = gz_vec(k);

    for i = 1:N_gth
        gamma_theta = gth_vec(i);

        for j = 1:N_zp
            zphat = zphat_vec(j);

            try
                [C(i,j,k), P(i,j,k), R(i,j,k)] = ND_2D_Function( ...
                    z0_hat, gamma_z, gamma_theta, zphat, bcType);
            catch ME
                warning('ND_2D_Function failed at (k,i,j)=(%d,%d,%d): %s', k, i, j, ME.message);
                C(i,j,k) = NaN;
                P(i,j,k) = NaN;
                R(i,j,k) = NaN;
            end
        end
    end

    fprintf('Finished gamma_z %d/%d (gamma_z=%g). Elapsed %.1f s\n', ...
        k, N_gz, gamma_z, toc(tStart));
end
fprintf('Total runtime: %.2f s\n', toc(tStart));

% -------------------------
% OPTIONAL: save
% -------------------------
save('Analysis_2D_3Dscan.mat','C','P','R','gth_vec','zphat_vec','gz_vec','L','z0_hat','bcType')

%% -------------------------
% Convert arrays into meshgrid ordering expected by isosurface(x,y,z,V,...)
% MATLAB expects V to be size: [length(y), length(x), length(z)]
% Your storage is size:        [length(x), length(y), length(z)]
% So permute: (x,y,z) -> (y,x,z)
% -------------------------
Cmesh = permute(C, [2 1 3]);   % size = [N_zp, N_gth, N_gz]
Pmesh = permute(P, [2 1 3]);
Rmesh = permute(R, [2 1 3]);

% Handle NaNs (isosurface/isonormals can fail with NaNs).
% With a tiny 3x3x3 grid, even one NaN can break plotting.
Cmesh(isnan(Cmesh)) = 0;
Pmesh(isnan(Pmesh)) = 0;
Rmesh(isnan(Rmesh)) = 0;

% -------------------------
% Iso-values (edit to match your targets)
% -------------------------
C_CON = [3000 3100];
P_CON = [800 1000];
R_CON = [4.9 5.1];

% -------------------------
% Plot: Total Cells isosurfaces
% -------------------------
figure('Color','w'); hold on
title('Isosurfaces: Total Cells $\mathcal{C}$','Interpreter','latex')

pC1 = patch(isosurface(x, y, z, Cmesh, C_CON(1)));
set(pC1,'FaceColor',[0.7 0.0 0.9],'EdgeColor','none','FaceAlpha',0.25);
% normals for lighting (safe now because x,y,z are vectors):
isonormals(x, y, z, Cmesh, pC1);

pC2 = patch(isosurface(x, y, z, Cmesh, C_CON(2)));
set(pC2,'FaceColor',[0.1 0.1 0.1],'EdgeColor','none','FaceAlpha',0.25);
isonormals(x, y, z, Cmesh, pC2);

xlabel('$\log_{10}(\gamma_\theta)$','Interpreter','latex')
ylabel('$\hat z_p$','Interpreter','latex')
zlabel('$\log_{10}(\gamma_z)$','Interpreter','latex')
legend({['$\mathcal{C}=' num2str(C_CON(1)) '$'], ['$\mathcal{C}=' num2str(C_CON(2)) '$']}, ...
    'Interpreter','latex','Location','best')
grid on; box on; view(3); camlight; lighting gouraud

%% -------------------------
% Plot: Proliferative Cells isosurfaces
% -------------------------
figure('Color','w'); hold on
title('Isosurfaces: Proliferative Cells $\mathcal{P}$','Interpreter','latex')

pP1 = patch(isosurface(x, y, z, Pmesh, P_CON(1)));
set(pP1,'FaceColor',[0.0 0.5 1.0],'EdgeColor','none','FaceAlpha',0.25);
isonormals(x, y, z, Pmesh, pP1);

pP2 = patch(isosurface(x, y, z, Pmesh, P_CON(2)));
set(pP2,'FaceColor',[0.1 0.1 0.1],'EdgeColor','none','FaceAlpha',0.25);
isonormals(x, y, z, Pmesh, pP2);

xlabel('$\log_{10}(\gamma_\theta)$','Interpreter','latex')
ylabel('$\hat z_p$','Interpreter','latex')
zlabel('$\log_{10}(\gamma_z)$','Interpreter','latex')
legend({['$\mathcal{P}=' num2str(P_CON(1)) '$'], ['$\mathcal{P}=' num2str(P_CON(2)) '$']}, ...
    'Interpreter','latex','Location','best')
grid on; box on; view(3); camlight; lighting gouraud

%% -------------------------
% Plot: Renewal Time isosurfaces
% -------------------------
figure('Color','w'); hold on
title('Isosurfaces: Crypt Renewal Time $\mathcal{R}$','Interpreter','latex')

pR1 = patch(isosurface(x, y, z, Rmesh, R_CON(1)));
set(pR1,'FaceColor',[1.0 0.3 0.0],'EdgeColor','none','FaceAlpha',0.25);
isonormals(x, y, z, Rmesh, pR1);

pR2 = patch(isosurface(x, y, z, Rmesh, R_CON(2)));
set(pR2,'FaceColor',[0.1 0.1 0.1],'EdgeColor','none','FaceAlpha',0.25);
isonormals(x, y, z, Rmesh, pR2);

xlabel('$\log_{10}(\gamma_\theta)$','Interpreter','latex')
ylabel('$\hat z_p$','Interpreter','latex')
zlabel('$\log_{10}(\gamma_z)$','Interpreter','latex')
legend({['$\mathcal{R}=' num2str(R_CON(1)) '$'], ['$\mathcal{R}=' num2str(R_CON(2)) '$']}, ...
    'Interpreter','latex','Location','best')
grid on; box on; view(3); camlight; lighting gouraud

%% -------------------------
% Plot: Feasible region as a 3D mask (intersection of constraints)
% -------------------------
C_OK = (Cmesh >= C_CON(1)) & (Cmesh <= C_CON(2));
P_OK = (Pmesh >= P_CON(1)) & (Pmesh <= P_CON(2));
R_OK = (Rmesh >= R_CON(1)) & (Rmesh <= R_CON(2));
FEAS = C_OK & P_OK & R_OK;

figure('Color','w'); hold on
title('Feasible Region (intersection of constraints)','Interpreter','latex')

F = double(FEAS);
% Only plot if there is something feasible
if any(FEAS(:))
    pf = patch(isosurface(x, y, z, F, 0.5));
    set(pf,'FaceColor',[0 0 0],'EdgeColor','none','FaceAlpha',0.25);
    isonormals(x, y, z, F, pf);
    camlight; lighting gouraud
else
    text(0.1,0.5,'No feasible points found on this coarse grid','Units','normalized');
end

xlabel('$\log_{10}(\gamma_\theta)$','Interpreter','latex')
ylabel('$\hat z_p$','Interpreter','latex')
zlabel('$\log_{10}(\gamma_z)$','Interpreter','latex')
grid on; box on; view(3)

fprintf('Feasible points: %d / %d\n', nnz(FEAS), numel(FEAS));