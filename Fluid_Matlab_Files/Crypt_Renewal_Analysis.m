%% Crypt Renewal Time Surface Analysis
close all; clear; clc

% Parameter grids
Nalpha = 5;
alpha_th_vec = logspace(2, 4, Nalpha); % Range for alpha_theta

Nzp = 5;
zp_vec = linspace(10, 70, Nzp);       % Range for prolif height zp

% Preallocate storage matrix
CRT_Storage = zeros(Nalpha, Nzp);

% Nested parameter sweep loop
for i = 1:Nalpha
    for j = 1:Nzp
        % Pass alpha_th (and optional zp) to function
        CRT_Storage(i, j) = nd_cell_movement_2D_analysis(alpha_th_vec(i), zp_vec(j));
    end
end

% Create meshgrid for surface plotting
[ZP_mesh, ALPHA_mesh] = meshgrid(zp_vec, alpha_th_vec);

% -------------------------
% Plot Renewal Time Surface
% -------------------------
figure(1)
set(gcf, 'Color', 'w')

surf(ZP_mesh, log10(ALPHA_mesh), CRT_Storage)

% Formatting to match your analysis style
title('Crypt Renewal Time Surface', 'FontSize', 18, 'Interpreter', 'latex')
xlabel('Proliferative Height $z_p$', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('$\log_{10}(\alpha_{\theta})$', 'FontSize', 14, 'Interpreter', 'latex')
zlabel('Crypt Renewal Time (Days)', 'FontSize', 14, 'Interpreter', 'latex')

shading interp
colormap turbo
colorbar
grid on
grid minor
box on
view([-37.5, 30])