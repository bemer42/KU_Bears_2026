%% Crypt Fission 2 
clear; close all; clc

% Discretize space:
N_u = 1e2;
u_north = linspace(0, pi, N_u);       
u_south = linspace(pi, 2*pi, N_u);   

N_v = 3e2;
L   = 78.8783;
v   = linspace(0, L, N_v)';

% Create Mesh Grids:
[U_mesh_n, V_mesh_n] = meshgrid(u_north, v);
[U_mesh_s, V_mesh_s] = meshgrid(u_south, v);

% Base Radius Parameters:
r_b = 41/2/pi;
r_t = 10/pi;   
a   = 0.3;   
b   = 0.15;
w   = 5;

R_v = r_b*(1 - exp(-a * V_mesh_n)) + r_t * exp(a * (V_mesh_n - L));

% Crypt Fission Animation Loop:
figure(1)
set(gcf, 'Color', 'w') 

for c = -10 : 2 : 90
    
    alpha    = 0.5 * (1 + tanh((c - V_mesh_n)/w));
    beta_val = 1 + b * exp(-(c - V_mesh_n).^2 / w^2);
    bR       = beta_val .* R_v;
    
    % North Lobe Surface Coordinates:
    rn0_x = bR .* cos(U_mesh_n);
    rn0_y = bR .* sin(U_mesh_n);
    rn1_x = 2 * bR .* sin(U_mesh_n) .* cos(U_mesh_n);
    rn1_y = 2 * bR .* sin(U_mesh_n) .* sin(U_mesh_n);
    
    Dn_X = (1 - alpha) .* rn0_x + alpha .* rn1_x;
    Dn_Y = (1 - alpha) .* rn0_y + alpha .* rn1_y;
    Dn_Z = V_mesh_n;
    
    % South Lobe Surface Coordinates:
    rs0_x = bR .* cos(U_mesh_s);
    rs0_y = bR .* sin(U_mesh_s);
    rs1_x = -2 * bR .* sin(U_mesh_s) .* cos(U_mesh_s);
    rs1_y = -2 * bR .* sin(U_mesh_s) .* sin(U_mesh_s);
    
    Ds_X = (1 - alpha) .* rs0_x + alpha .* rs1_x;
    Ds_Y = (1 - alpha) .* rs0_y + alpha .* rs1_y;
    Ds_Z = V_mesh_s;
    
    % Extract North Cell Chain (along apex at u = pi/2):
    [~, idx_n] = min(abs(u_north - pi/2));
    Rn_x = Dn_X(:, idx_n);
    Rn_y = Dn_Y(:, idx_n);
    Rn_z = Dn_Z(:, idx_n);
    
    % Extract South Cell Chain (along apex at u = 3*pi/2):
    [~, idx_s] = min(abs(u_south - 3*pi/2));
    Rs_x = Ds_X(:, idx_s);
    Rs_y = Ds_Y(:, idx_s);
    Rs_z = Ds_Z(:, idx_s);
   
    mesh(Dn_X, Dn_Y, Dn_Z, 'EdgeColor', [.5 .5 .5], 'FaceAlpha', 0.1)
    hold on
    mesh(Ds_X, Ds_Y, Ds_Z, 'EdgeColor', [.5 .5 .5], 'FaceAlpha', 0.1)
    
    % Plot 1D Cell Chains running along North and South walls:
    plot3(Rn_x, Rn_y, Rn_z, 'k', 'LineWidth', 4)
    plot3(Rs_x, Rs_y, Rs_z, 'k', 'LineWidth', 4)
    hold off
    
    % Formatting and LaTeX Styling:
    set(gca, 'FontSize', 14)
    title('Crypt Fission Geometry (Pinch Model)', 'FontSize', 20, 'Interpreter', 'latex')
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex')
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex')
    zlabel('$z$', 'FontSize', 14, 'Interpreter', 'latex')
    
    legend('$\vec{X}_{\mathrm{North}}(\theta,z)$', '$\vec{X}_{\mathrm{South}}(\theta,z)$', ...
           '$\vec{R}_{\mathrm{North}}(z)$', '$\vec{R}_{\mathrm{South}}(z)$', ...
           'FontSize', 12, 'Interpreter', 'latex', 'Location', 'northeast')
    
    grid on
    grid minor
    box on
    view([1 -1 1])
    xlim([-30 30])
    ylim([-30 30])
    zlim([-2 80])
    
    if c == -10
        pause(0.5)
    else
        drawnow
    end
end