%% Crypt fission 
% clear; close all; clc
% 
% % Discretize space:
% N_u = 1e2;
% u = linspace(0, 2*pi, N_u);
% du = u(2) - u(1);
% 
% N_v = 3e2;
% v = linspace(0,1,N_v)';
% dv = v(2)-v(1);
% 
% %Create mesh grid:
% [U_mesh,V_mesh] = meshgrid(u,v);
% 
% % Define radius function:
% L = 78;
% r_b = 41/2/pi;
% r_t = 10/pi;   
% a = 0.3;   
% beta = 0.15;
% 
% 
% max_separation = 13;
% 
% Z_mesh = L * V_mesh;
% R_z = r_b*(1 - exp(-a * Z_mesh)) + r_t * exp(a * (Z_mesh - L));
% 
% 
% 
% for c = -30 : 150
% 
%     A = 1./ (1+exp(beta * (Z_mesh-c)));
% 
%     if c < 0
%         A = A * exp(c); 
%     end
% 
%     %g = 1 + 1 ./ (1 + exp(beta * (Z_mesh - c)));
%     g = 1 + (sqrt(2) - 1) ./ (1 + exp(beta* (Z_mesh-c)));
% 
%     %p_u = sin(U_mesh);
%     %f = 1 - A + 2*A .* (p_u.^2);
%     %f = 1 - A .* (cos(U_mesh/2).^2);
% 
%     %lambda = f .* g .* R_z;
%     lambda = g .* R_z;
% 
%     shift = max_separation * A;
% 
%     X1 = lambda .* cos(U_mesh) - shift;
%     Y1 = lambda .* sin(U_mesh);
%     Z = Z_mesh;
% 
%     X2 = lambda .* cos(U_mesh) + shift;
%     Y2 = lambda .* sin(U_mesh);
% 
% 
% 
%     % Generate surface plot
%     figure(1)
%     surf(X1, Y2, Z, lambda)
%     hold on
% 
%     surf(X2, Y2, Z, lambda)
%     hold off
% 
%     shading interp
%     grid on
%     grid minor
%     colormap winter
%     axis equal
%     xlim([-45 45])
%     ylim([-45 45])
%     zlim([-2 80])
%     if c == -30
%         pause(1)
%     else
%          drawnow
%     end
% 
% 
% end

clear; close all; clc

% Discretize space
Nth    = 8e1;
th_0   = 0;
th_end = 2*pi;
th     = linspace(th_0, th_end, Nth);
dth    = th(2) - th(1);

Nz  = 8e1;
z_0 = 0;
L   = 78.8783; 
z   = linspace(z_0, L, Nz)';
dz  = z(2) - z(1); 

%Create Mesh Grid:
[Th, Z_mesh] = meshgrid(th, z);

%Base Radius Parameters:
rb   = 41/2/pi;
rt   = 10/pi;   
a    = 0.3;   
beta = 0.15;
max_separation = 13;


R_z = rb*(1 - exp(-a * Z_mesh)) + rt * exp(a * (Z_mesh - L));

% Crypt Fission Animation Loop:
figure(1)
set(gcf, 'Color', 'w') % Clean white background

for c = -30 : 2 : 40
   
    A = 1 ./ (1 + exp(beta * (Z_mesh - c)));
    if c < 0
        A = A * exp(c); 
    end

    % Geometry growth factor
    g = 1 + (sqrt(2) - 1) ./ (1 + exp(beta * (Z_mesh - c)));
    lambda = g .* R_z;
    shift  = max_separation * A;

    % Parametric coordinates for both splitting crypt branches
    X1 = lambda .* cos(Th) - shift;
    X2 = lambda .* cos(Th) + shift;
    Y  = lambda .* sin(Th);
    Z  = Z_mesh;

    % 1D Cell chain R(z) positioned along theta = 0 on Branch 1
    R1_x = lambda(:, 1) .* cos(0) - shift(:, 1);
    R1_y = lambda(:, 1) .* sin(0);
    R1_z = z;

    % 1D Cell chain R(z) positioned along theta = 0 on Branch 2
    R2_x = lambda(:, 1) .* cos(0) + shift(:, 1);
    R2_y = lambda(:, 1) .* sin(0);
    R2_z = z;

    % Plot Surfaces with grey mesh formatting
    mesh(X1, Y, Z, 'EdgeColor', [.5 .5 .5], 'FaceAlpha', 0.1)
    hold on
    mesh(X2, Y, Z, 'EdgeColor', [.5 .5 .5], 'FaceAlpha', 0.1)

    % Plot the 1D Cell Chains R(z)
    plot3(R1_x, R1_y, R1_z, 'k', 'LineWidth', 4)
    plot3(R2_x, R2_y, R2_z, 'k', 'LineWidth', 4)
    hold off

    % Styling and Formatting matching Figure 1
    set(gca, 'FontSize', 14)
    title('Crypt Fission Geometry', 'FontSize', 20, 'Interpreter', 'latex')
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex')
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex')
    zlabel('$z$', 'FontSize', 14, 'Interpreter', 'latex')

    legend('$\vec{X}_1(\theta,z)$', '$\vec{X}_2(\theta,z)$', '$\vec{R}_1(z)$', '$\vec{R}_2(z)$', ...
        'FontSize', 12, 'Interpreter', 'latex', 'Location', 'northeast')

    grid on 
    grid minor
    box on
    view([1 -1 1])
    xlim([-45 45])
    ylim([-45 45])
    zlim([-2 80])

    if c == -30
        pause(0.5)
    else
        drawnow
    end
end