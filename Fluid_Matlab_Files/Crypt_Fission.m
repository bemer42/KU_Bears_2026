clear; close all; clc

% Discretize space
Nth    = 8e1;
th_0   = 0;
th_end = 2*pi;
th     = linspace(th_0, th_end, Nth);

Nz  = 8e1;
z_0 = 0;
L   = 78.8783; 
z   = linspace(z_0, L, Nz)';

% Create Mesh Grid:
[Th, Z_mesh] = meshgrid(th, z);

% Base Radius Parameters:
rb   = 18/pi;
rt   = 10/pi;   
a    = 0.3;   
beta = 0.15;
max_separation = 13;

R_z = rb*(1 - exp(-a * Z_mesh)) + rt * exp(a * (Z_mesh - L));

%Plot
figure(1)
set(gcf, 'Color', 'w') 

for c = -30 : 2 : 35
   
    A = 1 ./ (1 + exp(beta * (Z_mesh - c)));
    if c < 0
        A = A * exp(c); 
    end

    % Geometry growth factor g(z, c) and shift function S(z, c)
    g = 1 + (sqrt(2) - 1) ./ (1 + exp(beta * (Z_mesh - c)));
    lambda = g .* R_z;
    shift  = max_separation * A;

    % Branch 1: flatten right side facing center
    pinch1 = 1 - 0.85 * (shift / max_separation) .* max(0, cos(Th));
    % Branch 2: flatten left side facing center
    pinch2 = 1 - 0.85 * (shift / max_separation) .* max(0, -cos(Th));

    % Parametric coordinates for non-overlapping crypt branches
    X1 = (lambda .* pinch1) .* cos(Th) - shift;
    X2 = (lambda .* pinch2) .* cos(Th) + shift;
    Y1 = (lambda .* pinch1) .* sin(Th);
    Y2 = (lambda .* pinch2) .* sin(Th);
    Z  = Z_mesh;

    % Single 1D Cell chain R1(z) along theta = 0 on Branch 1
    R1_x = X1(:, 1);
    R1_y = Y1(:, 1);
    R1_z = z;

    % Plot Surfaces with grey mesh formatting
    mesh(X1, Y1, Z, 'EdgeColor', [.5 .5 .5], 'FaceAlpha', 0.1)
    hold on
    mesh(X2, Y2, Z, 'EdgeColor', [.5 .5 .5], 'FaceAlpha', 0.1)

    % Plot the SINGLE 1D Cell Chain R1(z)
    plot3(R1_x, R1_y, R1_z, 'k', 'LineWidth', 4)
    hold off

    % Styling and Formatting
    set(gca, 'FontSize', 14)
    title('Crypt Fission Geometry', 'FontSize', 20, 'Interpreter', 'latex')
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex')
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex')
    zlabel('$z$', 'FontSize', 14, 'Interpreter', 'latex')
    legend('$\vec{X}_1(\theta,z)$', '$\vec{X}_2(\theta,z)$', '$\vec{R}_1(z)$', ...
        'FontSize', 12, 'Interpreter', 'latex', 'Location', 'northeast')
    grid on; grid minor; box on;
    view([1 -1 1])
    xlim([-45 45]); ylim([-45 45]); zlim([-2 80])

    if c == -30
        pause(0.5)
    else
        drawnow
    end
end


clf; % Clear the final frame with overlapping inner walls and R1

% Wrap angles into [0, 2*pi] to mask out inner-facing coordinates
Th_deg = mod(Th, 2*pi);

% Branch 1 (-x branch): Retain outer facing angles (pi/2 to 3*pi/2)
mask1 = (Th_deg >= pi/2) & (Th_deg <= 3*pi/2);
X1_outer = X1; Y1_outer = Y1;
X1_outer(~mask1) = NaN; 
Y1_outer(~mask1) = NaN;

% Branch 2 (+x branch): Retain outer facing angles (-pi/2 to pi/2)
mask2 = (Th_deg <= pi/2) | (Th_deg >= 3*pi/2);
X2_outer = X2; Y2_outer = Y2;
X2_outer(~mask2) = NaN; 
Y2_outer(~mask2) = NaN;

% Render non-overlapping outer envelope
mesh(X1_outer, Y1_outer, Z, 'EdgeColor', [.4 .4 .4], 'FaceAlpha', 0.12)
hold on
mesh(X2_outer, Y2_outer, Z, 'EdgeColor', [.4 .4 .4], 'FaceAlpha', 0.12)

% Re-apply formatting (R1 line and its legend handle are excluded)
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex')
title('Final Unified Crypt Fission Geometry', 'FontSize', 20, 'Interpreter', 'latex')
xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex')
zlabel('$z$', 'FontSize', 14, 'Interpreter', 'latex')
legend('$\vec{X}_1(\theta,z)$', '$\vec{X}_2(\theta,z)$', ...
    'FontSize', 12, 'Interpreter', 'latex', 'Location', 'northeast')
grid on; grid minor; box on;
view([1 -1 1])
xlim([-45 45]); ylim([-45 45]); zlim([-2 80])