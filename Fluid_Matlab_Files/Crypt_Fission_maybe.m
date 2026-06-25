% Crypt fission maybe
clear; close all; clc

% Discretize space:
N_u = 1e2;
u = linspace(0, 2*pi, N_u);
du = u(2) - u(1);

N_v = 3e2;
v = linspace(0,1,N_v)';
dv = v(2)-v(1);

%Create mesh grid:
[U_mesh,V_mesh] = meshgrid(u,v);

% Define radius function:
L = 78;
r_b = 41/2/pi;
r_t = 10/pi;   
a = 0.3;   
beta = 0.15;

Z_mesh = L * V_mesh;
R_z = r_b*(1 - exp(-a * Z_mesh)) + r_t * exp(a * (Z_mesh - L));

for c = 0 : 78

    A = 1./ (1+exp(beta * (Z_mesh-c)));

    p_u = sin(U_mesh);
    f = 1 - A + 2*A .* (p_u.^2);

    lambda = f .* R_z;

    X = lambda .* cos(U_mesh);
    Y = lambda .* sin(U_mesh);
    Z = Z_mesh;

    % Generate surface plot
    figure(1)
    surf(X, Y, Z, lambda)
    shading interp
    grid on
    grid minor
    colormap winter
    axis equal
    if c == 1
        pause
    else
        drawnow
    end
    xlim([-30 30])
    ylim([-30 30])
    zlim([-2 80])

end