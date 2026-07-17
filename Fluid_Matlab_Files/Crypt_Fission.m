% Crypt fission 
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


max_separation = 13;

Z_mesh = L * V_mesh;
R_z = r_b*(1 - exp(-a * Z_mesh)) + r_t * exp(a * (Z_mesh - L));



for c = -30 : 150

    A = 1./ (1+exp(beta * (Z_mesh-c)));

    if c < 0
        A = A * exp(c); 
    end

    %g = 1 + 1 ./ (1 + exp(beta * (Z_mesh - c)));
    g = 1 + (sqrt(2) - 1) ./ (1 + exp(beta* (Z_mesh-c)));

    %p_u = sin(U_mesh);
    %f = 1 - A + 2*A .* (p_u.^2);
    %f = 1 - A .* (cos(U_mesh/2).^2);

    %lambda = f .* g .* R_z;
    lambda = g .* R_z;

    shift = max_separation * A;

    X1 = lambda .* cos(U_mesh) - shift;
    Y1 = lambda .* sin(U_mesh);
    Z = Z_mesh;

    X2 = lambda .* cos(U_mesh) + shift;
    Y2 = lambda .* sin(U_mesh);
   
   

    % Generate surface plot
    figure(1)
    surf(X1, Y2, Z, lambda)
    hold on

    surf(X2, Y2, Z, lambda)
    hold off

    shading interp
    grid on
    grid minor
    colormap winter
    axis equal
    xlim([-45 45])
    ylim([-45 45])
    zlim([-2 80])
    if c == -30
        pause(1)
    else
         drawnow
    end
    

end