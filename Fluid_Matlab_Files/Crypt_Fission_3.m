%% Crypt Fission (Cassini?)
clear; close all; clc

% Discretize space
N_u = 120; 
u = linspace(0, 2*pi, N_u);
N_v = 120; 
v = linspace(0, 1, N_v)';
[U_mesh, V_mesh] = meshgrid(u, v);

% Crypt Geometry
L = 78;                
r_b = 41 / (2 * pi); 
r_t = 10 / pi;       
alpha = 0.3;       

Z = L * V_mesh;
R_z = r_b * (1 - exp(-alpha * Z)) + r_t * exp(alpha * (Z - L));


% Plot
figure(1)
for c = -10 : 120
   
    A = 1 ./ (1 + exp(0.15 * (Z - c)));
    if c < 0
        A = A * exp(c); 
    end 
    
   
    b = R_z; 
    a = 0.90 * b .* A; 
   % b = sqrt(R_z.^2 + R_z.*a);

    
    % Cassini
    r_sq = a.^2 .* cos(2 * U_mesh) + sqrt(b.^4 - a.^4 .* sin(2 * U_mesh).^2);
    lambda = sqrt(r_sq);
    
    % Cartesian coordinates
    X = lambda .* cos(U_mesh);
    Y = lambda .* sin(U_mesh);
    
    % Surface
    surf(X, Y, Z, lambda)
    shading interp
    colormap winter
    colorbar
    axis equal
    view(3)
    xlim([-20 20])
    ylim([-20 20])
    zlim([-2 80])
    drawnow
end