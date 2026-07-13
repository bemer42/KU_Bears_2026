%Crypt Fission 2
clear; close all; clc

% Discretize space 
N_u = 1e2;
u_north = linspace(0, pi, N_u);       
u_south = linspace(pi, 2*pi, N_u);   

N_v = 3e2;
v = linspace(0,78,N_v)';

% Create mesh grids
[U_mesh_n, V_mesh_n] = meshgrid(u_north, v);
[U_mesh_s, V_mesh_s] = meshgrid(u_south, v);

% Define radius function 
L = 78;
r_b = 41/2/pi;
r_t = 10/pi;   
a = 0.3;   
b = 0.15;
w = 5;
R_v = r_b*(1 - exp(-a * V_mesh_n)) + r_t * exp(a * (V_mesh_n - L));

for c = -10 : 90
    
    alpha = 0.5 * (1 + tanh((c - V_mesh_n)/w));
    beta_val = 1 + b * exp(-(c - V_mesh_n).^2 / w^2);
    bR = beta_val .* R_v;
    
    % North
    rn0_x = bR .* cos(U_mesh_n);
    rn0_y = bR .* sin(U_mesh_n);
    rn1_x = 2 * bR .* sin(U_mesh_n) .* cos(U_mesh_n);
    rn1_y = 2 * bR .* sin(U_mesh_n) .* sin(U_mesh_n);
    
    Dn_X = (1 - alpha) .* rn0_x + alpha .* rn1_x;
    Dn_Y = (1 - alpha) .* rn0_y + alpha .* rn1_y;
    Dn_Z = V_mesh_n;
    
    % South
    rs0_x = bR .* cos(U_mesh_s);
    rs0_y = bR .* sin(U_mesh_s);
    rs1_x = -2 * bR .* sin(U_mesh_s) .* cos(U_mesh_s);
    rs1_y = -2 * bR .* sin(U_mesh_s) .* sin(U_mesh_s);
    
    Ds_X = (1 - alpha) .* rs0_x + alpha .* rs1_x;
    Ds_Y = (1 - alpha) .* rs0_y + alpha .* rs1_y;
    Ds_Z = V_mesh_s;
    
    % Surface plot
    figure(1)
    surf(Dn_X, Dn_Y, Dn_Z, beta_val)
    hold on
    surf(Ds_X, Ds_Y, Ds_Z, beta_val)
    hold off
    
    shading interp
    grid on
    grid minor
    colormap winter
    axis equal
    xlim([-30 30])
    ylim([-30 30])
    zlim([-2 80])
    
    if c == -10
        pause(1)
    else
        drawnow
    end
end