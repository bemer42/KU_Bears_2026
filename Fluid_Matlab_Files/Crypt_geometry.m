% Crypt Geometry
clear; close all; clc

% Discretize theta space:
N_t = 1e2;
t_0 = 0;
t_end = 2*pi;
t = linspace(t_0, t_end, N_t)';
dt = t(2) - t(1);

% Discretize z space:
N_z = 3e2;
z_0 = 0;
L = 78;
z = linspace(z_0, L, N_z)';
dz = z(2) - z(1);

%Create mesh grid:
[T_mesh,Z_mesh] = meshgrid(t,z);

% Define radius function:
r_b = 41/2/pi;
r_t = 10/pi;   
a = 0.3;   


r = @(t,z) r_b*(1 - exp(-a * z)) + r_t * exp(a * (z - L));


% Define surface generating function:
X = @(t,z) r(t,z) .* cos(t);
Y = @(t,z) r(t,z) .* sin(t);
Z = @(t,z) z;


% Generate surface plot
figure(1)
surf(X(T_mesh,Z_mesh),Y(T_mesh,Z_mesh), Z(T_mesh,Z_mesh))
shading interp
grid on
grid minor
colormap winter
xlim([-30 30])
zlim([-2 80])