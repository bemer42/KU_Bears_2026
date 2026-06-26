%% Crypt Geometry
clear, close all, clc

% Discretize theta space:
N_t   = 1e2;
t_0   = 0;
t_end = 2*pi;
t     = linspace(t_0, t_end, N_t)';

% Discretize z space:
N_z = 3e2;
z_0 = 0;
L   = 78;
z   = linspace(z_0, L, N_z)';

% Create a Meshgrid:
[T_mesh,Z_mesh] = meshgrid(t,z);

% Define radius function: 
rb = 41/2/pi;
rt = 10/pi;
a  = 0.45;

r = @(t,z) rb*(1-exp(-a*z)) + rt*exp(a*(z-L)); 

% Define surface generating function: 
X = @(t,z) r(t,z).*cos(t);
Y = @(t,z) r(t,z).*sin(t); 
Z = @(t,z) z;

% Generate surface plot: 
figure(1)
surf(X(T_mesh,Z_mesh),Y(T_mesh,Z_mesh),Z(T_mesh,Z_mesh))
shading interp
grid on 
grid minor
colormap winter
xlim([-30 30])
ylim([-30 30])
zlim([-2 80])

% Surface Area: 




