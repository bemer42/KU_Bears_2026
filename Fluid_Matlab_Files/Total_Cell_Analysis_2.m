% Total Cell Analysis 2
close all; clear; clc

% gamma vec
Ngamma = 10;
gamma_vec = logspace(-3, -1, Ngamma); 

% zp vec
Nzp = 10;
zp_vec = linspace(0.1, 0.9, Nzp);   

TC_Storage = zeros(Ngamma, Nzp);
bcType = "DbDt";                     

% nested loop
for i = 1:Ngamma
    for j = 1:Nzp
        
        TC_Storage(i, j) = nd_cell_movement(gamma_vec(i), zp_vec(j), bcType);
    end
end

% Create meshgrid
[ZP_mesh, GAMMA_mesh] = meshgrid(zp_vec, gamma_vec);

% Plot
figure(1)

surf(ZP_mesh, log10(GAMMA_mesh), TC_Storage)

xlabel('ZP')
ylabel('gamma')
zlabel('Total Cells')
shading interp
colormap turbo
colorbar
grid on
grid minor