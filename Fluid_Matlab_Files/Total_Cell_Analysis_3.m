%% Total Cell Analysis 3
close all; clear; clc

Nalpha = 5;
alpha_z_vec = logspace(-6, -4, Nalpha); 

Nzp = 5;
zp_vec = linspace(10, 70, Nzp);        

TC_Storage = zeros(Nalpha, Nzp);
bcType = "NbDt";

% nested for
for i = 1:Nalpha
    for j = 1:Nzp
        TC_Storage(i, j) = stem_cell_movement(alpha_z_vec(i), zp_vec(j), bcType);
    end
end

[ZP_mesh, ALPHA_mesh] = meshgrid(zp_vec, alpha_z_vec);

% Plot
figure(1)
surf(ZP_mesh, log10(ALPHA_mesh), TC_Storage)
xlabel('ZP')
ylabel('alpha_z')
zlabel('Total Cells')
shading interp
colormap turbo
colorbar
grid on
grid minor;