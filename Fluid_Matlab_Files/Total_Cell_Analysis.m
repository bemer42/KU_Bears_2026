%% Relationship between total cells and k,zu

close all; clear; clc

% k vector:
Nk = 10;
k_vec = logspace(-10,-5, Nk);

% zu vector
Nzu = 10;
zu_vec = linspace(10, 80, Nzu);

TC_Storage = zeros(Nk,Nzu);


% nested loop


for i = 1: Nk

    for j = 1 : Nzu
        TC_Storage(i, j) = cell2D_fun(k_vec(i),zu_vec(j));

    end

end

% Mesh grid:

[ZU_mesh, K_mesh] = meshgrid(zu_vec, k_vec);

%%
figure(1)
surf(ZU_mesh, log10(K_mesh), TC_Storage)
shading interp
colormap turbo
colorbar
grid on
grid minor

xlabel('Upper proliferation boundary (zu)', 'Font', 'Times New Roman', 'FontSize',14)
ylabel('k', 'Font', 'Times New Roman','FontSize',14)
zlabel('Total cell count', 'Font', 'Times New Roman','FontSize',14)



