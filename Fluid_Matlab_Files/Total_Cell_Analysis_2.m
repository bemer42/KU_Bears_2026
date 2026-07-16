% Total Cell Analysis 2
close all; clear; clc

% gamma vec
Ngamma = 5;
gamma_vec = logspace(-3, -1, Ngamma); 

% zp vec
Nzp = 5;
zp_vec = linspace(0.1, 0.9, Nzp);   

TC_Storage = zeros(Ngamma, Nzp);
PR_Storage = zeros(Ngamma, Nzp);
CR_Storage = zeros(Ngamma, Nzp);

bcType = "DbNt";    



% nested loop
for i = 1:Ngamma
    for j = 1:Nzp
        
        [TC_Storage(i, j), PR_Storage(i,j), CR_Storage(i,j)] = nd_cell_movement(gamma_vec(i), zp_vec(j), bcType);
    end
end

% Create meshgrid
[ZP_mesh, GAMMA_mesh] = meshgrid(zp_vec, gamma_vec);

%% Plot

figure(1)
surf(ZP_mesh, log10(GAMMA_mesh), TC_Storage)
shading interp; colormap turbo; colorbar; grid on; grid minor;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
xlabel('Upper Proliferation Boundary (z_p)', 'FontSize', 14)
ylabel('log_{10}(\gamma)', 'FontSize', 14)
zlabel('Total Cell Count', 'FontSize', 14)
% C_z_Value = max(max(TC_Storage)) + 1;
% C_TC_3300 = contourc(zp_vec, log10(gamma_vec), TC_Storage, [3300 3300]);
% zp_TC_3300 = C_TC_3300(1,2:end);
% gamma_TC_3300 = C_TC_3300(2, 2:end);
% 
% C_PR_800 = contourc(zp_vec, log10(gamma_vec), PR_Storage, [800 800]);
% zp_PR_800 = C_PR_800(1,2:end);
% gamma_PR_800 = C_PR_800(2, 2:end);
% 
% C_CR_4 = contourc(zp_vec, log10(gamma_vec), TC_Storage, [4 4]);
% zp_CR_4 = C_CR_4(1,2:end);
% gamma_CR_4 = C_CR_4(2, 2:end);
% 
% 
% 
% figure(1)
% 
% surf(ZP_mesh, log10(GAMMA_mesh), TC_Storage)
% hold on
% plot3(zp_TC_3300,gamma_TC_3300, C_z_Value*ones(size(zp_TC_3300)), 'w-', 'LineWidth',3)
% plot3(zp_PR_800,gamma_PR_800, C_z_Value*ones(size(zp_PR_800)), 'w-', 'LineWidth',3)
% plot3(zp_CR_4,gamma_CR_4, C_z_Value*ones(size(zp_CR_4)), 'w-', 'LineWidth',3)
% 
% xlabel('ZP')
% ylabel('gamma')
% zlabel('Total Cells')
% shading interp
% colormap turbo
% colorbar
% grid on
% grid minor
% view(2)
%% Total Cells
figure(2)

[zp_TC, gamma_TC] = get_contour_coords(zp_vec, log10(gamma_vec), TC_Storage, 3300);


surf(ZP_mesh, log10(GAMMA_mesh), TC_Storage)
hold on
z_TC = max(TC_Storage(:)) + 10; % Clear the top of the surface
plot3(zp_TC, gamma_TC, z_TC * ones(size(zp_TC)), 'w-', 'LineWidth', 2.5)

% Formatting
shading interp; colormap turbo; colorbar; grid on; grid minor;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
xlabel('Upper Proliferation Boundary (z_p)', 'FontSize', 14)
ylabel('log_{10}(\gamma)', 'FontSize', 14)
zlabel('Total Cell Count', 'FontSize', 14)
view(2)

%% Prolif cells
figure(3)

[zp_PR, gamma_PR] = get_contour_coords(zp_vec, log10(gamma_vec), PR_Storage, 800);
%target_PR = (min(PR_Storage(:)) + max(PR_Storage(:))) / 2; 


%C_PR = contourc(zp_vec, log10(gamma_vec), PR_Storage, [target_PR, target_PR]);

%zp_PR = C_PR(1, 2 : end);
%gamma_PR = C_PR(2, 2 :end);

% 3. Plot surface and hover the contour line
surf(ZP_mesh, log10(GAMMA_mesh), PR_Storage)
hold on
z_PR = max(PR_Storage(:)) + 10;
plot3(zp_PR, gamma_PR, z_PR * ones(size(zp_PR)), 'w-', 'LineWidth', 2.5)

% Formatting
shading interp; colormap turbo; colorbar; grid on; grid minor;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
xlabel('Upper Proliferation Boundary (z_p)', 'FontSize', 14)
ylabel('log_{10}(\gamma)', 'FontSize', 14)
zlabel('Proliferating Cells', 'FontSize', 14)
view(2)

%% Crypt Renewal
figure(4)
[zp_CR, gamma_CR] = get_contour_coords(zp_vec, log10(gamma_vec), CR_Storage, 6);
%target_CR = (min(CR_Storage(:)) + max(CR_Storage(:))) / 2; 


%C_CR = contourc(zp_vec, log10(gamma_vec), CR_Storage, [target_CR, target_CR]);
%zp_CR = C_CR(1, 2 : end);
%gamma_CR = C_CR(2, 2 :end);

% surface
surf(ZP_mesh, log10(GAMMA_mesh), CR_Storage)
hold on
z_CR = max(CR_Storage(:)) + 1; 
plot3(zp_CR, gamma_CR, z_CR * ones(size(zp_CR)), 'w-', 'LineWidth', 2.5)
shading interp; colormap turbo; colorbar; grid on; grid minor;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
xlabel('Upper Proliferation Boundary (z_p)', 'FontSize', 14)
ylabel('log_{10}(\gamma)', 'FontSize', 14)
zlabel('Crypt Renewal Time (Days)', 'FontSize', 14)
view(2)

figure(5)


surf(ZP_mesh, log10(GAMMA_mesh), TC_Storage)
hold on


[x_TC, y_TC] = get_contour_coords(zp_vec, log10(gamma_vec), TC_Storage, 3300);
[x_PR, y_PR] = get_contour_coords(zp_vec, log10(gamma_vec), PR_Storage, 800);
[x_CR, y_CR] = get_contour_coords(zp_vec, log10(gamma_vec), CR_Storage, 4);




% 4. Plot all three lines with distinct colors and styles
h = []; % Collect handles for a neat legend
if ~isempty(x_TC); h(1) = plot3(x_TC, y_TC, z_TC*ones(size(x_TC)), 'w-', 'LineWidth', 3);   end 
if ~isempty(x_PR); h(2) = plot3(x_PR, y_PR, z_TC*ones(size(x_PR)), 'm-', 'LineWidth', 3);  end 
if ~isempty(x_CR); h(3) = plot3(x_CR, y_CR, z_TC*ones(size(x_CR)), 'c-', 'LineWidth', 3);  end 


shading interp
colormap turbo
colorbar
grid on
grid minor
view(2) 

% Axes & Labels
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
xlabel('Upper Proliferation Boundary (z_p)', 'FontSize', 14)
ylabel('log_{10}(\gamma)', 'FontSize', 14)


% Clear legend explaining each line
legend_labels = {'Total Cells = 3300', 'Prolif Cells = 800', 'Crypt Renewal = 4 Days'};
legend(h(h > 0), legend_labels(h > 0), 'Location', 'northeastoutside')
