% Total Cell Analysis 2
close all; clear; clc

% gamma vec
Ngamma = 2;
gamma_vec = logspace(log10(0.01), log10(0.25), Ngamma);

% zp vec
Nzp = 2;
zp_vec = linspace(0.10, 0.58, Nzp);   

TC_Storage = zeros(Ngamma, Nzp);
PR_Storage = zeros(Ngamma, Nzp);
CR_Storage = zeros(Ngamma, Nzp);

bcType = "NbDt";    



% nested loop
for i = 1:Ngamma
    for j = 1:Nzp
        
        [TC_Storage(i, j), PR_Storage(i,j), CR_Storage(i,j)] = nd_cell_movement(gamma_vec(i), zp_vec(j), bcType);
    end
end

% Create meshgrid
[ZP_mesh, GAMMA_mesh] = meshgrid(zp_vec, gamma_vec);

%% Plot


tc_low  = 2900;  tc_high = 3200;
pr_low  = 800;   pr_high = 1000;
cr_low  = 4.5;   cr_high = 6.5;


[zp_tc_L, gamma_tc_L] = get_contour_coords(zp_vec, log10(gamma_vec), TC_Storage, tc_low);
[zp_tc_H, gamma_tc_H] = get_contour_coords(zp_vec, log10(gamma_vec), TC_Storage, tc_high);

[zp_pr_L, gamma_pr_L] = get_contour_coords(zp_vec, log10(gamma_vec), PR_Storage, pr_low);
[zp_pr_H, gamma_pr_H] = get_contour_coords(zp_vec, log10(gamma_vec), PR_Storage, pr_high);

[zp_cr_L, gamma_cr_L] = get_contour_coords(zp_vec, log10(gamma_vec), CR_Storage, cr_low);
[zp_cr_H, gamma_cr_H] = get_contour_coords(zp_vec, log10(gamma_vec), CR_Storage, cr_high);


%% Total Cells
figure(1)
surf(ZP_mesh, log10(GAMMA_mesh), TC_Storage)
shading interp; colormap turbo; colorbar; grid on; grid minor;


set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)
xlabel('Upper Proliferation Boundary ($z_p$)', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\log_{10}(\gamma)$', 'Interpreter', 'latex', 'FontSize', 14)
zlabel('Total Cell Count', 'Interpreter', 'latex', 'FontSize', 14)
view(2)

%% Total Cells Contours
figure(2)
surf(ZP_mesh, log10(GAMMA_mesh), TC_Storage)
hold on
z_TC = max(TC_Storage(:)) + 10; 


h_tc(1) = plot3(zp_tc_L, gamma_tc_L, z_TC * ones(size(zp_tc_L)), 'w--', 'LineWidth', 2.5);
h_tc(2) = plot3(zp_tc_H, gamma_tc_H, z_TC * ones(size(zp_tc_H)), 'w-', 'LineWidth', 2.5);

shading interp; colormap turbo; colorbar; grid on; grid minor;
view(2)


set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)
xlabel('Upper Proliferation Boundary ($z_p$)', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\log_{10}(\gamma)$', 'Interpreter', 'latex', 'FontSize', 14)
zlabel('Total Cell Count', 'Interpreter', 'latex', 'FontSize', 14)

lgd2 = legend(h_tc, {['Low: ' num2str(tc_low)], ['High: ' num2str(tc_high)]}, 'Location', 'best');
set(lgd2, 'Interpreter', 'latex', 'FontSize', 11)


%% Prolif Cells
figure(3)


z_PR = max(PR_Storage(:)) + 10;


plot3(zp_pr_L, gamma_pr_L, z_PR * ones(size(zp_pr_L)), 'm--', 'LineWidth', 2.5);
hold on
plot3(zp_pr_H, gamma_pr_H, z_PR * ones(size(zp_pr_H)), 'm-', 'LineWidth', 2.5);

surf(ZP_mesh, log10(GAMMA_mesh), PR_Storage)

h_pr(3) = plot3(zp_tc_L, gamma_tc_L, z_PR * ones(size(zp_tc_L)), 'w--', 'LineWidth', 2.5);
h_pr(4) = plot3(zp_tc_H, gamma_tc_H, z_PR * ones(size(zp_tc_H)), 'w-', 'LineWidth', 2.5);

%h_pr(3) = [];
%h_pr(4) = [];
shading interp; colormap turbo; colorbar; grid on; grid minor;
view(2)


set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)
xlabel('Upper Proliferation Boundary ($z_p$)', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\log_{10}(\gamma)$', 'Interpreter', 'latex', 'FontSize', 14)
zlabel('Proliferating Cells', 'Interpreter', 'latex', 'FontSize', 14)

lgd3 = legend({...
    ['Prolif Low: ' num2str(pr_low)], ...
    ['Prolif High: ' num2str(pr_high)], ...
   
}, 'Location', 'best');
set(lgd3, 'Interpreter', 'latex', 'FontSize', 11)


%% Crypt Renewal
figure(4)

z_CR = max(CR_Storage(:)) + 1; 


plot3(zp_cr_L, gamma_cr_L, z_CR * ones(size(zp_cr_L)), 'c--', 'LineWidth', 2.5);

hold on
plot3(zp_cr_H, gamma_cr_H, z_CR * ones(size(zp_cr_H)), 'c-', 'LineWidth', 2.5);

surf(ZP_mesh, log10(GAMMA_mesh), CR_Storage)

h_cr(3) = plot3(zp_tc_L, gamma_tc_L, z_CR * ones(size(zp_tc_L)), 'w--', 'LineWidth', 2.5);
h_cr(4) = plot3(zp_tc_H, gamma_tc_H, z_CR * ones(size(zp_tc_H)), 'w-', 'LineWidth', 2.5);

shading interp; colormap turbo; colorbar; grid on; grid minor;
view(2)


set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)
xlabel('Upper Proliferation Boundary ($z_p$)', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\log_{10}(\gamma)$', 'Interpreter', 'latex', 'FontSize', 14)
zlabel('Crypt Renewal Time (Days)', 'Interpreter', 'latex', 'FontSize', 14)

lgd4 = legend( {...
    ['Renewal Low: ' num2str(cr_low) ' Days'], ...
    ['Renewal High: ' num2str(cr_high) ' Days'], ...
  
}, 'Location', 'best');
set(lgd4, 'Interpreter', 'latex', 'FontSize', 11)


%% Combined Plot
% figure(5)
% hold on
% z_comb = max(TC_Storage(:)) + 15; 
% 
% %% 1. Interpolate to a Fine Grid for Seamless Shading
% % Create a high-resolution mesh (200x200) for smooth boundary matching
% zp_fine = linspace(min(zp_vec), max(zp_vec), 200);
% gamma_log_fine = linspace(min(log10(gamma_vec)), max(log10(gamma_vec)), 200);
% [ZP_fine, GAMMA_log_fine] = meshgrid(zp_fine, gamma_log_fine);
% 
% % Interpolate original data onto the fine grid
% TC_fine = interp2(ZP_mesh, log10(GAMMA_mesh), TC_Storage, ZP_fine, GAMMA_log_fine, 'linear');
% PR_fine = interp2(ZP_mesh, log10(GAMMA_mesh), PR_Storage, ZP_fine, GAMMA_log_fine, 'linear');
% CR_fine = interp2(ZP_mesh, log10(GAMMA_mesh), CR_Storage, ZP_fine, GAMMA_log_fine, 'linear');
% 
% % Calculate the logical mask on the high-resolution grid
% feasible_region_fine = (TC_fine >= tc_low) & (TC_fine <= tc_high) & ...
%                        (PR_fine >= pr_low) & (PR_fine <= pr_high) & ...
%                        (CR_fine >= cr_low) & (CR_fine <= cr_high);
% 
% % Draw the smooth shaded region
% if any(feasible_region_fine(:))
%     [~, h_shade] = contourf(ZP_fine, GAMMA_log_fine, double(feasible_region_fine), [0.5 0.5]);
% 
%     % Soft green shading with no border
%     set(h_shade, 'FaceColor', [0.4660, 0.6740, 0.1880], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
% 
%     % Push the shading patch to the same Z-plane as the lines
%     drawnow;
%     try
%         h_shade.ContourParameters.ZData = z_comb * ones(size(h_shade.ContourParameters.XData));
%     catch
%     end
% end
% 
% %% 2. Plot the 6 Boundary Lines
% h_all(1) = plot3(zp_tc_L, gamma_tc_L, z_comb*ones(size(zp_tc_L)), 'k--', 'LineWidth', 3);
% h_all(2) = plot3(zp_tc_H, gamma_tc_H, z_comb*ones(size(zp_tc_H)), 'k-', 'LineWidth', 3);
% 
% h_all(3) = plot3(zp_pr_L, gamma_pr_L, z_comb*ones(size(zp_pr_L)), 'm--', 'LineWidth', 3);
% h_all(4) = plot3(zp_pr_H, gamma_pr_H, z_comb*ones(size(zp_pr_H)), 'm-', 'LineWidth', 3);
% 
% h_all(5) = plot3(zp_cr_L, gamma_cr_L, z_comb*ones(size(zp_cr_L)), 'c--', 'LineWidth', 3);
% h_all(6) = plot3(zp_cr_H, gamma_cr_H, z_comb*ones(size(zp_cr_H)), 'c-', 'LineWidth', 3);
% 
% grid on; grid minor;
% view(2) 
% 
% % LaTeX Styling
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)
% xlabel('Upper Proliferation Boundary ($z_p$)', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('$\log_{10}(\gamma)$', 'Interpreter', 'latex', 'FontSize', 14)
% 
% % Generate legend labels with math mode
% labels = {...
%     ['Total Cells (Low: ' num2str(tc_low) ')'], ...
%     ['Total Cells (High: ' num2str(tc_high) ')'], ...
%     ['Prolif Cells (Low: ' num2str(pr_low) ')'], ...
%     ['Prolif Cells (High: ' num2str(pr_high) ')'], ...
%     ['Renewal (Low: ' num2str(cr_low) ' Days)'], ...
%     ['Renewal (High: ' num2str(cr_high) ' Days)']...
% };
% 
% if any(feasible_region_fine(:))
%     lgd5 = legend([h_all, h_shade], [labels, {'Feasible Region (All Targets Met)'}], 'Location', 'northeastoutside');
% else
%     lgd5 = legend(h_all, labels, 'Location', 'northeast');
% end
% set(lgd5, 'Interpreter', 'latex', 'FontSize', 11)
%% Combined Plot (Figure 5)
figure(5)
hold on
z_comb = max(TC_Storage(:)) + 15; 

%% 1. High-Resolution Interpolation & Feasible Region Shading
zp_fine = linspace(min(zp_vec), max(zp_vec), 200);
gamma_log_fine = linspace(min(log10(gamma_vec)), max(log10(gamma_vec)), 200);
[ZP_fine, GAMMA_log_fine] = meshgrid(zp_fine, gamma_log_fine);

TC_fine = interp2(ZP_mesh, log10(GAMMA_mesh), TC_Storage, ZP_fine, GAMMA_log_fine, 'spline');
PR_fine = interp2(ZP_mesh, log10(GAMMA_mesh), PR_Storage, ZP_fine, GAMMA_log_fine, 'spline');
CR_fine = interp2(ZP_mesh, log10(GAMMA_mesh), CR_Storage, ZP_fine, GAMMA_log_fine, 'spline');

feasible_region_fine = (TC_fine <= tc_high) & ...
    (PR_fine >= pr_low) & (PR_fine <= pr_high) & ...
    (CR_fine >= cr_low) & (CR_fine <= cr_high);

h_shade = [];
if any(feasible_region_fine(:))
    [~, h_shade] = contourf(ZP_fine, GAMMA_log_fine, double(feasible_region_fine), [0.5 0.5]);
    set(h_shade, 'FaceColor', [0.4660, 0.6740, 0.1880], 'FaceAlpha', 0.35, 'EdgeColor', 'none');
end

%% 2. Plot Boundary Lines Dynamic Legend Setup
legend_handles = [];
legend_labels = {};

% Helper inline function to safely add lines and legend entries
add_line = @(zp_c, g_c, style, color, label) ...
    deal(plot3(zp_c, g_c, z_comb*ones(size(zp_c)), style, 'Color', color, 'LineWidth', 2.5), label);

if ~isempty(zp_tc_L)
    [h, l] = add_line(zp_tc_L, gamma_tc_L, '--', 'k', ['TC Low (' num2str(tc_low) ')']);
    legend_handles = [legend_handles, h]; legend_labels{end+1} = l;
end
if ~isempty(zp_tc_H)
    [h, l] = add_line(zp_tc_H, gamma_tc_H, '-', 'k', ['TC High (' num2str(tc_high) ')']);
    legend_handles = [legend_handles, h]; legend_labels{end+1} = l;
end
if ~isempty(zp_pr_L)
    [h, l] = add_line(zp_pr_L, gamma_pr_L, '--', 'm', ['Prolif Low (' num2str(pr_low) ')']);
    legend_handles = [legend_handles, h]; legend_labels{end+1} = l;
end
if ~isempty(zp_pr_H)
    [h, l] = add_line(zp_pr_H, gamma_pr_H, '-', 'm', ['Prolif High (' num2str(pr_high) ')']);
    legend_handles = [legend_handles, h]; legend_labels{end+1} = l;
end
if ~isempty(zp_cr_L)
    [h, l] = add_line(zp_cr_L, gamma_cr_L, '--', 'c', ['Renewal Low (' num2str(cr_low) ')']);
    legend_handles = [legend_handles, h]; legend_labels{end+1} = l;
end
if ~isempty(zp_cr_H)
    [h, l] = add_line(zp_cr_H, gamma_cr_H, '-', 'c', ['Renewal High (' num2str(cr_high) ')']);
    legend_handles = [legend_handles, h]; legend_labels{end+1} = l;
end

% Append shaded region handle if present
if ~isempty(h_shade)
    legend_handles = [legend_handles, h_shade];
    legend_labels{end+1} = 'Feasible Region';
end

%% 3. Apply Legend and Formatting
grid on; grid minor; view(2);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
xlabel('Upper Proliferation Boundary ($z_p$)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\log_{10}(\gamma)$', 'Interpreter', 'latex', 'FontSize', 14);

if ~isempty(legend_handles)
    lgd = legend(legend_handles, legend_labels, 'Location', 'northeastoutside');
    set(lgd, 'Interpreter', 'latex', 'FontSize', 11);
end