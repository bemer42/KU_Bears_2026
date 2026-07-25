function run_1D_nondimensional_cell_movement_pdepe()
clear; close all; clc

% 1. Discretize Spatial Axis z_hat in [0, 1]
Nz = 150;
z_hat = linspace(0, 1, Nz);

% 2. Time Mesh to Integration Steady State
Nt = 100;
t_hat = linspace(0, 10, Nt);

% 3. Model Parameters
gamma_z = 7;
Zp_hat  = 0.34; % Proliferative height (zp / L)

% 4. Solve using pdepe (m = 0 for 1D Cartesian)
m = 0;
sol = pdepe(m, @(z, t, q, dqdz) pdefun(z, t, q, dqdz, gamma_z, Zp_hat), ...
    @icfun, ...
    @bcfun, ...
    z_hat, t_hat);

% Extract Steady-State Cell Density profile (final time step)
q_ss = sol(end, :);

% 5. Plot Steady State Cell Density vs z_hat
figure('Color', 'w');
plot(z_hat, q_ss, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Cell Density');
hold on;
xline(Zp_hat, 'r--', 'LineWidth', 2, 'DisplayName', 'Proliferative Boundary ($z_p$)');

set(gca, 'FontSize', 14);
title('1D Nondimensional Steady-State Cell Density $\hat{q}(\hat{z})$', 'FontSize', 18, 'Interpreter', 'latex');
xlabel('Crypt Axis $\hat{z}$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Cell Density $\hat{q}$', 'FontSize', 14, 'Interpreter', 'latex');
grid on; grid minor; box on;
xlim([0 1]);
legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12);
end

% --- Helper Functions for pdepe ---

function [c, f, s] = pdefun(z, ~, q, dqdz, gamma_z, Zp_hat)
% Geometry metric function g_hat(z)
r_b = 41 / (2*pi);
r_t = 10 / pi;
a   = 0.3;
L   = 78.8783;

R_hat  = (1 - exp(-a*L*z)) + (r_t/r_b)*exp(a*L*(z - 1));
Rz_hat = a*L*exp(-a*L*z) + (r_t/r_b)*a*L*exp(a*L*(z - 1));
epsilon = r_b / L;

% Arc length element metric g_hat(z)
g_hat = sqrt(1 + (epsilon * Rz_hat)^2);

% Smooth Heaviside proliferation term H(Zp_hat - z)
H = 0.5 * (1 + tanh(50 * (Zp_hat - z)));

% c * dq/dt = d/dz(f) + s
c = gamma_z * g_hat;
f = (1 / (g_hat * q^2)) * dqdz; % Non-linear crowding diffusion flux
s = g_hat * H * q;             % Proliferation growth source term
end

function q0 = icfun(~)
% Initial uniform cell density q_hat = 1
q0 = 1.0;
end

function [pl, ql, pr, qr] = bcfun(~, ~, ~, ur, ~)
% Base (z_hat = 0): Zero flux -> f = 0 (Neumann)
pl = 0;
ql = 1;

% Top (z_hat = 1): Fixed density q_hat = 1 -> u - 1 = 0 (Dirichlet)
pr = ur - 1.0;
qr = 0;
end