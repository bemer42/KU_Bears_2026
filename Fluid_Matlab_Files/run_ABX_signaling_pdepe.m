function run_ABX_signaling_pdepe()
clear; close all; clc

% Discretize Space
L  = 78.8783; 
Nz = 100;
z  = linspace(0, L, Nz);

% Discretize time
Nt = 50;
t  = linspace(0, 100, Nt); 

% Solve 1D PDE System 
m = 0; 
sol = pdepe(m, @pdefun, @icfun, @bcfun, z, t);

% Extract Steady-State Profiles 
A_ss = sol(end, :, 1); % APC
B_ss = sol(end, :, 2); % Beta-catenin
X_ss = sol(end, :, 3); % Axin

% Plot 
figure('Color', 'w');
plot(z, A_ss, 'r-',  'LineWidth', 2.5, 'DisplayName', 'APC ($A$)'); hold on;
plot(z, B_ss, 'b-',  'LineWidth', 2.5, 'DisplayName', '$\beta$-catenin ($B$)');
plot(z, X_ss, 'g--', 'LineWidth', 2.5, 'DisplayName', 'Axin ($X$)');

set(gca, 'FontSize', 14);
title('Multi-Cell ABX Signaling Steady-State Profiles', 'FontSize', 18, 'Interpreter', 'latex');
xlabel('Crypt Axial Position $z$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Concentration', 'FontSize', 14, 'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 12);
grid on; grid minor; box on;
xlim([0 L]);
end

% --- Helper Functions for pdepe ---

function [c, f, s] = pdefun(z, t, u, Dudz)

A = u(1);
B = u(2);
X = u(3);

% Diffusion coefficients [d_a; d_b; d_x]
da = 0.5;
db = 1.0;
dx = 0.2;
c = [1; 1; 1]; % c*du/dt = d/dz(f) + s

% Flux terms (f = d * du/dz)
f = [da; db; dx] .* Dudz;

% Wnt
W = exp(-z / 10);

% Reaction dynamics f_a, f_b, f_x 
v0 = 0.1; 
v_deg = 0.5;

s1 = v0 - v_deg * A * B;             % APC dynamics
s2 = W * 2.0 - (v_deg / (1 + X)) * B; % Beta-catenin dynamics (Wnt upregulates)
s3 = 0.05 - 0.1 * X + 0.2 * A;       % Axin dynamics

s = [s1; s2; s3];
end

function u0 = icfun(z)
% Initial conditions at t = 0 (uniform initial concentrations)
u0 = [1.0; 0.5; 0.2];
end

function [pl, ql, pr, qr] = bcfun(xl, ul, xr, ur, t)
% Zero-flux Neumann boundary conditions
pl = [0; 0; 0];
ql = [1; 1; 1]; % Left boundary (z = 0)

pr = [0; 0; 0];
qr = [1; 1; 1]; % Right boundary (z = L)
end