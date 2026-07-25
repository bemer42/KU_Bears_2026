function run_2D_cell_movement_pdetoolbox()
% 1. Create a 2D PDE Model container
model = createpde();

% 2. Define Rectangular Domain: (theta, z_hat) in [0, 2*pi] x [0, 1]
% Geometry format: [3 (rect), 4 (sides), xmin, xmax, xmax, xmin, ymin, ymin, ymax, ymax]
R1 = [3, 4, 0, 2*pi, 2*pi, 0, 0, 0, 1, 1]';
g = decsg(R1);
geometryFromEdges(model, g);

% 3. Specify Boundary Conditions
% Edge 1: Bottom boundary (z_hat = 0) -> Zero flux (Neumann)
applyBoundaryCondition(model, 'neumann', 'Edge', 1, 'g', 0, 'q', 0);

% Edge 3: Top boundary (z_hat = 1) -> Fixed density (Dirichlet q_hat = 1)
applyBoundaryCondition(model, 'dirichlet', 'Edge', 3, 'u', 1);

% Edges 2 and 4: Periodic in Theta (theta = 0 and theta = 2*pi)
% PDE Toolbox handles periodicity via identifying coupled edges
applyBoundaryCondition(model, 'neumann', 'Edge', [2, 4], 'g', 0, 'q', 0);

% 4. Model Parameters
gamma_z  = 1.2;
gamma_th = 0.12;
Zp_hat   = 0.34; % Proliferative height (zp / L)

% 5. Define Coefficients: d*dq/dt - div(c*grad(q)) + a*q = f
d_coeff = gamma_z;

% Function handles for non-linear diffusion matrix 'c' and reaction term 'a'
c_coeff = @(location, state) diffusionTensor(location, state, gamma_z, gamma_th);
a_coeff = @(location, state) proliferationTerm(location, state, Zp_hat);
f_coeff = 0;

specifyCoefficients(model, 'm', 0, 'd', d_coeff, 'c', c_coeff, ...
    'a', a_coeff, 'f', f_coeff);

% 6. Initial Condition: Uniform density q_hat = 1
setInitialConditions(model, 1);

% 7. Generate Fine Mesh & Solve over Time
generateMesh(model, 'Hmax', 0.05);

tlist = linspace(0, 10, 21); % Nondimensional time
results = solvepde(model, tlist);

% 8. Plot Steady-State Density Surface
q_ss = results.NodalSolution(:, end);

figure('Color', 'w');
pdeplot(model, 'XYData', q_ss, 'Mesh', 'off');
title('Steady State Cell Density $\hat{q}(\theta, \hat{z})$', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\hat{z}$', 'Interpreter', 'latex', 'FontSize', 14);
colormap turbo;
colorbar;
end

% --- Helper Functions ---

function c = diffusionTensor(location, state, gamma_z, gamma_th)
    % Ensure location coordinates are row vectors (1 x N)
    x = location.x(:)';
    z_hat = location.y(:)';
    N = length(z_hat);
    
    % Radius profile R_hat(z_hat)
    r_b = 41/(2*pi);
    r_t = 10/pi;
    a   = 0.3;
    L   = 78.8783;
    R_hat = (1 - exp(-a*L*z_hat)) + (r_t/r_b)*exp(a*L*(z_hat - 1));
    
    % Extract current cell density q_hat safely as a 1 x N row vector
    if nargin < 2 || ~isstruct(state) || ~isfield(state, 'u') || isempty(state.u) || any(isnan(state.u))
        q_hat = ones(1, N);
    else
        q_hat = state.u(:)'; % Force row vector (1 x N)
    end
    
    % Compute tensor components explicitly as 1 x N row vectors
    c11 = gamma_z ./ (gamma_th * (R_hat.^2));
    c12 = zeros(1, N);
    c21 = zeros(1, N);
    c22 = 1 ./ (q_hat.^2);
    
    % Vertically concatenate into a 4 x N matrix for PDE Toolbox
    c = [c11; c12; c21; c22];
end

function a = proliferationTerm(location, ~, Zp_hat)
    % Ensure z_hat is a row vector
    z_hat = location.y(:)';
    
    % Smooth Heaviside step function around Zp_hat
    H = 0.5 * (1 + tanh(50 * (Zp_hat - z_hat)));
    
    % Reaction coefficient: a = -H
    a = -H; 
end