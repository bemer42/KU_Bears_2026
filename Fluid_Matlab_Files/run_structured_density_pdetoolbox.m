%% Stemness using PDE solver

function run_structured_density_pdetoolbox()
    % 1. Create a 2D PDE Model container
    model = createpde();
    
    % 2. Define the Rectangular Domain: (sigma, z) in [0,1] x [0,1]
    R1 = [3, 4, 0, 1, 1, 0, 0, 0, 1, 1]';
    g = decsg(R1);
    geometryFromEdges(model, g);
    
    % 3. Specify Boundary Conditions (No-flux)
    applyBoundaryCondition(model, 'neumann', 'Edge', 1:4, 'g', 0, 'q', 0);
    
    % 4. Define Coefficients
    d_coeff = 1;
    c_coeff = [0.01; 0.01]; % Diffusion: [alpha_sigma; alpha_z]
    
    a_coeff = @(location, state) reactionCoefficient(location, state);
    f_coeff = @(location, state) advectionSourceTerm(location, state);
    
    specifyCoefficients(model, 'm', 0, 'd', d_coeff, 'c', c_coeff, ...
                               'a', a_coeff, 'f', f_coeff);
                           
    % 5. Set Initial Condition p(sigma, z, t=0)
    setInitialConditions(model, @(location) initialDensity(location));
    
    % 6. Generate Mesh & Solve over Time
    generateMesh(model, 'Hmax', 0.03); 
    
    tlist = linspace(0, 5, 21);
    results = solvepde(model, tlist);
    
    % 7. Plot Selected Time Snapshots
    figure;
    snapshot_indices = [1, 6, 11, 21];
    for idx = 1:length(snapshot_indices)
        k = snapshot_indices(idx);
        subplot(2, 2, idx);
        pdeplot(model, 'XYData', results.NodalSolution(:, k));
        title(sprintf('t = %.2f', tlist(k)));
        xlabel('Stemness (\sigma)');
        ylabel('Crypt Position (z)');
        colorbar;
    end
end

% --- Fixed Helper Functions ---

function a = reactionCoefficient(location, ~)
    % Reaction term: a = -Gamma(sigma, z)
    sigma = location.x;
    z = location.y;
    gamma = 0.8 * sigma .* (1 - z);
    a = -gamma; 
end

function f = advectionSourceTerm(location, state)
    % Advection term: - (v_sigma * dp/d_sigma + v_z * dp/dz)
    v_sigma = 0.3; % Stemness degradation speed
    v_z = 0.8;     % Upward spatial speed
    
    % Safe check for lowercase gradient fields (ux, uy)
    if nargin < 2 || ~isstruct(state) || ~isfield(state, 'ux') || isnan(state.ux(1))
        f = zeros(size(location.x));
    else
        % state.ux is dp/d_sigma, state.uy is dp/dz
        f = -(v_sigma * state.ux + v_z * state.uy);
    end
end

function p0 = initialDensity(location)
    sigma = location.x;
    z = location.y;
    p0 = exp(-((sigma - 1).^2 + (z - 0).^2) / 0.02);
end