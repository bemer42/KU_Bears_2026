function solve_crypt_pde_system()
    % Spatial domain z in [0, L] and time domain t
    L = 1.0;
    z = linspace(0, L, 100);
    t = linspace(0, 10, 101);

    % Solve PDE system: [q; s; r]
    % m = 0 for slab (1D Cartesian) geometry
    m = 0; 
    sol = pdepe(m, @pdefun, @icfun, @bcfun, z, t);

    % Extract variables: sol(time_idx, space_idx, pde_idx)
    q = sol(:,:,1);
    s = sol(:,:,2);
    r = sol(:,:,3);

    % Plot final state
    figure;
    subplot(3,1,1); plot(z, q(end,:), 'b', 'LineWidth', 1.5); ylabel('q'); title('Solution at t_{final}'); grid on;
    subplot(3,1,2); plot(z, s(end,:), 'r', 'LineWidth', 1.5); ylabel('s'); grid on;
    subplot(3,1,3); plot(z, r(end,:), 'g', 'LineWidth', 1.5); ylabel('r'); xlabel('z'); grid on;
end

% PDE Definition: c(z,t,u,dudz) * du/dt = x^(-m) * d/dz(x^m * f(z,t,u,dudz)) + s(z,t,u,dudz)
function [c, f, s_term] = pdefun(z, t, u, dudz)
    q = u(1); s = u(2); r = u(3);
    dqdz = dudz(1); dsdz = dudz(2); drdz = dudz(3);

    alpha_z = 0.1;
    Ds = 0.05;
    Dr = 0.01;
    zp = 0.5;
    Tc = 1.0;
    gamma_uptake = 0.2;

    % Capacity (d u / d t coefficients)
    c = [1; 1; 1];

    % Flux terms: f = c_coeff * du/dz
    % Note: non-linear q diffusion flux uses alpha_z / ((1 + drdz^2) * q^2)
    f = [ (alpha_z / ((1 + drdz^2) * (q^2 + 1e-6))) * dqdz;
          Ds * dsdz;
          Dr * drdz ];

    % Source terms
    H = double(z <= zp);
    dt_r = 0; % Can be coupled or updated iteratively
    
    fq = (log(2)/Tc) * H * q - ((drdz * dt_r) / (1 + drdz^2)) * q;
    fs = 0.1 * exp(-z) - gamma_uptake * s * q;
    fr = 0.02 * (q - 1.0);

    s_term = [fq; fs; fr];
end

% Initial Conditions
function u0 = icfun(z)
    u0 = [1.0; 0.5; 0.2 + 0.05*z];
end

% Boundary Conditions: [pl, ql, pr, qr] such that p + q*f = 0
function [pl, ql, pr, qr] = bcfun(xl, ul, xr, ur, t)
    % Left boundary (z = 0): Zero flux -> f = 0
    pl = [0; 0; 0];
    ql = [1; 1; 1];

    % Right boundary (z = L): Dirichlet q=1, s=0.5, r=0.25 -> u - u_bound = 0
    pr = [ur(1) - 1.0; ur(2) - 0.5; ur(3) - 0.25];
    qr = [0; 0; 0];
end