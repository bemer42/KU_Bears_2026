function dqdt = dqdt_2D_snipsnap(t, q_int, Th, Z, r, rth, rz, Dth, Dz, par, bcType)

% Extract all parameters:
alpha_th = par.alpha_th;
alpha_z  = par.alpha_z;
k        = par.k;
zu       = par.zu;
zb       = par.zb;
thc      = par.thc;
thw      = par.thw;

% Sizes
Nth = size(Dth,1);
Nz  = size(Dz,1);
nz  = Nz - 2;

% Reshape interior state into matrix (Nth x (Nz-2))
Q_int = reshape(q_int, Nth, nz);

% Initialize full matrix and insert interior z-columns
Q = zeros(Nth, Nz);
Q(:,2:Nz-1) = Q_int;

% Convenience: boundary derivative stencil coefficients from Dz
b0 = Dz(1,1);       b1 = Dz(1,2);       b2 = Dz(1,3);          
t2 = Dz(end,end-2); t1 = Dz(end,end-1); t0 = Dz(end,end);  

% Helper that solves Flux_z(:,boundary)=0 for the boundary column qB = Q(:,jB)
% using one-sided stencil: a0*qB + a1*q1 + a2*q2 for q_z at boundary
    function qB = solveFluxNeumannAtBoundary(jB, a0, a1, a2, q1, q2)

        % Coefficients at that boundary column
        Cz  = spdiags(alpha_z*(rth(:,jB).^2 + r(:,jB).^2), 0, Nth, Nth);
        Cth = spdiags(alpha_th*(rth(:,jB).*rz(:,jB)),      0, Nth, Nth);

        % Enforce (numerator of Flux_z)=0:
        %   Cz*(a0*qB + a1*q1 + a2*q2) - Cth*(Dth*qB) = 0
        A = a0*Cz - Cth*Dth;
        b = -Cz*(a1*q1 + a2*q2);
        qB = A \ b;
    end

% Apply boundary conditions
switch bcType
    case "NbDt"
        % bottom Neumann 
        Q(:,1) = solveFluxNeumannAtBoundary(1, b0, b1, b2, Q(:,2), Q(:,3));
        % top Dirichlet
        Q(:,end) = 1;
    case "NbNt"
        % bottom Neumann
        Q(:,1) = solveFluxNeumannAtBoundary(1, b0, b1, b2, Q(:,2), Q(:,3));
        % top Neumann 
        Q(:,end) = solveFluxNeumannAtBoundary(Nz, t0, t1, t2, Q(:,end-1), Q(:,end-2));
    case "DbNt"
        % bottom Dirichlet
        Q(:,1) = 1;
        % top Neumann
        Q(:,end) = solveFluxNeumannAtBoundary(Nz, t0, t1, t2, Q(:,end-1), Q(:,end-2));
    case "DbDt"
        % bottom Dirichlet
        Q(:,1) = 1;
        % top Dirichlet
        Q(:,end) = 1;
    otherwise
        error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
end

% Metric tensor determinant:
det_g = rth.^2 + r.^2 .* (rz.^2 + 1);

% Fluxes:
Flux_th = (1 ./ (Q.^2) ./ sqrt(det_g)) .* ...
          ( alpha_th*(rz.^2+1).*(Dth*Q) - alpha_z*(rth.*rz).*(Q*Dz') );
Flux_z  = (1 ./ (Q.^2) ./ sqrt(det_g)) .* ...
          ( -alpha_th*(rth.*rz).*(Dth*Q) + alpha_z*(rth.^2 + r.^2).*(Q*Dz') );

% Cell proliferation:
P = k * Q .* (Z < zu) .* (Z > zb) .* (abs(Th - thc) < thw);

% PDE:
dQdt_full = (1 ./ sqrt(det_g)) .* (Dth*Flux_th) + ...
            (1 ./ sqrt(det_g)) .* (Flux_z*Dz') + P;

% Return only interior z-columns 
dQdt = dQdt_full(:,2:end-1);
dqdt = dQdt(:);

end