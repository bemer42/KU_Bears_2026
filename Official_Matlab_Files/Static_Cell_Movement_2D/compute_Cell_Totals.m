function [C, P] = compute_Cell_Totals(Q, S, th, z, maskP, scaleFactor)
% compute_CP_thetaSum_zTrapz
% Computes total cells C and proliferative cells P by:
%   - theta integral: periodic trapezoid on uniform grid -> dth * sum(...)
%   - z integral: composite trapezoid using trapz(z, ...)
%
% Inputs:
%   Q          (Nth x Nz) cell density on (theta,z)
%   S          (Nth x Nz) surface/geometry factor to integrate against
%              (dimensional: G; nondimensional: S_hat)
%   th         (Nth x 1) theta grid, assumed uniform and periodic (no 2*pi endpoint)
%   z          (Nz  x 1) z grid (can be uniform or nonuniform)
%   maskP      (Nth x Nz) logical/double mask for proliferative region (1 inside, 0 outside)
%              (e.g. Z < Zp or Z_hat < Zp_hat)
%   scaleFactor (scalar) optional multiplier (default = 1)
%              (dimensional: 1; nondimensional: r_b*L/ell to convert to physical count scale)
%
% Outputs:
%   Cval, Pval scalars (not rounded)

    if nargin < 6 || isempty(scaleFactor)
        scaleFactor = 1;
    end

    % Basic checks
    if ~isequal(size(Q), size(S))
        error('Q and S must have the same size.');
    end
    if ~isequal(size(Q), size(maskP))
        error('maskP must have the same size as Q.');
    end
    if size(Q,1) ~= numel(th)
        error('Length of th must match size(Q,1).');
    end
    if size(Q,2) ~= numel(z)
        error('Length of z must match size(Q,2).');
    end

    % Theta spacing (uniform)
    dth = th(2) - th(1);

    % Integrate over theta (periodic trapezoid rule on uniform grid)
    thetaIntC = dth * sum(Q .* S, 1);                 % 1 x Nz
    thetaIntP = dth * sum(Q .* S .* maskP, 1);        % 1 x Nz

    % Integrate over z with trapz
    C = scaleFactor * trapz(z, thetaIntC);
    P = scaleFactor * trapz(z, thetaIntP);
end