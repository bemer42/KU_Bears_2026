function dpdt = dpdt_1D_snipsnap(t, p_int, S, Z, rz, vs, Ds, Dz, par, bcType)

% Define limiter for upwind scheme:
if isfield(par,'limiter')
    limiterType = string(par.limiter);
else
    limiterType = "vanleer"; % "vanleer", "mc", "minmod"
end

% Collect parameters:
alpha_s = par.alpha_s;
alpha_z = par.alpha_z;
k       = par.k;
s_stem  = par.s_stem;
s_ta    = par.s_ta;

% Collect dimension sizes:
Ns = size(Ds,1);
Nz = size(Dz,1);

% Build full P from interior
P_int = reshape(p_int, Ns-2, Nz-2);
P = zeros(Ns, Nz);
P(2:Ns-1, 2:Nz-1) = P_int;

% Uniform grid spacing inferred from mesh
svec = S(:,1);   ds = svec(2) - svec(1);
zvec = Z(1,:).'; dz = zvec(2) - zvec(1);

% trapezoid weights for q(z)=∫ p ds
wts = ones(Ns,1); wts(1)=0.5; wts(end)=0.5;
qcol = @(pcol) ds * (wts.' * pcol);

% z one-sided stencil coeffs (for q_z=0 at boundaries)
b0 = Dz(1,1);       b1 = Dz(1,2);       b2 = Dz(1,3);
t2 = Dz(end,end-2); t1 = Dz(end,end-1); t0 = Dz(end,end);

% s boundary stencil coeffs (for Robin no-flux closure)
a0L = Ds(1,1);       a1L = Ds(1,2);       a2L = Ds(1,3);
a0R = Ds(end,end);   a1R = Ds(end,end-1); a2R = Ds(end,end-2);

% ------------------------------------------------------------
% 1) Enforce s no-flux boundary values on interior z-columns
% ------------------------------------------------------------
P = applyNoFlux_s(P, vs, 2:Nz-1);

% ------------------------------------------------------------
% 2) Apply z boundary conditions (q-based)
% ------------------------------------------------------------
switch bcType
    case "NbDt"
        % bottom: q_z(0)=0 => b0*q1 + b1*q2 + b2*q3 = 0
        q2 = qcol(P(:,2));
        q3 = qcol(P(:,3));
        q1 = -(b1*q2 + b2*q3)/b0;
        P(:,1) = scaleColumnToIntegral(P(:,2), q2, q1);

        % top: q(Lz)=1
        qNm1 = qcol(P(:,end-1));
        P(:,end) = P(:,end-1)/qNm1;

    case "NbNt"
        q2 = qcol(P(:,2));
        q3 = qcol(P(:,3));
        q1 = -(b1*q2 + b2*q3)/b0;
        P(:,1) = scaleColumnToIntegral(P(:,2), q2, q1);

        qNm1 = qcol(P(:,end-1));
        qNm2 = qcol(P(:,end-2));
        qN   = -(t1*qNm1 + t2*qNm2)/t0;
        P(:,end) = scaleColumnToIntegral(P(:,end-1), qNm1, qN);

    case "DbNt"
        q2 = qcol(P(:,2));
        P(:,1) = P(:,2)/q2;

        qNm1 = qcol(P(:,end-1));
        qNm2 = qcol(P(:,end-2));
        qN   = -(t1*qNm1 + t2*qNm2)/t0;
        P(:,end) = scaleColumnToIntegral(P(:,end-1), qNm1, qN);

    case "DbDt"
        q2 = qcol(P(:,2));
        P(:,1) = P(:,2)/q2;

        qNm1 = qcol(P(:,end-1));
        P(:,end) = P(:,end-1)/qNm1;

    otherwise
        error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
end

% enforce s no-flux also on z-boundary columns
P = applyNoFlux_s(P, vs, [1 Nz]);

% ------------------------------------------------------------
% 3) Build q(z), q_z(z), curvature factor g, and z-velocity
% ------------------------------------------------------------
qz = ds * (wts.' * P);               
if any(qz <= 0)
    error('Encountered nonpositive q(z,t); cannot form alpha_z*q_z/q^3 safely.');
end

qz_z = Dz * qz(:);                     
u    = alpha_z * (qz_z.' ./ (qz.^3));  

g    = sqrt(1 + rz.^2);                 % Ns x Nz
invG = 1 ./ g;                          % Ns x Nz

% w(z) = (1/g) * (alpha_z/q^3 * q_z)
Wz = invG .* repmat(u, Ns, 1);          % Ns x Nz

% ------------------------------------------------------------
% 4) Divergences
% ------------------------------------------------------------
div_s = div_s_MUSCL(P, vs, alpha_s, ds, limiterType);

% z-term: (1/g) * d/dz( (w(z)*p) )
dFdz  = div_z_MUSCL(P, Wz, dz, limiterType);  % returns d/dz( Wz * P )
div_z = invG .* dFdz;

% proliferation: k*p*I(zb<z<zu)
prolifMask = (S < s_stem) & (S > s_ta);
Source = k * P .* prolifMask;

dPdt_full = div_s + div_z + Source;

% return interior
dpdt = dPdt_full(2:Ns-1, 2:Nz-1);
dpdt = dpdt(:);

% ============================================================
% helpers
% ============================================================
    function Pout = applyNoFlux_s(Pin, VSin, zCols)
        % Enforce -vs*p + alpha_s*p_s = 0 at s boundaries, columnwise.
        Pout = Pin;

        vsL = VSin(1,   zCols);
        vsR = VSin(end, zCols);

        if alpha_s == 0
            maskL = abs(vsL) > 0;
            maskR = abs(vsR) > 0;
            Pout(1,   zCols(maskL)) = 0;
            Pout(end, zCols(maskR)) = 0;
            Pout(1,   zCols(~maskL)) = Pout(2,     zCols(~maskL));
            Pout(end, zCols(~maskR)) = Pout(end-1, zCols(~maskR));
            return
        end

        denomL = (-vsL + alpha_s*a0L);
        denomR = (-vsR + alpha_s*a0R);

        if any(abs(denomL) < 1e-14) || any(abs(denomR) < 1e-14)
            error('s no-flux enforcement ill-conditioned (check vs, alpha_s, Ds).');
        end

        Pout(1, zCols)   = -alpha_s*(a1L*Pout(2,zCols)     + a2L*Pout(3,zCols))     ./ denomL;
        Pout(end, zCols) = -alpha_s*(a1R*Pout(end-1,zCols) + a2R*Pout(end-2,zCols)) ./ denomR;
    end

    function pScaled = scaleColumnToIntegral(pRef, qRef, qTarget)
        if qRef <= 0 || qTarget <= 0
            error('scaleColumnToIntegral: qRef and qTarget must be positive.');
        end
        pScaled = (qTarget/qRef) * pRef;
    end

    function divs = div_s_MUSCL(Pin, VSin, alpha_s, ds, limiterType)
        % d/ds( -vs*p + alpha_s*p_s ) with MUSCL-TVD for advection
        eps0 = 1e-14;

        % slopes in s
        dP = Pin(2:Ns,:) - Pin(1:Ns-1,:);               % (Ns-1) x Nz
        slope = zeros(Ns, Nz);

        r = dP(1:Ns-2,:) ./ (dP(2:Ns-1,:) + eps0);      % (Ns-2) x Nz for i=2..Ns-1
        phi = limiter(r, limiterType);
        slope(2:Ns-1,:) = phi .* dP(2:Ns-1,:);

        % face reconstructions
        PL = Pin(1:Ns-1,:) + 0.5*slope(1:Ns-1,:);
        PR = Pin(2:Ns,:)   - 0.5*slope(2:Ns,:);

        vface = 0.5*(VSin(1:Ns-1,:) + VSin(2:Ns,:));
        Pup   = (vface >= 0).*PL + (vface < 0).*PR;

        Fadv  = -vface .* Pup;
        Fdiff = alpha_s * (Pin(2:Ns,:) - Pin(1:Ns-1,:)) / ds;

        F = Fadv + Fdiff;                               % (Ns-1) x Nz

        divs = zeros(Ns, Nz);
        divs(2:Ns-1,:) = (F(2:Ns-1,:) - F(1:Ns-2,:)) / ds;
    end

    function dFdz = div_z_MUSCL(Pin, VZin, dz, limiterType)
        % returns d/dz( Vz * p ) using MUSCL-TVD upwind across z
        eps0 = 1e-14;

        dPz = Pin(:,2:Nz) - Pin(:,1:Nz-1);             % Ns x (Nz-1)
        slope = zeros(Ns, Nz);

        r = dPz(:,1:Nz-2) ./ (dPz(:,2:Nz-1) + eps0);    % Ns x (Nz-2) for j=2..Nz-1
        phi = limiter(r, limiterType);
        slope(:,2:Nz-1) = phi .* dPz(:,2:Nz-1);

        PL = Pin(:,1:Nz-1) + 0.5*slope(:,1:Nz-1);
        PR = Pin(:,2:Nz)   - 0.5*slope(:,2:Nz);

        vface = 0.5*(VZin(:,1:Nz-1) + VZin(:,2:Nz));
        Pup   = (vface >= 0).*PL + (vface < 0).*PR;

        F = vface .* Pup;                               % Ns x (Nz-1)

        dFdz = zeros(Ns, Nz);
        dFdz(:,2:Nz-1) = (F(:,2:Nz-1) - F(:,1:Nz-2)) / dz;
    end

    function phi = limiter(r, limiterType)
        switch lower(char(limiterType))
            case 'minmod'
                phi = max(0, min(1, r));
            case 'mc'
                phi = max(0, min(min(2*r, (1+r)/2), 2));
            case 'vanleer'
                phi = (r + abs(r)) ./ (1 + abs(r));
            otherwise
                error('Unknown limiterType. Use "vanleer", "mc", or "minmod".');
        end
    end
end