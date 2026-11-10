function [R, Rth, Rz, G, R_s] = Crypt_Geom_fun(p, split, diffmat, geom)

% Note:
% p is the percentage of the cassini oval that it will max out at, i.e.
% p = 0 is no split at all, p = 1 is the complete spit;
% split is the level at which the splitting is taking place, i.e.
% split = 0 is no split, split = infinity is the full crypt fission,
% split = L/2 is about halfway split.

% Collect parameters:
Z   = geom.Z;
Th  = geom.Th;
r_b = geom.r_b;
r_t = geom.r_t;
a   = geom.a;
L   = geom.L;

% Collect differential matrices:
Dz  = diffmat.Dz;
Dth = diffmat.Dth;

% Define static radius:
R_s   = r_b*(1-exp(-a*Z))+r_t*exp(a*(Z-L));

% Output static crypt if p = 0:
if p == 0
    
    R = R_s;

else

    % Define current crypt structure:
    A = p./(1 + exp(-.2*(split-Z)));
    C = A .* (1 + sqrt(5))/2.*R_s;
    R = sqrt(C.^2.*cos(2*Th)+...
        sqrt((R_s.^2+R_s.*C).^2-C.^4.*sin(2*Th).^2));

end

% Define partial derivatives:
Rz = R*Dz';
Rth = Dth*R;

% Define surface area element:
G = sqrt(Rth.^2+R.^2.*(Rz.^2 + 1));