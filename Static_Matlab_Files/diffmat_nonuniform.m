function Dz = diffmat_nonuniform(z)
%FIRSTDERIVMATRIXNONUNIFORM  Sparse 1st-derivative matrix on nonuniform grid
%   fz ≈ Dz * f   for f as Nz×1

z = z(:);
Nz = numel(z);

Dz = spalloc(Nz, Nz, 3*Nz);

% ---- interior nodes: i = 2..Nz-1 (3-point, 2nd-order, nonuniform) ----
for i = 2:Nz-1
    hm = z(i)   - z(i-1);   % h_{i-1}
    hp = z(i+1) - z(i);     % h_i

    % weights for f'(z_i) using f_{i-1}, f_i, f_{i+1}
    w_im1 = -hp/(hm*(hm+hp));
    w_i   = (hp - hm)/(hm*hp);
    w_ip1 =  hm/(hp*(hm+hp));

    Dz(i,i-1) = w_im1;
    Dz(i,i)   = w_i;
    Dz(i,i+1) = w_ip1;
end

% ---- left boundary: i = 1 (one-sided, 2nd-order, nonuniform using z1,z2,z3) ----
h1 = z(2) - z(1);
h2 = z(3) - z(2);
% weights derived from Lagrange polynomial through (z1,z2,z3), evaluated at z1
Dz(1,1) = -(2*h1 + h2)/(h1*(h1+h2));
Dz(1,2) =  (h1 + h2)/(h1*h2);
Dz(1,3) = -h1/(h2*(h1+h2));

% ---- right boundary: i = Nz (one-sided, 2nd-order, nonuniform using z_{n-2},z_{n-1},z_n) ----
h1 = z(Nz)   - z(Nz-1);
h2 = z(Nz-1) - z(Nz-2);
% mirror of left boundary, evaluated at z_n
Dz(Nz,Nz)   =  (2*h1 + h2)/(h1*(h1+h2));
Dz(Nz,Nz-1) = -(h1 + h2)/(h1*h2);
Dz(Nz,Nz-2) =  h1/(h2*(h1+h2));
end