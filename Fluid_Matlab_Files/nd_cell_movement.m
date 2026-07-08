function total_cells = nd_cell_movement(gamma, zp, bcType)


Nz = 100;                 
z = linspace(0, 1, Nz)';  
dz = z(2) - z(1);


geom.z = z;
geom.rz = zeros(size(z)); 


Dz = zeros(Nz, Nz);
for i = 2:Nz-1
    Dz(i, i-1) = -1 / (2*dz);
    Dz(i, i+1) =  1 / (2*dz);
end

Dz(1, 1:3) = [-3/(2*dz), 4/(2*dz), -1/(2*dz)];
Dz(end, end-2:end) = [1/(2*dz), -4/(2*dz), 3/(2*dz)];

diffmat.Dz = Dz;


par.gamma = gamma;
par.zp    = zp;
par.ell   = 1.0;  


q_init_full = ones(Nz, 1);
q_init_int  = q_init_full(2:end-1); 

tspan = [0, 50]; 

options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, 'Stats', 'off');


[~, q_sol] = ode15s(@(t, q_int) dqdt_1D_snipsnap(t, q_int, diffmat, geom, par, bcType, "nonDim"), ...
    tspan, q_init_int, options);


q_final_int = q_sol(end, :)'; 


switch bcType
    case "NbDt"
        q_l = -(Dz(1,2)*q_final_int(1) + Dz(1,3)*q_final_int(2))/Dz(1,1);
        q_r = par.ell;
    case "NbNt"
        q_l = -(Dz(1,2)*q_final_int(1) + Dz(1,3)*q_final_int(2))/Dz(1,1);
        q_r = -(Dz(end,end-1)*q_final_int(end) + Dz(end,end-2)*q_final_int(end-1))/Dz(end,end);
    case "DbNt"
        q_l = par.ell;
        q_r = -(Dz(end,end-1)*q_final_int(end) + Dz(end,end-2)*q_final_int(end-1))/Dz(end,end);
    case "DbDt"
        q_l = par.ell;
        q_r = par.ell;
end

q_final_full = [q_l; q_final_int; q_r];


total_cells = trapz(z, q_final_full);
end