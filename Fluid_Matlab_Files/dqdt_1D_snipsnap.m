
function dqdt = dqdt_1D_snipsnap(t, q_int, z, rp, Dz, alpha, k)

% Note that input column vector q is only interior values

% Left boundary -- no flux;
q_l = (-2*q_int(1) + 1/2*q_int(2))/(-3/2);

% Right boundary -- dirichlet:
q_r = 1;

%Extend to full q:
q_full = [q_l; q_int; q_r];

%Apply PDE to q_full:
arc  = 1./sqrt(rp.^2+1);
dqdt = arc.*(Dz*(arc.*alpha./q_full.^2 .*(Dz*q_full)))  + ...
       k .* q_full .*(z<20); 

%Chop to interior q:
dqdt = dqdt(2:end-1);

end