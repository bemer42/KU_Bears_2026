
function dqdt = dqdt_snipsnap(t, Q,z, xp, Dz, alpha, k)

% Note that input column vector q is only interior values

%Left boundary -- no flux;

q_l = (-2*Q(1) + 1/2*Q(2))/(-3/2);

% Right boundary -- dirichlet:

q_r = 1;

%Extend to full q:
q_full = [q_l;Q;q_r];

%Apply PDE to q_full:
arc = 1./sqrt(xp.^2+1);
dqdt = arc.*(Dz*(arc.*alpha./q_full.^2 .*(Dz*q_full)))  + ...
       k .* q_full .*(z<20); 

%Chop to interior q:
dqdt = dqdt(2:end-1);


end