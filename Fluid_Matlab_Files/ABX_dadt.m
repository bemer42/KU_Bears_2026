%% ABX_dadt

function f = ABX_dadt(a, b)


% Parameters

k19 = 8.33e-4;
K16 = 30.0;

Kt = 39.2115;
Kb = 34.0445;

TCF0 = 15.0;
v18 = 0.1774;

%Right hand side
f = (v18 ./ (1 + (TCF0 .* b) ./ (Kt .* (K16 + b)) + a ./ Kb)) - k19 .* a;



end