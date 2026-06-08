%% ABX_dadt

function F = ABX_dadt(t, A, B,x)

% Discretize time: 
N_t = 1e3;
t_0 = 0;
t_end = 5e3;
t = linspace(t_0,t_end,N_t);


% Parameters

k19 = 8.33e-4;
K16 = 30.0;

Kt = 39.2115;
Kb = 34.0445;

TCF0 = 15.0;
v18 = 0.1774;

F = (v18 ./ (1 + (TCF0 .* B) ./ (Kt .* (K16 + B)) + B ./ Kb)) - k19 .* A;



end