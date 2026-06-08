%% ABX_dbdt

function G = ABX_dbdt(t, A, B,x, W)

% Discretize time: 
N_t = 1e3;
t_0 = 0;
t_end = 5e3;
t = linspace(t_0,t_end,N_t);


% Parameters

k13 = 2.57e-4; 
K7 = 50.0; 
K17 = 1200; 
K20 = 1; 
K21 = 1; 
v12 = 0.423;


end