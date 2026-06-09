%% ABX_dxdt

function h = ABX_dxdt(a, b,x, w)

% Discretize time: 
N_t = 1e3;
t_0 = 0;
t_end = 5e3;
t = linspace(t_0,t_end,N_t);


% Parameters

k1 = 0.182;
k2 = 1.82e-2;
k3 = 5.0e-2;
k4 = 0.267;
k5 = 0.133;
k6 = 9.09e-2;
k_6 = 0.909;
k9 = 206.0;
k10 = 206.0;
k11 = 0.417;
k13 = 2.57e-4;
k15 = 0.333;
k19 = 8.33e-4;
K7 = 50.0;
K8 = 120.0;
K16 = 30.0;
K17 = 1200;
K20 = 1;
K21 = 1;
Kt = 39.2115;
Kb = 34.0445;
GSK0 = 50.0;
Dsh0 = 100.0;
TCF0 = 15.0;
v18 = 0.1774;
v12 = 0.423;
v14 = 8.22e-5;
Km = 98;

f_fun = @(w) (k1*k4*k6*k9*K21/ (k5*K7*K8))*((w+k2/k1)./(k1*(Dsh0*k3+k_6)*w+k2*k_6));
g_fun = @(a, b, x , w) ((K21+x)./GSK0)+(((k9+k10).*a.*x)./(k9*k10*GSK0)).*(((k4+k5)*K8*k10)./(k4*(k9+k10))+b).*f_fun(w);

F_fun = @(a, b) (v18./(1+(TCF0.*b)./(Kt.*(K16+b))+b./Kb))-k19.*a;
G_fun = @(a, b, x, w) v12 - (k13 + a.*x.*(f_fun(w)./g_fun(a,b,x, w))+F_fun(a,b)./K17).*b ;
H_fun = @(a, b, x)  v14 - ((k15.*a)./(Km+a)+F_fun(a,b)./K7).*x;



h = ((-x/K20).*G_fun(a,b,x,w)+(1+a/K7+x/K20).*H_fun(a,b,x))./ ...
    ((1+a/K7+x/K20).*(1+ a/K7+ b/K20 + K21./((K21+x).*g_fun(a,b,x,w)))-(b.*x/(K20^2)));


end