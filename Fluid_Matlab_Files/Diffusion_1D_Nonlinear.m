%% Diffusion 1D Nonlinear

% Discretize time: 
N_t = 1e3;
t_0 = 0;
t_end = 5;
t = linspace(t_0,t_end,N_t);

% Discretize space:
N_z = 100;
z_0 = 0;
z_end = 82;
z = linspace(z_0, z_end, N_z)';
dz = z(2) - z(1);


% Define parameters

r_b = 5;
r_t = 10;   
a = 0.02;   
alpha = 0.05;
   

% Crypt width
x = r_b*(1 - exp(-a * z)) + r_t * exp(a * (z - z_end));

dx_dz = gradient(x, dz);

% Arc length
ds = sqrt(1+dx_dz.^2) * dz;

% positions s along wall
s = [0; cumsum(ds(1:end-1))]; %looked this one up to find total distance along curve
ds_size = mean(diff(s)); %also looked up it will supposedly calculate distance 
                         %between consecutive pairs of points along arc
                         %length coordinate vector s and average them


% Differentiation Matrix:
Dss = toeplitz([-2 1 zeros(1,N_z-2)]);
Dss(1,2) = 2;
Dss(end, end-1) = 2;
Dss = Dss/ds_size^2;



% Initial Condition
f = @(s) cos(s*2*pi) + 1;
U0 = f(s);

% Define the Right hand side:

dUdt = @(t,U) Dss * ((alpha ./ (U.^2 + 1e-6)) .* (Dss*U));

% Solve the heat equation:
tic
options = odeset('Stats', 'on');
[t,U] = ode23s(dUdt, t, U0, options);
toc


% Animation:

for i = 1: 5: N_t

    figure(1)
    plot(s,U(i,:), 'k', 'LineWidth',3); hold off;
    if i == 1
        pause
    else
        drawnow
    end
    grid on
    grid minor
    ylim([0,2.5])

end


