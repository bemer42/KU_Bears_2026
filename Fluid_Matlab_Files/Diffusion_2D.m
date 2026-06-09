%% Diffusion 2D

% Discretize time: 
N_t = 2e3;
t_0 = 0;
t_end = 5;
t = linspace(t_0,t_end,N_t);

% Discretize x space:
N_x = 1e1;
x_0 = 0;
x_end = 1;
x = linspace(x_0, x_end, N_x)';
dx = x(2) - x(1);

% Discretize y space:
N_y = 3e1;
y_0 = 0;
y_end = 1;
y = linspace(y_0, y_end, N_y)';
dy = y(2) - y(1);

%Create mesh grid:

[X,Y] = meshgrid(y,x);



% Parameters
kx = .1;
ky = .1;


% Differentiation Matrix:
Dxx = toeplitz([-2 1 zeros(1,N_x-2)]);
Dxx(1,2) = 2;
Dxx(end, end-1) = 2;
Dxx = kx*Dxx/dx^2;

% Differentiation Matrix:
Dyy = toeplitz([-2 1 zeros(1,N_y-2)]);
Dyy(1,2) = 2;
Dyy(end, end-1) = 2;
Dyy = ky*Dyy/dy^2;



% Initial Condition
f = @(x,y) 2*ones(size(X)).*(x>.25).*(x<.75).*(y>.25).*(y<.75);
U0 = f(X,Y);
U0 = U0(:);

% Define the Right hand side:

dUdt = @(t,U) reshape(Dxx*reshape(U, N_x,N_y) + reshape(U, N_x,N_y)*Dyy', [],1);


% Animation:

for i = 1: N_t

    U_plot = reshape(U(i,:)', N_x, N_y);
    figure(1)
    surf(X,Y,U_plot); hold off;

    if i == 1
        pause
    end

end


