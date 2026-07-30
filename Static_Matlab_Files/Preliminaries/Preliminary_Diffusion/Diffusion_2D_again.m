%% Spatial ODE with space varying Wnt Signal
close; clear all; clc

tic 

% Discretize time: 
N_t   = 2e3;
t_0   = 0;
t_end = 2;
t     = linspace(t_0,t_end,N_t);

% Discretize x space: 
N_x   = 1e1; 
x_0   = 0; 
x_end = 1; 
x     = linspace(x_0,x_end,N_x)';
dx    = x(2) - x(1);

% Discretize y space: 
N_y   = 2e1; 
y_0   = 0; 
y_end = 1; 
y     = linspace(y_0,y_end,N_y)';
dy    = y(2) - y(1);

% Create Mesh grid
[Y,X] = meshgrid(y,x);

% Parameters:
kx = .1;
ky = .1; 

% Differentiation Matrix: 
Dxx            = toeplitz([-2 1 zeros(1,N_x-2)]);
Dxx(1,2)       = 2;
Dxx(end,end-1) = 2;
Dxx            = kx*Dxx/dx^2;

% Differentiation Matrix: 
Dyy            = toeplitz([-2 1 zeros(1,N_y-2)]);
Dyy(1,2)       = 2;
Dyy(end,end-1) = 2;
Dyy            = ky*Dyy/dy^2;

% Initial condition: 
f = @(x,y) 2*ones(size(x)).*(x>.25).*(x<.75).*(y>.25).*(y<.75);

U0 = f(X,Y);
U0 = U0(:);

% Define the Right hand side: 
dUdt = @(t,U) reshape(Dxx*reshape(U,N_x,N_y) + reshape(U,N_x,N_y)*Dyy',[],1);

% Solve the Heat equation:
tic
options = odeset('Stats','on');
[t,U] = ode23s(dUdt,t,U0,options);
toc

% % Animation: 
for i = 1:N_t
   
    U_plot = reshape(U(i,:),N_x,N_y);

    figure(1)
    surf(X,Y,U_plot); hold off;
    if i ==1 
        pause
    end

end
% 
% 
% % Surface Plot: 
% [T,X] = meshgrid(t,x);
% 
% figure(2)
% surf(T,X,U')
% shading interp
% hold on
% 
% for i = 1:6:N_x
% 
%     plot3(t,x(i)*ones(size(t)),U(:,i),'k','linewidth',2)
%     if i == 1
%         plot3(zeros(size(x)),x,U(1,:),'k','linewidth',5)
%     end
%     hold on
% end