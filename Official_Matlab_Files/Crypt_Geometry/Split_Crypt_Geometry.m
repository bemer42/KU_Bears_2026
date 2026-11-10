%% A Static Split Crypt
close all; clear all; clc

% Discretize theta space:
Nth    = 1e2;
th_0   = 0;
th_end = 2*pi;
th     = linspace(th_0, th_end, Nth)';
dth    = th(2) - th(1);

% Discretize z space:
Nz  = 1e2;
z_0 = 0;
L   = 79.3515;
z   = linspace(z_0, L, Nz)';
dz  = z(2) - z(1); 

% Create a Meshgrid:
[Z,Th] = meshgrid(z,th);

% Define geometry parameters:
r_b   = 18/pi;
r_t   = 10/pi;
a     = 0.3;

% Splitting parameters:
p     = 0.99;
split = 50;

% Define static radius: 
R_s   = r_b*(1-exp(-a*Z))+r_t*exp(a*(Z-L));

% Define current crypt structure: 
A = p./(1 + exp(-.2*(split-Z)));
C = A .* (1 + sqrt(5))/2.*R_s; 
R = sqrt(C.^2.*cos(2*Th)+...
    sqrt((R_s.^2+R_s.*C).^2-C.^4.*sin(2*Th).^2));

figure(1)
mesh(R.*cos(Th),R.*sin(Th),Z,'EdgeColor',[.5 .5 .5])
set(gca,'fontsize',42)
title('Crypt Fission Geometry','fontsize',55,'interpreter','latex')
xlabel('$x$','fontsize',50,'interpreter','latex')
ylabel('$y$','fontsize',50,'interpreter','latex')
zlabel('$z$','fontsize',50,'interpreter','latex')
legend('$\vec{X}(\theta,z)$','fontsize',50,'interpreter','latex')
grid on 
grid minor
box on
% colormap summer
view([1 -1 1])
xlim([-40 40])
ylim([-40 40])
zlim([-2 80])

% Differentiation Matrix in theta
Dth = spdiags([-1/2*ones(Nth,1) 1/2*ones(Nth,1)],[-1 1],Nth,Nth);
Dth(1,end) = -1/2;
Dth(end,1) =  1/2;
Dth = Dth/dth;

% Define the differentiation matrix in z
Dz = diffmat_nonuniform(z);
% Define partial derivatives: 
Rz  = R*Dz';
Rth = Dth*R;

% Define surface area element: 
G = sqrt(Rth.^2+R.^2.*(Rz.^2 + 1));

% Output the total number of cells via surface area: 
Total_Cells = trapz(z, trapz(th, G))

