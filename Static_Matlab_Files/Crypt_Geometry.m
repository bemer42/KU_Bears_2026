%% Crypt Geometry
clear, close all, clc

% Discretize theta space:
Nth    = 8e1;
th_0   = 0;
th_end = 2*pi;
th     = linspace(th_0, th_end, Nth)';
dth    = th(2) - th(1);

% Discretize z space:
Nz  = 8e1;
z_0 = 0;
L   = 78.8783;
z   = linspace(z_0, L, Nz)';
dz  = z(2) - z(1); 

% Create a Meshgrid:
[Th,Z] = meshgrid(th,z);

% Define radius function: 
rb = 41/2/pi;
rt = 10/pi;
a  = 0.3;

r   = @(th,z) rb.*(1-exp(-a.*z)) + rt.*exp(a.*(z-L)); 
rz  = @(th,z) -a*rb.*exp(-a.*z) + a*rt.*exp(a.*(z-L));
rth = @(th,z) zeros(size(th)); 

% Define surface generating function: 
X = @(th,z) r(th,z).*cos(th);
Y = @(th,z) r(th,z).*sin(th); 

% Generate surface plot: 
figure(1)
mesh(X(Th,Z),Y(Th,Z),Z,'EdgeColor',[.5 .5 .5])
hold on
plot3(r(th,z),zeros(size(z)),z,'k','linewidth',10)
set(gca,'fontsize',42)
title('Crypt Geometry','fontsize',55,'interpreter','latex')
xlabel('$x$','fontsize',50,'interpreter','latex')
ylabel('$y$','fontsize',50,'interpreter','latex')
zlabel('$z$','fontsize',50,'interpreter','latex')
legend('$\vec{X}(\theta,z)$','$\vec{R}(z)$','fontsize',50,'interpreter','latex')
% shading interp
% lighting gouraud
grid on 
grid minor
box on
% colormap summer
view([1 -1 1])
xlim([-40 40])
ylim([-40 40])
zlim([-2 80])

% Arc-Length: 
Arc_Length = trapz(z,sqrt(rz(th,z).^2+1));

% Surface Area: 
Surface_Area = trapz(th,trapz(z,sqrt(rth(Th,Z).^2 + r(Th,Z).^2.*(rz(Th,Z).^2+1))));

% Optimizing the value of L: 
rz_fit = @(x) -a*rb.*exp(-a.*linspace(0,x,1e3)) + a*rt.*exp(a.*(linspace(0,x,1e3)-x));
f = @(x) trapz(linspace(0,x,1e3),sqrt(rz_fit(x).^2+1));
L = fzero(@(x) f(x) - 82,78.87);




