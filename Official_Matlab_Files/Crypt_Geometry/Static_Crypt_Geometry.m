%% Crypt Geometry
% This code will output the general static crypt geometry by first creating
% the one-dimensional curve that represents a chain of cells, then rotating
% it about the z-axis to create the crypt.  Based on arc-length and surface
% area, the optimal value of L is also determined at the end of the code.
% The parameters used in this code represent the baseline parameters for a
% fully formed single crypt. A figure if formed using a uniform grid in
% (theta,z) space.
clear, close all, clc; figure_on = 0; L_optimize = 0;

% Discretize theta space:
Nth    = 1e4;
th_0   = 0;
th_end = 2*pi;
th     = linspace(th_0, th_end, Nth)';
dth    = th(2) - th(1);

% Discretize z space:
Nz  = 1e4;
z_0 = 0;
L   = 79.3515;
z   = linspace(z_0, L, Nz)';
dz  = z(2) - z(1);

% Create a Meshgrid:
[Z,Th] = meshgrid(z,th);

% Define radius function:
rb = 18/pi;
rt = 10/pi;
a  = 0.3;

% Radius function with partial derivatives:
r   = @(th,z) rb.*(1-exp(-a.*z)) + rt.*exp(a.*(z-L));
rz  = @(th,z) -a*rb.*exp(-a.*z) + a*rt.*exp(a.*(z-L));
rth = @(th,z) zeros(size(th));

% Define surface generating function:
X = @(th,z) r(th,z).*cos(th);
Y = @(th,z) r(th,z).*sin(th);

% Surface plot for 1D and 2D generating functions:
while figure_on == 1

    fig = figure(1);
    fig.WindowState = 'maximized';

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
    grid on
    grid minor
    box on
    view([1 -1 1])
    xlim([-40 40])
    ylim([-40 40])
    zlim([-2 80])

    figure_on = 0;

end

% Arc-Length:
Arc_Length_Cum = cumtrapz(z,sqrt(rz(th,z).^2+1));
Arc_Length     = Arc_Length_Cum(end)

% Surface Area:
Surface_Area_Cum = cumtrapz(th,trapz(z,sqrt(rth(Th,Z).^2 + r(Th,Z).^2.*(rz(Th,Z).^2+1))))
Surface_Area     = Surface_Area_Cum(end)

% Optimizing the value of L:
while L_optimize == 1

    N_opt = 1e4;
    rz_fit = @(L) -a*rb.*exp(-a.*linspace(0,L,N_opt)) + a*rt.*exp(a.*(linspace(0,L,N_opt)-L));
    f = @(L) trapz(linspace(0,L,N_opt),sqrt(rz_fit(L).^2+1));
    L_opt = fzero(@(x) f(x) - 82,79)
    L_optimize = 0;

end

% Determining z_0 using Arc-Length:
pos_al = find(Arc_Length_Cum >= 1);
z_0_al = z(pos_al(1))

% Determining z_0 using Surface Area:
pos_sa = find(Surface_Area_Cum >= 32);
z_0_sa = z(pos_sa(1))

