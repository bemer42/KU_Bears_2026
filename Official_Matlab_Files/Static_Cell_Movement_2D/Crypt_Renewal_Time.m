%% Calculating Crypt renewal time given a velocity field.
close all; clear; clc;

% Load in Velocity Field and other data:
load('Velocity_Field_10-6.mat')

% -------------------------
% Initial condition vector:
% -------------------------
N_vec    = 10;
th_0_vec = linspace(th(2), th(end-1), N_vec);
z_0      = 1.64;

% -------------------------
% Initialize renewal time:
% -------------------------
Time = zeros(N_vec,1);

% -------------------------
% Loop through intial conditions:
% -------------------------
ex = @(x) [x; x(1,:)];
fig = figure(1);
fig.WindowState = 'maximized';
for i = 1:N_vec

    % -------------------------
    % Period parameters:
    % -------------------------
    thMin = 0;
    P     = 2*pi;

    % -------------------------
    % Extend theta grid to include 2*pi and extend fields:
    % -------------------------
    th_ext  = [th; 2*pi];
    Uth_ss = ex(Vth_ss) ./ ex(sqrt(G));
    Uz_ss  = ex(Vz_ss) ./ ex(sqrt(G));

    % -------------------------
    % Create gridded interpolants:
    % -------------------------
    Fth = griddedInterpolant({th_ext, z}, Uth_ss, "linear", "nearest");
    Fz  = griddedInterpolant({th_ext, z}, Uz_ss,  "linear", "nearest");

    % -------------------------
    % Define ODE:
    % -------------------------
    rhs = @(t,y) rhsPeriodic(t, y, Fth, Fz, thMin, P);

    % -------------------------
    % Event to stop when z reaches L
    % -------------------------
    opts = odeset("Events", @(t,y) eventTop(t,y,L), ...
        "RelTol", 1e-7, "AbsTol", 1e-9);

    % -------------------------
    % Integrate:
    % -------------------------
    tSpan = [0, 1e6];
    y0    = [th_0_vec(i); z_0];
    [t,y,te,ye,ie] = ode45(rhs, tSpan, y0, opts);

    % -------------------------
    % System check:
    % -------------------------
    if isempty(te)
        warning("Trajectory did not reach z = L within tSpan (or Vz never becomes positive enough).");
        travelTime = NaN;
    else
        travelTime = te(end);
    end

    % -------------------------
    % Define trajectory:
    % -------------------------
    th_traj = y(:,1);
    z_traj  = y(:,2);

    % -------------------------
    % wrap theta into [0,2pi):
    % -------------------------
    th_traj = thMin + mod(th_traj - thMin, P);
    r = @(z) r_b*(1-exp(-a*z)) + r_t*exp(a*(z-L));

    % -------------------------
    % Store renewal time:
    % -------------------------
    Time(i) = travelTime;

    % -------------------------
    % Plot:
    % -------------------------
    figure(1)
    subplot(1,2,1)
    surf([Th; 2*pi*ones(1,Nz)], ex(Z), ex(Q_ss));
    hold on
    quiver3(Th(2:dm:end,2:dm:end-dm), Z(2:dm:end,2:dm:end-dm), Q_ss(2:dm:end,2:dm:end-dm), ...
        Uth_ss(2:dm:end,2:dm:end-dm), Uz_ss(2:dm:end,2:dm:end-dm), zeros(size(Th(2:dm:end,2:dm:end-dm))), ...
        1.0,'color',[0.15 0.15 0.15],'Linewidth',2.5,'MaxHeadSize',.008);
    plot3([th; 2*pi], Zp([th; 2*pi], t_end), max(max(Q_ss))*ones(Nth+1,1), 'w:', 'linewidth', 5)
    plot3(th_traj, z_traj, max(max(Q_ss))*ones(size(z_traj)), 'm', 'linewidth', .5)
    plot3(th_0_vec(i), z_0, max(max(Q_ss)), 'go', 'linewidth', 6)
    plot3(th_traj(end), z_traj(end), max(max(Q_ss)), 'ro', 'linewidth', 6)
    set(gca,'fontsize',42)
    title('Cell Density and Velocity','fontsize',55,'interpreter','latex')
    xlabel('$\theta$','fontsize',50,'interpreter','latex')
    ylabel('$z$','fontsize',50,'interpreter','latex')
    zlabel('$q_s(\theta,z)$ and $\vec{v}_s(\theta,z)$','fontsize',50,'interpreter','latex')
    legend('$q_s(\theta,z)$','$\vec{v}_s(\theta,z)$','fontsize',50,'interpreter','latex','location','northwest')
    grid on; grid minor; box on
    colormap turbo
    shading interp
    colorbar
    caxis([1 max(max(Q_ss))])
    xlim([0 2*pi])
    ylim([0 L])
    view([0 0 1])
    hold on

    subplot(1,2,2)
    surf(ex(R).*cos([Th; 2*pi*ones(1,Nz)]), ...
        ex(R).*sin([Th; 2*pi*ones(1,Nz)]), ex(Z), ex(Q_ss))
    colormap turbo
    shading interp
    lighting gouraud
    hold on
    mesh(1.2*ex(R).*cos([Th; 2*pi*ones(1,Nz)]), ...
        1.2*ex(R).*sin([Th; 2*pi*ones(1,Nz)]), Zp_plot, 'edgecolor',[.5 .5 .5])
    plot3(r(z_traj).*cos(th_traj), r(z_traj).*sin(th_traj), z_traj, 'm', 'linewidth', .5)
    plot3(r(z_0).*cos(th_0_vec(i)), r(z_0).*sin(th_0_vec(i)), z_0, 'go', 'linewidth', 6)
    plot3(r(z_traj(end)).*cos(th_traj(end)), r(z_traj(end)).*sin(th_traj(end)), z_traj(end), 'ro', 'linewidth', 6)
    set(gca,'fontsize',42)
    title('Cell Density on Crypt','fontsize',55,'interpreter','latex')
    xlabel('$x$','fontsize',50,'interpreter','latex')
    ylabel('$y$','fontsize',50,'interpreter','latex')
    zlabel('$z$','fontsize',50,'interpreter','latex')
    legend('$q(\theta,z,t)$','fontsize',50,'interpreter','latex','location','northwest')
    grid on; grid minor; box on
    colorbar
    xlim([-40 40]); ylim([-40 40]); zlim([0 80])
    caxis([1 max(max(Q_full))])
    hold on
    if i == 1
        pause
    end

end

% -------------------------
% Print renewal times:
% -------------------------
fprintf("Travel time to reach z=L (converted to hours): %.6g\n", mean(Time));
fprintf("Travel time to reach z=L (converted to days):  %.6g\n", mean(Time)/24);

subplot(1,2,2)
text(.1,.5, ['$\mathcal{R} = $ ' num2str(mean(Time)/24)], ...
    'units','normalized','fontsize',40,'interpreter','latex')

%% Helpers (local functions)
function [value, isterminal, direction] = eventTop(~, y, L)
value      = y(2) - L;
isterminal = 1;
direction  = +1;
end

function dydt = rhsPeriodic(~, y, Fth, Fz, thMin, P)
th = thMin + mod(y(1) - thMin, P);
z  = y(2);
dydt = [Fth(th,z); Fz(th,z)];
end