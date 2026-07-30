%% ABX PDE
clear, close all, clc

% Discretize time:
N_t = 3e3;
t_0 = 0;
t_end = 3e5;
t = linspace(t_0,t_end,N_t);

% Discretize space:
N_z = 1e2;
z_0 = 0;
z_end = 82;
z = linspace(z_0, z_end, N_z)';
dz = z(2) - z(1);

% Parameters
da = 0e-9;
db = 0e-9;
dx = 0e-9;

% Differentiation Matrix:
Dzz = toeplitz([-2 1 zeros(1,N_z-2)]);
Dzz(1,2) = 2;
Dzz(end, end-1) = 2;
Dzz = Dzz/dz^2;

Da = da*Dzz;
Db = db*Dzz;
Dx = dx*Dzz;

% Define initial condition:
A0 = @(z) 25 + (114-25)/82*z;
B0 = @(z) 289 + (23-289)/82*z;
X0 = @(z) 9e-4 + (5e-4 - 9e-4)/82*z;

U0 = [A0(z); B0(z); X0(z)];

% Define a Wnt signal:
W = @(z) exp(-z/10);

% Define the Right hand side:
dUdt = @(t,U) [ABX_dadt(U(1:N_z),U(N_z+1:2*N_z)) + Da*U(1:N_z);
               ABX_dbdt(U(1:N_z),U(N_z+1:2*N_z),U(2*N_z+1:end),W(z)) + Db*U(N_z+1:2*N_z);
               ABX_dxdt(U(1:N_z),U(N_z+1:2*N_z),U(2*N_z+1:end),W(z)) + Dx*U(2*N_z+1:end)];

% Solve the heat equation:
tic
options = odeset('Stats', 'on');
[t,U] = ode23s(dUdt, t, U0, options);
toc

A = U(:,1:N_z);
B = U(:,N_z+1:2*N_z);
X = U(:,2*N_z+1:end);

% Plot of Average Values
N_avg = 300;
A_avg = mean(A(end-N_avg:end,:));
B_avg = mean(B(end-N_avg:end,:));
X_avg = mean(X(end-N_avg:end,:));

figure(1)
plot(z,max(A(end-N_avg:end,:))/max(max(A)), 'g', 'LineWidth',3); hold on;
plot(z,max(B(end-N_avg:end,:))/max(max(B)), 'r', 'LineWidth',3);
plot(z,max(X(end-N_avg:end,:))/max(max(X)), 'b', 'LineWidth',3);
plot(z,W(z),'m:','linewidth',3);
plot(z,min(A(end-N_avg:end,:))/max(max(A)), 'g', 'LineWidth',3); 
plot(z,min(B(end-N_avg:end,:))/max(max(B)), 'r', 'LineWidth',3);
plot(z,min(X(end-N_avg:end,:))/max(max(X)), 'b', 'LineWidth',3);
plot(z,A_avg/max(max(A)), 'g:', 'LineWidth',3); 
plot(z,B_avg/max(max(B)), 'r:', 'LineWidth',3);
plot(z,X_avg/max(max(X)), 'b:', 'LineWidth',3);
legend('APC','beta-catenin','Axin','Wnt','fontsize',16)
title('Average Concentrations along Crypt','fontsize',20)
xlabel('crypt level','fontsize',16)
ylabel('normalized concentration','fontsize',16)
set(gca,'fontsize',14)
grid on
grid minor

%% Concentration Animation:

dt = round(.01*N_t);

for i = 1:dt:N_t

    figure(2)
    plot(z,A(i,:)/max(max(A)), 'g', 'LineWidth',3); hold on;
    plot(z,B(i,:)/max(max(B)), 'r', 'LineWidth',3);
    plot(z,X(i,:)/max(max(X)), 'b', 'LineWidth',3);
    plot(z,W(z),'m:','LineWidth',3)
    legend('APC','beta-catenin','Axin','Wnt','fontsize',16)
    title('Protein Concentration along Crypt','fontsize',20)
    xlabel('crypt level','fontsize',16)
    ylabel('normalized concentration','fontsize',16)
    set(gca,'fontsize',14)
    hold off;
    if i == 1
        pause
    end
    grid on
    grid minor
    ylim([0,1])

end

%% Stemness Animation:
S = B.*sqrt(15)./A./sqrt(B + 30);
dt = round(.01*N_t);

for i = 1:dt:N_t

    figure(3)
    plot(z,S(i,:)/max(max(S)), 'k', 'LineWidth',3); hold on;
    plot(z,W(z),'m:','LineWidth',3)
    title('Stemness along Crypt','fontsize',20)
    xlabel('crypt level','fontsize',16)
    ylabel('normalized concentration','fontsize',16)
    set(gca,'fontsize',14)
    hold off;
    if i == 1
        pause
    end
    grid on
    grid minor
    ylim([0,1])

end

%%  Surface Plots: 
[T,Z] = meshgrid(t,z);

% APC:
figure(4)
surf(T,Z,A'); 
hold on;
shading interp
title('APC','fontsize',20, 'FontName', 'Times New Roman')
xlabel('time (t)','fontsize',16, 'FontName', 'Times New Roman')
ylabel('crypt level (z)','fontsize',16, 'FontName', 'Times New Roman')
zlabel('APC (nM)', 'FontSize',16,'FontName', 'Times New Roman')
set(gca,'fontsize',14)

% APC trajectories:
for i = 1:6:N_z

    figure(4)
    plot3(t,z(i)*ones(size(t)),A(:,i),'k','linewidth',3); hold on;

end

%  Beta Catenin: 
figure(5)
surf(T,Z,B'); 
hold on;
shading interp
title('Beta-Catenin','fontsize',20,'FontName', 'Times New Roman')
xlabel('time (t)','fontsize',16, 'FontName', 'Times New Roman')
ylabel('crypt level (z)','fontsize',16, 'FontName', 'Times New Roman')
zlabel('Beta-Catenin (nM)', 'FontSize',16, 'FontName', 'Times New Roman')
set(gca,'fontsize',14)

% Beta Catenin trajectories
for i = 1:6:N_z

    figure(5)
    plot3(t,z(i)*ones(size(t)),B(:,i),'k','linewidth',3); hold on;

end

%  Axin: 
figure(6)
surf(T,Z,X'); 
hold on;
shading interp
title('Axin','fontsize',20,'FontName', 'Times New Roman')
xlabel('time (t)','fontsize',16, 'FontName', 'Times New Roman')
ylabel('crypt level (z)','fontsize',16,'FontName', 'Times New Roman')
zlabel('Axin (nM)', 'FontSize',16, 'FontName', 'Times New Roman')
set(gca,'fontsize',14)

% Axin trajectories:
for i = 1:6:N_z

    figure(6)
    plot3(t,z(i)*ones(size(t)),X(:,i),'k','linewidth',3); hold on;

end
