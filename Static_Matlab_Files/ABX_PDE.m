%% ABX PDE
% Discretize time:
N_t = 2e3;
t_0 = 0;
t_end = 1e5;
t = linspace(t_0,t_end,N_t);

% Discretize space:
N_z = 3e2;
z_0 = 0;
z_end = 82;
z = linspace(z_0, z_end, N_z)';
dz = z(2) - z(1);

% Parameters
da = 1e-9;
db = 1e-9;
dx = 1e-9;

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

%% Concentration Animation:

for i = 1:3:N_t

    figure(1)
    plot(z,A(i,:)/max(max(A)), 'k', 'LineWidth',3); hold on;
    plot(z,B(i,:)/max(max(B)), 'r', 'LineWidth',3);
    plot(z,X(i,:)/max(max(X)), 'b', 'LineWidth',3);
    legend('APC','beta-catenin','Axin','fontsize',16)
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

%% Stemness:
S = B.*sqrt(15)./A./sqrt(B + 30);

% Animation:
for i = 1:3:N_t

    figure(2)
    plot(z,S(i,:)/max(max(S)), 'k', 'LineWidth',3); hold off;
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

%% Average Values
N_avg = 300;
A_avg = mean(A(end-N_avg:end,:));
B_avg = mean(B(end-N_avg:end,:));
X_avg = mean(X(end-N_avg:end,:));
S_avg = mean(S(end-N_avg:end,:));

figure(3)
plot(z,A_avg/max(max(A)), 'k:', 'LineWidth',3); hold on;
plot(z,B_avg/max(max(B)), 'r:', 'LineWidth',3);
plot(z,X_avg/max(max(X)), 'b:', 'LineWidth',3);
plot(z,S_avg/max(max(S)), 'm:', 'LineWidth',3);
plot(z,W(z),'g','linewidth',3);
plot(z,max(A(end-N_avg:end,:))/max(max(A)), 'k', 'LineWidth',3); 
plot(z,max(B(end-N_avg:end,:))/max(max(B)), 'r', 'LineWidth',3);
plot(z,max(X(end-N_avg:end,:))/max(max(X)), 'b', 'LineWidth',3);
plot(z,max(S(end-N_avg:end,:))/max(max(S)), 'm', 'LineWidth',3);
plot(z,min(A(end-N_avg:end,:))/max(max(A)), 'k', 'LineWidth',3); 
plot(z,min(B(end-N_avg:end,:))/max(max(B)), 'r', 'LineWidth',3);
plot(z,min(X(end-N_avg:end,:))/max(max(X)), 'b', 'LineWidth',3);
plot(z,min(S(end-N_avg:end,:))/max(max(S)), 'm', 'LineWidth',3);
plot(z,W(z),'g','linewidth',3);
legend('APC','beta-catenin','Axin','Stemness','Wnt','fontsize',16)
title('Average Concentrations along Crypt','fontsize',20)
xlabel('crypt level','fontsize',16)
ylabel('normalized concentration','fontsize',16)
set(gca,'fontsize',14)
grid on
grid minor
