%% Analysis of 2D Nondimensional model:
clear; close all; clc

% Define parameter vectors:
L         = 78.8783;
N         = 2e1;
% gz_vec    = logspace(-1.25,-.25,N);
gz        = .17;
gth_vec   = logspace(-1.25,-.25,N);
zp_vec    = linspace(10,30,N)/L;

% Define bc type:
bcType = "NbDt";

% Define Mesh:
[Gth,Zp] = meshgrid(gth_vec,zp_vec);

% Preallocate storage matrices:
C = zeros(N);
P = zeros(N);
% For loop
tic
for i = 1:N
    for j = 1:N

        [C(i,j), P(i,j)] = ND_2D_Function(gz,gth_vec(i),zp_vec(j),bcType);

    end
end
toc

%% Surface Plots:

% Total Cell Contours: 
MC  = max(max(C))+1; 
CON = [3280 3380];

C1 = contourc(gth_vec,zp_vec,C,[CON(1) CON(1)]);
gc1 = C1(1,2:end); 
zc1 = C1(2,2:end); 

C2 = contourc(gth_vec,zp_vec,C,[CON(2) CON(2)]);
gc2 = C2(1,2:end); 
zc2 = C2(2,2:end); 

% Total Cells:
figure(1)
surf(log10(Gth),Zp,C)
hold on
plot3(log10(gc1),zc1,MC*ones(size(zc1)),'w','linewidth',4)
plot3(log10(gc2),zc2,MC*ones(size(zc2)),'w','linewidth',4)
set(gca,'fontsize',42)
title('Total Number of Cells','fontsize',55,'interpreter','latex')
xlabel('$\gamma_z$','fontsize',50,'interpreter','latex')
ylabel('$z_p$','fontsize',50,'interpreter','latex')
zlabel('$\mathcal{C}$','fontsize',50,'interpreter','latex')
view([0 0 1])
grid on
grid minor
box on
shading interp
colormap turbo
colorbar
lighting gouraud
xlim([min(log10(gth_vec)) max(log10(gth_vec))])
ylim([min(zp_vec) max(zp_vec)])
zlim([min(min(C)) MC])


% Proliferative Cell Contours: 
MP  = max(max(P))+1; 
CON = [800 1000];

P1 = contourc(gth_vec,zp_vec,P,[CON(1) CON(1)]);
gp1 = P1(1,2:end); 
zp1 = P1(2,2:end); 

P2 = contourc(gth_vec,zp_vec,P,[CON(2) CON(2)]);
gp2 = P2(1,2:end); 
zp2 = P2(2,2:end); 

% Proliferative Cells:
figure(2)
surf(log10(Gth),Zp,P)
hold on
plot3(log10(gp1),zp1,MP*ones(size(zp1)),'w','linewidth',4)
plot3(log10(gp2),zp2,MP*ones(size(zp2)),'w','linewidth',4)
set(gca,'fontsize',42)
title('Total Number of Proliferative Cells','fontsize',55,'interpreter','latex')
xlabel('$\gamma_z$','fontsize',50,'interpreter','latex')
ylabel('$z_p$','fontsize',50,'interpreter','latex')
zlabel('$\mathcal{P}$','fontsize',50,'interpreter','latex')
view([0 0 1])
colorbar 
grid on
grid minor
box on
shading interp
colormap turbo
lighting gouraud
xlim([min(log10(gth_vec)) max(log10(gth_vec))])
ylim([min(zp_vec) max(zp_vec)])
zlim([min(min(P)) MP])


% Crypt Renewal Contours: 
% MR  = max(max(R))+1; 
% CON = [4 5];
% 
% R1 = contourc(gth_vec,zp_vec,R,[CON(1) CON(1)]);
% gr1 = R1(1,2:end); 
% zr1 = R1(2,2:end); 
% 
% R2 = contourc(gth_vec,zp_vec,R,[CON(2) CON(2)]);
% gr2 = R2(1,2:end); 
% zr2 = R2(2,2:end); 

% Crypt Renewal Time
% figure(3)
% surf(log10(Gth),Zp,R)
% hold on
% plot3(log10(gr1),zr1,MR*ones(size(zr1)),'w','linewidth',4)
% plot3(log10(gr2),zr2,MR*ones(size(zr2)),'w','linewidth',4)
% set(gca,'fontsize',42)
% title('Crypt Renewal Time','fontsize',55,'interpreter','latex')
% xlabel('$\gamma$','fontsize',50,'interpreter','latex')
% ylabel('$z_p$','fontsize',50,'interpreter','latex')
% zlabel('$\mathcal{R}$','fontsize',50,'interpreter','latex')
% view([0 0 1])
% colorbar 
% grid on
% grid minor
% box on
% shading interp
% colormap turbo
% lighting gouraud
% xlim([min(log10(gth_vec)) max(log10(gth_vec))])
% ylim([min(zp_vec) max(zp_vec)])
% zlim([min(min(R)) MR])

