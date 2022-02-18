clear; close all; clc;
% This script solves Model A and saves results into a .mat file + generates
% fun animated .gifs to show how the diffusion parameter affects model
% behavior

%%% parameters
D     = 80;   %apparent diffusion (ml/s)
Pair  = 150; %atmospheric oxygen partial pressure (mmHg)
Pin   = 45;  %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)
Vvasc = 1;   %volume of vascular space (ml)
Valv  = 1;   %alveolar volume (ml)

Vp = 5; %ventilation flow (ml/s)
Qp = 5; %blood flow (ml/s)

% load('ModelC_optimization_results.mat','D_opt')
% D = D_opt;

par = [D Pair Pin Vvasc Valv Vp  Qp]; %parameters vector for ODE system

%%% grid for contour plot
de = 0.5;
q  = 0:de:70;
v  = 0:de:70;
[Q,V] = meshgrid(q,v);

Dv = 0.1:0.1:10; %vector of diffusion parameters to explore
Dpvasc = zeros(length(v),length(q),length(Dv));
Dpalv  = zeros(length(v),length(q),length(Dv));
Ddpac  = zeros(length(v),length(q),length(Dv));

% computing contour plot
palv  = (Pin*Q.*(exp(-D./Q)-1)-Pair*V)./(Q.*(exp(-D./Q)-1)-V);
pvasc = palv + (Pin-palv).*exp(-D./Q);
dpac  = palv - pvasc;

% contour plots for vector of diffusion parameters to explore
for j = 1:length(Dv)
    Dpalv(:,:,j)  = (Pin*Q.*(exp(-Dv(j)./Q)-1)-V.*Pair)./(Q.*(exp(-Dv(j)./Q)-1)-V);
    Dpvasc(:,:,j) = Dpalv(:,:,j) + (Pin-Dpalv(:,:,j)).*exp(-Dv(j)./Q);
    Ddpac(:,:,j)  = Dpalv(:,:,j) - Dpvasc(:,:,j);
end

save('ModelC_results.mat')

figure;
surf(Q,V,pvasc, 'EdgeColor','none')
title('Pvasc')
xlabel('q')
ylabel('v')
% set(gca,'zscale','log')
% axis equal

figure;
surf(Q,V,palv, 'EdgeColor','none')
title('Palv')
xlabel('q')
ylabel('v')
% set(gca,'zscale','log')
% axis equal

figure;
surf(Q,V,dpac, 'EdgeColor','none')
title('\Delta P')
xlabel('q')
ylabel('v')
% set(gca,'zscale','log')
% axis equal

CON = -30:5:150;
x = 0:10;
figure; % Contour Plot
subplot(1,3,1)
contour(Q,V,pvasc,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Vascular (mmHg)')
xlabel('Blood Flow (ml/s)')
ylabel('Air Flow (ml/s)')
axis equal
grid on

subplot(1,3,2)
contour(Q,V,palv,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Alveolar Space (mmHg)')
xlabel('Blood Flow (ml/s)')
axis equal
grid on

subplot(1,3,3)
contour(Q,V,dpac,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Alv-Vasc Gradient (mmHg)')
xlabel('Blood Flow (ml/s)')
axis equal
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% big old plots
% figure; % Contour Plot
% subplot(3,3,1)
% contour(Q,V,squeeze(Dpvasc(:,:,1)),CON,'ShowText','on')
% hold on
% plot(x,x,'k--')
% set(gca,'fontsize',18)
% title('Vascular (mmHg)')
% ylabel('Air Flow (ml/s)')
% axis equal
% grid on
% 
% subplot(3,3,2)
% contour(Q,V,squeeze(Dpalv(:,:,1)),CON,'ShowText','on')
% hold on
% plot(x,x,'k--')
% set(gca,'fontsize',18)
% title('Alveolar Space (mmHg)')
% axis equal
% grid on
% 
% subplot(3,3,3)
% contour(Q,V,squeeze(Ddpac(:,:,1)),CON,'ShowText','on')
% hold on
% plot(x,x,'k--')
% set(gca,'fontsize',18)
% title('Alv-Vasc Gradient (mmHg)')
% axis equal
% grid on
% 
% subplot(3,3,4)
% contour(Q,V,squeeze(Dpvasc(:,:,2)),CON,'ShowText','on')
% hold on
% plot(x,x,'k--')
% set(gca,'fontsize',18)
% ylabel('Air Flow (ml/s)')
% axis equal
% grid on
% 
% subplot(3,3,5)
% contour(Q,V,squeeze(Dpalv(:,:,2)),CON,'ShowText','on')
% hold on
% plot(x,x,'k--')
% set(gca,'fontsize',18)
% axis equal
% grid on
% 
% subplot(3,3,6)
% contour(Q,V,squeeze(Ddpac(:,:,2)),CON,'ShowText','on')
% hold on
% plot(x,x,'k--')
% set(gca,'fontsize',18)
% axis equal
% grid on
% 
% subplot(3,3,7)
% contour(Q,V,squeeze(Dpvasc(:,:,3)),CON,'ShowText','on')
% hold on
% plot(x,x,'k--')
% set(gca,'fontsize',18)
% xlabel('Blood Flow (ml/s)')
% ylabel('Air Flow (ml/s)')
% axis equal
% grid on
% 
% subplot(3,3,8)
% contour(Q,V,squeeze(Dpalv(:,:,3)),CON,'ShowText','on')
% hold on
% plot(x,x,'k--')
% set(gca,'fontsize',18)
% xlabel('Blood Flow (ml/s)')
% axis equal
% grid on
% 
% subplot(3,3,9)
% contour(Q,V,squeeze(Ddpac(:,:,3)),CON,'ShowText','on')
% hold on
% plot(x,x,'k--')
% set(gca,'fontsize',18)
% xlabel('Blood Flow (ml/s)')
% axis equal
% grid on

%%% making gifs
CON = -30:5:150;
x = 0:10;
h = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(Dv) %iterate by size of Dv
    %%%% make contour plots
    cla(gca)
    
    contour(Q,V,squeeze(Dpvasc(:,:,i)),CON,'ShowText','on')
    hold on
    plot(x,x,'k--')
    set(gca,'fontsize',18)
    xlabel('Cardiac Output (ml/s)')
    ylabel('Ventilation Magnitude (ml/s)')
    text(9,1,['D = ', num2str(Dv(i))],'fontsize',18)
    axis equal
    grid on
    
    if i == 1
        gif('ModelC_vasc.gif','DelayTime',1/15) %make gif file
    else
        gif %append frame to gif
    end
    disp(i)
end

h2 = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(Dv) %iterate by size of Dv
    %%%% make contour plots
    cla(gca)
    
    contour(Q,V,squeeze(Dpalv(:,:,i)),CON,'ShowText','on')
    hold on
    plot(x,x,'k--')
    set(gca,'fontsize',18)
    xlabel('Cardiac Output (ml/s)')
    ylabel('Ventilation Magnitude (ml/s)')
    text(9,1,['D = ', num2str(Dv(i))],'fontsize',18)
    axis equal
    grid on
    
    if i == 1
        gif('ModelC_alv.gif','DelayTime',1/15) %make gif file
    else
        gif %append frame to gif
    end
    disp(i)
end

h3 = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(Dv) %iterate by size of Dv
    %%%% make contour plots
    cla(gca)
    
    contour(Q,V,squeeze(Ddpac(:,:,i)),CON,'ShowText','on')
    hold on
    plot(x,x,'k--')
    set(gca,'fontsize',18)
    xlabel('Cardiac Output (ml/s)')
    ylabel('Ventilation Magnitude (ml/s)')
    text(9,1,['D = ', num2str(Dv(i))],'fontsize',18)
    axis equal
    grid on
    
    if i == 1
        gif('ModelC_grad.gif','DelayTime',1/15) %make gif file
    else
        gif %append frame to gif
    end
    disp(i)
end
