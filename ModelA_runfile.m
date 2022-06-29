clear; close all; clc;
% This script solves Model A and saves results into a .mat file + generates
% fun animated .gifs to show how the diffusion parameter affects model
% behavior

%%% parameters
D     = 500;  %apparent diffusion (L/min)
Pair  = 150; %atmospheric oxygen partial pressure (mmHg)
Pin   = 45;  %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)
Vvasc = 1;   %volume of vascular space (ml)
Valv  = 1;   %alveolar volume (ml)

Vp = 5; %ventilation flow (L/min)
Qp = 5; %blood flow (L/min)


%load optimized diffusion (D) parameter
load('ModelA_optimization_v2_results.mat','JA','DA')
del = 30;% index to delete from - the first del-many indices are ~1e30 and screw up the minimum finding below
JA(1:del) = []; DA(1:del) = [];
[~, jApi] = min(JA);DAp = DA(jApi);
D = DAp;


%%% parameter vector
par = [D Pair Pin Vvasc Valv Vp  Qp]; %parameters vector for ODE system

%%% grids for contour plots
de = 0.5;
q  = 0:de:70;
v  = 0:de:70;
[Q,V] = meshgrid(q,v);

Dv = 1:100; %vector of diffusion parameters to explore
Dpvasc = zeros(length(q),length(v),length(Dv));
Dpalv  = zeros(length(q),length(v),length(Dv));
Ddpac  = zeros(length(q),length(v),length(Dv));

% computing contour plots
pvasc = (Pin*Q.*V+Pin*D*Q+Pair*D*V)./(Q.*V+D*Q+D*V);
palv  = (Pair*Q.*V+Pair*D*V+Pin*D*Q)./(Q.*V+D*Q+D*V);
dpac  = (Q.*V)*(Pair-Pin)./(Q.*V+D*Q+D*V);

% contour plots for vector of diffusion parameters to explore
for j = 1:length(Dv)
    Dpvasc(:,:,j) = (Pin*Q.*V+Pin*Dv(j)*Q+Pair*Dv(j)*V)./(Q.*V+Dv(j)*Q+Dv(j)*V);
    Dpalv(:,:,j)  = (Pair*Q.*V+Pair*Dv(j)*V+Pin*Dv(j)*Q)./(Q.*V+Dv(j)*Q+Dv(j)*V);
    Ddpac(:,:,j)  = Dpalv(:,:,j) - Dpvasc(:,:,j);
end

%%% solve system of ODEs
% X0 = [Pin Pair];
% [t,X] = ode45(@ModelA_RHS,[0 1],X0,[],par);





save('ModelA_results.mat')

%%% plots
% figure; %state time series
% plot(t,X)
% xlabel('Time (s)')
% ylabel('Oxygen Tension (mmHg)')
% legend('Vasc', 'Alv')

figure; %Pvasc surface
surf(Q,V,pvasc, 'EdgeColor','none')
title('Pvasc')
xlabel('q')
ylabel('v')
axis equal

figure; %Pvasc surface
surf(Q,V,palv, 'EdgeColor','none')
title('Palv')
xlabel('q')
ylabel('v')
axis equal

figure; %Alv-Cap O2 Gradient  surface
surf(Q,V,dpac, 'EdgeColor','none')
title('\Delta P')
xlabel('q')
ylabel('v')
axis equal

CON = -30:5:150;
x = 0:70;
figure; % Contour Plot
subplot(1,3,1)
contour(Q,V,pvasc,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Vascular (mmHg)')
xlabel('Blood Flow (L/min)')
ylabel('Air Flow (L/min)')
axis equal
grid on

subplot(1,3,2)
contour(Q,V,palv,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Alveolar Space (mmHg)')
xlabel('Blood Flow (L/min)')
axis equal
grid on

subplot(1,3,3)
contour(Q,V,dpac,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Alv-Vasc Gradient (mmHg)')
xlabel('Blood Flow (L/min)')
axis equal
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% big old plots
figure; % Contour Plot
subplot(3,3,1)
contour(Q,V,squeeze(Dpvasc(:,:,1)),CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Vascular (mmHg)')
ylabel('Air Flow (ml/s)')
axis equal
grid on

subplot(3,3,2)
contour(Q,V,squeeze(Dpalv(:,:,1)),CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Alveolar Space (mmHg)')
axis equal
grid on

subplot(3,3,3)
contour(Q,V,squeeze(Ddpac(:,:,1)),CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Alv-Vasc Gradient (mmHg)')
axis equal
grid on

subplot(3,3,4)
contour(Q,V,squeeze(Dpvasc(:,:,2)),CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
ylabel('Air Flow (ml/s)')
axis equal
grid on

subplot(3,3,5)
contour(Q,V,squeeze(Dpalv(:,:,2)),CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
axis equal
grid on

subplot(3,3,6)
contour(Q,V,squeeze(Ddpac(:,:,2)),CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
axis equal
grid on

subplot(3,3,7)
contour(Q,V,squeeze(Dpvasc(:,:,3)),CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
xlabel('Blood Flow (ml/s)')
ylabel('Air Flow (ml/s)')
axis equal
grid on

subplot(3,3,8)
contour(Q,V,squeeze(Dpalv(:,:,3)),CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
xlabel('Blood Flow (ml/s)')
axis equal
grid on

subplot(3,3,9)
contour(Q,V,squeeze(Ddpac(:,:,3)),CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
xlabel('Blood Flow (ml/s)')
axis equal
grid on

% %%% making gifs
% CON = -30:5:150;
% h = figure('units','normalized','outerposition',[0 0 1 1]);
% for i = 1:length(Dv) %iterate by size of Dv
%     %%%% make contour plots
%     cla(gca)
%     
%     contour(Q,V,squeeze(Dpvasc(:,:,i)),CON,'ShowText','on')
%     hold on
%     plot(x,x,'k--')
%     set(gca,'fontsize',18)
%     xlabel('Cardiac Output (L/min)')
%     ylabel('Ventilation Magnitude (L/min)')
%     text(9,1,['D = ', num2str(Dv(i))],'fontsize',18)
%     axis equal
%     grid on
%     
%     if i == 1
%         gif('ModelA_vasc.gif','DelayTime',1/15) %make gif file
%     else
%         gif %append frame to gif
%     end
%     disp(i)
% end
% 
% h2 = figure('units','normalized','outerposition',[0 0 1 1]);
% for i = 1:length(Dv) %iterate by size of Dv
%     %%%% make contour plots
%     cla(gca)
%     
%     contour(Q,V,squeeze(Dpalv(:,:,i)),CON,'ShowText','on')
%     hold on
%     plot(x,x,'k--')
%     set(gca,'fontsize',18)
%     xlabel('Cardiac Output (L/min)')
%     ylabel('Ventilation Magnitude (L/min)')
%     text(9,1,['D = ', num2str(Dv(i))],'fontsize',18)
%     axis equal
%     grid on
%     
%     if i == 1
%         gif('ModelA_alv.gif','DelayTime',1/15) %make gif file
%     else
%         gif %append frame to gif
%     end
%     disp(i)
% end
% 
% h3 = figure('units','normalized','outerposition',[0 0 1 1]);
% for i = 1:length(Dv) %iterate by size of Dv
%     %%%% make contour plots
%     cla(gca)
%     
%     contour(Q,V,squeeze(Ddpac(:,:,i)),CON,'ShowText','on')
%     hold on
%     plot(x,x,'k--')
%     set(gca,'fontsize',18)
%     xlabel('Cardiac Output (L/min)')
%     ylabel('Ventilation Magnitude (L/min)')
%     text(9,1,['D = ', num2str(Dv(i))],'fontsize',18)
%     axis equal
%     grid on
%     
%     if i == 1
%         gif('ModelA_grad.gif','DelayTime',1/15) %make gif file
%     else
%         gif %append frame to gif
%     end
%     disp(i)
% end
