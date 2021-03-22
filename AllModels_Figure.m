clear; close all; clc;
% this script loads the result of contour plots for each model and merges
% them all together into one common figure

load('ModelA_results.mat','Q','V','palv','pvasc','dpac')
QA = Q; VA = V; palvA = palv; pvascA = pvasc; dpacA = dpac;
clear Q V palv pvasc dpac
load('ModelB_results.mat','Q','V','palv','pvasc','dpac')
QB = Q; VB = V; palvB = palv; pvascB = pvasc; dpacB = dpac;
clear Q V palv pvasc dpac
load('ModelC_results.mat','Q','V','palv','pvasc','dpac')
QC = Q; VC = V; palvC = palv; pvascC = pvasc; dpacC = dpac;
clear Q V palv pvasc dpac
load('ModelD_results.mat','Q','V','palv','pvasc','dpac')
QD = Q; VD = V; palvD = palv; pvascD = pvasc; dpacD = dpac;
clear Q V palv pvasc dpac

CON = 5:5:150;
x = 0:70;

figure; % Contour Plot
subplot(4,3,1)
contour(QA,VA,palvA,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Alveolar Oxygen Tension (mmHg)')
ylabel('Air Flow (ml/s)')
axis equal
% grid on

subplot(4,3,2)
contour(QA,VA,pvascA,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Vascular Oxygen Tension (mmHg)')
axis equal
% grid on

subplot(4,3,3)
contour(QA,VA,dpacA,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Alv-Vasc Difference (mmHg)')
axis equal
% grid on

subplot(4,3,4)
contour(QB,VB,palvB,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
ylabel('Air Flow (ml/s)')
axis equal
% grid on

subplot(4,3,5)
contour(QB,VB,pvascB,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
axis equal
% grid on

subplot(4,3,6)
contour(QB,VB,dpacB,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
axis equal
% grid on

subplot(4,3,7)
contour(QC,VC,palvC,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
ylabel('Air Flow (ml/s)')
axis equal
% grid on

subplot(4,3,8)
contour(QC,VC,pvascC,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
axis equal
% grid on

subplot(4,3,9)
contour(QC,VC,dpacC,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
axis equal
% grid on

subplot(4,3,10)
contour(QD,VD,palvD,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
xlabel('Blood Flow (ml/s)')
ylabel('Air Flow (ml/s)')
axis equal
% grid on

subplot(4,3,11)
contour(QD,VD,pvascD,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
xlabel('Blood Flow (ml/s)')
axis equal
% grid on

subplot(4,3,12)
contour(QD,VD,dpacD,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
xlabel('Blood Flow (ml/s)')
axis equal
% grid on

