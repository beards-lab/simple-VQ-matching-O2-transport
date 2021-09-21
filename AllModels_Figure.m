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
load('ModelB2_results.mat')
QE = Q; VE = V; palvE = palv; pvascE = pvasc; dpacE = dpac;
clear Q V palv pvasc dpac
load('ModelD2_results.mat')
QF = Q; VF = V; palvF = palv; pvascF = pvasc; dpacF = dpac;
clear Q V palv pvasc dpac


CON = 10:10:150;
x = 0:70;

xl = 40;
yl = 70;
lw = 1.5;

%%% master contour plot
figure;
subplot(3,6,1)
contour(QA,VA,palvA,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
title('Model A')
ylabel('Air Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

subplot(3,6,7)
contour(QA,VA,pvascA,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
ylabel('Air Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

subplot(3,6,13)
contour(QA,VA,dpacA,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
ylabel('Air Flow (L/min)')
xlabel('Blood Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

subplot(3,6,2)
contour(QB,VB,palvB,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
title('Model B')
xlim([0 xl])
ylim([0 yl])

subplot(3,6,8)
contour(QB,VB,pvascB,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlim([0 xl])
ylim([0 yl])

subplot(3,6,14)
contour(QB,VB,dpacB,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlabel('Blood Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

subplot(3,6,3)
contour(QC,VC,palvC,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
title('Model C')
xlim([0 xl])
ylim([0 yl])

subplot(3,6,9)
contour(QC,VC,pvascC,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlim([0 xl])
ylim([0 yl])

subplot(3,6,15)
contour(QC,VC,dpacC,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlabel('Blood Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

subplot(3,6,4)
contour(QD,VD,palvD,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
title('Model D')
xlim([0 xl])
ylim([0 yl])

subplot(3,6,10)
contour(QD,VD,pvascD,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlim([0 xl])
ylim([0 yl])

subplot(3,6,16)
contour(QD,VD,dpacD,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlabel('Blood Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

subplot(3,6,5)
contour(QE,VE,palvE,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
title('Model E')
xlim([0 xl])
ylim([0 yl])

subplot(3,6,11)
contour(QE,VE,pvascE,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlim([0 xl])
ylim([0 yl])

subplot(3,6,17)
contour(QE,VE,dpacE,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlabel('Blood Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

subplot(3,6,6)
contour(QF,VF,palvF,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
title('Model F')
xlim([0 xl])
ylim([0 yl])

subplot(3,6,12)
contour(QF,VF,pvascF,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlim([0 xl])
ylim([0 yl])

subplot(3,6,18)
contour(QF,VF,dpacF,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlabel('Blood Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Models A-D
figure;
subplot(3,4,1)
contour(QA,VA,palvA,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
title('Model A')
ylabel('Air Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

subplot(3,4,5)
contour(QA,VA,pvascA,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
ylabel('Air Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

subplot(3,4,9)
contour(QA,VA,dpacA,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
ylabel('Air Flow (L/min)')
xlabel('Blood Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

subplot(3,4,2)
contour(QB,VB,palvB,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
title('Model B')
xlim([0 xl])
ylim([0 yl])

subplot(3,4,6)
contour(QB,VB,pvascB,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlim([0 xl])
ylim([0 yl])

subplot(3,4,10)
contour(QB,VB,dpacB,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlabel('Blood Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

subplot(3,4,3)
contour(QC,VC,palvC,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
title('Model C')
xlim([0 xl])
ylim([0 yl])

subplot(3,4,7)
contour(QC,VC,pvascC,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlim([0 xl])
ylim([0 yl])

subplot(3,4,11)
contour(QC,VC,dpacC,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlabel('Blood Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

subplot(3,4,4)
contour(QD,VD,palvD,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
title('Model D')
xlim([0 xl])
ylim([0 yl])

subplot(3,4,8)
contour(QD,VD,pvascD,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlim([0 xl])
ylim([0 yl])

subplot(3,4,12)
contour(QD,VD,dpacD,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlabel('Blood Flow (L/min)')
xlim([0 xl])
ylim([0 yl])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Models E-F

pCON = 140:2:150;

figure;
subplot(3,2,1)
contour(QE,VE,palvE,pCON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
title('Model E')
xlim([0 xl])
ylim([0 yl])

subplot(3,2,3)
contour(QE,VE,pvascE,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlim([0 xl])
ylim([0 yl])

subplot(3,2,5)
contour(QE,VE,dpacE,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlabel('Blood Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

subplot(3,2,2)
contour(QF,VF,palvF,pCON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
title('Model F')
xlim([0 xl])
ylim([0 yl])

subplot(3,2,4)
contour(QF,VF,pvascF,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlim([0 xl])
ylim([0 yl])

subplot(3,2,6)
contour(QF,VF,dpacF,CON,'ShowText','on','linewidth',lw)
hold on
plot(x,x,'k--','linewidth',lw)
set(gca,'fontsize',18)
xlabel('Blood Flow (L/min)')
xlim([0 xl])
ylim([0 yl])

