% clear; close all; clc;
% this script loads the result of optimization and plots all of the
% residuals

load('ModelA_optimization_results.mat','DATA','PO2_pred')
rA = DATA.PO2 - PO2_pred;
clear PO2_pred
load('ModelB_optimization_results.mat','PO2_pred')
rB = DATA.PO2 - PO2_pred;
clear PO2_pred
load('ModelC_optimization_results.mat','PO2_pred')
rC = DATA.PO2 - PO2_pred;
clear PO2_pred
load('ModelD_optimization_results.mat','PO2_pred')
rD = DATA.PO2 - PO2_pred;
clear PO2_pred
load('ModelB2_optimization_results.mat','PO2_pred')
rE = DATA.PO2 - PO2_pred;
clear PO2_pred
load('ModelD2_optimization_results.mat','PO2_pred')
rF = DATA.PO2 - PO2_pred;
clear PO2_pred

MS = 15; %markersize
lw = 1.5; %linewidth

plot_CO = DATA.CO;

x = 7:20;
figure;
plot(plot_CO,rA,'o-','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rB,'x:','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rC,'+--','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rD,'*-','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rE,'s:','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rF,'d--','markersize',MS,'linewidth',lw)
hold on
plot(x,0*x,'k','linewidth',2)
xlabel('Cardiac Output (L/min)')
ylabel('Residuals (mmHg)')
legend('A','B','C','D','E','F')
set(gca,'fontsize',18)
xlim([8 18])

figure;
subplot(1,2,1)

plot(plot_CO,rA,'o-','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rB,'x:','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rC,'+--','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rD,'*-','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rE,'s:','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rF,'d--','markersize',MS,'linewidth',lw)
hold on
plot(x,0*x,'k','linewidth',2)
xlabel('Cardiac Output (L/min)')
ylabel('Residuals (mmHg)')
legend('A','B','C','D','E','F')
set(gca,'fontsize',18)
xlim([8 18])


subplot(1,2,2)
FLAG = 1;
IND = [5 7];


if FLAG == 1
    rA(IND) = []; rB(IND) = []; rC(IND) = [];
    rD(IND) = []; rE(IND) = []; rF(IND) = [];
    plot_CO(IND) = [];
end


plot(plot_CO,rA,'o-','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rB,'x:','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rC,'+--','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rD,'*-','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rE,'s:','markersize',MS,'linewidth',lw)
hold on
plot(plot_CO,rF,'d--','markersize',MS,'linewidth',lw)
hold on
plot(x,0*x,'k','linewidth',2)
xlabel('Cardiac Output (L/min)')
ylabel('Residuals (mmHg)')
legend('A','B','C','D','E','F')
set(gca,'fontsize',18)
xlim([8 18])


