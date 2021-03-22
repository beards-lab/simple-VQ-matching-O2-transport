clear; close all; clc;
% This script optimizes the apparent diffusion parameter (D) to fit the
% CPET data processes in the CO_estimator script

%%% load data
load('ExerciseData','DATA')

%%% parameters
D     = 50;      %apparent diffusion (ml/s)
Pair  = 150;    %atmospheric oxygen partial pressure (mmHg)
Pin   = 45;     %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)

par = [D Pair Pin];

%%% optimization
N = 300; %number of grid points for numerical discretization
lb = 1e4;
ub = 1e12;

tic;
[D_opt,J,EXITFLAG,OUTPUT] = fminsearch(@OBJECTIVE_C, D,[],par,DATA);
RT = toc;
disp(['Optimizaiton runtime: ',num2str(RT),' s'])


[~,PO2_pred] = OBJECTIVE_C(D_opt,par,DATA);

%%% save and load plots
save('ModelC_optimization_results.mat')
load('ModelC_results.mat','Q','V','pvasc')

%%% plots
CON = -30:5:150;
x = 0:35;

figure;
plot(DATA.CO,DATA.V,'ko','linewidth',2,'markersize',10)
hold on
contour(Q,V,pvasc,CON,'ShowText','on')
hold on
plot(x,x,'k--')
xlabel('Cardiac Output (L/min)')
ylabel('Ventilation (L/min)')
set(gca,'fontsize',18)
axis equal

