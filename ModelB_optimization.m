clear; close all; clc;
% This script optimizes the apparent diffusion parameter (D) to fit the
% CPET data processes in the CO_estimator script

%%% load data
load('ExerciseData','DATA')

%%% parameters
D     = 30;      %apparent diffusion (ml/s)
Pair  = 150;    %atmospheric oxygen partial pressure (mmHg)
Pin   = 45;     %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)
Vvasc = 1;      %volume of vascular space (ml)
Valv  = 1;      %alveolar volume (ml)

alpha = 1.3e-6*1e3;  % O2 solubility  in water/plasma(mM/mmHg)
CHb   = 0.021*1e3;   % Hb binding site conc (mmol/L of RBC's)
Hct   = 0.40;    % hematocrit (unitless)
C0    = CHb*Hct; % blood oxygen binding capacity (mol/L)
n     = 2.7;     % Hill exponent
P50   = 30;      % half-max saturation of Hb
beta  = 16800*1e-3; % O2 solubility in air (mmHg/mM)

Vp = 5; %ventilation flow (ml/s)
Qp = 5; %blood flow (ml/s)

par = [D Pair Pin Vvasc Valv alpha beta C0 n P50 Vp Qp];

%%% optimization
N = 300; %number of grid points for numerical discretization
lb = 1;
ub = 1000;

opt = optimset('MaxFUnEvals',2000);
tic;
[D_opt,J,EXITFLAG,OUTPUT] = fminsearch(@OBJECTIVE_B, D,opt,par,DATA);
RT = toc;
disp(['Optimizaiton runtime: ',num2str(RT),' s'])


[~,PO2_pred] = OBJECTIVE_B(D_opt,par,DATA);

%%% save and load plots
save('ModelB_optimization_results.mat')
load('ModelB_results.mat','Q','V','pvasc')

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

