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

%%% testing this janky cost function works
% [JT, JV, JQ, vt, qt] = OBJECTIVE_B_v2(10000,par,DATA);

%%% proper optimization
% opt = optimset('MaxFUnEvals',2000);
% tic;
% [D_opt,J,EXITFLAG,OUTPUT] = fminsearch(@OBJECTIVE_B_v2, D,opt,par,DATA);
% RT = toc;
% disp(['Optimizaiton runtime: ',num2str(RT),' s'])

%%% brute force
DB = linspace(1e2,1e5,1000);
JB = zeros(length(DB),1);

par = [D Pair Pin Vvasc Valv alpha beta C0 n P50 Vp Qp];
parfor i = 1:length(DB)
    JB(i) = OBJECTIVE_B_v2(DB(i),par,DATA);
    disp(i)
end
figure;
loglog(DB,JB,'linewidth',2)
title('Model B Cost')
ylabel('J')
xlabel('D')
set(gca,'fontsize',18)

save('ModelB_optimization_v2_results.mat')


