clear; close all; clc;
% This script optimizes the apparent diffusion parameter (D) and maximum O2 consumption (VO2M) to fit the
% CPET data processes in the CO_estimator script

%%% load data
load('ExerciseData','DATA')

%%% parameters
%%% parameters
D     = 500;      %apparent diffusion (L/min)
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
l     = 1; %length of capillary

Vp = 70; %ventilation flow (L/min)
Qp = 70; %blood flow (L/min)
VO2M = 8; %O2 consumption rate (L/min)

k = 2.5;%conversion factor to convert O2 content (ml O2/100ml blood) to concentration (mM)

% empirical function of Cardiac Ouput as a function of work rate - from Stringer 1997
% Reference --> Stringer WW, Hansen JE, Wasserman K. Cardiac output estimated noninvasively from oxygen uptake during exercise. Journal of Applied Physiology. 1997 Mar 1;82(3):908-12.
a = 0.105;
b = 5.72;

% load('ModelD2_optimization_results.mat','D_opt')

par = [D Pair Pin alpha beta l a b k VO2M];

%%% look up table
load('Lookup.mat') %outputs LOOK
Plook = LOOK.Plookup;
Clook = LOOK.Clookup;

N = 1000; %number of grid points for numerical discretization

%%% testing that this janky cost function works
% [JT, JV, JQ, vt, qt] = OBJECTIVE_F_v2(D,par,DATA,Plook,Clook,N);

%%% proper optimization
% opt = optimset('MaxFUnEvals',2000,'Display','iter');
% tic;
% [D_opt,J,EXITFLAG,OUTPUT] = fminsearch(@OBJECTIVE_F_v2, D,opt,par,DATA,Plook,Clook,N);
% RT = toc;
% disp(['Optimizaiton runtime: ',num2str(RT),' s'])

%%% brute force optimization
DF = linspace(450,1000,100);
JF = zeros(length(DF),1);

par = [D Pair Pin alpha beta l a b k VO2M];
N = 1000;
for i = 1:length(DF)
    JF(i) = OBJECTIVE_F_v2(DF(i),par,DATA,Plook,Clook,N);
    disp(i)
end
figure;
loglog(DF,JF,'linewidth',2)
title('Model F Cost')
ylabel('J')
xlabel('D')
set(gca,'fontsize',18)

save('ModelF_optimization_v2_results.mat')





