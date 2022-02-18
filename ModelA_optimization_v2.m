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
opt = optimset('MaxFUnEvals',2000);


%%% testing this janky cost function works
[JT, JV, JQ, vt, qt] = OBJECTIVE_A_v2(100,par,DATA);

%%% proper optimization
tic;
[D_opt,J,EXITFLAG,OUTPUT] = fminsearch(@OBJECTIVE_A_v2, 100,opt,par,DATA);
RT = toc;
disp(['Optimizaiton runtime: ',num2str(RT),' s'])

%%% brute force
DA = linspace(0,1000,1000);
JA = zeros(length(DA),1);
for i = 1:length(DA)
    JA(i) = OBJECTIVE_A_v2(DA(i),par,DATA);
end
figure;
loglog(DA,JA,'linewidth',2)
title('Model A Cost')
ylabel('J')
xlabel('D')
set(gca,'fontsize',18)

save('ModelA_optimization_v2_results.mat')

