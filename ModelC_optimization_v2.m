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


%%% testing this janky cost function works
% [JT, JV, JQ, vt, qt] = OBJECTIVE_C_v2(100,par,DATA);

%%% proper optimization
tic;
[D_opt,J,EXITFLAG,OUTPUT] = fminsearch(@OBJECTIVE_C_v2, D,[],par,DATA);
RT = toc;
disp(['Optimizaiton runtime: ',num2str(RT),' s'])

%%% brute force optimization
DC = linspace(10,1000,1000);
JC = zeros(length(DC),1);

par = [D Pair Pin];
for i = 1:length(DC)
    JC(i) = OBJECTIVE_C_v2(DC(i),par,DATA);
end
figure;
loglog(DC,JC,'linewidth',2)
title('Model C Cost')
ylabel('J')
xlabel('D')
set(gca,'fontsize',18)

save('ModelC_optimization_v2_results.mat')