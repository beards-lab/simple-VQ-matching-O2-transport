clear; close all; clc;
% This script optimizes the apparent diffusion parameter (D) to fit the
% CPET data processes in the CO_estimator script

%%% load data
load('ExerciseData','DATA')

%%% parameters
D     = 1000;      %apparent diffusion (ml/s)
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

par = [D Pair Pin alpha beta l];

%%% look up table
load('Lookup.mat') %outputs LOOK
Plook = LOOK.Plookup;
Clook = LOOK.Clookup;

N = 300; %number of spatial nodes for discretization

%%% testing this janky cost function
[JT, JV, JQ, vt, qt] = OBJECTIVE_D_v2(D,par,DATA,Plook,Clook,N);

%%% actual optimization
opt = optimset('MaxFUnEvals',2000);
[D_opt,J] = fminsearch(@OBJECTIVE_D_v2,D,opt,par,DATA,Plook,Clook,N);

%%% brute force optimization
DD = linspace(100,1000,200);
JD = zeros(length(DD),1);

par = [D Pair Pin alpha beta l];
N = 300;
for i = 1:length(DD)
    JD(i) = OBJECTIVE_D_v2(DD(i),par,DATA,Plook,Clook,N);
    disp(i)
end
figure;
loglog(DD,JD,'linewidth',2)
title('Model D Cost')
ylabel('J')
xlabel('D')
set(gca,'fontsize',18)

save('ModelD_optimization_v2_results.mat')