clear; close all; clc;
% This script runs Model D with a variety of inputs to test it's numerical
% stability.

%%% parameters
D     = 15;      %apparent diffusion (L/min)
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

Vp = 5; %ventilation flow (L/min)
Qp = 5; %blood flow (L/min)
VO2M = 4; %O2 consumption rate (L/min)

k = 2.5;%conversion factor to convert O2 content (ml O2/100ml blood) to concentration (mM)

% empirical function of Cardiac Ouput as a function of work rate - from Stringer 1997
% Reference --> Stringer WW, Hansen JE, Wasserman K. Cardiac output estimated noninvasively from oxygen uptake during exercise. Journal of Applied Physiology. 1997 Mar 1;82(3):908-12.
a = 0.105;
b = 5.72;

% load('ModelD2_optimization_results.mat','D_opt')

par = [D Pair Pin alpha beta l a b k VO2M];

%%% look up table
load('Lookup.mat') %outputs LOOK
P = LOOK.Plookup;
C = LOOK.Clookup;

%%% Relaxing the system to numerically find the fixed points
% N = [1e2 1e3 1e4 1e5 1e6];  %number of discrete spatial nodes
N = [50 100 200 300 400 500 600];  %number of discrete spatial nodes
nv = length(N);

Pv{nv}   = [];
Pa{nv}   = [];
Cvi{nv}  = [];
Pvi{nv}  = [];
EPSv{nv} = [];
dx{nv}   = [];

%%% testing the effect of varing the granularity of the spatial discretization
for i = 1:nv
    tic;
    [Pv{i},Pa{i},Cvi{i},Pvi{i},EPSv{i},dx{i}] = modelD2_SS_relaxation(N(i),par,P,C,Vp,Qp);
    toc;
end


%%% plots and stuffs
figure;
for i = 1:nv
    clear x temp
    x  = dx{i}:dx{i}:l;
    hold on
    plot(x',Pvi{i})
end
legend(num2str(N'))

figure;
for i = 1:nv
    hold on
    plot(dx{i}:dx{i}:l,Cvi{i})
end
set(gca,'fontsize',18)
xlabel('Capillary Length (mm)')
ylabel('O_2 Concnetration (mM)')
legend(num2str(N'))

figure;
for i = 1:nv
    hold on
    plot(dx{i}:dx{i}:l,Pvi{i})
end
set(gca,'fontsize',18)
xlabel('Capillary Length (mm)')
ylabel('O_2 Tension (mmHg)')
legend(num2str(N'))

figure;
subplot(1,2,1)
for i = 1:nv
    hold on
    semilogy(EPSv{i},'linewidth',2)
end
set(gca,'fontsize',18,'yscale','log')
xlabel('Iteration')
ylabel('Error')
legend(num2str(N'))
ylim([1e-12 1e-2])

Iter = zeros(1,nv);
for i = 1:nv
    Iter(i) = length(EPSv{i});
end
% figure;
subplot(1,2,2)
plot(N,Iter)
xlabel('N - number of grid points')
ylabel('Iterations for convergence')
set(gca,'fontsize',18)

