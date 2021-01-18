clear; close all; clc;

%%% parameters
D     = 130;      %apparent diffusion (ml/s)
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

Vp = 5; %ventilation flow (ml/s)
Qp = 5; %blood flow (ml/s)

par = [D Pair Pin alpha beta l];

%%% look up table
load('Lookup.mat') %outputs LOOK
P = LOOK.Plookup;
C = LOOK.Clookup;

%%% Relaxing the system to numerically find the fixed points
N = [100 200 500 1000 2000 5000 1e4];  %number of discrete spatial nodes
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
    [Pv{i},Pa{i},Cvi{i},Pvi{i},EPSv{i},dx{i}] = modelD_SS_relaxation(N(i),par,P,C,Vp,Qp);
    toc;
end

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
for i = 1:nv
    hold on
    semilogy(EPSv{i},'linewidth',2)
end
set(gca,'fontsize',18,'yscale','log')
xlabel('Iteration')
ylabel('Error')
legend(num2str(N'))

Iter = zeros(1,nv);
for i = 1:nv
    Iter(i) = length(EPSv{i});
end
figure;
plot(N,Iter)
xlabel('N - number of grid points')
ylabel('Iterations for convergence')
set(gca,'fontsize',18)

