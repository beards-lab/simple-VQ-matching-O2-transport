clear; close all; clc;
% This script solves Model D and saves results into a .mat file

%%% parameters
D     = 100;      %apparent diffusion (L/min)
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

load('ModelD2_optimization_results.mat','D_opt')

par = [D_opt Pair Pin alpha beta l a b k VO2M];

%%% look up table
load('Lookup.mat') %outputs LOOK
P = LOOK.Plookup;
C = LOOK.Clookup;

N = 300; %number of grid points for numerical discretization

%%% grid for contour plot
de = 1;
q  = 0:de:Qp;
v  = 0:de:Vp;
[Q,V] = meshgrid(q,v);
np = length(q); %number of points for air and blood flow to use

%%% solve for fixed points
pvasc = zeros(np); %preallocate matricies to store fixed points
palv = zeros(np);
tic;
parfor i = 1:np
   for j = 1:np
       [pvasc(i,j),palv(i,j)] = modelD2_SS_relaxation(N,par,P,C,v(i),q(j)); %loop through q and v and store fixed points
       disp([i j])
   end
end
T = toc;
disp(['Runtime: ', T, ' s'])
dpac = palv - pvasc;

save('ModelD2_results.mat')

%%% pretty plots and stuff
figure;
surf(Q,V,pvasc, 'EdgeColor','none')
title('Pvasc')
xlabel('q')
ylabel('v')

figure;
surf(Q,V,palv, 'EdgeColor','none')
title('Palv')
xlabel('q')
ylabel('v')

figure;
surf(Q,V,dpac, 'EdgeColor','none')
title('\Delta P')
xlabel('q')
ylabel('v')

CON = -30:5:150;
x = 0:max(v);
figure; % Contour Plot
subplot(1,3,1)
contour(Q,V,palv,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Alveolar Oxygen Tension (mmHg)')
xlabel('Blood Flow (L/min)')
ylabel('Air Flow (L/min)')
axis equal
grid on
xlim([0 40])

subplot(1,3,2)
contour(Q,V,pvasc,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Vascular Oxygen Tension (mmHg)')
xlabel('Blood Flow (L/min)')
axis equal
grid on
xlim([0 40])

subplot(1,3,3)
contour(Q,V,dpac,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Alv-Vasc Oxygen Difference (mmHg)')
xlabel('Blood Flow (L/min)')
axis equal
grid on
xlim([0 40])
