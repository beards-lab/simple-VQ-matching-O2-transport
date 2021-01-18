clear; close all; clc;

%%% parameters
D     = 50;      %apparent diffusion (ml/s)
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

N = 300; %number of grid points for numerical discretization

%%% grid for contour plot
de = 0.1;
q  = 0:de:10;
v  = 0:de:10;
[Q,V] = meshgrid(q,v);
np = length(q); %number of points for air and blood flow to use

%%% solve for fixed points
pvasc = zeros(np); %preallocate matricies to store fixed points
palv = zeros(np);
parfor i = 1:np
   for j = 1:np
       [pvasc(i,j),palv(i,j)] = modelD_SS_relaxation(N,par,P,C,v(i),q(j)); %loop through q and v and store fixed points
       disp([i j])
   end
end
dpac = palv - pvasc;

save('ModelD_results.mat')

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
x = 0:10;
figure; % Contour Plot
subplot(1,3,1)
contour(Q,V,pvasc,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Vascular (mmHg)')
xlabel('Blood Flow (ml/s)')
ylabel('Air Flow (ml/s)')
axis equal
grid on

subplot(1,3,2)
contour(Q,V,palv,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Alveolar Space (mmHg)')
xlabel('Blood Flow (ml/s)')
axis equal
grid on

subplot(1,3,3)
contour(Q,V,dpac,CON,'ShowText','on')
hold on
plot(x,x,'k--')
set(gca,'fontsize',18)
title('Alv-Vasc Gradient (mmHg)')
xlabel('Blood Flow (ml/s)')
axis equal
grid on

