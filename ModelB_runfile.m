clear; close all; clc;

%%% parameters
D     = 30;      %apparent diffusion (ml/s)
Pair  = 150;    %atmospheric oxygen partial pressure (mmHg)
Pin   = 75;     %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)
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

%%% look up table
load('Lookup.mat') %outputs LOOK
P = LOOK.Plookup;
C = LOOK.Clookup;

%%% solve system of ODEs
X0 = [interp1(P,C,Pin) Pin];
[t,X] = ode15s(@ModelB_RHS,[0 10],X0,[],par,LOOK);

Cvasc = X(:,1);
Pvasc = interp1(C,P,Cvasc);
Palv  = X(:,2);

PSIM = fsolve(@ModelB_FixedPoint_Objective2,[100 100],[],par);
pvasc_sim = PSIM(1);
palv_sim  = PSIM(2);

%%% grid for contour plot
de = 0.1;
q  = 0:de:10;
v  = 0:de:10;
np = length(q);
[Q,V] = meshgrid(q,v);

% computing contour plot
pvasc = zeros(np);
palv  = zeros(np);
FLAG = zeros(np);
tic;
for i = 1:np
   for j = 1:np
       clear tpar P
       tpar = par; tpar(11) = v(i); tpar(12) = q(j);
       [P, ~, FLAG(i,j)] = fsolve(@ModelB_FixedPoint_Objective2,[100 100],[],tpar);
       pvasc(i,j) = P(1); palv(i,j)  = P(2);
       
       disp([i j])
   end
end
toc;

dpac  = palv - pvasc;


%%% plots
figure;
plot(t,Pvasc,t,Palv,'linewidth',2)
xlabel('Time (s)')
ylabel('Oxygen Tension (mmHg)')
legend('Vasc', 'Alv')

figure;
surf(Q,V,pvasc, 'EdgeColor','none')
title('Pvasc')
xlabel('q')
ylabel('v')
axis equal

figure;
surf(Q,V,palv, 'EdgeColor','none')
title('Palv')
xlabel('q')
ylabel('v')
axis equal

figure;
surf(Q,V,dpac, 'EdgeColor','none')
title('\Delta P')
xlabel('q')
ylabel('v')
axis equal



CON = -30:5:150;
x = 0:10;
figure;
subplot(1,3,1)
contour(Q,V,pvasc,CON,'ShowText','on')
hold on
plot(x,x,'k--')
hold on
plot(5,5,'rx')
set(gca,'fontsize',18)
title('Vascular (mmHg)')
xlabel('Cardiac Output (ml/s)')
ylabel('Ventilation Magnitude (ml/s)')
axis equal
grid on

% figure;
subplot(1,3,2)
contour(Q,V,palv,CON,'ShowText','on')
hold on
plot(x,x,'k--')
hold on
plot(5,5,'rx')
set(gca,'fontsize',18)
title('Alveolar Space (mmHg)')
xlabel('Cardiac Output (ml/s)')
axis equal
grid on

% figure;
subplot(1,3,3)
contour(Q,V,dpac,CON,'ShowText','on')
hold on
plot(x,x,'k--')
hold on
plot(5,5,'rx')
set(gca,'fontsize',18)
title('Gradient (mmHg)')
xlabel('Cardiac Output (ml/s)')
axis equal
grid on


