clear; close all; clc;

%%% gas exchange parameters
alpha   = 1.3e-6*1e3;  % O2 solubility by Henry's law(mM/mmHg)
CHb     = 0.021*1e3;   % Hb binding site conc (mmol/L of RBC's)
Hct     = 0.40;    % hematocrit (unitless)
C0      = CHb*Hct; % blood oxygen binding capacity (mol/L)
n       = 2.7;     % Hill exponent
P50     = 27;      % half-max saturation of Hb
DO      = 1;     % (apparent) O2 diffusion coefficient (ml/s)

% O2 look up table (concetration to partial pressure conversion)
Plookup = 0:1:200; %look up table
Clookup = alpha*Plookup + C0*(Plookup.^n)./(Plookup.^n + P50.^n);

SAT = (Plookup.^n)./(Plookup.^n + P50.^n);
FREE = alpha*Plookup;

LOOK.Plookup = Plookup;
LOOK.Clookup = Clookup;

figure;

plot(Plookup,Clookup,'k','linewidth',2)
hold on
plot(Plookup, FREE,'k:','linewidth',2)
hold on
plot(Plookup,C0*SAT,'k--','linewidth',2)
set(gca,'fontsize',18)
xlabel('Oxygen Partial Pressure (mmHg)')
ylabel('Oxygen Concentration (mM)')
legend('Total','Freely dissolved', 'Hb - bound','location','northwest')


save('Lookup.mat','LOOK')