%% load data from West and test against our model
% Data from Table 1, West 1962: Regional differences in gas exchange in the lung of erect man	https://journals.physiology.org/doi/epdf/10.1152/jappl.1962.17.6.893

% prep data
HbLookUp = load('Lookup.mat'); %outputs Hb dissociation curve lookup table
HbDisP = HbLookUp.LOOK.Plookup;
HbDisC = HbLookUp.LOOK.Clookup;

westData = readtable('westData');
westData.cO2 = interp1(HbDisP,HbDisC,westData.PO2 , "linear"); % Pulmonary venous distributed sum

Vols = westData.Vol;
qrs = westData.Q;
vrs = westData.Va;

%% parameters
D     = 285;      %apparent diffusion (L/min)
Pair  = 134;    %atmospheric oxygen partial pressure (mmHg) - 164?
Pin   = 45;     %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)
Vvasc = 1;      %volume of vascular space (ml)
Valv  = 1;      %alveolar volume (ml)
l     = 1; %length of capillary
alpha = 1.3e-6*1e3;  % O2 solubility  in water/plasma(mM/mmHg)
beta  = 16800*1e-3; % O2 solubility in air (mmHg/mM)

%load optimized diffusion (D) parameter
load('ModelD_optimization_v2_results.mat','JD','DD')
[~, jDpi] = min(JD); DDp = DD(jDpi);D = DDp;
par = [D Pair Pin alpha beta l];
%
figure(1);clf;

% tune down the diffusion
D = 385;

[pv, cv, p, c, validIds] = calculateDistributedAlveoliD(par, vrs, qrs, Vols, D.*Vols/100, 1);
%
% subplot(231);hold on;
% plot(qrs, Vols, 'o');
% % xlabel('Perfusion blood');
% % ylabel('Chunk size');
% 
% subplot(232);hold on;
% plot(qrs, vrs, 'o-');
% % xlabel('Perfusion blood');
% % ylabel('Perfusion air');
% 
% 
% subplot(233);hold on;
% plot(vrs, Vols, 'o');
% % xlabel('Perfusion air');
% % ylabel('Chunk size');

% distribution of concentrations
subplot(223);hold on;
plot(qrs, westData.cO2, 'o-');
xl = xlim();
cmart = interp1(HbDisP,HbDisC, 97, "linear"); % c mixed arterial from West's pO2
plot(xl, [cmart cmart], 'm--');
legend('Capillary cO2', '1 comp cO2', 'dist venous cO2', 'data West per chunk', 'West total', 'Location','southwest')

% distribution of partial pressure
subplot(224);hold on;
plot(qrs, westData.PO2, 'o-');
xl = xlim();
plot(xl, [97 97], 'm--');
legend('Capillary pO2', '1 comp pO2', 'dist venous pO2', 'data West p02', 'West total pO2', 'Location','northeast')
% PLot over our optimal VQ matching - this moght not really work!
% we compare exercise matching with gravitational matching

% Best fit Exponential Curve Fit (exp2)
a    = 2.3063; 
b    = 0.1459;  
c    = 0.0032;   
d    = 0.5296;
VfunCo = @(CO) a*exp(b*CO) + c*exp(d*CO);

% scale to full size
Qs = qrs./Vols*100;
Vs = vrs./Vols*100;

% scale back
vrsMatch = VfunCo(Qs).*Vols/100;

subplot(232);hold on;
plot(qrs, vrsMatch, 'o-');
xlabel('Perfusion blood');
ylabel('Perfusion air');
legend('West', 'Optimal VQ')

%% plot

figure(2);
chs = 1:numel(Vols);
bar(chs, Vols/100);
xlabel('# Chunk')
figure(3);
perf = qrs./Vols*100;
bar(chs, perf)
xlabel('# Chunk')