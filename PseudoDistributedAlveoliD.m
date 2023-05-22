clear;close all;
%% Test Hb dissociation look up table
alpha = 1.3e-6*1e3;  % O2 solubility  in water/plasma(mM/mmHg)
CHb   = 0.021*1e3;   % Hb binding site conc (mmol/L of RBC's)
Hct   = 0.40;    % hematocrit (unitless)
C0    = CHb*Hct; % blood oxygen binding capacity (mol/L)
n     = 2.7;     % Hill exponent
P50   = 27;      % half-max saturation of Hb
beta  = 16800*1e-3; % O2 solubility in air (mmHg/mM)

HbLookUp = load('Lookup.mat'); %outputs Hb dissociation curve lookup table
HbDisP = HbLookUp.LOOK.Plookup;
HbDisC = HbLookUp.LOOK.Clookup;
% Ceq = zeros(size(HbDisP));
% for iqs = 1:length(HbDisP)
%     Ceq(iqs) = alpha*HbDisP(iqs)+C0*((HbDisP(iqs)^n)/(HbDisP(iqs)^n+P50^n)); % transform to Pa
% end
% test interp1
% Cx = 1:1.7:12;
% Pvdist = interp1(HbDisC, HbDisP,Cx , "linear");
% plot(HbDisP, HbDisC, HbDisP, Ceq, '--', Pvdist, Cx, 'x-', 'Linewidth', 2);

%% Prepare pseudo-random continuous distribution
N = 10;
SD = 1;
% SD = 2.5e-1;

% N = 20;
% SD = 1e-1;

in = linspace(0.0, 2, N);
% normal probability distribution density
bell = normpdf(in, 1, SD);

% cut off low probabilities for compuational reasons
trashold = 1e-3;
in = in(bell>trashold);
bell = bell(bell>trashold);
N = numel(bell);

Vdist = in/ mean(in); % ventilation sweep with mean of 1
Pd = bell / mean(bell); % probability distribution with mean of one
% sum

% plot(Qdist, Qpd, Qdist, bell2)
% figure;plot(Vdist, '*');
% plot(Pd, '*');

figure(1);clf;subplot(221);
plot(Vdist, Pd, '*-');title('Probability distribution with mean = 1');
% bar(vrs, Pd);title(sprintf('Distribution of ventilation among %d elements', N));
% xlabel('Ventilated elements (L/)');
% 
% % distribution of flow
% figure(2);clf;
% bar(vrs, qrs);title('Distribution of flow');

%% parameters
D     = 285;      %apparent diffusion (L/min)
Pair  = 150;    %atmospheric oxygen partial pressure (mmHg)
Pin   = 45;     %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)
Vvasc = 1;      %volume of vascular space (ml)
Valv  = 1;      %alveolar volume (ml)
l     = 1; %length of capillary

%load optimized diffusion (D) parameter
load('ModelD_optimization_v2_results.mat','JD','DD')
[~, jDpi] = min(JD); DDp = DD(jDpi);D = DDp;
par = [D Pair Pin alpha beta l];

%% Run ModelD eval - test the distribution
NN = 500;
tic
Qt = 0.05;Vpt = 0.3;
par(6) = 1;
par(1) = D;
[Pv,Pa,Cvi,Pvi,EPSv,dx] = modelD_SS_relaxation(NN,par,HbDisP,HbDisC,5,5);
t = toc;
fprintf('At Q = %2.2f and Vp  = %2.2f, pO2_{dist} = %2.1f (alv. = %2.1f), in %2.0f ms \n', Qt, Vpt, Pv, Pa, t*1000);
% VQset = 0:0.01:0.2;
% Cvid = zeros(numel(VQset), 1);
% for i = 1:numel(VQset)
%     fprintf('Jedu %2.0f %% \n' , i/numel(VQset)*100);
%     [~,~,C] = modelD_SS_relaxation(NN,par,HbDisP,HbDisC,VQset(i),VQset(i));
%     Cvid(i) = C(end);
% end
% %
% cla;hold on;
% plot(VQset, Cvid, 'o-');
% yl = ylim;
% plot(VQset(isnan(Cvid)), yl(1) + (yl(2)-yl(1))*0.9*ones(size(find(isnan(Cvid)))), 'o-');
% t = toc;
% fprintf('Cvi %1.3f, CviSum %1.3f \n', Cvi(end), Csum);



%% init
% inputs
V = 5;
CO = 5;

vrs = V/N * Vdist; % pseudorandom ventilation set
% qrs = CO/N .* ones(1, N); % flow - start simple unmatched
qrs = vrs; % flow matches ventilation 1:1
% Best fit Exponential Curve Fit (exp2)
% a    = 2.3063; b    = 0.1459;  c    = 0.0032;   d    = 0.5296;
% VfunCo = @(CO) a*exp(b*CO) + c*exp(d*CO);
% a = 13.3119; b = 0.0044;c = -18.0414; d = -0.1465;
% COfunV = @(V) a*exp(b*V) + c*exp(d*V);
% qrs = COfunV(vrs*V)/V;

qrs = CO/N * Vdist;
vrs = V/N .* ones(1, N);

% qrs as a function of q
% A = struct();
figure(1)
A = calculateDistributedAlveoliD(par, vrs, qrs, Pd, true);

%% save results
% DMR = struct(); % distributed model results
% DMR.CO = CO; % cardiac output
% DMR.Vp = Vpd; % ventilation
% DMR.Pv = Pvdist; % PO2 in pulmonary veins
% DMR.Qpd = Qpd; % flow probability distribution - same for all data points
% DMR.pVascND = pVascND; % non distributed
% DMR.pAlvND = pAlvND;
% snapfile = sprintf('Data\\snapshot_V%d_Q%d.mat', Vps, Qs);
% DMR.snap = load(snapfile).snap; % paralel for snapshot
% 
% save(filename, "DMR");
%% load results
% clear;
% figure(1);clf;
% DMR = load(filename).DMR; % distributed model results
% 
% % DMR = load("Data/DistributedModel_SD_0.mat").DMR; % distributed model results
% % DMR = load("DistributedModel_SD_2e-1.mat").DMR; % distributed model results
% % DMR = load("DistributedModel_SD_1e-3.mat").DMR; % distributed model results
% % N = length(DMR.snap.Qds)*10;
% 
% % DMR.snap.PvNDs = DMR.pVascND(3, 3);
% % DMR.snap.Pvs = DMR.snap.Ps;
% % subplot(131); histogram(dmr.snap.Qds, N, "EdgeAlpha",0.2);
% % subplot(132); histogram(dmr.snap.Pds, N); 
% %%
% subplot(231); cla; hold on;
% % b = DMR.snap.Qds.*DMR.Qpd;
% bar(DMR.snap.Qds, DMR.Qpd/sum(DMR.Qpd)*100);
% title(sprintf('Flow pseudodistribution at N = %d', N));
% ylabel('Count (%)');xlabel('BLood flow');
% % plot(DMR.snap.Qds,DMR.snap.Vds, 'x');
% % title(['Blood flow distribution in alveoli, Q_{total} = ' num2str(round(sum(DMR.snap.Qds), 2)) ' L/min']);
% % xlabel('Q');ylabel('V');
% % Qsm = mean(DMR.snap.Qds);
% % Vsm = mean(DMR.snap.Vds);
% % plot([Qsm Qsm], [min(DMR.snap.Vds) max(DMR.snap.Vds)], 'r--')
% % plot([min(DMR.snap.Qds) max(DMR.snap.Qds)], [Vsm Vsm], 'r--')
% %%
% subplot(232); cla;hold on;
% % plot(DMR.snap.Pds,1:numel(DMR.snap.Qds), 'x'); 
% b = DMR.snap.Qds.*DMR.Qpd;
% bar(DMR.snap.Vds, DMR.Qpd/sum(DMR.Qpd)*100);
% title(sprintf('Flow pseudodistribution at N = %d', N));
% ylabel('Count (%)');xlabel('BLood flow');
% 
% 
% histogram(DMR.snap.Pds, N/10);
% title(['Partial pressure in distributed alveoli, Pa = ' num2str(round(DMR.snap.Pvs, 2))]);
% plot([DMR.snap.PvNds DMR.snap.PvNds], [0 max(ylim)], 'c:', 'LineWidth',1.5)
% plot([DMR.snap.Pvs DMR.snap.Pvs], [0 max(ylim)], 'r--')
% legend('Capillary pO2', '1 comp pO2', 'dist venous pO2', 'Location','northwest')
% xlabel('P_{O2}');ylabel('Count in category');
% 
% subplot(233); cla;hold on;
% [X Y] = meshgrid(DMR.CO, DMR.Vp);
% % surf(X, Y, DMR.Pv');
% contour(X, Y, DMR.pVascND', '--', 'ShowText','on', 'LineWidth',2);
% contour(X, Y, DMR.Pv', '-', 'ShowText','on', 'LineWidth',1);
% 
% xlabel('CO L/min');
% ylabel('V (L/min)');
% zlabel('P_O2 (mmHg)');
% legend('1 comp', 'Distributed')
% colorbar
% % view(2)
% 
% % PLot surf
% subplot(223);cla;hold on;
% surf(X, Y, DMR.pVascND', 'LineStyle','--', 'LineWidth',1);
% surf(X, Y, DMR.Pv', 'LineStyle','-', 'LineWidth',1);
% xlabel('CO L/min');
% ylabel('V (L/min)');
% zlabel('P_O2 (mmHg)');
% colorbar
% % view(2)
% legend('Distributed', 'Single comp')
% 
% % plot Px
% % figure(3);clf;hold on;
% subplot(224);cla;hold on;
% serco = 1:1:numel(DMR.CO); % selected range CO
% servp = 3:3:numel(DMR.Vp);% selected range Vp
% plot(DMR.CO(serco), DMR.pVascND(serco, servp), 'o-');
% set(gca, 'ColorOrderIndex', 1)
% plot(DMR.CO(serco), DMR.Pv(serco, servp), 'x--');
% xlabel('whole organ CO (L/min)');ylabel('pO2 (mmHg)')
% title('Distributed alveoli with normal distribution of flow at Vp rates')
% legend(string([DMR.Vp(servp) DMR.Vp(servp)]), 'NumColumns',2)
% % legend('Single compartment', 'Distributed')
% % legend('SD = 0.001', 'SD = 0.1','SD = 1');
% saveas(gcf, sprintf("figures/DistributedModelD_VQ_SD_%1.1e.fig", SD));
% saveas(gcf, sprintf("figures/DistributedModelD_VQ_SD_%1.1e.png", SD));