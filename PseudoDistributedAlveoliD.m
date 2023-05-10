clear;close all;
%% Constants and Hb dissociation look up table
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
Ceq = zeros(size(HbDisP));
for iqs = 1:length(HbDisP)
    Ceq(iqs) = alpha*HbDisP(iqs)+C0*((HbDisP(iqs)^n)/(HbDisP(iqs)^n+P50^n)); % transform to Pa
end

% test interp1
% Cx = 1:1.7:12;
% Pvdist = interp1(HbDisC, HbDisP,Cx , "linear");
% plot(HbDisP, HbDisC, HbDisP, Ceq, '--', Pvdist, Cx, 'x-', 'Linewidth', 2);

%% Prepare pseudo-random continuous distribution
N = 20;
SD = 2e-1;
in = linspace(0, 2, N);
% normal probability distribution density
bell = normpdf(in, 1, SD);

% cut off low probabilities for compuational reasons
trashold = 0.01;
in = in(bell>trashold);
bell = bell(bell>trashold);
N = numel(bell);
% bell2 = (bell) / sum(bell); % sums to one

Vdist = in/ mean(in); % ventilation sweep with mean of 1
Pd = bell / mean(bell); % probability distribution with mean of one
% sum

% plot(Qdist, Qpd, Qdist, bell2)
figure;plot(Vdist, '*');
plot(Pd, '*');

drawPlots = true;
if drawPlots
    % figure(1);
    % subplot(121);hold on;
    % subplot(132);hold on;
    % subplot(133);hold on;
end

%% Prepare inputs
CO = [5];
V = [5];

% CO = [10];
% Vp = [10];

% CO = 5;
% Vp = 5;
Qs = 5; % flow snapshot
Vps = 5; % Ventilation snapshot

% filename = sprintf("Data/DistributedModel_SD_%1.1e.mat", SD)
filename = sprintf("Data/PseudoDistributedModelQV_SD_%1.1e.mat", SD)

%% get all the param shift tigerther

%%% parameters
D     = 285;      %apparent diffusion (L/min)
Pair  = 150;    %atmospheric oxygen partial pressure (mmHg)
Pin   = 45;     %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)
Vvasc = 1;      %volume of vascular space (ml)
Valv  = 1;      %alveolar volume (ml)
l     = 1; %length of capillary

%load optimized diffusion (D) parameter
load('ModelD_optimization_v2_results.mat','JD','DD')
[~, jDpi] = min(JD); DDp = DD(jDpi);
D = DDp;

par = [D Pair Pin alpha beta l];

%% Run ModelD eval - test the distribution
NN = 500;
% tic
% Qt = 0.05;Vpt = 0.3;
% par(6) = 1;
% par(1) = D;
% [Pv,Pa,Cvi,Pvi,EPSv,dx] = modelD_SS_relaxation(NN,par,HbDisP,HbDisC,Vps,Qs);
% t = toc;
% fprintf('At Q = %2.2f and Vp  = %2.2f, pO2_{dist} = %2.1f (alv. = %2.1f), in %2.0f ms \n', Qt, Vpt, Pv, Pa, t*1000);
% 
% NNN = 3;
% Cvid = zeros(NNN, 1);
% for i = 1:NNN
%     par(1) = D/NNN;
%     [~,~,C] = modelD_SS_relaxation(NN,par,HbDisP,HbDisC,Vps/NNN,Qs/NNN);
%     Cvid(i) = C(end);
% end
% Csum = sum(Cvid*Qs/NNN)/Qs;
% 
% t = toc;
% fprintf('Cvi %1.3f, CviSum %1.3f \n', Cvi(end), Csum);
% 


%% init loop

pVascND = zeros(length(CO), length(V)); % non distributed
pAlvND = zeros(length(CO), length(V)); % non distributed
Pvdist = zeros(length(CO), length(V)); % distributed
Pds = 0;Qds = 0;Ps = 0;
params = par;
pfs = @(snap, V, Q) save(sprintf('Data\\snapshot_V%d_Q%d.mat', V, Q), "snap"); % paralel for snapshot
COEl = 1:length(CO);
VpEl = 1:length(V); % VpElements

for iv = VpEl % iterate ventilation
    par = params;

    disp("Jedu " + num2str(round(iv/length(V)*100)))
    vrs = V(iv)/N * Vdist; % pseudorandom ventilation set
    

    Cvi = 0;Pvi = 0;

    for iq = COEl % iterate flow
        qrs = CO(iq)/N .* ones(N, 1); % flow - start simple unmatched
        % qrs = CO(iq)/N .* Vpd(iv)/N ; % flow matches ventilation 1:1
        % qrs as a function of q
        
        % equal Q and Vp in all alveoli
        par(1) = D;
        [pv1,~, ~, ~] = modelD_SS_relaxation(NN,par,HbDisP,HbDisC,V(iv),CO(iq)); 
        pVascND(iq, iv) = pv1;
        
        tic
        % distribution within lungs
        c = zeros(size(vrs));
        p = zeros(size(vrs));
        for iqs = 1:length(vrs) % iterate submodel flows
            par(1) = D/N;
            [~,~, Cvi, Pvi] = modelD_SS_relaxation(NN,par,HbDisP,HbDisC,vrs(iqs),qrs(iqs)); %loop through q and v and store fixed points
            c(iqs) = Cvi(end); % pulmonary end-capillary concentration
            p(iqs) = Pvi(end); % pulmonary end-capillary pO2
        end
        t = toc;
        Peld = interp1(HbDisC, HbDisP,c , "linear"); % partial pressure per ventilation element
        Cvd = (c.*qrs'.*Pd); % concentration venous distributed, weighted by "probability distribution"
        % Pvd = interp1(HbDisC, HbDisP,Cvd , "linear"); 
        Cv = sum(Cvd)/sum(qrs); % weighted average of concentrations by flow
        % just another mean of calculation of the previous
        % Cv = sum(1/N*c.*qrs'.*Pd)/mean(qrs); 
        pv = interp1(HbDisC, HbDisP,Cv , "linear"); % Pulmonary venous distributed sum
        Pvdist(iq, iv) = pv;
        fprintf('At Q = %2.1f (%2.1f) and Vp  = %2.1f, pO2_{dist} = %2.1f (single comp. = %2.1f), in %2.0f ms (Ctrl %1.3e)\n', sum(qrs), CO(iq), Vpd(iv), Pvdist(iq, iv), pv1, t*1000, CO(iq)-sum(qrs));

        %% plot that
        figure(1);clf;subplot(211);
        bar(vrs, Vdist.*Pd/N);title(sprintf('Distribution of ventilation among %d elements', N));
        xlabel('Ventilated elements (L/)');

        % distribution of flow
        figure(2);clf;
        bar(vrs, qrs);title('Distribution of flow');

        % distribution of concentrations
        % figure(3);clf;
        % bar(vrs, c);title('Distribution of concentration (raw)');

        % distribution of partial pressure
        figure(3);clf;
        bar(vrs, Peld);title('Distribution of partial pressures (raw)');

        % Total partial pressure
        figure(4);clf;hold on;
        bar(vrs, Pvd);title(sprintf('Weighted partial pressures, total: %0.2f', Pvdist(iq, iv)));
        plot(xlim, [pv1 pv1], 'r--')
        plot(xlim, [pv pv], 'c:', 'LineWidth',1.5)
        legend('Capillary pO2', '1 comp pO2', 'dist venous pO2', 'Location','northwest')
        
    end
end
%% save results
DMR = struct(); % distributed model results
DMR.CO = CO; % cardiac output
DMR.Vp = Vpd; % ventilation
DMR.Pv = Pvdist; % PO2 in pulmonary veins
DMR.Qpd = Qpd; % flow probability distribution - same for all data points
DMR.pVascND = pVascND; % non distributed
DMR.pAlvND = pAlvND;
snapfile = sprintf('Data\\snapshot_V%d_Q%d.mat', Vps, Qs);
DMR.snap = load(snapfile).snap; % paralel for snapshot

save(filename, "DMR");
%% load results
% clear;
figure(1);clf;
DMR = load(filename).DMR; % distributed model results

% DMR = load("Data/DistributedModel_SD_0.mat").DMR; % distributed model results
% DMR = load("DistributedModel_SD_2e-1.mat").DMR; % distributed model results
% DMR = load("DistributedModel_SD_1e-3.mat").DMR; % distributed model results
% N = length(DMR.snap.Qds)*10;

% DMR.snap.PvNDs = DMR.pVascND(3, 3);
% DMR.snap.Pvs = DMR.snap.Ps;
% subplot(131); histogram(dmr.snap.Qds, N, "EdgeAlpha",0.2);
% subplot(132); histogram(dmr.snap.Pds, N); 
%%
subplot(231); cla; hold on;
% b = DMR.snap.Qds.*DMR.Qpd;
bar(DMR.snap.Qds, DMR.Qpd/sum(DMR.Qpd)*100);
title(sprintf('Flow pseudodistribution at N = %d', N));
ylabel('Count (%)');xlabel('BLood flow');
% plot(DMR.snap.Qds,DMR.snap.Vds, 'x');
% title(['Blood flow distribution in alveoli, Q_{total} = ' num2str(round(sum(DMR.snap.Qds), 2)) ' L/min']);
% xlabel('Q');ylabel('V');
% Qsm = mean(DMR.snap.Qds);
% Vsm = mean(DMR.snap.Vds);
% plot([Qsm Qsm], [min(DMR.snap.Vds) max(DMR.snap.Vds)], 'r--')
% plot([min(DMR.snap.Qds) max(DMR.snap.Qds)], [Vsm Vsm], 'r--')
%%
subplot(232); cla;hold on;
% plot(DMR.snap.Pds,1:numel(DMR.snap.Qds), 'x'); 
b = DMR.snap.Qds.*DMR.Qpd;
bar(DMR.snap.Vds, DMR.Qpd/sum(DMR.Qpd)*100);
title(sprintf('Flow pseudodistribution at N = %d', N));
ylabel('Count (%)');xlabel('BLood flow');


histogram(DMR.snap.Pds, N/10);
title(['Partial pressure in distributed alveoli, Pa = ' num2str(round(DMR.snap.Pvs, 2))]);
plot([DMR.snap.PvNds DMR.snap.PvNds], [0 max(ylim)], 'c:', 'LineWidth',1.5)
plot([DMR.snap.Pvs DMR.snap.Pvs], [0 max(ylim)], 'r--')
legend('Capillary pO2', '1 comp pO2', 'dist venous pO2', 'Location','northwest')
xlabel('P_{O2}');ylabel('Count in category');

subplot(233); cla;hold on;
[X Y] = meshgrid(DMR.CO, DMR.Vp);
% surf(X, Y, DMR.Pv');
contour(X, Y, DMR.pVascND', '--', 'ShowText','on', 'LineWidth',2);
contour(X, Y, DMR.Pv', '-', 'ShowText','on', 'LineWidth',1);

xlabel('CO L/min');
ylabel('V (L/min)');
zlabel('P_O2 (mmHg)');
legend('1 comp', 'Distributed')
colorbar
% view(2)

% PLot surf
subplot(223);cla;hold on;
surf(X, Y, DMR.pVascND', 'LineStyle','--', 'LineWidth',1);
surf(X, Y, DMR.Pv', 'LineStyle','-', 'LineWidth',1);
xlabel('CO L/min');
ylabel('V (L/min)');
zlabel('P_O2 (mmHg)');
colorbar
% view(2)
legend('Distributed', 'Single comp')

% plot Px
% figure(3);clf;hold on;
subplot(224);cla;hold on;
serco = 1:1:numel(DMR.CO); % selected range CO
servp = 3:3:numel(DMR.Vp);% selected range Vp
plot(DMR.CO(serco), DMR.pVascND(serco, servp), 'o-');
set(gca, 'ColorOrderIndex', 1)
plot(DMR.CO(serco), DMR.Pv(serco, servp), 'x--');
xlabel('whole organ CO (L/min)');ylabel('pO2 (mmHg)')
title('Distributed alveoli with normal distribution of flow at Vp rates')
legend(string([DMR.Vp(servp) DMR.Vp(servp)]), 'NumColumns',2)
% legend('Single compartment', 'Distributed')
% legend('SD = 0.001', 'SD = 0.1','SD = 1');
saveas(gcf, sprintf("figures/DistributedModelD_VQ_SD_%1.1e.fig", SD));
saveas(gcf, sprintf("figures/DistributedModelD_VQ_SD_%1.1e.png", SD));