%% get a random distribution
N = 100;
CO = 5; % Cardiac output mL/min

qrs = CO/N + CO/10*randn(N, 1)/N;% random flow set

histogram(qrs);

%% Test Hb dissociation look up table

HbLookUp = load('Lookup.mat'); %outputs Hb dissociation curve lookup table
HbDisP = HbLookUp.LOOK.Plookup;
HbDisC = HbLookUp.LOOK.Clookup;
Ceq = zeros(size(HbDisP));
for iqs = 1:length(HbDisP)
    Ceq(iqs) = alpha*HbDisP(iqs)+C0*((HbDisP(iqs)^n)/(HbDisP(iqs)^n+P50^n)); % transform to Pa
end

% test interp1
Cx = 1:1.7:12;
Padist = interp1(HbDisC, HbDisP,Cx , "linear");
plot(HbDisP, HbDisC, HbDisP, Ceq, 'o', Padist, Cx, 'x-', 'Linewidth', 2);

%% prepare random distribution
N = 500;
SD = 0.001;
rng(3); % set the random generator seed 
rnd = SD*randn(N, 1) + 1;
rnd = rnd(find(rnd > 0)); % select only non-zero flows
rnd = (rnd) / mean(rnd); % rescale to match total flow
% subplot(122);plot(rnd, '.')
% subplot(121);histogram(rnd, 'FaceAlpha',0.2);
% run CO x V scan for distributed alveoli
drawPlots = true;
if drawPlots
    % figure(1);
    % subplot(121);hold on;
    % subplot(132);hold on;
    % subplot(133);hold on;
end

%% Prepare inputs
CO = 3:0.5:20;
Vp = 3:0.5:20;
% CO = 5;
% Vp = 5;
Qs = 5; % flow snapshot
Vps = 5; % Ventilation snapshot
Padist = zeros(length(CO), length(Vp));
Pc = zeros(length(CO), length(Vp));

%% get all the param shift tigerther

%%% parameters
D     = 2200;      %apparent diffusion (ml/s)
Pair  = 150;    %atmospheric oxygen partial pressure (mmHg)
Pin   = 45;     %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)
Vvasc = 1;      %volume of vascular space (ml)
Valv  = 1;      %alveolar volume (ml)

alpha = 1.3e-6*1e3;  % O2 solubility  in water/plasma(mM/mmHg)
CHb   = 0.021*1e3;   % Hb binding site conc (mmol/L of RBC's)
Hct   = 0.40;    % hematocrit (unitless)
C0    = CHb*Hct; % blood oxygen binding capacity (mol/L)
n     = 2.7;     % Hill exponent
P50   = 27;      % half-max saturation of Hb
beta  = 16800*1e-3; % O2 solubility in air (mmHg/mM)

Vp = 5; %ventilation flow (ml/s)
Qp = 5; %blood flow (ml/s)

%load optimized diffusion (D) parameter
load('ModelB_optimization_v2_results.mat','JB','DB')
[~, jBpi] = min(JB); DBp = DB(jBpi);
D = DBp;

par = [D Pair Pin Vvasc Valv alpha beta C0 n P50 Vp Qp];

%% Run ModelB eval

opt = optimset('MaxFUnEvals',100, 'Display', 'None');
for iq = 1:length(CO) % iterate flow
    disp("Jedu " + num2str(round(iq/length(CO)*100)))
    qrs = CO(iq)/N * rnd; % random flow set

    for iv = 1:length(Vp) % iterate ventilation
        vrs = Vp(iv)/N *rnd; % ventilation random set
        vrs = ones(size(qrs))*Vp(iv)/N; % start simple - equal ventilation
        
        % equal Q and Vp in all alveoli
        [P, ~, FLAG(iqs)] = fsolve(@ModelB_FixedPoint_Objective2,[100 100],opt,par);
        Pc(iq, iv) = P(1);

        % distribution within lungs
        c = zeros(size(qrs));        
        for iqs = 1:length(qrs) % iterate submodel flows
            par(12) = qrs(iqs); % flow in 
            par(11) = vrs(iqs);
            [P, ~, FLAG(iqs)] = fsolve(@ModelB_FixedPoint_Objective2,[100 100],opt,par);
            [~, oc] = ModelB_FixedPoint_Objective2(P, par); % get the contrectration
            c(iqs) = oc(2); % arterial concentration
            p(iqs) = P(1);
        end
        Q = sum(qrs); % sum of flows
        Ca = sum(c.*qrs)/Q; % weighted average of concentrations by flow
        Padist(iq, iv) = interp1(HbDisC, HbDisP,Ca , "linear"); % Pa distributed
        if drawPlots && CO(iq) == Qs && Vp(iv) == Vps, 
            subplot(131); histogram(rnd, N/10, "EdgeAlpha",0.2);
            subplot(132); histogram(p); 
            Pds = p; % PO2 distribution snapshot
            Qds = qrs; % flow distribution snapshot
            Ps = Padist(iq, iv); %
        end
    end
end
%% save results
DistributedModelB.CO = CO; % cardiac output
DistributedModelB.Vp = Vp; % ventilation
DistributedModelB.Pv = Padist; % PO2 in pulmonary veins
DistributedModelB.Pvds = Pds; % PO2 in pulmonary veins - snapshot of distribution at Qs and Vps
DistributedModelB.Qds = Qds; % flow distribtiuon snapshot
DistributedModelB.Vps = Vps; % ventilation snapshot

save(DistributedModelB, 'DistributedModelB', '-mat');
%% load results
clear;
dmb = load(DistributedModelB);
CO = dbm.CO;
Vp = dbm.Vp;
Padist = dbm.Pv;


subplot(133); cla;hold on;
[X Y] = meshgrid(CO, Vp);
contour(X, Y, Padist', '--', 'ShowText','on');
contour(X, Y, Pc', 'ShowText','on');
% surf(X, Y, Px');
xlabel('CO L/min');
ylabel('V (L/min)');
zlabel('P_O2 (mmHg)');
colorbar
view(2)

% size(Px)
if drawPlots
    subplot(131); title(['Blood flow distribution in alveoli, Q_{total} = ' num2str(round(Qs, 2)) ' L/min']);
    xlabel('Q');ylabel('Count in category');
    subplot(132); title(['Partial pressure in distributed alveoli, Pa = ' num2str(round(Ps, 2))]);
    Ps = 121.96; hold on;
    plot([Ps Ps], [0 max(ylim)], 'r--')
    xlabel('P_{O2}');ylabel('Count in category');
end
disp(['SD ' num2str(SD) ', P_{O2, a} = ' num2str(round(Padist(end), 2)) ' mmHg']);
%% plot Px
figure(3);hold on; 
plot(CO, Padist);
xlabel('whole organ CO (L/min)');ylabel('pO2 (mmHg)')
title('Distributed alveoli with normal distribution of flow')
% legend('SD = 0.001', 'SD = 0.1','SD = 1');
%% plot contour
% clf;
hold on;
CON = 30:10:150;
[xq,yv] = meshgrid(CO,Vp);
contour(xq, yv, Padist, CON,'ShowText','on');
xlabel('Q');ylabel('V');