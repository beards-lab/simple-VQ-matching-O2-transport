%% get a random distribution
N = 100;
CO = 5; % Cardiac output mL/min

qrs = CO/N + CO/10*randn(N, 1)/N;% random flow set

histogram(qrs);

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
%%

HbLookUp = load('Lookup.mat'); %outputs Hb dissociation curve lookup table
HbDisP = HbLookUp.LOOK.Plookup;
HbDisC = HbLookUp.LOOK.Clookup;
Ceq = zeros(size(HbDisP));
for iqs = 1:length(HbDisP)
    Ceq(iqs) = alpha*HbDisP(iqs)+C0*((HbDisP(iqs)^n)/(HbDisP(iqs)^n+P50^n)); % transform to Pa
end

% test interp1
Cx = 1:1.7:12;
Px = interp1(HbDisC, HbDisP,Cx , "linear");
plot(HbDisP, HbDisC, HbDisP, Ceq, 'o', Px, Cx, 'x-', 'Linewidth', 2);



%%
opt = optimset('MaxFUnEvals',100, 'Display', 'None');
tic


N = 100;
SD = 1e-1;
rng(4); % set the random generator seed 
rnd = SD*randn(N, 1); % random normal distribution

CO = 3:1:20;
Vp = 3:1:20;
% CO = 5;
% Vp = 5;
Px = zeros(length(CO), length(Vp));
%
for iq = 1:length(CO) % iterate flow
    disp("Jedu " + num2str(round(iq/length(CO)*100)))
    qrs = CO(iq)/N * (1 + rnd); % random flow set

    for iv = 1:length(Vp) % iterate ventilation
        vrs = Vp(iv)/N *(1 + rnd); % ventilation random set
        vrs = ones(size(qrs))*Vp(iv)/N; % start simple
        
        % distribution within lungs
        c = zeros(size(qrs));        
        for iqs = 1:length(qrs) % iterate submodel flows
            par(12) = qrs(iqs); % flow in 
            par(11) = vrs(iqs);
            [P, ~, FLAG(iqs)] = fsolve(@ModelB_FixedPoint_Objective2,[100 100],opt,par);
            [~, oc] = ModelB_FixedPoint_Objective2(P, par); % get the contrectration
            c(iqs) = oc(3); % arterial concentration
        end
        Q = sum(qrs); % sum of flows
        Ca = sum(c.*qrs)/Q; % weighted average of concentrations by flow
        Px(iq, iv) = interp1(HbDisC, HbDisP,Ca , "linear");
    % Px
    end
end
%%
% clf;
hold on;
CON = 30:10:150;
[xq,yv] = meshgrid(CO,Vp);
contour(xq, yv, Px, CON,'ShowText','on');
xlabel('Q');ylabel('V');