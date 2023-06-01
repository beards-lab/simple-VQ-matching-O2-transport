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
%load optimized diffusion (D) parameter
% load('ModelD_optimization_v2_results.mat','JD','DD')
% [~, jDpi] = min(JD); DDp = DD(jDpi);D = DDp;


M = .445 + .395; % Mean lung weight in young normal adult man (kg)
D     = 385/M;      %apparent diffusion (L/min/kg)
% D     = 285/M;      %apparent diffusion (L/min/kg)
Pair  = 134;    %atmospheric oxygen partial pressure (mmHg) - 164?
% Pair  = 150;    %atmospheric oxygen partial pressure (mmHg) - 164?
Pin   = 45;     %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)
Vvasc = 1;      %volume of vascular space (ml)
Valv  = 1;      %alveolar volume (ml)
l     = 1; %length of capillary
alpha = 1.3e-6*1e3;  % O2 solubility  in water/plasma(mM/mmHg)
beta  = 16800*1e-3; % O2 solubility in air (mmHg/mM)


par = [D Pair Pin alpha beta l M];
%
figure(1);clf;

% tune down the diffusion
% D = 385;
Ms = M*Vols/100;
vu = vrs./Ms; % ventilation per mass
qu = qrs./Ms; % perfusion per mass
%%
[pv, cv, p, c_i, validIds] = calculateDistributedAlveoliD(par, vu, qu, Ms, D, 1);
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
%%
% distribution of concentrations
xn = 1:numel(qrs);
xl = [xn(1) xn(end)];

subplot(223);hold on;
plot(xn, westData.cO2, 'o-');

cmart = interp1(HbDisP,HbDisC, 97, "linear"); % c mixed arterial from West's pO2
plot(xl, [cmart cmart], 'm--');
legend('Capillary cO2', '1 comp cO2', 'dist venous cO2', 'data West per chunk', 'West total', 'Location','southwest')

% distribution of partial pressure
subplot(224);hold on;
plot(xn, westData.PO2, 'o-');
xn = xlim();
plot(xl, [97 97], 'm--');
legend('Capillary pO2', '1 comp pO2', 'dist venous pO2', 'data West p02', 'West total pO2', 'Location','northeast')
% PLot over our optimal VQ matching - this moght not really work!
% we compare exercise matching with gravitational matching

% Best fit Exponential Curve Fit (exp2)
a    = 2.3063; 
b    = 0.1459;  
c_i    = 0.0032;   
d    = 0.5296;
VfunCo = @(CO) a*exp(b*CO) + c_i*exp(d*CO);

% scale to normal size
Qs = qrs*M./Ms;
Vs = vrs.*M./Ms;

% ideal VQ matching
vrsMatch = VfunCo(Qs);

subplot(244);hold on;
plot(Qs./M, vrsMatch./M, 'o-');
xlabel('Perfusion blood');
ylabel('Perfusion air');
legend('West', 'Optimal VQ')

%% Match ventilation

pO2_target = 100;
NN = 100;
HbLookUp = load('Lookup.mat'); %outputs Hb dissociation curve lookup table
HbDisP = HbLookUp.LOOK.Plookup;
HbDisC = HbLookUp.LOOK.Clookup;
vu_i = zeros(size(qu));
c_i = zeros(size(qu));

CO0 = sum(qu.*Ms);
% CO_t

for i = 1:length(qu) % iterate submodel flows
    % i = 1
    sig = 1; % optimal VQ, initial guess is 1:1
    step = 0.5;
    for j = 1:11
    vu_i(i) = sig*qu(i); % ventilation iteration
    [~,~, Cvi, Pvi] = modelD_SS_relaxation(NN,par,HbDisP,HbDisC,vu_i(i),qu(i)); %loop through q and v and store fixed points
    E = Pvi(end) - pO2_target;
    c_i(i) = Cvi(end);
    if E < 0
        % increase
        sig = sig + step;
        step = step/2;
    else % decrease
        sig = sig - step;
        step = step/2;
    end
    % vu_i = vu_i - max(min(E, vu_i/3), -vu_i/3);
    
    end

    fprintf('For qu %1.2f the vu is found to be %1.2f (Err %1.3f) \n', qu(i), vu(i), E);
  
    % p(i) = Pvi(end); % pulmonary end-capillary pO2
end

cv = sum(c_i.*qu.*Ms)/sum(qu.*Ms);
% cv = sum(c.*q./Ms)/sum(q./Ms);
pv = interp1(HbDisC, HbDisP,cv , "linear"); % Pulmonary venous distributed sum


%% plot that
figure(2); clf; hold on;
plot(qu, vu_i, 'd-');
plot([0 max(qu(end), vu(end))], [0 max(qu(end), vu(end))], 'k--');
