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
D     = 385/M;      %apparent diffusion (L/min/kg) - hand-tuned
% D     = 285/M;      %apparent diffusion (L/min/kg)

% Pair  = 164;    %atmospheric oxygen partial pressure (mmHg) - realistic
% Pair  = 150;    %atmospheric oxygen partial pressure (mmHg) - original model
Pair  = 134;    %atmospheric oxygen partial pressure (mmHg) - hand tuned

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

%% Plot West data
% distribution of concentrations
xn = 1:numel(qrs);
xl = [xn(1) xn(end)];

clf; 
subplot(223);hold on;
plot(xn, westData.cO2, 'mo-', 'MarkerSize', 12, 'LineWidth',2);
cmart = interp1(HbDisP,HbDisC, 97, "linear"); % c mixed arterial from West's pO2
plot(xl, [cmart cmart], 'm--', 'LineWidth',3);
set(gca,'ColorOrderIndex',1)
% legend('Capillary cO2', '1 comp cO2', 'dist venous cO2', 'data West per chunk', 'West total', 'Location','southwest')

% distribution of partial pressure
subplot(224);hold on;
plot(xn, westData.PO2, 'mo-', 'MarkerSize', 12, 'LineWidth',2);
xn = xlim();
% plot(xl, [97 97], 'm--', 'LineWidth',3, 'Color', 'm');
legend('Capillary pO2', '1 comp pO2', 'dist venous pO2', 'data West p02', 'West total pO2', 'Location','northeast')
set(gca,'ColorOrderIndex',1)
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
plot([0,12], [0 12], 'k--', 'LineWidth',1); % identity line
% plot(Qs./M, vrsMatch./M, 'mo-', 'MarkerSize', 12, 'LineWidth',2);
xlabel('Perfusion blood');
ylabel('Perfusion air');
legend('Identity', 'Optimal VQ');
xlim([0 12])
set(gca,'ColorOrderIndex',1)

%% Run the West with equal ventilation
V_total = sum(vu.*Ms);
vu_eq = ones(size(vu))*V_total/sum(Ms); % equal distribution scaled to total ventilation

[pv_eq, cv, p, c_i, validIds] = calculateDistributedAlveoliD(par, vu_eq, qu, Ms, D, 1, true);

%% Calculate West distribution

% output pv1 too - single compartment pO2 pv1
[pv_west, cv, p, c_i, validIds, pv1] = calculateDistributedAlveoliD(par, vu, qu, Ms, D, 1, true);


%% optimal V distribution at given total ventilation
V_target = V_total; % target ventilation to keep

vu = vrs./Ms; % start from West's distribution - is there any better?

vuRed = vu(1:end-1); % reduced set missing the last element, because the sum is constrained
newVuFun = @(vuRed)[vuRed;(V_target - sum(vuRed.*Ms(1:end-1)))/Ms(end)]; % the last element is calculated so that the sum equals V_total
% newVuFun(vuRed)
costFun = @(vuRed) -calculateDistributedAlveoliD(par, newVuFun(vuRed), qu, Ms, D, 1, false);
% costFun(vuRed)

options = optimset('Display','iter', 'TolFun', 1e-3, 'TolX', 1e-2, 'PlotFcns', @optimplotfval, 'MaxIter', 500);

% [vuRedOpt,FVAL,EXITFLAG,OUTPUT] = fminsearch(costFun, vuRed, options);
vuOpt = newVuFun([0.9696    2.3050    3.2172    4.4220    5.3585    6.2412    7.4915    9.1125]'); % result of the optim

pv_opt = calculateDistributedAlveoliD(par, vuOpt, qu, Ms, D, 1, true);

%% legends and stuff
subplot(243);
legend('Equal ventilation', 'West', 'Optimal ventilation');

subplot(244);
legend('Identity', 'Equal ventilation', 'West', 'Optimal ventilation');

subplot(223); cla;
bar([pv_eq, 97, pv_west, pv_opt, pv1])
set(gca, 'XTick', 1:5, 'XTickLabel', {'Equal ventilation'; 'West (data)';'West (sim)';'Optimal';'Single compartment'});
ylim([90, 110])

subplot(224);
legend('West (Data)', 'Equal ventilation', 'West (Simulation)', 'Optimal ventilation');


%% Match ventiolation to optimal pO2 target during exercise

CO0 = sum(qrs); % West's baseline
CO_t = [6, 10, 15, 20, 25];
CO_t = [5, 10, 15, 20];
scaleCOproportionally = false;
pO2_target = 100; % Each compartment ventilation is optimized to have this optimal pO2

qu_i = zeros(numel(CO_t), numel(qu));
vu_i = zeros(numel(CO_t), numel(qu));
c_i = zeros(numel(CO_t), numel(qu));

for i_el = 1:numel(CO_t) % iterate Exercise Level

    vu = vrs./Ms; % ventilation per mass
    qu = qrs./Ms; % perfusion per mass
    
    if ~scaleCOproportionally
        % adding perfusion equally
        qu = qu + (CO_t(i_el) - sum(qu.*Ms))/sum(Ms);
    
    else
        % scaling proportionally
        qu = qu*CO_t(i_el)/CO0;
    end
    qu_i(i_el, :) = qu;
    % sum(qu.*Ms) % check
    %% matfch ventilation to ideal target
    HbLookUp = load('Lookup.mat'); %outputs Hb dissociation curve lookup table
    HbDisP = HbLookUp.LOOK.Plookup;
    HbDisC = HbLookUp.LOOK.Clookup;
    NN = 100;
    for i_chunk = 1:length(qu) % iterate submodel flows
        %%
        % i = 1
        sig = 1; % optimal VQ, initial guess is 1:1
        step = 0.5;
        for i = 1:11 % iterate to optimize
        vu_i(i_el, i_chunk) = sig*qu(i_chunk); % ventilation iteration
        [~,~, Cvi, Pvi] = modelD_SS_relaxation(NN,par,HbDisP,HbDisC,vu_i(i_el, i_chunk),qu(i_chunk)); %loop through q and v and store fixed points
        E = Pvi(end) - pO2_target;
        c_i(i_el, i_chunk) = Cvi(end);
        if E < 0
            % increase
            sig = sig + step;
            % step = step/2;
        else % decrease
            sig = sig - step;
            step = step/2;
        end
        % vu_i = vu_i - max(min(E, vu_i/3), -vu_i/3);
        
        end
    
        fprintf('For qu %1.2f the vu is found to be %1.2f (Err %1.3f) \n', qu(i_chunk), vu_i(i_el, i_chunk), E);
      
        % p(i) = Pvi(end); % pulmonary end-capillary pO2
    end
    
    cv(i_el) = sum(c_i(i_el, :)'.*qu.*Ms)/sum(qu.*Ms);
    % cv = sum(c.*q./Ms)/sum(q./Ms);
    pv(i_el) = interp1(HbDisC, HbDisP,cv(i_el) , "linear"); % Pulmonary venous distributed sum
    fprintf('For CO of %1.2f the total optimal ventilation is %1.2f \n', sum(qu.*Ms), sum(vu_i(i_el, :)'.*Ms));
end

%% get the concentrations
p_i = reshape(interp1(HbDisC, HbDisP,c_i(:) , "linear"), size(c_i)); % Pulmonary venous distributed sum

%% Load stored vals
% load env_adding.mat
% load env_scaling.mat

%% plot that
figure(2); clf; 
subplot(221); hold on;
bar(qu_i');
title('Perfusion per element (Qu)');
if scaleCOproportionally
    subtitle('Given by West and scaled by addition');
else
    % proportionally
    subtitle('Given by West and scaled proportionally');
end
xlabel('# element')
ylabel('Pefusion L/min/kg')

subplot(222); hold on;
bar(vu_i');
title('Ventilation per element Vu')
subtitle(sprintf('Optimized for pO2 %1.0f mmHg in each compartment', pO2_target));
legend('5 L/min', '10 L/min', '15 L/min', '20 L/min', 'Location', 'NorthWest');
xlabel('# element')
ylabel('Ventilation L/min/kg air')


subplot(223); hold on;
plot((vu_i./qu_i)', 'd-');
plot([0, 9], [1 1], 'k--');
% scatter(repmat(1:size(p_i, 1), [1 9])', (vu_i(:)./qu_i(:)), 20, p_i(:), 'filled')
title('V/Q ratio')
xlabel('# element')
ylabel('-')


subplot(224);hold on;
plot(CO_t, sum(vu_i'.*Ms), '-', 'LineWidth',2)
scatter(CO_t, sum(vu_i'.*Ms), 45, pv, 'filled')
scatter(qu_i(2:end), vu_i(2:end), 20, p_i(2:end)) % get rid of the mismatched pO2 for near zero flow, as it zooms the coloring out
colormap(jet)
colorbar
title('Total V/Q')
legend('Total V/Q', 'Total pO2', 'Chunk pO2', 'Location', 'northwest')
xlabel('CO')
ylabel('V')


fontsize(gcf, 18, 'points');
%%

plot(qu, vu_i, 'd-');
plot([0 max(qu(end), vu(end))], [0 max(qu(end), vu(end))], 'k--');
