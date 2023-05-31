function [pv, cv, p, c, validIds] = calculateDistributedAlveoliD(par, v, q, Ms, D, ExcludeNaNs)
% Calculates model values for ventilation VRS and perfusion QRS vectors.
% Usage:
%   [pv, cv, p, c, validIds] = calculateDistributedAlveoliD(par, vrs, qrs, Pd, Ds, ExcludeNaNs)
%
% where par is a parameter vector par = [D Pair Pin alpha beta l W]
%   vrs - ventilation rate set per compartment (L/min/kg)
%   qrs - perfusion rate set per compartment (L/min/kg)
%   Ms - mass set of compartments (Kg)
%   Ds - apparent diffusion coefficients
%   ExludeNaNs - boolean whether to exclude resulting NaN concentrations into output set 
% outputs


D = par(1);
N = length(q);
% if isempty(D)
%     % D = ones(N, 1)*D/N;
% end

NN = 100;
HbLookUp = load('Lookup.mat'); %outputs Hb dissociation curve lookup table
HbDisP = HbLookUp.LOOK.Plookup;
HbDisC = HbLookUp.LOOK.Clookup;

% equal Q and Vp in all alveoli
par(1) = D;
[pv1,~, cv1, ~] = modelD_SS_relaxation(NN,par,HbDisP,HbDisC,sum(v.*Ms),sum(q.*Ms)); 
cv1 = cv1(end);
        
tic
% distribution within lungs
c = zeros(size(v));
p = zeros(size(v));
for i = 1:length(v) % iterate submodel flows
    par(1) = D;
    [~,~, Cvi, Pvi] = modelD_SS_relaxation(NN,par,HbDisP,HbDisC,v(i),q(i)); %loop through q and v and store fixed points
    c(i) = Cvi(end); % pulmonary end-capillary concentration
    p(i) = Pvi(end); % pulmonary end-capillary pO2
end
t = toc;
if ExcludeNaNs
    validIds = ~isnan(c);
    c = c(validIds);
    p = p(validIds);
    % % store the nans
    qrsNans = q(~validIds);
    vrsNans = v(~validIds);
    PdNans = Ms(~validIds);
    % redeclare number of elements
    N = numel(c);
    sn = sum(~validIds);% sum of nans
    % Get rid of nans
    q = q(validIds);
    v = v(validIds);
    Ms = Ms(validIds);
    Ms = Ms/sum(Ms); % Sum of Pd = 1
else
    validIds = ones(size(c));
    sn = 0; qrsNans = [];vrsNans = [];PdNans = [];
    Ms = Ms/sum(Ms); % Sum of Pd = 1 - just in case here
end
    
%% c(isnan(c)) = 0; p(isnan(p)) = 0;
% Peld = interp1(HbDisC, HbDisP,c , "linear"); % partial pressure per ventilation element
Cvd = (c.*q.*Ms); % concentration venous distributed, weighted by "probability distribution"
% Pvd = interp1(HbDisC, HbDisP,Cvd , "linear"); 
cv = sum(Cvd)/sum(q); % weighted average of concentrations by flow
% just another mean of calculation of the previous
cv = sum(c.*q.*(Ms*1/N))/sum(q); 
cv = sum(c.*q.*Ms)/sum(q.*Ms);
% cv = sum(c.*q./Ms)/sum(q./Ms);
pv = interp1(HbDisC, HbDisP,cv , "linear"); % Pulmonary venous distributed sum

fprintf('At Q = %2.1f and Vp  = %2.1f, pO2_{dist} = %2.1f (single comp. = %2.1f), with %1.0f NaNs in %2.0f ms \n', ...
    sum(q), sum(v), pv, pv1, sn, t*1000);

%% plot that

% distribution of concentrations
% figure(3);clf;
% bar(vrs, c);title('Distribution of concentration (raw)');
clf;
% PdPrc = Ms/sum(Ms)*100;
subplot(241);
xn = 1:numel(Ms);
xl = [xn(1), xn(end)];
bar(xn, Ms*1000);title(sprintf('Size of compartments (M_T = %1.2f kg)', sum(Ms)));
% plot(qrs, PdPrc, '^');title('Pseudodistribution of perfused elements');
xlabel('Element # ');ylabel('Compartment size (g)')
% xl = xlim;

subplot(242);hold on;
bar(xn, q);
xlabel('Element # ');ylabel('Compartment flow (mL/min/kg)')

% plot(vrsNans, qrsNans, 'rx');
% title(sprintf('Flow (%1.1f L/min) to ventilation (%1.1f L/min) relation', sum(qrs), sum(vrs)));
% xlabel('Perfusion L/min');ylabel('Ventilation L/min');

subplot(243);hold on;
bar(xn, v);
xlabel('Element # ');ylabel('Compartment ventilation (mL/min/kg)')

subplot(244);hold on;
plot(q, v, 'd-')
plot(vrsNans*1000, qrsNans*1000, 'rx');
title(sprintf('Flow (Σ=%1.1f L/min) to ventilation (Σ=%1.1f L/min) relation', sum(q), sum(v)));
xlabel('Perfusion L/min/kg');ylabel('Ventilation mL/min/kg');

% % bar(vrs, Pd);title('Pseudodistribution of ventilated elements');
% plot(vrs, PdPrc, 'v');title('Pseudodistribution of ventilated elements');
% xlabel('Ventilation L/min');ylabel('Chunk size (%)');


% % distribution of concentrations
subplot(223);hold on;
plot(xn, c, 'x-');title(sprintf('Weighted concentration, total: %0.2f (single %0.2f, %0.2f%%)', cv, cv1, 100 - cv1/cv*100));
plot(xn(~validIds), max(c)*ones(size(vrsNans)), '*-');
plot(xl, [cv1 cv1], 'r--')
plot(xl, [cv cv], 'c:', 'LineWidth',1.5)
legend('Capillary cO2', '1 comp cO2', 'dist venous cO2', 'Location','northwest')
% xlim(xl);

% distribution of partial pressure
subplot(224);hold on;
plot(xn, p, 'x-');title(sprintf('Weighted partial pressures, total: %0.2f  (single %0.2f, %0.2f%%)', pv, pv1, 100 - pv1/pv*100));
plot(xn(~validIds), max(p)*ones(size(vrsNans)), '*-');
plot(xl, [pv1 pv1], 'r--')
plot(xl, [pv pv], 'c:', 'LineWidth',1.5)
legend('Capillary pO2', '1 comp pO2', 'dist venous pO2', 'Location','northwest')
%%
% Pdperc = Pd/sum(Pd)*100; % probability distribution in percent
% figure(2);clf
% subplot(221);
% plot(c, Pdperc, 's-');
% xlabel('Concentration');ylabel('Frequency %');
% subplot(222);
% plot(p, Pdperc, 'o-'); xlabel('pO2');ylabel('Frequency %');
% subplot(234);hold on;
% bar(c.*qrs);xlabel('# element');ylabel('c*Q');
% subplot(235);hold on;
% bar(c.*Pd);xlabel('# element');ylabel('c*Pd');
% subplot(236);
% bar(c.*qrs.*Pd); xlabel('# element');ylabel('c*Q*Pd');

disp('Been there, don det');