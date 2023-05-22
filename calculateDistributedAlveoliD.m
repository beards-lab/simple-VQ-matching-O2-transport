function [pv, cv, p, c, validIds] = calculateDistributedAlveoliD(par, vrs, qrs, Pd, ExcludeNaNs)
D = par(1);
N = length(qrs);
NN = 100;
HbLookUp = load('Lookup.mat'); %outputs Hb dissociation curve lookup table
HbDisP = HbLookUp.LOOK.Plookup;
HbDisC = HbLookUp.LOOK.Clookup;

% equal Q and Vp in all alveoli
par(1) = D;
[pv1,~, cv1, ~] = modelD_SS_relaxation(NN,par,HbDisP,HbDisC,sum(vrs),sum(qrs)); 
cv1 = cv1(end);
        
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
if ExcludeNaNs
    validIds = ~isnan(c);
    c = c(validIds);
    p = p(validIds);
    % % store the nans
    qrsNans = qrs(~validIds);
    vrsNans = vrs(~validIds);
    PdNans = Pd(~validIds);
    % redeclare number of elements
    N = numel(c);
    sn = sum(~validIds);% sum of nans
    % Get rid of nans
    qrs = qrs(validIds);
    vrs = vrs(validIds);
    Pd = Pd(validIds);
    Pd = Pd/sum(Pd); % Sum of Pd = 1
else
    validIds = ones(size(c));
    sn = 0; qrsNans = [];vrsNans = [];PdNans = [];
    Pd = Pd/sum(Pd); % Sum of Pd = 1 - just in case here
end
    
%% c(isnan(c)) = 0; p(isnan(p)) = 0;
% Peld = interp1(HbDisC, HbDisP,c , "linear"); % partial pressure per ventilation element
Cvd = (c.*qrs.*Pd); % concentration venous distributed, weighted by "probability distribution"
% Pvd = interp1(HbDisC, HbDisP,Cvd , "linear"); 
cv = sum(Cvd)/sum(qrs); % weighted average of concentrations by flow
% just another mean of calculation of the previous
cv = sum(c.*qrs.*(Pd*1/N))/sum(qrs); 
cv = sum(c.*qrs.*Pd)/sum(qrs.*Pd);
pv = interp1(HbDisC, HbDisP,cv , "linear"); % Pulmonary venous distributed sum

fprintf('At Q = %2.1f and Vp  = %2.1f, pO2_{dist} = %2.1f (single comp. = %2.1f), with %1.0f NaNs in %2.0f ms \n', ...
    sum(qrs), sum(vrs), pv, pv1, sn, t*1000);

% plot that

% distribution of concentrations
% figure(3);clf;
% bar(vrs, c);title('Distribution of concentration (raw)');
clf;
subplot(221);
bar(vrs, Pd);title('Pseudodistribution of ventilated elements');
xlabel('Ventilation L/min');ylabel('Weight (-)')

subplot(222);hold on;
plot(vrs, qrs, 'o');
plot(vrsNans, qrsNans, 'rx');
title(sprintf('Flow (%1.1f L/min) to ventilation (%1.1f L/min) relation', sum(qrs), sum(vrs)));
ylabel('Perfusion L/min');xlabel('Ventilation L/min')
xl = xlim;

% subplot(221);
% plot(qrs, vrs, 'x');title(fprintf('Flow (%1.1f L/min) to ventilation (%1.1f L/min) distribution', sum(qrs), sum(vrs)));
% xlim('Perfusion L/min');ylim('Ventilation L/min')

% distribution of concentrations
subplot(223);hold on;
plot(vrs, c, '*-');title(sprintf('Weighted concentration, total: %0.2f (single %0.2f, %0.2f%%)', cv, cv1, 100 - cv1/cv*100));
plot(vrsNans, max(c)*ones(size(vrsNans)), '*-');
plot(xl, [cv1 cv1], 'r--')
plot(xl, [cv cv], 'c:', 'LineWidth',1.5)
legend('Capillary cO2', '1 comp cO2', 'dist venous cO2', 'Location','northwest')
xlim(xl);

% distribution of partial pressure
subplot(224);hold on;
plot(vrs, p, '*-');title(sprintf('Weighted partial pressures, total: %0.2f  (single %0.2f, %0.2f%%)', pv, pv1, 100 - pv1/pv*100));
plot(vrsNans, max(p)*ones(size(vrsNans)), '*-');
plot(xlim, [pv1 pv1], 'r--')
plot(xlim, [pv pv], 'c:', 'LineWidth',1.5)
legend('Capillary pO2', '1 comp pO2', 'dist venous pO2', 'Location','northwest')
%%
Pdperc = Pd/sum(Pd)*100; % probability distribution in percent
figure(2);clf
subplot(221);
plot(c, Pdperc, 's-');
xlabel('Concentration');ylabel('Frequency %');
subplot(222);
plot(p, Pdperc, 'o-'); xlabel('pO2');ylabel('Frequency %');
subplot(234);hold on;
bar(c.*qrs);xlabel('# element');ylabel('c*Q');
subplot(235);hold on;
bar(c.*Pd);xlabel('# element');ylabel('c*Pd');
subplot(236);
bar(c.*qrs.*Pd); xlabel('# element');ylabel('c*Q*Pd');

disp('Been there, don det');