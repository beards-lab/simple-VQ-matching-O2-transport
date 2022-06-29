clear; close all; clc;
% This script solves Model A and saves results into a .mat file + generates
% fun animated .gifs to show how the diffusion parameter affects model
% behavior

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

%%% look up table
load('Lookup.mat') %outputs LOOK
P = LOOK.Plookup;
C = LOOK.Clookup;

%%% solve system of ODEs
% X0 = [interp1(P,C,Pin) Pin];
% [t,X] = ode15s(@ModelB_RHS,[0 10],X0,[],par,LOOK);
% 
% Cvasc = X(:,1); 
% Pvasc = interp1(C,P,Cvasc);
% Palv  = X(:,2);


opt = optimset('MaxFUnEvals',20000);
PSIM = fsolve(@ModelB_FixedPoint_Objective2,[100 100],opt,par);
pvasc_sim = PSIM(1);
palv_sim  = PSIM(2);

%%% grid for contour plot
de = 0.5;
q  = 0:de:70;
v  = 0:de:70;
np = length(q);
[Q,V] = meshgrid(q,v);

Dv = 50; %0:1:100; %vector of diffusion parameters to explore
Dpvasc = zeros(length(q),length(v),length(Dv));
Dpalv  = zeros(length(q),length(v),length(Dv));

%%% computing contour plot
% preallocating cell arrays of parameters to run for-loops in parallel
tpar{np,np} = [];
for i = 1:np
   for j = 1:np
       tpar{i,j} = par;  tpar{i,j}(11) = v(i); tpar{i,j}(12) = q(j);
   end
end

pvasc = zeros(np);
palv  = zeros(np);
FLAG = zeros(np);
tic;
parfor i = 1:np
   for j = 1:np
%        clear tpar P
%        tpar = par; tpar(11) = v(i); tpar(12) = q(j);
       [P, ~, FLAG(i,j)] = fsolve(@ModelB_FixedPoint_Objective2,[100 100],opt,tpar{i,j});
       pvasc(i,j) = P(1); palv(i,j)  = P(2);
       
       disp([i j])
   end
end
dpac  = palv - pvasc;
toc;


% preallocating cell arrays of parameters to run for-loops in parallel
Dtpar{np,np,length(Dv)} = [];
for k = 1:length(Dv)
    for i = 1:np
        for j = 1:np
            Dtpar{i,j,k} = par;  Dtpar{i,j,k}(11) = v(i); Dtpar{i,j,k}(12) = q(j); Dtpar{i,j,k}(1) = Dv(k);
        end
    end
end


tic;
DFLAG = zeros(np,np,np);
parfor k = 1:length(Dv)
    for i = 1:np
        for j = 1:np
%             clear tpar P
%             tpar = par; tpar(11) = v(i); tpar(12) = q(j); tpar(1) = Dv(k);
            [P, ~, DFLAG(i,j,k)] = fsolve(@ModelB_FixedPoint_Objective2,[100 100],[],Dtpar{i,j,k});
            Dpvasc(i,j,k) = P(1); Dpalv(i,j,k)  = P(2);
            
            disp([i j k])
        end
    end
end
Ddpac  = Dpalv - Dpvasc;
toc;

save('ModelB_results.mat')

%%% plots
% figure;
% plot(t,Pvasc,t,Palv,'linewidth',2)
% xlabel('Time (s)')
% ylabel('Oxygen Tension (mmHg)')
% legend('Vasc', 'Alv')

figure;
surf(Q,V,pvasc, 'EdgeColor','none')
title('Pvasc')
xlabel('q')
ylabel('v')
axis equal

figure;
surf(Q,V,palv, 'EdgeColor','none')
title('Palv')
xlabel('q')
ylabel('v')
axis equal

figure;
surf(Q,V,dpac, 'EdgeColor','none')
title('\Delta P')
xlabel('q')
ylabel('v')
axis equal


CON = -30:5:150;
x = 0:70;
figure;
subplot(1,3,1)
contour(Q,V,pvasc,CON,'ShowText','on')
hold on
plot(x,x,'k--')
hold on
plot(5,5,'rx')
set(gca,'fontsize',18)
title('Vascular (mmHg)')
xlabel('Cardiac Output (ml/s)')
ylabel('Ventilation Magnitude (ml/s)')
axis equal
grid on

% figure;
subplot(1,3,2)
contour(Q,V,palv,CON,'ShowText','on')
hold on
plot(x,x,'k--')
hold on
plot(5,5,'rx')
set(gca,'fontsize',18)
title('Alveolar Space (mmHg)')
xlabel('Cardiac Output (ml/s)')
axis equal
grid on

% figure;
subplot(1,3,3)
contour(Q,V,dpac,CON,'ShowText','on')
hold on
plot(x,x,'k--')
hold on
plot(5,5,'rx')
set(gca,'fontsize',18)
title('Gradient (mmHg)')
xlabel('Cardiac Output (ml/s)')
axis equal
grid on



% %%% making gifs
% CON = -30:5:150;
% x = 0:10;
% h = figure('units','normalized','outerposition',[0 0 1 1]);
% for i = 1:length(Dv) %iterate by size of Dv
%     %%%% make contour plots
%     cla(gca)
%     
%     contour(Q,V,squeeze(Dpvasc(:,:,i)),CON,'ShowText','on')
%     hold on
%     plot(x,x,'k--')
%     set(gca,'fontsize',18)
%     xlabel('Cardiac Output (ml/s)')
%     ylabel('Ventilation Magnitude (ml/s)')
%     text(9,1,['D = ', num2str(Dv(i))],'fontsize',18)
%     axis equal
%     grid on
%     
%     if i == 1
%         gif('ModelB_vasc.gif','DelayTime',1/15) %make gif file
%     else
%         gif %append frame to gif
%     end
%     disp(i)
% end
% 
% h2 = figure('units','normalized','outerposition',[0 0 1 1]);
% for i = 1:length(Dv) %iterate by size of Dv
%     %%%% make contour plots
%     cla(gca)
%     
%     contour(Q,V,squeeze(Dpalv(:,:,i)),CON,'ShowText','on')
%     hold on
%     plot(x,x,'k--')
%     set(gca,'fontsize',18)
%     xlabel('Cardiac Output (ml/s)')
%     ylabel('Ventilation Magnitude (ml/s)')
%     text(9,1,['D = ', num2str(Dv(i))],'fontsize',18)
%     axis equal
%     grid on
%     
%     if i == 1
%         gif('ModelB_alv.gif','DelayTime',1/15) %make gif file
%     else
%         gif %append frame to gif
%     end
%     disp(i)
% end
% 
% h3 = figure('units','normalized','outerposition',[0 0 1 1]);
% for i = 1:length(Dv) %iterate by size of Dv
%     %%%% make contour plots
%     cla(gca)
%     
%     contour(Q,V,squeeze(Ddpac(:,:,i)),CON,'ShowText','on')
%     hold on
%     plot(x,x,'k--')
%     set(gca,'fontsize',18)
%     xlabel('Cardiac Output (ml/s)')
%     ylabel('Ventilation Magnitude (ml/s)')
%     text(9,1,['D = ', num2str(Dv(i))],'fontsize',18)
%     axis equal
%     grid on
%     
%     if i == 1
%         gif('ModelB_grad.gif','DelayTime',1/15) %make gif file
%     else
%         gif %append frame to gif
%     end
%     disp(i)
% end

