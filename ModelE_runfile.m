% clear; close all; clc;
% This script solves Model E and saves results into a .mat file + generates
% fun animated .gifs to show how the diffusion parameter affects model
% behavior

%%% parameters
D     = 20;      %apparent diffusion (ml/s)
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

Vp = 70; %ventilation flow (L/min)
Qp = 70; %blood flow (L/min)
VO2M = 8; %O2 consumption rate (L/min)

k = 2.5;%conversion factor to convert O2 content (ml O2/100ml blood) to concentration (mM)

% empirical function of Cardiac Ouput as a function of work rate - from Stringer 1997
% Reference --> Stringer WW, Hansen JE, Wasserman K. Cardiac output estimated noninvasively from oxygen uptake during exercise. Journal of Applied Physiology. 1997 Mar 1;82(3):908-12.
a = 0.105;
b = 5.72;


%%% grid for contour plot
de = 0.25;
q  = 0:de:Vp;
v  = 0:de:Qp;
np = length(q);
[Q,V] = meshgrid(q,v);

% %%% Setting the diffusion parameter
% load('ModelE_optimization_results.mat','D_opt_0', 'D_opt_1')
% D_0 = D_opt_0;
% 
% % assuming some kind of vascular recruitment where the apparent diffusion
% % increases as blood flow increases - D_opt is a unitless parameter here
% D_1 = D_opt_1*Q; 
% 
% % estimation D algebraically to set Pv to a physiological value at baseline
% pv_T = 100; % target vascular oxygen tension
% v_T = 6; %target basal ventilation
% q_T = 5; %target basal blood flow
% VO2_T = -VO2M.*b.*q_T./(100.*(a.*q_T-VO2M)); %target oxygen consumption
% D_2 = 1/(alpha*(k*(Pair-pv_T)/VO2_T - beta/v_T));

% setting the parameter

%load optimized diffusion (D) parameter
load('ModelE_optimization_v2_results.mat','JE','DE')
[~, jEpi] = min(JE); DEp = DE(jEpi);
D = DEp;

%%% computing contour plot
VO2 = -VO2M.*b.*Q./(a.*Q-VO2M);
pvasc = Pair - VO2.*(beta./V + 1/alpha./D)./k;
palv  = Pair-beta.*VO2./V./k ;
dpac  = palv - pvasc;


% %%% look up table
% load('Lookup.mat') %outputs LOOK
% P = LOOK.Plookup;
% C = LOOK.Clookup;
% 
% cvasc = real(interp1(P,C,pvasc));
% cin = cvasc - VO2./Q./k;
% pin = real(interp1(P,C,cin));

% %%% explore affect of apparent diffusion parameters
% Dv = 0:0.5:30; %vector of diffusion parameters to explore
% Dpvasc = zeros(length(q),length(v),length(Dv));
% Dpalv  = zeros(length(q),length(v),length(Dv));
% 
% tic;
% for i = 1:length(Dv)
%     Dpvasc(:,:,i) = Pair - VO2.*(beta./V+1/alpha/Dv(i))./k;
%     Dpalv(:,:,i)  = Pair-beta.*VO2./V./k;
%     disp(i)
% end
% Ddpac  = Dpalv - Dpvasc;
% 
% Dcvasc = interp1(P,C,Dpvasc);
% Dcin = Dcvasc - VO2./Q./k;
% Dpin = real(interp1(P,C,Dcin));
% toc;

save('ModelE_results.mat','Q','V','palv','pvasc','dpac')

%%% plots
xl = 40;
yl = 70;


% figure;
% surf(Q,V,pvasc, 'EdgeColor','none')
% title('Pvasc')
% xlabel('q')
% ylabel('v')
% % axis equal
% 
% figure;
% surf(Q,V,palv, 'EdgeColor','none')
% title('Palv')
% xlabel('q')
% ylabel('v')
% % axis equal
% 
% figure;
% surf(Q,V,dpac, 'EdgeColor','none')
% title('\Delta P')
% xlabel('q')
% ylabel('v')
% 
% 
% figure;
% surf(Q,V,pin, 'EdgeColor','none')
% title('Pin')
% xlabel('q')
% ylabel('v')
% % axis equal

CON = -30:5:150;
x = 0:max(v);
figure;
subplot(1,3,1)
contour(Q,V,pvasc,CON,'ShowText','on')
hold on
plot(x,x,'k--')
hold on
set(gca,'fontsize',18)
title('Vascular (mmHg)')
xlabel('Cardiac Output (ml/s)')
ylabel('Ventilation Magnitude (ml/s)')
% axis equal
grid on
xlim([0 xl])
ylim([0 yl])

% figure;
subplot(1,3,2)
contour(Q,V,palv,CON,'ShowText','on')
hold on
plot(x,x,'k--')
hold on
set(gca,'fontsize',18)
title('Alveolar Space (mmHg)')
xlabel('Cardiac Output (ml/s)')
% axis equal
grid on
xlim([0 xl])
ylim([0 yl])

% figure;
subplot(1,3,3)
contour(Q,V,dpac,CON,'ShowText','on')
hold on
plot(x,x,'k--')
hold on
set(gca,'fontsize',18)
title('Gradient (mmHg)')
xlabel('Cardiac Output (ml/s)')
% axis equal
grid on
xlim([0 xl])
ylim([0 yl])


% %%% making gifs
% CON = -30:5:150;
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
%     text(50,10,['D = ', num2str(Dv(i))],'fontsize',18)
%     axis equal
%     grid on
%     
%     if i == 1
%         gif('ModelB2_vasc.gif','DelayTime',1/15) %make gif file
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
%     text(50,10,['D = ', num2str(Dv(i))],'fontsize',18)
%     axis equal
%     grid on
%     
%     if i == 1
%         gif('ModelB2_alv.gif','DelayTime',1/15) %make gif file
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
%     text(50,10,['D = ', num2str(Dv(i))],'fontsize',18)
%     axis equal
%     grid on
%     
%     if i == 1
%         gif('ModelB2_grad.gif','DelayTime',1/15) %make gif file
%     else
%         gif %append frame to gif
%     end
%     disp(i)
% end
