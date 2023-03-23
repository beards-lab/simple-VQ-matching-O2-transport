clear; close all; clc;

%%% load exercise data
load('ExerciseData','CO_data','V_data')

%%% model D parameters
Pair  = 150;    %atmospheric oxygen partial pressure (mmHg)
Pin   = 45;     %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)
Vvasc = 1;      %volume of vascular space (ml)
Valv  = 1;      %alveolar volume (ml)

alpha = 1.3e-6*1e3;  % O2 solubility  in water/plasma(mM/mmHg)
beta  = 16800*1e-3; % O2 solubility in air (mmHg/mM)
l     = 1; %length of capillary

Vp = 5; %ventilation flow (L/min)
Qp = 5; %blood flow (L/min)

%load optimized diffusion (D) parameter
load('ModelD_optimization_v2_results.mat','JD','DD')
[~, jDpi] = min(JD); DDp = DD(jDpi);
D = DDp;

par = [D Pair Pin alpha beta l];

%%% look up table
load('Lookup.mat') %outputs LOOK
P = LOOK.Plookup;
C = LOOK.Clookup;

N = 300; %number of grid points for numerical discretization


%%%%%% find empirical estimates for the coefficients for the optimal V/Q matching curve.
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0, 0],...
               'Upper',[Inf,Inf, Inf],...
               'StartPoint',100*[1 1 1]);
ft = fittype('a*x/(b*x+c)','options',fo); %this is the operative definition of the "optimal curve"
[opt_curve,metric] = fit(V_data,CO_data,ft); %this line does the fitting

%%% Set up ventilation and perfusion distributison 
n = 3;
V_test = linspace(5,50,n);%5:10:70; %grid of sample ventilations

Q_25  = mean(V_test)*ones(n,1); % uniform blood flows
Q_per = V_test;% perfect V/Q matching
Q_opt = opt_curve(V_test); % blood flow predicted by the optimality condition

%%%% Simulate blood flows using model D
P_per = zeros(n,1); %preallocate
P_25  = zeros(n,1);
P_opt = zeros(n,1);

for i = 1:n
    P_25(i)  = modelD_SS_relaxation(N,par,P,C,V_test(i),Q_25(i));
    P_per(i) = modelD_SS_relaxation(N,par,P,C,V_test(i),Q_per(i));
    P_opt(i) = modelD_SS_relaxation(N,par,P,C,V_test(i),Q_opt(i));
end

%%% calcualte total O2 in outlfow
C_30  = interp1(P,C,P_25);
C_per = interp1(P,C,P_per); %O2 conc (mM) per sample ventilation
C_opt = interp1(P,C,P_opt);

C_tot_25  = sum(C_30.*Q_25)/(sum(Q_25));
C_tot_per = sum(C_per.*Q_per')/(sum(Q_per)); %total O2 conc in washout (mM)
C_tot_opt = sum(C_opt.*Q_opt)/(sum(Q_opt));

O2_25  = interp1(C,P,C_tot_25);
O2_per = interp1(C,P,C_tot_per); %total O2 tension (mnmHg) in washout 
O2_opt = interp1(C,P,C_tot_opt);

save('optimality_test_results.mat')

%%% plotting
load('optimality_test_results.mat')
close all;

lw = 2; %linewidth
ms = 10; %markersize

figure;
subplot(2,2,1)
plot(V_test,Q_25,'bo','linewidth',lw,'markersize',ms)
hold on
plot(V_test,Q_per,'rd','linewidth',lw,'markersize',ms)
hold on
plot(V_test,Q_opt,'ks','linewidth',lw,'markersize',ms)
legend('Uniform', 'V=Q','V=AQ/(B-CQ)','location','northwest')
set(gca,'fontsize',20)
ylabel('Blood Flow (L/min)')
ylim([0 max(V_test)+5])
xlim([0 max(V_test)+5])

% figure;
subplot(2,2,3)
plot(V_test,P_25,'bo','linewidth',lw,'markersize',ms)
hold on
plot(V_test,P_per,'rd','linewidth',lw,'markersize',ms)
hold on
plot(V_test,P_opt,'ks','linewidth',lw,'markersize',ms)
% legend('Uniform - 15', 'Uniform - 30','Optimal','location','best')
set(gca,'fontsize',20)
xlabel('Ventilation (L/min)')
ylabel('Oxygen Tension Outflow (mmHg)')
xlim([0 max(V_test)+5])
ylim([0 130])

% figure;
M = 2; %multiplier to scale total maker size
subplot(2,2,2)
plot(1,O2_25,'bo','linewidth',lw,'markersize',M*ms)
hold on
plot(2,O2_per,'rd','linewidth',lw,'markersize',M*ms)
hold on
plot(3,O2_opt,'ks','linewidth',lw,'markersize',M*ms)
set(gca,'fontsize',20)
xticks([1 2 3])
xticklabels({'Uniform','V=Q','V=AQ/(B-CQ)'})
xtickangle(0)
ylabel('Total Oxygen Tension (mmHg)')
% legend('Uniform', 'V=Q','V=AQ/(B-CQ)', 'location','northwest')
xlim([0.5 3.5])
ylim([0 120])


% figure;
subplot(2,2,4)
plot(1,C_tot_25,'bo','linewidth',lw,'markersize',M*ms)
hold on
plot(2,C_tot_per,'rd','linewidth',lw,'markersize',M*ms)
hold on
plot(3,C_tot_opt,'ks','linewidth',lw,'markersize',M*ms)
set(gca,'fontsize',20)
xticks([1 2 3])
xticklabels({'Uniform','V=Q','V=AQ/(B-CQ)'})
xtickangle(0)
ylabel('Total Oxygen Conc. (mM)')
%legend('Uniform', 'V=Q','V=AQ/(B-CQ)', 'location','northwest')
xlim([0.5 3.5])
ylim([5 9])



