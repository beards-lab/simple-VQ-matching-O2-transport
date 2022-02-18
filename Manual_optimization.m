clear; close all; clc;
%%% brute forcing our way to find the optimal diffusion parameter

%%% load data
load('ExerciseData','DATA')

%%% look up table
load('Lookup.mat') %outputs LOOK
P = LOOK.Plookup;
C = LOOK.Clookup;

%%% parameters
D     = 100;      %apparent diffusion (L/min)
Pair  = 150;    %atmospheric oxygen partial pressure (mmHg)
Pin   = 45;     %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)
Vvasc = 1;      %volume of vascular space (ml)
Valv  = 1;      %alveolar volume (ml)

alpha = 1.3e-6*1e3;  % O2 solubility  in water/plasma(mM/mmHg)
CHb   = 0.021*1e3;   % Hb binding site conc (mmol/L of RBC's)
Hct   = 0.40;    % hematocrit (unitless)
C0    = CHb*Hct; % blood oxygen binding capacity (mol/L)
n     = 2.7;     % Hill exponent
P50   = 30;      % half-max saturation of Hb
beta  = 16800*1e-3; % O2 solubility in air (mmHg/mM)
l     = 1; %length of capillary

Vp = 70; %ventilation flow (L/min)
Qp = 40; %blood flow (L/min)
VO2M = 8; %O2 consumption rate (L/min)

k = 2.5;%conversion factor to convert O2 content (ml O2/100ml blood) to concentration (mM)

% empirical function of Cardiac Ouput as a function of work rate - from Stringer 1997
% Reference --> Stringer WW, Hansen JE, Wasserman K. Cardiac output estimated noninvasively from oxygen uptake during exercise. Journal of Applied Physiology. 1997 Mar 1;82(3):908-12.
a = 0.105;
b = 5.72;

%%%% flags
RFLAG = 0; %run flag (1 == run the stuff, 0 == don't run the stuff, but load and plot stuff)

if RFLAG == 0
    load('Manual_optimization_results.mat')

    %%% finding the values to plot
    DAp = 500; DBp = 6.1e3;%1e4;
    DCp = 1e2; DDp = 380;%436.54; 
    DEp = 1e4; DFp = 650; %best fit D
    
    JAp = interp1(DA, JA, DAp);
    JBp = interp1(DB, JB, DBp);
    JCp = interp1(DC, JC, DCp);
    JDp = interp1(DD, JD, DDp);%3.05e3;
    JEp = interp1(DE, JE, DEp);
    JFp = interp1(DF, JF, DFp);

    figure;
    subplot(2,3,1)
    loglog(DA,JA,DAp,JAp,'ro','linewidth',2,'markersize',7)
    title('Model A Cost')
    ylabel('J')
%     xlabel('D')
    set(gca,'fontsize',18)
    ylim([1e2 1e6])

    subplot(2,3,2)
    loglog(DB,JB,DBp,JBp,'ro','linewidth',2,'markersize',7)
    title('Model B Cost')
    %     ylabel('J')
    %     xlabel('D')
    set(gca,'fontsize',18)
    ylim([1e2 1e6])
    xlim([min(DB) max(DB)])

    subplot(2,3,3)
    loglog(DC,JC,DCp,JCp,'ro','linewidth',2,'markersize',7)
    title('Model C Cost')
    %     ylabel('J')
    %     xlabel('D')
    set(gca,'fontsize',18)
    ylim([1e2 1e6])

    subplot(2,3,4)
    loglog(DD,JD,DDp,JDp,'ro','linewidth',2,'markersize',7)
    title('Model D Cost')
    ylabel('J')
    xlabel('D')
    set(gca,'fontsize',18)
    ylim([1e2 1e6])
    xlim([min(DD) max(DD)])

    subplot(2,3,5)
    loglog(DE,JE,DEp,JEp,'ro','linewidth',2,'markersize',7)
    title('Model E Cost')
%     ylabel('J')
    xlabel('D')
    set(gca,'fontsize',18)
    ylim([1e2 1e6])

    subplot(2,3,6)
    loglog(DF,JF,DFp,JFp,'ro','linewidth',2,'markersize',7)
    title('Model F Cost')
%     ylabel('J')
    xlabel('D')
    set(gca,'fontsize',18)
    ylim([1e2 1e6])


elseif RFLAG == 1
    %%% model A
    DA = linspace(0,10000,1000);
    JA = zeros(length(DA),1);

    par = [D Pair Pin];
    for i = 1:length(DA)
        JA(i) = OBJECTIVE_A(DA(i),par,DATA);
    end
    figure;
    loglog(DA,JA,'linewidth',2)
    title('Model A Cost')
    ylabel('J')
    xlabel('D')
    set(gca,'fontsize',18)

    %%% model B
    DB = linspace(1e2,1e6,1000);
    JB = zeros(length(DB),1);

    par = [D Pair Pin Vvasc Valv alpha beta C0 n P50 Vp Qp];
    for i = 1:length(DB)
        JB(i) = OBJECTIVE_B(DB(i),par,DATA);
    end
    figure;
    loglog(DB,JB,'linewidth',2)
    title('Model B Cost')
    ylabel('J')
    xlabel('D')
    set(gca,'fontsize',18)

    %%% model C
    DC = linspace(10,1000,300);
    JC = zeros(length(DC),1);

    par = [D Pair Pin];
    for i = 1:length(DC)
        JC(i) = OBJECTIVE_C(DC(i),par,DATA);
    end
    figure;
    loglog(DC,JC,'linewidth',2)
    title('Model C Cost')
    ylabel('J')
    xlabel('D')
    set(gca,'fontsize',18)

    %%% model D
    DD = linspace(100,1000,200);
    JD = zeros(length(DD),1);

    par = [D Pair Pin alpha beta l];
    N = 300;
    for i = 1:length(DD)
        JD(i) = OBJECTIVE_D(DD(i),par,P,C,DATA,N);
        disp(i)
    end
    figure;
    loglog(DD,JD,'linewidth',2)
    title('Model D Cost')
    ylabel('J')
    xlabel('D')
    set(gca,'fontsize',18)

    %%% model E
    DE = linspace(0,1e5,1000);
    JE = zeros(length(DE),1);

    par = [Pair alpha beta a b k VO2M];
    for i = 1:length(DE)
        JE(i) = OBJECTIVE_E(DE(i),par,DATA,0);
    end
    figure;
    loglog(DE,JE,'linewidth',2)
    title('Model E Cost')
    ylabel('J')
    xlabel('D')
    set(gca,'fontsize',18)

    %%% model F
    DF = linspace(100,1000,100);
    JF = zeros(length(DF),1);

    par = [D Pair Pin alpha beta l a b k VO2M];
    N = 1000;
    for i = 1:length(DF)
        JF(i) = OBJECTIVE_F(DF(i),par,P,C,DATA,N);
        disp(i)
    end
    figure;
    loglog(DF,JF,'linewidth',2)
    title('Model F Cost')
    ylabel('J')
    xlabel('D')
    set(gca,'fontsize',18)

    save('Manual_optimization_results.mat')
end





