clear; close all; clc;
% This script uses data/information taken from a few sources to convert
% work rate to cardiac ouput. The ultimate purpose here is to prepare a
% data set from a CPET study to compare to our suite of simple pulmonary
% oxygen transport models.

%%% empirical function of Cardiac Ouput as a function of work rate - from Stringer 1997
% Reference --> Stringer WW, Hansen JE, Wasserman K. Cardiac output estimated noninvasively from oxygen uptake during exercise. Journal of Applied Physiology. 1997 Mar 1;82(3):908-12.
a = 0.105;
b = 5.72;

% plot of %VO2max and cardiac output from Stringer er al 
VM = 5;
vp = 0:100;
figure;
subplot(2,2,1)
plot(vp,(VM*vp)./(a*vp+b),'linewidth',2)
xlabel('%V_{O2,Max}')
ylabel('Cardiac Output (L/min)')
set(gca,'fontsize',18)

%%% Data from Riley 2000 - VO2 (L/min) as a funtion of work (Watts)
% Reference --> Riley MS, Pórszász J, Engelen MP, Brundage BH, Wasserman K. Gas exchange responses to continuous incremental cycle ergometry exercise in primary pulmonary hypertension in humans. European journal of applied physiology. 2000 Sep;83(1):63-70.
DATA1 = [10.176991150442483, 0.5978636615044246
7.227138643067853, 0.8575716721976401
22.566371681415944, 1.0590431415929202
12.24188790560472, 0.802982116519174
15.78171091445428, 0.7991450036873156
19.026548672566378, 0.8265521294247786
25.516224188790574, 0.6821476308997052
35.25073746312685, 0.6139783831120944
28.466076696165203, 0.8697052452064895
32.005899705014755, 0.9010243823746311
42.920353982300895, 0.8934250553097345
46.46017699115045, 0.8934941924778761
38.79056047197642, 0.9929537702802358
48.5250737463127, 1.0322063974926254
52.654867256637175, 1.1045526825221237
55.60471976401182, 1.0499227968289084
58.8495575221239, 0.9503767975663715
62.38938053097347, 1.2180240597345133
65.3392330383481, 1.1516754240412979
68.28908554572273, 1.1224361633480824
72.71386430678467, 1.1166632098082594
82.15339233038351, 1.1617694505899703
75.36873156342185, 1.2788244376843658
78.02359882005902, 1.3062200405604718
84.21828908554573, 1.2653254056047196
92.4778761061947, 1.2400961006637168
95.13274336283187, 1.2850698285398232
99.5575221238938, 1.33984375
101.62241887905607, 1.3769934550147491
105.45722713864309, 1.400505853613569
88.64306784660768, 1.4040837020648969
111.35693215339235, 1.5041367072271385
114.89675516224192, 1.5276433443952804
124.33628318584074, 1.5082964601769913
107.52212389380533, 1.5743743086283186
98.0825958702065, 1.5663774428466075
117.55162241887906, 1.6233983222713864
120.50147492625372, 1.629315311578171
127.58112094395283, 1.594297335914454
131.1209439528024, 1.6490539730825957
134.36578171091446, 1.697945473820059
142.0353982300885, 1.647314021017699
137.0206489675516, 1.7741692016961652
147.6401179941003, 1.78023598820059
143.21533923303838, 1.778196441740413
149.70501474926255, 1.7490263182153392
154.424778761062, 1.8077122511061947
160.91445427728615, 1.909401502581121
158.55457227138646, 1.8995897861356932
164.45427728613572, 2.0520487647492627
180.67846607669622, 2.0816625184365782
170.94395280235995, 1.9935817662241888
173.5988200589971, 1.888164869100295
167.69911504424783, 1.8548465154867255
177.13864306784666, 1.7983902562684366
183.62831858407083, 1.8239076327433628];

W   = DATA1(:,1);
VO2 = DATA1(:,2);

[F,D] = fit(W,VO2,'poly1'); %linear regression
c = F.p1; d = F.p2;         %extracting the regression parameters

% plot of work rate and VO2 data from Riley et al + regression line
wp1 = 0:1:200;
subplot(2,2,2)
plot(W,VO2,'o')
hold on
plot(wp1,c*wp1+d,'k')
xlabel('Work (Watts)')
ylabel('V_{O2} (L/min)')
set(gca,'fontsize',18)

%%% plotting Cardiac Ouput as a function of work rate by combinining info from Riley et al and Stringer et al.
par = [a b c d];
subplot(2,2,3)
plot(wp1,CO_W(wp1,4.5,par),'linewidth',2)
hold on
plot(wp1,CO_W(wp1,5.5,par),'linewidth',2)
hold on
plot(wp1,CO_W(wp1,8,par),'linewidth',2)
xlabel('Work (Watts)')
ylabel('Cardiac Output (L/min)')
set(gca,'fontsize',18)
legend('V_{O2,Max} = 4.5', 'V_{O2,Max} = 5.5', 'V_{O2,Max} = 8','location','northwest')

%%% loading exercise data from ODonnell 2017
% Reference --> O'Donnell DE, Elbehairy AF, Berton DC, Domnik NJ, Neder JA. Advances in the evaluation of respiratory pathophysiology during exercise in chronic lung diseases. Frontiers in physiology. 2017 Feb 22;8:82.
DATA2 = [0, 9.0032154340836;
10, 9.0032154340836;
20, 10.932475884244369;
29.83870967741936, 13.922829581993568;
40.161290322580655, 15.080385852090032;
60, 17.97427652733119;
80.00000000000003, 28.102893890675247;
100.00000000000003, 36.01286173633441;
120.16129032258067, 44.88745980707396;
140.1612903225807, 65.91639871382638];

DATA3 = [-7.105427357601002e-15, 98;
10, 97.97333333333333;
20, 98;
30.000000000000007, 98;
40.000000000000014, 98.98666666666666;
60.000000000000014, 98;
80.15625000000003, 98.98666666666666;
100.00000000000001, 98;
120.31250000000004, 98;
140.15625000000006, 98];

W_data = DATA2(:,1);
V_data = DATA2(:,2);
SO2_data = DATA3(:,2)/100;

%%% converting watts to cardiac output (L/min)
CO_data = CO_W(W_data,4.5,par);

%%% converting SO2 to partial pressure
n   = 2.7;     % Hill exponent
P50 = 27;      % half-max saturation of Hb

PO2_data = P50./nthroot((1./SO2_data)-1,n);

PO2_data = 110+0*PO2_data;

%%% stacking things into a nicely packaged struct
DATA.PO2 = PO2_data;
DATA.CO  = CO_data;
DATA.V   = V_data;

%%% ploting cardiac output and ventilation data from O'Donnell et al.
subplot(2,2,4)
plot(CO_data,V_data,'ko','linewidth',2,'markersize',10)
hold on
plot(0:70,0:70,'k--','linewidth',2)
xlabel('Cardiac Output (L/min)')
ylabel('Ventilation (L/min)')
set(gca,'fontsize',18)
xlim([5 20])
ylim([5, 70])

% figure;
% plot(CO_data,V_data,'ko','linewidth',2,'markersize',10)
% hold on
% plot(0:70,0:70,'k--','linewidth',2)
% xlabel('Cardiac Output (L/min)')
% ylabel('Ventilation (L/min)')
% set(gca,'fontsize',18)
% xlim([8 18])

% saving all of our hard work into a .mat file to use in other scripts
save('ExerciseData.mat','PO2_data','CO_data','V_data','DATA')

function Q = CO_W(W,VM,par) %function that converts work rate to cardiac output
    %VM is the maximum oxygen consumption rate
    a = par(1); b = par(2); c = par(3); d = par(4); %unpack empirical regression parameters
       
    VO2           = c*W+d;                       %VO2 as a function work rate
    VO2(VO2 > VM) = VM;                          %bound VO2 based on the maximum VO2
    Q             = 100*VO2./(a*(100*VO2/VM)+b); %cardiac output as a function of VO2
end

