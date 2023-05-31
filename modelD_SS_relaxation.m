function [Pv,Pa,Cvi,Pvi,EPSv,dx] = modelD_SS_relaxation(N,par,P,C,vp,qp)
% this function solves the steady state solutions of Model D by
% numerically relaxing the system.
% Usage:
%   modelD_SS_relaxation(N,par,P,C,Vp,Q
% 
% where   
% N - number of discretization elements
% par = [D Pair Pin alpha beta l M];
% P and C - O2 partial pressure and O2 concentration look-up table
% vp - ventilation (L/min/kg)
% vp - perfusion (L/min/kg)

%%% parameters
D     = par(1); % apparent diffusion coefficient L/min/kg
Pair  = par(2); % Air partial pressure (mmHg)
Pin   = par(3); % Systemic venous oxygen partial pressure - lung input (mmHg)
alpha = par(4); % O2 solubility  in water/plasma(mM/mmHg)
beta  = par(5); % O2 solubility in air (mmHg/mM)
l     = par(6); % length of the capillary (m) (value may be contained in D already)
% M     = par(7); % mass of the compartment (kg) (value may be contained in D already)
maxIter = 10000;


%%% set up stuff
dx = l/N; %spatial step size (m)
Cin = interp1(P,C,Pin); %inflow O2 concentration (mM/L)
Cvi = Cin*ones(N,1); % initialize O2 conc along the capillary

%%% while loop
EPS   = 1; %arbitrary initialization of error
EPSv  = EPS; %the eventual vector to store the error
for i = 1:maxIter 
    %iterating to update? not sure how many iterations will be needed
    clear Pvi Pa Pvm dCv Cvt EPS
    Pvi = interp1(C,P,Cvi); %O2 pp along the capillary
    Pvm = mean(Pvi); %average O2 pp along the capillary
    % PaT = (vp*M*Pair+alpha*beta*D*M*Pvm)/(vp*M+alpha*beta*D*M); % alveolar O2 (mmHg)
    PaT = (vp*Pair+alpha*beta*D*Pvm)/(vp+alpha*beta*D); % alveolar O2 (mmHg)
    try
    dCv = dx*alpha*D*(PaT-Pvi)/(qp*l); % change in Cv
    catch
        a=1;
    end
    Cvt        = zeros(N,1); %temporary vector for Cvi
    Cvt(2:end) = Cvi(1:end-1)+dCv(1:end-1); %apply iterative update
    Cvt(1)     = Cin; %reapplay inflow O2 BC
    
    EPS  = abs((Cvt(end)-Cvi(end))/(Cvi(end))); % percent change (error) in O2 leaving the capillary
    EPSv = [EPSv EPS]; %stack error into a vector
    
    Cvi = Cvt; %update Cvi vector

    if EPS < 1e-12 || isnan(EPS)
        % either we are under tolerance or in deep nans, aborting
        break;
    end
end
if i == maxIter
    % iteration not successful
    Cvi(end) = NaN;
end
EPSv(1) = []; %delete first entry, as this is just 1 and not a real result

Pvi = interp1(C,P,Cvi); %O2 pp along the capillary 
Pv  = Pvi(end); %steady state vascular O2

Pa = (vp*Pair+alpha*beta*D*mean(Pvi))/(vp+alpha*beta*D); %steady state alveolar O2
