function [Pv,Pa,Cvi,Pvi,EPSv,dx] = modelD2_SS_relaxation(N,par,P,C,Vp,Qp)
%%% this function solves the steady state solutions of Model D2 by
%%% numerically relaxing the system.

%%% parameters
D     = par(1);
Pair  = par(2);
Pin   = par(3);
alpha = par(4);
beta  = par(5);
l     = par(6);
a     = par(7);
b     = par(8);
k     = par(9);
VO2M  = par(10);

%%% set up stuff
dx = l/N; %spatial step size
C0 = interp1(P,C,100); %inflow O2 concentration
Cvi = C0*ones(N,1); %initialize O2 conc along the capillary

%%% while loop
EPS   = 1; %arbitrary initialization of error
EPSv  = EPS; %the eventual vector to store the error
while EPS > 1e-8 %iterating to update? not sure how many iterations will be needed
    clear Pvi Pa Pvm dCv Cvt EPS
    Pvi = interp1(C,P,Cvi); %O2 pp along the capillary
    Pvm = mean(Pvi); %average O2 pp along the capillary
    PaT = (Vp*Pair+alpha*beta*D*Pvm)/(Vp+alpha*beta*D); %alveolar O2
    dCv = dx*alpha*D*(PaT-Pvi)./(Qp*l); %change in Cv
    
    
    VO2 = -VO2M*b*Qp./(100*(a*Qp-VO2M)); %O2 consumption
    Cin = Cvi(end) - VO2/Qp/k; %oxygen concentration entering the pulmonary circulation
    
    Cvt        = zeros(N,1); %temporary vector for Cvi
    Cvt(2:end) = Cvi(1:end-1)+dCv(1:end-1); %apply iterative update
    Cvt(1)     = Cin; %reapplay inflow O2 BC
    
    EPS  = abs((Cvt(end)-Cvi(end))/(Cvi(end))); % percent change (error) in O2 leaving the capillary
    EPSv = [EPSv EPS]; %stack error into a vector
    
    Cvi = Cvt; %update Cvi vector
end
EPSv(1) = []; %delete first entry, as this is just 1 and not a real result

Pvi = interp1(C,P,Cvi); %O2 pp along the capillary 
Pv  = Pvi(end); %steady state vascular O2

Pa = (Vp*Pair+alpha*beta*D*mean(Pvi))/(Vp+alpha*beta*D); %steady state alveolar O2
