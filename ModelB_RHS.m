function dX = ModelB_RHS(~,X,par,LOOK)
% this function defines the kinetic right hand side of Model B to solved
% with an ode solver

%Paramaters
D     = par(1); %apparent diffusion (ml/s)
Pair  = par(2); %atmospheric oxygen (ml/s)
Pin   = par(3); %inlet oxygen (ml/s)
Vvasc = par(4); %vascular volume (ml)
Valv  = par(5); %alveolar volume (ml)
alpha = par(6);
beta  = par(7);

V = par(11); %airflow (ml/s)
Q = par(12); %blood flow (ml/s)

%%% states
Cvasc = X(1); %vascular oxygen (mM)
Palv  = X(2); %alveolar oxygen (mmol)

%%% O2 interpolation for hemoglobin solubility
P = LOOK.Plookup;
C = LOOK.Clookup;
Cin   = interp1(P,C,Pin);
Pvasc = interp1(C,P,Cvasc);

%%% RHS
dCvasc = (Q*(Cin - Cvasc) + alpha*D*(Palv - Pvasc))./Vvasc;
dPalv  = (V*(Pair-Palv) - alpha*beta*D*(Palv - Pvasc))./Valv;

dX = [dCvasc; dPalv];
