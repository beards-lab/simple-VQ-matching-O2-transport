function OUT = ModelB2_FixedPoint_Objective(X,par)
%%% this function defines the algebraic relationship to dolve for the the
%%% fixed points in ModelB - intended to be numerically solved with fsolve

%Paramaters
D     = par(1); %apparent diffusion (L/min)
Pair  = par(2); %atmospheric oxygen (mmHg)
% Pin   = par(3); %inlet oxygen (mmHg)
alpha = par(6);
beta  = par(7);
C0    = par(8);
n     = par(9);
P50   = par(10);

V = par(11); %airflow (L/min) 
Q = par(12); %blood flow (L/min)
VO2 = par(13); %O2 consumption (L/min)


%%% variables
Pv = X(1);
Pa = X(2);

%%% O2 interpolation for hemoglobin solubility
Cv = alpha*Pv+C0*((Pv^n)/(Pv^n+P50^n));

%%% circulatory loop
k = 2.5; %conversion factor for mM to ml/100 ml

Cin = Cv - 1e3*VO2/(k*Q);
if Cin < 0
   Cin = 0;
%    disp('ded')
end

%%% function to minimize
pvm = Q*(Cin - Cv) + alpha*D*(Pa - Pv);
pam = V*(Pair - Pa) - alpha*beta*D*(Pa - Pv);

OUT = [pvm; pam];