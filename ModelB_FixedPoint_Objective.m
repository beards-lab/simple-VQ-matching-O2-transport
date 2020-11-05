function OUT = ModelB_FixedPoint_Objective(Pv,par,LOOK)
%%% this function defines the algebraic relationship to dolve for the the
%%% fixed points in ModelB - intended to be numerically solved

%Paramaters
D     = par(1); %apparent diffusion (ml/s)
Pair  = par(2); %atmospheric oxygen (ml/s)
Pin   = par(3); %inlet oxygen (ml/s)
alpha = par(6);
beta  = par(7);
C0    = par(8);
n     = par(9);
P50   = par(10);

V = par(11); %airflow (ml/s)
Q = par(12); %blood flow (ml/s)


%%% O2 interpolation for hemoglobin solubility
P = LOOK.Plookup;
C = LOOK.Clookup;
Cin = interp1(P,C,Pin);

OUT = Q*(Cin - alpha*Pv - C0*((Pv^n)/(Pv^n + P50^n))) + alpha*D*(((V*Pair + alpha*beta*D*Pv)/(V + alpha*beta*D)) - Pv);