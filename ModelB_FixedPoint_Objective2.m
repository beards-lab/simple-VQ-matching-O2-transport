function [OUT c]     = ModelB_FixedPoint_Objective2(X,par)
%%% this function defines the algebraic relationship to dolve for the the
%%% fixed points in ModelB - intended to be numerically solved with fsolve

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

%%% variables
Pv = X(1);
Pa = X(2);

%%% O2 interpolation for hemoglobin solubility
Cin = alpha*Pin+C0*((Pin^n)/(Pin^n+P50^n));
Cv = alpha*Pv+C0*((Pv^n)/(Pv^n+P50^n));
% Ca = alpha*Pa+C0*((Pa^n)/(Pa^n+P50^n));

%%% function to minimize
pvm = Q*(Cin - Cv) + alpha*D*(Pa - Pv);
pam = V*(Pair - Pa) - alpha*beta*D*(Pa - Pv);

OUT = [pvm; pam];
% c = [Cin, Cv, Ca];
c = [Cin, Cv];
