function [J, Pv_m, r] = OBJECTIVE_A(D,par,DATA)
%%% this function definies the cost function to fit Model A to the exercise
%%% pacing data

%%% parameters
Pair = par(2);
Pin  = par(3);

%%% data
Q = DATA.CO;
V = DATA.V;

%%% calcualte vascular oxygen
Pv_m = (Pin*Q.*V+Pin*D*Q+Pair*D*V)./(Q.*V+D*Q+D*V);

%%% residuals and cost function
r = DATA.PO2 - Pv_m;
J = r'*r;