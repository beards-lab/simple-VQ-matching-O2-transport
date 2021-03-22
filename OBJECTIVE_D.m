function [J, Pv_m] = OBJECTIVE_D(D,par,P,C,DATA,N)
%%% this function definies the cost function to fit Model D to the exercise
%%% pacing data

% set diffusion into the parameter vector
par(1) = D;

% calculate model D vascular oxygen tension
n = length(DATA.CO);
Pv_m = zeros(n,1);
parfor i = 1:n
    Pv_m(i) = modelD_SS_relaxation(N,par,P,C,DATA.V(i),DATA.CO(i));
end

% calcualte residuals and cost
r = DATA.PO2 - Pv_m;
J = r'*r;