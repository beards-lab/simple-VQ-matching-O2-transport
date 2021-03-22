function [J, Pv_m, r] = OBJECTIVE_C(D,par,DATA)
%%% this function definies the cost function to fit Model C to the exercise
%%% pacing data

% parameters
Pair = par(2);
Pin  = par(3);

% data
Q = DATA.CO;
V = DATA.V;

% calculate model vascular oxygen
n = length(DATA.CO);
Pv_m = zeros(n,1);
for i = 1:n
    clear palv
    palv    = (Pin*Q(i).*(exp(-D./Q(i))-1)-Pair*V(i))./(Q(i).*(exp(-D./Q(i))-1)-V(i));
    Pv_m(i) = palv + (Pin-palv).*exp(-D./Q(i));
end

% calculate residuals and cost
r = DATA.PO2 - Pv_m;
J = r'*r;