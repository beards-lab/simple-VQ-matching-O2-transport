function [J, Pv_m, r] = OBJECTIVE_B(D,par,DATA)
opt = optimset('MaxFUnEvals',2000);
par(1) = D;

Q = DATA.CO;
V = DATA.V;


n = length(DATA.CO);
Pv_m = zeros(n,1);
for i = 1:n
    clear tpar PSIM
    tpar = par;
    tpar(end)   = Q(i);
    tpar(end-1) = V(i);
    
    PSIM = fsolve(@ModelB_FixedPoint_Objective2,[100 100],opt,tpar);
    Pv_m(i) = PSIM(1);
end

r = DATA.PO2 - Pv_m;
J = r'*r;