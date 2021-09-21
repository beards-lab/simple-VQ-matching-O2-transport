function [J, Pv_m, r] = OBJECTIVE_E(epar,par,DATA,FLAG)

Pair  = par(1);
alpha = par(2);
beta  = par(3);
a     = par(4);
b     = par(5);
k     = par(6);
VO2M  = par(7);

Q = DATA.CO;
V = DATA.V;

if FLAG == 0
    Dn = epar(1);
elseif FLAG == 1
    Dn = epar(1)*Q;
end

VO2 = -VO2M.*b.*Q./(100.*(a.*Q-VO2M));
Pv_m = Pair - VO2.*(beta./V + 1/alpha./Dn)./k;

r = DATA.PO2 - Pv_m;
J = r'*r;