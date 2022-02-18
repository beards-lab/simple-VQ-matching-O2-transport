function [JT, JV, JQ, vt, qt] = OBJECTIVE_A_v2(D,par,DATA)

par(1) = D;

%%% data
Q  = DATA.CO;
V  = DATA.V;
PO2 = DATA.PO2;

n = length(PO2); %number of data points

%%% find V as function of Q
opt = optimset('MaxFUnEvals',2000);

vt = zeros(n,1);
qt = zeros(n,1);
for i = 1:n
    vt(i) = fminsearch(@modA_eval,V(i),opt,Q(i),PO2(i),par,0);
    qt(i) = fminsearch(@modA_eval,Q(i),opt,V(i),PO2(i),par,1);
end

JV = sum((V-vt).^2);
JQ = sum((Q-qt).^2);

JT = JV+JQ;
