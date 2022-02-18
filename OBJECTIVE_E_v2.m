function [JT, JV, JQ, vt, qt] = OBJECTIVE_E_v2(D,par,DATA)

%%% data
Q  = DATA.CO;
V  = DATA.V;
PO2 = DATA.PO2;

n = length(PO2); %number of data points

%%% janky optimization wuz gud
opt = optimset('MaxFUnEvals',2000);

vt = zeros(n,1);
qt = zeros(n,1);
for i = 1:n
    vt(i) = fminsearch(@modE_eval,V(i),opt,Q(i),PO2(i),par,0,D);
    qt(i) = fminsearch(@modE_eval,Q(i),opt,V(i),PO2(i),par,1,D);
end

JV = sum((V-vt).^2);
JQ = sum((Q-qt).^2);

JT = JV+JQ;
