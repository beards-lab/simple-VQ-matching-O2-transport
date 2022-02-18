function [JT, JV, JQ, vt, qt] = OBJECTIVE_D_v2(D,par,DATA,Plook,Clook,N)

par(1) = D;

%%% data
Q  = DATA.CO;
V  = DATA.V;
PO2 = DATA.PO2;

n = length(PO2); %number of data points

%%% janky optimization wuz gud
opt = optimset('MaxFUnEvals',2000);

vt = zeros(n,1);
qt = zeros(n,1);
parfor i = 1:n
    vt(i) = fminsearch(@modD_eval,V(i),opt,Q(i),PO2(i),par,0,Plook,Clook,N);
end
parfor j = 1:n
    qt(j) = fminsearch(@modD_eval,Q(j),opt,V(j),PO2(j),par,1,Plook,Clook,N);
end

JV = sum((V-vt).^2);
JQ = sum((Q-qt).^2);

JT = JV+JQ;
