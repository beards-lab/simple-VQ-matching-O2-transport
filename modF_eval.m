function J = modF_eval(vq1,vq2,P,par,FLAG,D,Plook,Clook,N)

if FLAG == 0 %vertical errors
    v = vq1; q = vq2;
elseif FLAG == 1 %horizontal errors
    q = vq1; v = vq2;
end

if v < 0 || q < 0
    J = 1e12;
else
    %%% calcualte vascular oxygen
    par(1) = D;
    Pv_m = modelF_SS_relaxation(N,par,Plook,Clook,v,q);

    %%% error to minimize
    J = (P-Pv_m).^2;
end


