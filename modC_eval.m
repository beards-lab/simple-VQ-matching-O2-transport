function J = modC_eval(vq1,vq2,P,par,FLAG)

if FLAG == 0 %vertical errors
    v = vq1; q = vq2;
elseif FLAG == 1 %horizontal errors
    q = vq1; v = vq2;
end

if v < 0 || q < 0
    J = 1e12;
else
    %%% calcualte vascular oxygen
    D = par(1);
    Pair = par(2);
    Pin  = par(3);

    palv    = (Pin*q.*(exp(-D./q)-1)-Pair*v)./(q.*(exp(-D./q)-1)-v);
    Pv_m = palv + (Pin-palv).*exp(-D./q);

    %%% error to minimize
    J = (P-Pv_m).^2;
end


