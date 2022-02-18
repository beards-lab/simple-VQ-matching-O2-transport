function J = modA_eval(vq1,vq2,P,par,FLAG)

if FLAG == 0 %vertical errors
    v = vq1; q = vq2;
elseif FLAG == 1 %horizontal errors
    q = vq1; v = vq2;
end

if v < 0 || q < 0
    J = 1e12;
else
    %%% parameters
    D    = par(1);
    Pair = par(2);
    Pin  = par(3);

    %%% calcualte vascular oxygen
    Pv_m = (Pin*q.*v+Pin*D*q+Pair*D*v)./(q.*v+D*q+D*v);

    %%% error to minimize
    J = (P-Pv_m).^2;
end


