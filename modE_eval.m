function J = modE_eval(vq1,vq2,P,par,FLAG,D)

if FLAG == 0 %vertical errors
    v = vq1; q = vq2;
elseif FLAG == 1 %horizontal errors
    q = vq1; v = vq2;
end

if v < 0 || q < 0
    J = 1e12;
else
    %%% calcualte vascular oxygen
    Pair  = par(1);
    alpha = par(2);
    beta  = par(3);
    a     = par(4);
    b     = par(5);
    k     = par(6);
    VO2M  = par(7);


    VO2 = -VO2M.*b.*q./(a.*q-VO2M);
    Pv_m = Pair - VO2.*(beta./v + 1/alpha./D)./k;

    %%% error to minimize
    J = (P-Pv_m).^2;
end


