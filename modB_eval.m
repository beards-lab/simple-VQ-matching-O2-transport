function J = modB_eval(vq1,vq2,P,par,FLAG)

if FLAG == 0 %vertical errors
    v = vq1; q = vq2;
elseif FLAG == 1 %horizontal errors
    q = vq1; v = vq2;
end

if v < 0 || q < 0
    J = 1e12;
else
    %%% calcualte vascular oxygen
    par(11) = v; %airflow (ml/s)
    par(12) = q; %blood flow (ml/s)

    opt = optimset('MaxFUnEvals',2000);
    Pg = fsolve(@ModelB_FixedPoint_Objective2,[100 100],opt,par);
    Pv_m = Pg(1);

    %%% error to minimize
    J = (P-Pv_m).^2;
end


