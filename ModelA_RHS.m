function dX = ModelA_RHS(~,X,par)

%%% parameters
D     = par(1); %apparent diffusion (ml/s)
Pair  = par(2); %atmospheric oxygen (ml/s)
Pin   = par(3); %inlet oxygen (ml/s)
Vvasc = par(4); %vascular volume (ml)
Valv  = par(5); %alveolar volume (ml)

V = par(6); %airflow (ml/s)
Q = par(7); %blood flow (ml/s)

%%% states
Pvasc = X(1); %vascular oxygen (mmHg)
Palv  = X(2); %alveolar oxygen (mmHg)

%%% RHS
dPvasc = (Q*(Pin - Pvasc)+D*(Palv - Pvasc))./Vvasc;
dPalv  = (V*(Pair - Palv)-D*(Palv - Pvasc))./Valv;

dX = [dPvasc; dPalv];
