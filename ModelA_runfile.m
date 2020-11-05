clear; close all; clc;

%%% parameters
D     = 5;   %apparent diffusion (ml/s)
Pair  = 150; %atmospheric oxygen partial pressure (mmHg)
Pin   = 75;  %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)
Vvasc = 1;   %volume of vascular space (ml)
Valv  = 1;   %alveolar volume (ml)

Vp = 5; %ventilation flow (ml/s)
Qp = 5; %blood flow (ml/s)

par = [D Pair Pin Vvasc Valv Vp  Qp]; %parameters vector for ODE system

%%% grid for contour plot
de = 0.1;
q  = 0:de:10;
v  = 0:de:10;
[Q,V] = meshgrid(q,v);

% computing contour plot
pvasc = (Pin*Q.*V+Pin*D*Q+Pair*D*V)./(Q.*V+D*Q+D*V);
palv  = (Pair*Q.*V+Pair*D*V+Pin*D*Q)./(Q.*V+D*Q+D*V);
dpac  = (Q.*V)*(Pair-Pin)./(Q.*V+D*Q+D*V);

%%% solve system of ODEs
X0 = [Pin Pair];
[t,X] = ode45(@ModelA_RHS,[0 10],X0,[],par);

%%% plots
figure;
surf(Q,V,pvasc, 'EdgeColor','none')
title('Pvasc')
xlabel('q')
ylabel('v')
axis equal

figure;
surf(Q,V,palv, 'EdgeColor','none')
title('Palv')
xlabel('q')
ylabel('v')
axis equal

figure;
surf(Q,V,dpac, 'EdgeColor','none')
title('\Delta P')
xlabel('q')
ylabel('v')
axis equal

figure;
plot(t,X)

x = 0:10;
figure;
subplot(1,3,1)
contour(Q,V,pvasc,'ShowText','on')
hold on
plot(x,x,'k--')
hold on
plot(5,5,'rx')
set(gca,'fontsize',18)
title('Vascular (mmHg)')
xlabel('Cardiac Output (ml/s)')
ylabel('Ventilation Magnitude (ml/s)')
axis equal
grid on

% figure;
subplot(1,3,2)
contour(Q,V,palv,'ShowText','on')
hold on
plot(x,x,'k--')
hold on
plot(5,5,'rx')
set(gca,'fontsize',18)
title('Alveolar Space (mmHg)')
xlabel('Cardiac Output (ml/s)')
axis equal
grid on

% figure;
subplot(1,3,3)
contour(Q,V,dpac,'ShowText','on')
hold on
plot(x,x,'k--')
hold on
plot(5,5,'rx')
set(gca,'fontsize',18)
title('Gradient (mmHg)')
xlabel('Cardiac Output (ml/s)')
axis equal
grid on

