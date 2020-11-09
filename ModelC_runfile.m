clear; close all; clc;

%%% parameters
D     = 1;   %apparent diffusion (ml/s)
Pair  = 150; %atmospheric oxygen partial pressure (mmHg)
Pin   = 45;  %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)
Vvasc = 1;   %volume of vascular space (ml)
Valv  = 1;   %alveolar volume (ml)

Vp = 5; %ventilation flow (ml/s)
Qp = 5; %blood flow (ml/s)

par = [D Pair Pin Vvasc Valv Vp  Qp]; %parameters vector for ODE system

%%% grid for contour plot
de = 0.05;
q  = 0:de:10;
v  = 0:de:10;
[Q,V] = meshgrid(q,v);

% computing contour plot
palv  = (Pin*Q.*(1-exp(-D./Q))-V.*Pair)./(Q.*(1-exp(-D./Q))-V-2*D);
pvasc = palv + (Pin-palv).*exp(-D./Q);

dpac  = palv - pvasc;

%%% solve system of ODEs
% X0 = [Pin Pair];
% [t,X] = ode45(@ModelC_RHS,[0 10],X0,[],par);

%%% plots
figure;
surf(Q,V,pvasc, 'EdgeColor','none')
title('Pvasc')
xlabel('q')
ylabel('v')
% set(gca,'zscale','log')
% axis equal

figure;
surf(Q,V,palv, 'EdgeColor','none')
title('Palv')
xlabel('q')
ylabel('v')
% set(gca,'zscale','log')
% axis equal

figure;
surf(Q,V,dpac, 'EdgeColor','none')
title('\Delta P')
xlabel('q')
ylabel('v')
% set(gca,'zscale','log')
% axis equal

% figure;
% plot(t,X)


CON = -30:5:150;
x = 0:10;
figure;
subplot(1,3,1)
contour(Q,V,pvasc,CON,'ShowText','on')
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
contour(Q,V,palv,CON,'ShowText','on')
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
contour(Q,V,dpac,CON,'ShowText','on')
hold on
plot(x,x,'k--')
hold on
plot(5,5,'rx')
set(gca,'fontsize',18)
title('Gradient (mmHg)')
xlabel('Cardiac Output (ml/s)')
axis equal
grid on


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

