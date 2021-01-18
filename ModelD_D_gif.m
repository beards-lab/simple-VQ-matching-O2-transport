clear; close all; clc;

%%% parameters
D     = 50;      %apparent diffusion (ml/s)
Pair  = 150;    %atmospheric oxygen partial pressure (mmHg)
Pin   = 45;     %mixed venous oxygen partial pressure - pulmonary inlet (mmHg)
Vvasc = 1;      %volume of vascular space (ml)
Valv  = 1;      %alveolar volume (ml)

alpha = 1.3e-6*1e3;  % O2 solubility  in water/plasma(mM/mmHg)
CHb   = 0.021*1e3;   % Hb binding site conc (mmol/L of RBC's)
Hct   = 0.40;    % hematocrit (unitless)
C0    = CHb*Hct; % blood oxygen binding capacity (mol/L)
n     = 2.7;     % Hill exponent
P50   = 30;      % half-max saturation of Hb
beta  = 16800*1e-3; % O2 solubility in air (mmHg/mM)
l     = 1; %length of capillary

Vp = 5; %ventilation flow (ml/s)
Qp = 5; %blood flow (ml/s)

par = [D Pair Pin alpha beta l];

%%% look up table
load('Lookup.mat') %outputs LOOK
P = LOOK.Plookup;
C = LOOK.Clookup;

N = 300; %number of grid points for numerical discretization

%%% grid for contour plot
de = 0.1;
q  = 0:de:10;
v  = 0:de:10;
[Q,V] = meshgrid(q,v);
np = length(q); %number of points for air and blood flow to use

Dv = 70:10:150; %vector of diffusion parameters to explore
nD = length(Dv);

%%% make different vectors for diffusion
Dpar{nD} = [];
for i = 1:nD
    Dpar{i} = par;
    Dpar{i}(1) = Dv(i);
end

Dpvasc = zeros(np,np,nD); %preallocate matricies to store fixed points
Dpalv  = zeros(np,np,nD);

%%% make all of the things for diffusion
parfor k = 1:length(Dv)
    for i = 1:np
        for j = 1:np
            [Dpvasc(i,j,k),Dpalv(i,j,k)] = modelD_SS_relaxation(N,Dpar{k},P,C,0+i*de,0+j*de); %loop through q and v and store fixed points
            disp([k i j])
        end
    end
end
Ddpac = Dpalv - Dpvasc;

% save('ModelD_gif_results.mat')
load('ModelD_gif_results.mat')

%%% making gifs
CON = -30:5:150;
x = 0:10;
h = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(Dv) %iterate by size of Dv
    %%%% make contour plots
    cla(gca)
    
    contour(Q,V,squeeze(Dpvasc(:,:,i)),CON,'ShowText','on')
    hold on
    plot(x,x,'k--')
    set(gca,'fontsize',18)
    xlabel('Cardiac Output (ml/s)')
    ylabel('Ventilation Magnitude (ml/s)')
    text(9,1,['D = ', num2str(Dv(i))],'fontsize',18)
    axis equal
    grid on
    
    if i == 1
        gif('ModelD_vasc.gif','DelayTime',1/15) %make gif file
    else
        gif %append frame to gif
    end
    disp(i)
end

h2 = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(Dv) %iterate by size of Dv
    %%%% make contour plots
    cla(gca)
    
    contour(Q,V,squeeze(Dpalv(:,:,i)),CON,'ShowText','on')
    hold on
    plot(x,x,'k--')
    set(gca,'fontsize',18)
    xlabel('Cardiac Output (ml/s)')
    ylabel('Ventilation Magnitude (ml/s)')
    text(9,1,['D = ', num2str(Dv(i))],'fontsize',18)
    axis equal
    grid on
    
    if i == 1
        gif('ModelD_alv.gif','DelayTime',1/15) %make gif file
    else
        gif %append frame to gif
    end
    disp(i)
end

h3 = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(Dv) %iterate by size of Dv
    %%%% make contour plots
    cla(gca)
    
    contour(Q,V,squeeze(Ddpac(:,:,i)),CON,'ShowText','on')
    hold on
    plot(x,x,'k--')
    set(gca,'fontsize',18)
    xlabel('Cardiac Output (ml/s)')
    ylabel('Ventilation Magnitude (ml/s)')
    text(9,1,['D = ', num2str(Dv(i))],'fontsize',18)
    axis equal
    grid on
    
    if i == 1
        gif('ModelD_grad.gif','DelayTime',1/15) %make gif file
    else
        gif %append frame to gif
    end
    disp(i)
end
