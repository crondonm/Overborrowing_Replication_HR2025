%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% In this code: We plot all the figures included in the paper                   
% 
% Authors: Juan Herreño, jherrenolopera@ucsd.edu 
%              Carlos Rondón Moreno, crondon@bcentral.cl
%
% Last Update:  March 2025
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping

clearvars
clear global
close all

% Load databases


load('../Replication/Data/FICE.mat')
load('../Replication/Data/FIP.mat')
load('../Replication/Data/IIPCC.mat')
load('../Replication/Data/IICECC.mat')
FIPrec = load('../Replication/Data/FIP_rec.mat');
FIPrec = FIPrec.FIP;

% Figure parameters

Format.FontSize = 16;
Format.colors = {'k','r','k','r'};
Format.colors = {[0/255 0/255 102/255], [65/255 105/255 225/255],[123/255 104/255 238/255],[186/255 85/255 211/255]};
Format.widths = {3,3,3,3};
Format.styles = {'-',':','-.','-.'};
Format.figsize  = [460 200 700 525];
Format.figsize2 = [460 200 600 250];
Format.figsize3 = [460 200 700 250];
Format.figsize4 = [460 200 600 450];
Format.figsize6 = [460 200 875 875];
Format.FontSizeAxes = 12;
Format.fontweight = 'bold';
Format.XTick = [1975 1980 1985 1990 1995 2000 2005 2010 2015] ;
Format.Xlim =[1975 2015] ;
Format.Shade = 0.8 ;
Format.ScatterSize = [50,50] ;
Format.ScatterMarker = {'o','s'} ;

% Load Parameters

load('../Replication/Data/Param.mat')
fprintf("Parameters loaded... \n")

% Assign parameters
r = Param.r;
g =  Param.g;  % Mean growth rate of permanent component
burn = Param.burn;
nstd = Param.nstd;
init_debt = Param.init_debt;
horizn = Param.horizn;
ss = 21; % Position at which the IRF returns to origin
window = Param.window;


%% Figure 7: Optimal Macroprudential Policy: Optimal Tax Functions
% Optimal Tax policy functions
% We need to plot the functions for a particular realization of the states.
% To do this, first we pick which realization we are going to draw:

close all

% Standard Deviations
ZBBT = -nstd*std(IIPCC.S(:,1));
ZBBN = -nstd*std(IIPCC.S(:,2));
GBBT = -nstd*std(IIPCC.S(:,3));
GTBBT = -nstd*std(IIPCC.S(:,4));
GNBBT = -nstd*std(IIPCC.S(:,5));

% IIP Policy Function:
% Pick the realization that is closest to a 1 standard deviation on all the shocks

%Temp1 = [findClosest2(ZBBT,IIPCC.yt) findClosest2(ZBBN,IIPCC.yn) findClosest2(GBBT,IIPCC.gt+g) findClosest2(GTBBT,IIPCC.gT+g) findClosest2(GTBBT,IIPCC.gN+g)] ;
%Temp4 = [findClosest2(ZBBT,FIP.yt) findClosest2(ZBBN,FIP.yn) findClosest2(GBBT,FIP.gt+g)];

Temp1 = [10	10 6 6 6]; 
Temp4 = [10 10 6];

% Approximate the grid to our simulation:

Temp2 = [findClosest2(FIP.S(:,1),FIP.yt) findClosest2(FIP.S(:,2),FIP.yn) findClosest2(FIP.S(:,3),FIP.gt+g)] ;
Temp3 = [findClosest2(IIPCC.S(:,1),IIPCC.yt) findClosest2(IIPCC.S(:,2),IIPCC.yn) findClosest2(IIPCC.S(:,3),IIPCC.gt+g) findClosest2(IIPCC.S(:,4),IIPCC.gT+g) findClosest2(IIPCC.S(:,5),IIPCC.gN+g)] ;

% Compute average total income:
IIPmY = mean(IIPCC.b(IIPCC.SimB(2:end)).*exp(IIPCC.Sim(burn+2:end-1,1)')./IIPCC.DtoY(2:end));
FIPmY = mean(FIP.SimBhat(2:end).*exp(FIP.Sim(burn+1:end-1,1)')./FIP.DtoY(2:end));

% Pick the events where our simulation observes the desired states:

index = find((Temp3(:,1)==Temp1(:,1)).*(Temp3(:,2)==Temp1(:,2)).*(Temp3(:,3)==Temp1(:,3).*(Temp3(:,4)==Temp1(:,4)).*(Temp3(:,5)==Temp1(:,5))));
index2 = find((Temp2(:,1)==Temp4(:,1)).*(Temp2(:,2)==Temp4(:,2)).*(Temp2(:,3)==Temp4(:,3)));

% Extract the optimal tax corresponding to the realization
IIPCC.StdPol = IIPCC.TAO(index,:) ;
FIP.StdPol   = FIP.TAO(index2,:) ;

% Here we create the shaded areas for our figure:
% First we need to locate the values for which the tax is zero:
AA = IIPCC.StdPol;
AA(1:176) = NaN;
BB  = IIPCC.StdPol;
BB(1:176) = 0;
BB(177:end) = NaN;

% We repeat the same process for the full information model:
Zzero = find(FIP.StdPol==0);
CC = FIP.StdPol;
CC(1:278) = NaN;
DD  = FIP.StdPol;
DD(1:278) = 0;
DD(279:end) = NaN;
Format.styles = {'-','-','-.','+'};

f7 = figure('Position',Format.figsize2,'Color',[1 1 1]);
subplot(1,2,1)
p1 = plot(IIPCC.b/IIPmY*100, AA*100,IIPCC.b/IIPmY*100, BB);
for linei = 1:2
    set(p1(linei),'LineStyle',Format.styles{linei})
    set(p1(linei),'color',Format.colors{linei})
    set(p1(linei),'linewidth',Format.widths{linei})
end
title('Imperfect Information','FontSize',Format.FontSize,'FontWeight',Format.fontweight, 'Interpreter','Latex')
ylabel('Optimal Tax (\%)','FontSize',Format.FontSize, 'Interpreter','Latex')
xlabel('Bond Holdings (Percent of GDP)','FontSize',Format.FontSize, 'Interpreter','Latex')
xlim([-0.40*100 -0.1*100 ]) ;
ylim([0 50]) ;
set(gca,'FontSize',Format.FontSizeAxes)
set(f7,'PaperPositionMode', 'auto')
annotation(f7,'rectangle',...
    [0.24 0.1362 0.09 0.789],...
    'LineStyle','none',...
    'FaceColor',[0.800000011920929 0.800000011920929 0.800000011920929],...
    'FaceAlpha',0.4);
subplot(1,2,2)
p1 = plot(FIP.b/FIPmY*100,CC*100, FIP.b/FIPmY*100, DD);
for linei = 1:2
    set(p1(linei),'LineStyle',Format.styles{linei})
    set(p1(linei),'color',Format.colors{linei})
    set(p1(linei),'linewidth',Format.widths{linei})
end
title('Full Information','FontSize',Format.FontSize,'FontWeight',Format.fontweight, 'Interpreter','Latex')
ylabel('Optimal Tax (\%)','FontSize',Format.FontSize, 'Interpreter','Latex')
xlabel('Bond Holdings (Percent of GDP)','FontSize',Format.FontSize, 'Interpreter','Latex')
xlim([-0.3*100 -0.1*100 ]) ;
ylim([0 30]);
v=axis;
set(gca,'FontSize',Format.FontSizeAxes)
set(f7,'PaperPositionMode', 'auto')
annotation(f7,'rectangle',...
     [0.695 0.1357 0.1 0.789],...
     'LineStyle','none',...
     'FaceColor',[0.800000011920929 0.800000011920929 0.800000011920929],...
     'FaceAlpha',0.4);
% add a bit space to the figure
fig = gcf;
fig.Position(3) = fig.Position(3) + 125;
fig.Position(4) = fig.Position(4) + 30;
%saveas(f7,'Figure7_TaxFun_IIPvsFIP','png');
filename = 'Figure7.png';
resolution = 300; % DPI
% Export the graphics
exportgraphics(gcf, filename, 'Resolution', resolution);

Format.styles = {'-',':','-.','-.'};


%% Figure 8: Optimal Macroprudential Policy: Optimal Tax Functions
% Ergodic Distribution of Debt

[datay1, datax1] = ksdensity(FIP.TAOSim(FIP.TAOSim<10)*100, 'Support','nonnegative', 'BoundaryCorrection', 'reflection');
[datay2, datax2] = ksdensity(IIPCC.TAOSim   *100, 'Support','nonnegative', 'BoundaryCorrection', 'reflection');    

Format.styles = {'-.', '-', '-.','-.'};

f9 = figure('Color',[1 1 1], 'Position',Format.figsize);
p1 = plot(datax1, datay1, datax2, datay2) ;
for linei = 1:2
    set(p1(linei),'LineStyle',Format.styles{linei})
    set(p1(linei),'color',Format.colors{linei})
    set(p1(linei),'linewidth',Format.widths{linei})
end
%title('Ergodic Distribution of Asset Holdings: Incomplete Information','FontSize',Format.FontSize,'FontWeight',Format.fontweight)
ylabel('Frequency','FontSize',Format.FontSize, 'Interpreter','Latex')
xlabel('Optimal Tax Rate (\%) ','FontSize',Format.FontSize, 'Interpreter','Latex')
legend('Perfect Information','Imperfect Information','Location','NorthEast', 'Fontsize', 16)
xlim([0 50 ]) ;
%ylim([-0.5 .07 ]) ;
set(gca,'FontSize',Format.FontSizeAxes)
set(f9,'PaperPositionMode', 'auto')
%saveas(f9,'Figure8_ErgDist_TaxesIIvsFI','png');

filename = 'Figure8.png';
resolution = 300; % DPI
% Export the graphics
exportgraphics(gcf, filename, 'Resolution', resolution);
