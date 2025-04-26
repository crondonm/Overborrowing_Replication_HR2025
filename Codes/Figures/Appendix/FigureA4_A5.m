%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% In this code: Plot Figure A4 to A5
% 
% Authors:  Juan Herreño, jherrenolopera@ucsd.edu
%               Carlos Rondón Moreno, crondon@bcentral.cl
%
% Date: March 2025
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping

clearvars
clear global
close all

% Load databases


load('../../Replication/Data/FICE.mat')
load('../../Replication/Data/FIP.mat')
load('../../Replication/Data/IIPCC.mat')
load('../../Replication/Data/IICECC.mat')

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

load('../../Replication/Data/Param.mat')
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

%% Figure A4: Appendix - Policy Functions Imperfect Information vs Full Information

close all

Format.styles = {'-',':','-.','--'};

% Pick the realization that is closest to a 1 standard deviation on the
% permanent component

Temp1 = [10	10 6 6 6]; 
Temp4 = [10 10 6];

% Approximate the grid to our simulation:

Temp2 = [findClosest2(FICE.S(:,1),FICE.yt) findClosest2(FICE.S(:,2),FICE.yn) findClosest2(FICE.S(:,3),FICE.gt+g)] ;
Temp3 = [findClosest2(IIPCC.S(:,1),IIPCC.yt) findClosest2(IIPCC.S(:,2),IIPCC.yn) findClosest2(IIPCC.S(:,3),IIPCC.gt+g) findClosest2(IIPCC.S(:,4),IIPCC.gT+g) findClosest2(IIPCC.S(:,5),IIPCC.gN+g)] ;

% Pick the events where our simulation observes the desired states:

index = find((Temp3(:,1)==Temp1(:,1)).*(Temp3(:,2)==Temp1(:,2)).*(Temp3(:,3)==Temp1(:,3).*(Temp3(:,4)==Temp1(:,4)).*(Temp3(:,5)==Temp1(:,5))));
index2 = find((Temp2(:,1)==Temp4(:,1)).*(Temp2(:,2)==Temp4(:,2)).*(Temp2(:,3)==Temp4(:,3)));

f10 = figure('Position',Format.figsize,'Color',[1 1 1]);
p1 = plot(IIPCC.b, IIPCC.b(IIPCC.Pol(index,:)), IICECC.b, IICECC.b(IICECC.Pol(index,:)), FIP.b, FIP.b(FIP.Pol(index2,:)), FICE.b, FICE.b(FICE.Pol(index2,:)), FICE.b, FICE.b);
for linei = 1:5
    if linei < 5
        set(p1(linei),'LineStyle',Format.styles{linei})
        set(p1(linei),'color',Format.colors{linei})
        set(p1(linei),'linewidth',Format.widths{linei})
    else
        set(p1(linei),'color',[0,0,0])
    end
end
ylabel('Next Period Bond Holdings (Level)','FontSize',Format.FontSize, 'Interpreter','Latex')
xlabel('Current Bond Holdings (Level)','FontSize',Format.FontSize, 'Interpreter','Latex')
legend('Social Planner: Imperfect Information','Decentralized Equilibrium: Imperfect Information', 'Social Planner: Full Information', 'Decentralized Equilibrium: Full Information','Location','SouthEast')

xlim([-1.16 -0.2 ]) ;
ylim([-1.1 -0.2 ]) ;
set(gca,'FontSize',Format.FontSizeAxes)
set(f10,'PaperPositionMode', 'auto')
%saveas(f10,'FigureA4','png');


filename = 'FigureA4.png';
resolution = 300; % DPI

% Export the graphics
exportgraphics(gcf, filename, 'Resolution', resolution);

%% Figure A5a: Appendix - Policy Functions Imperfect Information vs Recalibrated

close all

Format.styles = {'-',':','-.','--'};

FICErec = load('../../Replication/Data/FICE_rec.mat');
FICErec = FICErec.FICE;
FIPrec = load('../../Replication/Data/FIP_rec.mat');
FIPrec = FIPrec.FIP;

% Pick the realization that is closest to a 1 standard deviation on the
% permanent component

Temp1 = [10	10 6 6 6]; 
Temp4 = [10 10 6];

% Approximate the grid to our simulation:

Temp2 = [findClosest2(FICErec.S(:,1),FICErec.yt) findClosest2(FICErec.S(:,2),FICErec.yn) findClosest2(FICErec.S(:,3),FICErec.gt+g)] ;
Temp3 = [findClosest2(IIPCC.S(:,1),IIPCC.yt) findClosest2(IIPCC.S(:,2),IIPCC.yn) findClosest2(IIPCC.S(:,3),IIPCC.gt+g) findClosest2(IIPCC.S(:,4),IIPCC.gT+g) findClosest2(IIPCC.S(:,5),IIPCC.gN+g)] ;

% Pick the events where our simulation observes the desired states:

index = find((Temp3(:,1)==Temp1(:,1)).*(Temp3(:,2)==Temp1(:,2)).*(Temp3(:,3)==Temp1(:,3).*(Temp3(:,4)==Temp1(:,4)).*(Temp3(:,5)==Temp1(:,5))));
index2 = find((Temp2(:,1)==Temp4(:,1)).*(Temp2(:,2)==Temp4(:,2)).*(Temp2(:,3)==Temp4(:,3)));

f11 = figure('Position',Format.figsize,'Color',[1 1 1]);
p1 = plot(IIPCC.b, IIPCC.b(IIPCC.Pol(index,:)), IICECC.b, IICECC.b(IICECC.Pol(index,:)), FIPrec.b, FIPrec.b(FIPrec.Pol(index2,:)), FICErec.b, FICErec.b(FICErec.Pol(index2,:)), FICErec.b, FICErec.b);
for linei = 1:5
    if linei < 5
        set(p1(linei),'LineStyle',Format.styles{linei})
        set(p1(linei),'color',Format.colors{linei})
        set(p1(linei),'linewidth',Format.widths{linei})
    else
        set(p1(linei),'color',[0,0,0])
    end
end
%title('Imperfect Information','FontSize',Format.FontSize,'FontWeight',Format.fontweight)
ylabel('Next Period Bond Holdings (Level)','FontSize',Format.FontSize, 'Interpreter','Latex')
xlabel('Current Bond Holdings (Level)','FontSize',Format.FontSize, 'Interpreter','Latex')
legend('Social Planner: Baseline','Decentralized Equilibrium: Baseline', 'Social Planner: Recalibrated', 'Decentralized Equilibrium: Recalibrated','Location','SouthEast')

xlim([-1.16 -0.2 ]) ;
ylim([-1.1 -0.2 ]) ;
set(gca,'FontSize',Format.FontSizeAxes)
set(f11,'PaperPositionMode', 'auto')
%saveas(f11,'FigureA5a','png');

filename = 'FigureA5a.png';
resolution = 300; % DPI

% Export the graphics
exportgraphics(gcf, filename, 'Resolution', resolution);

%% Figure A5b: Appendix - Policy Functions Recalibrated Economy - Full Information

close all

% Pick the realization that is closest to a 1 negative standard deviation on the permanent component

Temp4 = [10 10 6];

% Approximate the grid to our simulation:

Temp2 = [findClosest2(FICErec.S(:,1),FICErec.yt) findClosest2(FICErec.S(:,2),FICErec.yn) findClosest2(FICErec.S(:,3),FICErec.gt+g)] ;

% Pick the events where our simulation observes the desired states:

index = find((Temp2(:,1)==Temp4(:,1)).*(Temp2(:,2)==Temp4(:,2)).*(Temp2(:,3)==Temp4(:,3)));

f12 = figure('Position',Format.figsize,'Color',[1 1 1]);
p1 = plot(FIP.b, FIP.b(FIP.Pol(index,:)), FICE.b, FICE.b(FICE.Pol(index,:)), FIPrec.b, FIPrec.b(FIPrec.Pol(index,:)), FICErec.b, FICErec.b(FICErec.Pol(index,:)), FICE.b, FICE.b);
for linei = 1:5
    if linei < 5
        set(p1(linei),'LineStyle',Format.styles{linei})
        set(p1(linei),'color',Format.colors{linei})
        set(p1(linei),'linewidth',Format.widths{linei})
    else
        set(p1(linei),'color',[0,0,0])
    end
end
%title('Full Information','FontSize',Format.FontSize,'FontWeight',Format.fontweight)
ylabel('Next Period Bond Holdings (Level)','FontSize',Format.FontSize, 'Interpreter','Latex')
xlabel('Current Bond Holdings (Level)','FontSize',Format.FontSize, 'Interpreter','Latex')
legend('Social Planner: Baseline','Decentralized Equilibrium: Baseline', 'Social Planner: Recalibrated', 'Decentralized Equilibrium: Recalibrated','Location','SouthEast')
xlim([-1.16 -0.2 ]) ;
ylim([-1.16 -0.2 ]) ;
set(gca,'FontSize',Format.FontSizeAxes)
set(f12,'PaperPositionMode', 'auto')
%saveas(f12,'FigureA5b','png');

filename = 'FigureA5b.png';
resolution = 300; % DPI

% Export the graphics
exportgraphics(gcf, filename, 'Resolution', resolution);