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

load('../Replication/Data/Param.mat')
fprintf("Parameters loaded... \n")

load('../Replication/Data/IIPCCsim.mat')
fprintf("Database loaded... \n")

g =  Param.g;  % Mean growth rate of permanent component

Tsim = Param.Tsim;   % Simulation points
burn = Param.burn; % Burn-in period for simulation
nstd = Param.nstd;
window = Param.window;

% Settings for figures

Format.FontSize = 18;
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

%% Simulation

fprintf("Starting crises analysis: Imperfect Information ... \n")

AAA = b(SimB) > BCSim + (b(2) - b(1))/2 ; 
CCC =CA;
CCCT = nstd*std(CCC);
Crisis = (CCC > CCCT).*(1 - AAA) ;
Freq = sum(Crisis)/(length(CCC)) ;
CrInd = find(Crisis == 1) ;
CrInd = CrInd(CrInd > window + 1) ; 
CrInd = CrInd(CrInd < Tsim - burn - window) ;  

% Another definition Boom-Bust 

clear BBII 

BBII.mP = mP;
BBII.mCT = mCT;

% BB in the permanent component

simg = Sim(burn+2:end,3);
stdGt = nstd*std(simg);

CrIndG = find(simg<=-stdGt);
CrIndG = CrIndG(CrIndG > window + 1); 
CrIndG = CrIndG(CrIndG < Tsim - burn - window);

CrIndG_pos = find(simg>=stdGt);
CrIndG_pos = CrIndG_pos(CrIndG_pos > window + 1);
CrIndG_pos = CrIndG_pos(CrIndG_pos < Tsim - burn - window);

% BB in the transitory component tradable

simz = Sim(burn+2:end,1);
stdZt = nstd*std(simz);

CrIndZ = find(simz<-stdZt);
CrIndZ = CrIndZ(CrIndZ > window + 1);
CrIndZ = CrIndZ(CrIndZ < Tsim - burn - window);

CrIndZ_pos = find(simz>stdZt);
CrIndZ_pos = CrIndZ_pos(CrIndZ_pos > window + 1);
CrIndZ_pos = CrIndZ_pos(CrIndZ_pos < Tsim - burn - window);

% BB in the transitory component nontradable

simzn = Sim(burn+2:end,2);
stdZn = nstd*std(simzn);

CrIndZn = find(simzn<-stdZn);
CrIndZn = CrIndZn(CrIndZn > window + 1);
CrIndZn = CrIndZn(CrIndZn < Tsim - burn - window);

CrIndZn_pos = find(simzn>stdZn);
CrIndZn_pos = CrIndZn_pos(CrIndZn_pos > window + 1);
CrIndZn_pos = CrIndZn_pos(CrIndZn_pos < Tsim - burn - window);

PSim = PSim + mP;
Ytot = ((exp(simz(2:end)'+simg(2:end)'+g) + PSim(2:end).*exp(simzn(2:end)'+simg(2:end)'+g)))./exp(simz(1:end- 1)');


for i=-window:window
    BBII.TAOSim(i + window + 1,:) = TAOSim(CrInd + i);
    BBII.Ytot(i + window + 1,:) = Ytot(CrInd + i - 1);
    BBII.TAOSimG(i + window + 1,:) = TAOSim(CrIndG + i);
    BBII.TAOSimZ(i + window + 1,:) = TAOSim(CrIndZ + i);
    BBII.TAOSimZn(i + window + 1,:) = TAOSim(CrIndZn + i);
    BBII.GT(i + window + 1,:) = gt(Posteriorg(CrIndG + i));
    BBII.ZT(i + window + 1,:) = yt(Posterioryt(CrIndZ + i));
    BBII.ZN(i + window + 1,:) = yn(Posterioryn(CrIndZn + i));

    BBII.TAOSimG_pos(i + window + 1,:) = TAOSim(CrIndG_pos + i);
    BBII.TAOSimZ_pos(i + window + 1,:) = TAOSim(CrIndZ_pos + i);
    BBII.TAOSimZn_pos(i + window + 1,:) = TAOSim(CrIndZn_pos + i);
end 

BBII.IRTAOSim = mean(BBII.TAOSim, 2);
BBII.IRTAOSimG = mean(BBII.TAOSimG, 2);
BBII.IRTAOSimZ = mean(BBII.TAOSimZ, 2);
BBII.IRTAOSimZn = mean(BBII.TAOSimZn, 2);
BBII.IRTAOSimG_pos = mean(BBII.TAOSimG_pos, 2);
BBII.IRTAOSimZ_pos = mean(BBII.TAOSimZ_pos, 2);
BBII.IRTAOSimZn_pos = mean(BBII.TAOSimZn_pos, 2);
BBII.GT = mean(BBII.GT, 2);
BBII.ZT = mean(BBII.ZT, 2);
BBII.ZN = mean(BBII.ZN, 2);
BBII.IRYtot = mean(BBII.Ytot, 2);


%% Full Information

clearvars -except BBII g Tsim burn nstd window Format 

load('../Replication/Data/FIPsim.mat')


fprintf("Starting crises analysis: Full Information \n")

AAA = b(SimB) > (BCSim2 + (b(2)-b(1))/2) ; 
CCC = CA ; 
CCCT = nstd*std(CCC);
Crisis = (CCC>CCCT).*(1-AAA) ;
Crisis2 = sum(Crisis)/(length(CCC)) ;

CrInd = find(Crisis==1) ;
CrInd = CrInd(CrInd>window + 1) ; 
CrInd = CrInd(CrInd< Tsim - burn - window) ;   

% BB in the permanent component
Simg = Sim(burn + 1:end, 3);
stdGt = nstd*std(Simg);

CrIndG = find(Simg<-stdGt);
CrIndG = CrIndG(CrIndG > window + 1); 
CrIndG = CrIndG(CrIndG < Tsim - burn - window);

CrIndG_pos = find(Simg>stdGt);
CrIndG_pos = CrIndG_pos(CrIndG_pos > window + 1);
CrIndG_pos = CrIndG_pos(CrIndG_pos < Tsim - burn - window);

% BB in the transitory component Tradable
Simz = Sim(burn + 1:end, 1);
stdZt = nstd*std(Simz);

CrIndZ = find(Simz<-stdZt);
CrIndZ = CrIndZ(CrIndZ > window + 1);
CrIndZ = CrIndZ(CrIndZ < Tsim - burn - window);

CrIndZ_pos = find(Simz>stdZt);
CrIndZ_pos = CrIndZ_pos(CrIndZ_pos > window + 1);
CrIndZ_pos = CrIndZ_pos(CrIndZ_pos < Tsim - burn - window);

% BB in the transitory component Non-Tradable

Simzn = Sim(burn + 1:end, 2);
stdZn = nstd*std(Simzn);

CrIndZn = find(Simzn<-stdZn);
CrIndZn = CrIndZn(CrIndZn > window + 1);
CrIndZn = CrIndZn(CrIndZn < Tsim - burn - window);

CrIndZn_pos = find(Simzn>stdZn);
CrIndZn_pos = CrIndZn_pos(CrIndZn_pos > window + 1);
CrIndZn_pos = CrIndZn_pos(CrIndZn_pos < Tsim - burn - window);

PSim = PSim + mP;
Ytot = ((exp(Simz(2:end)'+Simg(2:end)'+g) + PSim(2:end).*exp(Simzn(2:end)'+Simg(2:end)'+g)))./exp(Simz(1:end- 1)');

for i=-window:window
    BBFI.Ytot(i + window + 1,:) = Ytot(CrInd + i - 1);
    BBFI.TAOSim(i + window + 1,:) = TAOSim(CrInd + i);
    BBFI.TAOSimG(i + window + 1,:) = TAOSim(CrIndG + i );
    BBFI.TAOSimZ(i + window + 1,:) = TAOSim(CrIndZ + i);
    BBFI.TAOSimZn(i + window + 1,:) = TAOSim(CrIndZn + i);
    BBFI.TAOSimG_pos(i + window + 1,:) = TAOSim(CrIndG_pos + i);
    BBFI.TAOSimZ_pos(i + window + 1,:) = TAOSim(CrIndZ_pos + i);
    BBFI.TAOSimZn_pos(i + window + 1,:) = TAOSim(CrIndZn_pos + i);
    BBFI.GT(i + window + 1,:) = Simg(CrIndG + i);
    BBFI.GT_pos(i + window + 1,:) = Simg(CrIndG_pos + i);
    BBFI.IRB(i + window + 1,:) = b(SimB(CrIndG_pos + i));
    BBFI.IRBC(i + window + 1,:) = BCSim2(CrIndG_pos + i) ;
end 

BBFI.IRYtot = mean(BBFI.Ytot, 2);
BBFI.IRTAOSim = mean(BBFI.TAOSim, 2, "omitmissing");
BBFI.IRTAOSimG = mean(BBFI.TAOSimG, 2, "omitmissing");
BBFI.IRTAOSimZ = mean(BBFI.TAOSimZ, 2);
BBFI.IRTAOSimZn = mean(BBFI.TAOSimZn, 2);
BBFI.IRTAOSimG_pos = mean(BBFI.TAOSimG_pos, 2);
BBFI.IRTAOSimZ_pos = mean(BBFI.TAOSimZ_pos, 2);
BBFI.IRTAOSimZn_pos = mean(BBFI.TAOSimZn_pos, 2);
BBFI.IRGT = mean(BBFI.GT, 2);
BBFI.IRGT_pos = mean(BBFI.GT_pos, 2);

%% Figure 12

close all

ff= figure('Position',Format.figsize2,'Color',[1 1 1]);
tt = tiledlayout(1,2);
ax=nexttile;
p1 = plot(-window:window, BBFI.IRTAOSim*100);
linei = 1;
set(p1,'LineStyle',Format.styles{linei})
set(p1,'color',Format.colors{linei})
set(p1,'linewidth',Format.widths{linei})
set(gca,'ycolor',Format.colors{linei})
ylim([0,18])
xticks([-5,-4,-3,-2,-1,0,1,2,3,4,5 ]) ;
xlabel('Years','FontSize',14, 'Interpreter','Latex')
ylabel('Tax Rate (\%)','FontSize',14, 'Interpreter','Latex')
title('Perfect Information', 'FontSize',16, 'Interpreter','Latex')
yyaxis right
p1 = plot(-window:window, BBFI.IRYtot);
linei = 2;
set(p1,'LineStyle',Format.styles{linei})
set(p1,'color',Format.colors{linei})
set(p1,'linewidth',Format.widths{linei})
set(gca,'ycolor',Format.colors{linei})
ylabel('Total Income (Level) ','FontSize',14, 'Interpreter','Latex')
ylim([4,6])
ax=nexttile;
p1 = plot(-window:window, BBII.IRTAOSim * 100);
linei = 1;
set(p1,'LineStyle',Format.styles{linei})
set(p1,'color',Format.colors{linei})
set(p1,'linewidth',Format.widths{linei})
set(gca,'ycolor',Format.colors{linei})
ylim([0,18])
xticks([-5,-4,-3,-2,-1,0,1,2,3,4,5 ]) ;
xlabel('Years','FontSize',14, 'Interpreter','Latex')
ylabel('Tax Rate (\%)','FontSize',14, 'Interpreter','Latex')
title('Imperfect Information', 'FontSize',16, 'Interpreter','Latex')

yyaxis right
p2 = plot(-window:window, BBII.IRYtot);
linei = 2;
set(p2,'LineStyle',Format.styles{linei})
set(p2,'color',Format.colors{linei})
set(p2,'linewidth',Format.widths{linei})
set(gca,'ycolor',Format.colors{linei})
ylim([4,6])
xticks([-5,-4,-3,-2,-1,0,1,2,3,4,5 ]) ;
ylabel('Total Income (Level)','FontSize',14, 'Interpreter','Latex')
lg  = legend(ax,'Optimal Tax (\%, Left-hand Axis)', 'Total Income (Right-hand Axis)', ...
             'Orientation','Horizontal','NumColumns',2,'FontSize', Format.FontSize,'Interpreter','Latex');
lg.Layout.Tile = 'South';
% add a bit space to the figure
fig = gcf;
fig.Position(3) = fig.Position(3) + 375;
fig.Position(4) = fig.Position(4) + 75;

% Save figure

saveas(ff,'Figure12','png');

corr(BBII.IRYtot, BBII.IRTAOSim)
corr(BBFI.IRYtot, BBFI.IRTAOSim)