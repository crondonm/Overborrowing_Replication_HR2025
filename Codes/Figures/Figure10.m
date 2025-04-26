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
load("../Replication/Data/IIPCCsim.mat") ; % load imperfect info CE solution
fprintf("Database loaded... \n")

g =  Param.g;  % Mean growth rate of permanent component
Tsim = Param.Tsim;   % Simulation points
burn = Param.burn; % Burn-in period for simulation
nstd = Param.nstd;
window = 2;

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

%% Imperfect Information

fprintf("Starting Boom-Bust: Imperfect Information ... \n")

% Total Output

Simyt = Sim(burn+2:end, 1);
Simyn = Sim(burn+2:end, 2);
Simgt = Sim(burn+2:end, 3);
Ytot = ((exp(Simyt(2:end)'+Simgt(2:end)' + g) + PSim(2:end).*exp(Simyn(2:end)'+Simgt(2:end)' + g)))./exp(Simyt(1:end- 1)');
Ytot_fwd = Ytot(window+1:end);
Ytot_bwd = Ytot(1:end-window);

% Income Components

Simyt = yt(Posterioryt);
Simyt_fwd = Simyt(window+1:end);
Simyt_bwd = Simyt(1:end-window);

Simgt = gt(Posteriorg);
Simgt_fwd = Simgt(window+1:end);
Simgt_bwd = Simgt(1:end-window);

Simyn = yn(Posterioryn);
Simyn_fwd = Simyn(window+1:end);
Simyn_bwd = Simyn(1:end-window);


Crisis = (Simyt_fwd < -nstd*std(Simyt)).*(Simyt_bwd > nstd*std(Simyt)) ;
CrInd = find(Crisis == 1) ;
CrInd = CrInd(CrInd > window + 1) ; 
CrInd = CrInd(CrInd < Tsim - burn - window - 5);  

Crisis2 = (Simgt_fwd < -nstd*std(Simgt)).*(Simgt_bwd > nstd*std(Simgt)) ;
CrIndg = find(Crisis2 == 1) ;
CrIndg = CrIndg(CrIndg > window + 1) ; 
CrIndg = CrIndg(CrIndg < Tsim - burn - window - 5) ;  

Crisis3 = (Simyn_fwd < -nstd*std(Simyn)).*(Simyn_bwd > nstd*std(Simyn)) ;
CrIndn = find(Crisis3 == 1) ;
CrIndn = CrIndn(CrIndn > window + 1) ;
CrIndn = CrIndn(CrIndn < Tsim - burn - window - 5) ;

bbsInd = find(Ytot_fwd<(mean(Ytot)-std(Ytot)) & Ytot_bwd>(mean(Ytot) + std(Ytot)));
bbsInd = bbsInd(bbsInd > window + 1) ; 
bbsInd = bbsInd(bbsInd < Tsim - burn - window - 5);  

% Boom-Bust 

clear BBII
for i=-window:window+2
    BBII.TAOSimZt(i + window + 1,:) = TAOSim(CrInd + i);
    BBII.SimZt(i + window + 1,:) = Simyt(CrInd + i);
    BBII.TAOSimGt(i + window + 1,:) = TAOSim(CrIndg + i);
    BBII.SimGt(i + window + 1,:) = Simgt(CrIndg + i);
    BBII.TAOSimNt(i + window + 1,:) = TAOSim(CrIndn + i);
    BBII.SimNt(i + window + 1,:) = Simyn(CrIndn + i);
    BBII.TAOSimYt(i + window + 1,:) = TAOSim(bbsInd + i + 1);
    BBII.SimYt(i + window + 1,:) = Ytot(bbsInd + i);
end

BBII.IRTAOSimZt = mean(BBII.TAOSimZt, 2);
BBII.IRSimZt = mean(BBII.SimZt, 2);
BBII.IRTAOSimGt = mean(BBII.TAOSimGt, 2);
BBII.IRSimGt = mean(BBII.SimGt, 2);
BBII.IRTAOSimNt = mean(BBII.TAOSimNt, 2);
BBII.IRSimNt = mean(BBII.SimNt, 2);
BBII.IRTAOYt = mean(BBII.TAOSimYt, 2);
BBII.IRYt = mean(BBII.SimYt, 2);




%% Full Information

clearvars -except BBII g Tsim burn nstd window Format 

load('../Replication/Data/FIPsim.mat')

idx = find(TAOSim>10);
TAOSim(idx) = nan;

% Total Output

Simyt = Sim(burn+1:end, 1);
Simyn = Sim(burn+1:end, 2);
Simgt = Sim(burn+1:end, 3);
Ytot = ((exp(Simyt(2:end)'+Simgt(2:end)' + g) +(PSim(2:end)+mP).*exp(Simyn(2:end)'+Simgt(2:end)' + g)))./exp(Simyt(1:end- 1)');
Ytot_fwd = Ytot(window+1:end);
Ytot_bwd = Ytot(1:end-window);

bbsInd = find(Ytot_fwd<(mean(Ytot)-std(Ytot)) & Ytot_bwd>(mean(Ytot) + std(Ytot)));
bbsInd = bbsInd(bbsInd > window + 1) ; 
bbsInd = bbsInd(bbsInd < Tsim - burn - window - 5);  


% Income components:

Posterioryt = findClosest2(Sim(burn+1:end,1), yt);
Posterioryn = findClosest2(Sim(burn+1:end,2), yn);
Posteriorg = findClosest2(Sim(burn+1:end,3), gt);

Simyt = yt(Posterioryt);
Simyt_fwd = Simyt(window+1:end);
Simyt_bwd = Simyt(1:end-window);

Simgt =  gt(Posteriorg);
Simgt_fwd = Simgt(window+1:end);
Simgt_bwd = Simgt(1:end-window);

Simyn = yn(Posterioryn);
Simyn_fwd = Simyn(window+1:end);
Simyn_bwd = Simyn(1:end-window);

Crisis = (Simyt_fwd < -nstd*std(Simyt)).*(Simyt_bwd > nstd*std(Simyt)) ;

CrInd = find(Crisis == 1) ;
CrInd = CrInd(CrInd > window + 1) ;
CrInd = CrInd(CrInd < Tsim - burn - window - 5) ;

Crisis2 = (Simgt_fwd < -nstd*std(Simgt)).*(Simgt_bwd > nstd*std(Simgt)) ;
CrIndg = find(Crisis2 == 1) ;
CrIndg = CrIndg(CrIndg > window + 1) ;
CrIndg = CrIndg(CrIndg < Tsim - burn - window - 5) ;

Crisis3 = (Simyn_fwd < -nstd*std(Simyn)).*(Simyn_bwd > nstd*std(Simyn)) ;
CrIndn = find(Crisis3 == 1) ;
CrIndn = CrIndn(CrIndn > window + 1) ;
CrIndn = CrIndn(CrIndn < Tsim - burn - window - 5) ;

% Boom-Bust 

clear BBFI

for i=-window:window+2
    BBFI.TAOSimZt(i + window + 1,:) = TAOSim(CrInd + i) ;
    BBFI.SimZt(i + window + 1,:) = Simyt(CrInd + i);
    BBFI.SimZt(i + window + 1,:) = Simyt(CrInd + i);
    BBFI.TAOSimGt(i + window + 1,:) = TAOSim(CrIndg + i);
    BBFI.SimGt(i + window + 1,:) = Simgt(CrIndg + i);
    BBFI.TAOSimNt(i + window + 1,:) = TAOSim(CrIndn + i);
    BBFI.SimNt(i + window + 1,:) = Simyn(CrIndn + i);
    BBFI.TAOSimYt(i + window + 1,:) = TAOSim(bbsInd + i + 1);
    BBFI.SimYt(i + window + 1,:) = Ytot(bbsInd + i);
end

BBFI.IRTAOSimZt = mean(BBFI.TAOSimZt, 2, "omitmissing");
BBFI.IRSimZt = mean(BBFI.SimZt, 2, "omitmissing");
BBFI.IRTAOSimGt = mean(BBFI.TAOSimGt, 2, "omitmissing");
BBFI.IRSimGt = mean(BBFI.SimGt, 2, "omitmissing");
BBFI.IRTAOSimNt = mean(BBFI.TAOSimNt, 2, "omitmissing");
BBFI.IRSimNt = mean(BBFI.SimNt, 2, "omitmissing");
BBFI.IRTAOSimYt = mean(BBFI.TAOSimYt, 2, "omitmissing");
BBFI.IRSimYt = mean(BBFI.SimYt, 2, "omitmissing");


%% Total Income: Boombust

close all

ff= figure('Position',Format.figsize2,'Color',[1 1 1]);
tt = tiledlayout(1,2);
ax=nexttile;
p1 = plot(-window:window+2, BBFI.IRTAOSimYt*100);
linei = 1;
set(p1,'LineStyle',Format.styles{linei})
set(p1,'color',Format.colors{linei})
set(p1,'linewidth',Format.widths{linei})
set(gca,'ycolor',Format.colors{linei})
%ylim([3,18])
xticks([-2,-1,0,1,2,3,4 ]) ;
xlabel('Years','FontSize',14, 'Interpreter','Latex')
ylabel('Tax Rate (\%)','FontSize',14, 'Interpreter','Latex')
title('Perfect Information', 'FontSize',16, 'Interpreter','Latex')
yyaxis right
p1 = plot(-window:window+2, BBFI.IRSimYt);
linei = 2;
set(p1,'LineStyle',Format.styles{linei})
set(p1,'color',Format.colors{linei})
set(p1,'linewidth',Format.widths{linei})
set(gca,'ycolor',Format.colors{linei})
ylabel('Total Income (Level) ','FontSize',14, 'Interpreter','Latex')
%ylim([0.85,1.2])
ax=nexttile;
p1 = plot(-window:window+2, BBII.IRTAOYt*100);
linei = 1;
set(p1,'LineStyle',Format.styles{linei})
set(p1,'color',Format.colors{linei})
set(p1,'linewidth',Format.widths{linei})
set(gca,'ycolor',Format.colors{linei})
ylim([3,18])
xticks([-2,-1,0,1,2,3,4 ]) ;
xlabel('Years','FontSize',14, 'Interpreter','Latex')
ylabel('Tax Rate (\%)','FontSize',14, 'Interpreter','Latex')
title('Imperfect Information', 'FontSize',16, 'Interpreter','Latex')

yyaxis right
p2 = plot(-window:window+2, BBII.IRYt);
linei = 2;
set(p2,'LineStyle',Format.styles{linei})
set(p2,'color',Format.colors{linei})
set(p2,'linewidth',Format.widths{linei})
set(gca,'ycolor',Format.colors{linei})
%ylim([0.85,1.2])
xticks([-2,-1,0,1,2,3,4 ]) ;
ylabel('Total Income (Level)','FontSize',14, 'Interpreter','Latex')
lg  = legend(ax,'Optimal Tax (\%, Left-hand Axis)', 'Underlying Income Component (Right-hand Axis)', ...
             'Orientation','Horizontal','NumColumns',2,'FontSize', Format.FontSize,'Interpreter','Latex');
lg.Layout.Tile = 'South';
% add a bit space to the figure
fig = gcf;
fig.Position(3) = fig.Position(3) + 375;
fig.Position(4) = fig.Position(4) + 75;

% Save figure
saveas(ff,'Figure10','png');
