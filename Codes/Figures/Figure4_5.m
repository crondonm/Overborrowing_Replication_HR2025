%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% In this code: Replicate Figures 4 and 5 in the paper
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
load('../Replication/Data/IICEsim.mat')

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

fprintf("Starting crises analysis ... \n")

Simyt = Sim(burn+2:end,1)' ; Simyn = Sim(burn+2:end,2)' ; Simg = Sim(burn+2:end,3)' ; 

AAA = b(SimB) > BCSim + (b(2) - b(1))/2 ; 
CCC =CA;
CCCT = nstd*std(CCC);
Crisis = (CCC > CCCT).*(1 - AAA) ;
Crisis2 = sum(Crisis)/(length(CCC)) ;
CrInd = find(Crisis == 1) ;
CrInd = CrInd(CrInd > window + 1) ; 
CrInd = CrInd(CrInd < Tsim - burn - window) ;  

clear BBII 

BBII.CE.mP = mP;
BBII.CE.mCT = mCT;
BBII.CE.mC = mean(CSim);

for i=-window:window
    BBII.CE.IRB(i + window + 1,:) = b(SimB(CrInd + i));
    BBII.CE.IRCA(i + window + 1,:) = CA(CrInd + i);
    BBII.CE.IRBC(i + window + 1,:) = BCSim(CrInd + i);
    BBII.CE.IRC(i + window + 1,:) = CSim(CrInd + i);
    BBII.CE.IRCT(i + window + 1,:) = CTSim(CrInd + i);
    BBII.CE.IRCN(i + window + 1,:) = CNSim(CrInd + i - 1);
    BBII.CE.IRZ(i + window + 1,:) = yt(Posterioryt(CrInd + i));
    BBII.CE.IRZN(i + window + 1,:) = yn(Posterioryn(CrInd + i));
    BBII.CE.IRG(i + window + 1,:) = gt(Posteriorg(CrInd + i)) + g;    
    BBII.CE.IRZSim(i + window + 1,:) = Simyt(CrInd + i);
    BBII.CE.IRZNSim(i + window + 1,:) = Simyn(CrInd + i - 1);
    BBII.CE.IRGSim(i + window + 1,:) = Simg(CrInd + i) + g;    
    BBII.CE.IRP(i + window + 1,:) = PSim(CrInd + i);
    BBII.CE.IRDtoY(i + window + 1,:) = DtoY(CrInd + i); 
    BBII.CE.IRCtoY(i + window + 1,:) = CtoY(CrInd + i + 1);
    BBII.CE.IRCTtoY(i + window + 1,:) = CTtoY(CrInd + i);
    BBII.CE.IRCNtoY(i + window + 1,:) = CNtoY(CrInd + i);
    BBII.CE.IRCAtoY(i + window + 1,:) = CAtoY(CrInd + i);
    BBII.CE.IRYtot(i + window + 1,:) = Ytot(CrInd + i);
end 
 
BBII.CE.IRBMean = mean(BBII.CE.IRB, 2);
BBII.CE.IRCAMean = mean(BBII.CE.IRCA,2);
BBII.CE.IRBCMean = mean(BBII.CE.IRBC,2);
BBII.CE.IRCMean = mean(BBII.CE.IRC,2);
BBII.CE.IRCTMean = mean(BBII.CE.IRCT,2);
BBII.CE.IRCNMean = mean(BBII.CE.IRCN,2);
BBII.CE.IRPMean = mean(BBII.CE.IRP, 2);
BBII.CE.IRDtoYMean = mean(BBII.CE.IRDtoY,2);
BBII.CE.IRCtoYMean = mean(BBII.CE.IRCtoY,2);
BBII.CE.IRCTtoYMean = mean(BBII.CE.IRCTtoY,2);
BBII.CE.IRCNtoYMean = mean(BBII.CE.IRCNtoY,2);
BBII.CE.IRCAtoYMean = mean(BBII.CE.IRCAtoY,2);
BBII.CE.IRYtotMean = mean(BBII.CE.IRYtot,2);
BBII.CE.IRZMean = mean(BBII.CE.IRZ, 2);
BBII.CE.IRZNMean = mean(BBII.CE.IRZN, 2);
BBII.CE.IRGMean = mean(BBII.CE.IRG, 2);
BBII.CE.IRZSimMean = mean(BBII.CE.IRZSim, 2);
BBII.CE.IRZNSimMean = mean(BBII.CE.IRZNSim, 2);
BBII.CE.IRGSimMean = mean(BBII.CE.IRGSim, 2);

%% II: Planner 

clearvars -except BBII g Tsim burn nstd window Format 

load('../Replication/Data/IIPCCsim.mat')

fprintf("Starting crises analysis ... \n")

Simyt = Sim(burn+2:end,1)' ; Simyn = Sim(burn+2:end,2)' ; Simg = Sim(burn+2:end,3)' ; 

AAA = b(SimB) > BCSim + (b(2) - b(1))/2 ; 
CCC =CA;
CCCT = nstd*std(CCC);
Crisis = (CCC > CCCT).*(1 - AAA) ;
Crisis2 = sum(Crisis)/(length(CCC)) ;
CrInd = find(Crisis == 1) ;
CrInd = CrInd(CrInd > window + 1) ; 
CrInd = CrInd(CrInd < Tsim - burn - window) ;  

BBII.PP.mP = mP;
BBII.PP.mCT = mCT;
BBII.PP.mC = mean(CSim);

for i=-window:window
    BBII.PP.IRB(i + window + 1,:) = b(SimB(CrInd + i));
    BBII.PP.IRCA(i + window + 1,:) = CA(CrInd + i);
    BBII.PP.IRBC(i + window + 1,:) = BCSim(CrInd + i);
    BBII.PP.IRC(i + window + 1,:) = CSim(CrInd + i);
    BBII.PP.IRCT(i + window + 1,:) = CTSim(CrInd + i);
    BBII.PP.IRZSim(i + window + 1,:) = Simyt(CrInd + i);
    BBII.PP.IRZNSim(i + window + 1,:) = Simyn(CrInd + i - 1);
    BBII.PP.IRGSim(i + window + 1,:) = Simg(CrInd + i) + g;    
    BBII.PP.IRP(i + window + 1,:) = PSim(CrInd + i);
    BBII.PP.IRDtoY(i + window + 1,:) = DtoY(CrInd + i); 
    BBII.PP.IRCtoY(i + window + 1,:) = CtoY(CrInd + i + 1);
    BBII.PP.IRCTtoY(i + window + 1,:) = CTtoY(CrInd + i);
    BBII.PP.IRCNtoY(i + window + 1,:) = CNtoY(CrInd + i);
    BBII.PP.IRCAtoY(i + window + 1,:) = CAtoY(CrInd + i);
end 
 
BBII.PP.IRBMean = mean(BBII.PP.IRB, 2);
BBII.PP.IRCAMean = mean(BBII.PP.IRCA,2);
BBII.PP.IRBCMean = mean(BBII.PP.IRBC,2);
BBII.PP.IRCMean = mean(BBII.PP.IRC,2);
BBII.PP.IRCTMean = mean(BBII.PP.IRCT,2);
BBII.PP.IRPMean = mean(BBII.PP.IRP, 2);
BBII.PP.IRDtoYMean = mean(BBII.PP.IRDtoY,2);
BBII.PP.IRCtoYMean = mean(BBII.PP.IRCtoY,2);
BBII.PP.IRCTtoYMean = mean(BBII.PP.IRCTtoY,2);
BBII.PP.IRCNtoYMean = mean(BBII.PP.IRCNtoY,2);
BBII.PP.IRCAtoYMean = mean(BBII.PP.IRCAtoY,2);
BBII.PP.IRZSimMean = mean(BBII.PP.IRZSim, 2);
BBII.PP.IRZNSimMean = mean(BBII.PP.IRZNSim, 2);
BBII.PP.IRGSimMean = mean(BBII.PP.IRGSim, 2);

%% Full Information

clearvars -except BBII g Tsim burn nstd window Format 

load('../Replication/Data/FICEsim.mat')

fprintf("Starting crises analysis Full Information \n")

AAA = b(SimB) > BCSim + (b(2) - b(1))/2 ; 
CCC =CA;
CCCT = nstd*std(CCC);
Crisis = (CCC > CCCT).*(1 - AAA) ;
Crisis2 = sum(Crisis)/(length(CCC)) ;
%Crisis = [0 Crisis]; 
CrInd = find(Crisis == 1) ;
CrInd = CrInd(CrInd > window + 1) ; 
CrInd = CrInd(CrInd < Tsim - burn - window) ;  

clear BBFI 

BBFI.CE.mP = mP;
BBFI.CE.mCT = mCT;
BBFI.CE.mC = mean(CSim);

for i=-window:window
    BBFI.CE.IRB(i + window + 1,:) = SimBhat(CrInd + i);
    BBFI.CE.IRCA(i + window + 1,:) = CA(CrInd + i) ;
    BBFI.CE.IRBC(i + window + 1,:) = BCSim(CrInd + i) ;
    BBFI.CE.IRC(i + window + 1,:) = CSim(CrInd + i) ;
    BBFI.CE.IRCT(i + window + 1,:) = CTSim(CrInd + i) ;
    BBFI.CE.IRP(i + window + 1,:) = PSim(CrInd + i) ;
    BBFI.CE.IRDtoY(i + window + 1,:) = DtoY(CrInd + i + 1); 
    BBFI.CE.IRLambda(i + window + 1,:) = LambdaSim(CrInd + i);
    BBFI.CE.IRCtoY(i + window + 1,:) = CtoY(CrInd + i + 1);
    BBFI.CE.IRCAtoY(i + window + 1,:) = CAtoY(CrInd + i);
    BBFI.CE.IRCTtoY(i + window + 1,:) = CTtoY(CrInd + i);
    BBFI.CE.IRCNtoY(i + window + 1,:) = CNtoY(CrInd + i);
    BBFI.CE.IRYtot(i + window + 1,:) = Ytot(CrInd + i);
end 
 
BBFI.CE.IRBMean = mean(BBFI.CE.IRB, 2);
BBFI.CE.IRCAMean = mean(BBFI.CE.IRCA,2);
BBFI.CE.IRBCMean = mean(BBFI.CE.IRBC,2);
BBFI.CE.IRCMean = mean(BBFI.CE.IRC,2);
BBFI.CE.IRCTMean = mean(BBFI.CE.IRCT,2);
BBFI.CE.IRPMean = mean(BBFI.CE.IRP, 2);
BBFI.CE.IRDtoYMean = mean(BBFI.CE.IRDtoY,2);
BBFI.CE.IRCtoYMean = mean(BBFI.CE.IRCtoY,2);
BBFI.CE.IRCTtoYMean = mean(BBFI.CE.IRCTtoY,2);
BBFI.CE.IRCNtoYMean = mean(BBFI.CE.IRCNtoY,2);
BBFI.CE.IRCAtoYMean = mean(BBFI.CE.IRCAtoY,2);
BBFI.CE.IRYtotMean = mean(BBFI.CE.IRYtot,2);

%% FI: Planner

clearvars -except BBII BBFI g Tsim burn nstd window Format 

load('../Replication/Data/FIPsim.mat')

fprintf("Starting crises analysis Full Information \n")

AAA = b(SimB) > BCSim + (b(2) - b(1))/2 ; 
CCC =CA;
CCCT = nstd*std(CCC);
Crisis = (CCC > CCCT).*(1 - AAA) ;
Crisis2 = sum(Crisis)/(length(CCC)) ;
%Crisis = [0 Crisis]; 
CrInd = find(Crisis == 1) ;
CrInd = CrInd(CrInd > window + 1) ; 
CrInd = CrInd(CrInd < Tsim - burn - window) ;  

BBFI.PP.mP = mP;
BBFI.PP.mCT = mCT;
BBFI.PP.mC = mean(CSim);

for i=-window:window
    BBFI.PP.IRB(i + window + 1,:) = SimBhat(CrInd + i);
    BBFI.PP.IRCA(i + window + 1,:) = CA(CrInd + i) ;
    BBFI.PP.IRBC(i + window + 1,:) = BCSim(CrInd + i) ;
    BBFI.PP.IRC(i + window + 1,:) = CSim(CrInd + i) ;
    BBFI.PP.IRCT(i + window + 1,:) = CTSim(CrInd + i) ;
    BBFI.PP.IRP(i + window + 1,:) = PSim(CrInd + i) ;
    BBFI.PP.IRDtoY(i + window + 1,:) = DtoY(CrInd + i + 1); 
    BBFI.PP.IRCtoY(i + window + 1,:) = CtoY(CrInd + i + 1);
    BBFI.PP.IRCAtoY(i + window + 1,:) = CAtoY(CrInd + i);
    BBFI.PP.IRCTtoY(i + window + 1,:) = CTtoY(CrInd + i);
    BBFI.PP.IRCNtoY(i + window + 1,:) = CNtoY(CrInd + i);
end 
 
BBFI.PP.IRBMean = mean(BBFI.PP.IRB, 2);
BBFI.PP.IRCAMean = mean(BBFI.PP.IRCA,2);
BBFI.PP.IRBCMean = mean(BBFI.PP.IRBC,2);
BBFI.PP.IRCMean = mean(BBFI.PP.IRC,2);
BBFI.PP.IRCTMean = mean(BBFI.PP.IRCT,2);
BBFI.PP.IRPMean = mean(BBFI.PP.IRP, 2);
BBFI.PP.IRDtoYMean = mean(BBFI.PP.IRDtoY,2);
BBFI.PP.IRCtoYMean = mean(BBFI.PP.IRCtoY,2);
BBFI.PP.IRCTtoYMean = mean(BBFI.PP.IRCTtoY,2);
BBFI.PP.IRCNtoYMean = mean(BBFI.PP.IRCNtoY,2);
BBFI.PP.IRCAtoYMean = mean(BBFI.PP.IRCAtoY,2);

%% Figure 4

% Main Drivers: Imperfect Information

BBII.CE.IRZMean = mean(BBII.CE.IRZ, 2);
BBII.CE.IRZNMean = mean(BBII.CE.IRZN, 2);
BBII.CE.IRGMean = mean(BBII.CE.IRG, 2);
BBII.CE.IRZSimMean = mean(BBII.CE.IRZSim, 2);
BBII.CE.IRZNSimMean = mean(BBII.CE.IRZNSim, 2);
BBII.CE.IRGSimMean = mean(BBII.CE.IRGSim, 2);

f4= figure('Position',Format.figsize3,'Color',[1 1 1]);
tt = tiledlayout(1,3);
ax = nexttile;
p1 = plot(-window:window, BBII.CE.IRZMean, -window:window, BBII.CE.IRZSimMean);
   for linei = 1:2
       set(p1(linei),'LineStyle',Format.styles{linei})
       set(p1(linei),'color',Format.colors{linei})
       set(p1(linei),'linewidth',Format.widths{linei})
   end
%ylim([-11 1]) ;
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
ylabel('Level')
title('$Z_t^T$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
ax = nexttile;
p2 = plot(-window:window, BBII.CE.IRZNMean, -window:window, BBII.CE.IRZNSimMean);
   for linei = 1:2
       set(p2(linei),'LineStyle',Format.styles{linei})
       set(p2(linei),'color',Format.colors{linei})
       set(p2(linei),'linewidth',Format.widths{linei})
   end
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
%ylim([-1 3]) ;
ylabel('Level')
title('$Z_t^N$', 'FontSize', Format.FontSize, 'FontWeight', Format.fontweight,'Interpreter','latex')
ax = nexttile;
p3 = plot(-window:window, (BBII.CE.IRGMean)*100, -window:window, (BBII.CE.IRGSimMean)*100);
for linei = 1:2
       set(p3(linei),'LineStyle',Format.styles{linei})
       set(p3(linei),'color',Format.colors{linei})
       set(p3(linei),'linewidth',Format.widths{linei})
end
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
%ylim([-8 1]) ;
title('$g_t$', 'FontSize', Format.FontSize, 'FontWeight', Format.fontweight,'Interpreter','latex')
ylabel('Percent')
lg  = legend(ax,'Posterior', 'True', 'Orientation','Horizontal','NumColumns',2);
lg.Layout.Tile = 'South';
% add a bit space to the figure
fig = gcf;
fig.Position(3) = fig.Position(3) + 175;
fig.Position(4) = fig.Position(4) + 50;

%saveas(f4,'Figure4_bb_drivers','png');

filename = 'Figure4.png';
resolution = 300; % DPI
% Export the graphics
exportgraphics(gcf, filename, 'Resolution', resolution);


%% Figure 5

close all

f5 = figure('Position',Format.figsize4,'Color',[1 1 1]);
subplot(2,2,1)
p1 = plot(-window:window, 100*(BBFI.CE.IRCTMean/BBFI.CE.mCT-1), -window:window, 100*((BBII.CE.IRCTMean+BBII.CE.mCT)/BBII.CE.mCT-1)); 
    for linei = 1:2
        set(p1(linei),'LineStyle',Format.styles{linei})
        set(p1(linei),'color',Format.colors{linei})
        set(p1(linei),'linewidth',Format.widths{linei})
    end
title('$C_t^T$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
ylabel('Percent', 'FontSize',14,'FontWeight',  Format.fontweight,'Interpreter','latex')
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;

subplot(2,2,2)
p2 = plot(-window:window, 100*(BBFI.CE.IRPMean/BBFI.PP.mP-1),-window:window, 100*((BBII.CE.IRPMean+BBII.CE.mP)/BBII.CE.mP-1)); 
    for linei = 1:2
        set(p2(linei),'LineStyle',Format.styles{linei})
        set(p2(linei),'color',Format.colors{linei})
        set(p2(linei),'linewidth',Format.widths{linei})
    end
title('$P_t$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
ylabel('Percent', 'FontSize',14,'FontWeight',  Format.fontweight,'Interpreter','latex')
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')

subplot(2,2,3)
p3 = plot(-window:window, BBFI.CE.IRBMean,-window:window, BBII.CE.IRBMean); 
    for linei = 1:2
        set(p3(linei),'LineStyle',Format.styles{linei})
        set(p3(linei),'color',Format.colors{linei})
        set(p3(linei),'linewidth',Format.widths{linei})
    end
title('$B_{t+1}$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xlim
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
ylabel('Level', 'FontSize',14,'FontWeight',  Format.fontweight,'Interpreter','latex')
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')
subplot(2,2,4)
p3 = plot(-window:window, BBFI.CE.IRBCMean ,-window:window, BBII.CE.IRBCMean); 
    for linei = 1:2
        set(p3(linei),'LineStyle',Format.styles{linei})
        set(p3(linei),'color',Format.colors{linei})
        set(p3(linei),'linewidth',Format.widths{linei})
    end
title('Borrowing Limit','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
ylabel('Level', 'FontSize',14,'FontWeight',  Format.fontweight,'Interpreter','latex')
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')
legend('Perfect Information','Imperfect Information', 'Location', 'SouthEast', 'Orientation','horizontal', 'FontSize', 14,'Interpreter','latex');
% add a bit space to the figure
fig = gcf;
fig.Position(3) = fig.Position(3) + 50;
fig.Position(4) = fig.Position(4) + 100;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.2405  ;
Lgnd.Position(2) = 0.0101;

%saveas(f5, 'Figure5_responses_crises', 'png');

filename = 'Figure5.png';
resolution = 300; % DPI
% Export the graphics
exportgraphics(gcf, filename, 'Resolution', resolution);

