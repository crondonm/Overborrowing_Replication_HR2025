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

load('../Replication/Data/IICEsim.mat')
fprintf("Database loaded... \n")

g =  Param.g;  % Mean growth rate of permanent component

Tsim = Param.Tsim;   % Simulation points
burn = Param.burn; % Burn-in period for simulation
nstd = Param.nstd;
window = Param.window;

% Settings for figures

Format.FontSize = 19;
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

Simyt = Sim(burn+2:end,1)' ; Simyn = Sim(burn+2:end,2)' ; Simg = Sim(burn+2:end,3)' ; 

ZBBT = -nstd*std(Simyt); % Standard Deviation to consider a Boom Bust
ZBBN = -nstd*std(Simyn);
GBBT = -nstd*std(Simg);

Crisis = (Simyt < ZBBT);
FreqZT = sum(Crisis)/(length(Simyt));
CrIndZT = find(Crisis == 1) ;
CrIndZT = CrIndZT(CrIndZT > window) ; 
CrIndZT = CrIndZT(CrIndZT < Tsim - burn - window) ;

CrisisZN = (Simyn < ZBBN);
FreqZN = sum(CrisisZN)/(length(Simyn));
CrIndZN = find(CrisisZN == 1) ;
CrIndZN = CrIndZN(CrIndZN > window) ; 
CrIndZN = CrIndZN(CrIndZN < Tsim - burn - window) ;

CrisisGT = (Simg < GBBT);
FreqGT = sum(CrisisGT)/(length(Simg));
CrIndGT = find(CrisisGT == 1) ;
CrIndGT = CrIndGT(CrIndGT > window) ; 
CrIndGT = CrIndGT(CrIndGT < Tsim - burn - window) ;


clear BB

for i=-window:window
    BB.ZT.IRZ(i + window + 1,:) = yt(Posterioryt(CrIndZT + i));
    BB.ZT.IRZN(i + window + 1,:) = yn(Posterioryn(CrIndZT + i));
    BB.ZT.IRG(i + window + 1,:) =  gt(Posteriorg(CrIndZT + i)) + g;       
    BB.ZT.IRZSim(i + window + 1,:) = Simyt(CrIndZT + i);
    BB.ZT.IRZNSim(i + window + 1,:) = Simyn(CrIndZT + i);
    BB.ZT.IRGSim(i + window + 1,:) = Simg(CrIndZT + i) + g;

    BB.ZN.IRZ(i + window + 1,:) = yt(Posterioryt(CrIndZN + i));
    BB.ZN.IRZN(i + window + 1,:) = yn(Posterioryn(CrIndZN + i));
    BB.ZN.IRG(i + window + 1,:) =  gt(Posteriorg(CrIndZN + i)) + g;       
    BB.ZN.IRZSim(i + window + 1,:) = Simyt(CrIndZN + i);
    BB.ZN.IRZNSim(i + window + 1,:) = Simyn(CrIndZN + i);
    BB.ZN.IRGSim(i + window + 1,:) = Simg(CrIndZN + i) + g;

    BB.GT.IRZ(i + window + 1,:) = yt(Posterioryt(CrIndGT + i));
    BB.GT.IRZN(i + window + 1,:) = yn(Posterioryn(CrIndGT + i));
    BB.GT.IRG(i + window + 1,:) =  gt(Posteriorg(CrIndGT + i)) + g;       
    BB.GT.IRZSim(i + window + 1,:) = Simyt(CrIndGT + i);
    BB.GT.IRZNSim(i + window + 1,:) = Simyn(CrIndGT + i);
    BB.GT.IRGSim(i + window + 1,:) = Simg(CrIndGT + i) + g;

end 


BB.ZT.IRZMean = mean(BB.ZT.IRZ, 2);
BB.ZT.IRZNMean = mean(BB.ZT.IRZN, 2);
BB.ZT.IRGMean = mean(BB.ZT.IRG, 2);
BB.ZT.IRZSimMean = mean(BB.ZT.IRZSim,2);
BB.ZT.IRZNSimMean = mean(BB.ZT.IRZNSim,2);
BB.ZT.IRGSimMean = mean(BB.ZT.IRGSim,2);

BB.ZN.IRZMean = mean(BB.ZN.IRZ, 2);
BB.ZN.IRZNMean = mean(BB.ZN.IRZN, 2);
BB.ZN.IRGMean = mean(BB.ZN.IRG, 2);
BB.ZN.IRZSimMean = mean(BB.ZN.IRZSim,2);
BB.ZN.IRZNSimMean = mean(BB.ZN.IRZNSim,2);
BB.ZN.IRGSimMean = mean(BB.ZN.IRGSim,2);

BB.GT.IRZMean = mean(BB.GT.IRZ, 2);
BB.GT.IRZNMean = mean(BB.GT.IRZN, 2);
BB.GT.IRGMean = mean(BB.GT.IRG, 2);
BB.GT.IRZSimMean = mean(BB.GT.IRZSim,2);
BB.GT.IRZNSimMean = mean(BB.GT.IRZNSim,2);
BB.GT.IRGSimMean = mean(BB.GT.IRGSim,2);


%% Figure 1:
% Response of Beliefs to a negative shock to the permanent growth component

f1 = figure('Position',Format.figsize6,'Color',[1 1 1]);
subplot(3,3,1)
p1 = plot(-window:window, BB.ZT.IRZMean,-window:window, BB.ZT.IRZSimMean); 
    for linei = 1:2
        set(p1(linei),'LineStyle',Format.styles{linei})
        set(p1(linei),'color',Format.colors{linei})
        set(p1(linei),'linewidth',Format.widths{linei})
    end
title('$Z_t^T$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;

subplot(3,3,2)
p2 = plot(-window:window, BB.ZT.IRZNMean,-window:window, BB.ZT.IRZNSimMean); 
    for linei = 1:2
        set(p2(linei),'LineStyle',Format.styles{linei})
        set(p2(linei),'color',Format.colors{linei})
        set(p2(linei),'linewidth',Format.widths{linei})
    end
title('$Z_t^N$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;

subplot(3,3,3)
p3 = plot(-window:window, BB.ZT.IRGMean,-window:window, (g +0*BB.ZT.IRGSimMean)); 
    for linei = 1:2
        set(p3(linei),'LineStyle',Format.styles{linei})
        set(p3(linei),'color',Format.colors{linei})
        set(p3(linei),'linewidth',Format.widths{linei})
    end
title('$g_t$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
subplot(3,3,4)
p4 = plot(-window:window, BB.ZN.IRZMean,-window:window, BB.ZN.IRZSimMean); 
    for linei = 1:2
        set(p4(linei),'LineStyle',Format.styles{linei})
        set(p4(linei),'color',Format.colors{linei})
        set(p4(linei),'linewidth',Format.widths{linei})
    end
title('$Z_t^T$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;

subplot(3,3,5)
p5 = plot(-window:window, BB.ZN.IRZNMean,-window:window, BB.ZN.IRZNSimMean); 
    for linei = 1:2
        set(p5(linei),'LineStyle',Format.styles{linei})
        set(p5(linei),'color',Format.colors{linei})
        set(p5(linei),'linewidth',Format.widths{linei})
    end
title('$Z_t^N$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
subplot(3,3,6)
p6 = plot(-window:window, BB.ZN.IRGMean,-window:window, (g+0.*BB.ZN.IRGSimMean)); 
    for linei = 1:2
        set(p6(linei),'LineStyle',Format.styles{linei})
        set(p6(linei),'color',Format.colors{linei})
        set(p6(linei),'linewidth',Format.widths{linei})
    end
title('$g_t$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
subplot(3,3,7)
p7 = plot(-window:window, BB.GT.IRZMean ,-window:window, 0*BB.GT.IRZSimMean); 
    for linei = 1:2
        set(p7(linei),'LineStyle',Format.styles{linei})
        set(p7(linei),'color',Format.colors{linei})
        set(p7(linei),'linewidth',Format.widths{linei})
    end
title('$Z_t^T$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
subplot(3,3,8)
p8 = plot(-window:window, BB.GT.IRZNMean ,-window:window, 0.*BB.GT.IRZNSimMean); 
    for linei = 1:2
        set(p8(linei),'LineStyle',Format.styles{linei})
        set(p8(linei),'color',Format.colors{linei})
        set(p8(linei),'linewidth',Format.widths{linei})
    end
title('$Z_t^N$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
subplot(3,3,9)
p9 = plot(-window:window, BB.GT.IRGMean ,-window:window, (BB.GT.IRGSimMean)); 
    for linei = 1:2
        set(p9(linei),'LineStyle',Format.styles{linei})
        set(p9(linei),'color',Format.colors{linei})
        set(p9(linei),'linewidth',Format.widths{linei})
    end
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
title('$g_t$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
legend('Posterior','True', 'Location', 'SouthEast', 'Orientation','horizontal', 'FontSize', 20,'Interpreter','latex');

% Create textbox
annotation(f1,'textbox',...
    [0.0910137363010006 0.735 0.17 0.0443690906343448],...
    'String','Shock to $Z_t^T$',...
    'Rotation',90,...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',22,...
    'EdgeColor','none');

annotation(f1,'textbox',...
    [0.0910137363010006 0.43 0.18 0.0443690906343448],...
    'String','Shock to $Z_t^N$',...
    'Rotation',90,...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',22,...
    'EdgeColor','none');

annotation(f1,'textbox',...
    [0.0910137363010006 0.135 0.17 0.0443690906343448],...
    'String','Shock to $g_t$',...
    'Rotation',90,...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',22,...
    'EdgeColor','none');

% add a bit space to the figure
fig = gcf;
fig.Position(3) = fig.Position(3) - 100;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.384;
Lgnd.Position(2) = 0.02651;

saveas(f1, 'Figure1', 'png');

