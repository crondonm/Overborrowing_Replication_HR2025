%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% In this code: Plot figure A2
% 
% Authors:  Juan Herreño, jherrenolopera@ucsd.edu
%               Carlos Rondón Moreno, crondon@bcentral.cl
%
% Last:  March 2025
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping

clearvars
clear global
close all

% Load databases

load('../../Replication/Data/Param.mat')
fprintf("Parameters loaded... \n")
load('../../Replication/Data/FICESim.mat')
fprintf("Database loaded... \n")

g =  Param.g;  % Mean growth rate of permanent component

Tsim = Param.Tsim;   % Simulation points
burn = Param.burn; % Burn-in period for simulation
nstd = Param.nstd;
window = Param.window;
FICEmP = mP;

% Settings for figures

Format.FontSize = 22;
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

ZBBT = nstd*std(Simyt); % Standard Deviation to consider a Boom Bust
ZBBN = nstd*std(Simyn);
GBBT = nstd*std(Simg);

Crisis = (Simyt > ZBBT);
FreqZT = sum(Crisis)/(length(Simyt));
CrIndZT = find(Crisis == 1) ;
CrIndZT = CrIndZT(CrIndZT > window) ; 
CrIndZT = CrIndZT(CrIndZT < Tsim - burn - window-5) ;

CrisisZN = (Simyn > ZBBN);
FreqZN = sum(CrisisZN)/(length(Simyn));
CrIndZN = find(CrisisZN == 1) ;
CrIndZN = CrIndZN(CrIndZN > window) ; 
CrIndZN = CrIndZN(CrIndZN < Tsim - burn - window-5) ;

CrisisGT = (Simg > GBBT);
FreqGT = sum(CrisisGT)/(length(Simg));
CrIndGT = find(CrisisGT == 1) ;
CrIndGT = CrIndGT(CrIndGT > window) ; 
CrIndGT = CrIndGT(CrIndGT < Tsim - burn - window-5) ;


clear BBFI

BBFI.mP = mP;
BBFI.mCT = mCT;


for i=-window:window

    BBFI.ZT.IRB(i + window + 1,:) = SimBhat(CrIndZT + i + 2);
    BBFI.ZT.IRCA(i + window + 1,:) = CA(CrIndZT + i + 2) ;
    BBFI.ZT.IRBC(i + window + 1,:) = BCSim(CrIndZT + i + 1) ;
    BBFI.ZT.IRC(i + window + 1,:) = CSim(CrIndZT + i + 2) ;
    BBFI.ZT.IRCT(i + window + 1,:) = CTSim(CrIndZT+ i + 1) ;
    BBFI.ZT.IRP(i + window + 1,:) = PSim(CrIndZT + i + 1) ;
    BBFI.ZT.IRDtoY(i + window + 1,:) = DtoY(CrIndZT + i + 2); 
    BBFI.ZT.IRCtoY(i + window + 1,:) = CtoY(CrIndZT + i + 2);
    BBFI.ZT.IRCAtoY(i + window + 1,:) = CAtoY(CrIndZT + i + 2);
    BBFI.ZT.IRCTtoY(i + window + 1,:) = CTtoY(CrIndZT + i + 2);
    BBFI.ZT.IRCNtoY(i + window + 1,:) = CNtoY(CrIndZT + i + 2);
    BBFI.ZT.IRYtot(i + window + 1,:) = Ytot(CrIndZT + i + 2);

    BBFI.ZN.IRB(i + window + 1,:) = SimBhat(CrIndZN + i + 1);
    BBFI.ZN.IRCA(i + window + 1,:) = CA(CrIndZN + i + 2) ;
    BBFI.ZN.IRBC(i + window + 1,:) = BCSim(CrIndZN + i + 2) ;
    BBFI.ZN.IRC(i + window + 1,:) = CSim(CrIndZN + i + 2) ;
    BBFI.ZN.IRCT(i + window + 1,:) = CTSim(CrIndZN+ i + 2) ;
    BBFI.ZN.IRP(i + window + 1,:) = PSim(CrIndZN + i + 1) ;
    BBFI.ZN.IRDtoY(i + window + 1,:) = DtoY(CrIndZN + i + 2); 
    BBFI.ZN.IRCtoY(i + window + 1,:) = CtoY(CrIndZN + i + 2);
    BBFI.ZN.IRCAtoY(i + window + 1,:) = CAtoY(CrIndZN + i + 2);
    BBFI.ZN.IRCTtoY(i + window + 1,:) = CTtoY(CrIndZN + i + 2);
    BBFI.ZN.IRCNtoY(i + window + 1,:) = CNtoY(CrIndZN + i + 2);
    BBFI.ZN.IRYtot(i + window + 1,:) = Ytot(CrIndZN + i + 2);

    BBFI.GT.IRB(i + window + 1,:) = SimBhat(CrIndGT + i + 1);
    BBFI.GT.IRCA(i + window + 1,:) = CA(CrIndGT + i + 2) ;
    BBFI.GT.IRBC(i + window + 1,:) = BCSim(CrIndGT + i + 1) ;
    BBFI.GT.IRC(i + window + 1,:) = CSim(CrIndGT + i + 2) ;
    BBFI.GT.IRCT(i + window + 1,:) = CTSim(CrIndGT+ i + 1) ;
    BBFI.GT.IRP(i + window + 1,:) = PSim(CrIndGT + i + 1) ;
    BBFI.GT.IRDtoY(i + window + 1,:) = DtoY(CrIndGT + i + 2); 
    BBFI.GT.IRCtoY(i + window + 1,:) = CtoY(CrIndGT + i + 2);
    BBFI.GT.IRCAtoY(i + window + 1,:) = CAtoY(CrIndGT + i + 2);
    BBFI.GT.IRCTtoY(i + window + 1,:) = CTtoY(CrIndGT + i + 2);
    BBFI.GT.IRCNtoY(i + window + 1,:) = CNtoY(CrIndGT + i + 2);
    BBFI.GT.IRYtot(i + window + 1,:) = Ytot(CrIndGT + i + 2);

end 
 
BBFI.ZT.IRBMean = mean(BBFI.ZT.IRB, 2);
BBFI.ZT.IRCAMean = mean(BBFI.ZT.IRCA,2);
BBFI.ZT.IRBCMean = mean(BBFI.ZT.IRBC,2);
BBFI.ZT.IRCMean = mean(BBFI.ZT.IRC,2);
BBFI.ZT.IRCTMean = mean(BBFI.ZT.IRCT,2);
BBFI.ZT.IRPMean = mean(BBFI.ZT.IRP, 2);
BBFI.ZT.IRDtoYMean = mean(BBFI.ZT.IRDtoY,2);
BBFI.ZT.IRCtoYMean = mean(BBFI.ZT.IRCtoY,2);
BBFI.ZT.IRCAtoYMean = mean(BBFI.ZT.IRCAtoY,2);
BBFI.ZT.IRCTtoYMean = mean(BBFI.ZT.IRCTtoY,2);
BBFI.ZT.IRCNtoYMean = mean(BBFI.ZT.IRCNtoY,2);
BBFI.ZT.IRYtotMean = mean(BBFI.ZT.IRYtot,2);

BBFI.ZN.IRBMean = mean(BBFI.ZN.IRB, 2);
BBFI.ZN.IRCAMean = mean(BBFI.ZN.IRCA,2);
BBFI.ZN.IRBCMean = mean(BBFI.ZN.IRBC,2);
BBFI.ZN.IRCMean = mean(BBFI.ZN.IRC,2);
BBFI.ZN.IRCTMean = mean(BBFI.ZN.IRCT,2);
BBFI.ZN.IRPMean = mean(BBFI.ZN.IRP, 2);
BBFI.ZN.IRDtoYMean = mean(BBFI.ZN.IRDtoY,2);
BBFI.ZN.IRCtoYMean = mean(BBFI.ZN.IRCtoY,2);
BBFI.ZN.IRCAtoYMean = mean(BBFI.ZN.IRCAtoY,2);
BBFI.ZN.IRCTtoYMean = mean(BBFI.ZN.IRCTtoY,2);
BBFI.ZN.IRCNtoYMean = mean(BBFI.ZN.IRCNtoY,2);
BBFI.ZN.IRYtotMean = mean(BBFI.ZN.IRYtot,2);

BBFI.GT.IRBMean = mean(BBFI.GT.IRB, 2);
BBFI.GT.IRCAMean = mean(BBFI.GT.IRCA,2);
BBFI.GT.IRBCMean = mean(BBFI.GT.IRBC,2);
BBFI.GT.IRCMean = mean(BBFI.GT.IRC,2);
BBFI.GT.IRCTMean = mean(BBFI.GT.IRCT,2);
BBFI.GT.IRPMean = mean(BBFI.GT.IRP, 2);
BBFI.GT.IRDtoYMean = mean(BBFI.GT.IRDtoY,2);
BBFI.GT.IRCtoYMean = mean(BBFI.GT.IRCtoY,2);
BBFI.GT.IRCAtoYMean = mean(BBFI.GT.IRCAtoY,2);
BBFI.GT.IRCTtoYMean = mean(BBFI.GT.IRCTtoY,2);
BBFI.GT.IRCNtoYMean = mean(BBFI.GT.IRCNtoY,2);
BBFI.GT.IRYtotMean = mean(BBFI.GT.IRYtot,2);

%% 

clearvars -except BBFI g Tsim burn nstd window CrIndGT CrIndZN CrIndZT Format FICEmP

load("../../Replication/Data/IICEsim.mat")

clear BBII

BBII.mP = mP;
BBII.mCT = mCT;

for i=-window:window

    BBII.ZT.IRB(i + window + 1,:) = b(SimB(CrIndZT + i + 1));
    BBII.ZT.IRCA(i + window + 1,:) = CA(CrIndZT + i + 2) ;
    BBII.ZT.IRBC(i + window + 1,:) = BCSim(CrIndZT + i + 1) ;
    BBII.ZT.IRC(i + window + 1,:) = CSim(CrIndZT + i + 2) ;
    BBII.ZT.IRCT(i + window + 1,:) = CTSim(CrIndZT+ i + 1) ;
    BBII.ZT.IRP(i + window + 1,:) = PSim(CrIndZT + i + 1) ;
    BBII.ZT.IRDtoY(i + window + 1,:) = DtoY(CrIndZT + i + 2); 
    BBII.ZT.IRCtoY(i + window + 1,:) = CtoY(CrIndZT + i + 2);
    BBII.ZT.IRCAtoY(i + window + 1,:) = CAtoY(CrIndZT + i + 2);
    BBII.ZT.IRCTtoY(i + window + 1,:) = CTtoY(CrIndZT + i + 2);
    BBII.ZT.IRCNtoY(i + window + 1,:) = CNtoY(CrIndZT + i + 2);
    BBII.ZT.IRYtot(i + window + 1,:) = Ytot(CrIndZT + i + 2);

    BBII.ZN.IRB(i + window + 1,:) = b(SimB(CrIndZN + i + 1));
    BBII.ZN.IRCA(i + window + 1,:) = CA(CrIndZN + i + 2) ;
    BBII.ZN.IRBC(i + window + 1,:) = BCSim(CrIndZN + i + 2) ;
    BBII.ZN.IRC(i + window + 1,:) = CSim(CrIndZN + i + 2) ;
    BBII.ZN.IRCT(i + window + 1,:) = CTSim(CrIndZN+ i + 2) ;
    BBII.ZN.IRP(i + window + 1,:) = PSim(CrIndZN + i + 2) ;
    BBII.ZN.IRDtoY(i + window + 1,:) = DtoY(CrIndZN + i + 2); 
    BBII.ZN.IRCtoY(i + window + 1,:) = CtoY(CrIndZN + i + 2);
    BBII.ZN.IRCAtoY(i + window + 1,:) = CAtoY(CrIndZN + i + 2);
    BBII.ZN.IRCTtoY(i + window + 1,:) = CTtoY(CrIndZN + i + 2);
    BBII.ZN.IRCNtoY(i + window + 1,:) = CNtoY(CrIndZN + i + 2);
    BBII.ZN.IRYtot(i + window + 1,:) = Ytot(CrIndZN + i + 2);

    BBII.GT.IRB(i + window + 1,:) = b(SimB(CrIndGT + i + 1));
    BBII.GT.IRCA(i + window + 1,:) = CA(CrIndGT + i + 2) ;
    BBII.GT.IRBC(i + window + 1,:) = BCSim(CrIndGT + i + 1) ;
    BBII.GT.IRC(i + window + 1,:) = CSim(CrIndGT + i + 2) ;
    BBII.GT.IRCT(i + window + 1,:) = CTSim(CrIndGT+ i + 1) ;
    BBII.GT.IRP(i + window + 1,:) = PSim(CrIndGT + i + 1) ;
    BBII.GT.IRDtoY(i + window + 1,:) = DtoY(CrIndGT + i + 2); 
    BBII.GT.IRCtoY(i + window + 1,:) = CtoY(CrIndGT + i + 2);
    BBII.GT.IRCAtoY(i + window + 1,:) = CAtoY(CrIndGT + i + 2);
    BBII.GT.IRCTtoY(i + window + 1,:) = CTtoY(CrIndGT + i + 2);
    BBII.GT.IRCNtoY(i + window + 1,:) = CNtoY(CrIndGT + i + 2);
    BBII.GT.IRYtot(i + window + 1,:) = Ytot(CrIndGT + i + 2);

end 

BBII.ZT.IRBMean = mean(BBII.ZT.IRB, 2);
BBII.ZT.IRCAMean = mean(BBII.ZT.IRCA,2);
BBII.ZT.IRBCMean = mean(BBII.ZT.IRBC,2);
BBII.ZT.IRCMean = mean(BBII.ZT.IRC,2);
BBII.ZT.IRCTMean = mean(BBII.ZT.IRCT,2);
BBII.ZT.IRPMean = mean(BBII.ZT.IRP, 2);
BBII.ZT.IRDtoYMean = mean(BBII.ZT.IRDtoY,2);
BBII.ZT.IRCtoYMean = mean(BBII.ZT.IRCtoY,2);
BBII.ZT.IRCAtoYMean = mean(BBII.ZT.IRCAtoY,2);
BBII.ZT.IRCTtoYMean = mean(BBII.ZT.IRCTtoY,2);
BBII.ZT.IRCNtoYMean = mean(BBII.ZT.IRCNtoY,2);
BBII.ZT.IRYtotMean = mean(BBII.ZT.IRYtot,2);

BBII.ZN.IRBMean = mean(BBII.ZN.IRB, 2);
BBII.ZN.IRCAMean = mean(BBII.ZN.IRCA,2);
BBII.ZN.IRBCMean = mean(BBII.ZN.IRBC,2);
BBII.ZN.IRCMean = mean(BBII.ZN.IRC,2);
BBII.ZN.IRCTMean = mean(BBII.ZN.IRCT,2);
BBII.ZN.IRPMean = mean(BBII.ZN.IRP, 2);
BBII.ZN.IRDtoYMean = mean(BBII.ZN.IRDtoY,2);
BBII.ZN.IRCtoYMean = mean(BBII.ZN.IRCtoY,2);
BBII.ZN.IRCAtoYMean = mean(BBII.ZN.IRCAtoY,2);
BBII.ZN.IRCTtoYMean = mean(BBII.ZN.IRCTtoY,2);
BBII.ZN.IRCNtoYMean = mean(BBII.ZN.IRCNtoY,2);
BBII.ZN.IRYtotMean = mean(BBII.ZN.IRYtot,2);

BBII.GT.IRBMean = mean(BBII.GT.IRB, 2);
BBII.GT.IRCAMean = mean(BBII.GT.IRCA,2);
BBII.GT.IRBCMean = mean(BBII.GT.IRBC,2);
BBII.GT.IRCMean = mean(BBII.GT.IRC,2);
BBII.GT.IRCTMean = mean(BBII.GT.IRCT,2);
BBII.GT.IRPMean = mean(BBII.GT.IRP, 2);
BBII.GT.IRDtoYMean = mean(BBII.GT.IRDtoY,2);
BBII.GT.IRCtoYMean = mean(BBII.GT.IRCtoY,2);
BBII.GT.IRCAtoYMean = mean(BBII.GT.IRCAtoY,2);
BBII.GT.IRCTtoYMean = mean(BBII.GT.IRCTtoY,2);
BBII.GT.IRCNtoYMean = mean(BBII.GT.IRCNtoY,2);
BBII.GT.IRYtotMean = mean(BBII.GT.IRYtot,2);



%% Figure 1 Panel A:
% Response of Beliefs to a negative shock to the permanent growth component

f1 = figure('Position',Format.figsize6,'Color',[1 1 1]);
subplot(3,4,1)
p1 = plot(-window:window, 100*((BBFI.ZT.IRCTMean)/BBFI.mCT - 1),-window:window, ((BBII.ZT.IRCTMean+BBII.mCT)/BBII.mCT - 1)*100); 
    for linei = 1:2
        set(p1(linei),'LineStyle',Format.styles{linei})
        set(p1(linei),'color',Format.colors{linei})
        set(p1(linei),'linewidth',Format.widths{linei})
    end
title('$C_t^T$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
ylabel('Percent', 'FontSize', Format.FontSizeAxes, 'FontWeight', Format.fontweight,'Interpreter','latex')
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]);

subplot(3,4,2)
p2 = plot(-window:window, (BBFI.ZT.IRPMean/BBFI.mP-1)*100,-window:window, ((BBII.ZT.IRPMean+BBII.mP)/BBII.mP-1)*100); 
    for linei = 1:2
        set(p2(linei),'LineStyle',Format.styles{linei})
        set(p2(linei),'color',Format.colors{linei})
        set(p2(linei),'linewidth',Format.widths{linei})
    end
title('$P_t$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')
ylabel('Percent', 'FontSize', Format.FontSizeAxes, 'FontWeight', Format.fontweight,'Interpreter','latex')

subplot(3,4,3)
p3 = plot(-window:window, BBFI.ZT.IRBMean,-window:window, BBII.ZT.IRBMean); 
    for linei = 1:2
        set(p3(linei),'LineStyle',Format.styles{linei})
        set(p3(linei),'color',Format.colors{linei})
        set(p3(linei),'linewidth',Format.widths{linei})
    end
title('$B_{t+1}$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')
ylabel('Level', 'FontSize', Format.FontSizeAxes, 'FontWeight', Format.fontweight,'Interpreter','latex')

subplot(3,4,4)
p3 = plot(-window:window, BBFI.ZT.IRBCMean ,-window:window, BBII.ZT.IRBCMean); 
    for linei = 1:2
        set(p3(linei),'LineStyle',Format.styles{linei})
        set(p3(linei),'color',Format.colors{linei})
        set(p3(linei),'linewidth',Format.widths{linei})
    end
title('Borrowing Limit','FontSize',20,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')
ylabel('Level', 'FontSize', Format.FontSizeAxes, 'FontWeight', Format.fontweight,'Interpreter','latex')

subplot(3,4,5)
p5 = plot(-window:window, (BBFI.ZN.IRCTMean - mCT)*100,-window:window, BBII.ZN.IRCTMean*100); 
    for linei = 1:2
        set(p5(linei),'LineStyle',Format.styles{linei})
        set(p5(linei),'color',Format.colors{linei})
        set(p5(linei),'linewidth',Format.widths{linei})
    end
title('$C_t^T$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')
ylabel('Percent', 'FontSize', Format.FontSizeAxes, 'FontWeight', Format.fontweight,'Interpreter','latex')

subplot(3,4,6)
p6 = plot(-window:window, (BBFI.ZN.IRPMean-FICEmP)*100,-window:window, 100*BBII.ZN.IRPMean); 
    for linei = 1:2
        set(p6(linei),'LineStyle',Format.styles{linei})
        set(p6(linei),'color',Format.colors{linei})
        set(p6(linei),'linewidth',Format.widths{linei})
    end
title('$P_t$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')
ylabel('Percent', 'FontSize', Format.FontSizeAxes, 'FontWeight', Format.fontweight,'Interpreter','latex')

subplot(3,4,7)
p7 =  plot(-window:window, BBFI.ZN.IRBMean,-window:window, BBII.ZN.IRBMean);  
    for linei = 1:2
        set(p7(linei),'LineStyle',Format.styles{linei})
        set(p7(linei),'color',Format.colors{linei})
        set(p7(linei),'linewidth',Format.widths{linei})
    end
title('$B_{t+1}$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')
ylabel('Level', 'FontSize', Format.FontSizeAxes, 'FontWeight', Format.fontweight,'Interpreter','latex')
subplot(3,4,8)
p8 = plot(-window:window, BBFI.ZN.IRBCMean ,-window:window, BBII.ZN.IRBCMean); 
    for linei = 1:2
        set(p8(linei),'LineStyle',Format.styles{linei})
        set(p8(linei),'color',Format.colors{linei})
        set(p8(linei),'linewidth',Format.widths{linei})
    end
title('Borrowing Limit','FontSize',20,'FontWeight',Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')
ylabel('Level', 'FontSize', Format.FontSizeAxes, 'FontWeight', Format.fontweight,'Interpreter','latex')
subplot(3,4,9)
p9 = plot(-window:window, (BBFI.GT.IRCTMean - mCT)*100,-window:window, 100*BBII.GT.IRCTMean); 
    for linei = 1:2
        set(p9(linei),'LineStyle',Format.styles{linei})
        set(p9(linei),'color',Format.colors{linei})
        set(p9(linei),'linewidth',Format.widths{linei})
    end
title('$C_t^T$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
ylabel('Percent', 'FontSize',Format.FontSizeAxes, 'FontWeight', Format.fontweight,'Interpreter','latex')
xlabel('Years', 'FontSize',Format.FontSizeAxes, 'FontWeight', Format.fontweight,'Interpreter','latex')
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
subplot(3,4,10)
p10 = plot(-window:window, (BBFI.GT.IRPMean-FICEmP)*100,-window:window, BBII.GT.IRPMean*100); 
    for linei = 1:2
        set(p10(linei),'LineStyle',Format.styles{linei})
        set(p10(linei),'color',Format.colors{linei})
        set(p10(linei),'linewidth',Format.widths{linei})
    end
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
title('$P_t$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')
ylabel('Percent', 'FontSize',Format.FontSizeAxes, 'FontWeight', Format.fontweight,'Interpreter','latex')
subplot(3,4,11)
p11 = plot(-window:window, BBFI.GT.IRBMean,-window:window, BBII.GT.IRBMean);  
    for linei = 1:2
        set(p11(linei),'LineStyle',Format.styles{linei})
        set(p11(linei),'color',Format.colors{linei})
        set(p11(linei),'linewidth',Format.widths{linei})
    end
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
title('$B_{t+1}$','FontSize',Format.FontSize,'FontWeight',Format.fontweight,'Interpreter','latex')
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')
ylabel('Level', 'Interpreter', 'latex');
subplot(3,4,12)
p12 = plot(-window:window, BBFI.GT.IRBCMean ,-window:window, BBII.GT.IRBCMean); 
    for linei = 1:2
        set(p12(linei),'LineStyle',Format.styles{linei})
        set(p12(linei),'color',Format.colors{linei})
        set(p12(linei),'linewidth',Format.widths{linei})
    end
xticks([-5,-4,-3,-2,-1,0,1,2,3,4, 5 ]) ;
xlabel('Years', 'FontSize',12, 'FontWeight', Format.fontweight,'Interpreter','latex')
ylabel('Level', 'Interpreter', 'latex');
title('Borrowing Limit','FontSize',20,'FontWeight',Format.fontweight,'Interpreter','latex')
legend('Perfect Information','Imperfect Information', 'Location', 'SouthEast', 'Orientation','horizontal', 'FontSize', 20,'Interpreter','latex');
% add a bit space to the figure
fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.3113  ;
Lgnd.Position(2) = 0.02651;

% Create textbox
annotation(f1,'textbox',...
    [0.0910137363010006 0.721769764420462 0.131750305175781 0.0443690906343448],...
    'String','Shock to $Z_t^T$',...
    'Rotation',90,...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',22,...
    'EdgeColor','none');

annotation(f1,'textbox',...
    [0.0910137363010006 0.43 0.131750305175781 0.0443690906343448],...
    'String','Shock to $Z_t^N$',...
    'Rotation',90,...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',22,...
    'EdgeColor','none');

annotation(f1,'textbox',...
    [0.0910137363010006 0.12 0.131750305175781 0.0443690906343448],...
    'String','Shock to $g_t$',...
    'Rotation',90,...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',22,...
    'EdgeColor','none');

saveas(f1, 'FigureA2', 'png');



