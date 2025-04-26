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

load('../Replication/Data/IIPCC.mat')
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


%% Figure 13: Optimal Macroprudential Policy: Relative to the Recalibrated Economy

[datay1, datax1] = ksdensity(FIPrec.TAOSim(FIPrec.TAOSim<10)*100, 'Support','nonnegative', 'BoundaryCorrection', 'reflection');
[datay2, datax2] = ksdensity(IIPCC.TAOSim   *100, 'Support','nonnegative', 'BoundaryCorrection', 'reflection');    

Format.styles = {'-.', '-', '-.','+'};

f13 = figure('Color',[1 1 1], 'Position',Format.figsize);
p1 = plot(datax1, datay1, datax2, datay2) ;
for linei = 1:2
    set(p1(linei),'LineStyle',Format.styles{linei})
    set(p1(linei),'color',Format.colors{linei})
    set(p1(linei),'linewidth',Format.widths{linei})
end
%title('Ergodic Distribution of Asset Holdings: Incomplete Information','FontSize',Format.FontSize,'FontWeight',Format.fontweight)
ylabel('Frequency','FontSize',Format.FontSize, 'Interpreter','Latex')
xlabel('Optimal Tax Rate (\%) ','FontSize',Format.FontSize, 'Interpreter','Latex')
legend('Perfect Information: Recalibrated Economy','Imperfect Information','Location','NorthEast')
xlim([0*100 1*100 ]) ;
set(gca,'FontSize',Format.FontSizeAxes)
set(f13,'PaperPositionMode', 'auto')
saveas(f13,'Figure13_ErgDist_TaxesIIvsRecalibrated','png');



