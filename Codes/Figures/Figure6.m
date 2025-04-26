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

%% Figure 6: Welfare Distribution
% Welfare Costs of the Pecunary Externality Under Different Information Sets

Format.styles = {'-',':','-.','-.'};
[datay1, datax1] = ksdensity(FICE.Welfcost_fi, 'Support','unbounded', 'BoundaryCorrection','reflection');
[datay2, datax2] = ksdensity(IICECC.Welfcost_ii, 'Support','unbounded', 'BoundaryCorrection','reflection');

f6 = figure('Color',[1 1 1]);
p1 = plot(datax1,datay1,datax2,datay2) ;
for linei = 1:2
    set(p1(linei),'LineStyle',Format.styles{linei})
    set(p1(linei),'color',Format.colors{linei})
    set(p1(linei),'linewidth',Format.widths{linei})
end
xlim([0 .5 ]) ;
%ylim([0 1 ]) ;
%title('Ergodic Distribution of Welfare Costs','FontSize',Format.FontSize,'FontWeight',Format.fontweight)
ylabel('Probability','FontSize',Format.FontSize, 'Interpreter','Latex')
xlabel('Welfare Gain (\%)','FontSize',Format.FontSize, 'Interpreter','Latex')
legend('Full Information','Imperfect Information','Location','NorthEast','FontSize',14)
set(gca,'FontSize',Format.FontSizeAxes)
%saveas(f6,'Figure6_Dist_Welfare','png');
filename = 'Figure6.png';
resolution = 300; % DPI
% Export the graphics
exportgraphics(gcf, filename, 'Resolution', resolution);

