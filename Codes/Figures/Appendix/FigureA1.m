%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% In this code: Plot figure A1
% 
% Authors:  Juan Herreño, jherrenolopera@ucsd.edu
%               Carlos Rondón Moreno, crondon@bcentral.cl
%
% Last:  March 2025
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping

clear all
clc
warning off

% Set path

cd('../../../Calibration')


% Import Data

opts = spreadsheetImportOptions("NumVariables", 9);
% Specify sheet and rangedelta
opts.Sheet = "Data";
opts.DataRange = "A30:I145";
% Specify column names and types
opts.VariableNames = ["Year", "Agriculture", "Manufacturing", "Services", "Tradable", "SignalYT", "SignalYN", "Total", "Ratio"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double"];
filename = "Data.xlsx";

% Import the dataco
Data = readtable(filename, opts, "UseExcel", false);

lnsignalYT = (Data.SignalYT) - mean(Data.SignalYT);
lnsignalYN = (Data.SignalYN) - mean(Data.SignalYN);
T = size(lnsignalYN,1);


% Fit a quadratic polynomial (order = 2) to the data
order = 2; 
t = 1:T;
%p = polyfit(t, lnsignalYT, order);
%quadratic_trend = polyval(p, t');    % Evaluate the fitted quadratic polynomial at the time points
%lnsignalYT = lnsignalYT - quadratic_trend; % Remove the quadratic trend from the time series
p = polyfit(t, lnsignalYN, order);
quadratic_trend = polyval(p, t');
lnsignalYN = lnsignalYN - quadratic_trend; % Remove the quadratic trend from the time series

lntotalY = (Data.Ratio);
num=10;
std_yt = num*std(lnsignalYT)^2;
std_yn = num*std(lnsignalYN)^2;
std_y = num*std(lntotalY)^2;
cov_yTN = num*cov(lnsignalYT', lnsignalYN');
grY = 100*1.01/100;


x_0   = [ 0.734679129925539; -0.255320869154147; 0.0336525070039930; 0.417047371095781; 0.496779121757470;  0.00461817366335938; 0.000379607811872954; 0.00137117483970497; 0.00327510105056510]';

lb = 0.01 ;
ub = 1.8 ;
ntry = 100 ;

array = linspace(lb,ub,ntry) ;

for jj = 1:ntry
    for ll = 1:ntry
        x_0_temp   = x_0 ;
        x_0_temp(5) = min(x_0(5)*array(jj),0.99) ;
        x_0_temp(end) = x_0(end)*array(ll) ;
        loglik(ll,jj) = LL_klm(x_0_temp, [lnsignalYT';lnsignalYN']) ;
    end
end


fig = figure; 
contour(array, array, loglik, 1000);
filename = '../Figures/Appendix/FigureA1.png';
resolution = 300; % DPI

% Export the graphics
exportgraphics(gcf, filename, 'Resolution', resolution);


