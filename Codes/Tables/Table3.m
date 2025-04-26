%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% In this code: We create all the tables included in the paper                   
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

% Base Model
load('../Replication/Data/Param.mat')
load('../Replication/Data/FICE.mat')
load('../Replication/Data/FIP.mat')
load('../Replication/Data/IIPCC.mat')
load('../Replication/Data/IICECC.mat')

% Recalibrated Model
fice_rec = load('../Replication/Data/FICE_rec.mat');
fip_rec = load('../Replication/Data/FIP_rec.mat');
fice_rec = fice_rec.FICE;
fip_rec = fip_rec.FIP;

% Deep parameters

r = Param.r;
eta = Param.eta;
rho = Param.rho;
beta = Param.beta;
kappa = Param.kappa;
omega = Param.omega;
burn = Param.burn;
g = Param.g;

 %% Tables


% Debt-to-GDP Ratios

fice_mdtoy = FICE.mDtoY*100;
iice_mdtoy = IICECC.mDtoY*100;
fip_mdtoy = FIP.mDtoY*100;
iip_mdtoy = IIPCC.mDtoY*100;

% Frequency of crises
fice_freq = FICE.Freq*100;
iice_freq = IICECC.Freq*100;
fip_freq = FIP.Freq*100;
iip_freq = IIPCC.Freq*100;

% Mean Taxes
fip_tao = mean(FIP.TAOSim)*100;
iip_tao = mean(IIPCC.TAOSim)*100;

% Mean welfare costs
fip_welfare = FICE.mWelfcost_fi;
iip_welfare = IICECC.mWelfcost_ii;


%% Long-Run Moments:  IICECC

Simyt = IICECC.Sim(burn:end, 1);
Simyn = IICECC.Sim(burn:end, 2);
Simg = IICECC.Sim(burn:end, 3);
PSim = IICECC.PSim+IICECC.mP;

IICECC.Ytot = ((exp(Simyt(3:end)'+Simg(3:end)'+g) + PSim.*exp(Simyn(3:end)'+Simg(3:end)'+g)))./exp(Simyt(2:end- 1)');
p = prctile(IICECC.CAtoY, [0.1, 99.9]);
IICECC_CAtoY = IICECC.CAtoY(IICECC.CAtoY >= p(1) & IICECC.CAtoY <= p(2));
std_iice_catoy = std(IICECC_CAtoY)*100;

p = prctile(IICECC.CA, [0.1, 99.9]);
IICECC_CA = IICECC.CA(IICECC.CA >= p(1) & IICECC.CA <= p(2));
temp = find((IICECC.CA>=p(1) & IICECC.CA<=p(2)));
Ytot_IICECCca = IICECC.Ytot(temp);
corr_iice_cay = corr(IICECC_CA', Ytot_IICECCca');

p = prctile(IICECC.CtoY, [0.1, 99.9]);
IICE_CSim = IICECC.CtoY(IICECC.CtoY >= p(1) & IICECC.CtoY <= p(2));
std_iice_csim = std(IICE_CSim)*100;



%% Long-Run Moments:  IIPCC

Simyt = IIPCC.Sim(burn:end, 1);
Simyn = IIPCC.Sim(burn:end, 2);
Simg = IIPCC.Sim(burn:end, 3);
PSim = IIPCC.PSim;

IIPCC.Ytot = ((exp(Simyt(3:end)'+Simg(3:end)'+g) + PSim.*exp(Simyn(3:end)'+Simg(3:end)'+g)))./exp(Simyt(2:end- 1)');
p = prctile(IIPCC.CAtoY, [0.1, 99.9]);
IIPCC_CAtoY = IIPCC.CAtoY(IIPCC.CAtoY >= p(1) & IIPCC.CAtoY <= p(2));
std_iip_catoy = std(IIPCC_CAtoY);

p = prctile(IIPCC.CA, [0.1, 99.9]);
IIPCC_CA = IIPCC.CA(IIPCC.CA >= p(1) & IIPCC.CA <= p(2));
temp = find((IIPCC.CA>=p(1) & IIPCC.CA<=p(2)));
Ytot_IIPCCca = IIPCC.Ytot(temp);
corr_iip_cay = corr(IIPCC_CA', Ytot_IIPCCca');


p = prctile(IIPCC.CtoY, [0.1, 99.9]);
IIP_CSim = IIPCC.CtoY(IIPCC.CtoY >= p(1) & IIPCC.CtoY <= p(2));
std_iip_csim = std(IIP_CSim);

%% Long-Run Moments:  FICE

Simyt = FICE.Sim(burn:end, 1);
Simyn = FICE.Sim(burn:end, 2);
Simg = FICE.Sim(burn:end, 3);
PSim = FICE.PSim;

FICE.Ytot = ((exp(Simyt(2:end)'+Simg(2:end)'+g) + PSim.*exp(Simyn(2:end)'+Simg(2:end)'+g)))./exp(Simyt(1:end- 1)');
p = prctile(FICE.CAtoY, [0.1, 99.9]);
FICE_CAtoY = FICE.CAtoY(FICE.CAtoY >= p(1) & FICE.CAtoY <= p(2));
std_fice_catoy = std(FICE_CAtoY);

p = prctile(FICE.CA, [0.1, 99.9]);
FICE_CA = FICE.CA(FICE.CA >= p(1) & FICE.CA <= p(2));
temp = find((FICE.CA>=p(1) & FICE.CA<=p(2)));
Ytot_FICEca = FICE.Ytot(temp);
corr_fice_cay = corr(FICE_CA', Ytot_FICEca');

p = prctile(FICE.CtoY, [0.1, 99.9]);
FICE_CSim = FICE.CtoY(FICE.CtoY >= p(1) & FICE.CtoY <= p(2));
std_fice_csim = std(FICE_CSim);

%% Long-Run Moments:  FIP

Simyt = FIP.Sim(burn:end, 1);
Simyn = FIP.Sim(burn:end, 2);
Simg = FIP.Sim(burn:end, 3);
PSim = FIP.PSim;

FIP.Ytot = ((exp(Simyt(2:end)'+Simg(2:end)'+g) + PSim.*exp(Simyn(2:end)'+Simg(2:end)'+g)))./exp(Simyt(1:end- 1)');
p = prctile(FIP.CAtoY, [0.1, 99.9]);
FIP_CAtoY = FIP.CAtoY(FIP.CAtoY >= p(1) & FIP.CAtoY <= p(2));
std_fip_catoy = std(FIP_CAtoY);

p = prctile(FIP.CA, [0.1, 99.9]);
FIP_CA = FIP.CA(FIP.CA >= p(1) & FIP.CA <= p(2));
temp = find((FIP.CA>=p(1) & FIP.CA<=p(2)));
Ytot_fipca = FIP.Ytot(temp);
corr_fip_cay = corr(FIP_CA', Ytot_fipca');

p = prctile(FIP.CtoY, [0.1, 99.9]);
FIP_CSim = FIP.CtoY(FIP.CtoY >= p(1) & FIP.CtoY <= p(2));
std_fip_csim = std(FIP_CSim);


%% Recalibration: FIP


fiprec_mdtoy = fip_rec.mDtoY;
fiprec_freq = fip_rec.Freq*100;
c_base_fiprec = min((fip_rec.Cycles.IRCMean/mean(fip_rec.CSim)-1)*100);
fiprec_welfare = fice_rec.mWelfcost_fi;
fiprec_tao = mean(fip_rec.TAOSim);

Simyt = fip_rec.Sim(burn:end, 1);
Simyn = fip_rec.Sim(burn:end, 2);
Simg = fip_rec.Sim(burn:end, 3);
PSim = fip_rec.PSim;

fip_rec.Ytot = ((exp(Simyt(2:end)'+Simg(2:end)'+g) + PSim.*exp(Simyn(2:end)'+Simg(2:end)'+g)))./exp(Simyt(1:end- 1)');

p = prctile(fip_rec.CA, [0.1, 99.9]);
FIPrec_CA = fip_rec.CA(fip_rec.CA >= p(1) & fip_rec.CA <= p(2));
temp = find((fip_rec.CA>=p(1) & fip_rec.CA<=p(2)));
Ytot_FIPca = fip_rec.Ytot(temp);
corr_fiprec_cay = corr(FIPrec_CA', Ytot_FIPca');

p = prctile(fip_rec.CAtoY, [0.1,99.9]);
FIPrec.CAtoY = fip_rec.CAtoY(fip_rec.CAtoY >= p(1) & fip_rec.CAtoY <= p(2));
std_fiprec_catoy = std(FIPrec.CAtoY);

fip_rec.Ytot = ((exp(Simyt(2:end)'+Simg(2:end)'+g) + PSim.*exp(Simyn(2:end)'+Simg(2:end)'+g)))./exp(Simyt(1:end- 1)');
p = prctile(fip_rec.CtoY, [0.1, 99.9]);
FIPrec_CSim = fip_rec.CtoY(fip_rec.CtoY >= p(1) & fip_rec.CtoY <= p(2));
std_fiprec_csim = std(FIPrec_CSim);

%% Recalibration: FICE

ficerec_mdtoy = fice_rec.mDtoY;
ficerec_freq = fice_rec.Freq*100;
c_base_ficerec = min((fice_rec.Cycles.IRCMean/mean(fice_rec.CSim)-1)*100);

Simyt = fice_rec.Sim(burn:end, 1);
Simyn = fice_rec.Sim(burn:end, 2);
Simg = fice_rec.Sim(burn:end, 3);
PSim = fice_rec.PSim;

fice_rec.Ytot = ((exp(Simyt(2:end)'+Simg(2:end)'+g) + PSim.*exp(Simyn(2:end)'+Simg(2:end)'+g)))./exp(Simyt(1:end- 1)');
p = prctile(fice_rec.CA, [0.1, 99.9]);
FICErec_CA = fice_rec.CA(fice_rec.CA >= p(1) & fice_rec.CA <= p(2));
temp = find((fice_rec.CA>=p(1) & fice_rec.CA<=p(2)));
Ytot_FICEca = fice_rec.Ytot(temp);
corr_ficerec_cay = corr(FICErec_CA', Ytot_FICEca');

p = prctile(fice_rec.CAtoY, [0.1, 99.9]);
FICErec.CAtoY = fice_rec.CAtoY(fice_rec.CAtoY >= p(1) & fice_rec.CAtoY <= p(2));
std_ficerec_catoy = std(FICErec.CAtoY)*100

p = prctile(fice_rec.CtoY, [0.1, 99.9]);
FICErec_CSim = fice_rec.CtoY(fice_rec.CtoY >= p(1) & fice_rec.CtoY <= p(2));
std_ficerec_csim = std(FICErec_CSim)*100;

%% Table 4

c_base_iice = min((IICECC.Cycles.IRCMean/mean(IICECC.CSim)-1)*100);
c_base_iip = min((IIPCC.Cycles.IRCMean/mean(IIPCC.CSim)-1)*100);
c_base_fice = min((FICE.Cycles.IRCMean/mean(FICE.CSim)-1)*100);
c_base_fip = min((FIP.Cycles.IRCMean/mean(FIP.CSim)-1)*100);


names = ["Debt-to-GDP Ratio", "Frequency of Fin. Crises", "C_t drop Fin. Crises"  ...
               "std(C/Y)", "corr(CA, Y)", "std(CA, Y)", "Welfare Gain", "Avg tax"];
fice = [fice_mdtoy, fice_freq,c_base_fice, std_fice_csim*100, corr_fice_cay, 100*std_fice_catoy, nan,  nan];
fip = [fip_mdtoy, fip_freq, c_base_fip, std_fip_csim*100, corr_fip_cay, 100*std_fip_catoy, fip_welfare,  fip_tao];
iice = [iice_mdtoy, iice_freq, c_base_iice, std_iice_csim, corr_iice_cay, std_iice_catoy, nan,  nan];
iip = [iip_mdtoy, iip_freq, c_base_iip, std_iip_csim*100, corr_iip_cay, 100*std_iip_catoy, iip_welfare,  iip_tao];

clc;

% Build raw table
T = table( names', fice', fip', iice', iip', ...
    'VariableNames', {'Variable', 'P_I_D_E', 'P_I_SP', 'I_I_D_E', 'I_I_SP'});

% Format numeric columns to strings with 2 decimals
for col = 2:width(T)
    T.(col) = cellfun(@(x) sprintf('%.2f', x), num2cell(T{:,col}), 'UniformOutput', false);
end

% Print header
fprintf('+-------------------------+------------+------------+------------+\n');
fprintf('+--------------------      BASELINE MODEL            ---------+------------+\n');
fprintf('+-------------------------+------------+------------+------------+\n');
fprintf('| %-23s | %-10s | %-10s | %-10s | %-10s |\n', ...
    'Variable', 'P.I D.E', 'P.I SP', 'I.I D.E', 'I.I SP');
fprintf('+-------------------------+------------+------------+------------+\n');

% Print rows
for i = 1:height(T)
    fprintf('| %-23s | %-10s | %-10s | %-10s | %-10s |\n', ...
        T.Variable{i}, T.P_I_D_E{i}, T.P_I_SP{i}, T.I_I_D_E{i}, T.I_I_SP{i});
end

% Footer
fprintf('+-------------------------+------------+------------+------------+\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   Recalibrated Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ficerec = [ficerec_mdtoy*100, ficerec_freq, c_base_ficerec, std_ficerec_csim, corr_ficerec_cay, std_ficerec_catoy, nan,  nan];
fiprec = [fiprec_mdtoy*100, fiprec_freq, c_base_fiprec, std_fiprec_csim*100, corr_fiprec_cay, 100*std_fiprec_catoy, fiprec_welfare,  100*fiprec_tao];

% Build raw table
T = table( names', ficerec', fiprec', ...
    'VariableNames', {'Variable', 'P_I_D_E', 'P_I_SP'});

% Format numeric columns to strings with 2 decimals
for col = 2:width(T)
    T.(col) = cellfun(@(x) sprintf('%.2f', x), num2cell(T{:,col}), 'UniformOutput', false);
end

fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf('\n');

fprintf('+-----------------+-----------+------------+\n');
fprintf('+-----      RECALIBRATED MODEL    --+------------+\n');
fprintf('+-----------------+-----------+------------+\n');
fprintf('| %-23s | %-10s | %-10s  |\n', ...
    'Variable', 'P.I D.E', 'P.I SP');
fprintf('+-----------------+-----------+------------+\n');

% Print rows
for i = 1:height(T)
    fprintf('| %-23s | %-10s | %-10s |\n', ...
        T.Variable{i}, T.P_I_D_E{i}, T.P_I_SP{i});
end

% Footer
fprintf('+-----------------+-----------+------------+\n');

