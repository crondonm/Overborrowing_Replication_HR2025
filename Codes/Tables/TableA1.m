%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% This code replicates Table A1 in the Online Appendix using a pre-existing database.
%
% To generate this database, you need to solve the baseline model and the sensitivity scenarios. Code that replicates the descriptive statistics for the baseline model can be found at the end of this file.
% 
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

load('../Replication/Data/Param.mat')
load('../Replication/Data/Descriptive_stats_appendix.mat')

% Deep parameters

r = Param.r;
eta = Param.eta;
rho = Param.rho;
beta = Param.beta;
kappa = Param.kappa;
omega = Param.omega;
burn = Param.burn;
g = Param.g;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sensitivity Analysis
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'TableA1.tex';
% Open the file for writing
fileID = fopen(filename, 'w');
% Print the table rows using fprintf
fprintf(fileID, '\\begin{landscape}\n');
fprintf(fileID, '\\begin{table}[]\n');
fprintf(fileID, '\\centering\n');
fprintf(fileID, '\\begin{adjustbox}{max width=1.5\\textwidth} \n');
fprintf(fileID, '\\begin{threeparttable} \n');
fprintf(fileID, '\\caption{Sensitivity Analysis}\n');
fprintf(fileID, '\\label{tab:senstable}\n');
fprintf(fileID, '\\begin{tabular}{@{}lccccccccccccccccccccccccccccc@{}}\n');
fprintf(fileID, '\\toprule \n');
fprintf(fileID, '\\toprule \n');
fprintf(fileID, ' &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  & \\multicolumn{12}{c}{\\textit{Severity of Financial Crises}} \\\\ \\cmidrule(l){19-30} \n');
fprintf(fileID, ' &  &  &  &  &  &  &  & \\multicolumn{4}{c}{\\textit{\\begin{tabular}[c]{@{}c@{}}Debt-to-Output\\\\ Ratio\\end{tabular}}} &  & \\multicolumn{4}{c}{\\textit{Probability of Crises}} &  & \\multicolumn{4}{c}{\\textit{Consumption}} & \\multicolumn{4}{c}{\\textit{RER}} & \\multicolumn{4}{c}{\\textit{Current Account}} \\\\ \\cmidrule(lr){9-12} \\cmidrule(lr){14-17} \\cmidrule(l){19-30} \n');
fprintf(fileID, ' &  & \\multicolumn{2}{c}{\\textit{\\begin{tabular}[c]{@{}c@{}}Welfare\\\\ Costs\\end{tabular}}} &  & \\multicolumn{2}{c}{\\textit{\\begin{tabular}[c]{@{}c@{}}Tax\\\\ on Debt\\end{tabular}}} &  & \\multicolumn{2}{c}{\\textit{Perf. Info}} & \\multicolumn{2}{c}{\\textit{Imp. Info}} &  & \\multicolumn{2}{c}{\\textit{Perf. Info}} & \\multicolumn{2}{c}{\\textit{Imp. Info}} &  & \\multicolumn{2}{c}{\\textit{Perf. Info}} & \\multicolumn{2}{c}{\\textit{Imp. Info}} & \\multicolumn{2}{c}{\\textit{Perf. Info}} & \\multicolumn{2}{c}{\\textit{Imp. Info}} & \\multicolumn{2}{c}{\\textit{Perf. Info}} & \\multicolumn{2}{c}{\\textit{Imp. Info}} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){6-7} \\cmidrule(lr){9-12} \\cmidrule(lr){14-17} \\cmidrule(l){19-30} \n');
fprintf(fileID, ' &  & \\textit{Perf. Info} & \\textit{Imp. Info} &  & \\textit{Perf. Info} & \\textit{Imp. Info} &  & \\textit{D.E} & \\textit{S.P} & \\textit{D.E} & \\textit{S.P} &  & \\textit{D.E} & \\textit{S.P} & \\textit{D.E} & \\textit{S.P} &  & \\textit{D.E} & \\textit{S.P} & \\textit{D.E} & \\textit{S.P} & \\textit{D.E} & \\textit{S.P} & \\textit{D.E} & \\textit{S.P} & \\textit{D.E} & \\textit{S.P} & \\textit{D.E} & \\textit{S.P} \\\\ \\midrule\n');
% Fill the column values with placeholder %8.2f starting from "Baseline"
fprintf(fileID, '\\textit{Baseline $(\\beta = 0.83,\\ \\kappa = 0.335)$}  &  & %8.2f & %8.2f &  & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\\n', welf_base_fice, welf_base_iice, tau_base_fip, tau_base_iip, dtoy_base_fice, dtoy_base_fip, dtoy_base_iice, dtoy_base_iip, freq_base_fice, freq_base_fip, freq_base_iice, freq_base_iip, c_base_fice, c_base_fip, c_base_iice, c_base_iip, rer_base_fice, rer_base_fip, rer_base_iice, rer_base_iip, ca_base_fice, ca_base_fip, ca_base_iice, ca_base_iip);
fprintf(fileID, '\\textit{$\\beta = 0.90$}     &  & %8.2f & %8.2f &  & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', welf_b90_fice, welf_b90_iice, tau_b90_fip, tau_b90_iip, dtoy_b90_fice, dtoy_b90_fip, dtoy_b90_iice, dtoy_b90_iip, freq_b90_fice, freq_b90_fip, freq_b90_iice, freq_b90_iip, c_b90_fice, c_b90_fip, c_b90_iice, c_b90_iip, rer_b90_fice, rer_b90_fip, rer_b90_iice, rer_b90_iip, ca_b90_fice, ca_b90_fip, ca_b90_iice, ca_b90_iip);
fprintf(fileID, '\\textit{$\\epsilon = 0.50$} &  & %8.2f & %8.2f &  & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', welf_s05_fice, welf_s05_iice, tau_s05_fip, tau_s05_iip, dtoy_s05_fice, dtoy_s05_fip, dtoy_s05_iice, dtoy_s05_iip, freq_s05_fice, freq_s05_fip, freq_s05_iice, freq_s05_iip, c_s05_fice, c_s05_fip, c_s05_iice, c_s05_iip, rer_s05_fice, rer_s05_fip, rer_s05_iice, rer_s05_iip, ca_s05_fice, ca_s05_fip, ca_s05_iice, ca_s05_iip);
fprintf(fileID, '\\textit{Recalibrated F.I Economy $(\\beta = 0.53, \\ \\kappa = 0.3525)$} &  & %8.2f & %8.2f &  & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \\midrule \n', welf_reca_fice, welf_base_iice, tau_reca_fip, tau_base_iip, dtoy_reca_fice, dtoy_reca_fip, dtoy_base_iice, dtoy_base_iip, freq_reca_fice, freq_reca_fip, freq_base_iice, freq_base_iip, c_reca_fice, c_reca_fip, c_base_iice, c_base_iip, rer_reca_fice, rer_reca_fip, rer_base_iice, rer_base_iip, ca_reca_fice, ca_reca_fip, ca_base_iice, ca_base_iip);
fprintf(fileID, '\\textit{Autocorrelation $\\rho_g$ (15 \\%% less)}   &  & %8.2f & %8.2f &  & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f\\\\\n', welf_rlow_fice, welf_rlow_iice, tau_rlow_fip, tau_rlow_iip, dtoy_rlow_fice, dtoy_rlow_fip, dtoy_rlow_iice, dtoy_rlow_iip, freq_rlow_fice, freq_rlow_fip, freq_rlow_iice, freq_rlow_iip, c_rlow_fice, c_rlow_fip, c_rlow_iice, c_rlow_iip, rer_rlow_fice, rer_rlow_fip, rer_rlow_iice, rer_rlow_iip, ca_rlow_fice, ca_rlow_fip, ca_rlow_iice, ca_rlow_iip);
fprintf(fileID, '\\textit{Autocorrelation $\\rho_g$ (15 \\%% more)} &  & %8.2f & %8.2f &  & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f\\\\\n', welf_rhigh_fice, welf_rhigh_iice, tau_rhigh_fip, tau_rhigh_iip, dtoy_rhigh_fice, dtoy_rhigh_fip, dtoy_rhigh_iice, dtoy_rhigh_iip, freq_rhigh_fice, freq_rhigh_fip, freq_rhigh_iice, freq_rhigh_iip, c_rhigh_fice, c_rhigh_fip, c_rhigh_iice, c_rhigh_iip, rer_rhigh_fice, rer_rhigh_fip, rer_rhigh_iice, rer_rhigh_iip, ca_rhigh_fice, ca_rhigh_fip, ca_rhigh_iice, ca_rhigh_iip);
fprintf(fileID, '\\textit{Volatility $\\sigma_g$ (15 \\%% less)}         &  & %8.2f & %8.2f &  & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f\\\\\n', welf_slow_fice, welf_slow_iice, tau_slow_fip, tau_slow_iip, dtoy_slow_fice, dtoy_slow_fip, dtoy_slow_iice, dtoy_slow_iip, freq_slow_fice, freq_slow_fip, freq_slow_iice, freq_slow_iip, c_slow_fice, c_slow_fip, c_slow_iice, c_slow_iip, rer_slow_fice, rer_slow_fip, rer_slow_iice, rer_slow_iip, ca_slow_fice, ca_slow_fip, ca_slow_iice, ca_slow_iip);
fprintf(fileID, '\\textit{Volatility $\\sigma_g$ (15 \\%% more)}       &  & %8.2f & %8.2f &  & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f &  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f\\\\\n', welf_shigh_fice, welf_shigh_iice, tau_shigh_fip, tau_shigh_iip, dtoy_shigh_fice, dtoy_shigh_fip, dtoy_shigh_iice, dtoy_shigh_iip, freq_shigh_fice, freq_shigh_fip, freq_shigh_iice, freq_shigh_iip, c_shigh_fice, c_shigh_fip, c_shigh_iice, c_shigh_iip, rer_shigh_fice, rer_shigh_fip, rer_shigh_iice, rer_shigh_iip, ca_shigh_fice, ca_shigh_fip, ca_shigh_iice, ca_shigh_iip);

% Close the environments and the file
fprintf(fileID, '\\bottomrule\n');
fprintf(fileID, '\\bottomrule\n');
fprintf(fileID, '\\end{tabular}%%\n');
fprintf(fileID, '\\begin{tablenotes} \n');
fprintf(fileID, '\\footnotesize \\baselineskip=14pt \n');
fprintf(fileID, '\\item Note: This table summarizes the key descriptive statistics for each alternative calibration tested. We measure the welfare costs of the pecuniary externality as the percentage of lifetime consumption lost relative to a Social Planner sharing the same information set. The probability of a financial crisis is defined as the frequency of crises occurring every 100 years. The consumption and RER changes during a financial crisis are the percentage changes relative to their respective ergodic means when the collateral constraint binds and the economy experiences a one-standard deviation current account reversal. \n');
fprintf(fileID, '\\end{tablenotes} \n');
fprintf(fileID, '\\end{threeparttable} \n');
fprintf(fileID, '\\end{adjustbox}\n');
fprintf(fileID, '\\end{table}\n');
fprintf(fileID, '\\end{landscape}\n');

fclose(fileID);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replication
% Descriptive statistics for the baseline model (as shown in Table A1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

replicate = 0;

if replicate == 1
    
    clearvars 
    close all

    load('../Replication/Data/Param.mat')
    load('../Replication/Data/FICE.mat')
    load('../Replication/Data/FIP.mat')
    load('../Replication/Data/IIPCC.mat')
    load('../Replication/Data/IICECC.mat')

    r = Param.r;
    eta = Param.eta;
    rho = Param.rho;
    beta = Param.beta;
    kappa = Param.kappa;
    omega = Param.omega;
    burn = Param.burn;
    g = Param.g;

   % IICECC
    rer_base_iice = 100*max(((omega^eta + (1 - omega)^eta*(IICECC.Cycles.IRPMean+IICECC.mP).^(1 - eta)).^(-1/(1-eta)))./((omega^eta + (1 - omega)^eta*(IICECC.mP).^(1 - eta)).^(-1/(1-eta)))-1);
    c_base_iice = min((IICECC.Cycles.IRCMean/mean(IICECC.CSim)-1)*100);
    cT_base_iice=min(IICECC.Cycles.IRCTMean)*100;
    ca_base_iice = max((IICECC.Cycles.IRCAtoYMean + IICECC.mCAtoY))*100;
    welf_base_iice = IICECC.mWelfcost_ii;
    freq_base_iice = IICECC.Freq*100;
    dtoy_base_iice = IICECC.mDtoY*100;
    
    % IIPCC
    rer_base_iip = 100*max(((omega^eta + (1 - omega)^eta*(IIPCC.Cycles.IRPMean).^(1 - eta)).^(-1/(1-eta)))./((omega^eta + (1 - omega)^eta*(IIPCC.mP).^(1 - eta)).^(-1/(1-eta)))-1);
    c_base_iip = min((IIPCC.Cycles.IRCMean/mean(IIPCC.CSim)-1)*100);
    cT_base_iip = min(IIPCC.Cycles.IRCTMean./IIPCC.mCT-1)*100;
    ca_base_iip = max((IIPCC.Cycles.IRCAtoYMean + IIPCC.mCAtoY))*100;
    freq_base_iip = IIPCC.Freq*100;
    dtoy_base_iip = IIPCC.mDtoY*100;
    tau_base_iip = mean(IIPCC.TAOSim)*100;
    
    % FICE
    rer_base_fice = 100*max(((omega^eta + (1 - omega)^eta*(FICE.Cycles.IRPMean).^(1 - eta)).^(-1/(1-eta)))./((omega^eta + (1 - omega)^eta*(FICE.mP).^(1 - eta)).^(-1/(1-eta)))-1);
    c_base_fice = min((FICE.Cycles.IRCMean/mean(FICE.CSim)-1)*100);
    cT_base_fice = min((FICE.Cycles.IRCTMean/mean(FICE.CTSim)-1)*100);
    ca_base_fice = max((FICE.Cycles.IRCAtoYMean))*100;
    welf_base_fice = FICE.mWelfcost_fi;
    freq_base_fice = FICE.Freq*100;
    dtoy_base_fice = FICE.mDtoY*100;
    
    % FIP
    tau_base_fip = mean(FIP.TAOSim)*100;
    freq_base_fip = FIP.Freq*100;
    dtoy_base_fip = FIP.mDtoY*100;
    rer_base_fip = 100*max(((omega^eta + (1 - omega)^eta*(FIP.Cycles.IRPMean).^(1 - eta)).^(-1/(1-eta)))./((omega^eta + (1 - omega)^eta*(FIP.mP).^(1 - eta)).^(-1/(1-eta)))-1);
    c_base_fip = min((FIP.Cycles.IRCMean/mean(FIP.CSim)-1)*100);
    ca_base_fip = max((FIP.Cycles.IRCAtoYMean))*100;
    clear p FICE FIP IICECC IIPCC 
    
end

