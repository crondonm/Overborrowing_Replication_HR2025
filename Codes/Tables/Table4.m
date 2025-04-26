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

load('../Replication/Data/Cyclicality_taxes.mat')


%% Baseline Model

names = ["Total Output", "Permanent Component", "Transitory Tradable", "Transitory Non-Tradable"];
perf = [corr_fip_tau_ytot, corr_fip_tau_g, corr_fip_tau_zt, corr_fip_tau_zn ];
imp = [corr_iip_tau_ytot, corr_iip_tau_g, corr_iip_tau_zt, corr_iip_tau_zn ];

perf = round(perf, 2);
imp = round(imp, 2);

fprintf('+------------+------------+-------------------------+\n');
fprintf('+-------------+BASELINE MODEL+----------------------+\n');
fprintf('+------------+------------+-------------------------+\n');

T = table( names', perf', imp', 'VariableNames', {'Variable', 'Perfect Info', 'Imp. Info'});
disp(T)

%% Elasticity of Substitution

perf = [corr_s05_fip_tau_ytot, corr_s05_fip_tau_g, corr_s05_fip_tau_zt, corr_s05_fip_tau_zn ];
imp = [corr_s05_iip_tau_ytot, corr_s05_iip_tau_g, corr_s05_iip_tau_zt, corr_s05_iip_tau_zn ];

perf = round(perf, 2);
imp = round(imp, 2);
fprintf('+------------+------------+-------------------------------------------+\n');
fprintf('+-------------+ELASTICITY OF SUBSTITUTION epsilon = 0.5+-------------------+\n');
fprintf('+------------+------------+------------------------------------------+\n');

T = table( names', perf', imp', 'VariableNames', {'Variable', 'Perfect Info', 'Imp. Info'});
disp(T)


%% Rho_g LOW

fprintf('+------------+------------+-------------------------+\n');
fprintf('+-------------rho_g LOW-------------------+\n');
fprintf('+------------+------------+-------------------------+\n');



perf = [corr_rlow_fip_tau_ytot, corr_rlow_fip_tau_g, corr_rlow_fip_tau_zt, corr_rlow_fip_tau_zn ];
imp = [corr_rlow_iip_tau_ytot, corr_rlow_iip_tau_g, corr_rlow_iip_tau_zt, corr_rlow_iip_tau_zn ];
perf = round(perf, 2);
imp = round(imp, 2);

T = table( names', perf', imp', 'VariableNames', {'Variable', 'Perfect Info', 'Imp. Info'});
disp(T)

%% Rho_g HIGH

fprintf('+------------+------------+-------------------------+\n');
fprintf('+-------------rho_g HIGH-------------------+\n');
fprintf('+------------+------------+-------------------------+\n');

perf = [corr_rhigh_fip_tau_ytot, corr_rhigh_fip_tau_g, corr_rhigh_fip_tau_zt, corr_rhigh_fip_tau_zn ];
imp = [corr_rhigh_iip_tau_ytot, corr_rhigh_iip_tau_g, corr_rhigh_iip_tau_zt, corr_rhigh_iip_tau_zn ];

perf = round(perf, 2);
imp = round(imp, 2);

T = table( names', perf', imp', 'VariableNames', {'Variable', 'Perfect Info', 'Imp. Info'});
disp(T)

%% Sigma_g LOW

fprintf('+------------+------------+-------------------------+\n');
fprintf('+-------------sigma_g LOW-------------------+\n');
fprintf('+------------+------------+-------------------------+\n');

perf = [corr_slow_fip_tau_ytot, corr_slow_fip_tau_g, corr_slow_fip_tau_zt, corr_slow_fip_tau_zn ];
imp = [corr_slow_iip_tau_ytot, corr_slow_iip_tau_g, corr_slow_iip_tau_zt, corr_slow_iip_tau_zn ];

perf = round(perf, 2);
imp = round(imp, 2);

T = table( names', perf', imp', 'VariableNames', {'Variable', 'Perfect Info', 'Imp. Info'});
disp(T)

%% Sigma_g HIGH

fprintf('+------------+------------+-------------------------+\n');
fprintf('+-------------sigma_g HIGH-------------------+\n');
fprintf('+------------+------------+-------------------------+\n');

perf = [corr_shigh_fip_tau_ytot, corr_shigh_fip_tau_g, corr_shigh_fip_tau_zt, corr_shigh_fip_tau_zn ];
imp = [corr_shigh_iip_tau_ytot, corr_shigh_iip_tau_g, corr_shigh_iip_tau_zt, corr_shigh_iip_tau_zn ];

perf = round(perf, 2);
imp = round(imp, 2);

T = table( names', perf', imp', 'VariableNames', {'Variable', 'Perfect Info', 'Imp. Info'});
disp(T)

