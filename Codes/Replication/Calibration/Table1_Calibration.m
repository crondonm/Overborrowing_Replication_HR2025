%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% IMPORTANT: Make sure you have downloaded the neccesary files from the repository. See README.m for more details
%
% In this code: 
%               1. Produce table 1 in the paper.
%               2. Generates point estimates used for stochastic processes
%                   in the model.
%          
% Authors:  Juan Herreño, jherrenolopera@ucsd.edu
%               Carlos Rondón Moreno, crondon@bcentral.cl
% Date:      March 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
warning('off', 'all');
clc

data = load("../Data/Calibration.mat");
load("../Data/Calibration.mat")
disp("************************************************")
disp("********TABLE 1*************")
disp("************************************************")

names = ["rho_zt_zt", "rho_zt_zn", "rho_zn_zt", "rho_zn_zn", "rho_g", "sigma_zt", "covar_zt_zn", "sigma_zn", "sigma_g"  ]

param = [xopt(1), xopt(2), xopt(3), xopt(4), xopt(5), sqrt(xopt(6)), xopt(7), sqrt(xopt(8)), sqrt(xopt(9))  ]
std_dev = sqrt(abs(diag(inv(hessiancsd(opt_LL, xopt)))))

T = table( names', param', std_dev, 'VariableNames', {'Parameter', 'Values', 'Std. Dev'});
disp(T)

%% Replication

replicate = 1;

if replicate == 1

    % Import Data

    % lnsignalYT and lnsignalYN contains the output data constructed as
    % detailed in the online appendix

    neg_LL = @(x) LL_klm(x, [lnsignalYT';lnsignalYN']);
    opt_LL = @(x) LL_klm_opt(x, [lnsignalYT';lnsignalYN']);
    
    lb   = [ -0.99; -0.99; -0.99;-0.99;  -0.99;  1e-6; -cov_yTN(1,2); 1e-6; 1e-6]';
    ub  = [ 0.99;  0.99;  0.99;  0.99;  0.99;  std_yt; cov_yTN(1,2); std_yn;  std_y]';
    
    %% Pattern Search + Simulated Annealing
    
    rng(10,'twister') % for reproducibility
    
    % x0 is computed as detailes in the online appendix

    options = optimoptions('patternsearch','Display','iter','PlotFcn',{@psplotbestf, @psplotbestx}, 'MaxIterations', 6e4, 'MaxFunctionEvaluations', 10e10);
    [sol1, fval2] = patternsearch(neg_LL, x0, [],[],[],[],lb,ub, options);
    options = optimoptions('fminunc','Display','iter','MaxIterations', 6e4, 'MaxFunctionEvaluations', 10e10, 'HessianApproximation', 'bfgs', 'PlotFcn',{@optimplotfval,@optimplotx,@optimplotfirstorderopt});
    [xopt, fval3,~,~,~,~] = fminunc(neg_LL, sol1, options);
    
    names = ["rho_zt_zt", "rho_zt_zn", "rho_zn_zt", "rho_zn_zn", "rho_g", "sigma_zt", "covar_zt_zn", "sigma_zn", "sigma_g"  ]
    param = [xopt(1), xopt(2), xopt(3), xopt(4), xopt(5), sqrt(xopt(6)), xopt(7), sqrt(xopt(8)), sqrt(xopt(9))  ]
    std_dev = sqrt(abs(diag(inv(hessiancsd(opt_LL,xopt)))))

    
    T = table( names', param', std_dev, 'VariableNames', {'Parameter', 'Values', 'Std. Dev'});
    disp(T)
     
end

