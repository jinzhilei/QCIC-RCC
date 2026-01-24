%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script file is used to perform Sensitivity of parameters

% runs: the number of samples to perform the random simulation
% num_vars: the number of parameters for sensitivity analysis
% num_results: the number of output results
% time = the time required for sensitivity analysis comparison

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runs = 100;
num_vars = 10;
num_results = 4;
time = 50;

Parameter_settings_LHS;
theta_D0_LHS = LHS_Call(1e7, theta_D0, 1e9, 0 ,runs,'unif');
theta_TN4_LHS = LHS_Call(5e6, theta_TN4, 8e9, 0 ,runs,'unif');
theta_TN8_LHS = LHS_Call(5e6, theta_TN8, 8e9, 0 ,runs,'unif');
theta_M_LHS = LHS_Call(1e7, theta_M, 8e8, 0 ,runs,'unif');
kappa_Th_LHS = LHS_Call(0, kappa_Th, 1, 0 ,runs,'unif');
kappa_Tr_LHS = LHS_Call(0, kappa_Tr, 1, 0 ,runs,'unif');
kappa_Tc_LHS = LHS_Call(0, kappa_Tc, 1, 0 ,runs,'unif');
beta_Th_LHS = LHS_Call(0, beta_Th, 0.45, 0 ,runs,'unif');
beta_Tr_LHS = LHS_Call(0, beta_Tr, 0.6, 0 ,runs,'unif');
beta_Tc_LHS = LHS_Call(0, beta_Tc, 1, 0 ,runs,'unif');

X=[theta_D0_LHS,theta_TN4_LHS,theta_TN8_LHS,theta_M_LHS,kappa_Th_LHS,...
    kappa_Tr_LHS,kappa_Tc_LHS,beta_Th_LHS,beta_Tr_LHS,beta_Tc_LHS];
D = zeros(runs,num_results);

for i=1:runs
    [t,y]=ode15s(@(t,x) ODE_LHS(t,x,X(i,:)),tspan,y0);
    C=[t y]; % [time y]

    % Save the outputs at ALL time points
    D(i,1)=C(time+1,13); % Th_D
    D(i,2)=C(time+1,16); % Tr_D
    D(i,3)=C(time+1,22); % Tc_D
    D(i,4)=C(time+1,25); % TAM
end
% Save the workspace
save ./data_prcc/Model_LHS.mat;

for var_idx = 1:num_results
    y = D(:,var_idx);
    mdl = fitlm(X, y);
    coefficients = mdl.Coefficients.Estimate(2:end);
    std_X = std(X);
    std_y = std(y);
    standardized_coefficients(:, var_idx) = coefficients .* (std_X' / std_y);
end

file=strcat('./data_prcc/model_prcc.mat');
save(file,'standardized_coefficients');
