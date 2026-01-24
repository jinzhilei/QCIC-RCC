
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This script file is used to perform parameter calibration for virtual patients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variable Name Set

sample  =  30000;                   % sample: Virtual patients number
pare    =  10;                      % pare: Model heterogeneity parameters number
t_start =  1;                       % t_start: Start time
t_end   =  250;                     % t_end: End time
var     =  28;                      % var: Variable number
cell    =  zeros(sample,var);       % cell: Virtual patient prestorage matrix

%% Parameter Sample

M = model_parameter(sample,pare);  
file_M='./VP_immune/VP_model_parameter.mat';
save(file_M,'M');

%% Parallel Computing and Simulation Calculation

delete (gcp('nocreate'))            % Close the current parallel computing pool
parpool(7);                         % Create a parallel pool

parfor i = 1:sample
    disp(i)

    theta_TN4 = M(i,1); theta_TN8 = M(i,2); theta_D0  = M(i,3); theta_M   = M(i,4);
    kappa_Th  = M(i,5); kappa_Tr  = M(i,6); kappa_Tc  = M(i,7);
    beta_Th   = M(i,8); beta_Tr   = M(i,9); beta_Tc   = M(i,10);

    main_parallel(i,t_start,t_end,theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc)

end

for i = 1:sample
    file = strcat('./VP_parallel/VP',num2str(i),'.mat');
    x=cell2mat(struct2cell(load(file)));
    cell(i,:) = x(t_end,:);
end

file=['./VP_immune/','VP_sample_immune.mat'];
save(file,'cell');

%% Data processing

Data_processing(sample);

%% VP Sample

select(sample);

%% Calculate Divergence

KL_JS_div();

%% Plot

picture_select();
picture_line();


