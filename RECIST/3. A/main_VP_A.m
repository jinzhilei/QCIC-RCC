
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script file is used to perform cell dynamics simulations of atezolizumab

% a: start time (therapy)
% b: end time (therapy)
% c: therapy period 
% d: simulation time
% m: the number of therapy periods
% z: simulate the number of virtual patients
% Drug_Ate: therapy dose of atezolizumab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_start = 1;
t_calibration = 84;
t_end = 360;
dim = 10; 

param_bounds =[1,10;1,10;
    1,10;1,10;
    0.1,10;1,10;
    0.1,10;1,10;
    0.1,10;1,10];

dose_A = 1200;
v = 30;

%% Parameter sampling

file = '../1、VP_Immunogenomic/VP_immune/Number_VP_patients.mat';
O_order=cell2mat(struct2cell(load(file)));
L = length(O_order);

Tumor_0 = 3.2279e10+rand(L,1)*(2.0675e11 - 3.2279e10);
file_Tumor_0 = './VP_immune/Tumor_0';
save(file_Tumor_0,"Tumor_0");

Par = model_parameter_select(L,7);
file = './VP_immune/immune_par.mat';
save(file,"Par");

file = '../2. NoDrug/VP_immune/P';
P =cell2mat(struct2cell(load(file)));

Ate_par=PK_model_atezolizumab_parameter(L);
file_C='./VP_A/PK_model_atezolizumab_paramete.mat';
save(file_C,'Ate_par');

pool = gcp('nocreate');
if ~isempty(pool)
    delete(pool);
end
parpool(20);  % 使用20个worker
pool = gcp;
parfor i =1:pool.NumWorkers
    rng(1,'twister')
end


%% Calculation of tumor parameters
[best_params] = PSO(Par,P,Ate_par,dose_A,dim,param_bounds,t_start,t_calibration,pool);

file = './VP_immune/P_X2';
save(file,'best_params');

%% Calculation of tumor parameters

S = tumor_parameter(L, P,best_params);
    futures = parallel.Future.empty();
    for j = 1:L
        futures(j) = parfeval(pool, @solve_one_patient, 1, ...
            j, ...  
            Par,S,...
            Ate_par,dose_A,...
            t_start,t_calibration, t_end); 
    end
    you = zeros(L, 28);  
    
    for completed = 1:L
        [idx, y_t] = fetchNext(futures);
        file = ['./VP_A/','VPsample_',num2str(idx),'.mat'];
        save(file, 'y_t');
    end


e(L)


%% This function is used to describe the parameter sampling of the pharmacokinetics of atezolizumab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C=PK_model_atezolizumab_parameter(sample)

%% Generating sampling matrix

C = zeros(sample,6);

C(:,1)=Beta_distribution(2.83,3.36,5,3,sample);     
C(:,2)=Beta_distribution(2.49,3.63,2,1,sample);      
C(:,3)=Beta_distribution(0.546,0.714,8,2,sample);  
C(:,4)=Beta_distribution(0.2,0.274,9,1,sample);   
C(:,5)=Beta_distribution(50,7000,7,2,sample);       
C(:,6)=Beta_distribution(1000,13000,9,2,sample);      

end

