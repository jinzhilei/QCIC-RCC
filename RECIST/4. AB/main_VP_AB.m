
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
dim = 22; 


param_bounds =[0.1,10;0.1,10;
    0.1,10;0.1,10;
    0.1,10;0.1,10;
    0.1,10;0.1,10;
    0.1,10;0.1,10;
    0.1,10;0.1,10;
    0.1,10;0.1,10;
    0.1,10;0.1,10;
    0.1,10;0.1,10;
    0.1,10;0.1,10;
    0.1,10;0.1,10];

dose_A = 1200;
dose_B = 1050;

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
P_X1 =cell2mat(struct2cell(load(file)));

file = '../3. A/VP_immune/P_X2';
P_X2 =cell2mat(struct2cell(load(file)));

Dose_par=PK_model_parameter(L);
file_C='./VP_immune/PK_model_paramete.mat';
save(file_C,'Dose_par');


pool = gcp('nocreate');
if ~isempty(pool)
    delete(pool);
end
parpool(20);  
pool = gcp;
parfor i =1:pool.NumWorkers
    rng(1,'twister')
end


%% Calculation of tumor parameters
[best_params] = PSO(Par,P_X1,P_X2,Dose_par,dose_A,dose_B,dim,param_bounds,t_start,t_calibration,pool);

file = './VP_immune/P_X34';%shiyan  _2  woziji meigang  monijieguo _3zuiyou
save(file,'best_params');
% file = './VP_immune/P_X34';
% best_params=cell2mat(struct2cell(load(file)));

%% Calculation of tumor parameters

S = tumor_parameter(L, P_X1,P_X2,best_params);
    futures = parallel.Future.empty();
    for j = 1:L
        futures(j) = parfeval(pool, @solve_one_patient, 1, ...
            j, ...  
            Par,S,...
            Dose_par,dose_A,dose_B,...
            t_start,t_calibration, t_end);  
    end
    
    you = zeros(L, 28);  
    
    for completed = 1:L
        [idx, y_t] = fetchNext(futures);
        
        file = ['./VP_AB/','VPsample_',num2str(idx),'.mat'];
        save(file, 'y_t');
    end


e(L)


%% This function is used to describe the parameter sampling of the pharmacokinetics of atezolizumab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C=PK_model_parameter(sample)

%% Generating sampling matrix

C = zeros(sample,12);

C(:,1)=Beta_distribution(2.83,3.36,5,3,sample);       
C(:,2)=Beta_distribution(2.49,3.63,2,1,sample);       
C(:,3)=Beta_distribution(0.546,0.714,8,2,sample);  
C(:,4)=Beta_distribution(0.2,0.274,9,1,sample);      
C(:,5)=Beta_distribution(50,7000,7,2,sample);       
C(:,6)=Beta_distribution(1000,13000,9,2,sample);       

C(:,7)=Beta_distribution(2.57,4.20,2,2,sample);  
C(:,8)=Beta_distribution(2.29,7.70,2,6,sample);      
C(:,9)=Beta_distribution(0.27,0.69,6,2,sample);     
C(:,10)=Beta_distribution(0.15,0.30,1,9,sample);      
C(:,11)=Beta_distribution(120,800,8,2,sample);      
C(:,12)=Beta_distribution(300,1300,4,5,sample);     


end

