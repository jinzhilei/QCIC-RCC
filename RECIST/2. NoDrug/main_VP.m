%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script file is used to perform cell dynamics simulations of no drug

% t_start: start time (therapy)
% t_calibration: simulation time
% z: simulate the number of virtual patients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_start = 1;
t_calibration = 84;
t_end = 360;
dim = 6; 

param_bounds = [0.1,50;0.1,50;0.1,50;0.1,50;50,100;50,100];
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
[best_params] = PSO(Par,dim,param_bounds,t_start,t_calibration,pool);

file = './VP_immune/P';
save(file,'best_params');


    
    S = tumor_parameter(L, best_params);
    %% Simulation calculation
    futures = parallel.Future.empty();
    for j = 1:L
       
        futures(j) = parfeval(pool, @solve_one_patient, 1, ...
            j, ...
            Par,S,...
            t_start, t_end);  %time 
    end
    
    you = zeros(L, 28);  
    
    for completed = 1:L
        [idx, y_t] = fetchNext(futures);
        file = ['./VP_NoDrug/','VPsample_',num2str(idx),'.mat'];
        save(file, 'y_t');
    end
    for i =1:L
    file = ['./VP_NoDrug/','VPsample_',num2str(i),'.mat'];
    y=cell2mat(struct2cell(load(file)));
    you(i,:) = y(t_calibration,:);
    end

    e(L);


