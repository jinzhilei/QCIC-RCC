function treat_AB()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script file is used to perform the cell dynamics simulation of different drug treatments

% t_start: start time 
% t_calibration: simulation time
% t_end: end time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t_start = 1;
t_calibration = 84;
t_end = 360;

Drug = [1200;1050];
Drug_ratio = [5/10,1,15/10,2];
programme = [21,7,7];
ratio = [1,1,1;1/3,1/3,1/3;0.6,0.3,0.1];
[~,number_scheme] = size(programme);
number_drug = length(Drug_ratio);

L = 100;

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

file = '../4. AB/VP_immune/P_X34';
P_X34 =cell2mat(struct2cell(load(file)));

Dose_par=PK_model_parameter(L);

pool = gcp('nocreate');
if ~isempty(pool)
    delete(pool);
end
parpool(20);
pool = gcp;

for p = 1:number_scheme
    c = programme(1,p);
    for m = 1:number_drug
        for n = 1:number_drug
            Drug_Ate_init =  Drug(1)*Drug_ratio(m);
            Drug_Bev_init = Drug(2)*Drug_ratio(n);

            S = tumor_parameter(L, P_X1,P_X2,P_X34);
            file_S=strcat('treat_AB/treat',num2str(p),'/(',num2str(m),',',num2str(n),')','_VP_tumor_parameter.mat');
            save(file_S,'S');

            futures = parallel.Future.empty();
            for j = 1:L

                futures(j) = parfeval(pool, @solve_one_patient, 1, ...
                    j, ratio,c,p,...  
                    Par,S,...
                    Dose_par,Drug_Ate_init,Drug_Bev_init,...
                    t_start,t_calibration, t_end);  
            end

            for completed = 1:L
                [idx, y_t] = fetchNext(futures);
                file=strcat('treat_AB/treat',num2str(p),'/(',num2str(m),',',num2str(n),')','_VPsample_AB',num2str(idx),'.mat');
                save(file,'y_t');
            end
        end
    end
end

end

%%
function C=PK_model_parameter(sample)

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