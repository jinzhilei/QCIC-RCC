function Adaptive_therapy()

a = 1;
d = 360;
c = [7,21];
z = 100;


%% Parameter sampling

file='./VP_patient/VP_model_parameter.mat';
A_immune_par=cell2mat(struct2cell(load(file)));

file = './VP_patient/Number_VP_patients.mat';  % VP_order_patient
O_order=cell2mat(struct2cell(load(file)));
L  = length(O_order);

%% Calculation of tumor parameters

S = tumor_parameter(z);
file_S='./VP_AB/VP_tumor_parameter.mat';
save(file_S,'S');

B_immune_par= A_immune_par(randi([1 L], z, 1),:);

%% Calculation of PK  parameters

C=PK_model_parameter(z);
file_C='./VP_AB/VP_PK_model_parameter.mat';
save(file_C,'C');

%% Initialized tumor cell count

r = 3.2279e10+rand(z,1) *(2.0675e11-3.2279e10);
file_R='./VP_AB/VP_tumor_number.mat';
save(file_R,'r');

%% Simulation calculation

delete (gcp('nocreate'))            % Close the current parallel computing pool
parpool(7);                         % Create a parallel pool
for j = 2
    Drug_day = c(j);
    for m = 1:3
        parfor i=1:z
            disp(i)
            % PK
            V_1_A = C(i,1);V_2_A = C(i,2);Q_A = C(i,3);CI_A = C(i,4);V_max_A = C(i,5);K_A_A = C(i,6);
            V_1_B = C(i,7);V_2_B = C(i,8);Q_B = C(i,9);CI_B = C(i,10);V_max_B = C(i,11);K_B_A = C(i,12);
            % immune
            theta_TN4 = B_immune_par(i,1);theta_TN8 = B_immune_par(i,2);theta_D0 = B_immune_par(i,3);
            theta_M = B_immune_par(i,4);kappa_Th = B_immune_par(i,5);kappa_Tr = B_immune_par(i,6);
            kappa_Tc = B_immune_par(i,7);beta_Th = B_immune_par(i,8);beta_Tr = B_immune_par(i,9);
            beta_Tc = B_immune_par(i,10);
            % tumor1
            G_X1 = S(i,1);beta_X1 = S(i,2);d_X1 = S(i,3);G_X2 = S(i,4);beta_X2 = S(i,5);
            d_X2 = S(i,6);Lambda_A = S(i,7);rho_X2X1 = S(i,8);K_A = S(i,9);
            % tumor2
            G_X3 = S(i,10);G_X4 = S(i,11);beta_X3 = S(i,12);beta_X4 = S(i,13);
            d_X3 = S(i,14);d_X4 = S(i,15);Lambda_B = S(i,16);rho_X3X1 = S(i,17);
            rho_X4X1 = S(i,18);rho_X4X2 = S(i,19);rho_X4X3 = S(i,20);K_B = S(i,21);

            treat_therapy(a,Drug_day,m,d,i,r,V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A,V_1_B,V_2_B,Q_B,CI_B,V_max_B,K_B_A,...
                theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,...
                d_X2,Lambda_A,rho_X2X1,K_A,G_X3,G_X4,beta_X3,beta_X4,d_X3,d_X4,Lambda_B,rho_X3X1,rho_X4X1,rho_X4X2,rho_X4X3,K_B)

        end
    end
end
end


function S = tumor_parameter(z)
S = zeros(z,21);

S(:,1) = Beta_distribution(5e12,5e13,9,3,z);       
S(:,2) = Beta_distribution(0.1,0.2,90,90,z);     
S(:,3) = Beta_distribution(0.03,0.15,50,20,z);      

S(:,4) = Beta_distribution(2.5e12,2.5e13,9,3,z);     
S(:,5) = Beta_distribution(0.05,0.15,0.71,0.7,z);    
S(:,6) = Beta_distribution(0.008,0.1,2,3.5,z);     
S(:,7) = Beta_distribution(0.02,5,0.9,1,z);                
S(:,8) = Beta_distribution(0.01,0.3,2,4,z);     
S(:,9) = Beta_distribution(3e2,1e3,3,8,z);    


S(:,10) = Beta_distribution(2.5e12,2.5e13,1,3,z);     
S(:,11) = Beta_distribution(2e12,2e13,1,3,z);      
S(:,12) = Beta_distribution(0.07,0.15,0.55,0.9,z);     
S(:,13) = Beta_distribution(0.03,0.13,1.1,0.9,z);    
S(:,14) = Beta_distribution(0.008,0.1,1.6,1,z);   
S(:,15) = Beta_distribution(0.005,0.08,4,2,z);   
S(:,16) = Beta_distribution(0.01,5,2,1,z);     
S(:,17) = Beta_distribution(0.01,0.3,1,4,z);      
S(:,18) = Beta_distribution(0.01,0.3,1,7,z);     
S(:,19) = Beta_distribution(0.01,0.3,1,7,z);     
S(:,20) = Beta_distribution(0.01,0.3,1,7,z);      
S(:,21) = Beta_distribution(1e2,8e2,8,3,z);     


end

function C=PK_model_parameter(sample)

%% Generating sampling matrix

C = zeros(sample,12);

C(:,1)=Beta_distribution(2.83,3.36,5,3,sample);       % V_1_A
C(:,2)=Beta_distribution(2.49,3.63,2,1,sample);       % V_2_A
C(:,3)=Beta_distribution(0.546,0.714,8,2,sample);     % Q_A
C(:,4)=Beta_distribution(0.2,0.274,9,1,sample);       % CI_A
C(:,5)=Beta_distribution(50,7000,7,2,sample);        % V_max_A
C(:,6)=Beta_distribution(1000,13000,9,2,sample);         % K_A

C(:,7)=Beta_distribution(2.57,4.20,2,2,sample);       % V_1_B
C(:,8)=Beta_distribution(2.29,7.70,2,6,sample);       % V_2_B
C(:,9)=Beta_distribution(0.27,0.69,6,2,sample);       % Q_B
C(:,10)=Beta_distribution(0.15,0.30,1,9,sample);       % CI_B
C(:,11)=Beta_distribution(120,800,8,2,sample);         % V_max_B
C(:,12)=Beta_distribution(300,1300,4,5,sample);        % K_B

end


%% This function is used to generate random numbers that follow a Beta distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=Beta_distribution(a,b,c,d,e)

y=a+(b-a)*betarnd(c,d,[1,e]);

end