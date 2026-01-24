
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script file is used to perform cell dynamics simulations of atezolizumab

% a: start time (therapy)
% d: simulation time
% m: the number of therapy periods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 1;
b = 360;
c = 21;
d = 360;
Drug_Ate = 1200;
z = 10000;
v = 30;

%% Parameter sampling

file='./VP_patient/VP_model_parameter.mat';
A_immune_par=cell2mat(struct2cell(load(file)));

file = './VP_patient/Number_VP_patients.mat';  % VP_order_patient
O_order=cell2mat(struct2cell(load(file)));
L  = length(O_order);

tout = zeros(1,d);
yout = zeros(z,v);

S = tumor_parameter(z);
file_S='./VP_A_bar/VP_tumor_parameter.mat';
save(file_S,'S');

B_immune_par= A_immune_par(randi([1 L], z, 1),:);

%% Simulation calculation

Ate_par=PK_model_atezolizumab_parameter(z);
file_C='./VP_A_bar/PK_model_atezolizumab_paramete.mat';
save(file_C,'Ate_par');

r = randi([3.2279e10, 2.0675e11], 1, z);
file_R='./VP_A_bar/VP_tumor_initialzation.mat';
save(file_R,'r');

for i=1:z
    disp(i)
    theta_TN4  =  B_immune_par(i,1);
    theta_TN8  =  B_immune_par(i,2);
    theta_D0     =  B_immune_par(i,3);
    theta_M       =  B_immune_par(i,4);
    kappa_Th    =  B_immune_par(i,5);
    kappa_Tr     =  B_immune_par(i,6);
    kappa_Tc    =  B_immune_par(i,7);
    beta_Th       =  B_immune_par(i,8);
    beta_Tr        =  B_immune_par(i,9);
    beta_Tc       =  B_immune_par(i,10);

    G_X1               =  S(i,1);
    beta_X1          =  S(i,2);
    d_X1                =  S(i,3);
    G_X2               =  S(i,4);
    beta_X2          =  S(i,5);
    d_X2                =  S(i,6);
    Lambda_A     =  S(i,7);
    rho_X2X1       =  S(i,8);
    K_A                 =  S(i,9);

    V_1_A       =  Ate_par(i,1);
    V_2_A       =  Ate_par(i,2);
    Q_A           =  Ate_par(i,3);
    CI_A          =  Ate_par(i,4);
    V_max_A  =  Ate_par(i,5);
    K_A_A      =  Ate_par(i,6);


    t_start = a;
    t_tau   = c;
    t_end  = t_start+t_tau;
    t_calibration = b;
    t_simend = d;
    m=ceil(t_calibration/t_tau);

    Q_Ate_start = [Drug_Ate 0];
    Q_Ate = zeros(1,30);
    Q_Ate(29) = Drug_Ate;

    x0=initialzation_VP_A(Q_Ate_start,r(i));

    for j = 1:m
        if j <m
            [t,y] = ode45(@(t,x) QCIC_VP_A(t,x,theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,d_X2,...
                Lambda_A,rho_X2X1,K_A, V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A),t_start:1:t_end, x0);
            nt = length(t);
            tout(t_start:1:t_end) = t_start:1:t_end;
            you(t_start:1:t_end,:) = y(1:nt,:);
            t_start = t_end;
            t_end = t_start + t_tau;
            x0 = Q_Ate+y(end,:);
        else
            [t,y] = ode45(@(t,x) QCIC_VP_A(t,x,theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,d_X2,...
                Lambda_A,rho_X2X1,K_A, V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A),t_start:1:t_simend, x0);
            nt = length(t);
            tout(t_start:1:t_simend) = t_start:1:t_simend;
            you(t_start:1:t_simend,:) = y(1:nt,:);
        end
    end
    [~,position] = min(you(43:360,25)+you(43:360,26)+you(43:360,27)+you(43:360,28));
    yout(i,:) = you(position+42,:);
end

file=['./VP_A_bar/','VPsample.mat'];
save(file,'yout');

counts();

%% This function is used to describe the parameter sampling of VP


%% This function is used to describe the parameter sampling of VP

function S = tumor_parameter(z)

S = zeros(z,9);

S(:,1) = Beta_distribution(5e12,5e13,9,3,z);       % G_X1         %1.5e9
S(:,2) = Beta_distribution(0.1,0.2,90,90,z);      % \beta_X1     %0.15
S(:,3) = Beta_distribution(0.03,0.15,50,20,z);      % \d_X1       %0.04

S(:,4) = Beta_distribution(2.5e12,2.5e13,9,3,z);       % G_X2         %8e8
S(:,5) = Beta_distribution(0.05,0.15,0.71,0.7,z);     % \beta_X2     %0.1
S(:,6) = Beta_distribution(0.008,0.1,2,3.5,z);     % \d_X2        %0.02
S(:,7) = Beta_distribution(0.02,5,0.9,1,z);                  % \Lambda_A    %8e-3
S(:,8) = Beta_distribution(0.01,0.3,2,4,z);      % \rho_X2X1    % 0.05
S(:,9) = Beta_distribution(3e2,1e3,3,8,z);       % K_A         % 转化 8e2

end

%% This function is used to describe the parameter sampling of the pharmacokinetics of atezolizumab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C=PK_model_atezolizumab_parameter(sample)

%% Generating sampling matrix

C = zeros(sample,6);

C(:,1)=Beta_distribution(2.83,3.36,5,3,sample);       % V_1_A
C(:,2)=Beta_distribution(2.49,3.63,2,1,sample);       % V_2_A
C(:,3)=Beta_distribution(0.546,0.714,8,2,sample);     % Q_A
C(:,4)=Beta_distribution(0.2,0.274,9,1,sample);       % CI_A
C(:,5)=Beta_distribution(50,7000,7,2,sample);        % V_max_A
C(:,6)=Beta_distribution(1000,13000,9,2,sample);         % K_A

end

%% This function is used to generate random numbers that follow a Beta distribution

function y=Beta_distribution(a,b,c,d,e)

y=a+(b-a)*betarnd(c,d,[1,e]);

end